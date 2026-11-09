"""Run regions with their own caller and aligner in parallel worker processes.

When regions of a TOML bring their own ``caller_settings``/``mapper_settings``
(see :class:`readfish._config.Region`), each non-control region becomes an
*arm*: a worker process with its own caller, aligner and decision state, so
that several selection algorithms can be compared on one flow cell. The main
process keeps the single Read Until connection:

#. :class:`ArmManager` starts the workers and waits until their plugins are
   initialised.
#. Every iteration, the main process forwards each arm the newest chunk of
   each of its channels (:meth:`ArmManager.dispatch`). Control regions are
   answered with ``stop_receiving`` directly and never base-called.
#. A worker (:func:`arm_worker`) collects everything forwarded since its last
   batch, keeps the newest chunk per channel, base-calls, aligns and decides
   exactly as the single-process loop does, and sends its actions back.
#. A thread in the main process passes the actions on to the Read Until client
   as they arrive.

So an arm never waits for another arm's base-calls or alignments: a slow
aligner only delays the decisions on its own channels.

Each worker writes, into the working directory:
    - ``<run_id>_readfish_<region>.tsv``: the per-read debug log, as ``readfish targets``
      writes for the whole flow cell (if the debug log is enabled)
    - ``<run_id>_readfish_<region>_batches.tsv``: per batch the number of reads, the
      time the oldest chunk waited before the worker picked it up, and the time spent
      base-calling, aligning and deciding
    - ``<run_id>_readfish_<region>.log``: the worker's log
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import queue
import threading
import time
import traceback
from pathlib import Path
from typing import Any, Optional

import attrs
import rtoml

from readfish._config import Conf
from readfish._loggers import flush_queued_logs, setup_logger

# How long the main process waits for workers to load their plugins (seconds)
STARTUP_TIMEOUT = 1800
BATCH_LOG_FIELDS = (
    "batch",
    "n_reads",
    "n_chunks_forwarded",
    "wait_s",
    "basecall_s",
    "align_s",
    "decide_s",
    "total_s",
    "timestamp",
)


@attrs.define
class ArmSpec:
    """Everything a worker process needs, picklable for ``spawn``.

    :param name: Name of the region
    :param region_index: Index of the region in ``Conf.regions``
    :param toml_text: The TOML text of the experiment
    :param channels: Number of channels on the flow cell
    :param padding: Target padding given on the command line
    :param run_information: Serialised ``ProtocolRunInfo`` of the run, if connected
    :param sample_rate: Sample rate of the device
    :param signal_dtype: Numpy dtype of the raw signal
    :param calibration_values: Channel to (scaling, offset), as ``ReadUntilClient.calibration_values``
    :param dry_run: Replace unblocks with ``stop_receiving``
    :param chemistry: Value of the :class:`readfish._cli_args.Chemistry` enum
    :param run_id: Run id, used to name the worker's files
    :param debug_log: Write the per-read debug log
    :param break_reads_after_seconds: MinKNOW's chunk length, for the statistics
    """

    name: str
    region_index: int
    toml_text: str
    channels: int
    padding: int
    run_information: Optional[bytes]
    sample_rate: float
    signal_dtype: Any
    calibration_values: dict
    dry_run: bool
    chemistry: str
    run_id: str
    debug_log: bool
    break_reads_after_seconds: float


def load_conf(toml_text: str, channels: int, padding: int = 0) -> Conf:
    """Parse the experiment TOML as ``readfish targets`` does, including padding."""
    conf = Conf.from_dict(rtoml.loads(toml_text), channels)
    if padding:
        for region in conf.regions:
            region.targets.padding = padding
    return conf


def _newest_per_channel(messages: list) -> tuple[list, float, int]:
    """Merge forwarded chunk lists, keeping the newest chunk per channel.

    Chunks of one read accumulate in the Read Until cache, so a newer chunk of a
    read contains the older ones. A newer read on a channel replaces the older read.

    :return: (channel, serialised ReadData) pairs, the forwarding time of the oldest
        message, and the number of chunks forwarded
    """
    newest, oldest, n = {}, None, 0
    for sent_at, chunks in messages:
        oldest = sent_at if oldest is None else min(oldest, sent_at)
        for channel, data in chunks:
            newest[channel] = data
            n += 1
    return list(newest.items()), oldest, n


def arm_worker(spec: ArmSpec, inbox: mp.Queue, outbox: mp.Queue) -> None:
    """Entry point of a worker process; see the module documentation.

    Messages received: ``("chunks", started_during_sequencing, sent_at, [(channel, bytes)])``
    or ``None`` to stop. Messages sent: ``("ready", name, description)``,
    ``("error", name, traceback)``, ``("actions", name, unblock, stop_receiving)``
    and ``("stopped", name, summary)``.
    """
    # Imported here: the module is also imported by the main process, which
    #   does not need minknow_api's protobufs or the decision logic
    from minknow_api import data_pb2, protocol_pb2
    from readfish.read_until.base import CALIBRATION
    from readfish.entry_points.targets import DecisionMaker
    from readfish._cli_args import Chemistry
    from readfish._statistics import ReadfishStatistics
    from readfish._utils import ChunkTracker
    from readfish.plugins.utils import PreviouslySentActionTracker, DuplexTracker

    prefix = f"{spec.run_id}_readfish_{spec.name}"
    logger = setup_logger(
        f"readfish.arm.{spec.name}",
        log_file=f"{prefix}.log",
        log_format="%(asctime)s %(name)s %(message)s",
        level=logging.INFO,
    )
    caller = mapper = None
    try:
        conf = load_conf(spec.toml_text, spec.channels, spec.padding)
        caller_settings, mapper_settings = conf.plugin_settings()[spec.name]
        run_information = None
        if spec.run_information is not None:
            run_information = protocol_pb2.ProtocolRunInfo.FromString(
                spec.run_information
            )
        logger.info(f"Initialising caller {caller_settings.name}")
        caller = caller_settings.load_object(
            "Caller", run_information=run_information, sample_rate=spec.sample_rate
        )
        logger.info(f"Initialising aligner {mapper_settings.name}")
        mapper = mapper_settings.load_object("Aligner")
        description = "\n".join(
            [caller.describe(), mapper.describe(conf.regions, conf.barcodes)]
        )
        logger.info(description)

        state = DecisionMaker()
        state.conf = conf
        state.logger = logger
        state.chemistry = Chemistry(spec.chemistry)
        state.dry_run = spec.dry_run
        state.readfish_started_during_sequencing = True
        state.chunk_tracker = ChunkTracker(spec.channels)
        state.previous_action_tracker = PreviouslySentActionTracker()
        state.duplex_tracker = DuplexTracker()
        state.loop_statistics = ReadfishStatistics(
            f"{prefix}.tsv" if spec.debug_log else None,
            spec.break_reads_after_seconds,
        )
        calibration = {
            ch: CALIBRATION(*values) for ch, values in spec.calibration_values.items()
        }
    except Exception:
        outbox.put(("error", spec.name, traceback.format_exc()))
        for plugin in (caller, mapper):
            if plugin is not None:
                plugin.disconnect()
        flush_queued_logs()
        return
    outbox.put(("ready", spec.name, description))

    batch_log = open(f"{prefix}_batches.tsv", "w", buffering=1)
    batch_log.write("\t".join(BATCH_LOG_FIELDS) + "\n")
    loop_counter = 0
    try:
        while True:
            message = inbox.get()
            if message is None:
                break
            messages = [message]
            # Everything forwarded while the last batch was processed
            stop = False
            while True:
                try:
                    message = inbox.get_nowait()
                except queue.Empty:
                    break
                if message is None:
                    stop = True
                    break
                messages.append(message)
            # The flag of the newest message is the current one
            state.readfish_started_during_sequencing = messages[-1][1]
            chunks, sent_at, n_forwarded = _newest_per_channel(
                [(m[2], m[3]) for m in messages]
            )
            if chunks:
                t0 = time.time()
                reads = [
                    (channel, data_pb2.GetLiveReadsResponse.ReadData.FromString(data))
                    for channel, data in chunks
                ]
                calls = list(caller.basecall(reads, spec.signal_dtype, calibration))
                t1 = time.time()
                aligns = list(mapper.map_reads(calls))
                t2 = time.time()
                loop_counter += 1
                unblock, stop_receiving = [], []
                for n, result in enumerate(aligns, start=1):
                    state.process_result(
                        result, loop_counter, n, stop_receiving, unblock
                    )
                outbox.put(("actions", spec.name, unblock, stop_receiving))
                t3 = time.time()
                state.loop_statistics.add_batch_performance(
                    number_of_reads=len(aligns), batch_time=t3 - t0
                )
                batch_log.write(
                    "\t".join(
                        map(
                            str,
                            (
                                loop_counter,
                                len(aligns),
                                n_forwarded,
                                f"{t0 - sent_at:.4f}",
                                f"{t1 - t0:.4f}",
                                f"{t2 - t1:.4f}",
                                f"{t3 - t2:.4f}",
                                f"{t3 - t0:.4f}",
                                f"{t3:.3f}",
                            ),
                        )
                    )
                    + "\n"
                )
                if loop_counter % 100 == 0:
                    logger.info(state.loop_statistics.get_batch_performance())
            if stop:
                break
    except Exception:
        outbox.put(("error", spec.name, traceback.format_exc()))
    finally:
        batch_log.close()
        summary = state.loop_statistics.get_batch_performance()
        logger.info(summary)
        caller.disconnect()
        mapper.disconnect()
        flush_queued_logs()
        outbox.put(("stopped", spec.name, summary))


class ArmManager:
    """Starts and feeds the worker processes, and forwards their actions.

    :param client: The Read Until client (an ``RUClient``)
    :param conf: The experiment configuration, with ``conf.parallel`` True
    :param toml_text: The TOML text the configuration was read from
    :param logger: Logger of the main process
    :param unblock_duration: Time, in seconds, to apply unblock voltage
    :param specs: Keyword arguments for :class:`ArmSpec` shared by all arms
    """

    def __init__(
        self,
        client,
        conf: Conf,
        toml_text: str,
        logger: logging.Logger,
        unblock_duration: float,
        **specs,
    ):
        self.client = client
        self.conf = conf
        self.logger = logger
        self.unblock_duration = unblock_duration
        # spawn, not fork: the main process runs gRPC threads
        ctx = mp.get_context("spawn")
        self.outbox = ctx.Queue()
        self.inboxes: dict[str, mp.Queue] = {}
        self.processes: dict[str, mp.Process] = {}
        # region index of each arm; control regions have no arm
        self.arm_of_region: dict[int, str] = {}
        for index, region in enumerate(conf.regions):
            if region.control:
                continue
            spec = ArmSpec(
                name=region.name, region_index=index, toml_text=toml_text, **specs
            )
            self.inboxes[region.name] = ctx.Queue()
            self.processes[region.name] = ctx.Process(
                target=arm_worker,
                args=(spec, self.inboxes[region.name], self.outbox),
                name=f"readfish-arm-{region.name}",
                daemon=True,
            )
            self.arm_of_region[index] = region.name
        self.failed: set[str] = set()
        self._forwarder: Optional[threading.Thread] = None

    def start(self) -> None:
        """Start the workers and wait until all have initialised their plugins.

        :raises RuntimeError: If a worker fails to initialise
        """
        for process in self.processes.values():
            process.start()
        waiting = set(self.processes)
        deadline = time.time() + STARTUP_TIMEOUT
        while waiting:
            try:
                kind, name, payload = self.outbox.get(timeout=1)[:3]
            except queue.Empty:
                for name in list(waiting):
                    if not self.processes[name].is_alive():
                        self.stop()
                        raise RuntimeError(
                            f"Worker for region {name!r} exited during start-up"
                        )
                if time.time() > deadline:
                    self.stop()
                    raise RuntimeError(
                        f"Workers {sorted(waiting)} did not start in time"
                    )
                continue
            if kind == "error":
                self.stop()
                raise RuntimeError(
                    f"Worker for region {name!r} failed to start:\n{payload}"
                )
            if kind == "ready":
                self.logger.info(f"Region {name} worker ready:\n{payload}")
                waiting.discard(name)
        self._forwarder = threading.Thread(
            target=self._forward_actions, name="readfish-arm-actions", daemon=True
        )
        self._forwarder.start()

    def _forward_actions(self) -> None:
        """Pass the workers' actions to the Read Until client as they arrive."""
        while True:
            message = self.outbox.get()
            if message is None:
                return
            kind, name = message[:2]
            if kind == "actions":
                unblock, stop_receiving = message[2:]
                if unblock:
                    self.client.unblock_read_batch(
                        unblock, duration=self.unblock_duration
                    )
                if stop_receiving:
                    self.client.stop_receiving_batch(stop_receiving)
            elif kind == "error":
                self.failed.add(name)
                self.logger.error(
                    f"Worker for region {name!r} failed, its channels get no more decisions:\n{message[2]}"
                )
            elif kind == "stopped":
                self.logger.info(f"Region {name} worker stopped: {message[2]}")

    def dispatch(self, chunks: list, started_during_sequencing: bool) -> list:
        """Forward chunks to the workers of their regions.

        :param chunks: ``(channel, ReadData)`` pairs from the Read Until client
        :param started_during_sequencing: Whether readfish started during sequencing
        :return: The ``(channel, ReadData)`` pairs of control regions, which are
            not forwarded
        """
        per_arm: dict[str, list] = {}
        control = []
        for channel, read in chunks:
            name = self.arm_of_region.get(self.conf._channel_map[channel])
            if name is None:
                control.append((channel, read))
            elif name not in self.failed:
                per_arm.setdefault(name, []).append((channel, read.SerializeToString()))
        sent_at = time.time()
        for name, arm_chunks in per_arm.items():
            self.inboxes[name].put(
                ("chunks", started_during_sequencing, sent_at, arm_chunks)
            )
        return control

    def stop(self, timeout: float = 30) -> None:
        """Stop the workers and the action forwarding."""
        for name, inbox in self.inboxes.items():
            if self.processes[name].is_alive():
                inbox.put(None)
        for name, process in self.processes.items():
            if process.pid is not None:
                process.join(timeout)
                if process.is_alive():
                    self.logger.warning(
                        f"Worker for region {name!r} did not stop, terminating"
                    )
                    process.terminate()
        if self._forwarder is not None:
            self.outbox.put(None)
            self._forwarder.join(timeout)
            self._forwarder = None
