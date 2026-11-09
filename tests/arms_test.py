"""Tests for regions with their own plugins: channel assignment, configuration
and the worker processes of readfish._arms (with the plugins in arms_plugins.py
and a fake Read Until client)."""

import logging
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pytest
import rtoml
from minknow_api import data_pb2

from readfish._arms import ArmManager
from readfish._config import Conf

if sys.version_info < (3, 11):
    from exceptiongroup import BaseExceptionGroup

TESTS = Path(__file__).parent
# the spawned workers import arms_plugins from here
sys.path.insert(0, str(TESTS))

CHANNELS = 12
DECISIONS = """
targets = ["target"]
single_on = "stop_receiving"
multi_on = "stop_receiving"
single_off = "unblock"
multi_off = "unblock"
no_seq = "proceed"
no_map = "unblock"
min_chunks = 1
max_chunks = 2
"""
CONTROL = """
control = true
targets = []
single_on = "stop_receiving"
multi_on = "stop_receiving"
single_off = "stop_receiving"
multi_off = "stop_receiving"
no_seq = "stop_receiving"
no_map = "stop_receiving"
"""


def make_toml(fast_delay=0.0, slow_delay=0.0, fail=False):
    return f"""
channel_assignment = "modulo"

[caller_settings.arms_plugins]

[[regions]]
name = "fast"
{DECISIONS}
mapper_settings.arms_plugins = {{ motif = "ACGT", delay = {fast_delay}, fail = {str(fail).lower()} }}

[[regions]]
name = "slow"
{DECISIONS}
mapper_settings.arms_plugins = {{ motif = "TTTT", delay = {slow_delay} }}

[[regions]]
name = "control"
{CONTROL}
"""


def test_modulo_assignment():
    conf = Conf.from_dict(rtoml.loads(make_toml()), CHANNELS)
    assert conf.parallel
    assert [conf.get_region(c).name for c in range(1, 7)] == [
        "fast",
        "slow",
        "control",
        "fast",
        "slow",
        "control",
    ]
    assert conf.region_channels(0) == [1, 4, 7, 10]
    assert sorted(conf.plugin_settings()) == ["fast", "slow"]
    assert "every 3th starting at channel 2" in conf.describe_experiment()


def test_split_assignment_unchanged():
    # without channel_assignment, regions are blocks of the flow cell as before
    toml = make_toml().replace('channel_assignment = "modulo"', "")
    toml = toml[: toml.index('[[regions]]\nname = "control"')]
    conf = Conf.from_dict(rtoml.loads(toml), 512)
    assert conf.channel_assignment == "split"
    assert len(conf.region_channels(0)) == len(conf.region_channels(1)) == 256
    assert conf.region_channels(0) != list(range(1, 513, 2))


def test_region_plugins_fall_back_to_top_level():
    toml = make_toml().replace(
        "[caller_settings.arms_plugins]",
        "[caller_settings.arms_plugins]\n[mapper_settings.arms_plugins]\nmotif = 'GGGG'",
    )
    toml = toml.replace(
        'mapper_settings.arms_plugins = { motif = "TTTT", delay = 0.0 }', ""
    )
    conf = Conf.from_dict(rtoml.loads(toml), CHANNELS)
    settings = conf.plugin_settings()
    assert settings["slow"][1].parameters == {"motif": "GGGG"}
    assert settings["fast"][1].parameters["motif"] == "ACGT"


def test_missing_plugins_rejected():
    toml = make_toml().replace("[caller_settings.arms_plugins]", "")
    # cattrs wraps errors raised while building the Conf in an exception group
    with pytest.raises(BaseExceptionGroup) as error:
        Conf.from_dict(rtoml.loads(toml), CHANNELS)
    assert "No caller_settings or mapper_settings for 'fast'" in str(
        error.value.exceptions
    )


def test_invalid_channel_assignment():
    toml = make_toml().replace('"modulo"', '"zigzag"')
    with pytest.raises(Exception):
        Conf.from_dict(rtoml.loads(toml), CHANNELS)


def test_to_file_roundtrip(tmp_path):
    conf = Conf.from_dict(rtoml.loads(make_toml()), CHANNELS)
    conf.to_file(tmp_path / "out.toml")
    assert Conf.from_file(tmp_path / "out.toml", CHANNELS) == conf


class FakeClient:
    """Records the actions the ArmManager forwards, with their arrival times."""

    def __init__(self):
        self.lock = threading.Lock()
        self.actions = []

    def unblock_read_batch(self, reads, duration=0.1):
        with self.lock:
            self.actions += [("unblock", ch, rid, time.time()) for ch, rid in reads]

    def stop_receiving_batch(self, reads):
        with self.lock:
            self.actions += [("stop", ch, rid, time.time()) for ch, rid in reads]

    def wait_for(self, n, timeout=60):
        deadline = time.time() + timeout
        while time.time() < deadline:
            with self.lock:
                if len(self.actions) >= n:
                    return list(self.actions)
            time.sleep(0.01)
        raise TimeoutError(f"got {len(self.actions)} of {n} actions")


def make_manager(toml_text, client, dry_run=False):
    conf = Conf.from_dict(rtoml.loads(toml_text), CHANNELS)
    return ArmManager(
        client,
        conf,
        toml_text,
        logging.getLogger("arms_test"),
        0.1,
        channels=CHANNELS,
        padding=0,
        run_information=None,
        sample_rate=4000.0,
        signal_dtype=np.dtype("int16"),
        calibration_values={c: (1.0, 0.0) for c in range(1, CHANNELS + 1)},
        dry_run=dry_run,
        chemistry="simplex",
        run_id="test",
        debug_log=True,
        break_reads_after_seconds=1.0,
    )


def chunk(channel, seq, read_id=None):
    return (
        channel,
        data_pb2.GetLiveReadsResponse.ReadData(
            id=read_id or f"read{channel}", raw_data=seq.encode()
        ),
    )


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_arms_decide_independently(workdir):
    client = FakeClient()
    manager = make_manager(make_toml(), client)
    manager.start()
    try:
        # channels 1, 4, 7, 10 -> fast (ACGT); 2, 5, 8, 11 -> slow (TTTT); 3, 6, ... control
        chunks = [chunk(c, "ACGTAAAA" if c % 2 else "TTTTAAAA") for c in range(1, 13)]
        control = manager.dispatch(chunks, started_during_sequencing=False)
        assert [c for c, _ in control] == [3, 6, 9, 12]
        actions = client.wait_for(8)
    finally:
        manager.stop()
    decided = {ch: kind for kind, ch, _, _ in actions}
    # fast keeps the ACGT reads (odd channels), slow the TTTT reads (even channels)
    assert decided == {
        1: "stop",
        4: "unblock",
        7: "stop",
        10: "unblock",
        2: "stop",
        5: "unblock",
        8: "stop",
        11: "unblock",
    }
    # each worker writes its debug log and batch log
    for name in ("fast", "slow"):
        assert (workdir / f"test_readfish_{name}.tsv").exists()
        batches = (
            (workdir / f"test_readfish_{name}_batches.tsv").read_text().splitlines()
        )
        assert batches[0].startswith("batch\tn_reads")
        assert batches[1].split("\t")[1] == "4"


def test_slow_arm_does_not_delay_fast_arm(workdir):
    client = FakeClient()
    manager = make_manager(make_toml(slow_delay=2.0), client)
    manager.start()
    try:
        t0 = time.time()
        manager.dispatch(
            [chunk(c, "ACGT") for c in range(1, 13)], started_during_sequencing=False
        )
        actions = client.wait_for(8, timeout=30)
    finally:
        manager.stop()
    fast = [t - t0 for _, ch, _, t in actions if ch % 3 == 1]
    slow = [t - t0 for _, ch, _, t in actions if ch % 3 == 2]
    assert max(fast) < 1.0
    assert min(slow) >= 2.0


def test_newest_chunk_per_channel_and_first_read_override(workdir):
    client = FakeClient()
    manager = make_manager(make_toml(slow_delay=1.0), client)
    manager.start()
    try:
        # the slow worker is busy with the first batch while the next two arrive;
        #   it then decides on the newest chunk of channel 2 only
        manager.dispatch([chunk(2, "TTTT", "r1")], started_during_sequencing=True)
        time.sleep(0.2)
        manager.dispatch([chunk(2, "AAAA", "r2")], started_during_sequencing=False)
        manager.dispatch([chunk(2, "AAAATTTT", "r2")], started_during_sequencing=False)
        client.wait_for(2, timeout=30)
        time.sleep(1.5)
    finally:
        manager.stop()
    assert [(kind, rid) for kind, _, rid, _ in client.actions] == [
        # first read on the channel while started during sequencing: sequenced
        ("stop", "r1"),
        # newest chunk of r2 contains TTTT
        ("stop", "r2"),
    ]
    batches = (workdir / "test_readfish_slow_batches.tsv").read_text().splitlines()
    assert [line.split("\t")[1:3] for line in batches[1:]] == [["1", "1"], ["1", "2"]]


def test_dry_run_sends_no_unblocks(workdir):
    client = FakeClient()
    manager = make_manager(make_toml(), client, dry_run=True)
    manager.start()
    try:
        manager.dispatch([chunk(1, "AAAA")], started_during_sequencing=False)
        actions = client.wait_for(1)
    finally:
        manager.stop()
    assert [kind for kind, *_ in actions] == ["stop"]
    assert "unblock" in (workdir / "test_readfish_fast.tsv").read_text()


def test_worker_start_failure(workdir):
    client = FakeClient()
    manager = make_manager(make_toml(fail=True), client)
    with pytest.raises(RuntimeError, match="aligner failed on purpose"):
        manager.start()
