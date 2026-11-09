"""Validate experiment configuration :doc:`TOML <toml>` files.

This script is used to check that an experiment configuration file can be loaded.
By default, this will attempt to load the ``Caller`` and ``Aligner`` plugins as specified in the
:doc:`TOML <toml>` file. The ``--no-check-plugins`` flag can be used to skip this step.
The ``--no-describe`` flag can be used to skip printing a description of the configuration to terminal. Description will error if there a re contigs listed as targets which are not found in the reference.

These basic checks are for compatibility and do not indicate that a configuration/plugins will work efficiently with readfish.

Any errors are passed back through and printed to terminal.
If you require help understanding or resolving an error you can check the :doc:`TOML <toml>` documentation pages or `open an issue`_.

Example run command::

   readfish validate my_exp.toml

.. _`open an issue`: https://github.com/LooseLab/readfish/issues/new
"""

import sys
import logging
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor

if sys.version_info < (3, 11):
    from exceptiongroup import BaseExceptionGroup

from readfish._config import Conf
from readfish._cli_args import BASE_ARGS
from readfish._utils import iter_exception_group

_help = "Readfish TOML Validator"
_cli = BASE_ARGS + (
    ("toml", dict(help="TOML file to validate")),
    (
        "--prom",
        dict(
            help="Use this flag if the target platform is a PromethION",
            action="store_true",
        ),
    ),
    (
        "--no-check-plugins",
        dict(
            help="Do not attempt to load the plugins with the supplied toml file configuration.",
            action="store_true",
        ),
    ),
    (
        "--no-describe",
        dict(
            help="Do not describe the experimental and plugin configuration. Always disabled if --no-check-plugins is passed.",
            action="store_true",
        ),
    ),
)


def _check_plugins(
    toml: str, channels: int, name: str
) -> tuple[str | None, str | None, str | None]:
    """Initialise the caller and aligner of one condition.

    :param toml: Path to the TOML file
    :param channels: Number of channels on the flow cell
    :param name: Key of the condition in :meth:`Conf.plugin_settings`
    :return: The caller's error, the aligner's error, and the aligner's
        description (None where there is none)
    """
    conf = Conf.from_file(toml, channels)
    caller_settings, mapper_settings = conf.plugin_settings()[name]
    caller_error = aligner_error = description = None
    try:
        _ = caller_settings.load_object("Caller")
    except Exception as exc:
        caller_error = str(exc)
    try:
        al = mapper_settings.load_object("Aligner")
    except Exception as exc:
        aligner_error = str(exc)
    else:
        try:
            description = al.describe(conf.regions, conf.barcodes)
        except Exception as exc:
            description = f"Aligner description failed: {exc}"
    return caller_error, aligner_error, description


def run(parser, args, extras) -> int:
    """
    Runs the validate entry point

    :parser: The parser. Unused, but must be included as targets.py requires it
    :type:
    """
    logger = logging.getLogger(f"readfish.{args.command}")
    channels = 3000 if args.prom else 512
    try:
        conf = Conf.from_file(args.toml, channels, logger=logger)
    except BaseExceptionGroup as exc:
        msg = f"Could not load TOML config ({args.toml}), see below for details:\n"
        msg += "\n".join(iter_exception_group(exc))
        logger.error(msg)
        return 1
    logger.info("Loaded TOML config without error")

    errors = 0
    if args.no_check_plugins:
        logger.info("Plugin initialisation and testing was skipped.")
    else:
        descriptions = []
        for name in conf.plugin_settings():
            label = "" if name == "default" else f" for region {name}"
            logger.info(f"Initialising Caller and Aligner{label}")
            if conf.parallel:
                # Each region runs in its own process during `readfish targets`,
                #   so check it in one too: some plugins cannot share a process
                with ProcessPoolExecutor(
                    max_workers=1, mp_context=mp.get_context("spawn")
                ) as pool:
                    result = pool.submit(_check_plugins, args.toml, channels, name)
                    caller_error, aligner_error, description = result.result()
            else:
                caller_error, aligner_error, description = _check_plugins(
                    args.toml, channels, name
                )
            if caller_error is not None:
                logger.error(f"Caller{label} could not be initialised.")
                logger.error(
                    "Possible reasons for this include a mismatch between the ont basecaller client and the version dorado you are connecting to."
                )
                logger.error(
                    "Additional information is available here: https://looselab.github.io/readfish/FAQ.html#connection-error-bad-reply-could-not-interpret-message-from-server-for-request-load-config-reply-invalid-protocol"
                )
                logger.error("See below for further details on this specific error.")
                logger.error(caller_error)
                errors += 1
            else:
                logger.info(f"Caller{label} initialised")
            if aligner_error is not None:
                logger.error(
                    f"Aligner{label} could not be initialised, see below for details"
                )
                logger.error(aligner_error)
                errors += 1
            else:
                logger.info(f"Aligner{label} initialised")
                descriptions.append(description)

        if args.no_describe:
            logger.info("Skipping descriptions of Config and Plugins.")
        elif not errors:
            logger.info(conf.describe_experiment())
            for description in descriptions:
                logger.info(description)
        else:
            logger.info("Skipping descriptions due to errors.")

    return errors
