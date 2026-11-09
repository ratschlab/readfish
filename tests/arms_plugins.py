"""Test plugins for arms_test.py, importable by the spawned worker processes.

The caller "base-calls" by decoding the raw data bytes as the sequence. The
aligner reports a read as mapped to ``contig`` if its sequence contains
``motif``, after sleeping ``delay`` seconds per batch, so that tests can make
one region slow.
"""

import time
from collections import namedtuple

from readfish.plugins.abc import AlignerABC, CallerABC
from readfish.plugins.utils import Result

Alignment = namedtuple("Alignment", "ctg r_st r_en strand")


class Caller(CallerABC):
    def __init__(self, *args, **kwargs):
        pass

    def validate(self):
        return None

    def basecall(self, reads, signal_dtype, daq_values=None):
        for channel, read in reads:
            yield Result(
                channel=channel,
                read_id=read.id,
                seq=read.raw_data.decode(),
            )

    def describe(self):
        return "test caller"

    def disconnect(self):
        return None


class Aligner(AlignerABC):
    def __init__(self, motif="ACGT", contig="target", delay=0.0, fail=False, **kwargs):
        if fail:
            raise RuntimeError("aligner failed on purpose")
        self.motif = motif
        self.contig = contig
        self.delay = float(delay)

    def validate(self):
        return None

    @property
    def initialised(self):
        return True

    def map_reads(self, basecall_results):
        time.sleep(self.delay)
        for result in basecall_results:
            if self.motif in result.seq:
                result.alignment_data = [Alignment(self.contig, 0, len(result.seq), 1)]
            else:
                result.alignment_data = []
            yield result

    def describe(self, regions, barcodes):
        return f"test aligner {self.motif} delay {self.delay}"

    def disconnect(self):
        return None
