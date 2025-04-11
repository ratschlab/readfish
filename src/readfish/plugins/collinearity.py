from functools import partial
import os
from typing import Optional, Iterable

from readfish._config import Region, Barcode
from readfish._loggers import setup_logger
from readfish.plugins.abc import AlignerABC
from readfish.plugins.utils import Result
from dataclasses import dataclass

from pycollinearity import Index, Alignment

import numpy as np
from fast_edit_distance import sub_edit_distance
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

class _Aligner(AlignerABC):
    def __init__(self, debug_log: Optional[str] = None, **kwargs):
        self.logger = setup_logger(__name__, log_file=debug_log)
        self.aligner_params = kwargs
        self.validate()
        self.aligner = Index(**kwargs)
        self._initialised = True

    def validate(self) -> None:
        if "n_threads" in self.aligner_params:
            os.environ['PARLAY_NUM_THREADS'] = self.aligner_params['n_threads']
        else:
            os.environ['PARLAY_NUM_THREADS'] = "1"
        if "input" not in self.aligner_params:
            raise RuntimeError(f"Input fasta file not provided")
        elif not (self.aligner_params['input'].endswith(".fasta") or self.aligner_params['input'].endswith(".fasta")):
            raise RuntimeError(f"Input file must have .fa or .fasta extension")

    @property
    def initialised(self) -> bool:
        return self._initialised

    def describe(self, regions: list[Region], barcodes: dict[Barcode]) -> str:
        return "Collinearity aligner"

    def map_reads(self, basecall_results: Iterable[Result]) -> Iterable[Result]:
        # queries = [query.seq for query in basecall_results]
        # alignments = self.aligner.query_batch(queries)
        # for reads, alignment in zip(basecall_results, alignments):
        #     reads.alignment_data = [alignment]
        # yield from basecall_results

        # alignment_results = []
        # for query, alignment in zip(basecall_results, alignments):
        #     alignment_results.append(query)
        #     alignment_results[-1].alignment_data = [alignment]
        # yield basecall_results

        for result in basecall_results:
            result.alignment_data = [self.aligner.query(result.seq)]
            yield result

    def disconnect(self) -> None:
        pass

Aligner = partial(_Aligner)