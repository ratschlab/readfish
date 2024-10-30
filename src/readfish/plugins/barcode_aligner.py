from functools import partial
from typing import Optional, Iterable

from readfish._config import Region, Barcode
from readfish._loggers import setup_logger
from readfish.plugins.abc import AlignerABC
from readfish.plugins.utils import Result
from dataclasses import dataclass

import numpy as np
from fast_edit_distance import sub_edit_distance
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor


@dataclass
class DummyAlignment:
    ctg: str
    r_st: int
    r_en: int
    strand: int  # +1 or -1


def get_most_probable_positions(seq: str):
    barcodes = ["AAAAAAAAAAGAUUCAGCAG", "AUACGGUCUGGAUCGUUGAC", "UAGCACUGAGGAAUCAGUCC"]
    n = len(seq)
    k = 20
    positions = np.full(len(barcodes), -1, dtype=np.int32)
    start = 0
    for i in range(len(barcodes)):
        distance, end = sub_edit_distance(barcodes[i], seq, max_ed=10)
        positions[i] = start + end - k
        start = end
        remaining = n - end
        if remaining < 15: break
    return positions


def has_barcode(query: Result) -> Result:
    n = len(query.seq)
    if 'AAAAAAAAAA' not in query.seq:  # not yet seen poly(A)
        query.alignment_data = []  # (proceed)
    else:
        positions = get_most_probable_positions(query.seq)
        if positions[0] >= 0:  # found poly(A)+Linker1
            if n - positions[0] < 100:  # not yet enough sequence to find other linkers
                query.alignment_data = []  # (proceed)
            else:
                if positions[1] > positions[0]:  # found Linker 2
                    if 25 < positions[1] - positions[0] < 45:  # Linker 2 is at a good distance from poly(A)+Linker1
                        if positions[2] > positions[1]:  # found Linker 3, check dist. from Linker 2
                            if 25 < positions[2] - positions[1] < 45:  # good distance
                                query.alignment_data = [DummyAlignment("yes", 0, 1, 1)]  # aligned (stop_receiving)
                            else:
                                query.alignment_data = [
                                    DummyAlignment("no", 0, 1, 1)]  # Linker 3 found but not in a good place (eject)
                        else:
                            query.alignment_data = [DummyAlignment("no", 0, 1, 1)]  # Linker 3 not found (eject)
                    else:
                        query.alignment_data = [
                            DummyAlignment("no", 0, 1, 1)]  # Linker 2 found but not in a good place (eject)
                else:
                    query.alignment_data = [DummyAlignment("no", 0, 1, 1)]  # Linker 2 not found. (eject)
        else:
            query.alignment_data = []  # haven't found poly(A)+Linker1 yet todo - low long do we keep checking?
    return query


class BarcodeAligner(AlignerABC):
    def __init__(self, debug_log: Optional[str] = None, **kwargs) -> None:
        self.logger = setup_logger(__name__, log_file=debug_log)
        self.aligner_params = kwargs
        self.validate()
        self.executor = ThreadPoolExecutor(max_workers=kwargs['n_threads'])

    def validate(self) -> None:
        if 'n_threads' not in self.aligner_params:
            raise RuntimeError("Need the number of threads for alignment.")

    @property
    def initialised(self) -> bool:
        return True

    def describe(self, regions: list[Region], barcodes: dict[Barcode]) -> str:
        return "This is where I will return a proper description someday"

    def map_reads(self, basecall_results: Iterable[Result]) -> Iterable[Result]:
        return list(self.executor.map(has_barcode, basecall_results))

    def disconnect(self) -> None:
        return


Aligner = partial(BarcodeAligner)

