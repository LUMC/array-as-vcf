"""
aav.readers
~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from pathlib import Path
from typing import Optional

from .variation import Variant, InfoFieldNumber, InfoField, Genotype


class Reader(object):
    """
    Generic reader object

    Readers are iterators that produce variants
    """
    def __init__(self, path: Path, n_header_lines: int = 0):
        self.path = path
        self.handle = self.path.open()

        for _ in range(n_header_lines):
            next(self.handle)  # discard header lines

    def __next__(self) -> Variant:
        raise NotImplementedError

    def __iter__(self):
        return self


class AffyReader(Reader):
    """
    Affymetrix files are expected to conform to the follow spec:

    ID AffymetrixSNPsID rsID Chromosome Position log2ratio_AB N_AB Call_test LOH_likeihood  # noqa

    Call_test follows the following scheme:
    0: unknown
    1: hom_ref
    2: het
    3: hom_alt
    """

    def __init__(self, path: Path, qual: int = 100,
                 prefix_chr: Optional[str] = None):
        super().__init__(path, n_header_lines=1)
        self.qual = qual
        self.prefix_chr = prefix_chr

    def __next__(self) -> Variant:
        line = next(self.handle).strip().split("\t")
        chrom = self.get_chrom(line[3])
        pos = int(line[4])
        rs_id = line[2]
        gt = self.get_gt(int(line[7]))

        infos = [
            InfoField("ID", line[0], InfoFieldNumber.one),
            InfoField("AffymetrixSNPsID", line[1], InfoFieldNumber.one),
            InfoField("log2ratio_AB", line[5], InfoFieldNumber.one),
            InfoField("N_AB", line[6], InfoFieldNumber.one),
            InfoField("LOH_likelihood", line[8], InfoFieldNumber.one)
        ]

        # TODO: get ref and alt alleles!!
        return Variant(chrom=chrom, pos=pos, ref="NA", alt=["NA"],
                       qual=self.qual, id=rs_id, info_fields=infos,
                       genotype=gt)

    def get_gt(self, val: int) -> Genotype:
        if val == 0:
            return Genotype.unknown
        if val == 1:
            return Genotype.hom_ref
        if val == 2:
            return Genotype.het
        if val == 3:
            return Genotype.hom_alt
        return Genotype.unknown

    def get_chrom(self, val):
        """23 = X"""
        if val == "23":
            val = "X"

        if self.prefix_chr is not None:
            return "{0}{1}".format(self.prefix_chr, val)

        return val
