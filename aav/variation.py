"""
aav.variation
~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from typing import List, Optional


class Variant(object):
    """
    Representations of a variant

    TODO: add info and genotypes
    """
    def __init__(self, chrom: str, pos: int, ref: str, alt: List[str],
                 qual: float, filters: List[str] = list(),
                 id: Optional[str] = None):
        self.chrom = chrom
        self.pos = pos
        if id is not None:
            self.id = id
        else:
            self.id = "."
        self.qual = qual
        self.ref = ref
        self.alt = alt
        self.filters = filters

    @property
    def vcf_line(self):
        fmt = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(
            self.chrom, self.pos, self.id,
            self.ref, ",".join(self.alt),
            self.qual, ",".join(self.filters)
        )
        return fmt
