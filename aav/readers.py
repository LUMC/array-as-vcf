"""
aav.readers
~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from math import log10
from pathlib import Path
from typing import Optional, List

from werkzeug.exceptions import NotFound

from .variation import Variant, InfoFieldNumber, InfoField, Genotype
from .lookup import RSLookup
from .utils import comma_float


GRCH37_LOOKUP = RSLookup("GRCh37")
GRCH38_LOOKUP = RSLookup("GRCh38")


class Reader(object):
    """
    Generic reader object

    Readers are iterators that produce variants
    """
    def __init__(self, path: Path, n_header_lines: int = 0):
        self.path = path
        self.handle = self.path.open()
        self.header_lines = []

        for _ in range(n_header_lines):
            self.header_lines.append(next(self.handle))

    def __next__(self) -> Variant:
        raise NotImplementedError

    def __iter__(self):
        return self

    @property
    def vcf_header(self) -> str:
        raise NotImplementedError


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

    def __init__(self, path: Path,
                 lookup_table: RSLookup,
                 qual: int = 100,
                 prefix_chr: Optional[str] = None):
        super().__init__(path, n_header_lines=1)
        self.qual = qual
        self.prefix_chr = prefix_chr
        self.lookup_table = lookup_table

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

        try:
            ref, alt = self.lookup_table[rs_id]
        except NotFound:
            ref = '.'
            alt = '.'

        if ref is None:
            ref = '.'
        if alt is None:
            alt = '.'

        return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt,
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


class CytoScanReader(Reader):
    """
    Cytoscan files are expected to conform to the following spec
    They have 12 header lines, with the following columns:

    Probe Set ID    Call Codes      Confidence      Signal A        Signal B        Forward Strand Base Calls       dbSNP RS ID     Chromosome      Chromosomal Position  # noqa

    """

    def __init__(self, path,
                 lookup_table: RSLookup,
                 prefix_chr: Optional[str] = None):
        super().__init__(path, 12)
        self.prefix_chr = prefix_chr
        self.lookup_table = lookup_table

    def __next__(self) -> Variant:
        line = next(self.handle).strip().split("\t")
        chrom = self.get_chrom(line[7])
        pos = int(line[8])
        id = line[6]

        try:
            ref, alt = self.lookup_table[id]
        except NotFound:
            ref = "."
            alt = "."

        gt = self.get_genotype(line[1])
        qual = self.get_qual(float(line[2]))

        infos = [
            InfoField("Probe_Set_ID", line[0], InfoFieldNumber.one),
            InfoField("Signal_A", line[3], InfoFieldNumber.one),
            InfoField("Signal_B", line[4], InfoFieldNumber.one)
        ]

        return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt, id=id,
                       qual=qual, info_fields=infos, genotype=gt)

    def get_genotype(self, call_code: str) -> Genotype:
        if len(set(call_code)) == 2:
            return Genotype.het
        elif call_code.upper() == "AA":
            return Genotype.hom_ref
        elif call_code.upper() == "BB":
            return Genotype.hom_alt
        else:
            return Genotype.unknown

    def get_chrom(self, chrom: str) -> str:
        if self.prefix_chr is None:
            return chrom
        return "{0}{1}".format(self.prefix_chr, chrom)

    def get_qual(self, confidence: float) -> float:
        if confidence == 0:
            return 0
        return -10 * log10(confidence)


class LumiReader(Reader):
    """
    Lumi readers have one header line.
    They required rsID lookups.

    The first two columns (rs id and chr) may be switched around
    """

    def __init__(self, path: Path,
                 lookup_table: RSLookup,
                 chr_prefix: Optional[str] = None,
                 qual=100):
        super().__init__(path, n_header_lines=1)
        self.lookup_table = lookup_table
        self.chr_prefix = chr_prefix
        self.qual = qual

    def __next__(self) -> Variant:
        line_items = next(self.handle).strip().split("\t")
        rs_id = self.get_rs_id(line_items)
        raw_chrom = self.get_raw_chrom(line_items)
        chrom = self.get_chrom(raw_chrom)
        pos = int(line_items[2])
        g_type = line_items[3]
        gt = self.get_genotype(g_type)

        infos = [
            InfoField(
                "Log_R_Ratio", comma_float(line_items[4]), InfoFieldNumber.one
            ),
            InfoField(
                "CNV_Value", int(line_items[5]), InfoFieldNumber.one
            ),
            InfoField(
                "Allele_Freq", comma_float(line_items[6]), InfoFieldNumber.one
            )
        ]

        try:
            ref, alt = self.lookup_table[rs_id]
        except NotFound:
            ref = "."
            alt = "."

        return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt,
                       qual=self.qual, id=rs_id, info_fields=infos,
                       genotype=gt)

    def get_chrom(self, chrom: str) -> str:
        if self.chr_prefix is None:
            return chrom
        return "{0}{1}".format(self.chr_prefix, chrom)

    def get_rs_id(self, line_items: List[str]) -> str:
        raise NotImplementedError

    def get_raw_chrom(self, line_items: List[str]) -> str:
        raise NotImplementedError

    def get_genotype(self, g_type: str) -> Genotype:
        if g_type == "NC":
            return Genotype.unknown
        elif len(set(g_type)) == 2:
            return Genotype.het
        elif g_type == "AA":
            return Genotype.hom_ref
        elif g_type == "BB":
            return Genotype.hom_alt
        else:
            return Genotype.unknown


class Lumi370kReader(LumiReader):

    def get_rs_id(self, line_items: List[str]) -> str:
        return line_items[0]

    def get_raw_chrom(self, line_items: List[str]) -> str:
        return line_items[1]


class Lumi317Reader(LumiReader):

    def get_rs_id(self, line_items: List[str]) -> str:
        return line_items[1]

    def get_raw_chrom(self, line_items: List[str]) -> str:
        return line_items[0]
