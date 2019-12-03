"""
aav.readers
~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from functools import reduce
from math import log10
from typing import Optional, List, Type, Tuple, Set
import logging

from .variation import (Variant, InfoFieldNumber, InfoField, Genotype,
                        InfoHeaderLine, GT_FORMAT, VCF_v_4_2,
                        program_header, date_header, chrom_header,
                        InfoFieldType)
from .lookup import RSLookup
from .utils import comma_float, empty_string


GRCH37_LOOKUP = RSLookup("GRCh37")
GRCH38_LOOKUP = RSLookup("GRCh38")


logger = logging.getLogger('ArrayReader')


class Reader(object):
    """
    Generic reader object

    Readers are iterators that produce variants
    """
    def __init__(self, path: str, n_header_lines: int = 0,
                 encoding: Optional[str] = None):
        self.path = path
        self.handle = open(path, mode="r", encoding=encoding)
        self.header_lines = []

        self.header_fields = [
            VCF_v_4_2, date_header(), program_header(), GT_FORMAT
        ]

        for _ in range(n_header_lines):
            self.header_lines.append(next(self.handle))

    def __next__(self) -> Variant:
        raise NotImplementedError

    def __iter__(self):
        return self

    def vcf_header(self, sample_name: str) -> str:
        s = reduce(lambda x, y: x + str(y) + "\n", self.header_fields, "")
        return s + chrom_header(sample_name) + '\n'


class OpenArrayReader(Reader):
    def __init__(self, path: str, lookup_table: RSLookup, sample: str,
                 qual: int = 100, prefix_chr: Optional[str] = None,
                 encoding: Optional[str] = None,
                 exclude_assays: Optional[Set[str]] = None):
        super().__init__(path, n_header_lines=18, encoding=encoding)
        self.qual = qual
        self.sample = sample
        self.lookup_table = lookup_table
        self.prefix_chr = prefix_chr
        self.linecount = 18  # n_header_lines
        self.unknown_call = {'INV', 'NOAMP', 'UND', '-/-'}
        if exclude_assays is not None:
            self.exclude_assays = exclude_assays
        else:
            self.exclude_assays = set()

        self.header_fields += [
            InfoHeaderLine("Assay_Name", InfoFieldNumber.one,
                           InfoFieldType.STRING),
            InfoHeaderLine("Assay_ID", InfoFieldNumber.one,
                           InfoFieldType.STRING),
            InfoHeaderLine("Gene_Symbol", InfoFieldNumber.unknown,
                           InfoFieldType.STRING)
        ]

        self._header_splitted = self.header_lines[-1].strip().split("\t")

    @property
    def chromsome_col_idx(self) -> int:
        return self._header_splitted.index("Chromosome #")

    @property
    def position_col_idx(self) -> int:
        return self._header_splitted.index("Position")

    @property
    def sample_col_idx(self) -> int:
        return self._header_splitted.index("Sample ID")

    @property
    def rsid_col_idx(self) -> int:
        return self._header_splitted.index("NCBI SNP Reference")

    @property
    def assay_name_col_idx(self) -> int:
        return self._header_splitted.index("Assay Name")

    @property
    def assay_id_col_idx(self) -> int:
        return self._header_splitted.index("Assay ID")

    @property
    def gene_symbol_col_idx(self) -> int:
        return self._header_splitted.index("Gene Symbol")

    @property
    def call_col_idx(self) -> int:
        return self._header_splitted.index("Call")

    def __next__(self):
        for raw_line in self.handle:
            self.linecount += 1
            if empty_string(raw_line):
                raise StopIteration  # end of initial list
            line = raw_line.strip('\n').split("\t")
            if len(line) < 8:  # may occur if assay design is dumped in file
                logger.debug(f"Skipping line {self.linecount}, to few columns")
                continue
            assay_id = line[self.assay_id_col_idx]
            if assay_id in self.exclude_assays:
                logger.debug("Skipping excluded assay {assay_id}")
                continue
            line_sample = line[self.sample_col_idx]
            if line_sample != self.sample:
                logger.debug(f"Skipping line {self.linecount}, wrong sample "
                             f"({line_sample} is not {self.sample})")
                continue
            rs_id = line[self.rsid_col_idx].strip()  # may have spaces :cry:
            try:
                raw_chrom = line[self.chromsome_col_idx]
            except IndexError:  # sometimes the entire row is truncated
                logger.debug((f"Skipping line {self.linecount}, entire row "
                              "truncated"))
                continue
            pos = line[self.position_col_idx]

            # Skip if fields we need are missing
            if empty_string(raw_chrom):
                logger.debug((f"Skipping line {self.linecount}, missing "
                              "chromosome"))
                continue
            if empty_string(pos):
                logger.debug((f"Skipping line {self.linecount}, missing "
                              "position"))
                continue
            if empty_string(rs_id):
                logger.debug((f"Skipping line {self.linecount}, missing "
                              "rs_id"))
                continue

            # Also skip if the rs_id is not in the lookup_table
            try:
                q_res = self.lookup_table[rs_id]
            except KeyError:
                logger.debug(f"Skipping {rs_id}, transcript not found")
                continue
            else:
                call = line[self.call_col_idx]
                ref = q_res.ref
                genotype, alt = self.get_genotype_and_alt(call, ref, q_res.alt)

            assay_name = line[self.assay_name_col_idx]
            raw_gene_symbol = line[self.gene_symbol_col_idx]

            infos = [
                InfoField("Assay_Name", assay_name, InfoFieldNumber.one),
                InfoField("Assay_ID", assay_id, InfoFieldNumber.one)
            ]

            if not empty_string(raw_gene_symbol):
                infos.append(InfoField("Gene_Symbol",
                                       raw_gene_symbol.split(";"),
                                       InfoFieldNumber.unknown))

            chrom = self.get_chrom(raw_chrom)
            return Variant(chrom=chrom, pos=int(pos), id=rs_id, ref=ref,
                           alt=alt, info_fields=infos, qual=self.qual,
                           genotype=genotype)
        else:
            raise StopIteration

    def get_chrom(self, chrom: str) -> str:
        if self.prefix_chr is None:
            return chrom
        return "{0}{1}".format(self.prefix_chr, chrom)

    def get_genotype_and_alt(self, call: str, ref: str,
                             fallback_alt: str) -> Tuple[Genotype, str]:
        # These are calls that indicate that no genotype could be determined
        if call in self.unknown_call:
            msg = f"Recognised {call}, which has no genotype information"
            logger.debug(msg)
            return Genotype.unknown, "."

        alleles = set(call.split("/"))
        # These are alleles that are not handled by this tool, so we print an
        # error when we encounter them
        for allele in alleles:
            if not set(allele).issubset({"A", "T", "C", "G", "N"}):
                logger.error(f"Skipping Unknown call {call}")
                return Genotype.unknown, "."

        if len(alleles) > 1:
            return Genotype.het, (alleles - {ref}).pop()

        homozygous_allele = alleles.pop()
        if homozygous_allele == ref:
            return Genotype.hom_ref, fallback_alt
        else:
            return Genotype.hom_alt, homozygous_allele


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

    def __init__(self, path: str,
                 lookup_table: RSLookup,
                 qual: int = 100,
                 prefix_chr: Optional[str] = None,
                 encoding: Optional[str] = None):
        super().__init__(path, n_header_lines=1, encoding=encoding)
        self.qual = qual
        self.prefix_chr = prefix_chr
        self.lookup_table = lookup_table

        self.header_fields += [
            InfoHeaderLine("ID", InfoFieldNumber.one, InfoFieldType.STRING),
            InfoHeaderLine("AffymetrixSNPsID", InfoFieldNumber.one,
                           InfoFieldType.STRING),
            InfoHeaderLine("log2ratio_AB", InfoFieldNumber.one,
                           InfoFieldType.FLOAT),
            InfoHeaderLine("N_AB", InfoFieldNumber.one, InfoFieldType.INT),
            InfoHeaderLine("LOH_likelihood", InfoFieldNumber.one,
                           InfoFieldType.FLOAT)
        ]

    def __next__(self) -> Variant:
        for raw_line in self.handle:
            line = raw_line.strip('\n').split("\t")
            chrom = self.get_chrom(line[3])
            pos = int(line[4])
            rs_id = line[2]
            try:
                q_res = self.lookup_table[rs_id]
            except KeyError:
                logger.info(f"Skipping {rs_id}, transcript not found")
                continue
            else:
                if q_res is None or q_res.ref_is_minor is None:
                    logger.info(f"Skipping {rs_id}, incomplete data: {q_res}")
                    continue
                else:
                    ref = q_res.ref
                    alt = q_res.alt
                    ref_is_minor = q_res.ref_is_minor
                    gt = self.get_gt(int(line[7]), ref_is_minor)

            infos = [
                InfoField("ID", line[0], InfoFieldNumber.one),
                InfoField("AffymetrixSNPsID", line[1], InfoFieldNumber.one),
                InfoField("log2ratio_AB", line[5], InfoFieldNumber.one),
                InfoField("N_AB", line[6], InfoFieldNumber.one),
                InfoField("LOH_likelihood", line[8], InfoFieldNumber.one)
            ]

            return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt,
                           qual=self.qual, id=rs_id, info_fields=infos,
                           genotype=gt)
        else:
            raise StopIteration

    def get_gt(self, val: int, ref_is_minor: bool) -> Genotype:
        if val == 0:
            return Genotype.unknown
        elif val == 1:
            return Genotype.unknown
        elif val == 2:
            return Genotype.het
        elif val == 3:
            return Genotype.unknown
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
                 prefix_chr: Optional[str] = None,
                 encoding: Optional[str] = None):
        super().__init__(path, 12, encoding=encoding)
        self.prefix_chr = prefix_chr
        self.lookup_table = lookup_table

        self.header_fields += [
            InfoHeaderLine("Probe_Set_ID", InfoFieldNumber.one,
                           InfoFieldType.STRING),
            InfoHeaderLine("Signal_A", InfoFieldNumber.one,
                           InfoFieldType.FLOAT),
            InfoHeaderLine("Signal_B", InfoFieldNumber.one,
                           InfoFieldType.FLOAT)
        ]

    def __next__(self) -> Variant:
        for raw_line in self.handle:
            line = raw_line.strip('\n').split("\t")
            chrom = self.get_chrom(line[7])
            pos = int(line[8])
            rs_id = line[6]

            try:
                q_res = self.lookup_table[rs_id]
            except KeyError:
                logger.debug(f"Skipping {rs_id}, transcript not found")
                continue
            else:
                if q_res is None or q_res.ref_is_minor is None:
                    logger.debug(f"Skipping {rs_id}, incomplete data: {q_res}")
                    continue
                else:
                    ref = q_res.ref
                    alt = q_res.alt
                    gt = self.get_genotype(ref, alt, line[5])

            qual = self.get_qual(float(line[2]))

            infos = [
                InfoField("Probe_Set_ID", line[0], InfoFieldNumber.one),
                InfoField("Signal_A", line[3], InfoFieldNumber.one),
                InfoField("Signal_B", line[4], InfoFieldNumber.one)
            ]

            return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt, id=rs_id,
                           qual=qual, info_fields=infos, genotype=gt)
        else:
            raise StopIteration

    def get_genotype(self, ref: str, alt: List[str], calls: str) -> Genotype:
        if calls is None or calls == "":
            return Genotype.unknown

        alleles = list(calls)
        if len(set(alleles)) > 1:
            return Genotype.het
        elif alleles[0] == ref:
            return Genotype.hom_ref
        elif any([alleles[0] == x for x in alt]):
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

    def __init__(self, path: str,
                 lookup_table: RSLookup,
                 prefix_chr: Optional[str] = None,
                 qual=100,
                 encoding: Optional[str] = None):
        super().__init__(path, n_header_lines=1, encoding=encoding)
        self.lookup_table = lookup_table
        self.chr_prefix = prefix_chr
        self.qual = qual

        self.header_fields += [
            InfoHeaderLine("Log_R_Ratio", InfoFieldNumber.one,
                           InfoFieldType.FLOAT),
            InfoHeaderLine("CNV_Value", InfoFieldNumber.one,
                           InfoFieldType.INT),
            InfoHeaderLine("Allele_Freq", InfoFieldNumber.one,
                           InfoFieldType.FLOAT)
        ]

    def __next__(self) -> Variant:
        for raw_line in self.handle:
            line = raw_line.strip('\n').split("\t")
            rs_id = self.get_rs_id(line)
            raw_chrom = self.get_raw_chrom(line)
            chrom = self.get_chrom(raw_chrom)
            pos = int(line[2])
            g_type = line[3]

            try:
                q_res = self.lookup_table[rs_id]
            except KeyError:
                logger.info(f"Skipping {rs_id}, transcript not found")
                continue
            else:
                if q_res is None or q_res.ref_is_minor is None:
                    logger.info(f"Skipping {rs_id}, incomplete data: {q_res}")
                    continue
                else:
                    ref = q_res.ref
                    alt = q_res.alt
                    ref_is_minor = q_res.ref_is_minor
                    gt = self.get_genotype(g_type, ref_is_minor)

            infos = [
                InfoField(
                    "Log_R_Ratio", comma_float(line[4]), InfoFieldNumber.one
                ),
                InfoField(
                    "CNV_Value", int(line[5]), InfoFieldNumber.one
                ),
                InfoField(
                    "Allele_Freq", comma_float(line[6]), InfoFieldNumber.one
                )
            ]

            return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt,
                           qual=self.qual, id=rs_id, info_fields=infos,
                           genotype=gt)
        else:
            raise StopIteration

    def get_chrom(self, chrom: str) -> str:
        if self.chr_prefix is None:
            return chrom
        return "{0}{1}".format(self.chr_prefix, chrom)

    def get_rs_id(self, line: List[str]) -> str:
        raise NotImplementedError

    def get_raw_chrom(self, line: List[str]) -> str:
        raise NotImplementedError

    def get_genotype(self, g_type: str, ref_is_minor: bool) -> Genotype:
        if g_type == "NC":
            return Genotype.unknown
        elif len(set(g_type)) == 2:
            return Genotype.het
        elif g_type == "AA" and not ref_is_minor:
            return Genotype.hom_ref
        elif g_type == "AA" and ref_is_minor:
            return Genotype.hom_alt
        elif g_type == "BB" and not ref_is_minor:
            return Genotype.hom_alt
        elif g_type == "BB" and ref_is_minor:
            return Genotype.hom_ref
        else:
            return Genotype.unknown


class Lumi370kReader(LumiReader):

    def get_rs_id(self, line: List[str]) -> str:
        return line[1]

    def get_raw_chrom(self, line: List[str]) -> str:
        return line[0]


class Lumi317kReader(LumiReader):

    def get_rs_id(self, line: List[str]) -> str:
        return line[0]

    def get_raw_chrom(self, line: List[str]) -> str:
        return line[1]


def autodetect_reader(path: str,
                      encoding: Optional[str] = None) -> Type[Reader]:
    """
    Detect type of reader for a certain array path
    :param path: instance of string pointing to path
    :param encoding: optional encoding of file
    :return: Reader class (NOT instance)
    :raises: NotImplementedError for unknown types.
    """
    pot_affy = None
    with open(path, encoding=encoding, mode="r") as handle:
        for i, line in enumerate(handle):
            if i == 0 and "Affymetrix" in line:
                pot_affy = line
            elif i == 0 and line.startswith("Name"):
                return Lumi317kReader
            elif i == 0 and line.startswith("Chr"):
                return Lumi370kReader
            elif i == 11 and line.startswith("Probe"):
                return CytoScanReader
            elif i > 11 and pot_affy is not None:
                return AffyReader
            elif i == 17 and line.startswith("Assay Name"):
                return OpenArrayReader
            elif i >= 18:
                raise NotImplementedError("Could not detect type of array")

    raise NotImplementedError
