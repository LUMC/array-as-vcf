"""
test_readers.py
~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center

:license: MIT
"""
from datetime import date
from pathlib import Path
import pytest
from aav import __version__
from aav.readers import (AffyReader, CytoScanReader,
                         Lumi317kReader, Lumi370kReader,
                         autodetect_reader, Reader,
                         OpenArrayReader)
from aav.lookup import RSLookup
from aav.variation import Genotype


@pytest.fixture(scope="module")
def grch37_lookup():
    return RSLookup(build="GRCh37", request_tries=3)


@pytest.fixture(scope="module")
def test_lookup_table():
    lookup_path = (Path(__file__).parent / Path("data")
                   / Path("lookup_table_test.json"))
    return RSLookup.from_path(lookup_path, build="GRCh37")


_lumi_317_path = (
    Path(__file__).parent / Path("data") / Path("lumi_317_test.txt")
)
_lumi_370_path = (
    Path(__file__).parent / Path("data") / Path("lumi_370_test.txt")
)
_affy_path = (
    Path(__file__).parent / Path("data") / Path("affy_test.txt")
)
_cytoscan_path = (
    Path(__file__).parent / Path("data") / Path("cytoscan_test.txt")
)
_open_array_path = (
    Path(__file__).parent / Path("data") / Path("open_array_test.txt")
)


@pytest.fixture
def lumi_317_reader(grch37_lookup):
    return Lumi317kReader(_lumi_317_path, grch37_lookup)


@pytest.fixture
def lumi_370_reader(grch37_lookup):
    return Lumi370kReader(_lumi_370_path, grch37_lookup)


@pytest.fixture
def affy_reader(grch37_lookup):
    return AffyReader(_affy_path, grch37_lookup)


@pytest.fixture
def cytoscan_reader(grch37_lookup):
    return CytoScanReader(_cytoscan_path, grch37_lookup)


@pytest.fixture
def open_array_reader(test_lookup_table):
    return OpenArrayReader(_open_array_path, test_lookup_table,
                           "e31a0a96465a", encoding="windows-1252")


_grch37_lookup = grch37_lookup()
_test_lookup = test_lookup_table()
genotype_test_data = [
    (
        affy_reader(_grch37_lookup),
        [Genotype.unknown, Genotype.unknown, Genotype.het, Genotype.unknown,
         Genotype.unknown, Genotype.unknown, Genotype.unknown,
         Genotype.unknown, Genotype.het, Genotype.unknown,
         Genotype.unknown, Genotype.unknown]
    ),
    (
        cytoscan_reader(_grch37_lookup),
        [Genotype.hom_alt, Genotype.hom_alt, Genotype.het, Genotype.unknown,
         Genotype.unknown]
    ),
    (
        lumi_317_reader(_grch37_lookup),
        [Genotype.hom_alt, Genotype.hom_alt, Genotype.het, Genotype.unknown,
         Genotype.unknown, Genotype.unknown]
    ),
    (
        lumi_370_reader(_grch37_lookup),
        [Genotype.hom_alt, Genotype.hom_alt, Genotype.het, Genotype.unknown,
         Genotype.unknown, Genotype.unknown]
    ),
    (
        open_array_reader(_test_lookup),
        [Genotype.unknown, Genotype.het, Genotype.hom_alt, Genotype.hom_ref,
         Genotype.hom_ref, Genotype.hom_ref, Genotype.hom_ref,
         Genotype.hom_ref, Genotype.hom_alt, Genotype.hom_alt,
         Genotype.hom_alt, Genotype.hom_ref, Genotype.hom_alt,
         Genotype.hom_alt, Genotype.het, Genotype.het, Genotype.hom_alt,
         Genotype.het, Genotype.hom_alt, Genotype.unknown, Genotype.het,
         Genotype.het, Genotype.het, Genotype.hom_alt, Genotype.unknown,
         Genotype.het, Genotype.het, Genotype.hom_ref, Genotype.het,
         Genotype.het, Genotype.hom_alt, Genotype.unknown, Genotype.het,
         Genotype.het, Genotype.hom_ref, Genotype.unknown, Genotype.hom_ref,
         Genotype.hom_ref, Genotype.het, Genotype.het, Genotype.hom_alt,
         Genotype.hom_alt, Genotype.het, Genotype.hom_alt, Genotype.het,
         Genotype.hom_ref, Genotype.hom_alt, Genotype.hom_alt,
         Genotype.hom_ref, Genotype.het, Genotype.hom_ref, Genotype.het,
         Genotype.het, Genotype.unknown, Genotype.hom_ref]
    )
]

chrom_test_data = [
    (
        AffyReader(_affy_path, _grch37_lookup),
        ["1", "1", "1", "1", "X", "X", "1", "1", "1", "1", "X", "X"]
    ),
    (
        AffyReader(_affy_path, _grch37_lookup, prefix_chr="chr"),
        ["chr1", "chr1", "chr1", "chr1", "chrX", "chrX", "chr1", "chr1",
         "chr1", "chr1", "chrX", "chrX"]
    ),
    (
        CytoScanReader(_cytoscan_path, _grch37_lookup),
        ["1", "1", "1", "X", "X"]
    ),
    (
        CytoScanReader(_cytoscan_path, _grch37_lookup, prefix_chr="chr"),
        ["chr1", "chr1", "chr1", "chrX", "chrX"]
    ),
    (
        Lumi317kReader(_lumi_317_path, _grch37_lookup),
        ["1", "1", "1", "1", "X", "X"]
    ),
    (
        Lumi317kReader(_lumi_317_path, _grch37_lookup, prefix_chr="chr"),
        ["chr1", "chr1", "chr1", "chr1", "chrX", "chrX"]
    ),
    (
        Lumi370kReader(_lumi_370_path, _grch37_lookup),
        ["1", "1", "1", "1", "X", "X"]
    ),
    (
        Lumi370kReader(_lumi_370_path, _grch37_lookup, prefix_chr="chr"),
        ["chr1", "chr1", "chr1", "chr1", "chrX", "chrX"]
    ),
    (
        OpenArrayReader(_open_array_path, _test_lookup, prefix_chr="chr",
                        sample="e31a0a96465a", encoding="windows-1252"),
        ["chr3", "chr19", "chr15", "chr4", "chr13", "chr21", "chr19",
         "chr12", "chr4", "chr18", "chr14", "chrY", "chr4", "chr19",
         "chr6", "chr13", "chr18", "chr11", "chr1", "chr8", "chr15", "chr1",
         "chr20", "chr15", "chrY", "chr18", "chr8", "chr18", "chr17",
         "chr13", "chr14", "chr6", "chr17", "chr13", "chr10", "chr19",
         "chr2", "chrY", "chr6", "chr21", "chr18", "chr14", "chr6", "chr8",
         "chr18", "chr2", "chr3", "chr22", "chr17", "chr1", "chr10",
         "chr10", "chr3", "chrY", "chr18"]
    )
]

ref_test_data = [
    (
        AffyReader(_affy_path, _grch37_lookup),
        ["T", "A", "T", "T", ".", ".", "T", "A", "T", "T", ".", "."]
    ),
    (
        CytoScanReader(_cytoscan_path, _grch37_lookup),
        ["A", "T", "C", ".", "."]
    ),
    (
        Lumi317kReader(_lumi_317_path, _grch37_lookup),
        ["C", "A", "C", "A", ".", "."]
    ),
    (
        Lumi370kReader(_lumi_370_path, _grch37_lookup),
        ["C", "A", "C", "A", ".", "."]
    ),
    (
        OpenArrayReader(_open_array_path, _test_lookup,
                        sample="e31a0a96465a", encoding="windows-1252"),
        ["A", "G", "G", "T", "A", "G", "G", "C", "T", "T", "T", "C", "T",
         "G", "G", "A", "A", "G", "C", "T", "G", "C", "G", "A", "T", "T",
         "T", "G", "G", "A", "T", "G", "G", "A", "C", "G", "T", "C", "T",
         "G", "A", "T", "T", "G", "G", "A", "A", "T", "T", "T", "C", "G",
         "A", "C", "A"]
    )
]

alt_test_data = [
    (
        AffyReader(_affy_path, _grch37_lookup),
        [["C"], ["C"], ["C"], ["A", "G"], ".", ".",
         ["C"], ["C"], ["C"], ["A", "G"], ".", "."]
    ),
    (
        CytoScanReader(_cytoscan_path, _grch37_lookup),
        [["G"], ["C"], ["T"], ".", "."]
    ),
    (
        Lumi317kReader(_lumi_317_path, _grch37_lookup),
        [["T"], ["G", "T"], ["T"], ["G"], ".", "."]
    ),
    (
        Lumi370kReader(_lumi_370_path, _grch37_lookup),
        [["T"], ["G", "T"], ["T"], ["G"], ".", "."]
    )
]

autodetect_reader_data = [
    (_lumi_317_path, "Lumi317kReader", None),
    (_lumi_370_path, "Lumi370kReader", None),
    (_affy_path, "AffyReader", None),
    (_cytoscan_path, "CytoScanReader", None),
    (_open_array_path, "OpenArrayReader", 'windows-1252')
]


def test_affy_reader_amount(affy_reader):
    assert len(list(affy_reader)) == 12


def test_cytoscan_reader_amount(cytoscan_reader):
    assert len(list(cytoscan_reader)) == 5


def test_lumi317_reader_amount(lumi_317_reader):
    assert len(list(lumi_317_reader)) == 6


def test_lumi370_reader_amount(lumi_370_reader):
    assert len(list(lumi_370_reader)) == 6


def test_open_array_reader_amount(open_array_reader):
    assert len(list(open_array_reader)) == 55


@pytest.mark.parametrize("reader, genotypes", genotype_test_data)
def test_reader_genotypes(reader, genotypes):
    found_gt = [x.genotype for x in reader]
    assert found_gt == genotypes


@pytest.mark.parametrize("path, class_name, encoding", autodetect_reader_data)
def test_autodetect_readers(path, class_name, encoding):
    if encoding is None:
        found_cls = autodetect_reader(path)
    else:
        found_cls = autodetect_reader(path, encoding=encoding)
    assert found_cls.__name__ == class_name


@pytest.mark.parametrize("reader, chroms", chrom_test_data)
def test_reader_chromosomes(reader, chroms):
    found_chroms = [x.chrom for x in reader]
    assert found_chroms == chroms


@pytest.mark.parametrize("reader, refs", ref_test_data)
def test_reader_refs(reader, refs):
    found_ref = [x.ref for x in reader]
    assert found_ref == refs


@pytest.mark.parametrize("reader, alts", alt_test_data)
def test_reader_alts(reader, alts):
    for i, rec in enumerate(reader):
        assert rec.alt == alts[i]


def test_base_reader_header():
    reader = Reader(_lumi_370_path)
    header = reader.vcf_header("sample_01")
    date_str = date.today().isoformat()
    exp = (
        "##fileformat=VCFv4.2\n"
        "##fileDate={0}\n"
        "##source=aav_v{1}\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_01\n"
           ).format(date_str, __version__)
    assert header == exp


def test_affy_reader_header(affy_reader):
    header = affy_reader.vcf_header("sample_01")
    date_str = date.today().isoformat()
    exp = (
        "##fileformat=VCFv4.2\n"
        "##fileDate={0}\n"
        "##source=aav_v{1}\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##INFO=<ID=ID,Number=1,Type=String,Description="A field">\n'
        '##INFO=<ID=AffymetrixSNPsID,Number=1,Type=String,Description="A field">\n'  # noqa
        '##INFO=<ID=log2ratio_AB,Number=1,Type=Float,Description="A field">\n'
        '##INFO=<ID=N_AB,Number=1,Type=Integer,Description="A field">\n'
        '##INFO=<ID=LOH_likelihood,Number=1,Type=Float,Description="A field">\n'  # noqa
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_01\n"
           ).format(date_str, __version__)
    assert header == exp


def test_cytoscan_reader_header(cytoscan_reader):
    header = cytoscan_reader.vcf_header("sample_01")
    date_str = date.today().isoformat()
    exp = (
        "##fileformat=VCFv4.2\n"
        "##fileDate={0}\n"
        "##source=aav_v{1}\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##INFO=<ID=Probe_Set_ID,Number=1,Type=String,Description="A field">\n'
        '##INFO=<ID=Signal_A,Number=1,Type=Float,Description="A field">\n'
        '##INFO=<ID=Signal_B,Number=1,Type=Float,Description="A field">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_01\n"
           ).format(date_str, __version__)
    assert header == exp


def test_lumi_reader_header(lumi_317_reader):
    header = lumi_317_reader.vcf_header("sample_01")
    date_str = date.today().isoformat()
    exp = (
        "##fileformat=VCFv4.2\n"
        "##fileDate={0}\n"
        "##source=aav_v{1}\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##INFO=<ID=Log_R_Ratio,Number=1,Type=Float,Description="A field">\n'
        '##INFO=<ID=CNV_Value,Number=1,Type=Integer,Description="A field">\n'
        '##INFO=<ID=Allele_Freq,Number=1,Type=Float,Description="A field">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_01\n"
           ).format(date_str, __version__)
    assert header == exp


def test_open_array_reader_header(open_array_reader):
    header = open_array_reader.vcf_header("e31a0a96465a")
    date_str = date.today().isoformat()
    exp = (
        "##fileformat=VCFv4.2\n"
        "##fileDate={0}\n"
        "##source=aav_v{1}\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##INFO=<ID=Assay_Name,Number=1,Type=String,Description="A field">\n'
        '##INFO=<ID=Assay_ID,Number=1,Type=String,Description="A field">\n'
        '##INFO=<ID=Gene_Symbol,Number=.,Type=String,Description="A field">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\te31a0a96465a\n"
    ).format(date_str, __version__)
    assert header == exp
