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
                         autodetect_reader, Reader)
from aav.lookup import RSLookup
from aav.variation import Genotype


@pytest.fixture(scope="module")
def grch37_lookup():
    return RSLookup(build="GRCh37", request_tries=3)


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


_grch37_lookup = grch37_lookup()
genotype_test_data = [
    (
        affy_reader(_grch37_lookup),
        [Genotype.unknown, Genotype.hom_ref, Genotype.het, Genotype.hom_ref,
         Genotype.unknown, Genotype.unknown, Genotype.unknown,
         Genotype.hom_ref, Genotype.het, Genotype.hom_ref,
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
    (_lumi_317_path, "Lumi317kReader"),
    (_lumi_370_path, "Lumi370kReader"),
    (_affy_path, "AffyReader"),
    (_cytoscan_path, "CytoScanReader")
]


def test_affy_reader_amount(affy_reader):
    assert len(list(affy_reader)) == 12


def test_cytoscan_reader_amount(cytoscan_reader):
    assert len(list(cytoscan_reader)) == 5


def test_lumi317_reader_amount(lumi_317_reader):
    assert len(list(lumi_317_reader)) == 6


def test_lumi370_reader_amount(lumi_370_reader):
    assert len(list(lumi_370_reader)) == 6


@pytest.mark.parametrize("reader, genotypes", genotype_test_data)
def test_reader_genotypes(reader, genotypes):
    found_gt = [x.genotype for x in reader]
    assert found_gt == genotypes


@pytest.mark.parametrize("path, class_name", autodetect_reader_data)
def test_autodetect_readers(path, class_name):
    found_cls = autodetect_reader(path)
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
