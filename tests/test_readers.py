"""
test_readers.py
~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center

:license: MIT
"""
from pathlib import Path
import pytest
from aav.readers import (AffyReader, CytoScanReader,
                         Lumi317kReader, Lumi370kReader,
                         autodetect_reader)
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
        [Genotype.unknown, Genotype.hom_ref, Genotype.het, Genotype.hom_alt,
         Genotype.unknown, Genotype.unknown]
    ),
    (
        cytoscan_reader(_grch37_lookup),
        [Genotype.hom_alt, Genotype.hom_ref, Genotype.het, Genotype.unknown,
         Genotype.unknown]
    ),
    (
        lumi_317_reader(_grch37_lookup),
        [Genotype.hom_alt, Genotype.hom_ref, Genotype.het, Genotype.unknown,
         Genotype.unknown, Genotype.unknown]
    ),
    (
        lumi_370_reader(_grch37_lookup),
        [Genotype.hom_alt, Genotype.hom_ref, Genotype.het, Genotype.unknown,
         Genotype.unknown, Genotype.unknown]
    )
]

chrom_test_data = [
    (
        AffyReader(_affy_path, _grch37_lookup),
        ["1", "1", "1", "1", "X", "X"]
    ),
    (
        AffyReader(_affy_path, _grch37_lookup, prefix_chr="chr"),
        ["chr1", "chr1", "chr1", "chr1", "chrX", "chrX"]
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

autodetect_reader_data = [
    (_lumi_317_path, "Lumi317kReader"),
    (_lumi_370_path, "Lumi370kReader"),
    (_affy_path, "AffyReader"),
    (_cytoscan_path, "CytoScanReader")
]


def test_affy_reader_amount(affy_reader):
    assert len(list(affy_reader)) == 6


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
