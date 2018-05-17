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
                         Lumi317Reader, Lumi370kReader)
from aav.lookup import RSLookup
from aav.variation import Genotype


@pytest.fixture(scope="module")
def grch37_lookup():
    return RSLookup(build="GRCh37", request_tries=3)


@pytest.fixture
def lumi_317_reader(grch37_lookup):
    p = Path(__file__).parent / Path("data") / Path("lumi_317_test.txt")
    return Lumi317Reader(p, grch37_lookup)


@pytest.fixture
def lumi_370_reader(grch37_lookup):
    p = Path(__file__).parent / Path("data") / Path("lumi_370_test.txt")
    return Lumi370kReader(p, grch37_lookup)


@pytest.fixture
def affy_reader(grch37_lookup):
    p = Path(__file__).parent / Path("data") / Path("affy_test.txt")
    return AffyReader(p, grch37_lookup)


@pytest.fixture
def cytoscan_reader(grch37_lookup):
    p = Path(__file__).parent / Path("data") / Path("cytoscan_test.txt")
    return CytoScanReader(p, grch37_lookup)


_grch37_lookup = grch37_lookup()
genotype_test_data = [
    (
        affy_reader(_grch37_lookup),
        [Genotype.unknown, Genotype.hom_ref, Genotype.het, Genotype.hom_alt]
    ),
    (
        cytoscan_reader(_grch37_lookup),
        [Genotype.hom_alt, Genotype.hom_ref, Genotype.het]
    ),
    (
        lumi_317_reader(_grch37_lookup),
        [Genotype.hom_alt, Genotype.hom_ref, Genotype.het, Genotype.unknown]
    ),
    (
        lumi_370_reader(_grch37_lookup),
        [Genotype.hom_alt, Genotype.hom_ref, Genotype.het, Genotype.unknown]
    )
]


def test_affy_reader_amount(affy_reader):
    assert len(list(affy_reader)) == 4


def test_cytoscan_reader_amount(cytoscan_reader):
    assert len(list(cytoscan_reader)) == 3


def test_lumi317_reader_amount(lumi_317_reader):
    assert len(list(lumi_317_reader)) == 4


def test_lumi370_reader_amount(lumi_370_reader):
    assert len(list(lumi_370_reader)) == 4


@pytest.mark.parametrize("reader, genotypes", genotype_test_data)
def test_reader_genotypes(reader, genotypes):
    found_gt = [x.genotype for x in reader]
    assert found_gt == genotypes
