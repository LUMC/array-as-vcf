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


@pytest.fixture
def lumi_317_reader():
    p = Path(__file__).parent / Path("data") / Path("lumi_317_test.txt")
    return Lumi317Reader(p, RSLookup(build="GRCh37", request_tries=3))


@pytest.fixture
def lumi_370_reader():
    p = Path(__file__).parent / Path("data") / Path("lumi_370_test.txt")
    return Lumi370kReader(p, RSLookup(build="GRCh37", request_tries=3))


@pytest.fixture
def affy_reader():
    p = Path(__file__).parent / Path("data") / Path("affy_test.txt")
    return AffyReader(p, RSLookup(build="GRCh37", request_tries=3))


@pytest.fixture
def cytoscan_reader():
    p = Path(__file__).parent / Path("data") / Path("cytoscan_test.txt")
    return CytoScanReader(p, RSLookup(build="GRCh37", request_tries=3))


def test_affy_reader_amount(affy_reader):
    assert len(list(affy_reader)) == 4


def test_cytoscan_reader_amount(cytoscan_reader):
    assert len(list(cytoscan_reader)) == 3


def test_lumi317_reader_amount(lumi_317_reader):
    assert len(list(lumi_317_reader)) == 4


def test_lumi370_reader_amount(lumi_370_reader):
    assert len(list(lumi_370_reader)) == 4
