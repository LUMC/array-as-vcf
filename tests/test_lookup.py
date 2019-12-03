"""
test_lookup.py
~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""

import os
import json
from datetime import datetime

import pytest
from requests.exceptions import Timeout

from array_as_vcf.lookup import query_ensembl, RSLookup, QueryResult


@pytest.fixture
def grch37_lookup():
    return RSLookup("GRCh37")


@pytest.fixture
def grch38_lookup():
    return RSLookup("GRCh38")


@pytest.fixture
def lookup_table():
    return os.path.join("tests", "data", "lookup_table_test.json")


@pytest.mark.xfail
def test_query_ensembl_known_rsid_grch38():
    assert query_ensembl("rs56116432", "GRCh38") == QueryResult(
        "C", ["T"], False
    )


@pytest.mark.xfail
def test_query_ensembl_known_rsid_grch37():
    assert query_ensembl("rs56116432", "GRCh37") == QueryResult(
        "C", ["T"], False
    )


@pytest.mark.xfail
def test_query_ensembl_unknown_minor_grch37():
    assert query_ensembl("rs3913290", "GRCh37") == QueryResult(
        "C", ["T"], None
    )


@pytest.mark.xfail
def test_query_multi_alt_grch38():
    assert query_ensembl("rs60", "GRCh38") == QueryResult(
        "A", ["G", "T"], True
    )


@pytest.mark.xfail
def test_query_multi_alt_grch37():
    assert query_ensembl("rs60", "GRCh37") == QueryResult(
        "A", ["G", "T"], True
    )


@pytest.mark.xfail
def test_query_ensembl_unknown_build():
    with pytest.raises(NotImplementedError):
        query_ensembl("rs56116432", "unknown")


@pytest.mark.xfail
def test_unknown_rsid():
    with pytest.raises(RuntimeError) as exc:
        query_ensembl("rs5611644432", "GRCh37")
        assert "rsID not found for human" in str(exc)


def test_timeout():
    with pytest.raises(Timeout):
        query_ensembl("rs60", "GRCh37", 0.00001)


@pytest.mark.xfail
def test_lookup_succeed_grch37(grch37_lookup):
    assert grch37_lookup['rs56116432'] == QueryResult(
        "C", ["T"], False
    )


@pytest.mark.xfail
def test_lookup_fail_grch37(grch37_lookup):
    with pytest.raises(KeyError, match='rs5611644432'):
        grch37_lookup['rs5611644432']


@pytest.mark.xfail
def test_lookup_succeed_grhc38(grch38_lookup):
    assert grch38_lookup['rs56116432'] == QueryResult(
        "C", ["T"], False
    )


@pytest.mark.xfail
def test_lookup_fail_grch38(grch38_lookup):
    with pytest.raises(KeyError, match='rs5611644432'):
        grch38_lookup['rs5611644432']


@pytest.mark.xfail
def test_lookup_speed_grch37(grch37_lookup):
    now = datetime.utcnow()
    _ = grch37_lookup['rs56']
    after = datetime.utcnow()
    init_delta = after - now

    _ = grch37_lookup['rs56']  # noqa
    after2 = datetime.utcnow()
    second_delta = after2 - after

    assert second_delta < init_delta


@pytest.mark.xfail
def test_lookup_speed_grch38(grch38_lookup):
    now = datetime.utcnow()
    _ = grch38_lookup['rs56']
    after = datetime.utcnow()
    init_delta = after - now

    _ = grch38_lookup['rs56']  # noqa
    after2 = datetime.utcnow()
    second_delta = after2 - after

    assert second_delta < init_delta


def test_lookup_timeout():
    look = RSLookup(build="GRCh37", request_timeout=0.0001)
    with pytest.raises(KeyError, match='rs56'):
        look['rs56']


@pytest.mark.xfail
def test_lookup_timeout_tries():
    look = RSLookup(build="GRCh37", request_timeout=0.0001, request_tries=5)
    with pytest.raises(KeyError, match='r56'):
        look['r56']


def test_lookup_table_load(lookup_table):
    lookup = RSLookup.from_path(lookup_table, "GRCh37")
    assert lookup._RSLookup__rsids == {
        "rs776746": QueryResult("C", ["T"], False),
        "rs497692": QueryResult("T", ["C"], False),
        "rs12248560": QueryResult("C", ["A", "T"], False),
        "rs768983": QueryResult("C", ["T"], False),
        "rs2032598": QueryResult("T", ["C"], False),
        "rs1884213": QueryResult("T", ["C"], True),
        "rs1633021": QueryResult("T", ["C"], False),
        "rs4735258": QueryResult("T", ["C"], False),
        "rs10373": QueryResult("A", ["G"], True),
        "rs2819561": QueryResult("A", ["G"], True),
        "rs7514030": QueryResult("T", ["C"], False),
        "rs7300444": QueryResult("C", ["T"], False),
        "rs622272": QueryResult("T", ["G"], True),
        "rs629562": QueryResult("A", ["G"], True),
        "rs2249028": QueryResult("G", ["A"], False),
        "rs2711823": QueryResult("G", ["A"], False),
        "rs1381532": QueryResult("A", ["G"], True),
        "rs4148973": QueryResult("T", ["G"], True),
        "rs2942": QueryResult("G", ["A"], True),
        "rs1019433": QueryResult("G", ["A"], False),
        "rs1037256": QueryResult("G", ["A", "C", "T"], True),
        "rs6124288": QueryResult("T", ["C"], False),
        "rs10203363": QueryResult("C", ["T"], False),
        "rs2229546": QueryResult("C", ["A", "G", 'T'], True),
        "rs4577050": QueryResult("G", ["A"], True),
        "rs9962023": QueryResult("T", ["C"], True),
        "rs9532292": QueryResult("A", ["G"], False),
        "rs10421632": QueryResult("G", ["A"], False),
        "rs4688963": QueryResult("T", ["C"], False),
        "rs2844682": QueryResult("G", ["A"], False),
        "rs1410592": QueryResult("G", ["A", "T"], True),
        "rs2228611": QueryResult("T", ["C"], True),
        "rs5993935": QueryResult("T", ["C"], False),
        "rs309557": QueryResult("T", ["C"], False),
        "rs2070203": QueryResult("G", ["A"], False),
        "rs8024825": QueryResult("A", ["G"], False),
        "rs9528543": QueryResult("A", ['G'], False),
        "rs9361875": QueryResult("C", ["T"], False),
        "rs9837496": QueryResult("C", ["A", "T"], False),
        "rs10965655": QueryResult("T", ["C"], False),
        "rs357004": QueryResult("A", ["G"], True),
        "rs998132": QueryResult("G", ["A"], False),
        "rs4675": QueryResult("T", ["C"], True),
        "rs17548783": QueryResult("T", ["C"], False),
        "rs17686195": QueryResult("T", ["C"], False),
        "rs17476242": QueryResult("G", ["A"], False),
        "rs2297995": QueryResult("G", ["A"], True),
        "rs349047": QueryResult("G", ["A"], False),
        "rs3912984": QueryResult("T", ["C"], True),
        "rs4617548": QueryResult("A", ["G"], True),
        "rs3913290": QueryResult("C", ["T"], None),
        "rs1147504": QueryResult("G", ["A"], True),
        "rs10883099": QueryResult("G", ["A"], True),
        "rs2395029": QueryResult("T", ["G"], False),
        "rs2980300": QueryResult("T", ["C"], True),
        "rs10907175": QueryResult("A", ["C"], False),
        "rs2887286": QueryResult("T", ["C"], True),
        "rs307378": QueryResult("T", ["A", "G"], True),
        "rs6696609": QueryResult("C", ["G", "T"], False),
        "rs3737728": QueryResult("A", ["C", "G", "T"], True),
        "rs12939215": QueryResult("A", ["C"], None)
    }


def test_lookup_table_dump(lookup_table):
    with open(lookup_table, "r") as handle:
        original = json.load(handle)

    rs_lookup = RSLookup.from_path(lookup_table, "GRCh37")
    generated = json.loads(rs_lookup.dumps())
    assert generated == original


@pytest.mark.xfail
def test_lookup_online():
    look = RSLookup(build="GRCh37", ensembl_lookup=True)
    look['rs3934834']


def test_lookup_offline():
    look = RSLookup(build="GRCh37", ensembl_lookup=False)
    with pytest.raises(KeyError):
        look['rs3934834']
