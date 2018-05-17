"""
test_lookup.py
~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""

from aav.lookup import query_ensembl, RSLookup
from datetime import datetime

import pytest
from werkzeug.exceptions import NotFound
from requests.exceptions import Timeout


@pytest.fixture
def grch37_lookup():
    return RSLookup("GRCh37")


@pytest.fixture
def grch38_lookup():
    return RSLookup("GRCh38")


def test_query_ensembl_known_rsid_grch38():
    assert query_ensembl("rs56116432", "GRCh38") == ("C", ["T"])


def test_query_ensembl_known_rsid_grch37():
    assert query_ensembl("rs56116432", "GRCh37") == ("C", ["T"])


def test_query_multi_alt_grch38():
    assert query_ensembl("rs60", "GRCh38") == ("A", ["G", "T"])


def test_query_multi_alt_grch37():
    assert query_ensembl("rs60", "GRCh37") == ("A", ["G", "T"])


def test_query_ensembl_unknown_build():
    with pytest.raises(NotImplementedError):
        query_ensembl("rs56116432", "unknown")


def test_unknown_rsid():
    with pytest.raises(NotFound) as exc:
        query_ensembl("rs5611644432", "GRCh37")
        assert "rsID not found for human" in str(exc)


def test_timeout():
    with pytest.raises(Timeout):
        query_ensembl("rs60", "GRCh37", 0.00001)


def test_lookup_succeed_grch37(grch37_lookup):
    assert grch37_lookup['rs56116432'] == ("C", ["T"])


def test_lookup_fail_grch37(grch37_lookup):
    with pytest.raises(NotFound):
        grch37_lookup['rs5611644432']


def test_lookup_succeed_grhc38(grch38_lookup):
    assert grch38_lookup['rs56116432'] == ("C", ["T"])


def test_lookup_fail_grch38(grch38_lookup):
    with pytest.raises(NotFound):
        grch38_lookup['rs5611644432']


def test_lookup_speed_grch37(grch37_lookup):
    now = datetime.utcnow()
    _ = grch37_lookup['rs56']
    after = datetime.utcnow()
    init_delta = after - now

    _ = grch37_lookup['rs56']  # noqa
    after2 = datetime.utcnow()
    second_delta = after2 - after

    assert second_delta < init_delta


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
    with pytest.raises(ValueError):
        look['rs56']


def test_lookup_timeout_tries():
    look = RSLookup(build="GRCh37", request_timeout=0.0001, request_tries=5)
    with pytest.raises(ValueError):
        look['r56']
