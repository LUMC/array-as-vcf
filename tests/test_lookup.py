"""
test_lookup.py
~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""

from aav.lookup import query_ensembl

import pytest
from werkzeug.exceptions import NotFound


def test_query_ensembl_known_rsid_grch38():
    assert query_ensembl("rs56116432", "GRCh38") == ("C", "T")


def test_query_ensembl_known_rsid_grch37():
    assert query_ensembl("rs56116432", "GRCh37") == ("C", "T")


def test_query_ensembl_unknown_build():
    with pytest.raises(NotImplementedError):
        query_ensembl("rs56116432", "unknown")


def test_unknown_rsid():
    with pytest.raises(NotFound) as exc:
        query_ensembl("rs5611644432", "GRCh37")
        assert "rsID not found for human" in str(exc)
