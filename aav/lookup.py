"""
aav.lookup
~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
import requests
from werkzeug.exceptions import NotFound
from typing import Tuple


def query_ensembl(rs_id: str, build: str) -> Tuple[str, str]:
    """Get ref and alt alleles for an rs id from ensembl"""

    if build.upper() == "GRCH38":
        url_prefix = "https://"
    elif build.upper() == "GRCH37":
        url_prefix = "https://grch37."
    else:
        raise NotImplementedError

    url = "{0}rest.ensembl.org/variation/human/{1}?{2}".format(
        url_prefix,
        rs_id,
        "content-type=application/json"
    )

    response = requests.get(url)

    if not 200 <= response.status_code < 300:
        try:
            error = response.json().get('error')
        except ValueError:
            pass
        else:
            if "not found for human" in error:
                raise NotFound("rsID not found for human")
        raise ValueError("Request failed with code {0}".format(
            response.status_code
        ))

    j = response.json()
    return j.get("ancestral_allele", ""), j.get("minor_allele", "")


class RSLookup(object):
    """
    Object to look up ref and alt positions for rs ids
    Only performs requests to ensembl when rs ids has not been
    accessed before.

    Behaves like dict.
    """

    def __init__(self, build: str,
                 init_d: dict = None):
        """
        Create lookup table.
        :param build: genome build. Either GRCH37 or GRCH38
        :param init_d: Optional dict with known rs ids
        """
        self.build = build
        if init_d:
            self.__rsids = init_d
        else:
            self.__rsids = {}

    def __getitem__(self, rs_id: str):
        if rs_id not in self.__rsids:
            self.__rsids[rs_id] = query_ensembl(rs_id, self.build)

        return self.__rsids[rs_id]
