"""
aav.lookup
~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
import requests
from werkzeug.exceptions import NotFound
from typing import Tuple, List

import json


def query_ensembl(rs_id: str, build: str,
                  timeout: float = 120) -> Tuple[str, List[str]]:
    """
    Get ref and alt alleles for an rs id from ensembl

    :param rs_id: The rsID to query
    :param build: genome build of interest. Either Grch37 or GRCh38
    :param timeout: request timeout in seconds (default = 120)
    """

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

    response = requests.get(url, timeout=timeout)

    if not 200 <= response.status_code < 300:
        try:
            error = response.json().get('error')
        except ValueError:
            pass
        else:
            if "not found for human" in error:
                raise NotFound("rsID not found for human")
        raise requests.HTTPError("Request failed with code {0}".format(
            response.status_code
        ))

    j = response.json()
    try:
        allele_string = j.get("mappings", [{}])[0].get("allele_string", "")
    except IndexError:  # mapping may be `mapping: []` when it does not map to genome # noqa
        raise NotFound("rsID does not map to genome")
    parts = allele_string.split("/")
    ref = parts[0]
    if len(parts) > 0:
        return ref, parts[1:]
    else:
        return ref, []


class RSLookup(object):
    """
    Object to look up ref and alt positions for rs ids
    Only performs requests to ensembl when rs ids has not been
    accessed before.

    Behaves like dict.
    """

    def __init__(self, build: str,
                 init_d: dict = None,
                 request_timeout: float = 120,
                 request_tries: int = 1):
        """
        Create lookup table.
        :param build: genome build. Either GRCH37 or GRCH38
        :param init_d: Optional dict with known rs ids
        :param request_timeout: timeout in seconds for requests
        :param request_tries: number of tries for a timed-out request
        """
        self.build = build
        self.request_tries = request_tries
        self.request_timeout = request_timeout
        if init_d:
            self.__rsids = init_d
        else:
            self.__rsids = {}

    def __getitem__(self, rs_id: str):
        if rs_id not in self.__rsids:
            self.__rsids[rs_id] = self._get_ensembl(rs_id)

        return self.__rsids[rs_id]

    def _get_ensembl(self, rs_id):
        for _ in range(self.request_tries):
            try:
                return query_ensembl(rs_id, self.build, self.request_timeout)
            except (requests.Timeout, requests.HTTPError):
                continue
        raise ValueError("Too many tries for request")

    def dumps(self):
        """Dump table to json-formatted string"""
        return json.dumps(self.__rsids)
