"""
test_variation.py
~~~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
import pytest
from aav.variation import Variant


vcf_test_data = [
    (
        ["1", 100, "A", ["T"], 500],
        "1\t100\t.\tA\tT\t500\tPASS"
    ),
    (
        ["1", 100, "A", ["T"], 500, ["DIDNOTPASS"]],
        "1\t100\t.\tA\tT\t500\tDIDNOTPASS"
    ),
    (
        ["1", 100, "A", ["T"], 500, ["DIDNOTPASS"], "rsUnknown"],
        "1\t100\trsUnknown\tA\tT\t500\tDIDNOTPASS"
    ),
    (
        ["1", 100, "A", ["T,C"], 500],
        "1\t100\t.\tA\tT,C\t500\tPASS"
    ),
    (
        ["1", 100, "A", ["T,C"], 500, ["DN1", "DN2"]],
        "1\t100\t.\tA\tT,C\t500\tDN1,DN2"
    )
]


@pytest.mark.parametrize("vcf_args, expected_line", vcf_test_data)
def test_vcf_line(vcf_args, expected_line):
    a = Variant(*vcf_args)
    assert a.vcf_line == expected_line
