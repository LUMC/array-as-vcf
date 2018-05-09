"""
test_variation.py
~~~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""

from aav.variation import Variant


def test_vcf_line_no_id_no_filter():
    a = Variant("1", 100, "A", ["T"], 500)
    assert a.vcf_line == "1\t100\t.\tA\tT\t500\t"


def test_vcf_line_no_id_with_filter():
    a = Variant("1", 100, "A", ["T"], 500, ["DIDNOTPASS"])
    assert a.vcf_line == "1\t100\t.\tA\tT\t500\tDIDNOTPASS"


def test_vcf_line_with_id_with_filter():
    a = Variant("1", 100, "A", ["T"], 500, ["DIDNOTPASS"], "rsUnknown")
    assert a.vcf_line == "1\t100\trsUnknown\tA\tT\t500\tDIDNOTPASS"
