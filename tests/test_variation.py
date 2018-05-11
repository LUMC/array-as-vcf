"""
test_variation.py
~~~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
import pytest
from aav.variation import Variant, InfoField, InfoFieldNumber, Genotype


info_test_data = [
    (
        ["INTFIELD", 100, InfoFieldNumber.one],
        "INTFIELD=100"
    ),
    (
        ["FLOATFIELD", 100.0, InfoFieldNumber.one],
        "FLOATFIELD=100.0"
    ),
    (
        ["CHARFIELD", "C", InfoFieldNumber.one],
        "CHARFIELD=C"
    ),
    (
        ["STRINGFIELD", "hello", InfoFieldNumber.one],
        "STRINGFIELD=hello"
    ),
    (
        ["FLAG", True, InfoFieldNumber.one, True],
        "FLAG"
    ),
    (
        ["FLAG", False, InfoFieldNumber.one, True],
        ""
    ),
    (
        ["MULTIINTFIELD", [100, 100], InfoFieldNumber.A],
        "MULTIINTFIELD=100,100"
    ),
    (
        ["MULTIFLOATFIELD", [100.0, 100.0], InfoFieldNumber.A],
        "MULTIFLOATFIELD=100.0,100.0"
    ),
    (
        ["MULTICHARFIELD", ["C", "C"], InfoFieldNumber.A],
        "MULTICHARFIELD=C,C"
    ),
    (
        ["MULTISTRINGFIELD", ["foo", "bar"], InfoFieldNumber.A],
        "MULTISTRINGFIELD=foo,bar"
    )
]


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
    ),
    (
        ["1", 100, "A", ["T,C"], 500, ["DN1", "DN2"], "rsUnknown",
         [InfoField("FOO", "bar", InfoFieldNumber.one)]],
        "1\t100\trsUnknown\tA\tT,C\t500\tDN1,DN2\tFOO=bar"
    ),
    (
        ["1", 100, "A", ["T,C"], 500, ["DN1", "DN2"], "rsUnknown",
         [
             InfoField("FOO", "bar", InfoFieldNumber.one),
             InfoField("BAZ", True, InfoFieldNumber.one, True)
         ]],
        "1\t100\trsUnknown\tA\tT,C\t500\tDN1,DN2\tFOO=bar;BAZ"
    ),
    (
        ["1", 100, "A", ["T,C"], 500, ["DN1", "DN2"], "rsUnknown",
         [
             InfoField("FOO", "bar", InfoFieldNumber.one),
             InfoField("BAZ", True, InfoFieldNumber.one, True)
         ],
         Genotype.hom_ref
         ],
        "1\t100\trsUnknown\tA\tT,C\t500\tDN1,DN2\tFOO=bar;BAZ\tGT\t0/0"
    ),
    (
        ["1", 100, "A", ["T,C"], 500, ["DN1", "DN2"], "rsUnknown",
         [
             InfoField("FOO", "bar", InfoFieldNumber.one),
             InfoField("BAZ", True, InfoFieldNumber.one, True)
         ],
         Genotype.het
         ],
        "1\t100\trsUnknown\tA\tT,C\t500\tDN1,DN2\tFOO=bar;BAZ\tGT\t0/1"
    ),
    (
        ["1", 100, "A", ["T,C"], 500, ["DN1", "DN2"], "rsUnknown",
         [
             InfoField("FOO", "bar", InfoFieldNumber.one),
             InfoField("BAZ", True, InfoFieldNumber.one, True)
         ],
         Genotype.hom_alt
         ],
        "1\t100\trsUnknown\tA\tT,C\t500\tDN1,DN2\tFOO=bar;BAZ\tGT\t1/1"
    ),
    (
        ["1", 100, "A", ["T,C"], 500, ["DN1", "DN2"], "rsUnknown",
         [
             InfoField("FOO", "bar", InfoFieldNumber.one),
             InfoField("BAZ", True, InfoFieldNumber.one, True)
         ],
         Genotype.unknown
         ],
        "1\t100\trsUnknown\tA\tT,C\t500\tDN1,DN2\tFOO=bar;BAZ\tGT\t./."
    )
]


@pytest.mark.parametrize("info_args, expected_str", info_test_data)
def test_info_str(info_args, expected_str):
    a = InfoField(*info_args)
    assert str(a) == expected_str


def test_info_flag_exc():
    with pytest.raises(ValueError):
        InfoField("FLAG", "notaboolean", InfoFieldNumber.one, True)


@pytest.mark.parametrize("vcf_args, expected_line", vcf_test_data)
def test_vcf_line(vcf_args, expected_line):
    a = Variant(*vcf_args)
    assert a.vcf_line == expected_line
