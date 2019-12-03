"""
test_variation.py
~~~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from datetime import date
import pytest
from array_as_vcf import __version__
from array_as_vcf.variation import (
    Variant, InfoField, InfoFieldNumber, Genotype,
    MetaLine, InfoFieldType, InfoHeaderLine,
    FormatHeaderLine, date_header, program_header,
    chrom_header)


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

header_line_data = [
    (
        MetaLine("key", "val"),
        "##key=val"
    ),
    (
        InfoHeaderLine("FOO", InfoFieldNumber.one,
                       InfoFieldType.STRING, "FOO value"),
        '##INFO=<ID=FOO,Number=1,Type=String,Description="FOO value">'  # noqa
    ),
    (
        InfoHeaderLine("FOO", InfoFieldNumber.one,
                       InfoFieldType.FLAG, "FOO value"),
        '##INFO=<ID=FOO,Number=0,Type=Flag,Description="FOO value">'  # noqa
    ),
    (
        InfoHeaderLine("FOO", InfoFieldNumber.one,
                       InfoFieldType.INT, "FOO value"),
        '##INFO=<ID=FOO,Number=1,Type=Integer,Description="FOO value">'  # noqa
    ),
    (
        InfoHeaderLine("FOO", InfoFieldNumber.one,
                       InfoFieldType.FLOAT, "FOO value"),
        '##INFO=<ID=FOO,Number=1,Type=Float,Description="FOO value">'  # noqa
    ),
    (
        InfoHeaderLine("FOO", InfoFieldNumber.A,
                       InfoFieldType.STRING, "FOO value"),
        '##INFO=<ID=FOO,Number=A,Type=String,Description="FOO value">'  # noqa
    ),
    (
        InfoHeaderLine("FOO", InfoFieldNumber.G,
                       InfoFieldType.STRING, "FOO value"),
        '##INFO=<ID=FOO,Number=G,Type=String,Description="FOO value">'  # noqa
    ),
    (
        InfoHeaderLine("FOO", InfoFieldNumber.R,
                       InfoFieldType.STRING, "FOO value"),
        '##INFO=<ID=FOO,Number=R,Type=String,Description="FOO value">'  # noqa
    ),
    (
        InfoHeaderLine("FOO", InfoFieldNumber.one,
                       InfoFieldType.STRING),
        '##INFO=<ID=FOO,Number=1,Type=String,Description="A field">'  # noqa
    ),
    (
        InfoHeaderLine("FOO", InfoFieldNumber.one,
                       InfoFieldType.STRING, '"FOO\\BAR"'),
        '##INFO=<ID=FOO,Number=1,Type=String,Description="\\"FOO\\\\BAR\\"">'  # noqa
    ),
    (
        FormatHeaderLine("FOO", InfoFieldNumber.one,
                         InfoFieldType.STRING),
        '##FORMAT=<ID=FOO,Number=1,Type=String,Description="A field">'
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


@pytest.mark.parametrize("header, expected_str", header_line_data)
def test_header_line(header, expected_str):
    assert str(header) == expected_str


def test_program_header():
    assert str(program_header()) == "##source=aav_v{0}".format(__version__)


def test_file_date():
    date_str = date.today().isoformat()
    assert str(date_header()) == "##fileDate={0}".format(date_str)


def test_chrom_header():
    e = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_01"
    assert chrom_header("sample_01") == e
