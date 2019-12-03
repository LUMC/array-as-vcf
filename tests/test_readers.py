"""
test_readers.py
~~~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center

:license: MIT
"""
from datetime import date
from pathlib import Path

import pytest

from array_as_vcf import __version__
from array_as_vcf.readers import (AffyReader, CytoScanReader,
                                  Lumi317kReader, Lumi370kReader,
                                  autodetect_reader, Reader,
                                  OpenArrayReader)
from array_as_vcf.lookup import RSLookup
from array_as_vcf.variation import Genotype
from array_as_vcf.variation import Variant


def grch37_lookup():
    return RSLookup(build="GRCh37", request_tries=3)


def grch37_lookup_no_ensembl():
    return RSLookup(build="GRCh37", ensembl_lookup=False)


def test_lookup_table():
    lookup_path = str(Path(__file__).parent / Path("data") /
                      Path("lookup_table_test.json"))
    return RSLookup.from_path(lookup_path, build="GRCh37",
                              ensembl_lookup=False)


_lumi_317_path = str(Path(__file__).parent / Path("data") /
                     Path("lumi_317_test.txt"))
_lumi_370_path = str(Path(__file__).parent / Path("data") /
                     Path("lumi_370_test.txt"))
_affy_path = str(Path(__file__).parent / Path("data") /
                 Path("affy_test.txt"))
_cytoscan_path = str(Path(__file__).parent / Path("data") /
                     Path("cytoscan_test.txt"))
_open_array_path = str(Path(__file__).parent / Path("data") /
                       Path("open_array_test.txt"))
_open_array_all_calls = str(Path(__file__).parent / Path("data") /
                            Path("open_array_all_calls.txt"))


@pytest.fixture
def lumi_317_reader():
    return Lumi317kReader(_lumi_317_path, grch37_lookup())


@pytest.fixture
def lumi_370_reader():
    return Lumi370kReader(_lumi_370_path, grch37_lookup())


@pytest.fixture
def affy_reader():
    return AffyReader(_affy_path, grch37_lookup())


@pytest.fixture
def cytoscan_reader():
    return CytoScanReader(_cytoscan_path, grch37_lookup())


@pytest.fixture
def open_array_reader():
    return OpenArrayReader(_open_array_path, test_lookup_table(),
                           "e31a0a96465a", encoding="windows-1252")


@pytest.fixture
def lumi_317_reader_no_ensembl():
    lookup = test_lookup_table()
    return Lumi317kReader(_lumi_317_path, lookup)


@pytest.fixture
def lumi_370_reader_no_ensembl():
    lookup = test_lookup_table()
    return Lumi370kReader(_lumi_370_path, lookup)


@pytest.fixture
def affy_reader_no_ensembl():
    lookup = test_lookup_table()
    return AffyReader(_affy_path, lookup)


@pytest.fixture
def cytoscan_reader_no_ensembl():
    lookup = test_lookup_table()
    return CytoScanReader(_cytoscan_path, lookup)


@pytest.fixture
def open_array_reader_no_ensembl():
    lookup = test_lookup_table()
    return OpenArrayReader(_open_array_path, lookup,
                           "e31a0a96465a", encoding="windows-1252")


@pytest.fixture
def open_array_reader_no_ensembl_lookup():
    lookup = test_lookup_table()
    return OpenArrayReader(_open_array_path, lookup,
                           "e31a0a96465a", encoding="windows-1252")


@pytest.fixture
def open_array_reader_all_calls():
    lookup = test_lookup_table()
    return OpenArrayReader(_open_array_all_calls, lookup, 'all_calls')


chrom_test_data = [
    (
        AffyReader(_affy_path, grch37_lookup()),
        ["1"]*8
    ),
    (
        AffyReader(_affy_path, grch37_lookup(), prefix_chr="chr"),
        ["chr1"]*8
    ),
    (
        CytoScanReader(_cytoscan_path, grch37_lookup()),
        ["1", "1", "1"]
    ),
    (
        CytoScanReader(_cytoscan_path, grch37_lookup(), prefix_chr="chr"),
        ["chr1", "chr1", "chr1"]
    ),
    (
        Lumi317kReader(_lumi_317_path, grch37_lookup()),
        ["1", "1", "1", "1"]
    ),
    (
        Lumi317kReader(_lumi_317_path, grch37_lookup(), prefix_chr="chr"),
        ["chr1", "chr1", "chr1", "chr1"]
    ),
    (
        Lumi370kReader(_lumi_370_path, grch37_lookup()),
        ["1", "1", "1", "1"]
    ),
    (
        Lumi370kReader(_lumi_370_path, grch37_lookup(), prefix_chr="chr"),
        ["chr1", "chr1", "chr1", "chr1"]
    ),
    (
        OpenArrayReader(_open_array_path, test_lookup_table(),
                        prefix_chr="chr", sample="e31a0a96465a",
                        encoding="windows-1252"),
        ["chr3", "chr19", "chr15", "chr4", "chr13", "chr21", "chr19",
         "chr12", "chr4", "chr18", "chr14", "chrY", "chr4", "chr19",
         "chr6", "chr13", "chr18", "chr11", "chr1", "chr8", "chr15", "chr1",
         "chr20", "chr15", "chrY", "chr18", "chr8", "chr18", "chr17",
         "chr13", "chr14", "chr6", "chr17", "chr13", "chr10", "chr19",
         "chr2", "chrY", "chr6", "chr21", "chr18", "chr14", "chr6", "chr8",
         "chr18", "chr2", "chr3", "chr22", "chr17", "chr1", "chr10",
         "chr10", "chr3", "chrY", "chr18"]
    )
]

ref_test_data = [
    (
        AffyReader(_affy_path, grch37_lookup()),
        ["T", "A", "T", "T", "T", "A", "T", "T"]
    ),
    (
        CytoScanReader(_cytoscan_path, grch37_lookup()),
        ["A", "T", "C"]
    ),
    (
        Lumi317kReader(_lumi_317_path, grch37_lookup()),
        ["C", "A", "C", "A"]
    ),
    (
        Lumi370kReader(_lumi_370_path, grch37_lookup()),
        ["C", "A", "C", "A"]
    ),
    (
        OpenArrayReader(_open_array_path, test_lookup_table(),
                        sample="e31a0a96465a", encoding="windows-1252"),
        ["A", "G", "G", "T", "A", "G", "G", "C", "T", "T", "T", "C", "T",
         "G", "G", "A", "A", "G", "C", "T", "G", "C", "G", "A", "T", "T",
         "T", "G", "G", "A", "T", "G", "G", "A", "C", "G", "T", "C", "T",
         "G", "A", "T", "T", "G", "G", "A", "A", "T", "T", "T", "C", "G",
         "A", "C", "A"]
    )
]

alt_test_data = [
    (
        AffyReader(_affy_path, grch37_lookup()),
        [["C"], ["C"], ["C"], ["A", "G"], ["C"],
         ["C"], ["C"], ["A", "G"]]
    ),
    (
        CytoScanReader(_cytoscan_path, grch37_lookup()),
        [["G"], ["C"], ["G", "T"]]
    ),
    (
        Lumi317kReader(_lumi_317_path, grch37_lookup()),
        [["T"], ["C", "G", "T"], ["T"], ["G"]]
    ),
    (
        Lumi370kReader(_lumi_370_path, grch37_lookup()),
        [["T"], ["C", "G", "T"], ["T"], ["G"]]
    )
]

autodetect_reader_data = [
    (_lumi_317_path, "Lumi317kReader", None),
    (_lumi_370_path, "Lumi370kReader", None),
    (_affy_path, "AffyReader", None),
    (_cytoscan_path, "CytoScanReader", None),
    (_open_array_path, "OpenArrayReader", 'windows-1252')
]


def test_affy_reader_amount(affy_reader):
    assert len(list(affy_reader)) == 8


def test_cytoscan_reader_amount(cytoscan_reader):
    assert len(list(cytoscan_reader)) == 3


def test_lumi317_reader_amount(lumi_317_reader):
    assert len(list(lumi_317_reader)) == 4


def test_lumi370_reader_amount(lumi_370_reader):
    assert len(list(lumi_370_reader)) == 4


def test_open_array_reader_amount(open_array_reader):
    assert len(list(open_array_reader)) == 55


def test_reader_genotypes(affy_reader, cytoscan_reader, lumi_317_reader,
                          lumi_370_reader, open_array_reader):
    genotype_test_data = [
        (
            affy_reader,
            [Genotype.unknown, Genotype.unknown, Genotype.het,
             Genotype.unknown, Genotype.unknown, Genotype.unknown,
             Genotype.het, Genotype.unknown]
        ),
        (
            cytoscan_reader,
            [Genotype.hom_alt, Genotype.hom_alt, Genotype.het]
        ),
        (
            lumi_317_reader,
            [Genotype.hom_alt, Genotype.hom_alt, Genotype.het,
             Genotype.unknown]
        ),
        (
            lumi_370_reader,
            [Genotype.hom_alt, Genotype.hom_alt, Genotype.het,
             Genotype.unknown]
        ),
        (
            open_array_reader,
            [Genotype.unknown, Genotype.het, Genotype.hom_alt,
             Genotype.hom_ref, Genotype.hom_ref, Genotype.hom_ref,
             Genotype.hom_ref, Genotype.hom_ref, Genotype.hom_alt,
             Genotype.hom_alt, Genotype.hom_alt, Genotype.hom_ref,
             Genotype.hom_alt, Genotype.hom_alt, Genotype.het,
             Genotype.het, Genotype.hom_alt, Genotype.het, Genotype.hom_alt,
             Genotype.unknown, Genotype.het, Genotype.het, Genotype.het,
             Genotype.hom_alt, Genotype.unknown, Genotype.het, Genotype.het,
             Genotype.hom_ref, Genotype.het, Genotype.het, Genotype.hom_alt,
             Genotype.unknown, Genotype.het, Genotype.het, Genotype.hom_ref,
             Genotype.unknown, Genotype.hom_ref, Genotype.hom_ref,
             Genotype.het, Genotype.het, Genotype.hom_alt, Genotype.hom_alt,
             Genotype.het, Genotype.hom_alt, Genotype.het, Genotype.hom_ref,
             Genotype.hom_alt, Genotype.hom_alt, Genotype.hom_ref,
             Genotype.het, Genotype.hom_ref, Genotype.het,
             Genotype.het, Genotype.unknown, Genotype.hom_ref]
        )
    ]
    for reader, genotypes in genotype_test_data:
        found_gt = [x.genotype for x in reader]
        assert found_gt == genotypes


@pytest.mark.parametrize("path, class_name, encoding", autodetect_reader_data)
def test_autodetect_readers(path, class_name, encoding):
    if encoding is None:
        found_cls = autodetect_reader(path)
    else:
        found_cls = autodetect_reader(path, encoding=encoding)
    assert found_cls.__name__ == class_name


@pytest.mark.parametrize("reader, chroms", chrom_test_data)
def test_reader_chromosomes(reader, chroms):
    found_chroms = [x.chrom for x in reader]
    assert found_chroms == chroms


@pytest.mark.parametrize("reader, refs", ref_test_data)
def test_reader_refs(reader, refs):
    found_ref = [x.ref for x in reader]
    assert found_ref == refs


@pytest.mark.parametrize("reader, alts", alt_test_data)
def test_reader_alts(reader, alts):
    found_alt = [x.alt for x in reader]
    assert found_alt == alts


def test_base_reader_header():
    reader = Reader(_lumi_370_path)
    header = reader.vcf_header("sample_01")
    date_str = date.today().isoformat()
    exp = (
        "##fileformat=VCFv4.2\n"
        "##fileDate={0}\n"
        "##source=aav_v{1}\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_01\n"
           ).format(date_str, __version__)
    assert header == exp


def test_affy_reader_header(affy_reader):
    header = affy_reader.vcf_header("sample_01")
    date_str = date.today().isoformat()
    exp = (
        "##fileformat=VCFv4.2\n"
        "##fileDate={0}\n"
        "##source=aav_v{1}\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##INFO=<ID=ID,Number=1,Type=String,Description="A field">\n'
        '##INFO=<ID=AffymetrixSNPsID,Number=1,Type=String,Description="A field">\n'  # noqa
        '##INFO=<ID=log2ratio_AB,Number=1,Type=Float,Description="A field">\n'
        '##INFO=<ID=N_AB,Number=1,Type=Integer,Description="A field">\n'
        '##INFO=<ID=LOH_likelihood,Number=1,Type=Float,Description="A field">\n'  # noqa
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_01\n"
           ).format(date_str, __version__)
    assert header == exp


def test_cytoscan_reader_header(cytoscan_reader):
    header = cytoscan_reader.vcf_header("sample_01")
    date_str = date.today().isoformat()
    exp = (
        "##fileformat=VCFv4.2\n"
        "##fileDate={0}\n"
        "##source=aav_v{1}\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##INFO=<ID=Probe_Set_ID,Number=1,Type=String,Description="A field">\n'
        '##INFO=<ID=Signal_A,Number=1,Type=Float,Description="A field">\n'
        '##INFO=<ID=Signal_B,Number=1,Type=Float,Description="A field">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_01\n"
           ).format(date_str, __version__)
    assert header == exp


def test_lumi_reader_header(lumi_317_reader):
    header = lumi_317_reader.vcf_header("sample_01")
    date_str = date.today().isoformat()
    exp = (
        "##fileformat=VCFv4.2\n"
        "##fileDate={0}\n"
        "##source=aav_v{1}\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##INFO=<ID=Log_R_Ratio,Number=1,Type=Float,Description="A field">\n'
        '##INFO=<ID=CNV_Value,Number=1,Type=Integer,Description="A field">\n'
        '##INFO=<ID=Allele_Freq,Number=1,Type=Float,Description="A field">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_01\n"
           ).format(date_str, __version__)
    assert header == exp


def test_open_array_reader_header(open_array_reader):
    header = open_array_reader.vcf_header("e31a0a96465a")
    date_str = date.today().isoformat()
    exp = (
        "##fileformat=VCFv4.2\n"
        "##fileDate={0}\n"
        "##source=aav_v{1}\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##INFO=<ID=Assay_Name,Number=1,Type=String,Description="A field">\n'
        '##INFO=<ID=Assay_ID,Number=1,Type=String,Description="A field">\n'
        '##INFO=<ID=Gene_Symbol,Number=.,Type=String,Description="A field">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\te31a0a96465a\n"
    ).format(date_str, __version__)
    assert header == exp


def test_lumi317_unknown_rsid(lumi_317_reader_no_ensembl):
    """ Variants with missing REF fields are not allowed according to the VCF
    standard
    """
    for variant in lumi_317_reader_no_ensembl:
        assert variant.ref != '.'


def test_lumi370_unknown_rsid(lumi_370_reader_no_ensembl):
    """ Variants with missing REF fields are not allowed according to the VCF
    standard
    """
    for variant in lumi_370_reader_no_ensembl:
        assert variant.ref != '.'


def test_openarray_unknown_rsid(open_array_reader_no_ensembl):
    """ Variants with missing REF fields are not allowed according to the VCF
    standard
    """
    for variant in open_array_reader_no_ensembl:
        assert variant.ref != '.'


def test_cytoscan_unknown_rsid(cytoscan_reader_no_ensembl):
    """ Variants with missing REF fields are not allowed according to the VCF
    standard
    """
    for variant in cytoscan_reader_no_ensembl:
        assert variant.ref != '.'


def test_affy_unknown_rsid(affy_reader_no_ensembl):
    """ Variants with missing REF fields are not allowed according to the VCF
    standard
    """
    for variant in affy_reader_no_ensembl:
        assert variant.ref != '.'


def variants_are_ordered(iterator):
    """ Return whether the Variants in iterator are sorted by CHROM and
    POS, so that they can be written to a valid vcf file
    """
    chroms = set()
    # Create an empty variant
    prev_variant = Variant('None', -1, 'A', ['C'], -1.0, ['filters'])
    for variant in iterator:
        if variant.chrom != prev_variant.chrom:
            # If we see this chromosome for the second time, the ordering is
            # wrong
            if variant.chrom in chroms:
                return False
            # Otherwise, add the chromosome to the set
            else:
                chroms.add(variant.chrom)
        # If the CHROM is the same, POS can not be smaller than the previous
        # position
        elif variant.pos < prev_variant.pos:
            return False
        prev_variant = variant
    else:
        return True


def test_open_array_is_sorted(open_array_reader_no_ensembl):
    """ The OpenArray file is not sorted, unlike the other array files. So we
    can use the OpenArray file to make sure we sort the variants properly

    Note that this is only required if you want to write the Variants from an
    array into a valid VCF file
    """
    variants = sorted(open_array_reader_no_ensembl)
    assert variants_are_ordered(variants) is True


def test_open_array_all_calls(open_array_reader_all_calls):
    """ All these calls should have './.' as genotype since we cannot determine
    the actual genotype from the "call" column
    """
    genotypes = [var.genotype for var in open_array_reader_all_calls]
    assert genotypes == [Genotype.unknown]*6
