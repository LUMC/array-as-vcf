"""
aav.variation
~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
import enum
from typing import List, Optional, Any
from datetime import date

from . import __version__


class InfoFieldNumber(enum.Enum):
    one = "1"
    A = "A"
    R = "R"
    G = "G"
    unknown = "."


class InfoFieldType(enum.Enum):
    STRING = "String"
    FLOAT = "Float"
    INT = "Integer"
    FLAG = "Flag"


class Genotype(enum.Enum):
    hom_ref = "0/0"
    het = "0/1"
    hom_alt = "1/1"
    unknown = "./."


class InfoField(object):
    """Info field"""

    def __init__(self, name: str, value: Any,
                 number: InfoFieldNumber, flag: bool = False):
        if flag and not isinstance(value, bool):
            raise ValueError("Value must be boolean if a flag")

        self.name = name
        self.value = value
        self.number = number
        self.flag = flag

    def __str__(self) -> Optional[str]:
        if self.flag and self.value:
            return self.name
        elif self.flag:
            return ""
        elif self.number == InfoFieldNumber.one:
            return "{0}={1}".format(self.name, self.value)
        else:
            return "{0}={1}".format(self.name, ",".join(map(str, self.value)))


class Variant(object):
    """
    Representations of a variant

    This currently only supports _one_ sample,
    with _one_ FORMAT field entry (GT).
    """
    def __init__(self, chrom: str, pos: int, ref: str, alt: List[str],
                 qual: float, filters: List[str] = list(),
                 id: Optional[str] = None,
                 info_fields: List[InfoField] = list(),
                 genotype: Optional[Genotype] = None):
        self.chrom = chrom
        self.pos = pos
        if id is not None:
            self.id = id
        else:
            self.id = "."
        self.qual = qual
        self.ref = ref
        self.alt = alt
        self.filters = filters
        self.info_fields = info_fields
        self.genotype = genotype

    @property
    def vcf_line(self) -> str:
        fmt = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(
            self.chrom, self.pos, self.id,
            self.ref, ",".join(self.alt),
            self.qual,
            ",".join(self.filters) if len(self.filters) > 0 else "PASS"
        )
        if len(self.info_fields) > 0:
            info_str = ";".join(
                str(x) for x in self.info_fields if str(x) != ""
            )
            fmt += "\t{0}".format(info_str)

        if self.genotype is not None:
            fmt += "\tGT\t{0}".format(self.genotype.value)

        return fmt

    def __lt__(self, other) -> bool:
        return (self.chrom, self.pos) < (other.chrom, other.pos)


class HeaderLine(object):

    def __str__(self) -> str:
        return "##"


class MetaLine(HeaderLine):
    def __init__(self, key, value):
        """
        Simple meta line with key-value pairs
        """
        self.key = key
        self.value = value

    def __str__(self) -> str:
        return super().__str__() + "{0}={1}".format(self.key, self.value)


class BracketHeaderLine(HeaderLine):
    header_type = None

    def __init__(self, id: str, number: InfoFieldNumber,
                 type: InfoFieldType, description: Optional[str] = None):
        """
        An INFO or FORMAT header line
        :param id: id of the field
        :param number: number of info field. Defaults to 0 if type is FLAG
        :param type: type of info field
        :param description: Optional description of field
        """
        self.id = id
        if type == InfoFieldType.FLAG:
            self.number = "0"
        else:
            self.number = number.value
        self.type = type.value
        if description is None:
            self.description = "A field"
        else:
            self.description = description.replace(
                '\\', '\\\\'
            ).replace('"', '\\"')

    def __str__(self) -> str:
        field = '{0}=<ID={1},Number={2},Type={3},Description="{4}">'.format(
            self.header_type, self.id, self.number,
            self.type, self.description
        )
        return super().__str__() + field


class InfoHeaderLine(BracketHeaderLine):
    header_type = "INFO"


class FormatHeaderLine(BracketHeaderLine):
    header_type = "FORMAT"


def program_header() -> MetaLine:
    return MetaLine("source", "aav_v{0}".format(__version__))


def date_header() -> MetaLine:
    d = date.today().isoformat()
    return MetaLine("fileDate", d)


def chrom_header(sample_name: str) -> str:
    """Create final header line. Single sample only"""
    s = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{0}".format(
        sample_name
    )
    return s


VCF_v_4_2 = MetaLine("fileformat", "VCFv4.2")
GT_FORMAT = FormatHeaderLine("GT", InfoFieldNumber.one,
                             InfoFieldType.STRING, "Genotype")
