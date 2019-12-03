"""
aav.utils
~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from typing import Optional


def comma_float(val: str) -> float:
    """
    Get float for a string that may contain commas in stead of dots
    :param val: the value to be casted to float
    :return: float
    """
    if "," in val and "." not in val:
        return float(val.replace(",", "."))
    elif "." in val and "," not in val:
        return float(val)
    elif "," not in val and "." not in val:
        return float(val)
    else:
        raise ValueError("Cannot parse string with both commas and dots")


def empty_string(val: Optional[str]) -> bool:
    """Check whether a field is an empty or non-existing string"""
    return val is None or val == ""
