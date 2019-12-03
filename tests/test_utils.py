"""
test_utils.py
~~~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center

:license: MIT
"""
import pytest
from array_as_vcf.utils import comma_float


comma_float_params = [
    ("0,001", 0.001),
    ("0.001", 0.001),
    ("1,34", 1.34),
    ("1.34", 1.34),
    ("500", 500),
    ("500,0", 500),
    ("500.0", 500)
]


@pytest.mark.parametrize("str_val, float_val", comma_float_params)
def test_comma_float(str_val, float_val):
    assert comma_float(str_val) == float_val


def test_comma_float_err():
    with pytest.raises(ValueError):
        comma_float("5,6.0")
