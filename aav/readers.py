"""
aav.readers
~~~~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""
from pathlib import Path

from .variation import Variant


class Reader(object):
    """
    Generic reader object

    Readers are iterators that produce variants
    """
    def __init__(self, path: Path):
        self.path = path

    def __next__(self) -> Variant:
        raise NotImplementedError

    def __iter__(self):
        return self
