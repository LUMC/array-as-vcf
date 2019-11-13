"""
setup.py
~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""

from os.path import abspath, dirname, join

from setuptools import setup

from aav import __version__, __author__

readme_file = join(abspath(dirname(__file__)), "README.md")
with open(readme_file) as desc_handle:
    long_desc = desc_handle.read()

setup(
    name="aav",
    version=__version__,
    description="Array to VCF",
    long_description=long_desc,
    author=__author__,
    author_email="a.h.b.bollen@lumc.nl",
    license="MIT",
    packages=["aav"],
    install_requires=[],
    entry_points={
        "console_scripts": [
            "aav = aav.cli:convert"
        ]
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
