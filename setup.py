"""
setup.py
~~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""

from os.path import abspath, dirname, join

from setuptools import setup

readme_file = join(abspath(dirname(__file__)), "README.md")
with open(readme_file) as desc_handle:
    long_desc = desc_handle.read()

setup(
    name="aav",
    version="0.0.1",
    description="Array to VCF",
    long_description=long_desc,
    author="Sander Bollen",
    author_email="a.h.b.bollen@lumc.nl",
    license="MIT",
    packages=["aav"],
    install_requires=[
        "click"
    ],
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)