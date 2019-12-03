"""
setup.py
~~~~~~~~

:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""

from os.path import abspath, dirname, join

from setuptools import setup, find_packages

readme_file = join(abspath(dirname(__file__)), "README.md")
with open(readme_file) as desc_handle:
    long_desc = desc_handle.read()

setup(
    name="array_as_vcf",
    version="1.0.0-dev",
    description="Array to VCF",
    long_description=long_desc,
    author="Sander Bollen, Redmar van den Berg",
    author_email="KG_bioinf@lumc.nl",
    license="MIT",
    package_dir={'': 'src'},
    packages=find_packages(),
    install_requires=["requests", "setuptools"],
    entry_points={
        "console_scripts": [
            "array-as-vcf = array_as_vcf.cli:convert",
            "aav = array_as_vcf.cli:convert"
        ]
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
