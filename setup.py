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
    name="array-as-vcf",
    version="1.0.0",
    description="Convert SNP array to VCF",
    long_description=long_desc,
    long_description_content_type="text/markdown",
    author="Sander Bollen, Redmar van den Berg",
    author_email="KG_bioinf@lumc.nl",
    license="MIT",
    keywords="array vcf SNP convert OpenArray Cytoscan lumi370k lumi317k Affymetrix",
    zip_safe=False,
    url="https://github.com/LUMC/array-as-vcf",
    package_dir={"": "src"},
    packages=find_packages("src"),
    install_requires=["requests", "setuptools"],
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "array-as-vcf = array_as_vcf.cli:convert",
            "aav = array_as_vcf.cli:convert",
        ]
    },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
