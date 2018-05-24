"""
aav.cli
~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""

import click
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
import json

from .readers import autodetect_reader
from .lookup import RSLookup


def green_message(msg: str) -> None:
    pre_msg = "[ {0} ]  ".format(
        datetime.utcnow().replace(tzinfo=timezone.utc).isoformat()
    )
    click.echo(click.style(pre_msg + msg, fg="green"), err=True)


@click.command()
@click.option("-p", "--path",
              type=click.Path(exists=True, readable=True),
              required=True,
              help="Path to array file")
@click.option("-b", "--build",
              type=click.Choice(["GRCh37", "GRCh38"]),
              help="")
@click.option("-s", "--sample-name",
              type=click.STRING,
              help="Name of sample in VCF file")
@click.option("-c", "--chr-prefix",
              type=click.STRING,
              required=False,
              help="Optional prefix to chromosome names")
@click.option("-l", "--lookup-table",
              type=click.Path(exists=True, readable=True),
              required=False,
              help="Optional path to existing lookup table for rsIDs.")
@click.option("-d", "--dump",
              type=click.Path(writable=True),
              required=False,
              help="Optional path to write generated lookup table")
def convert(path: str, build: str, sample_name: str,
            chr_prefix: Optional[str],
            lookup_table: Optional[str], dump: Optional[str]):
    true_path = Path(path)
    try:
        reader_cls = autodetect_reader(true_path)
    except NotImplementedError:
        raise click.FileError("Could not detect type of array")
    else:
        green_message(
            "Detected array file with type: {0}".format(reader_cls.__name__)
        )

    if lookup_table is None:
        rs_look = RSLookup(build=build)
        init_el = 0
    else:
        with Path(lookup_table).open() as lhandle:
            init_d = json.load(lhandle)
        init_el = len(init_d.keys())
        rs_look = RSLookup(build=build, init_d=init_d)

    green_message(
        "Initialized lookup table with {0} initial elements".format(init_el)
    )

    green_message("Start conversion")

    if chr_prefix is None:
        reader = reader_cls(true_path, lookup_table=rs_look)
    else:
        reader = reader_cls(true_path, lookup_table=rs_look,
                            prefix_chr=chr_prefix)

    print(reader.vcf_header(sample_name), end='')

    i = 0

    for record in reader:
        print(record.vcf_line)
        i += 1

    green_message("Converted {0} records".format(i))

    if dump is not None:
        click.echo(click.style("Dumping lookup table"), err=True)
        with Path(dump).open("w") as dhandle:
            dhandle.write(rs_look.dumps())
