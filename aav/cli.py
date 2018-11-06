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
from typing import Optional, Set

from .readers import autodetect_reader, OpenArrayReader
from .lookup import RSLookup


def green_message(msg: str) -> None:
    pre_msg = "[ {0} ]  ".format(
        datetime.utcnow().replace(tzinfo=timezone.utc).isoformat()
    )
    click.echo(click.style(pre_msg + msg, fg="green"), err=True)


def exlude_assays_callback(ctx, param, value) -> Optional[Set[str]]:
    if value is None:
        return None

    return set(value.split(","))


@click.command()
@click.option("-p", "--path",
              type=click.Path(exists=True, readable=True),
              required=True,
              help="Path to array file")
@click.option("-b", "--build",
              type=click.Choice(["GRCh37", "GRCh38"]),
              help="Genome build. Default = GRCh37", default="GRCh37")
@click.option("-s", "--sample-name",
              type=click.STRING,
              help="Name of sample in VCF file",
              required=True)
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
@click.option("--encoding", type=click.STRING, required=False,
              help="Optional encoding of array file. "
                   "Encoding defaults to UTF-8 if not given")
@click.option("--exclude-assays", type=click.STRING, required=False,
              callback=exlude_assays_callback,
              help="Optional comma-separated list of assay IDs "
                   "for OpenArray to ignore")
def convert(path: str, build: str, sample_name: str,
            chr_prefix: Optional[str],
            lookup_table: Optional[str], dump: Optional[str],
            encoding: Optional[str],
            exclude_assays: Optional[Set[str]]):
    true_path = Path(path)
    try:
        reader_cls = autodetect_reader(true_path, encoding=encoding)
    except NotImplementedError:
        raise click.FileError("Could not detect type of array.")
    else:
        green_message(
            "Detected array file with type: {0}.".format(reader_cls.__name__)
        )

    if lookup_table is None:
        rs_look = RSLookup(build=build)
    else:
        rs_look = RSLookup.from_path(Path(lookup_table), build=build)

    green_message(
        f"Initialized lookup table with {len(rs_look)} initial elements."
    )

    green_message("Start conversion.")

    if reader_cls == OpenArrayReader:
        reader = reader_cls(true_path, lookup_table=rs_look,
                            sample=sample_name, prefix_chr=chr_prefix,
                            encoding=encoding,
                            exclude_assays=exclude_assays)
    else:
        reader = reader_cls(true_path, lookup_table=rs_look,
                            prefix_chr=chr_prefix, encoding=encoding)

    print(reader.vcf_header(sample_name), end='')

    i = 0

    for record in reader:
        print(record.vcf_line)
        i += 1

    green_message("Converted {0} records.".format(i))

    if dump is not None:
        green_message("Dumping lookup table.")
        with Path(dump).open("w") as dhandle:
            dhandle.write(rs_look.dumps())
