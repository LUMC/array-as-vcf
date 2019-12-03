"""
aav.cli
~~~~~~~

:copyright: (c) 2018 Sander Bollen
:copyright: (c) 2018 Leiden University Medical Center
:license: MIT
"""

import argparse
import logging

from .readers import autodetect_reader, OpenArrayReader
from .lookup import RSLookup


def get_parser():
    """ Argument parsing """
    parser = argparse.ArgumentParser(
        description="Convert an array file to VCF format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--path", "-p", required=True,
                        help="Path to array file")
    parser.add_argument("--build", "-b", choices=["GRCh37", "GRCh38"],
                        default="GRCh37", help="Genome build")
    parser.add_argument("--sample-name", "-s", required=True,
                        help="Name of sample in VCF file")
    parser.add_argument("--chr-prefix", "-c", required=False,
                        help="Prefix to chromosome names")
    parser.add_argument("--lookup-table", "-l", required=False,
                        help="Path to existing lookup table for rsIDs")
    parser.add_argument("--dump", "-d", required=False,
                        help="Path to write generated lookup table")
    parser.add_argument("--encoding", default="UTF-8",
                        help="Encoding of the array file")
    parser.add_argument("--exclude-assays", nargs='+', required=False,
                        help="Assay IDs for OpenArray to ignore")
    parser.add_argument("--no-ensembl-lookup", action="store_true",
                        help="Lookup missing rsIDs on Ensembl")
    parser.add_argument("--log-level", default="INFO", required=False,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Set the verbosity of the logger")
    return parser


def convert():
    parser = get_parser()
    args = parser.parse_args()
    ensembl_lookup = not args.no_ensembl_lookup

    # Set up logging
    num_level = getattr(logging, args.log_level)
    log = logging.getLogger()
    log.setLevel(num_level)

    reader_cls = autodetect_reader(args.path, encoding=args.encoding)
    logging.info(f"Detected array file with type: {reader_cls.__name__}")

    if args.lookup_table is None:
        rs_look = RSLookup(build=args.build, ensembl_lookup=ensembl_lookup)
    else:
        rs_look = RSLookup.from_path(args.lookup_table, build=args.build,
                                     ensembl_lookup=ensembl_lookup)

    logging.info(f"Initialized lookup table with {len(rs_look)} elements.")

    logging.info(f"Start conversion.")

    if reader_cls == OpenArrayReader:
        reader = reader_cls(args.path, lookup_table=rs_look,
                            sample=args.sample_name,
                            prefix_chr=args.chr_prefix,
                            encoding=args.encoding,
                            exclude_assays=args.exclude_assays)
    else:
        reader = reader_cls(args.path, lookup_table=rs_look,
                            prefix_chr=args.chr_prefix, encoding=args.encoding)

    print(reader.vcf_header(args.sample_name), end='')

    # To print a valid vcf file, the Variants have to be sorted
    variants = sorted(reader)
    for i, record in enumerate(variants, 1):
        print(record.vcf_line)

    try:
        logging.info("Converted {0} records.".format(i))
    except UnboundLocalError:  # if there were 0 records, i is unset
        logging.info("Converted 0 records.")

    if args.dump is not None:
        logging.info("Dumping lookup table.")
        with open(args.dump, "w") as dhandle:
            dhandle.write(rs_look.dumps())
