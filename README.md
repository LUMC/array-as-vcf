# Array As VCF

`array-as-vcf` is a small library and tool to 
convert common SNP array formats to VCF format.

There are four currently supported array formats:

* Affymetrix (TSV export)
* Cytoscan HD Array (TSV export)
* Lumi 317k array (TSV export)
* Lumi 370k array (TSV export)
* Multi-sample OpenArray (TSV export)

Binary formats are not (yet) supported.
 

# Requirements

* Python 3.6
* requests

# CLI usage

The `array-as-vcf` tool will convert array files to VCF format.
It will auto-detect the type of array file, and throw an error if it can't
determine it. 

The generated VCF file is printed to stdout.

A sample name to be used in the VCF file _must_ be supplied.

The REF and ALT alleles will be queried from Ensembl if no `lookup-table` is
supplied. This requires a working internet connection, and can be quite slow
due the amount of HTTP requests that are necessary.

When supplied with `lookup-table`, no requests are made for the rsIDs 
which exist within the lookup table. The lookup table is a JSON file,
containing a single large object of shape:

```json
{
  "rs0": "{ref_allele}:{alt_alleles}:{ref_is_minor_allele}"
}
``` 

E.g. 

```json
{
  "rs1000003": "A:G:F"
}
```

If you have never run `array-as-vcf` before , you can run `array-as-vcf` sans lookup table
and `dump` the generated internal lookup table to a file for next iterations.

```bash
Usage: array-as-vcf [OPTIONS]

Options:
  -p, --path PATH              Path to array file  [required]
  -b, --build [GRCh37|GRCh38]
  -s, --sample-name TEXT       Name of sample in VCF file
  -c, --chr-prefix TEXT        Optional prefix to chromosome names
  -l, --lookup-table PATH      Optional path to existing lookup table for
                               rsIDs.
  -d, --dump PATH              Optional path to write generated lookup table
  --help                       Show this message and exit.

```
