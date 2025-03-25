#!/usr/bin/env python3

import click
import taxdmp_tools
from typing import TextIO
import sys

option_taxonomy = click.option(
   "--taxonomy",
   "taxonomy",
   required=True,
   type=click.Path(file_okay=False),
   help="Path to directory containing NCBI taxdump files"
)

@click.group(help="CLI tool to extract information from NCBI taxa")
def cli():
    pass

@cli.command(help="Get ancestor taxon at a given taxonomic rank")
@option_taxonomy
@click.option(
   "--input",
   "input_fp",
   type=click.File('r'),
   default=sys.stdin,
   help="File with Taxonomic IDs to get ancestors from"
)
@click.option(
    "--name/--no-name",
    default=False,
    help="Print ancestor taxon scientific name"
)
@click.option(
   "--target-rank",
   "target_rank",
   required=True,
   type=click.STRING,
   help="Name of ancestor taxonomic rank"
)
@click.option(
   "--rank/--no-rank",
   default=False,
   help="Print ancestor taxon rank"
)
def taxid2ancestor(input_fp: str, target_rank: str, name: bool, rank: bool, taxonomy: str):
    taxids = [line.strip() for line in input_fp]
    taxa = taxdmp_tools.create_taxa(taxonomy=taxonomy)

    for taxid in taxids:
        try:
            ancestor_taxid = taxdmp_tools.get_ancestor_at_rank(
                first_taxid=int(taxid), target_rank=target_rank, taxa=taxa
            )
            print(format_taxon_output(taxid=ancestor_taxid, name=name, rank=rank, taxa=taxa))
        except ValueError:
            print(taxid)

def format_taxon_output(taxid: int, name: bool, rank: bool, taxa: dict):
    fields = [str(taxid)]
    if name:
        fields.append(taxa[taxid]["name"])
    if rank:
        fields.append(taxa[taxid]["rank"])
    return("\t".join(fields))

def parse_taxids(input_fp: TextIO):
    taxids = []
    for line in input_fp:
        if not line.strip():
            continue
        taxids.append(line.strip())
    return(taxids)

if __name__ == "__main__":
    cli()
