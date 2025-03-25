#!/usr/bin/env python3

import click
import taxdmp_tools

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
   "--taxid",
   "taxid",
   required=True,
   type=click.INT,
   help="Taxonomic ID to get ancestor from"
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
def taxid2ancestor(taxid: int, target_rank: str, name: bool, rank: bool, taxonomy: str):
    taxa = taxdmp_tools.create_taxa(taxonomy=taxonomy)
    new_taxid = taxdmp_tools.get_ancestor_at_rank(
        first_taxid=taxid, target_rank=target_rank, taxa=taxa
    )
    fields = [str(new_taxid)]
    if name:
        fields.append(taxa[new_taxid]["name"])
    if rank:
        fields.append(taxa[new_taxid]["rank"])
    print("\t".join(fields))

if __name__ == "__main__":
    cli()
