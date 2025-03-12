#!/usr/bin/env python3

import click
import os
import re
from collections import OrderedDict
from typing import TextIO
import taxdmp_tools
from pathlib import Path


@click.command()
@click.option(
    "--samplesheet",
    "samplesheet_fp",
    required=True,
    type=click.File("r"),
    help="samplesheet file in csv format",
)
@click.option(
    "--nodes_dmp",
    "nodes_dmp_fp",
    required=True,
    type=click.File("r"),
    help="path to nodes.dmp file",
)
@click.option(
    "--names_dmp",
    "names_dmp_fp",
    required=True,
    type=click.File("r"),
    help="path to names.dmp file",
)
@click.option(
    "--output",
    "output_path",
    required=True,
    type=click.Path(dir_okay=False),
    help="path to output file",
)
@click.option(
    "--use-taxids",
    "use_taxids",
    is_flag=True,
    default=False,
    help="generate a lineage of taxonomic IDs instead of names"
)

def build_sylph_taxonomy(
    samplesheet_fp: TextIO,
    nodes_dmp_fp: TextIO,
    names_dmp_fp: TextIO,
    output_path: str,
    use_taxids: bool
):
    
    output_file = Path(output_path)
    output_dir = output_file.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(output_dir)

    taxa = taxdmp_tools.build_taxa_dict(nodes_dmp_fp, names_dmp_fp)

    wanted_ranks = (
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    )

    with open(output_file, "w") as f_taxonomy:
        next(samplesheet_fp)  # skip header
        for line in samplesheet_fp:
            if not line.strip():
                continue
            fields = line.split(sep=",")
            taxid = int(fields[1])
            fasta_dna = format_sylph_fasta_filename(fields[2]).split("/")[-1]
            taxonomy = format_sylph_taxonomy(
                taxdmp_tools.get_lineage(
                    first_taxid=taxid,
                    taxa=taxa,
                    wanted_ranks=wanted_ranks,
                    use_taxids=use_taxids)
            )
            f_taxonomy.write("\t".join([fasta_dna, taxonomy]) + "\n")


def format_sylph_taxonomy(taxonomy: OrderedDict) -> str:
    # Use the preserved order of insertion to substitute in taxonomic ranks
    sylph_str = "d__{};p__{};c__{};o__{};f__{};g__{};s__{}".format(
        *list(taxonomy.values())
    )
    return sylph_str


def format_sylph_fasta_filename(fname: str) -> str:
    r = re.search("(^.*?)_(?:ASM|genomic).*", fname)
    if not r:
        return fname
    return r.groups()[0]

if __name__ == "__main__":
    build_sylph_taxonomy()
