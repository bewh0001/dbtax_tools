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
    "--outdir",
    default="./results/sylph",
    help="path to output directory for taxonomy.tsv file",
)
def build_sylph_taxonomy(
    samplesheet_fp: TextIO, nodes_dmp_fp: TextIO, names_dmp_fp: TextIO, outdir: str
):
    
    Path(outdir).mkdir(parents=True, exist_ok=True)
    os.chdir(outdir)
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

    with open("taxonomy.tsv", "w") as f_taxonomy:
        f_taxonomy.write("")

    with open("taxonomy.tsv", "a") as f_taxonomy:
        next(samplesheet_fp)  # skip header
        for line in samplesheet_fp:
            if not line.strip():
                continue
            fields = line.split(sep=",")
            taxid = int(fields[1])
            fasta_dna = format_sylph_fasta_filename(fields[2]).split("/")[-1]
            taxonomy = format_sylph_taxonomy(
                taxdmp_tools.get_lineage(taxid, taxa, wanted_ranks)
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
