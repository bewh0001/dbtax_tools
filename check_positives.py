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
    "--input",
    "input_fp",
    required=True,
    type=click.File("r"),
    help="input tsv file with taxonomic hit per sample read",
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
    "--neg_ctrls",
    "neg_ctrls",
    required=True,
    type=click.STRING,
    help="comma-separated list of negative control sample identifiers",
)
@click.option(
    "--output",
    "output_path",
    required=True,
    type=click.Path(dir_okay=False),
    help="path to output file",
)

def main(
        input_fp: TextIO,
        nodes_dmp_fp: TextIO,
        names_dmp_fp: TextIO,
        neg_ctrls: str,
        output_path: str
):
    
    output_file = Path(output_path)
    output_directory = output_file.parent
    output_directory.mkdir(parents=True, exist_ok=True)

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

    ctrl_sample_ids = neg_ctrls.split(sep=',')
    
    sample_taxids = {}
    ctrl_taxids = set()

    # skip input file header
    next(input_fp)
    
    for line in input_fp:
        if not line.strip():
            continue
        sample, taxid = line.split(sep='\t')
        taxid = int(taxid)

        if sample in ctrl_sample_ids:
            ctrl_taxids.add(taxid)
        else:
            sample_taxids.setdefault(sample, set()).add(taxid)

    ctrl_lineages = []
    for taxid in ctrl_taxids:
        ctrl_lineages.append(taxdmp_tools.get_lineage(
            first_taxid = taxid,
            taxa = taxa,
            wanted_ranks = wanted_ranks
        ))

    ctrl_taxa = extract_lineage_ranks(ctrl_lineages)

    checked_lineages = {}

    for sample in sample_taxids:
        lineages = list(map(
            lambda taxid: taxdmp_tools.get_lineage(
                first_taxid = taxid,
                taxa = taxa,
                wanted_ranks = wanted_ranks
            ), sample_taxids[sample]
        ))

        checked_lineages[sample] = check_lineages(
            lineages = lineages,
            ctrl_taxa = ctrl_taxa 
        )
    
    with open(Path(output_path), "a") as f_output:
        header = "\t".join(["sample"] + list(wanted_ranks))
        f_output.write(header + "\n")
        for sample in checked_lineages:
            for lineage in checked_lineages[sample]:
                f_output.write("\t".join([sample] + list(lineage.values())) + "\n")

def extract_lineage_ranks(lineages: tuple[OrderedDict]):
    seen_taxa = {}
    for lineage in lineages:
        for rank, taxon in lineage.items():
            if not taxon:
                continue
            seen_taxa.setdefault(rank, set()).add(taxon)
    return(seen_taxa)

def check_lineages(
        lineages: tuple[OrderedDict],
        ctrl_taxa: dict
):
    def lineage_in_ctrl(lineage: OrderedDict):
        lowest_rank = list(lineage.keys())[-1]
        return(lineage[lowest_rank] in ctrl_taxa[lowest_rank])

    new_lineages = []
    for lineage in lineages:
        pruned_lineage = taxdmp_tools.prune_lineage_empty_ranks(lineage)
        if len(pruned_lineage.keys()) == 0:
            continue
        elif not lineage_in_ctrl(lineage = pruned_lineage):
            new_lineages.append(pruned_lineage)
    return(tuple(new_lineages))

if __name__ == "__main__":
    main()
