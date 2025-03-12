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
    help="input tsv samplesheet",
)
@click.option(
    "--taxonomy",
    "taxonomy",
    required=True,
    type=click.Path(file_okay=False),
    help="path to directory with NCBI taxdump files",
)
@click.option(
    "--output",
    "output_path",
    required=True,
    type=click.Path(dir_okay=False),
    help="path to output file",
)
@click.option(
    "--tool",
    "tool",
    required=True,
    type=click.STRING,
    help="classifier tool used to generate sample profiles"
)

def main(
        samplesheet_fp: TextIO,
        taxonomy: str,
        output_path: str,
        tool: str
):
    
    output_file = Path(output_path)
    output_directory = output_file.parent
    output_directory.mkdir(parents=True, exist_ok=True)
    nodes_dmp = Path(taxonomy + "/nodes.dmp")
    names_dmp = Path(taxonomy + "/names.dmp")

    with open(nodes_dmp, "r") as nodes_dmp_fp, open(names_dmp, "r") as names_dmp_fp:
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
    
    # skip samplesheet file header
    next(samplesheet_fp)

    samples = OrderedDict()

    for line in samplesheet_fp:
        if not line.strip():
            continue
        fields = line.split(sep="\t")
        sample = fields[0].strip()
        profile_path = Path(fields[1].strip())
        samples[sample] = {"profile": profile_path}
    
    taxid_counts = {}

    for sample in samples:
        with open(samples[sample]["profile"], "r") as profile_fp:
            def get_taxid(tool: str):
                match tool:
                    case "kraken2":
                        return(get_kraken2_taxid)
                    case "diamond":
                        return(get_diamond_taxid)
                    case "metabuli":
                        return(get_metabuli_taxid)
                    case "metacache":
                        return(get_metacache_taxid)
                    case "sylph":
                        return(get_sylph_taxid)

            for line in profile_fp:
                if not line.strip():
                    continue
                taxid = get_taxid(tool = tool)(line = line)
                if taxid > 1:
                    # initialize all sample read counts to 0
                    if taxid not in taxid_counts:
                        taxid_counts[taxid] = {s: 0 for s in samples}
                    taxid_counts[taxid][sample] += 1

    with open(output_file, "w") as output_fp:
        header = "\t".join(["taxonomy_id", "lineage"] + list(samples.keys()))
        output_fp.write(header + "\n")
        for taxid in taxid_counts:
            lineage = taxdmp_tools.get_lineage(
                first_taxid = taxid,
                taxa = taxa,
                wanted_ranks = wanted_ranks
            )
            line = "\t".join(
                [str(taxid), format_lineage(lineage)]
                + list(map(str,taxid_counts[taxid].values()))
            )
            output_fp.write(line + "\n")

def format_lineage(lineage: OrderedDict):
    pruned_lineage = taxdmp_tools.prune_lineage_empty_ranks(lineage)
    return(";".join(pruned_lineage.values()))

def get_kraken2_taxid(line: str):
    fields = line.split(sep="\t")
    classified = fields[0] == "C"
    if classified:
        taxid = int(fields[2])
        return(taxid)
    else:
        return(0)

def get_diamond_taxid(line: str):
    fields = line.split(setp="\t")
    taxid = int(fields[1])
    return(taxid)

def get_metabuli_taxid(line: str):
    fields = line.split(sep="\t")
    classified = fields[0] == "1"
    if classified:
        taxid = int(fields[2])
        return(taxid)
    else:
        return(0)

def get_metacache_taxid(line: str):
    # skip comments
    if line.strip()[0] == "#":
        return(0)
    fields = line.split(sep="|")
    taxid = int(fields[3].strip())
    return(taxid)

def get_sylph_taxid(line: str):
    fields = line.split(sep="\t")
    # skip first two lines
    if fields[0] in ("#SampleID", "clade_name"):
        return(0)
    
    clade = fields[0].split(sep="|")
    tax_string = clade.pop() # get last level in clade

    # only consider full clades down to species level (strain level is ignored)
    if tax_string[0:3] != "s__":
        return(0)
    while True:
        try:
            # pick the lowest level valid taxid in clade
            taxid = int(tax_string[3:])
            break
        except ValueError:
            # if no valid taxid in clade, return 28384 (other sequences)
            if len(clade) == 0:
                return(28384)
            tax_string = clade.pop()
    return(taxid)

if __name__ == "__main__":
    main()
