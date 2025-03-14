#!/usr/bin/env python3

import click
from collections import OrderedDict
from typing import TextIO
import taxdmp_tools
from pathlib import Path
import pandas


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
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    profiles = get_sample_profiles(samplesheet_fp = samplesheet_fp)
    taxonomic_counts_data = parse_profiles(profiles = profiles, classifier = tool)
    taxonomic_counts_data["lineage"] = get_lineages_from_taxids(
        taxdump = taxonomy, taxids = taxonomic_counts_data["taxid"]
    )

    wide_data = format_tax_data(taxonomic_counts_data)
    wide_data.to_csv(output_file, sep="\t")

def get_lineages_from_taxids(taxids: list[int], taxdump: str):
    nodes_dmp = Path(taxdump + "/nodes.dmp")
    names_dmp = Path(taxdump + "/names.dmp")
    with open(nodes_dmp, "r") as nodes_dmp_fp, open(names_dmp, "r") as names_dmp_fp:
        taxa = taxdmp_tools.build_taxa_dict(nodes_dmp_fp, names_dmp_fp)

    wanted_ranks = ("superkingdom", "phylum","class","order","family","genus","species")

    def taxid_to_lineage(taxid: int):
        return(format_lineage(taxdmp_tools.get_lineage(
            first_taxid=taxid, taxa=taxa, wanted_ranks=wanted_ranks)))
    
    return(list(map(taxid_to_lineage, taxids)))

def get_sample_profiles(samplesheet_fp: TextIO):
    # skip samplesheet file header
    next(samplesheet_fp)

    profiles = {}
    for line in samplesheet_fp:
        if not line.strip():
            continue
        fields = line.split(sep="\t")
        sample = fields[0].strip()
        profile_path = Path(fields[1].strip())
        profiles[sample] = profile_path
    return(profiles)

def parse_profiles(profiles: dict, classifier: str):        
    counts_data = {"sample": [], "taxid": [], "count": []}
    for sample in profiles:
        with open(profiles[sample], "r") as profile_fp:
            profile = profile_fp.readlines()
        match classifier:
            case "sylph":
                taxid_counts = parse_sylph_profile(profile = profile)
            case _:
                taxid_counts = parse_profile_reads(profile = profile, classifier = classifier)
        for taxid, reads_count in taxid_counts.items():
            counts_data["sample"].append(sample)
            counts_data["taxid"].append(taxid)
            counts_data["count"].append(reads_count)
    return(counts_data)

def format_tax_data(long_data: dict):
    wide_data = pandas.DataFrame(long_data).pivot_table(
        index=["taxid","lineage"],
        columns="sample",
        values="count",
        fill_value=0)
    return(wide_data)

def parse_profile_reads(profile: list, classifier: str):
    def get_fields(line: str): return(line.split(sep="\t"))
    def is_classified(taxid): return(taxid > 0)
    match classifier:
        case "kraken2" | "metabuli":
            def get_taxid(fields: list): return(int(fields[2]))
        case "metacache":
            def get_fields(line: str): return(line.split(sep="|"))
            def get_taxid(fields: list): return(int(fields[3]))
        case "diamond":
            def get_taxid(fields: list): return(int(fields[1]))

    taxid_counts = {}
    for line in profile:
        # skip metacache header
        if classifier == "metacache" and line.strip()[0] == "#":
            continue
        fields = get_fields(line)
        taxid = get_taxid(fields)
        if not is_classified(taxid):
            continue
        taxid_counts[taxid] = taxid_counts.get(taxid, 0) + 1
    return(taxid_counts)

def parse_sylph_profile(profile: list):
    taxid_counts = {}
    # skip profile header
    for line in profile[2:]:
        fields = line.split(sep="\t")
        clade = tuple(fields[0].split(sep="|"))
        # only consider full clades down to species level (strain level is not handled)
        if clade[-1][0:3] != "s__":
            continue
        taxid = get_sylph_taxid(clade = list(clade))
        read_count = int(fields[3])
        taxid_counts[taxid] = taxid_counts.get(taxid, 0) + read_count
    return(taxid_counts)

def get_sylph_taxid(clade: list):
    while True:
        tax_string = clade.pop()
        try:
            # pick the lowest level valid taxid in clade
            taxid = int(tax_string[3:])
            if taxid > 0:
                break
        except ValueError:
            # if no valid taxid in clade, return 238384 (other sequences)
            if len(clade) == 0:
                return(28384)
    return(taxid)

def format_lineage(lineage: OrderedDict):
    pruned_lineage = taxdmp_tools.prune_lineage_empty_ranks(lineage)
    return(";".join(pruned_lineage.values()))

if __name__ == "__main__":
    main()
