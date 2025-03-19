#!/usr/bin/env python3

import click
from collections import OrderedDict
from typing import TextIO
import taxdmp_tools
from pathlib import Path
import pandas
from functools import reduce


@click.command()
@click.option(
    "--samplesheet",
    "samplesheet_fp",
    required=True,
    type=click.File("r"),
    help="input tsv samplesheet"
)
@click.option(
    "--taxonomy",
    "taxonomy",
    required=True,
    type=click.Path(file_okay=False),
    help="path to directory with NCBI taxdump files"
)
@click.option(
    "--output",
    "output_path",
    required=True,
    type=click.Path(dir_okay=False),
    help="path to output file"
)
@click.option(
    "--tool",
    "tool",
    required=True,
    type=click.STRING,
    help="classifier tool used to generate sample profiles"
)
@click.option(
    "--summarise-at",
    "summarise_at",
    default="",
    type=click.STRING,
    help="summarise abundance profiles up to the given taxonomic rank and ignore abundances at higher ranks"
)

def main(
        samplesheet_fp: TextIO,
        taxonomy: str,
        output_path: str,
        tool: str,
        summarise_at: str
):

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    taxa = taxdmp_tools.create_taxa(taxonomy = taxonomy)
    profiles = get_sample_profiles(samplesheet_fp = samplesheet_fp)
    raw_profile_data = parse_profiles(profiles=profiles, classifier=tool)
    standardised_data = standardise_profiles(data=raw_profile_data,
                        classifier=tool, summarise_at=summarise_at, taxa=taxa)

    if summarise_at:
        taxid_map = taxdmp_tools.map_taxids_to_higher(taxids=standardised_data["taxonomy_id"],
                                         target_rank=summarise_at, taxa=taxa)
        summarised_data = summarise_data_at(
            data=standardised_data, taxa=taxa, taxid_map=taxid_map)
        wide_summarised_data = format_tax_data(summarised_data)
        taxid_map_file = Path(output_path + ".summarised_to_" + summarise_at + "_taxid_map")
        wide_summarised_data.to_csv(Path(output_path + ".summarised_to_" + summarise_at), sep="\t")
        taxid_map_to_df(taxid_map).to_csv(taxid_map_file, sep="\t", index=False)
        
    wide_data = format_tax_data(standardised_data)
    wide_data.to_csv(output_file, sep="\t")

def taxid_map_to_df(taxid_map: dict):
    df = pandas.DataFrame(
        {"from_taxid": list(taxid_map.keys()),
         "to_taxid": list(taxid_map.values())}
    )
    return(df)

def summarise_data_at(data: pandas.DataFrame, taxa: dict, taxid_map: dict):
    summarised_data = data[data.taxonomy_id.isin(list(taxid_map.keys()))].copy()
    summarised_data.loc[:,"taxonomy_id"] = [taxid_map[x] for x in summarised_data["taxonomy_id"]]
    summarised_data = summarised_data.groupby(
        ["taxonomy_id", "sample"], as_index=False).aggregate({"num_reads": "sum"})
    summarised_data["name"] = [taxa[taxid]["name"]
                               for taxid in summarised_data["taxonomy_id"]]
    summarised_data["rank"] = [taxa[taxid]["rank"]
                               for taxid in summarised_data["taxonomy_id"]]
    summarised_data["lineage"] = get_lineages_from_taxids(
        taxids=summarised_data["taxonomy_id"], taxa=taxa)
    return(summarised_data)

def standardise_profiles(data: dict, classifier: str, taxa: dict, summarise_at: str):
    std_data = {"sample": [], "taxonomy_id": [], "name": [], "rank": [],
                "num_reads": [], "lineage": []}
    for sample in data:
        taxid_counts = get_taxid_counts(classifier)(data=data[sample])
        std_data["sample"].extend([sample for i in range(len(taxid_counts))])
        std_data["taxonomy_id"].extend([taxid for taxid in taxid_counts.keys()])
        std_data["num_reads"].extend([count for count in taxid_counts.values()])    
        std_data["name"] = list(map(
            lambda taxid: taxdmp_tools.get_taxon_name(taxid, taxa),
            std_data["taxonomy_id"]))
        std_data["rank"] = list(map(
            lambda taxid: taxdmp_tools.get_taxon_rank(taxid, taxa),
            std_data["taxonomy_id"]))
        std_data["lineage"] = get_lineages_from_taxids(
            taxids=std_data["taxonomy_id"], taxa=taxa)
    return(pandas.DataFrame(std_data).astype({"taxonomy_id": int, "num_reads": int}))

def get_taxid_counts(classifier):
    match classifier:
        case "kraken2" | "metabuli" | "metacache":
            return(get_k2_counts)
        case "diamond":
            return(get_diamond_counts)
        case "sylph":
            return(get_sylph_counts)

def get_k2_counts(data: dict):
    taxid_counts = {
        taxid: count for taxid, count
        in zip(data["taxonomy_id"], data["taxon_reads"])
        if (count > 0 and taxid > 0)
    }
    return(taxid_counts)

def get_diamond_counts(data: dict):
    taxids = list({taxid for taxid in data["taxonomy_id"] if taxid > 0})
    taxid_counts = {}
    for taxid in taxids:
        num_reads = len([None for n in data["taxonomy_id"] if n == taxid])
        taxid_counts[taxid] = num_reads
    return(taxid_counts)

def get_sylph_counts(data: dict):
    taxid_counts = {
        taxid: count for taxid, count
        in zip(data["taxonomy_id"], data["num_reads"])
        if taxid > 0
    }
    return(taxid_counts)

def get_lineages_from_taxids(taxids: list[int], taxa: dict):
    wanted_ranks = ("superkingdom", "phylum","class","order","family","genus","species")
    lineages = []
    for taxid in taxids:
        if taxid == 0:
            return("")
        lineages.append(format_lineage(taxdmp_tools.get_lineage(
            first_taxid=taxid, taxa=taxa, wanted_ranks=wanted_ranks
        )))
    return(lineages)

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
    profile_data = {}
    for sample in profiles:
        with open(profiles[sample], "r") as profile_fp:
            profile = profile_fp.readlines()
        match classifier:
            case "kraken2" | "metabuli":
                profile_data[sample] = parse_k2_style_report(report=profile)
            case "sylph":
                profile_data[sample] = parse_sylph_profile(profile=profile)
            case "diamond":
                profile_data[sample] = parse_diamond_report(report=profile)
            case "metacache":
                profile_data[sample] = parse_metacache_report(report=profile)
    return(profile_data)

def format_tax_data(long_data: dict):
    wide_data = pandas.DataFrame(long_data).pivot_table(
        index=["taxonomy_id","name","rank","lineage"],
        columns="sample",
        values="num_reads",
        aggfunc="sum",
        fill_value=0)
    return(wide_data)

def parse_k2_style_report(report: list):
    columns = ("percent", "clade_reads", "taxon_reads", "taxonomy_lvl",
               "taxonomy_id", "name")
    data = {col: [] for col in columns}
    for line in report:
        fields = line.split(sep="\t")
        data["percent"].append(float(fields[0]))
        data["clade_reads"].append(int(fields[1]))
        data["taxon_reads"].append(int(fields[2]))
        data["taxonomy_lvl"].append(fields[3])
        data["taxonomy_id"].append(int(fields[4]))
        data["name"].append(fields[5])
    return(data)

def parse_metacache_report(report: list):
    columns = ("rank", "name", "taxonomy_id", "taxon_reads", "percent")
    data = {col: [] for col in columns}
    for line in report[2:]:
        fields = line.split(sep="|")
        data["rank"].append(fields[0].strip())
        data["name"].append(fields[1].strip())
        data["taxonomy_id"].append(int(fields[2]))
        data["taxon_reads"].append(int(fields[3]))
        data["percent"].append(fields[4].strip())
    return(data)

def parse_diamond_report(report: list):
    columns = ("query_id", "taxonomy_id", "e_value")
    data = {col: [] for col in columns}
    for line in report:
        fields = line.split()
        data["query_id"].append(fields[0])
        data["taxonomy_id"].append(int(fields[1]))
        data["e_value"].append(float(fields[2]))
    return(data)

def parse_sylph_profile(profile: list):
    data = {"clade_name": [], "rel_abundance": [], "seq_abundance": [],
             "taxonomy_id": [], "num_reads": []}
    # skip profile header
    for line in profile[2:]:
        fields = line.split()
        clade = tuple(fields[0].split(sep="|"))
        # only consider full clades down to species level (strain level is not handled)
        if clade[-1][0:3] != "s__":
            continue
        data["clade_name"].append(fields[0])
        data["rel_abundance"].append(float(fields[1]))
        data["seq_abundance"].append(float(fields[2]))
        data["num_reads"].append(int(fields[3]))
        data["taxonomy_id"].append(get_sylph_taxid(clade = list(clade)))
    return(data)

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
