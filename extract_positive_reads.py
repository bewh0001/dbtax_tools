#!/usr/bin/env python3

"""
For an associated sample fastq file, extract positive read IDs from
classifier profiles and expected positive taxa. Reads with the
same ancestor taxon at a given taxonomic rank are aggregated. Any
hit at a higher rank will be ignored. The samples to include are 
specified in a samplesheet and an output tsv file is generated with the format

sample  fastq   taxid   reads
sample1    /path/to/sample1.fq  taxid1  read1,read2,...,readN
...

where taxid1 is the taxonomic identifier of an expected taxon 
"""

import click
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
    help="input tsv samplesheet of sample profiles",
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
@click.option(
    "--taxonomy",
    "taxonomy",
    required=True,
    type=click.Path(file_okay=False),
    help="path to directory with NCBI taxdump files",
)
@click.option(
    "--summarise_at",
    "summarise_at",
    type=click.STRING,
    required=True,
    help="taxonomy level up to which lower level positive reads are aggregated"
)
@click.option(
    "--expected",
    "expected_fp",
    type=click.File("r"),
    help="text file with expected positive taxids, one per line"
)

def main(
        samplesheet_fp: TextIO,
        output_path: str,
        tool: str,
        taxonomy: str,
        summarise_at: str,
        expected_fp: TextIO
):

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    taxa = taxdmp_tools.create_taxa(taxonomy=taxonomy)

    samplesheet = parse_samplesheet(samplesheet_fp, classifier=tool)
    expected_taxa = parse_expected_taxa(expected_fp)
    
    classifier_profiles = parse_profiles(samplesheet=samplesheet, classifier=tool)
    std_profiles = standardise_profiles(profiles=classifier_profiles)
    summarised_profiles = summarise_profiles(profiles=std_profiles,
                                             summarise_at=summarise_at,
                                             taxa=taxa)
    filtered_profiles = filter_profiles(profiles=summarised_profiles, expected_taxa=expected_taxa)

    output_data = get_output_data(profiles=filtered_profiles, samplesheet=samplesheet, expected_taxa=expected_taxa)
    format_output(output_data).to_csv(output_file, sep="\t", index=False)

def parse_samplesheet(samplesheet_fp: TextIO, classifier: str):
    match classifier:
        case "kraken2":
            profile_idx = 2
        case "diamond":
            profile_idx = 3
        case "metabuli":
            profile_idx = 4
        case "metacache":
            profile_idx = 5
        case "sylph":
            profile_idx = 6
    data = {"sample": [], "fastq": [], "profile": []}
    next(samplesheet_fp)
    for line in samplesheet_fp:
        if not line.strip():
            continue
        fields = line.split("\t")
        data["sample"].append(fields[0])
        data["fastq"].append(fields[1].strip())
        data["profile"].append(fields[profile_idx].strip())
    return(pandas.DataFrame(data))

def get_output_data(profiles: dict, samplesheet: pandas.DataFrame, expected_taxa: tuple):
    df = pandas.DataFrame({"sample": [], "fastq": [], "taxid": [], "reads": []})
    for sample, profile in profiles.items():
        for taxid in expected_taxa:
            reads = list(profile.loc[profile["taxid"] == taxid, "read_id"])
            fastq = list(samplesheet.loc[samplesheet["sample"] == sample, "fastq"])[0]
            row = pandas.DataFrame(
                {"sample": [sample], "fastq": [fastq],
                 "taxid": [taxid], "reads": [reads]}
            )
            df = pandas.concat([df,row])
    return(df.astype({"taxid": int}))

def format_output(data: pandas.DataFrame):
    formatted_data = data.copy()
    formatted_data["reads"] = data["reads"].apply(lambda reads_list: ",".join(reads_list))
    return(formatted_data)

def standardise_profiles(profiles: dict):
    def standardise_profile(profile: pandas.DataFrame):
        return(profile.loc[profile["taxid"] > 0, ["taxid","read_id"]])
    std_profiles = {}
    for sample in profiles:
        std_profiles[sample] = standardise_profile(profile=profiles[sample])
    return(std_profiles)

def summarise_profiles(profiles: dict, summarise_at: str, taxa: dict):
    def summarise_taxid(taxid: int):
        return(taxdmp_tools.get_ancestor_at_rank(
            first_taxid=taxid, target_rank=summarise_at, taxa=taxa
        ))
    def summarise_profile(profile: pandas.DataFrame):
        summarised_profile = profile.copy()
        summarised_profile["taxid"] = list(map(summarise_taxid,
                                    profile["taxid"]))
        summarised_profile["tax_rank"] = []
        for taxid in summarised_profile["taxid"]:
            try:
                summarised_profile["tax_rank"].append(taxa[taxid]["rank"])
            except KeyError:
                summarised_profile["tax_rank"].append("unknown")
        return(summarised_profile.loc[
            summarised_profile["tax_rank"] == summarise_at, ["taxid", "read_id"]])
            
    summarised_profiles = {}
    for sample in profiles:
        summarised_profiles[sample] = summarise_profile(profile=profiles[sample])
    return(summarised_profiles)

def parse_expected_taxa(expected_fp: TextIO):
    expected_taxa = []
    for line in expected_fp:
        if not line.strip():
            continue
        expected_taxa.append(int(line.strip()))
    return(tuple(expected_taxa))

def filter_profiles(profiles: dict, expected_taxa: tuple):
    def filter_profile(profile: pandas.DataFrame):
        return(profile.loc[profile["taxid"].isin(expected_taxa), ["taxid", "read_id"]])
    filtered_profiles = {}
    for sample in profiles:
        filtered_profiles[sample] = filter_profile(profile=profiles[sample])
    return(filtered_profiles)

def parse_profiles(samplesheet: pandas.DataFrame, classifier: str):
    profile_data = {}
    for sample in samplesheet["sample"]:
        profile_path = list(samplesheet.loc[samplesheet["sample"] == sample, "profile"])[0]
        if not profile_path:
            profile_data[sample] = pandas.DataFrame(
                {"read_id": [], "taxid": []}
            )
            continue
        with open(profile_path, "r") as profile_fp:
            profile = profile_fp.readlines()
        match classifier:
            case "kraken2":
                profile_data[sample] = parse_k2_profile(profile=profile)
            case "metabuli":
                profile_data[sample] = parse_metabuli_profile(profile=profile)
            case "diamond":
                profile_data[sample] = parse_diamond_profile(profile=profile)
            case "metacache":
                profile_data[sample] = parse_metacache_profile(profile=profile)
            case "sylph":
                profile_data[sample] = parse_sylph_mapped_reads(profile=profile)
    return(profile_data)

def parse_diamond_profile(profile: list):
    columns = ("read_id", "taxid", "e-value")
    data = {col: [] for col in columns}
    for line in profile:
        if not line.strip():
            continue
        fields = line.split(sep="\t")
        data["read_id"].append(fields[0])
        data["taxid"].append(int(fields[1]))
        data["e-value"].append(float(fields[2]))
    return(pandas.DataFrame(data))

def parse_metabuli_profile(profile: list):
    columns = ("status", "read_id", "taxid", "read_len", "DNA_ident", "rank", "match_count")
    data = {col: [] for col in columns}
    for line in profile:
        if not line.strip():
            continue
        fields = line.split(sep="\t")
        data["status"].append(int(fields[0]))
        data["read_id"].append(fields[1])
        data["taxid"].append(int(fields[2]))
        data["read_len"].append(int(fields[3]))
        data["DNA_ident"].append(float(fields[4]))
        data["rank"].append(fields[5])
        if len(fields) == 7:
            data["match_count"].append(fields[6])
        else:
            data["match_count"].append(match_count)
    return(pandas.DataFrame(data))

def parse_metacache_profile(profile: list):
    columns = {"read_id", "rank", "taxname", "taxid"}
    data = {col: [] for col in columns}
    for line in profile:
        if not line.strip() or line[0] == "#":
            continue
        fields = line.split(sep="|")
        data["read_id"].append(fields[0].strip())
        data["rank"].append(fields[1].strip())
        data["taxname"].append(fields[2].strip())
        data["taxid"].append(int(fields[3].strip()))
    return(pandas.DataFrame(data))

def parse_k2_profile(profile: list):
    columns = ("status", "read_id", "taxid", "read_len", "LCA_mapping")
    data = {col: [] for col in columns}
    for line in profile:
        if not line.strip():
            continue
        fields = line.split(sep="\t")
        data["status"].append(fields[0])
        data["read_id"].append(fields[1])
        data["taxid"].append(int(fields[2]))
        data["read_len"].append(int(fields[3]))
        data["LCA_mapping"].append(fields[4])
    return(pandas.DataFrame(data))

def parse_sylph_mapped_reads(profile: list):
    columns = ("read_id", "taxid")
    data = {col: [] for col in columns}
    for line in profile:
        if not line.strip():
            continue
        fields = line.split(sep="\t")
        data["read_id"].append(fields[0])
        data["taxid"].append(int(fields[1]))
    return(pandas.DataFrame(data))

if __name__ == "__main__":
    main()
