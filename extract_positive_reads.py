#!/usr/bin/env python3

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
    "--taxnoodle",
    "taxnoodle_fp",
    required=True,
    type=click.File("r"),
    help="input taxnoodle tsv file",
)
@click.option(
    "--taxid-map",
    "taxid_map_fp",
    required=False,
    type=click.File("r"),
    help="tsv file mapping lower rank taxids to summarised rank taxids"
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
        taxnoodle_fp: TextIO,
        output_path: str,
        tool: str,
        taxid_map_fp: TextIO = None,
):

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if taxid_map_fp:
        taxid_map = parse_taxid_map(taxid_map_fp)
        taxnoodle = parse_taxnoodle(taxnoodle_fp=taxnoodle_fp, taxid_map=taxid_map)
    else:
        taxnoodle = parse_taxnoodle(taxnoodle_fp)
    samplesheet = parse_samplesheet(samplesheet_fp)
    raw_profiles = parse_profiles(samplesheet=samplesheet, classifier=tool)
    std_profiles = standardise_profiles(profiles=raw_profiles, classifier=tool)
    filtered_profiles = filter_profiles(profiles=std_profiles, taxnoodle=taxnoodle)

    output = get_output(profiles=filtered_profiles, samplesheet=samplesheet)
    output.to_csv(output_file, sep="\t", index=False)

def parse_taxid_map(taxid_map_fp: TextIO):
    taxid_map = {}
    next(taxid_map_fp)
    for line in taxid_map_fp:
        if not line.strip():
            continue
        fields = line.split("\t")
        original_taxid = int(fields[0])
        higher_taxid = int(fields[1])
        taxid_map[higher_taxid] = original_taxid
    return(taxid_map)

def filter_profiles(profiles: dict, taxnoodle: pandas.DataFrame):
    filtered_profiles = {}
    for sample,profile in profiles.items():
        taxids = list(taxnoodle.loc[taxnoodle["sample"] == sample, "taxid"])
        filtered_profile = profile.loc[profile["taxid"].isin(taxids), :].copy()
        filtered_profiles[sample] = filtered_profile
    return(filtered_profiles)

def parse_taxnoodle(taxnoodle_fp: TextIO, taxid_map: dict = None):
    taxnoodle = taxnoodle_fp.readlines()
    data = {"taxid": [], "sample": []}
    
    samples = taxnoodle[0].strip().split("\t")[4:]
    for line in taxnoodle[1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        if taxid_map:
            taxid = taxid_map[int(fields[0])]
        else:
            taxid = int(fields[0])
        counts = list(map(int,fields[4:]))
        pos_samples = [samples[i] for i in range(len(samples)) if counts[i] > 0]
        for sample in pos_samples:
            data["taxid"].append(taxid)
            data["sample"].append(sample)
    return(pandas.DataFrame(data))

def parse_samplesheet(samplesheet_fp: TextIO):
    data = {"sample": [], "fastq": [], "profile": []}
    next(samplesheet_fp)
    for line in samplesheet_fp:
        if not line.strip():
            continue
        fields = line.split("\t")
        data["sample"].append(fields[0])
        data["fastq"].append(fields[1].strip())
        data["profile"].append(fields[2].strip())
    return(pandas.DataFrame(data))

def get_output(profiles: dict, samplesheet: pandas.DataFrame, taxid_map: dict = None):
    df = pandas.DataFrame({"sample": [], "fastq": [], "taxid": [], "reads": []})
    for sample, profile in profiles.items():
        for taxid in set(profile["taxid"]):
            reads = list(profile.loc[profile["taxid"] == taxid, "read_id"])
            fastq = list(samplesheet.loc[samplesheet["sample"] == sample, "fastq"])[0]
            row = pandas.DataFrame(
                {"sample": [sample], "fastq": [fastq],
                 "taxid": [taxid], "reads": [",".join(map(str,reads))]}
            )
            df = pandas.concat([df,row])
    return(df)

def standardise_k2_profile(profile: pandas.DataFrame):
    std_profile = profile.loc[profile["status"] == "C", ["taxid","read_id"]].copy()
    return(std_profile)

def standardise_profiles(profiles: dict, classifier: str):
    def standardise_profile(profile: pandas.DataFrame):
        match classifier:
            case "kraken2":
                return(standardise_k2_profile)
    std_profiles = {}
    for sample in profiles:
        std_profiles[sample] = standardise_profile(classifier)(profile=profiles[sample])
    return(std_profiles)

def parse_profiles(samplesheet: pandas.DataFrame, classifier: str):
    profile_data = {}
    for sample in samplesheet["sample"]:
        profile_path = list(samplesheet.loc[samplesheet["sample"] == sample, "profile"])[0]
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
    return(profile_data)

def parse_k2_profile(profile: list):
    columns = ("status", "read_id", "taxid", "read_len", "LCA_mapping")
    data = {col: [] for col in columns}
    for line in profile:
        if not line:
            continue
        fields = line.split(sep="\t")
        data["status"].append(fields[0])
        data["read_id"].append(fields[1])
        data["taxid"].append(int(fields[2]))
        data["read_len"].append(int(fields[3]))
        data["LCA_mapping"].append(fields[4])
    return(pandas.DataFrame(data))

if __name__ == "__main__":
    main()
