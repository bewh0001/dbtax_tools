#!/usr/bin/env python3

"""
For an associated sample fastq file, extract positive read IDs from
classifier profiles and aggregate reads with the same predicted taxonomic ID.
The samples to include are specified in a samplesheet and the reads are
matched against an input file of expected taxa. A read will be considered
positive if its assigned taxon matches an expected taxon or a descended taxa
at a lower taxonomic rank. An output tsv file is generated with the format

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
    "--expected_taxa",
    "expected_taxa_fp",
    required=True,
    type=click.File("r"),
    help="input taxnoodle tsv file",
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
    required=False,
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

def main(
        samplesheet_fp: TextIO,
        expected_taxa_fp: TextIO,
        output_path: str,
        tool: str,
        taxonomy: str,
        summarise_at: str
):

    output_file = Path(output_path)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    taxa = taxdmp_tools.create_taxa(taxonomy=taxonomy)


    expected_taxids = parse_expected_taxa(expected_taxa_fp)
    samplesheet = parse_samplesheet(samplesheet_fp, classifier=tool)
    
    classifier_profiles = parse_profiles(samplesheet=samplesheet, classifier=tool)
    std_profiles = standardise_profiles(profiles=classifier_profiles)
    filtered_profiles = filter_profiles(
        profiles=std_profiles, expected_taxids=expected_taxids, taxa=taxa
    )

    output_data = get_output_data(
        profiles=filtered_profiles, samplesheet=samplesheet,
        summarise_at=summarise_at, taxa=taxa
    )
    format_output(output_data).to_csv(output_file, sep="\t", index=False)

def filter_profiles(profiles: dict, expected_taxids: list, taxa: dict):
    def is_taxid_expected(taxid: int):
        return(taxdmp_tools.ancestor_is_in(
            first_taxid=taxid, ancestors=expected_taxids, taxa=taxa
        ))
    filtered_profiles = {}
    for sample,profile in profiles.items():
        taxid_is_expected = list(map(is_taxid_expected, profile["taxid"]))
        filtered_profiles[sample] = profile.loc[taxid_is_expected, :].copy()
    return(filtered_profiles)

def parse_expected_taxa(expected_taxa_fp: TextIO):
    expected_taxa = expected_taxa_fp.readlines()
    taxids = [int(line.split("\t")[0]) for line in expected_taxa[1:]]
    return(taxids)

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

def get_output_data(
        profiles: dict, samplesheet: pandas.DataFrame, summarise_at: str, taxa:dict
):
    df = pandas.DataFrame({"sample": [], "fastq": [], "taxid": [], "reads": []})
    for sample, profile in profiles.items():
        for taxid in set(profile["taxid"]):
            reads = list(profile.loc[profile["taxid"] == taxid, "read_id"])
            fastq = list(samplesheet.loc[samplesheet["sample"] == sample, "fastq"])[0]
            summarised_taxid = taxdmp_tools.get_ancestor_at_rank(
                first_taxid=taxid, target_rank=summarise_at, taxa=taxa
            )
            row = pandas.DataFrame(
                {"sample": [sample], "fastq": [fastq],
                 "taxid": [summarised_taxid], "reads": [reads]}
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
        if not line.strip():
            continue
        fields = line.split(sep="|")
        data["read_id"] = fields[0].strip()
        data["rank"] = fields[1].strip()
        data["taxname"] = fields[2].strip()
        data["taxid"] = int(fields[3].strip())
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

if __name__ == "__main__":
    main()
