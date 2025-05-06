#!/usr/bin/env python3

from collections import OrderedDict
from typing import TextIO
from pathlib import Path


def create_taxa(taxonomy: str):
    nodes_dmp = Path(taxonomy + "/nodes.dmp")
    names_dmp = Path(taxonomy + "/names.dmp")
    with open(nodes_dmp, "r") as nodes_dmp_fp, open(names_dmp, "r") as names_dmp_fp:
        taxa = build_taxa_dict(nodes_dmp_fp, names_dmp_fp)
    return(taxa)

def build_taxa_dict(nodes_dmp_fp: TextIO, names_dmp_fp: TextIO) -> dict:
    taxa = {}
    for line in nodes_dmp_fp:
        if not line.strip():
            continue
        fields = line.split(sep="|")
        taxid = int(fields[0])
        taxa[taxid] = {"parent": int(fields[1]), "rank": fields[2].strip()}

    for line in names_dmp_fp:
        if not line.strip():
            continue
        fields = line.split(sep="|")
        name_class = fields[3].strip()
        if name_class == "scientific name":
            taxid = int(fields[0])
            taxa[taxid]["name"] = "_".join(
                fields[1].strip().split()
            )
    return(taxa)

def get_taxon_name(taxid: int, taxa: dict):
    if taxid > 0:
        return(taxa[taxid]["name"])
    elif taxid == 0:
        return("Unclassified")

def get_taxon_rank(taxid: int, taxa: dict):
    if taxid > 0:
        return(taxa[taxid]["rank"])
    elif taxid == 0:
        return("no rank")

def get_ancestor_at_rank(first_taxid: int, target_rank: str, taxa: dict):
    taxid = first_taxid
    while taxid > 1:
        try:
            if taxa[taxid]["rank"] == target_rank:
                # return ancestor with target rank
                return(taxid)
            taxid = taxa[taxid]["parent"]
        except KeyError:
            print("Unknown taxid. Unable to get ancestor.")
            break
    # return original taxid if no ancestor of target rank found
    return(first_taxid)

def ancestor_is_in(first_taxid: int, ancestors: list, taxa: dict):
    taxid = first_taxid
    while taxid > 1:
        if taxid in ancestors:
            return(True)
        taxid = taxa[taxid]["parent"]
    return(False)

def get_lineage(
        first_taxid: int,
        taxa: dict,
        wanted_ranks: tuple[str],
        output_taxids: bool = False
    ) -> OrderedDict:
    taxonomy = OrderedDict()
    for rank in wanted_ranks:
        taxonomy[rank] = None # OrderedDict preserves order of insertion
    taxid = first_taxid
    while True:
        if taxid == 1:
            break
        if taxa[taxid]["rank"] in wanted_ranks:
            if output_taxids:
                taxonomy[taxa[taxid]["rank"]] = taxid
            else:
                taxonomy[taxa[taxid]["rank"]] = taxa[taxid]["name"]
        taxid = taxa[taxid]["parent"]
    return taxonomy

def prune_lineage_empty_ranks(lineage: OrderedDict):
    return(
        OrderedDict(
            (rank, name)
            for rank, name in lineage.items()
            if not name == None
        )
    )

def map_taxids_to_higher(taxids: list, target_rank: str, taxa: dict):
    taxid_map = {}
    for taxid in taxids:
        new_taxid = taxid
        while new_taxid != 1:
            if taxa[new_taxid]["rank"] == target_rank:
                #taxid_map.setdefault(new_taxid, []).append(taxid)
                taxid_map[taxid] = new_taxid
                break
            new_taxid = taxa[new_taxid]["parent"]
    return(taxid_map)