#!/usr/bin/env python3

from collections import OrderedDict
from typing import TextIO


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

def get_lineage(first_taxid: int, taxa: dict, wanted_ranks: tuple[str]) -> OrderedDict:
    taxonomy = OrderedDict()
    for rank in wanted_ranks:
        taxonomy[rank] = None # OrderedDict preserves order of insertion
    taxid = first_taxid
    while True:
        if taxid == 1:
            break
        if taxa[taxid]["rank"] in wanted_ranks:
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