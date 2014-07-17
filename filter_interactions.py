#!/usr/bin/env python2.7
# encoding: utf-8

from __future__ import print_function

import sys
import os
import re

from collections import defaultdict

protein_blacklist = list( line.strip() for line in open("protein.blacklist") ) 


def main():
    
    protein_names = {}

    threshold = 5

    if "cutoff2" in os.getcwd():
        threshold = 2

    for line in open("protein_names"):
        (protein_id, protein_name, desc) = line.strip("\n").split("\t")
        protein_names[protein_id] = (protein_name, desc)

    drugs = set( int(line.split()[0] ) for line in open("full_interactions.tsv") )
    drugs_to_remove = set( -int(line.split()[0] ) for line in open("hobohm_remove_list.txt") )

    drugs = drugs.difference( drugs_to_remove )

    proteins = defaultdict(list)

    for line in open("full_interactions.tsv"):
        (cid, protein) = line.strip("\n").split("\t")

        cid = int(cid)

        if cid not in drugs:
            continue

        if protein.split("@")[0] in protein_blacklist:
            continue

        proteins[protein].append(cid)

    cids2protein = defaultdict(list)

    fh_out = os.popen("sort > target.counts", "w")

    for protein, cids in proteins.items():

        count = len(cids)

        if count < threshold:
            continue

        cids2protein[ tuple(sorted(cids)) ].append(protein)

        if "@" not in protein: 
            name, desc = protein_names[protein]
            print(name, count, protein, desc, sep="\t", file=fh_out)


    for cids, proteins in cids2protein.items():
        protein = "_".join(proteins)
        for cid in cids:
            print(cid, protein, sep="\t")




if __name__ == '__main__':
    main()
