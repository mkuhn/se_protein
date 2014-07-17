#!/usr/bin/env python2.7
# encoding: utf-8

from __future__ import print_function

import sys
import os
import re

from collections import defaultdict


def main():
    
    known_proteins = set()

    known_drugs = set()

    protein_map = {}

    for line in open("../interactions.tsv"):
        (drug, proteins) = line.strip("\n").split("\t")
        known_proteins.update( re.findall("ENSP\d+", proteins) )
        known_drugs.add(int(drug))

        for protein in proteins.split("_"):
            protein_map[ protein ] = proteins

    chemical_sources = defaultdict(list)

    for line in open("filtered_actann_interact.tsv"):
        (action, _, protein, chemical, source) = line.strip("\n").split("\t")

        if protein not in known_proteins:
            continue

        if int(chemical[3:]) not in known_drugs:
            continue

        chemical_sources[chemical, source].append( (protein, action) )

    for (chemical, source), proteins in chemical_sources.items():

        if len(proteins) >= 10:
            continue

        for protein, action in proteins:

            act_protein = protein + "@" + action

            if act_protein in protein_map:
                print(protein_map[ act_protein ], chemical, sep="\t")

            print(protein_map[ protein ], chemical, sep="\t")



if __name__ == '__main__':
    main()
