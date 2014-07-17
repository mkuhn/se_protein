#!/usr/bin/env python2.7
# encoding: utf-8

from __future__ import print_function

import sys
import os
import re

from collections import defaultdict

q_threshold = 0.01

def readSE():

    metabolizing_proteins = set( s.split()[0] for s in open("metabolizing_proteins") )

    fh_in = open("protein_se_pv.tsv")
    fh_in.next()

    current_se = None
    best_protein = None
    best_q = q_threshold

    for line in fh_in:
        (se, protein, p, q) = line.strip("\n").split("\t")

        if se != current_se:

            if best_protein:
                yield current_se, best_protein, best_q

            current_se = se
            best_protein = None
            best_q = q_threshold

        q = float(q)

        if q <= best_q:

            for _protein in re.findall(r"ENSP\d+", protein):
                if _protein not in metabolizing_proteins: 
                    best_q = q
                    best_protein = protein
                    break

    yield current_se, best_protein, best_q


def main():

    go_classification = {}

    for line in open("go_classification.tsv"):
        (protein, go_id, go_name) = line.strip("\n").split("\t")

        # by accident, sorting alphabetically is a good priority list
        if protein not in go_classification or go_name < go_classification[protein][1]:
            go_classification[protein] = (go_id, go_name)

    go_classification["ENSP00000231509"] = ("GO:0004879", "nuclear receptor")

    for se, protein, q in readSE():

        go_id, go_name = "?", "?"

        for _protein in re.findall(r"ENSP\d+", protein):
            if _protein in go_classification:
                go_id, go_name = go_classification[_protein]
                break                

        print(se, protein, q, go_id, go_name, sep="\t")

if __name__ == '__main__':
    main()
