#!/usr/bin/env python2.7
# encoding: utf-8

from __future__ import print_function

import sys
import os
import re

from collections import defaultdict


def main():
    
    valid_drugs = set( int(line ) for line in open("valid_drugs") )

    ses = defaultdict(set)

    se2name = {}

    threshold = 5

    if "cutoff2" in os.getcwd():
        threshold = 2

    for line in sys.stdin:
        (cid, _, concept, drug, name) = line.strip("\n").split("\t")

        cid = int(cid[2:])

        if cid not in valid_drugs:
            continue

        ses[concept].add( (cid, drug) )

        se2name.setdefault(concept, name)

    cids2concepts = defaultdict(list)

    for concept, drugs in ses.items():

        if len(drugs) < threshold:
            continue

        cids2concepts[ tuple(sorted(drugs)) ].append(concept)

    # merge side effects that occur for the same drug, take the one with the longest name
    for drugs, concepts in cids2concepts.items():

        concept = None
        name = ""

        for _concept in concepts:
            _name = se2name[_concept]
            if len(_name) > len(name):
                concept = _concept
                name = _name

        for (cid, drug) in drugs:
            print(cid, concept, drug, name, sep="\t")


if __name__ == '__main__':
    main()
