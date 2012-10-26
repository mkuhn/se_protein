#!/usr/bin/env python2.7
# encoding: utf-8

from __future__ import print_function

import sys
import os
import re

from collections import defaultdict

def readActions():

    current_cid = None
    actions = defaultdict(list)

    for line in os.popen("zcat actions.tsv.gz"):
        (cid, protein, action, _, _, score) = line.strip("\n").split("\t")

        if int(score) < 500: continue

        if cid != current_cid:
            if actions:
                yield current_cid, actions
            current_cid = cid
            actions = defaultdict(list)

        protein = protein.split(".", 1)[1]

        actions[protein].append(action)

    yield current_cid, actions

def main():

    for (cid, protein_actions) in readActions():
        cid = int(cid[4:])

        binding_proteins = sum( "binding" in actions for actions in protein_actions.values() )

        if binding_proteins > 100:
            continue

        for protein, actions in protein_actions.items():
            if "binding" not in actions and "catalysis" not in actions: continue

            print(cid, protein, sep="\t")

            if "inhibition" in actions:
                print(cid, protein + "@inh", sep="\t")

            if "activation" in actions:
                print(cid, protein + "@act", sep="\t")



if __name__ == '__main__':
    main()
