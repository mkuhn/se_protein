#!/usr/bin/env python2.7
# encoding: utf-8

from __future__ import print_function

import sys
import os
import re

from collections import defaultdict


def main():
    
    mapping = {}

    for line in open("string7to9.tsv"):
        (s7, s9) = line.strip("\n").split("\t")
        mapping[s7] = s9

    for line in sys.stdin:

        if line.startswith("#") or "ENSP" not in line: continue

        fields = line.strip("\n").split("\t")

        if re.search(r"(ENSP|D)[0-9]{5}", fields[2]):
            fields = fields[:2] + fields[3:]

        for protein in fields[1].split("_"):
            if protein.startswith("D"): continue

            if "@" in protein:
                protein, suffix = protein.split("@")
                suffix = "@" + suffix
            else:
                suffix = ""                

            if protein not in mapping:
                print("Not found:", protein, "in\n", line, file=sys.stderr)
                continue

            fields[1] = mapping[protein] + suffix

            print("\t".join(fields))


if __name__ == '__main__':
    main()
