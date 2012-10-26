#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

import sys
import os

from collections import defaultdict

import re
import math

import itertools

import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def main():

    d = {
        "GO:0016301" : "kinase",
        "GO:0004879" : "nuclear receptor",
        "GO:0004930" : "GPCR",
        "GO:0005216" : "ion channel",
        # "GO:0004888"
        "GO:0003824" : "non-kinase enzyme"
    }
    
    annotations = defaultdict(lambda : defaultdict(int))
    direct_annotations = defaultdict(lambda : defaultdict(int))
    
    for line in os.popen("zcat /g/bork8/mkuhn/data/se_protein/go/human_geneontology.tsv.gz | egrep -v 'GO:0000155|GO:0000156|GO:0000160|GO:0007234|GO:0070297|GO:0070298|GO:0070299|GO:0008020'"):
        (species, protein, ontology_class, go_term, evidence, direct_go_term) = line.strip("\n").split("\t") 
        
        if go_term in d:
            annotations[protein][go_term] += 1

        if direct_go_term in d:
            direct_annotations[protein][direct_go_term] += 1
            

    fh_out = open("/home/mkuhn/data/se_protein/go_classification.tsv", "w")

    targets = set( line.strip() for line in os.popen("cut -f 2 /home/mkuhn/data/se_protein/interactions.tsv | sort -u | egrep -o 'ENSP[0-9]+'") )
    
    for (protein, go_terms) in annotations.iteritems():
        
        if "GO:0003824" in go_terms and go_terms.get("GO:0016301", 0) >= go_terms.get("GO:0003824", 0):
            go_terms.pop("GO:0003824")

        direct_go_terms = direct_annotations.get(protein, ())

        if "GO:0003824" in direct_go_terms and direct_go_terms.get("GO:0016301", 0) >= direct_go_terms.get("GO:0003824", 0):
            direct_go_terms.pop("GO:0003824")
        
        if direct_go_terms:
            
            # if protein not in ("ENSP00000320239",):
            if protein in targets and protein not in ("ENSP00000345004",):
                if len(direct_go_terms) > 1:
                    print >> sys.stderr, "D", protein, direct_go_terms
                # assert len(direct_go_terms) == 1
                
            for go_term in direct_go_terms:
                print >> fh_out, "\t".join( (protein, go_term, d[go_term]) )
        else:
            
            max_count = max( go_terms.values() )
            
            max_go_terms = [ k for (k,v) in go_terms.items() if v == max_count ]
                        
            if protein in targets:
                if len(max_go_terms) > 1:
                    print >> sys.stderr, protein, go_terms,
                    if "GO:0004930" in max_go_terms: max_go_terms = ["GO:0004930",]
                    elif "GO:0004930" in max_go_terms: max_go_terms = ["GO:0004930",]
                    print >> sys.stderr, "resolved to", max_go_terms[0]
                # assert len(max_go_terms) == 1
            
            for go_term in max_go_terms:
                print >> fh_out, "\t".join( (protein, go_term, d[go_term]) )

        
if __name__ == '__main__':
    main()
