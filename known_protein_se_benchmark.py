#!/usr/bin/env python2.5
# encoding: utf-8
"""
untitled.py

Created by Michael Kuhn on 2008-04-07.
Copyright (c) 2008 __MyCompanyName__. All rights reserved.
"""

from __future__ import division

import sys
import os
import itertools
import re
import glob

from collections import defaultdict

protein_blacklist = list( line.strip() for line in open("protein.blacklist") ) 

priority_list = (
    "no"                ,
    "unknown"                ,
    "wrong"                ,
    "hypersensitivity_reaction"                ,
    "false_positive"    ,
    "anecdotal_evidence"    ,

    "opposite_phenotype"          ,
    "yes"    ,
    "best"    ,
    "literature_support",
    "related_phenotype" ,
    "same_phenotype"    ,
)

priorities = {}

for p in priority_list:
    priorities[p] = len(priorities)

over_under = "over"
# over_under = "under"

valid_se = set( line.strip() for line in open("valid_side_effects") )
valid_targets = set( line.strip() for line in open("valid_targets") )
for t in tuple(valid_targets): valid_targets.update(t.split("_"))



def cmpString(a,b):
    """take the string with fewer special chars, or the longer one
    >>> cmpString("(a)","a")
    -1
    >>> cmpString("(a)",",a")
    -1
    >>> cmpString("(a)",",b,")
    0
    >>> cmpString("abc","ab")
    1
    >>> sorted(("(a)","a"), cmp = cmpString)[-1]
    'a'
    >>> sorted(("abc","ab"), cmp = cmpString)[-1]
    'abc'
    """
    def c(s):
        return len(re.findall(r"[,()]", s))
        
    return cmp( c(b), c(a) ) or cmp( len(a), len(b) ) 
    

def print_pairs(se_protein_pairs, prefix, target2drugs, se2drugs, p2pn, s2sn, print_full):
    
    filename = "pred."+prefix+".se_protein_known.all.txt"
    
    if print_full:
        for line in open("protein_names"):
            (p, pn, _) = line.strip("\n").split("\t")
            p2pn.setdefault(p, pn)

        for line in open("se_names"):
            (s, sn) = line.strip("\n").split("\t")
            s2sn.setdefault(s, sn)
    
    known_action = "known_action" in prefix 
    mouse = "mouse" in prefix or "actinh" in prefix
    combined = prefix == "combined"
    actinh = "actinh" in prefix
    
    sp_pvalues = defaultdict(lambda : defaultdict(dict))
    fh_in = open("protein_se_pv.tsv")
    fh_in.next()
    for line in fh_in:
        (se, protein, p, q) = line.strip("\n").split("\t")

        if not print_full and se not in s2sn:
            continue

        for _protein in protein.split("_"):
            if _protein.split("@")[0] not in p2pn: continue
            sp_pvalues[se][_protein] = q 

    fh_complete = open(prefix+".se_protein_known.complete.txt", "w")
    fh_best = open(prefix+".se_protein_known.best.txt", "w")
    fh_all = open("plot."+prefix+".se_protein_known.all.txt", "w")
    fh_pred = os.popen("sort -t '\t' -k 5,6 > " + filename, "w")

    sp_known = {}

    # go through annotation of se-protein pairs
    for fields in se_protein_pairs:
        
        s,p,action,kind = fields[:4]
        reason = fields[4] if len(fields) == 5 else ""
    
        assert "_" not in p

        _p = p
    
        drugs_with_se = set(se2drugs.get(s, ()))

        results = []

        for suffix in ("", "@act", "@inh"):

            # print "\t".join( map(str, (t2t.get(p,"-"), p, s, suffix, pn, sn)) )
            # continue
            
            p = _p + suffix
            if (s,p) in sp_known: continue

            drugs_with_protein = set(target2drugs.get(p, ()))

            if not drugs_with_se or not drugs_with_protein:
                continue

            drugs_with_protein_and_se = drugs_with_protein.intersection( drugs_with_se )

            if drugs_with_protein_and_se:
                n_protein_with_se = "%.3g" % (100 * len(drugs_with_protein_and_se) / len(drugs_with_protein))
                n_se_with_protein = "%.3g" % (100 * len(drugs_with_protein_and_se) / len(drugs_with_se))
            else:
                n_protein_with_se = n_se_with_protein = "-"

            pvalue = sp_pvalues[s].get(p,"?")
            output = "\t".join( map(str, (kind, p, s, suffix, p2pn.get(p.split("@")[0], ""), s2sn.get(s, ""), len(drugs_with_protein), len(drugs_with_se), n_protein_with_se, n_se_with_protein, pvalue)) )

            if pvalue != "?":
                results.append( (pvalue, p, s, output) ) 
                print >> fh_complete, output

        
        if results:
            
            if known_action:
                
                for (_, p, s, _) in results:
                    
                    if "@" not in p:
                        known = "unspecified_action"
                    elif p.endswith(action):
                        known = "same_action"
                    else:
                        known = "opposite_action"
                    
                    sp_known[s,p] = known

            elif combined:
                
                for (_, p, s, _) in results:
                    
                    if kind == "not relevant":
                        known = "not_relevant"
                    elif "@" in p and p.endswith(action):
                        known = kind
                    else:
                        known = "other"

                    sp_known[s,p] = known, reason
                    
            elif mouse:        
            
                for (_, p, s, _) in results:
                    
                    sp_known[s,p] = kind, reason
            
            else:
                results.sort()
                print >> fh_best, results[0][-1]
            
                known = "best"
                for (_, p, s, _) in results:
                    sp_known[s,p] = known, reason
                    known = "yes"
        
    
    print >> fh_all, "\t".join(("known", "pvalue"))

    seen = set()

    for s, p_pvalues in sp_pvalues.items():

        
        for p, pvalue in p_pvalues.items():

            if (s,p) in seen: continue

            seen.add((s,p))

            known = "no"
            reason = ""


            for _p in p.split("_"):
                if (s, p) not in sp_known: continue

                _known, _reason = sp_known[s,p]
                if priorities[_known] > priorities[known]:
                    known, reason = _known, _reason

            print >> fh_all, "\t".join( map(str, (known, pvalue) ) ) 

            pn = p2pn.get(p.split("@")[0], "") 
            if "@" in p:
                pn = pn + "@" + p.split("@")[1]
            
            print >> fh_pred, "\t".join( map(str, (known, pvalue, s2sn.get(s, ""), pn, s, p, reason ) ) )

    return filename
    
    
    
def load_whitebread_data(print_full):

    # protein name to protein id
    pn2p = {}
    p2pn = {}
    for line in open("/home/mkuhn/src/se_protein/side_effect_protein_mapping.txt"):
        fields = line.strip().split("\t")
        
        pn = fields[0]
        p = fields[-1]
        # pn2p[pn] = t2t.get(p,p)
        for _p in p.split():
            p2pn[p] = pn
        pn2p[pn] = p

    # side effect name to side effect id
    sn2s = {}
    s2sn = {}
    for line in open("/home/mkuhn/src/se_protein/side_effect_costart_annotation.txt"):
        fields = line.strip().split("\t")

        sn = fields[0].lower()
        s = fields[1]
        sn2s[sn] = s
        s2sn[s] = sn


    se_protein_pairs = set()

    for line in open("/home/mkuhn/src/se_protein/side_effect_protein_annotation.txt"):
        (pn, _, _, sn) = line.strip().split("\t")
    
        p = pn2p[pn]
        s = sn2s[sn.lower()]

        if s not in valid_se: continue

        if p == "?" or s == "?": continue 

        for _p in p.split():
            se_protein_pairs.add( (s,_p,"","","Whitebread et al.") )

    se_protein_pairs = list(se_protein_pairs)

    print >> sys.stderr, "Loaded", len(se_protein_pairs), "se-protein-pairs for Whitebread data"

    return se_protein_pairs, p2pn, s2sn


def load_annotated_data(rx):
    
    r = re.compile(rx).search
    
    s2sn = {}
    p2pn = {}
    se_protein_pairs = set()
    
    for line in open("/home/mkuhn/src/se_protein/se_annotation.tsv"):
        (s, p, sn, pn, ann, pmid) = line.strip("\n").split("\t")[:6]

        if s not in valid_se: continue
        if not r(ann): continue
        
        s2sn[s] = sn
        p2pn[p] = pn
        
        if ann.endswith("activation"):  
            action = "act"
        elif ann.endswith("inhibition"): 
            action = "inh"
        else:
            action = ""
        
        kind = ann.split("_")[0].strip()
        
        reason = "annotated from http://www.ncbi.nlm.nih.gov/pubmed/%d" % int(pmid[4:])
        
        se_protein_pairs.add( (s,p,action,kind,reason) )

    se_protein_pairs = list(se_protein_pairs)

    print >> sys.stderr, "Loaded", len(se_protein_pairs), "se-protein-pairs for regex", rx
    
    
    return se_protein_pairs, p2pn, s2sn



    
def load_mouse_data(print_full):
    
    s2sn = {}
    p2pn = {}
    se_protein_pairs = set()
    seen = set()

    action = "inh"
    kind = "same_phenotype"
    
    for line in open("mouse_se.tsv"):
        (s, p, sn, pn) = line.strip("\n").split("\t")[:4]

        if s not in valid_se: continue
        
        se_protein_pairs.add( (s,p,action,kind,"") )
        seen.add((s,p))

        s2sn[s] = sn
        p2pn[p] = pn
            
    print >> sys.stderr, "Direct phenotypes:", len(se_protein_pairs)

    equivalent_se = defaultdict(list)
    
    for line in open("/home/mkuhn/src/se_protein/equivalent_se.tsv"):
        (a, b) = line.strip("\n").split("\t")
        equivalent_se[b].append(a)
        equivalent_se[a].append(b)

    kind = "related_phenotype"

    # for gathering related phenotypes, go through the whole unfiltered set of mouse SE
    for line in open("/home/mkuhn/data/mgi/mouse_se.tsv"):
        (ms, p, sn, pn, mp) = line.strip("\n").split("\t")[:5]

        l = []
        l.extend( equivalent_se.get(ms, ()))
        l.extend( equivalent_se.get(mp, ()))
        
        for s in l:
            if (s,p) in seen: continue
            if p not in valid_targets: continue
            if not print_full:
                if s not in s2sn: continue
                if p.split("@")[0] not in p2pn: continue

            reason = sn
            
            se_protein_pairs.add( (s,p,action,kind, reason) )
            seen.add((s,p))
    
    relevant_se = set( s for (s,p) in seen )

    for line in open("src/ann/ann"):
        if not line.strip(): continue
    
        fields = line.strip().split("\t")
        
        assert not re.match("(D|ENSP)\d+", fields[2])
            
        (s, pp) = fields[:2]
        reason = " ".join(fields[2:])
    
        assert "_" not in pp, line
    
        kind = "literature_support"
        
        if s.startswith("#"):
            continue
        elif s.startswith("!"):
            s = s[1:]
            kind = "opposite_phenotype"
        elif s.startswith("~"):
            s = s[1:]
            kind = "false_positive"
        elif s.startswith("?"):
            s = s[1:]
            kind = "anecdotal_evidence"
        elif s.startswith("*"):
            s = s[1:]
            kind = "hypersensitivity_reaction"
        elif s.startswith("."):
            s = s[1:]
            kind = "unknown"
        elif s.startswith("-"):
            s = s[1:]
            kind = "wrong"
            
        if not print_full and s not in relevant_se: continue
        
        proteins = pp.split(" ")

        _proteins = set( protein.split("@")[0] for protein in proteins if (protein.startswith("ENSP") and protein in valid_targets) )

        for p in proteins:   

            if (s,p) in seen: continue

            if p not in valid_targets: continue
        
            if not print_full:
                if s not in s2sn: continue
                if p.split("@")[0] not in p2pn: continue
    

            if kind == "literature_support":
                assert reason

            se_protein_pairs.add( (s,p,action,kind,reason) )
            seen.add((s,p))
    
    se_protein_pairs = list(se_protein_pairs)
    print >> sys.stderr, "Loaded", len(se_protein_pairs), "se-protein-pairs for mouse"

    return se_protein_pairs, p2pn, s2sn



def load_general_info():

    target2drugs = defaultdict(list)
    target2target = {}

    for line in open("interactions.tsv"):
        (drug, target) = line.strip("\n").split("\t")
        drug = int(drug)
        
        for _target in target.split("_"):

            if target.split("@")[0] in protein_blacklist:
                continue

            target2drugs[_target].append(drug)
            target2target[_target] = target

    assert target2drugs

    se2drugs = defaultdict(list)
    for line in open("side_effects.tsv"):
        (drug, concept, _, _) = line.strip("\n").split("\t")

        se2drugs[concept].append(int(drug))
    
    assert se2drugs
    
    return target2drugs, se2drugs, target2target
    

def main():

    import doctest
    doctest.testmod()
    
    target2drugs, se2drugs, _ = load_general_info()
    
    load_known_action = lambda : load_annotated_data(rx = r"^caused_by_" )
    load_treated_known_action = lambda : load_annotated_data(rx = r"^treated_by_" )
    load_all_annotated = lambda _ : load_annotated_data(rx = r"^caused" )
    load_combined = lambda : load_annotated_data(rx = r"_by_|not relevant" )

    for print_full in (True, False):
        for actinh in (True, False):

            filenames = []
            # for prefix, load in (("whitebread", load_whitebread_data), ("caused", load_all_annotated), ("known_action", load_known_action), ("treated_known_action", load_treated_known_action),("combined", load_combined),):
            for prefix, load in (("whitebread", load_whitebread_data), ("caused", load_all_annotated), ("actinh" if actinh else "mouse", load_mouse_data),):
                se_protein_pairs, p2pn, s2sn = load(print_full)
            
                full_prefix = ("full_" if print_full else "") + ("under_" if over_under == "under" else "") + prefix
                    
                filenames.append( print_pairs(se_protein_pairs, full_prefix, target2drugs, se2drugs, p2pn, s2sn, print_full) )
            
            if actinh: continue
            if not print_full: continue
            
            print >> sys.stderr, "Merging results"
            
            print >> sys.stderr, " ".join(filenames)
            
            
            prefix = ("full_" if print_full else "") + ("under_" if over_under == "under" else "") + "merged"
            filename = "pred."+prefix+".se_protein_known.all.txt"
            fh_merged = open(filename, "w")
            
            for lines in zip( *[ open(fn) for fn in filenames ] ):
                
                marker = "no"
                remainder = None
                
                for line in lines:
                    _marker, _remainder = line.split("\t", 1)

                    if priorities[marker] > priorities[_marker]: continue
                    
                    marker = _marker

                    if remainder is not None:
                        if _remainder != remainder:
                            r = remainder.split("\t")
                            _r = _remainder.split("\t")
                            
                            assert r[0] == _r[0]
                            assert r[4] == _r[4]

                    remainder = _remainder

                print >> fh_merged, "%s\t%s" % (marker, remainder),
                    
                
    
    



    

if __name__ == '__main__':
    main()

