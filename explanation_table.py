#!/usr/bin/env python2.5
# encoding: utf-8

from __future__ import division
from __future__ import with_statement

import sys
import os
import re
import itertools

from collections import defaultdict

re_clean_targets = re.compile(r"(ENSP\d+)")

def load_se2protein(cutoff):
    
    se2protein = defaultdict(list)
    
    fh_in = open("protein_se_pv.tsv")
    fh_in.next()

    for line in fh_in:

        se, protein, _, qv = line.strip("\n").split("\t")

        qv = float(qv)

        if qv >= cutoff:
            continue
        
        se2protein[se].append(protein)

    print >> sys.stderr, "Number of SE with q-values:", len(se2protein)

    return se2protein 

def get_go(protein2go, target):
    
    go = set()
    
    for _target in re.split(r"(?:@...)?_?", target):
        go.update(protein2go.get(_target, ()))

    return go


def load_go_classifications(protein_drug_proteins):
    
    protein2go = defaultdict(list)
    ds = set()
    
    # first, load GO classifications for each target
    for line in open("go_classification.tsv"):
        (protein, go_term, description) = line.strip("\n").split("\t")
        protein2go[protein].append(description)
        ds.add(description)
        ds.add("main " + description)
        
    for d in ds:
        protein_drug_proteins[d] = defaultdict(list)
        
    main_targets = protein_drug_proteins["main targets"]
        
    # then, go through all drug-target relations and distribute the targets into the GO categories
    for (drug, targets) in protein_drug_proteins["all targets"].iteritems():
        for target in targets:
            is_main_target = drug in main_targets and target in main_targets[drug] 
            for go in get_go(protein2go, target):
                protein_drug_proteins[go][drug].append(target)
                if is_main_target:
                    protein_drug_proteins["main "+go][drug].append(target)
                    
    
    
def load_drug_targets():

    drug2targets = defaultdict(list)

    for line in open("interactions.tsv"):
        (drug, target) = line.strip("\n").split("\t")
        drug = int(drug)
        drug2targets[drug].append(target)

    return drug2targets
    

def load_main_targets():

    drug2main_targets = defaultdict(list)

    for line in open("actions/main_targets.tsv"):
        (protein, _chemical) = line.strip("\n").split("\t")
        chemical = int(_chemical[3:])
        drug2main_targets[chemical].append(protein)
            
    return drug2main_targets
            
            
def load_protein_names():
    
    protein_names = {}
    protein_annotations = {}

    for line in open("protein_names"):
        (protein, name, annotation) = line.strip("\n").split("\t")
        protein_names[protein] = name
        protein_annotations[protein] = annotation
      
    assert protein_annotations
    assert protein_names 
      
    return protein_names, protein_annotations


def filter_by_re(d2t, search, parents, protein_annotations):
    
    matched_proteins = set()
    
    for (pp, parents) in parents.iteritems():
        for p in pp.split("_"):
            
            p = p.split("@")[0]
            
            if search( protein_annotations.get(p, "") ):
                matched_proteins.update(parents)
    
    print >> sys.stderr, "Number of matched proteins", len(matched_proteins)
    
    
    drugs2targets = defaultdict(list)
    
    for (d, ts) in d2t.iteritems():
        result = matched_proteins.intersection(ts)
        
        if result:
            drugs2targets[d] = result

    return drugs2targets


def dict_map(f, d1, d2, cast = set):
    
    assert d1
    assert d2
    
    d = {}
    
    for k in set( d1.keys() ).union( d2.keys() ):
        v1 = d1.get(k, ())
        v2 = d2.get(k, ())
        result = f( cast(v1), cast(v2) )
        
        if result:
            d[k] = result
        
    return d
    
def consolidate(protein_drug_proteins):
    
    ## do a majority vote among main targets to only keep the most important one
    
    keys = ("main GPCR", "main non-kinase non-metabolizing enzyme", "main kinase", "main ion channel", "main nuclear receptor")
    
    drugs = set()
    
    for k in keys: drugs.update( protein_drug_proteins[k].keys() )
    
    for drug in drugs:
        
        m = max( len(protein_drug_proteins[k].get(drug, ())) for k in keys )
        
        for k in keys:
            if drug not in protein_drug_proteins[k]: continue
            if len(protein_drug_proteins[k][drug]) < m:
                protein_drug_proteins[k].pop(drug)
        


def load_drug_protein_mapping():
    
    protein_drug_proteins = defaultdict(lambda : defaultdict(list))

    protein_names, protein_annotations = load_protein_names()

    protein_drug_proteins["all targets"] = protein_drug_proteins["targets"] = load_drug_targets()

    # need to have all targets to load further subclasses
    protein_drug_proteins["main targets"] = load_main_targets()
    
    load_go_classifications(protein_drug_proteins)    

    # re_metabolizing = re.compile(r"Multidrug|^CYP[0-9]|ATP.binding.cassette|Thromboxane.A|Glutathione.S.transferase", re.I).search
    # re.compile(r"Multidrug|Cytochrome|ATP.binding.cassette|Thromboxane.A|arachidonate.5.lipoxygenase|Glutathione.S.transferase", re.I).search
    # protein_drug_proteins["metabolizing proteins"] = filter_by_re( protein_drug_proteins["targets"], re_metabolizing, parents, protein_annotations )
    metabolizing_proteins = set(line.strip().split("\t")[0] for line in open("metabolizing_proteins"))
    
    for drug, targets in protein_drug_proteins["targets"].items():
        l = []
        for target in targets:
            if any( t in metabolizing_proteins for t in re_clean_targets.findall(target) ): 
                l.append(target)
        
        if l:
            protein_drug_proteins["metabolizing proteins"][drug] = l

    protein_drug_proteins["main targets"] = dict_map( set.difference, protein_drug_proteins["main targets"], protein_drug_proteins["metabolizing proteins"]  ) # remove metabolizing proteins from main targets
    protein_drug_proteins["non-metabolizing proteins"] = dict_map( set.difference, protein_drug_proteins["targets"], protein_drug_proteins["metabolizing proteins"]  )
    protein_drug_proteins["main non-metabolizing proteins"] = dict_map( set.difference, protein_drug_proteins["main targets"], protein_drug_proteins["metabolizing proteins"]  )
    
    x = dict_map( set.intersection, protein_drug_proteins["GPCR"], protein_drug_proteins["metabolizing proteins"]  )
    assert not x

    x = dict_map( set.difference, protein_drug_proteins["GPCR"], protein_drug_proteins["non-metabolizing proteins"]  )
    assert not x
    
    protein_drug_proteins["off-targets"] = dict_map( set.difference, protein_drug_proteins["non-metabolizing proteins"], protein_drug_proteins["main targets"] )

    protein_drug_proteins["non-kinase non-metabolizing enzyme"] = dict_map( set.difference, protein_drug_proteins["non-kinase enzyme"], protein_drug_proteins["metabolizing proteins"]  )
    protein_drug_proteins["main non-kinase non-metabolizing enzyme"] = dict_map( set.difference, protein_drug_proteins["main non-kinase enzyme"], protein_drug_proteins["metabolizing proteins"]  )

    for s in ("GPCR", "non-kinase non-metabolizing enzyme", "kinase", "ion channel", "nuclear receptor", ):
        protein_drug_proteins["!"+s] = dict_map( set.difference, protein_drug_proteins["non-metabolizing proteins"], protein_drug_proteins[s] )
        
        
    protein_drug_proteins = dict(protein_drug_proteins)

    for (s, d2t) in protein_drug_proteins.items():
        assert all(d2t.values())
        # print >> sys.stderr, "   ", s, len(d2t)
        # print >> sys.stderr, d2t.items()[:3]
        
    consolidate(protein_drug_proteins)
    
    return protein_drug_proteins
    

def merge_redundant_proteins(l, merging_cache = {}):
    """
    >>> merge_redundant_proteins(["ENSP1", "ENSP2", "ENSP3"])
    3
    >>> merge_redundant_proteins(["ENSP1_ENSP2", "ENSP2", "ENSP3"])
    2
    >>> merge_redundant_proteins(["ENSP1_ENSP2_ENSP3", "ENSP2", "ENSP3"])
    1
    >>> merge_redundant_proteins(["ENSP1_ENSP2", "ENSP2_ENSP3", "ENSP3"])
    1
    >>> merge_redundant_proteins(['ENSP00000000442_ENSP00000252542_ENSP00000261532_ENSP00000292123_ENSP00000339587_ENSP00000354584', 'D047629@act_ENSP00000343925@act', 'ENSP00000000442@act_ENSP00000252542@act_ENSP00000261532@act_ENSP00000292123@act_ENSP00000339587@act_ENSP00000354584@act'])
    2
    """

    c = merging_cache.get(tuple(l))
    if c is not None: return c
    
    clusters = [ set(re_clean_targets.findall(s)) for s in l ]
    
    n = 0
    
    while len(clusters) != n:
        
        n = len(clusters)
        new_clusters = []
        
        for i in range(len(clusters)):
            current = clusters[i]
            if current is None: continue
            new_clusters.append(current)
            for j in range(i+1,len(clusters)):
                if clusters[j] is not None and current.intersection(clusters[j]):
                    # print >> sys.stderr, "Merge", current, clusters[j]
                    
                    current.update(clusters[j])
                    clusters[j] = None
        
        clusters = new_clusters
    
    merging_cache[tuple(l)] = n

    return n


def calc(se2drugs, se2proteins, drug2proteins, drug_set, se_set, prefix, fh_out, se_with_explanations):
    
    store_se_with_explanations = se_with_explanations is not None and len(se_with_explanations) == 0

    total_se_drug_pairs = 0
    explained_se_drug_pairs = []
    
    n_explained_se_one = 0
    n_explained_se_half = 0
    n_se = 0
    
    
    # for all SE, count how many SE-drug pairs can be explained
    for (se, drugs) in se2drugs.items(): 

        if se not in se_set: 
            continue


        assert drugs

        n_se += 1

        n_explained = 0
        n_drugs_in_set = 0
        
        se_proteins = se2proteins[se]
        
        l = []
        
        for drug in drugs:
            if drug not in drug_set: continue

            n_drugs_in_set += 1
            total_se_drug_pairs += 1

            explaining_proteins = [ protein for protein in drug2proteins.get(drug, ()) if protein in se_proteins ]
            n_explaining_proteins = merge_redundant_proteins( explaining_proteins )

            if n_explaining_proteins:
                explained_se_drug_pairs.append( (se, drug) )
                n_explained += 1
                
                l.append(drug)
                
                if store_se_with_explanations:
                    se_with_explanations.add(se)
        
        if prefix and n_drugs_in_set:
            print >> fh_out, "%s\t%d\t%s\t%d\t%d\t%s" % (prefix, se in se_with_explanations, se, n_drugs_in_set, n_explained, " ".join(map(str,l)))
            # print >> fh_out, "%s\t%d\t%s\t%d\t%d" % (prefix, se in se_with_explanations, se, n_drugs_in_set, n_explained)
        
        if n_explained:
            n_explained_se_one += 1
            
            if n_explained * 2 >= n_drugs_in_set:
                n_explained_se_half += 1
        
            
            
                
                
    return (explained_se_drug_pairs, total_se_drug_pairs, n_explained_se_one, n_explained_se_half, n_se)
        

def load_se2drugs():

    se2drugs = defaultdict(list)

    # collect all SE-drug relationships (from the labels)
    for line in open("side_effects.tsv"):
        (_cid, se) = line.strip().split("\t")[:2]
        cid = int(_cid)
        se2drugs[se].append(cid)

    assert se2drugs
        
    return se2drugs


def main():

    n_jobs = 1
    this_job = 0

    if "submit" in sys.argv:
        submit()
        return
    elif len(sys.argv) == 3:
        n_jobs = int(sys.argv[2])
        this_job = int(sys.argv[1])


    bins = (
        # ("all", "LOW BIN1 BIN2 BIN3 HIGH"),
        ("all", "BIN1 BIN2 BIN3 HIGH"),
        ("rare", "BIN1"),
        ("average", "BIN2"),
        ("frequent", "BIN3"),
        # ("excl low/high", "BIN1 BIN2 BIN3"),
    )

    # the first line for each protein / drug requirement should be the one to base "explained SE" on
    # -> : use protein requirement as drug requirement
    
    to_analyze = """
main non-metabolizing proteins	non-metabolizing proteins
main non-metabolizing proteins	GPCR
main non-metabolizing proteins	non-kinase non-metabolizing enzyme
# non-metabolizing proteins	metabolizing proteins
main non-metabolizing proteins	ion channel
main non-metabolizing proteins	kinase
main non-metabolizing proteins	nuclear receptor
main GPCR	non-metabolizing proteins
main GPCR	GPCR
main GPCR	!GPCR
main ion channel	non-metabolizing proteins
main ion channel	ion channel
main ion channel	!ion channel
main non-kinase non-metabolizing enzyme	non-metabolizing proteins
main non-kinase non-metabolizing enzyme	non-kinase non-metabolizing enzyme
main non-kinase non-metabolizing enzyme	!non-kinase non-metabolizing enzyme
# metabolizing proteins	all targets
# main kinase	non-metabolizing proteins
# main kinase	kinase
# main kinase	!kinase
main nuclear receptor	non-metabolizing proteins
main nuclear receptor	nuclear receptor
main nuclear receptor	!nuclear receptor
main targets	non-metabolizing proteins
main targets	main targets
# main targets	metabolizing proteins
# main targets	non-metabolizing proteins
main targets	off-targets
"""
    
    
    # check if the anaysis options are correctly formated
    for l in to_analyze.strip().split("\n"):
        if l.startswith("#"): continue
        assert len(l.split("\t")) == 2


    # for cutoff in (0.01, 0.05):
    for cutoff in (0.01,):

        se2proteins = load_se2protein(cutoff)
        se_set = set( se2proteins.keys() )

        fh_results = open("explanation_table.%g.tsv" % cutoff, "w")
        fh_all_se = open("explained_se.%g.tsv" % cutoff, "w")
    
        se_with_explanations = {}
    
        last_drug_requirement = None

        se2drugs = load_se2drugs()

        fh_drug_count = open("/g/bork8/mkuhn/data/se_protein/drug_counts.%g.tsv" % cutoff, "w") if this_job == 0 else None
        # go through each type of analysis
        for l in to_analyze.strip().split("\n"):

            if l.startswith("#"): continue

            (drug_requirement, protein_requirement) = l.split("\t")

            if drug_requirement == "->":
                drug_requirement = protein_requirement

            protein_drug_proteins = load_drug_protein_mapping()

            if protein_requirement not in protein_drug_proteins:
                print >> sys.stderr, "don't know drug-protein mapping for: %s, skipping!" % protein_requirement
                continue


            drug2proteins = protein_drug_proteins[protein_requirement]

            drug_set = set( protein_drug_proteins[drug_requirement] )

            drug_blacklist = set( int(s) for s in open("hobohm_remove_list.txt") )
            drug_set = drug_set.difference(drug_blacklist) 

            if last_drug_requirement != drug_requirement:
                last_drug_requirement = drug_requirement

            print >> sys.stderr, "number of drugs with", drug_requirement, "=", len(drug_set)
        
            if fh_drug_count:
                print >> fh_drug_count, "%s\t%d" % (drug_requirement, len(drug_set))

            prefix = ""
            current_se_with_explanations = set()
            
            prefix = "\t".join( (drug_requirement, protein_requirement) )
            


            explained_se_drug_pairs, total_se_drug_pairs, n_explained_se_one, n_explained_se_half, n_se = calc(se2drugs, se2proteins, drug2proteins, drug_set, se_set, prefix, fh_all_se, current_se_with_explanations)

            assert "C0018681" not in current_se_with_explanations
            assert current_se_with_explanations

            print >> fh_results, "\t".join( map(str, (drug_requirement, protein_requirement, len(explained_se_drug_pairs), total_se_drug_pairs, n_explained_se_one, n_explained_se_half, n_se, ) ) )
    


      
        

        



if __name__ == '__main__':
    main()
