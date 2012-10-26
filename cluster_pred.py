#!/g/bork3/bin/python

# fisher.test( matrix(c(12,77,132,1901), nrow=2), alternative="g")

from __future__ import division
from __future__ import with_statement

import sys, os

import time
import math
import functools, itertools

from collections import defaultdict

from Bio.Cluster import *

import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)
                                         
import re

from known_protein_se_benchmark import load_general_info, priorities

cutoff = 1/2

metabolizing_proteins = set(line.strip().split("\t")[0] for line in open("metabolizing_proteins"))
protein_blacklist = list( line.strip() for line in open("protein.blacklist") ) 


def se_target_corr(drugs, proteins, target2drugs, target2target, p2pn):
    
    drugs = set(drugs)
    corr = {}

    for p1 in proteins:
        for p2 in proteins:
            sim = 0
    
            if target2target[p1] == target2target[p2]:
                sim = 2 
            elif p1.split("@")[0] == p2.split("@")[0]:
                sim = 1.5
            else:
                d1 = drugs.intersection(target2drugs[p1])
                # d1 = set(target2drugs[p1])
                n1 = len(d1)

                d2 = drugs.intersection(target2drugs[p2])
                # d2 = set(target2drugs[p2])
                n2 = len(d2)

                intersection = len( d1.intersection(d2) )

                sim = max( (intersection / n1) if n1 else 0, (intersection / n2) if n2 else 0 )
        
                # print p2pn[p1], p2pn[p2], intersection / n1
                # print p2pn[p1], p2pn[p2], sim, intersection / n1, intersection / n2
            
            corr[p1, p2] = sim
    
    return corr


def gather(tree, root):
    
    if root >= 0:
        return [root]
    else:
        n = tree[-root-1]
        c = []
        c.extend( gather(tree, n.left) )
        c.extend( gather(tree, n.right) )
        return c

def getclusters(tree, root):
    
    if root >= 0:
        return [[root]]
    
    n = tree[-root-1]
    
    # print "x", n
    
    if n.distance < 2 - cutoff:
        c = []
        c.extend( gather(tree, n.left) )
        c.extend( gather(tree, n.right) )
        return (c,)
    else:
        c = []
        c.extend( getclusters(tree, n.left))
        c.extend( getclusters(tree, n.right))
        return c
        
    
def cluster_groups(proteins, corr):

    if not proteins: return ()
    
    proteins = list(proteins)
    
    # print proteins
    
    if len(proteins) == 1:
        return (proteins,)
    
    n_nodes = len(proteins)
    distancematrix = []
    
    for i in range(n_nodes):
        p1 = proteins[i]
        for j in range(0, i):
            p2 = proteins[j]
            
            sim = corr[p1,p2]
            distancematrix.append(2-sim)

    tree = treecluster( distancematrix = distancematrix, method = "m" )    

    # print tree

    clusters = getclusters(tree, -len(tree))
    pclusters = []
        
    for c in clusters:
        pclusters.append( [proteins[i] for i in c] )

    # print pclusters

    return pclusters





def main():
    
    target2drugs, se2drugs, target2target = load_general_info()
    
    proteins_per_se = defaultdict(set)
    
    s2sn = {}
    p2pn = {}
    
    known_sp = {}
    
    pvalues = {}

    counts = defaultdict(int)

    do_only = ""
    # do_only = "@inh"
    # do_only = "@act"

    pvalue_cutoff = 0.01

    # fn = "full_under_mouse"
    fn = "full_merged"
    # fn = "mouse"
    # fn = "full_actinh"

    if len(sys.argv) == 4:
        (fn, do_only, _pvalue_cutoff) = sys.argv[1:]
        pvalue_cutoff = float(_pvalue_cutoff)
        
        if not do_only.startswith("@"):
            do_only = ""
        
    
    fn_descr = fn + do_only.replace("@",".") + "." + str(pvalue_cutoff)
    
    fh_out = open("pred_lines.%s.tsv" % fn_descr, "w")

    mouse_se = set( x.strip() for x in os.popen("cut -f 1 mouse_se.tsv | sort -u") ) 

    for line in os.popen("cat pred.%s.se_protein_known.all.txt" % fn):
        (known, pvalue, sn, pn, s, p, reason)  = line.strip("\n").split("\t")
        
        if p.split("@")[0] in protein_blacklist:
            continue

        # if do_only and "inh" not in p: continue
        # if s != "C0039231": continue
        
        # if sn != "priapism": continue
        # if "SCN5A" not in pn: continue
        
        if float(pvalue) > pvalue_cutoff: continue
        # if float(contribution) < 20: continue
        # if float(pvalue) == 1: continue

        print >> fh_out, line,

        proteins_per_se[s].add(p)
        
        s2sn.setdefault(s, sn)
        p2pn.setdefault(p, pn)
        
        pvalues[s,p] = float(pvalue)
        
        if known == "literature_support":
            assert reason

        if known != "no":
            known_sp[s,p] = known, reason
    
    # print known_sp
    
    fh_out = open("pred.%s.tsv" % fn_descr, "w")
    
    se_outcome = {}
    
    first = True
    
    for (se, proteins) in sorted(proteins_per_se.items()):
        
        drugs = se2drugs[se]
        corr = se_target_corr(drugs, proteins, target2drugs, target2target, p2pn)
        
        explaining_proteins = [ p for p in proteins if (se,p) in known_sp ]
        explaining_clusters = [ tuple(x) for x in cluster_groups(explaining_proteins, corr) ]
        
        # print "explaining_clusters", explaining_clusters
        
        protein_cluster_dict = defaultdict(list)
        
        remaining_proteins = []
        
        for protein in proteins:
            m = 0
            mc = None
            
            for explaining_cluster in explaining_clusters:
                for explaining_protein in explaining_cluster:
                    sim = corr[explaining_protein, protein]
                    if sim > m:
                        m = sim
                        mc = explaining_cluster
                        
            if m > cutoff:
                protein_cluster_dict[mc].append(protein)
            else:
                remaining_proteins.append(protein)
                
        protein_clusters = protein_cluster_dict.values()

        if remaining_proteins:
            protein_clusters.extend( cluster_groups(remaining_proteins, corr) )
        
        for cluster in protein_clusters:
            
            known = "no"
            reason = ""
            
            rproteins = defaultdict(set)
            
            for p in cluster:
                k, r = known_sp.get((se,p), (None,""))
                if k:
                    # don't overwrite "same_phenotype" with "related_phenotype"
                    if priorities[known] < priorities[k]:
                        known = k
                        reason = r
                    rproteins[k].add(p)

            if known == "literature_support":
                assert reason
            
            cluster = sorted(cluster, key = lambda p : pvalues[se,p] )
            
            _epns = sorted(p2pn[p] for p in cluster if p in rproteins[known])
            _base_epns = set(p.split("@")[0] for p in cluster if p in rproteins[known])

            epns = " ".join(_epns)
            pns = " ".join(sorted(p2pn[p] for p in cluster))
            
            only_metabolzing = all( p.split("@")[0] in metabolizing_proteins for p in cluster )
            _only_metabolizing = "metabolizing" if only_metabolzing else ""
            
            # for require_known in (True, False):
            #     for require_at in (True, False):
            #         for p in cluster:
            #             if require_at and "@" not in p: continue
            #             if require_known and p not in rproteins[known]: continue
            
            pv = min( pvalues[se,p] for p in cluster )
            
            assert pv == pvalues[se, cluster[0]]
            
            lpv = log10(pv)
            
            if do_only and do_only not in pns + epns: continue
            
            print "\t".join((known, "%.3g" % lpv, se, s2sn[se], epns, pns, reason, " ".join(cluster), _only_metabolizing))
            
            m = "mouse" if se in mouse_se else "other" 
            ai = " ".join( set(re.findall(r"@(act|inh)", epns + pns ) ) )
            
            
            counts[known] += 1
            
            if not only_metabolzing:
                if se not in se_outcome or priorities[ se_outcome[se] ] < priorities[known]:
                    se_outcome[se] = known

            if first:
                print >> fh_out, "\t".join(("log10(min qvalue)", "found in mouse?", "act/inh", "status", "side effect", "explaining proteins", "comment", "UMLS code of SE", "all proteins", "only metabolizing?"))
            
            if known == "no":
                known = "not verified"
            elif known == "best":
                known = "literature support"
            
            print >> fh_out, "\t".join(("%.3g" % lpv, m, ai, known.replace("_", " "), s2sn[se], epns, reason, se, pns, _only_metabolizing))
            
            first = False
            
    print >> sys.stderr, "Same phenotype:", counts["same_phenotype"] 
    print >> sys.stderr, "Total SE / protein-cluster pairs:", sum(counts.values())
    print >> sys.stderr, "Number of different SE:", len(proteins_per_se)
    print >> sys.stderr
    
    
    se_count = defaultdict(int)
    for known in se_outcome.values():
        se_count[known] += 1
     
    single_supporting_evidence = r"literature_support|best|related_phenotype|opposite_phenotype|same_phenotype|yes".split("|")
    
    n_single_strong_support = sum( [ se_count[s] for s in single_supporting_evidence ] )
    
    n_some_evidence = se_count["anecdotal_evidence"]
    
    n_wrong = sum( [se_count[s] for s in ("wrong", "hypersensitivity_reaction")] )
    n_tm = sum( [se_count[s] for s in ("false_positive",)] )
    
    n_unknown =  sum( [se_count[s] for s in ("unknown",)] )
    n_not_checked =  sum( [se_count[s] for s in ("no",)] )
    
    n_total = sum( se_count.values() )
    
    check_total = sum( ( n_single_strong_support, n_some_evidence, n_tm, n_wrong, n_unknown, n_not_checked) )
    
    assert n_total == check_total, n_total - check_total
    
    print >> sys.stderr, "\t".join( map(str, (n_single_strong_support, n_some_evidence, n_wrong, n_unknown, n_total) ) )
    print >> sys.stderr, "all checked?", se_count["no"] == 0
    
    all_checked = "yes" if se_count["no"] == 0 else "no"
    
    se_kind = "all" if fn.startswith("full") else fn
    
    act_inh = do_only[1:] if do_only.startswith("@") else "-"
    
    fh_summary = open("summary.tsv", "a")
    
    print >> fh_summary, "\t".join( map(str, (pvalue_cutoff, se_kind, act_inh, all_checked, n_single_strong_support, n_some_evidence, n_unknown, n_not_checked, n_tm, n_wrong, n_total) ) )
    
    
    
    

if __name__ == '__main__':
    
    main()

    
