"""
    trim_gaps.py
    
    This script will remove sites (i.e. columns) from an alignment
    that lack homologous sequence (i.e. they have indels) compared
    to a "seed" sequence.
    
    COMMANDS:
    --fastapath P, where P is a file path to a fasta-formatted alignment
    --seed S, where S is a taxa name found in P.
    --outpath R, where R is the filepath to which the trimmed alignment will be written.
    --indelmax D, (optional) where D is the upper limit proportion of indels allowed in a sequence.
                    Taxa exceeding this limit will be culled.

"""

import os, sys, re
from argParser import *
ap = ArgParser(sys.argv)

fastapath = ap.getArg("--fastapath")
seedseq = ap.getArg("--seed")
outpath = ap.getArg("--outpath")
indel_max = ap.getOptionalArg("--indelmax") # any species with this proportion of indels, or higher, will be culled from the analysis.

fin = open(fastapath, "r")
lines = fin.readlines()
fin.close()

taxa_seq = {}
lastseq = None
for l in lines:
    l = l.strip()
    if l.startswith(">"):
        lastseq = re.sub("\>", "", l)
        taxa_seq[lastseq] = ""
    elif lastseq != None and l.__len__() > 10:
        taxa_seq[lastseq] += l

if seedseq not in taxa_seq:
    print "\n. Hmm, I cannot find", seedseq, "in your alignment."
    exit()

# Place trim sequences in taxa_newseq
taxa_newseq = {}
for taxa in taxa_seq:
    taxa_newseq[taxa] = ""
for site in range(0, taxa_seq[seedseq].__len__()):
    #print site
    if taxa_seq[seedseq][site] != "-":
        for taxa in taxa_seq:
            taxa_newseq[taxa] += taxa_seq[taxa][site]
    
# Cull taxa with too many indels:
taxa_finalseq = {}
if indel_max != False:
    print "\n. Culling indel-rich sequences from", fastapath
    indel_max = float(indel_max)
    for taxa in taxa_newseq:
        count = taxa_newseq[taxa].count("-")
        p = float(count) / taxa_newseq[taxa].__len__()
        #print taxa, p
        if p >= indel_max:
            print ". Removing", taxa, "with %.1f"%(p * 100), "% indels."
        else:
            taxa_finalseq[taxa] = taxa_newseq[taxa]

fout = open(outpath, "w")
for taxa in taxa_finalseq:
    fout.write(">" + taxa + "\n")
    fout.write(taxa_finalseq[taxa] + "\n")
fout.close()