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

    Optional commands:
    --verbose
    --auto_thresh

"""

import os, sys, re
from argParser import *
ap = ArgParser(sys.argv)

v = ap.getOptionalToggle("--verbose")
if v != False:
    v = True
fastapath = ap.getArg("--fastapath")
seedseq = ap.getOptionalArg("--seed")
cthresh = ap.getOptionalArg("--cthreshold") # sites with N indels >= this value will be culled.
outpath = ap.getArg("--outpath")
indel_max = ap.getOptionalArg("--indelmax") # any species with this proportion of indels, or higher, will be culled from the analysis.
if indel_max == False:
    indel_max = 1.0
if cthresh == False:
    cthresh = 1.0

if seedseq == False and cthresh == False:
    print "\n. You need to specify either --seedseq or --cthreshold"
    exit()

if seedseq != False and cthresh != False:
    print "\n. You cannot specify both --seedseq and --cthreshold"
    exit()


def read_taxa_seq(msapath):
    fin = open(msapath, "r")
    lines = fin.readlines()
    fin.close()
    
    # Read the sequences from the original MSA. . .
    taxa_seq = {}
    lastseq = None
    for l in lines:
        l = l.strip()
        if l.startswith(">"):
            lastseq = re.sub("\>", "", l)
            taxa_seq[lastseq] = ""
        elif lastseq != None and l.__len__() > 10:
            taxa_seq[lastseq] += l
    return taxa_seq

def write_fasta(taxa_seq):
    fout = open(outpath, "w")
    for taxa in taxa_seq:
        fout.write(">" + taxa + "\n")
        fout.write(taxa_seq[taxa] + "\n")
    fout.close()

def trim_nonseed_sites(taxa_seq, seedseq):
    # Place trimmed sequences in taxa_newseq. . .
    taxa_newseq = {}
    if seedseq not in taxa_seq:
        print "\n. Hmm, I cannot find", seedseq, "in your alignment."
        exit()
    for taxa in taxa_seq:
        taxa_newseq[taxa] = ""
    for site in range(0, taxa_seq[seedseq].__len__()):
        #print site
        if taxa_seq[seedseq][site] != "-":
            for taxa in taxa_seq:
                taxa_newseq[taxa] += taxa_seq[taxa][site]
    return taxa_newseq

def trim_indel_sites(taxa_seq, cthresh, verbose=False):
    #print "84", taxa_seq
    taxa_newseq = {}
    if verbose:
        print "\n. Culling indel-rich sites from", fastapath
    for taxa in taxa_seq:
        taxa_newseq[taxa] = ""
    for site in range(0, taxa_seq[ taxa_seq.keys()[0] ].__len__()):
        count_indels = 0
        for taxa in taxa_seq:
            #print taxa_seq[taxa][site]
            if taxa_seq[taxa][site] == "-":
                count_indels += 1
        cprop = float(count_indels) / taxa_seq.keys().__len__()
        if cprop >= cthresh:
            if verbose:
                print "--> Culling site", site, "with %.2f"%cprop,"indels."
        else:
            if verbose:
                print "--> Keeping site", site, "with %.2f"%cprop,"indels."
            for taxa in taxa_seq:
                taxa_newseq[taxa] += taxa_seq[taxa][site]         
    if taxa_newseq[ taxa_newseq.keys()[0] ].__len__() == 0:
        print "\n. Hmmm... your trimming criteria are too stringent.  No sequence data remains."
        return None    
    if taxa_newseq.__len__() == 0:
        print "\n. Hmmm, your trimming criteria are too stringent.  No sequence data remains."
        return None
    #print "108", taxa_newseq
    return taxa_newseq

def trim_taxa(taxa_seq, indel_max, verbose=False):
    taxa_newseq = {}
    if verbose:
        print "\n. Culling indel-rich sequences from", fastapath
    for taxa in taxa_seq:
        count = taxa_seq[taxa].count("-")
        p = float(count) / taxa_seq[taxa].__len__()
        if p >= indel_max:
            if verbose:
                print "--> Removing", taxa, "with %.1f"%(p * 100), "% indels."
        else:
            taxa_newseq[taxa] = taxa_seq[taxa]
    return taxa_newseq

def get_stats_taxa_seq(taxa_seq):
    ntaxa = taxa_seq.keys().__len__()
    nsites = taxa_seq[ taxa_seq.keys()[0] ].__len__()
    countchars = 0
    countindels = 0
    for taxa in taxa_seq:
        for char in taxa_seq[taxa]:
            countchars += 1
            if char == "-":
                countindels += 1
    return [ntaxa, nsites, countindels, countchars]
    

taxa_seq = read_taxa_seq(fastapath)
if seedseq != False:
    taxa_seq = trim_nonseed_sites(taxa_seq, seedseq)
    
if taxa_seq.__len__() == 0:
    print "\n. Hmmm, something is wrong with your input."
    exit()

if ap.getOptionalToggle("--auto_thresh"):
    last = 0.0
    if cthresh != False:
        cstep = 0.01
        c = 0.0
        print "threshold ntaxa nsites nindels nchars propindels"
        while(c < float(cthresh) ):
            c += cstep
            taxa_newseq = trim_indel_sites(taxa_seq, c )
            if taxa_newseq == None:
                continue
            print taxa_newseq
            taxa_newseq = trim_taxa(taxa_newseq, float(indel_max))
            x = get_stats_taxa_seq(taxa_newseq)
            print c, x[0], x[1], x[2], x[3], "%.3f"%(float(x[2])/x[3])
            this = x[3]
            if this >= last:
                last = this
            else:
                cthresh = c - cstep
                print ". Setting per-site conservation threshold to", cthresh
                c = cthresh + 1
        if c == float(cthresh):
            print ". Setting per-site conservation threshold to", c
            print ". The search for an optimal threshold stopped early because you specified"
            print "  the upper limit to be", cthresh, "using the --cthreshold command."


c = float(cthresh)

taxa_newseq = trim_indel_sites(taxa_seq, c, verbose=v)
taxa_newseq = trim_taxa(taxa_newseq, float(indel_max), verbose=v)
x = get_stats_taxa_seq(taxa_newseq)
print "\n. Writing", x[0], "taxa with", x[1], "sites."
write_fasta(taxa_newseq)