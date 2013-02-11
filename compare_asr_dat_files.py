#
# USAGE:
#
# python compare_asr_dat_files.py [<id> <dat filepath> . . .]
#
# CREATES THESE FILES:
unaligned_fasta_path = "ancestors.unaligned.fasta"
#

import os
import sys
import re

def getprobs(inpath):
    fin = open(inpath, "r")
    lines = fin.readlines()
    fin.close()

    site_states_probs = {}
    for l in lines:
        tokens = l.split()
        site = int(tokens[0])
        site_states_probs[ site ] = {}
        i = 1
        while i < tokens.__len__():
            s = tokens[i]
            p = float(tokens[i+1])
            foundgap = False
            if p > 1.0:
                p = 1.0
                foundgap = True
            site_states_probs[site][s] = p
            i += 2
            if foundgap:
                i = tokens.__len__() # stop early
    return site_states_probs
    
#
# Capture command-line parameters
#
dat_files = {} # key = unique name, value = path to .dat file
i = 1
while i < sys.argv.__len__():
    name = sys.argv[i]
    i += 1
    path = sys.argv[i]
    dat_files[name] = path
    if False == os.path.exists( path ):
        print "\n\nERROR: the following path does not seem to exist:", path
        exit()
    i += 1
    print "I found a file:", path, "with ID", name

# sanity check:
#if (i <= 2):
#    print "\n\n! ERROR: Hey, you need to specify at least two *.dat files"
#    exit()
    
# parse DAT files:
anc_data = {}
for anc in dat_files.keys():
    anc_data[anc] = getprobs( dat_files[anc] )

#
def get_ml_state(states_probs):
    maxp = None
    maxc = None
    for state in states_probs:
        if maxp == None:
            maxp = states_probs[state]
            maxc = state
        elif maxp < states_probs[state]:
            maxp = states_probs[state]
            maxc = state
    return maxc

# compare dat files:
countindel = 0
countred = 0 # different states
countorange = 0 # diff. states, but pairwise secondary states
countgreen = 0 # same state, but adds uncertainty
redsites = []
orangesites = []
greensites = []
nsites = anc_data[ anc_data.keys()[0] ].__len__()
print "\n"
for site in range(1, nsites+1):
    indel = False
    red = False
    orange = False
    green = False
    compare_to_this_state = get_ml_state( anc_data[anc_data.keys()[0]][site] )
    for anc in anc_data.keys():
        this_ml_state = get_ml_state( anc_data[anc][site] )
        # test for indel mismatch:
        if this_ml_state != compare_to_this_state:
            if this_ml_state == "-" or compare_to_this_state == "-":
                indel = True
        # test for red:
        if this_ml_state != compare_to_this_state and indel == False:
            if False == anc_data[anc][site].keys().__contains__(compare_to_this_state):
                red = True
            elif False == anc_data[ anc_data.keys()[0] ][site].keys().__contains__(this_ml_state):
                red = True
            elif anc_data[anc_data.keys()[0]][site][this_ml_state] < 0.05 and anc_data[anc][site][compare_to_this_state] < 0.05:
                red = True
        # test for orange:
        if this_ml_state != compare_to_this_state and red == False and indel == False:
            orange = True
        # test for green
        if this_ml_state == compare_to_this_state:
            if (anc_data[anc][site][this_ml_state] < 0.8 and anc_data[anc_data.keys()[0]][site][this_ml_state] > 0.8) or (anc_data[anc][site][this_ml_state] > 0.8 and anc_data[anc_data.keys()[0]][site][this_ml_state] < 0.8):
                green = True
    if indel:
        countindel += 1
    if red:
        countred += 1
        redsites.append(site)
        print "site", site, " - case 1"
        print "\t",get_ml_state( anc_data[anc_data.keys()[0]][site] ), anc_data[anc_data.keys()[0]][site]
        print "\t",get_ml_state( anc_data[anc_data.keys()[1]][site] ), anc_data[anc_data.keys()[1]][site]
    if orange:
        countorange +=1
        orangesites.append(site)
        print "site", site, " - case 2"
        print "\t", get_ml_state( anc_data[anc_data.keys()[0]][site] ), anc_data[anc_data.keys()[0]][site]
        print "\t",get_ml_state( anc_data[anc_data.keys()[1]][site] ), anc_data[anc_data.keys()[1]][site]
    if green:
        countgreen += 1
        greensites.append(site)
        print "site", site, " - case 3"
        print "\t",get_ml_state( anc_data[anc_data.keys()[0]][site] ), anc_data[anc_data.keys()[0]][site]
        print "\t",get_ml_state( anc_data[anc_data.keys()[1]][site] ), anc_data[anc_data.keys()[1]][site]
print "============================================================"
print "Legend:"
print "- case 1: sites with disagreeing ML states, and neither ancestral vector has strong support for the others' state."
print "- case 2: sites with disagreeing ML states, but one or both of the ancestral vectors has PP < 0.05 for the others' ML state."
print "- case 3: sites with the same ML state, but one of the ancestral vectors strongly supports (PP > 0.8) the state, while the other vector poorly supports (PP < 0.8) the state."
print "============================================================"
print "\n"

print "============================================================"
print "Summary:"
print "nsites=", nsites
print "indel_mismatch=", countindel
print "case 1=", countred, "sites: ",  redsites
print "case 2=", countorange, "sites: ",  orangesites
print "case 3=", countgreen, "sites: ", greensites
print "============================================================"
print "\n\n"

"""
#
# Write a FASTA file with all the sequences
#
fasta_sequences = {} # key = unique name, value = string with sequence
for dat in dat_files:
    fasta_sequences[dat] = ""
    fin = open(dat_files[dat], "r")
    for line in fin.readlines():
        line = line.strip()
        tokens = line.split()
        fasta_sequences[ dat ] += tokens[1]
    fin.close()
fout = open(unaligned_fasta_path, "w")
for f in fasta_sequences:
    fout.write( ">" + f + "\n")
    fout.write( fasta_sequences[f] + "\n" )
fout.close()
"""