#
# clean_trees.py
# 
# Victor Hanson-Smith
# victorhansonsmith@gmail.com
#
# This script will remove branch support values (e.g. bootstraps) from a Newick-formatted tree.
# This is useful if you're using PAML, which seems to be unable to deal with
# branch labels.
#
# USAGE: 
#     python clean_trees.py FILEIN FILEOUT
#
# . . .where FILEIN is a text file containing one or more Newick-formatted trees,
# with one tree on each line.  FILEOUT is the desired output filename.  
# After the script finished, FILEOUT will contain the support-less versions 
# of each tree in FILEIN. If FILEOUT does not exist, then it will be created de novo.
#

import re,sys,os
from argParser import *
ap = ArgParser(sys.argv)

def strip_support(line):
    line = re.sub("\)\d+\.*\d+\:", "):", line)
    return line

def clean_mixtrees(fin):
    #cleantrees = []
    clean = ""
    for line in fin.xreadlines():
        line = re.sub("\[\d+\.\d+\,\d+\.*\d+\]", "", line)
        line = re.sub("\:", "", line)
        #line = re.sub("\]", "", line)        
        clean += line
    return clean

def clean_singletree(fin):
    clean = ""
    for line in fin.xreadlines():
        if ap.getOptionalArg("--strip_support"):
            line = strip_support(line)
        line = re.sub("\[", "", line)
        line = re.sub("\]", "", line)
        clean += line
    return clean


fin = open(sys.argv[1], "r")
if ap.getOptionalArg("--split_mixtrees"):
    cleantree = clean_mixtrees(fin)
else:
    cleantree = clean_singletree(fin)
fin.close()

fout = open(sys.argv[2], "w")
fout.write(cleantree)
fout.close()