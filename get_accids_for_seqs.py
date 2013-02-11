

#
# --inpath
# --dbname
# --outfile
#

import cgi, os
import cgitb
import Cookie
import sqlite3
import math
import os
import pickle
import Queue
import random
import re
import string
import sys
import time

from argParser import *
ap = ArgParser(sys.argv)

BIOPYTHON_DIRECTORY = "/Users/victor/Applications/MacBioPython100a4" # on Victor's laptop
sys.path.append(BIOPYTHON_DIRECTORY)

#from sshDispatcher import *
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

taxa_seq = {}
taxa_to_accid = {}

#
# Parse a PHYLIP file
#
def get_Sequences():
    fin = open( ap.getArg("--inpath"), "r")
    for line in fin.readlines()[1:]:
        print line
        tokens = line.split()
        taxa = tokens[0]
        seq = re.sub("-", "", tokens[1])
        print seq
        taxa_seq[taxa] = seq

# input: an array of NCBI accession numbers for full genomes
# post-condition: this directory will contain a collection of *.gbk files, one for each NCBI ID.
def fetchGenbankData():    
    Entrez.email = "victorhs@cs.uoregon.edu"
    count = 0
    try:
        for taxa in taxa_seq:
            seq = taxa_seq[taxa]
            retry_count = 0
            # fetch the GenBank record; retry up to 3 times if the connection is problematic.
            while retry_count < 3:
                try:
                    qpath = taxa.__str__() + ".gbk"
                    if False == os.path.exists(qpath):
                        print "Searching Genbank for sequence: ", taxa, seq
                        dbname = ap.getOptionalArg("--dbname")
                        if dbname == False:
                            dbname = "nucleotide"
                        print "."
                        handle = NCBIWWW.qblast("blastp", "nr", sequence=seq, hitlist_size=1)
                        print "|"

                        f = open(qpath, "w")
                        f.write( handle.read() )
                        f.close()
                        handle.close()
                    
                    accid = None
                    f = open( qpath, "r")
                    for l in f.readlines():
                        if l.__contains__("Hit_accession"):
                            l = re.sub("\<", ">", l)
                            tokens = l.split(">")
                            accid = tokens[2]
                            print "Found accid", accid, "for taxa", taxa
                    f.close()

                    taxa_to_accid[ taxa ] = accid
                    retry_count = 3
                except ValueError:
                    time.sleep(3)
                    print "Something went wrong, my GenBank query for NCBI #" + seq.__str__() + " returned no records."
                    print "I'm trying again. . ."
                    retry_count += 1
    except ValueError:
        print "Something went wrong, my GenBank query for NCBI #" + taxa.__str__() + " returned no records."
        print "I'm not going retry anymore.  Sorry."
        exit(1) 

get_Sequences()
fetchGenbankData()

#
# Now write the accids and taxa names to a FASTA file
#
fout = open( ap.getArg("--outpath"), "w" )
for taxa in taxa_seq:
    this_taxa = re.sub("\_", ".", taxa)
    this_accid = re.sub("\_", ".", taxa_to_accid[taxa])
    fout.write(">" + taxa + "." + this_accid + "\n")
    fout.write(taxa_seq[taxa] + "\n")
    print "Added sequence for", this_taxa, this_accid
fout.close()