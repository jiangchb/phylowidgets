#
# Usage:
# python blast_multi.py --phylip_path <path>
#
# Input: a sequential PHYLIP alignment.
#
# This script will reverse-BLAST every sequence, select the top-hit for each BLAST search,
# and then report the official name of the sequence.
#
# This is useful, for example, if someone gave you a sequence alignment with poorly-annotated names
# and you would like to know more about the provenance of the data.
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
from Bio.Blast import NCBIWWW
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML

#
#
#
def fetchGenbankData(seq_list):
    Entrez.email = "victorhs@cs.uoregon.edu"
    try:
        for taxa in seq_list.keys():
            seq = seq_list[taxa]  
            print "BLAST-ing NCBI for sequence ID: " + taxa.__str__()
            retry_count = 0
            # fetch the GenBank record; retry up to 3 times if the connection is problematic.
            while retry_count < 3:
                try:
                    blast_handle = NCBIWWW.qblast('blastp', 'nr', seq)
    
                    blast_handle.seek(0)
                    blast_file = open( taxa.__str__() + '.xml', 'w' )
                    blast_file.write( blast_handle.read() )
                    blast_file.close()
                    blast_handle.close()
                    print ". . . results written to " + taxa.__str__() + '.xml'
                    break # if we get the handle OK, then break out of the loop
                except ValueError:
                    sleep(3)
                    print "Something went wrong, my GenBank query for taxa " + taxa.__str__() + " returned no records."
                    print "I'm trying again. . ."
                    retry_count += 1
            time.sleep(2)
    except ValueError:
        print "Something went wrong, my GenBank query for taxa " + taxa.__str__() + " returned no records."
        print "I'm not going retry anymore.  Sorry."
        exit(1)
        
        
def Genbank_Genes_ToFasta(ncbi_list):
    ncbi_to_name = {} 
    ncbi_to_seq = {}
    
    for taxa in seq_list:
        seq = seq_list[taxa]  
        file_handle = open(taxa.__str__() + '.xml', "r")
        blast_records = NCBIXML.read(file_handle)
        
        max_score = 0
        max_query = []
        max_name = None
        max_acc = None
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                if hsp.score > max_score:
                    max_score = hsp.score
                    max_match = hsp.match
                    max_name = alignment.title
                    max_name = re.sub(".*\[", "", max_name)
                    max_name = re.sub("\].*", "", max_name)
                    max_query = hsp.query
        print taxa, max_name
        file_handle.close()



seq_list = {} # key = name, value = sequence
fin = open( ap.getArg("--phylip_path"), "r")
for line in fin.readlines()[1:]:
    line = line.strip()
    if line.__len__() > 2:
        tokens = line.split()
        seq_list[ tokens[0] ] = tokens[1]
        #print tokens[0], tokens[1]
fin.close()

#if True == ap.getOptionalArg("--phylip_path"):
fetchGenbankData(seq_list)
Genbank_Genes_ToFasta(seq_list)