#
# Usage:
# python get_genbank_data --inpath <filename1> --outpath <filename2> --accids --speciesnames --dbname <GenBank db name> --dbaction <fetch/search>
#
# Optional Arg:
#     --skip_query True : skip the Genbank query and assume all the *.gbk files already exist.
#
#     --accid_for_name True : use accession IDs for the taxa names in the output
#
#     --taxonmy_for_name N : use the first N elements of the taxonomy for the output name
#
#     --skip_fasta True : don't write the sequences to disk.
#
#     --translation_key True : write a text file with accids and species names
#
# where <filename1> is a text file with one GenBank accession ID per line.
# This script will query Genbank for each acc-ID, and grab the first search result.
# The GenBank files will be written to the current directory as *.gbk files.
# The sequences from all the GenBank files will be collected into a non-aligned FASTA
# file at <filename2>.
#
# if --skip_query True, then this script will skip the query and assume that several *.gbk files
#     already exist in the local directory.
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
from urllib2 import HTTPError

from argParser import *
ap = ArgParser(sys.argv)

BIOPYTHON_DIRECTORY = "/Users/victor/Applications/MacBioPython100a4" # on Victor's laptop
sys.path.append(BIOPYTHON_DIRECTORY)

#from sshDispatcher import *
from Bio import Entrez, SeqIO

#
# Returns a list of accession ids or species names, parsed from the file at --inpath
#
def Fill_helper():
    ncbi_list = []
    fin = open( ap.getArg("--inpath"), "r")
    last_line = ""
    for line in fin.readlines():
        if True == line.__contains__(">") or True == line.__contains__("_"):
            ncbi_list.append(last_line)
        else:
            last_line += line
    ncbi_list.append(last_line)        
    return ncbi_list    

def Fill_ncbi_list():
    if False != ap.getOptionalArg("--sequences"):
        return Fill_helper()
    ncbi_list = []
    fin = open( ap.getArg("--inpath"), "r")
    for line in fin.readlines():
        line = line.strip()
        
        if False != ap.getOptionalArg("--accids"): 
            if line.__len__() > 2:
                tokens = line.split()
                for t in tokens:
                    ncbi_list.append( t )
        elif False != ap.getOptionalArg("--speciesnames"):
            if line.__len__() > 2:
                line = re.sub("\ ", "_", line)
                ncbi_list.append(line)                
    fin.close()
    print ncbi_list
    return ncbi_list


# input: an array of NCBI accession numbers for full genomes
# post-condition: this directory will contain a collection of *.gbk files, one for each NCBI ID.
def fetchGenbankData(ncbi_list):
    Entrez.email = "victorhs@cs.uoregon.edu"
    try:
        for ncbi in ncbi_list:  
            print "fetching genbank data for acc. ID: " + ncbi.__str__()
            retry_count = 0
            handle = None
            # fetch the GenBank record; retry up to 3 times if the connection is problematic.
            while retry_count < 3:
                try:
                    dbname = ap.getOptionalArg("--dbname")
                    if dbname == False:
                        dbname = "nucleotide"
                    dbaction = ap.getOptionalArg("--dbaction")
                    if dbaction == "search":
                        pre_query = re.sub("\_", "\ ", ncbi)
                        pre_handle = Entrez.esearch(db=dbname, retmax=1, term=pre_query)
                        pre_record = Entrez.read(pre_handle)
                        if pre_record["IdList"].__len__() < 1:
                            print "skipping this ID. . .\n"
                            retry_count = 3
                            break
                        print "ID= ", pre_record["IdList"][0].__str__()
                        handle = Entrez.efetch(db=dbname, id = pre_record["IdList"][0].__str__(), rettype="gb")
                    elif dbaction == "fetch":
                        handle = Entrez.efetch(db=dbname, id = ncbi.__str__(), rettype="gb")

                    break # if we get the handle OK, then break out of the loop
                except ValueError:
                    time.sleep(1)
                    print "Something went wrong, my GenBank query for NCBI #" + ncbi.__str__() + " returned no records."
                    print "I'm trying again. . ."
                    retry_count += 1
                except HTTPError:
                    time.sleep(1)
                    print "HTTP error with my GenBank query for NCBI #" + ncbi.__str__() + " returned no records."
                    print "I'm trying again. . ."
                    retry_count += 1                    
            if retry_count < 3 and handle != None:
            # write the GBK file:
                print ". . . ok, found", ncbi
                f = open(ncbi.__str__() + ".gbk", "w")
                f.write( handle.read() )
                f.close()
                
                f = open(ncbi.__str__() + ".gbk", "r")
                for l in f.readlines():
                    if l.__contains__("Cannot connect to the server"):
                        print "WARNING: The GenBank file for " + ncbi.__str__() + " might be incomplete."
                f.close()
            
            handle.close()
            time.sleep(1)
    except ValueError:
        print "Something went wrong, my GenBank query for NCBI #" + ncbi.__str__() + " returned no records."
        print "I'm not going retry anymore.  Sorry."
        exit(1)

#
# Returns an array [names, seqs]
# where names is a hashtable with keys = accession IDS, values = species name read from the GBK file
# and seqs is a hashtable with keys = accession IDS, values = sequences read from the GBK file
#
def Read_GBK_Files(ncbi_list):
    ncbi_to_name = {} 
    ncbi_to_seq = {}
    for ncbi in ncbi_list:
        print "Opening", ncbi.__str__() + ".gbk"
        try:
            handle = open(ncbi.__str__() + ".gbk", "r")
            record = SeqIO.read(handle, "genbank")
            #continue
            if False == ap.getOptionalArg("--taxonomy_for_name"):
                ncbi_to_name[ ncbi ] = ncbi
            else:
                n = int( ap.getOptionalArg("--taxonomy_for_name") )
                taxonomy = record.annotations["taxonomy"]
                name = ""
                for i in range(0,n):
                    name += taxonomy[i] + "."
                name = name[0:name.__len__()-1]
                ncbi_to_name[ ncbi ] = name
            ncbi_to_seq[ ncbi ] = record.seq.__str__()
            #handle.close()
        except ValueError:
            print "Hmmm, I had a problem with", ncbi
    return [ncbi_to_name, ncbi_to_seq]

#
# pre-condition: a collection of *.gbk files exists in the DATA_DIRECTORY
# post-condition: an unaligned FASTA file will exist for each gene, and a file named gene_list.txt will be written.
#
def Genbank_Genes_ToFasta(ncbi_list):
    [ncbi_to_name, ncbi_to_seq] = Read_GBK_Files(ncbi_list)
    f = open(ap.getArg("--outpath"), "w")
    for ncbi in ncbi_list:
        if False != ap.getOptionalArg("--accid_for_name"):
            f.write(">" + ncbi + "\n")
        else:
            f.write(">" + ncbi_to_name[ncbi] + "\n")
        f.write(ncbi_to_seq[ncbi] + "\n")
    f.close()

#
#
#
def Print_Translation_Key(ncbi_list):
    [ncbi_to_name, ncbi_to_seq] = Read_GBK_Files(ncbi_list)
    for ncbi in ncbi_list:
        name = ncbi_to_name[ncbi]
        tokens = name.split("_")
        shortname = tokens[0]
        if tokens.__len__() > 1 and False == tokens[1].startswith("sp."):
            shortname += "." + tokens[1]
        shortname = re.sub("\(", "", shortname)
        shortname = re.sub("\)", "", shortname)
        print ncbi + "   " + shortname+ "-" + ncbi

###########################################
# main
#
ncbi_list = Fill_ncbi_list()
if False == ap.getOptionalArg("--skip_query"):
    fetchGenbankData(ncbi_list)
if False == ap.getOptionalArg("--skip_fasta"):
    Genbank_Genes_ToFasta(ncbi_list)
if False != ap.getOptionalArg("--translation_key"):
    Print_Translation_Key(ncbi_list)