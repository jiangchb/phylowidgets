################################
#
# input: a Newick-formatted tree
#
# method: bins all branch lengths into 0.05-sized bins and counts
# the proportion of BLs in each bin
#
# output: a CRAN script and a PDF plot (generated from the CRAN script)
#
#

import math
import os
import re
import sys
from dendropy import Tree

from argParser import ArgParser
argParser = ArgParser(sys.argv)

inpath = sys.argv[1]
id = sys.argv[2]

import numpy as np
import cogent.maths.stats.test as stats

colors = ["orangered","royalblue", "black"]
pch = [20,5,3]
COLORS = {}
COLORS["unimax"] = "tomato"
COLORS["multimax"] = "deepskyblue"
COLORS["true"] = "black"
COLORS["bionj"] = "snow3"

# set is an array of floats
def mean(set):
    if set.__len__() == 0:
        return None
    sum = 0.0
    for x in set:
        sum += x
    return sum / float( set.__len__() )

# standard deviation
def sd(set):
    avg = mean(set)
    if avg == None:
        return None
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() ) )

# calculates variance
def var(set):
    avg = mean(set)
    if avg == None:
        return None
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() - 1 ) ) 

def stderr(set):
    return (sd(set) / math.sqrt( set.__len__() ) )

def get_bls(inpath):
    bls = []
    fin = open(inpath, "r")
    for l in fin.readlines():
        if l.__len__() > 1:
            tokens = l.split(":")
            for t in tokens[1:]:
                bl = t.split(")")[0]
                bl = bl.split(",")[0]
                bl = float(bl)
                bls.append( bl )
    fin.close()
    return bls

def stats_about_bls(path, bls):
    # calculate mean
    sum = 0.0
    for b in bls:
        sum += b
    mean = sum / bls.__len__()
    mean_str = "%.3f"%mean
    
    # median
    bls.sort()
    median = bls[len(bls)/2]
    median_str = "%.3f"%median
    
    print path, "mean=", mean_str, "median=", median_str, "sem=", stderr(bls)

bls = get_bls(inpath)
#
# returns a hashtable, key = bin numbers, rounded to the nearest 0.01, value = proportion of total branches in that bin
#
minbl = bls[0]
for b in bls:
    if b < minbl:
        minbl = b
maxbl = bls[0]
for b in bls:
    if b > maxbl:
        maxbl = b

#MAX_BIN_LENGTH = float("%.5f"%maxbl)
#MIN_BIN_LENGTH = float("%.5f"%minbl)
#BL_BIN_SIZE = (MAX_BIN_LENGTH - MIN_BIN_LENGTH) / 20
#MAX_BIN_LENGTH += BL_BIN_SIZE
#MIN_BIN_LENGTH -= BL_BIN_SIZE

MIN_BIN_LENGTH = 0.100
MAX_BIN_LENGTH = 0.165
BL_BIN_SIZE =    0.0025

print MIN_BIN_LENGTH, MAX_BIN_LENGTH, BL_BIN_SIZE

def get_bin_count():
    return int( (MAX_BIN_LENGTH-MIN_BIN_LENGTH)/BL_BIN_SIZE )


def calculate_bins(bls):    
    max_bin = None
    bins = {} # key = floor, value = count
    
    i = MIN_BIN_LENGTH
    while(i < MAX_BIN_LENGTH):
        bins[i] = 0.0
        i += BL_BIN_SIZE
        print i
    
    for b in bls:
        if b < MIN_BIN_LENGTH or b > MAX_BIN_LENGTH:
            continue
        
        bin = MIN_BIN_LENGTH
        while (bin <= b):
            bin += BL_BIN_SIZE
        bins[bin - BL_BIN_SIZE] += 1
    
    print bins
    
    normalized_bins = {}
    for b in bins.keys():
        normalized_bins[b] = float(bins[b]) / bls.__len__()
    return normalized_bins

#
# data[xgroup][series] = value
#
def barplot1(data, xlab, ylab, filekeyword):
    
    pointsets = data.keys()
    pointsets.sort()
    
    finalset = pointsets[ pointsets.__len__()-1 ]
    tablepath = "barplot.table." + filekeyword + ".txt"
    fout = open(tablepath, "w")
    for p in pointsets:
        if p != finalset:
            fout.write(p.__str__() + "\t")
        else:
            fout.write(p.__str__() )
    fout.write("\n")
    for p in pointsets:
        if p != finalset:
            fout.write( data[p].__str__() + "\t")
        else:
            fout.write( data[p].__str__() )            
    fout.write("\n")
    fout.close()
    
            
    pdfpath = "barplot." + filekeyword + id.__str__() + ".pdf"
    cranstr = "pdf(\"" + pdfpath + "\", width=8, height=4);\n"    
    cranstr += "bars <- read.table(\"" + tablepath + "\", header=T, sep=\"\\t\")\n"
    
    print pointsets


    cranstr += "pointsets <- c("
    for p in pointsets:
        cranstr += (p).__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"
    
    cranstr += "barx = barplot(as.matrix(bars), beside=TRUE, col=c(\"" + COLORS["unimax"] + "\", \"" + COLORS["multimax"] + "\"), ylim=range(0,1.0), names.arg=pointsets);\n"
    cranpath = "barplot." + filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save < " + cranpath)

#
# Fill the bins with data
#
bins = calculate_bins( bls )
barplot1(bins, "BL bins", "proportion", "bl_distribution" + id.__str__())
print bls
print bins

