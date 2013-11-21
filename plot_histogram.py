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

from argParser import ArgParser
argParser = ArgParser(sys.argv)

inpath = sys.argv[1]

import cogent.maths.stats.test as stats

MIN_BIN_LENGTH = 0.0
MAX_BIN_LENGTH = 1.0
BIN_SIZE =    1.0


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


def read_data(inpath):
    fin = open(inpath, "r")

    lines = fin.readlines()
    nsets = lines[0].split().__len__()
    
    set_data = {}
    for i in range(0, nsets):
        set_data[i] = []
    
    for l in lines:
        if l.__len__() > 5:
            tokens = l.split()
            for i in range(0, nsets):
                set_data[i].append( float(tokens[i]) )
    return set_data

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

data = read_data(inpath)

#
# returns a hashtable, key = bin numbers, rounded to the nearest 0.01, value = proportion of total branches in that bin
#
maxval = None
minval = None
for set in data:
    if maxval == None:
        maxval = max(data[set])
    if max(data[set]) > maxval:
        maxval = max(data[set])
    if minval == None:
        minval = min(data[set])
    if min(data[set]) < minval:
        minval = min(data[set])

MAX_BIN_LENGTH = float("%.5f"%maxval)
MIN_BIN_LENGTH = float("%.5f"%minval)
BIN_SIZE = (MAX_BIN_LENGTH - MIN_BIN_LENGTH) / 20
MAX_BIN_LENGTH += BIN_SIZE
MIN_BIN_LENGTH -= BIN_SIZE

print MIN_BIN_LENGTH, MAX_BIN_LENGTH, BIN_SIZE

def get_bin_count():
    return int( (MAX_BIN_LENGTH-MIN_BIN_LENGTH)/BIN_SIZE )


def calculate_bins(bls):    
    bins = {} # key = floor, value = count
    bins[MAX_BIN_LENGTH] = 0.0    
    
    for b in bls:        
        if b > MAX_BIN_LENGTH:
            bins[MAX_BIN_LENGTH] += 1
            print bin, b
        else:    
            bin = MIN_BIN_LENGTH
            while (bin <= b):
                if bin not in bins:
                    bins[bin] = 0.0
                if bin+BIN_SIZE > b and bin <= b:
                    bins[bin] += 1
                    print "130", bin, b
                bin += BIN_SIZE
                
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
    
    print "153"
    
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

    print "pointsets=", pointsets

    cranstr += "pointsets <- c("
    for p in pointsets:
        cranstr += (p).__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"
    
    cranstr += "barx = barplot(as.matrix(bars), beside=TRUE, ylim=range(0,1.0), names.arg=pointsets);\n"
    cranpath = "barplot." + filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save < " + cranpath)

#
# Fill the bins with data
#
for set in data:
    bins = calculate_bins( data[set] )
    print "bins=", bins
    barplot1(bins, "bins", "proportion", "histogram" + set.__str__())