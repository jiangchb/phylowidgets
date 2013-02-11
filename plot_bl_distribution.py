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

t1path = sys.argv[1]
t2path = sys.argv[2]
id = argParser.getOptionalArg("--id")
if id == False:
    id = ""

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

def get_bls(tree_path):
    # clean the tree of any support values, so we're left only with BLs
    bls = []
    t = Tree()
    t.read_from_path( tree_path, "newick" )
    
    i = t.level_order_edge_iter()
    while True:
        try:
            e = i.next() # in Python 2.x
            len = e.length
            if len != None:
                bls.append( len )
        except StopIteration:
            break
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

filepath_bls = {}
paths = []
for path in [t1path, t2path]:
    these_bls = get_bls(path)
    if these_bls.__len__() > 0:
        paths.append(path)
        filepath_bls[path] = these_bls
#
# returns a hashtable, key = bin numbers, rounded to the nearest 0.01, value = proportion of total branches in that bin
#
MAX_BIN_LENGTH = 0.3
BL_BIN_SIZE = 0.01

def get_bin_count():
    return int( MAX_BIN_LENGTH/BL_BIN_SIZE )

def get_bin_for_length(bl):
    return  (int((bl / BL_BIN_SIZE) % get_bin_count() ) )

def get_length_for_bin(b):
    return (int(b) * BL_BIN_SIZE)

def calculate_bins(bls):    
    max_bin = None
    bins = {}
    for b in bls:
        rounded_bl = get_bin_for_length( b )
        if max_bin == None:
            max_bin = rounded_bl
        elif max_bin < rounded_bl:
            max_bin = rounded_bl
        if bins.__contains__( rounded_bl ):
            bins[ rounded_bl ] += 1
        else:
            bins[ rounded_bl ] = 1
    # zero-out any empty bins
    for j in range(0, get_bin_count()):
        if False == bins.__contains__(j):
            bins[j] = 0
    print bins
    
    #normalize to the number of branches
    #max = 0
    #for b in bins.keys():
    #    if bins[b] > max:
    #        max = bins[b]
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
            fout.write(get_length_for_bin(p).__str__() + "\t")
        else:
            fout.write(get_length_for_bin(p).__str__() )
    fout.write("\n")
    for path in paths:
        for p in pointsets:
            if p != finalset:
                fout.write( data[p][path].__str__() + "\t")
            else:
                fout.write( data[p][path].__str__() )            
        fout.write("\n")
    fout.close()
    
    pdfpath = "barplot." + filekeyword + id.__str__() + ".pdf"
    cranstr = "pdf(\"" + pdfpath + "\", width=8, height=4);\n"    
    cranstr += "bars <- read.table(\"" + tablepath + "\", header=T, sep=\"\\t\")\n"
    
    cranstr += "pointsets <- c("
    for p in pointsets:
        cranstr += get_length_for_bin(p).__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"
    
    cranstr += "barx = barplot(as.matrix(bars), beside=TRUE, col=c(\"" + COLORS["unimax"] + "\", \"" + COLORS["multimax"] + "\"), legend = c(\"Unimax ML\", \"Multimax ML\", ylim=range(0.0, 1.0)), names.arg=pointsets);\n"
    cranpath = "barplot." + filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save < " + cranpath)

#
# Fill the bins with data
#
bins = {}
for path in paths:
    bins[path] = calculate_bins( filepath_bls[path] ) 

"""
# print CRAN script
scriptpath = "plot.bl_distribution.cran"
fout = open(scriptpath, "w")
varcount = 0
for path in paths:
    normalized_bins = bins[path]
    lstr = "lengths" + varcount.__str__() + " <- c("
    keys = normalized_bins.keys()
    keys.sort()
    for x in keys:
        this_length = "%.2f"%x
        lstr += this_length + ","
    lstr = re.sub(",$", "", lstr)
    lstr += ")"
    fout.write(lstr + "\n")
    
    pstr = "proportions" + varcount.__str__() + " <- c("
    for x in keys:
        this_prop = "%.2f"%normalized_bins[x]
        pstr += this_prop + ","
    pstr = re.sub(",$", "", pstr)
    pstr += ")"
    fout.write(pstr + "\n")
    varcount += 1
fout.write("pdf(\"plot.bl_distribution.pdf" + "\", width=8, height=4)\n")
fout.write("plot(lengths0, proportions0, type='l',xlab=\"BLs, binned\", ylab=\"proportion\", col=\"" + colors[0].__str__() + "\", lwd='2', pch=" + pch[0].__str__() + ", main=\"BL Distribution\");\n")
for i in range(1, varcount):
    fout.write("points(lengths" + i.__str__() + ", proportions" + i.__str__() + ", type='l', col=\"" + colors[i].__str__() + "\", lwd='2', pch=" + pch[i].__str__() + ")\n")
fout.write("dev.off()\n")
fout.close()
os.system("r --no-save < " + scriptpath)
"""

#
# barplot
#
dataseries = {}
for binid in range(0, get_bin_count()):
    dataseries[binid] = {}
    for path in paths:
        if bins[path].__contains__( binid ):
            dataseries[binid][path] = bins[path][binid]
        else:
            dataseries[binid][path] = 0.0
barplot1(dataseries, "BL bins", "proportion", "bl_distribution" + id.__str__())

# print stats
print "\n.\n. Stats about these ML branch length distributions:\n."
for path in paths:
    stats_about_bls( path, filepath_bls[path])
[t, p] = stats.t_two_sample(filepath_bls[ paths[0] ], filepath_bls[ paths[1] ])
print "T = ", t, "P=", p
