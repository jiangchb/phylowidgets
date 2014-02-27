#!/usr/bin/env python
#################################################################
#
# Estimate ASR error.
#
# Victor Hanson-Smith
# victorhs@cs.uoregon.edu
#
# This script estimates reconstructed ancestral sequence error due
# to the phylogenetic model only.
#
# INSTALLATION:
# This script requires that you already installed the following programs:
# 1. Lazarus
# 2. Seq-Gen
# Scroll down in this script, and modify the variable RUNLAZARUS and
# SEQGEN to point to locations that appropriate for your computer.
#
# USAGE:
# The following command-line parameters are required
# --treepath
# --ingroup
# --outgroup
# --modelpath
# --outputdir
# --nreps
# --seqlen
# --seqgennode
#
#################################################################

__author__ = "Victor Hanson-Smith"
__date__ = "08192012"
__usage__ = "estimate_asr_error.py"

import os
import random
import re
import sys
import math
from argParser import *
ap = ArgParser(sys.argv)

#
# Point to necessary software. . .
#
RUNLAZARUS = "python /Users/victor/Documents/SourceCode/Lazarus/lazarus.py"
SEQGEN = "~/Applications/Seq-Gen.v1.3.2/source/seq-gen"

mltreepath = ap.getArg("--treepath")
workspaceDirectory = ap.getArg("--outputdir")
reps = int(ap.getArg("--nreps"))
sequenceLength = int(ap.getArg("--seqlen"))
MODELPATH = ap.getArg("--modelpath")


def floatToString(f):
    return "%.3f"%f

# set is an array of floats
def calculateAverage(set):
    sum = 0.0
    for x in set:
        sum += x
    return sum / float( set.__len__() )

def calculateStandardDeviation(set):
    avg = calculateAverage(set)
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() ) )

def calculateVariance(set):
    avg = calculateAverage(set)
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return float( sumofsquares ) / float( set.__len__() ) 

# returns a "bin" number for this P value
def binForProb(p):
    return int(p / 0.05)

# return the P value of the floor of this bin
def probForBin(b):
    x = float(b*5) / float(100)
    return x
    
    #if x == 1.00:
    #    return x
    
    # Correct the X mark location
    # For example, the bin with X range 0.050 to 0.100 needs to be plotted at 0.075
    #return x #+ 0.025

def getRepName(rep):
    repName = "rep-" + rep.__str__()
    return repName

def buildOutputDirectory():
    if os.path.exists(workspaceDirectory) == False:
        os.system("mkdir " + workspaceDirectory)
    for i in range(1, reps+1):
        repName = getRepName(i)
        if os.path.exists(workspaceDirectory + "/" + repName) == False:
            print ". I'm building the output directory for " + repName
            os.system("mkdir " + workspaceDirectory + "/" + repName)
            if os.path.exists(workspaceDirectory + "/" + repName):
                print "... OK."
            else:
                print "\n. Hmmm, something is wrong.  I was unable to create the directory ", (workspaceDirectory + "/" + repName)
                print ". Perhaps you do not have privileges to write?"
                exit()

def runSeqGen():
    for i in range(1, reps+1):
        repName = getRepName(i)
        if os.path.exists(workspaceDirectory + "/" + repName + "/SEQGEN") == False:
            os.system("mkdir " + workspaceDirectory + "/" + repName + "/SEQGEN")    
    for i in range(1, reps+1):
        repName = getRepName(i)
        treePath = mltreepath
        outputPath = workspaceDirectory + "/" + repName + "/SEQGEN/" + "seqgen_output.txt"
        model = ap.getArg("--modelname")
        rn = random.randint(1,100000000)
        command = SEQGEN + " -z " + rn.__str__() + " -g4 -a2.893 -m" + model + " -n1 -l" + sequenceLength.__str__() + " -wa < " + treePath + " > " + outputPath;
        print command
        os.system(command)

def getExtantFromSeqGen():
    for i in range(1, reps+1):
        repName = getRepName(i)
        if os.path.exists(workspaceDirectory + "/" + repName + "/LAZARUS") == False:
            os.system("mkdir " + workspaceDirectory + "/" + repName + "/LAZARUS")    
            
    for i in range(1, reps+1):
        repName = getRepName(i)
    
        # consult the method named 'runSeqGen' for the output path:
        print "PARSING SEQ-GEN OUTPUT for " + repName
        fsq = open(workspaceDirectory + "/" + repName + "/SEQGEN/" + "seqgen_output.txt", "r")
        falign = open(workspaceDirectory + "/" + repName + "/LAZARUS/alignment.fasta", "w")
        lines = fsq.readlines()
        for l in lines:
            l = l.strip()
            tokens = l.split()
            tokens[0] = re.sub(">", "", tokens[0])
            if tokens[0][0].isalpha():
                print tokens[0]
                falign.write( ">" + tokens[0] + "\n")
                falign.write( tokens[1] + "\n" )
        falign.close()
        fsq.close()

def runLazarus():
    for i in range(1, reps+1):
        repName = getRepName(i)
        dir = workspaceDirectory + "/" + repName + "/LAZARUS"                
        command = RUNLAZARUS + " --codeml --outputdir " + dir + " --verbose 9 --alignment " + dir + "/alignment.fasta --tree " + mltreepath + " --model " + MODELPATH
        print command
        os.system(command)

# given a file mlAnc.node.X.txt, this method returns the filepath to the *.dat file containing
# the posterior distribution of ancestral states.
def get_mlAncFilePath( mlAnc_path ):
    fin = open( mlAnc_path, "r" )
    lines = fin.readlines()
    fin.close()
    tokens = lines[10].strip().split()
    
    tree_num = re.sub("\(", "", tokens[4])
    tree_num = re.sub("\)", "", tree_num)
    tree_num = re.sub("\,", "", tree_num)
    tree_num = re.sub("\#", "", tree_num)
    
    node_num = re.sub("\#", "", tokens[8])

    return "tree" + tree_num + "/node" + node_num + ".dat"

#
# This is a helper method for 'getAncestors'
#
def getAnc(repName, ingroup, outgroup, nodeid):
    dir = workspaceDirectory + "/" + repName + "/LAZARUS"
    command = RUNLAZARUS + " --getanc --ingroup " + ingroup + " --outgroup " + outgroup + " --outputdir " + dir
    #print command
    os.system(command)   
    command = "mv " + dir + "/ancestor.out.txt " + dir + "/mlAnc.node." + nodeid + ".txt"
    #print command
    os.system(command)

def getAncestors():
    for i in range(1, reps+1):
        repName = getRepName(i)
        outgroup = ap.getArg("--outgroup")
        ingroup = ap.getArg("--ingroup")
        getAnc(repName, ingroup, outgroup, "a")

#
# plot PP (binned into 5% groups) versus ASR accuracy
#
#
def plotPPvsACCAll():
    mlppbins = {} # key = bin number, value = an array of pps
    mlccbins = {} # key = bin number, value = integer number of correct inferences within the bin
   
    for i in range(0,21):
        mlppbins[i] = []
        mlccbins[i] = 0
   
    for i in range(1, reps+1):
        repName = getRepName(i) 
        dir = workspaceDirectory + "/" + repName + "/LAZARUS"
        seqgennode = ap.getArg("--seqgennode")
                
        #
        # parse the SEQGEN output
        #
        truesequence = ""
        f = open(workspaceDirectory + "/" + repName + "/SEQGEN/seqgen_output.txt", "r")
        lines = f.readlines()
        for l in lines:
            matchString = "^" + seqgennode.__str__() + "\t"
            if None != re.match(matchString, l):
                l = l.strip()
                tokens = l.split()
                truesequence = tokens[1]
                #print "227:", truesequence
        f.close()
        if os.path.exists(dir + "/ancestor-ml.dat") == False:
            print "ERROR! I cannot find " + dir + "/ancestor-ml.dat"
            exit(1)
            continue
        
        #
        # parse the ML reconstruction
        #
        f = open(dir + "/ancestor-ml.dat" , "r")
        lines = f.readlines()
        site = 0
        for l in lines:
            l = l.strip()
            tokens = l.split()
            thisstate = tokens[1]
            thispp = float( tokens[2] )
                
            mlppbins[ binForProb(thispp)].append( thispp )
            if truesequence[site] == thisstate:
                mlccbins[ binForProb(thispp) ] += 1
            site += 1
        f.close()
    
    #print "mlppbins:", mlppbins
    
    plotdir = workspaceDirectory + "/PLOTS"
    if os.path.exists(plotdir) == False:
        os.system("mkdir " + plotdir)
    
    #
    # Plot the number of points in each bin
    #
    data = {}
    for i in range(0,21):
        data[i] = mlccbins[i] / (float(ap.getArg("--nreps"))*float(ap.getArg("--seqlen")))
        if data[i] == 0.0:
            data[i] = 0.0001
    barplot2(data, "ASR probability", "Proportion of Sites", "binsize")
    
    
    #
    # normalize the counts into P values:
    # remove bins with fewer than 20 data points
    #
    mlpoints = {} # key = x value (mean PP for this bin), value = y value (accuracy of bin)
    averageml = {} # key = bin, value = X average for this bin
    #mlerrorindex = 0.0
    
    for i in range(0, 21):
        if mlppbins[i].__len__() < 5:
            continue # skip bins with less than 20 inferences.
        
        avgmlpp = float(sum(mlppbins[i])) / float(mlppbins[i].__len__())
        fractioncorrect = float(mlccbins[i]) / float(mlppbins[i].__len__())
        mlpoints[avgmlpp] = fractioncorrect        
        averageml[i] = avgmlpp

            
    #print "mlpoints:", mlpoints
    #print "averageml:", averageml

    gerr = 0.0 # geometric error from ideal x=y data
    for i in averageml:
        actual = mlpoints[averageml[i]]
        ideal = probForBin(i) + 0.025
        if i == 21:
            ideal = 1.0
        gerr += abs(actual - ideal)
    error = float(gerr/averageml.keys().__len__())
    print "geo. error = ", error   
    plotdir = workspaceDirectory + "/PLOTS"
    f = open(plotdir + "/summary.pp.txt", "w")
    f.write("geometric error of PP from y=x: " + error.__str__() + "\n")
    f.close()


    # create a CRAN script:
    f = open(plotdir + "/ppacc.all.cran", "w")

    # add ML points
    xvals = mlpoints.keys()
    xvals.sort()
    x = "mlx <- c("
    for k in xvals:
        x += floatToString(k) + ","
    x = re.sub(",$", "", x)
    x += ")"
    f.write( x + "\n")
    mly = "mly <- c("
    for k in xvals:
        mly += floatToString(mlpoints[k])  + ","
    mly = re.sub(",$", "", mly)
    mly += ")"
    f.write( mly + "\n")


    f.write("pdf(\"" + plotdir + "/ppacc.all.pdf" + "\", height=4, width=4) \n")
    f.write("plot(mlx, mly, xlab=\"ASR probability\", ylab=\"Proportion of Correct Inferences\", xlim=c(0.0, 1.0), ylim=c(0.0, 1.0), pch=20,cex=2, col=\"red\")\n")  
    f.write("abline(a=0.0, b=1)\n")
    f.write("dev.off()\n")
    f.close()
    
    os.system("r --no-save < " + plotdir + "/ppacc.all.cran > " + plotdir + "/ppacc.all.cran.log") 

# returns:
# bars: bars[data.keys()] = mean of data[i]
# meanvals: meanvals = [bars[0], bars[1], ...]
def prepare_data_for_barplot2( data ):
    lengths = {}
    bars = {}
    meanvals = []
    semupper = []
    semlower = []    
            
    keys = data.keys()
    keys.sort()        
    
    # init data structs
    for i in keys:
        bars[i] = {}
        bars[i] = 0.0
        lengths[i] = []
    #
    for i in keys:
        #print "data[%d]"%i, data[i]
        # sum all the tree length errors for this pointset
        bars[i] = mean( data[i] )
        meanvals.append( bars[i] )
        semupper.append( bars[i] + stderr( data[i] ) )     
        semlower.append( bars[i] - stderr( data[i] ) )            
    return [ bars, meanvals, semupper, semlower ]

#
# data[xgroup][series] = value
#
def barplot2(data, xlab, ylab, filekeyword):
        
    pointsets = data.keys()
    pointsets.sort()
    finalset = pointsets[ pointsets.__len__()-1 ]
    
    plotdir = workspaceDirectory + "/PLOTS"
    tablepath = plotdir + "/barplot2.txt"
    fout = open(tablepath, "w")
    for p in pointsets:
        if p != finalset:
            fout.write( data[p].__str__() + "\t")
        else:
            fout.write( data[p].__str__() )            
    fout.write("\n")
    fout.close()
    pdfpath = plotdir + "/barplot2.pdf"
    cranstr = "pdf(\"" + pdfpath + "\", height=4, width=4)\n"    
    cranstr += "bars <- read.table(\"" + tablepath + "\", header=F, sep=\"\\t\")\n"
    cranstr += "barx = barplot(as.matrix(bars),axes = FALSE, axisnames = FALSE, beside=TRUE, col=c(\"red\"), xlab=\"" + xlab + "\", ylab=\"" + ylab + "\");\n"
    
    xlabels = "labels <- c("
    for p in pointsets:
        if p == 20:
            xlabels += "\"" + probForBin(p).__str__() + "\","
        else:
            xlabels += "\"" + probForBin(p).__str__() + " - " + (probForBin(p)+0.05).__str__() + "\","
    xlabels = re.sub(",$", "", xlabels)
    #cranstr += "labels <- paste(\"This is bar #\", 1:16, sep ='');\n"
    cranstr += xlabels + ");\n"
    cranstr += "text(barx, par(\"usr\")[3], labels = labels, srt = 45, adj = 1.15, cex = 0.5,xpd = TRUE);\n"
    cranstr += "axis(2);\n"
    
    cranpath = plotdir + "/barplot2.cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    os.system("r --no-save < " + cranpath + " > " + plotdir + "/barplot2.cran.log")


#
# report the number of inferences which are strongly-supported and incorrect
#
def reportError():
    seqgennodes = [ap.getArg("--seqgennode")]
    nodes = ["a"]

    rep_mlaccuracy = {} # key = rep number, value = proportion correct for all nodes
   
    for i in range(1, reps+1):
        print ". analyzing replicate ", i
        repName = getRepName(i) 
        dir = workspaceDirectory + "/" + repName + "/LAZARUS"

        rep_mlcount = 0 # the number of ML inferences for this rep
        rep_mlcorrect = 0 # the number of ML *correct* inferences for this rep
        
        for j in range(0, seqgennodes.__len__()):
            seqgennode = seqgennodes[j]
            node = nodes[j]
        
            #
            # parse the SEQGEN output
            #
            truesequence = ""
            f = open(workspaceDirectory + "/" + repName + "/SEQGEN/seqgen_output.txt", "r")
            #print "seqgen out at: " + workspaceDirectory + "/" + repName + "/SEQGEN/seqgen_output.txt", "r"
            lines = f.readlines()
            for l in lines:
                matchString = "^" + seqgennode.__str__() + "\t"
                if None != re.match(matchString, l):
                    l = l.strip()
                    tokens = l.split()
                    truesequence = tokens[1]
            f.close()
            
            #
            # parse the ML reconstruction
            #
            datpath = get_mlAncFilePath( dir + "/mlAnc.node." + node.__str__() + ".txt" )
            f = open(dir + "/" + datpath , "r")
            
            lines = f.readlines()
            site = 0
            for l in lines:
                l = l.strip()
                tokens = l.split()
                thisstate = tokens[1]
                thispp = float( tokens[2] )
                
                rep_mlcount += 1
                if truesequence[site] == thisstate:
                    rep_mlcorrect += 1
                site += 1
            f.close()

        
        # sum stats for this node:
        ml_accuracy = float(rep_mlcorrect) / float(rep_mlcount)
        rep_mlaccuracy[i] = ml_accuracy
        
    for i in rep_mlaccuracy:
        print "rep", i, "ML accuracy= %.2f"%rep_mlaccuracy[i]
    
    plotdir = workspaceDirectory + "/PLOTS"
    f = open(plotdir + "/summary.txt", "w")
    f.write("mean proportion of sites correctly inferred = " + calculateAverage(rep_mlaccuracy.values()).__str__() + "\n" )
    f.write("s.d. proportion of sites correctly inferred = " + calculateStandardDeviation(rep_mlaccuracy.values()).__str__() + "\n" )       
    f.close()


########################################
#
# main()...
#
# 1. build a directory for each replicate
# 2. use Seq-Gen to simulate replicates
# 3. gather the descendant seqeunces from Seq-Gen
# 4. invoke Lazarus on each replicate
# 5. plot the results
#
########################################
if False == ap.getOptionalToggle("--skip_asr"):
    buildOutputDirectory()
    runSeqGen()
    getExtantFromSeqGen()
    runLazarus()
print "\n. Retrieving ancestors. . ."
getAncestors()
plotPPvsACCAll()
reportError()
