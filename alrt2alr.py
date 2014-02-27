#
# alrt2.alr.py
# 
# Written by Victor Hanson-Smith
# Email: victorhansonsmith@gmail.com
#
# This script will convert approximate likelihood ratio test statistic values (aLRTs)
# into approximate likelihood ratios (aLRs).
#
# aLRTs can be computed by PhyML, using the command-line option "-b -1"
#
# USING THIS SCRIPT
# %> python alrt2alr.py FILEIN > FILEOUT
#
# . . .where FILEIN is a text file containing one or more Newick-formatted trees,
# with one tree on each line, and with each tree containing aLRT values on branches.  
# FILEOUT is the desired output filename that will contain the aLR trees. 
#
import math,re,sys,os
fin = open(sys.argv[1], "r")
newline = ""
for line in fin.readlines():
    if line.__len__() > 10:
        tokens = line.split(":")
        for t in tokens:
            #print t
            if t.__contains__(")") and False == t.__contains__(";"):
                #print t
                ts = t.split(")")
                newline += ts[0] + ")"
                alrt = float(ts[1])
                if alrt > 1418:
                    alrt = 1418.0
                alr = math.exp(alrt/2.0)
                #print alrt, alr
                newline += "%.4e"%alr
                newline += ":"
            else:
                newline += t
                newline += ":"
        print newline
fin.close()
