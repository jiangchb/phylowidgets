import math
import os
import re
import sys

from argParser import ArgParser
argp = ArgParser(sys.argv)

spath = argp.getArg("--input")
fin = open( spath, "r" )
lines = fin.readlines()
fin.close()

output = argp.getArg("--output")

def plot_comb(cat):
    cranscript = "plot." + output + "." + cat.__str__() + ".cran"
    f = open(cranscript, "w")
 
    y = "y <- c("
    x = "x <- c("
    
    for site in site_pps:
        y += site_pps[site][cat].__str__() + ","
        x += site.__str__() + ","

    y = re.sub(",$", "", y)
    y += ")"
    f.write( y + "\n")  
    
    x = re.sub(",$", "", x)
    x += ")"
    f.write( x + "\n")
    
    f.write("pdf('plot." + output.__str__() + "." + cat.__str__() + ".pdf', width=8, height=3);\n")
    f.write("barplot(y, x, xlab='site', ylab='PP(site | cat)', ylim=range(0,1.0), lwd=2, col=\"black\", main ='" + output.__str__() + "." + cat.__str__() + "');\n")
    f.write("dev.off();\n")
    f.close()
    os.system("r --no-save < " + cranscript)


site_pps = {}
found_first_line = False
for line in lines:
    if line.startswith("0"):
        found_first_line = True
        
    if found_first_line == False:
        continue
    else:
        tokens = line.strip().split()
        site = int(tokens[0])
        print site
        lks = tokens[2:]
        sum = 0.0
        for i in range(0, lks.__len__() ):
            lks[i] = float(lks[i])
            sum += math.exp( lks[i] )

        pps = []
        for i in range(0, lks.__len__() ):
            pp = math.exp(lks[i]) / sum
            if pp < 0.01:
                pp = 0.0
            pps.append( math.exp(lks[i]) / sum )
        
        site_pps[site] = pps
        #print site, pps

for cat in range(0, site_pps[0].__len__()):
    plot_comb(cat)
    
