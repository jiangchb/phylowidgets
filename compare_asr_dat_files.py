#################################################
#
# USAGE:
#
# python compare_asr_dat_files.py [<id> <dat filepath> . . .]
#

import os
import sys
import re
from map_anc_2_anc import *
from argParser import ArgParser
ap = ArgParser(sys.argv)

anc1 = ap.getArg("--anc1")
nick1 = ap.getArg("--nick1")
msa1 = ap.getArg("--msa1")
seed1 = ap.getArg("--seed1")

anc2 = ap.getArg("--anc2")
nick2 = ap.getArg("--nick2")
msa2 = ap.getArg("--msa2")
seed2 = ap.getArg("--seed2")

runid = ap.getOptionalArg("--runid")
if runid == False:
    exit()

rsitesa = []
x = ap.getOptionalList("--restrict_sites_1")
if x != None:
    for s in x:
        rsitesa.append( int(s) )
rsitesb = []
x = ap.getOptionalList("--restrict_sites_2")
if x != None:
    for s in x:
        rsitesb.append( int(s) )
        
colors = ["blue", "red", "orange", "green", "black"]
names = ["indel change", "case 1", "case 2", "case 3", "no change"]

########################################################
#
# Do the paths exist?

if False == os.path.exists( anc1 ):
    print "\n\nERROR: the following path does not seem to exist:", anc1
    exit()
if False == os.path.exists( anc2 ):
    print "\n\nERROR: the following path does not seem to exist:", anc2
    exit()
if False == os.path.exists( msa1 ):
    print "\n\nERROR: the following path does not seem to exist:", msa1
    exit()
if False == os.path.exists( msa2 ):
    print "\n\nERROR: the following path does not seem to exist:", msa2
    exit()

#############################################################

#
# Function defs. go here. . .
#
def getprobs(ancpath, rsites):
    """Returns site_states_probs"""
    fin = open(ancpath, "r")
    lines = fin.readlines()
    fin.close()

    site_states_probs = {}
    for l in lines:
        tokens = l.split()
        site = int(tokens[0])
        if site not in rsites:
            continue
        
        # does this site contain a gap?
        if l.__contains__("-"):
            #print "75: I found a gap at site", site
            site_states_probs[ site ] = {}
            site_states_probs[site]["-"] = 0.0
        elif site not in rsites:
            continue
        else:
            site_states_probs[ site ] = {}
            i = 1
            while i < tokens.__len__():
                s = tokens[i]
                p = float(tokens[i+1])
                site_states_probs[site][s] = p
                i += 2
        
    return site_states_probs

#
def get_ml_state(states_probs):
    maxp = None
    maxc = None
    for state in states_probs:
        if maxp == None:
            maxp = states_probs[state]
            maxc = state
        elif maxp < states_probs[state]:
            maxp = states_probs[state]
            maxc = state
    return maxc
    
def renumber_sites(site_states_probs):
    sites = site_states_probs.keys()
    sites.sort()
    
    ret = {}
    for ii in range(0,sites.__len__() ):
        ret[ii+1] = site_states_probs[ sites[ii] ]    

    return ret

##########################################################

#
# "main" begins here. . .
#

# Generate a site map between the two ancestors
t = build_tuples(seed1, msa1, anc1, seed2, msa2, anc2)

#print t[0]
#print t[1]
#exit()

ra = []
rb = []
for ii in t:
    if ii[0] in rsitesa and ii[1] in rsitesb:
        ra.append(ii[0])
        rb.append(ii[1])

anc_data = {}
anc_data[nick1] = getprobs( anc1, ra )
anc_data[nick2] = getprobs( anc2, rb )

#
# Remove sites where both ancestors have "-"
#
for ii in t:
    if ii[0] in rsitesa and ii[1] in rsitesb:
        if anc_data[nick1][ ii[0] ].__contains__("-") and anc_data[nick2][ ii[1] ].__contains__("-"):
            print "Removing tuple", ii
            anc_data[nick1].__delitem__(ii[0])
            anc_data[nick2].__delitem__(ii[1])
        else:
            pass

nsites = anc_data[ nick1 ].__len__()

#
# Write a brief summary of the two ancestors
#
fout = open(runid + ".summary.txt", "w")
fout.write( "Comparing " + msa1 + " (" + nick1 + ") ")
fout.write("to " + msa2 + " (" + nick2 + ")\n")
fout.write("-------------------------------------------------------\n")

sum1 = 0
sum2 = 0

outl = ""
for ii in t:
    if ii[0] in anc_data[ nick1 ] and ii[1] in anc_data[ nick2 ]:
        state1 = get_ml_state( anc_data[ nick1 ][ii[0]] )
        state2 = get_ml_state( anc_data[ nick2 ][ii[1]] )
        sum1 += anc_data[ nick1 ][ii[0]][state1]
        sum2 += anc_data[ nick2 ][ii[1]][state2]
        outl += nick1 + ":" + ii[0].__str__() + "\t" + state1 + " (" + anc_data[ nick1 ][ii[0]][state1].__str__() + ")"
        outl += "\t\t" + nick2 + "\t" + ii[1].__str__() + "\t" + state2 + " (" + anc_data[ nick2 ][ii[1]][state2].__str__() + ")"
        if state1 != state2:
            outl += "\t\t"
            if state1 == "-" or state2 == "-":
                outl += "(indel)"
            elif anc_data[ nick1 ][ii[0]][state2] < 0.05 and anc_data[ nick2 ][ii[1]][state1] < 0.05:
                outl += "(case 1)"
            else:
                outl += "(case 2)"
            #elif anc_data[ nick1 ][ii[0]][state1] >= 0.8 and anc_data[ nick2 ][ii[1]][state2] < 0.8:
            #    outl += "(case 2)"
            #elif anc_data[ nick1 ][ii[0]][state1] < 0.8 and anc_data[ nick2 ][ii[1]][state2] >= 0.8:
            #    outl += "(case 2)"            
        else:
            outl += "\t"
            if anc_data[ nick1 ][ii[0]][state1] < 0.8 and anc_data[ nick2 ][ii[1]][state2] >= 0.8:
                outl += "(case 3)" 
            elif anc_data[ nick1 ][ii[0]][state1] >= 0.8 and anc_data[ nick2 ][ii[1]][state2] < 0.8:
                outl += "(case 3)" 
        outl += "\n"

fout.write("mean PP " + nick1 + " = %.3f"%(sum1/nsites) + "\n")
fout.write("mean PP " + nick2 + " = %.3f"%(sum2/nsites) + "\n")
fout.write("-------------------------------------------------------\n")
fout.write(outl)
fout.close()

#
# Renumber to 1-based site counting
#
anc_data[nick1] = renumber_sites( anc_data[nick1] )
anc_data[nick2] = renumber_sites( anc_data[nick2] )


#
# Sanity Check!
#
if anc_data[nick1].__len__() != anc_data[nick2].__len__():
    print "Ancestors have different lengths!"
    for site in anc_data[ nick1 ].keys():
        state1 = get_ml_state( anc_data[ nick1 ][site] )
        state2 = get_ml_state( anc_data[ nick2 ][site] )
        print site, state1, state2
    exit()


#
# Now compare the distributions. . . 
#

indelsites = []
redsites = []
orangesites = []
greensites = []
blacksites = []

print "\n"
for site in anc_data[ nick1 ].keys():
    indel = False
    red = False
    orange = False
    green = False
    state1 = get_ml_state( anc_data[ nick1 ][site] )
    state2 = get_ml_state( anc_data[ nick2 ][site] )
  
    # test for indel mismatch:
    if state2 != state1:
        if state2 == "-" or state1 == "-":
            indel = True
    
    # test for red:
    if state2 != state1 and indel == False:
        if False == anc_data[ nick2 ][site].keys().__contains__(state1):
            red = True
        elif False == anc_data[ nick1 ][site].keys().__contains__(state2):
            red = True
        elif anc_data[nick1][site][state2] < 0.05 and anc_data[nick2][site][state1] < 0.05:
            red = True
    
    # test for orange:
    if state2 != state1 and red == False and indel == False:
        orange = True
    
    # test for green
    if state2 == state1 and indel == False:
        if (anc_data[nick1][site][state1] < 0.8 and anc_data[nick2][site][state2] > 0.8) or (anc_data[nick1][site][state1] > 0.8 and anc_data[nick2][site][state2] < 0.8):
            green = True
    
    if indel:
        indelsites.append(site)
    if red:
        redsites.append(site)
    if orange:
        orangesites.append(site)
    if green:
        greensites.append(site)
    if indel == False and red == False and orange == False and green == False:
        blacksites.append(site)



print "============================================================"
print nick1 + " to " + nick2
print "nsites=", nsites
print "indel_mismatch=", indelsites.__len__(), indelsites
print "case 1:", redsites.__len__(), "sites:",  redsites
print "case 2:", orangesites.__len__(), "sites:",  orangesites
print "case 3:", greensites.__len__(), "sites:", greensites
print "case 4:", blacksites.__len__(), "sites:", blacksites
print ""


#
# Plot scatterplot: orange, red, and green sites
#

ix = "ix<-c("
iy = "iy<-c("
rx = "rx<-c("
ry = "ry<-c("
ox = "ox<-c("
oy = "oy<-c("
gx = "gx<-c("
gy = "gy<-c("
bx = "bx<-c("
by = "by<-c("

for i in indelsites:
    mlstate = get_ml_state( anc_data[nick1][i] )
    pa = anc_data[nick1][i][mlstate]
    mlstate = get_ml_state( anc_data[nick2][i] )
    pb = anc_data[nick2][i][mlstate]    
    ix += pa.__str__() + ","
    iy += pb.__str__() + ","

for s in redsites:
    mlstate = get_ml_state( anc_data[nick1][s] )
    pa = anc_data[nick1][s][mlstate]
    mlstate = get_ml_state( anc_data[nick2][s] )
    pb = anc_data[nick2][s][mlstate]    
    rx += pa.__str__() + ","
    ry += pb.__str__() + ","
for s in orangesites:
    mlstate = get_ml_state( anc_data[nick1][s] )
    pa = anc_data[nick1][s][mlstate]
    mlstate = get_ml_state( anc_data[nick2][s] )
    pb = anc_data[nick2][s][mlstate]    
    ox += pa.__str__() + ","
    oy += pb.__str__() + ","
for s in greensites:
    mlstate = get_ml_state( anc_data[nick1][s] )
    pa = anc_data[nick1][s][mlstate]
    mlstate = get_ml_state( anc_data[nick2][s] )
    pb = anc_data[nick2][s][mlstate]    
    gx += pa.__str__() + ","
    gy += pb.__str__() + ","
for s in blacksites: # black and green get plotted in the same bin
    mlstate = get_ml_state( anc_data[nick1][s] )
    pa = anc_data[nick1][s][mlstate]
    mlstate = get_ml_state( anc_data[nick2][s] )
    pb = anc_data[nick2][s][mlstate]    
    bx += pa.__str__() + ","
    by += pb.__str__() + ","

    
rx = re.sub(",$", "", rx)
rx += ");\n"
ry = re.sub(",$", "", ry)
ry += ");\n"
ox = re.sub(",$", "", ox)
ox += ");\n"
oy = re.sub(",$", "", oy)
oy += ");\n"
gx = re.sub(",$", "", gx)
gx += ");\n"

gy = re.sub(",$", "", gy)
gy += ");\n"

bx = re.sub(",$", "", bx)
bx += ");\n"

by = re.sub(",$", "", by)
by += ");\n"

ix = re.sub(",$", "", ix)
ix += ");\n"

iy = re.sub(",$", "", iy)
iy += ");\n"

cranstr = "colors <-c("
for c in colors:
    cranstr += "\""  + c.__str__() + "\","
cranstr = re.sub(",$", "", cranstr)
cranstr += ");\n"


#names = ["indel change", "case 1", "case 2", "case 3", "no change"]
cranstr += "n <-c("
cranstr += "\"indel change (" + indelsites.__len__().__str__() + ")\","
cranstr += "\"case 1 (" + redsites.__len__().__str__() + ")\","
cranstr += "\"case 2 (" + orangesites.__len__().__str__() + ")\","
cranstr += "\"case 3 (" + greensites.__len__().__str__() + ")\","
cranstr += "\"no change (" + blacksites.__len__().__str__() + ")\");\n"

fout = open(runid + ".scatter.rscript", "w")
cranstr += "pdf(\"" + runid + ".scatter.pdf\", width=4, height=4);\n" 
fout.write(cranstr)

fout.write("plot(c(0.0,1.0),c(0.0,1.0),type='n',xlim=range(0,1.0),ylim=range(0,1.0),xlab=\"" + nick1 + " PPs\",ylab=\"" + nick2 + " PPs\",main=\"" + runid + "\");\n")
fout.write(rx)
fout.write(ry)
fout.write(ox)
fout.write(oy)
fout.write(gx)
fout.write(gy)
fout.write(bx)
fout.write(by)
fout.write(ix)
fout.write(iy)
fout.write("rect(0.8,0.8,1.0,1.0,lty=\"blank\",col=\"lightyellow\");\n")
fout.write("points(bx,by,col='black',pch=1);\n")
fout.write("points(ix,iy,col='blue',pch=2);\n")
fout.write("points(gx,gy,col='green',pch=5);\n")
fout.write("points(ox,oy,col='orange',pch=20);\n")
fout.write("points(rx,ry,col='red',pch=18);\n")
fout.write("abline(0,1.0,col='black',lwd=0.5,lty=\"dashed\");\n")

fout.write("legend(0.1,1.0,n,col=colors,pch=c(2,18,20,5,1),pt.cex=0.5,cex=0.5,box.lwd=0.5,box.lty=\"dashed\");\n")
fout.write("dev.off();\n")

fout.close()
os.system("r --no-save < " + runid + ".scatter.rscript")


#
# Plot histogram: proportion of orange, red, and green sites, relative to total
#

#
# data[xgroup][series] = value
#
def pieplot(bars, names, colors, xlab, ylab, filekeyword):
                
    pdfpath = filekeyword + ".pie.pdf"
    #cranstr = "pdf(\"" + pdfpath + "\", width=6, height=3);\n"    
    cranstr = "pdf(\"" + pdfpath + "\", width=6, height=6);\n"
    
    cranstr += "y <-c("
    for b in bars:
        cranstr += b.__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"
    
    cranstr += "colors <-c("
    for c in colors:
        cranstr += "\""  + c.__str__() + "\","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"

    cranstr += "n <-c("
    for n in names:
        cranstr += "\"" + n.__str__() + "\","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"  

    cranstr += "pie(y, labels=n, col=colors, main=\"" + filekeyword + "\");\n"    
    cranstr += "dev.off();\n"  

    fout = open(filekeyword + ".pie.rscript", "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save < " + filekeyword + ".pie.rscript")

bd = []
nsites = 1.0
bd.append( float(indelsites.__len__()) / nsites )
bd.append( float(redsites.__len__()) / nsites )
bd.append( float(orangesites.__len__()) / nsites )
bd.append( float(greensites.__len__()) / nsites )
bd.append( float(blacksites.__len__()) / nsites )
pieplot(bd,names,colors,"Sets", "Proportion of Sites", runid)

