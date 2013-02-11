#
# This script converts *.trprobs tree file (from Mr. Bayes) into a text file
# containing Newick trees.
#

import sys
import re

input = sys.argv[1] # our nexus tree file
output = sys.argv[2] # our desired newick tree file

f = open(input, "r")
g = open(output, "w")

translate = {} # key = number, value = taxon name
trees = [] # each cell in this array contains an untranslated newick tree

stateTranslate = False # are we parsing through the translation section?

lines = f.readlines()
for l in lines:
    #print "line = " + l
    if re.match("^\s+\n", l):
        continue
    
    if l.__contains__("translate") or l.__contains__("TRANSLATE") or l.__contains__("Translate"):
        stateTranslate = True
        continue
    elif stateTranslate == True:
        if re.match("^\s+\d+\s+\w+", l):
            l = l.strip()
            l = re.sub(",", "", l)
            l = re.sub(";", "", l)
            tokens = l.split()
            number = int(tokens[0])
            name = tokens[1]
            translate[number] = name
            #print "found a translation: " + number.__str__() + ", " + name
    
    if re.match("^\s+tree", l):
        #print "line = " + l
        # this line contains a tree:
        stateTranslate = False
        l = l.strip()
        tokens = l.split()
        tree = tokens[ tokens.__len__() -1 ]
        trees.append(tree)
        #print "found a tree: " + tree

print "I found " + trees.__len__().__str__() + " trees in " + sys.argv[1]


sys.stdout.write("Translating your trees")
#
# This three-level loop could use some optimization
#
nums = translate.keys()
nums.reverse()
#print "numbers = " + nums.__str__()
translatedTrees = []
for t in trees:
    sys.stdout.write(".")
    sys.stdout.flush()
    #print "erg tree = " + t
    tokens = t.split(",")
    translatedTokens = []
    for k in tokens:
        for s in nums:
            if k.__contains__(s.__str__()):
                k = re.sub(s.__str__(), translate[s], k)
                translatedTokens.append(k)
                break
    
    #print translatedTokens.__str__()
    
    transTree = ""
    for x in translatedTokens:
        transTree += x.__str__() + ","
    #print transTree
    transTree = re.sub(";,", ";", transTree)
    #print "trans tree = " + transTree
    translatedTrees.append(transTree)
    
    g.write(transTree + "\n")

print ""
    
f.close()
g.close()