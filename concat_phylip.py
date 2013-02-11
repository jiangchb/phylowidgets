#
# Concats two PHYLIP alignments
#
# python concat_phylip <phylip_path1> <phylip_path2>
#
# The concat-ed alignement is printed to STDOUT.
#

import sys

def read_alignment( path ):
    fin = open(path, "r")
    lines = fin.readlines()
    fin.close()
    a = {}
    for i in range(1,lines.__len__()): #skip the first line of the PHYLIP alignment
        line = lines[i].strip()
        tokens = line.split()
        taxa = tokens[0]
        sequence = tokens[1]
        a[taxa] = sequence
    return a


def write_alignment(a1, a2):
    ntax = a1.__len__()
    length = a1[ a1.keys()[0] ].__len__() + a2[ a2.keys()[0] ].__len__()
    print " " + ntax.__str__() + " " + length.__str__() 
    for taxa in a1.keys():
        if False == a2.__contains__(taxa):
            pass
            #print "WARNING: ignoring", taxa, "because it does not exist in file", sys.argv[2]
        else:
            print taxa + " " + a1[taxa] + a2[taxa] 


a1 = read_alignment( sys.argv[1] )
a2 = read_alignment( sys.argv[2] )
write_alignment(a1, a2)