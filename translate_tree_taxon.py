#
# Usage: translate_tree_taxa <newick tree path> <translate key>
#
# Where <translate key> is a text file with the following format:
# <taxon name in tree> <new taxon name>
#
# The reformatted tree will be printed to the screen.

import re, sys

ftree = open(sys.argv[2])
tree = ftree.readline()
ftree.close()

ftrans = open(sys.argv[1])
translines = ftrans.readlines()
ftrans.close()

for l in translines:
    tokens = l.split()
    old = tokens[0]
    new = tokens[1]
    tree = re.sub(old, new, tree)

print tree