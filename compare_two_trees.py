import os
import sys
from dendropy import Tree

t1path = sys.argv[1]
t2path = sys.argv[2]

t1 = Tree()
t1.read_from_path(t1path, "newick")
t2 = Tree()
t2.read_from_path(t2path, "newick")


s = t1.symmetric_difference(t2)
s = t2.symmetric_difference(t1)
print "symmetric diff. = ", s

print t1.length()
print t2.length()

