#
# Alphabetize an alignment
#
# Usage:
# Python alphabetize_alignment.py <filepath>
#
# Prints the results to STDOUT.
#

from argParser import *

from cogent import *
import math
import os
import sys

alignment = sys.argv[1]
print "Alignment = ", alignment

seq = LoadSeqs( alignment )

names = seq.getSeqNames()
names.sort()

for name in names:
    s = seq.getGappedSeq(name)
    print ">", name
    print s


