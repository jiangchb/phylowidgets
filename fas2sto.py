#
# Converts FASTA-formatted alignment into Stockholm format.
#

from Bio import AlignIO
import sys
al = AlignIO.read( open(sys.argv[1]), "fasta" )
print al.format("stockholm")