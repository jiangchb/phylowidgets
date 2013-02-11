from Bio import AlignIO
import sys
al = AlignIO.read( open(sys.argv[1]), "stockholm" )
print al.format("fasta")