#!/usr/bin/env python
import sys
firstarg=sys.argv[1]
import Bio
from Bio import AlignIO
alignment = AlignIO.read(open(firstarg),"fasta")
print("Alignment length %i" % alignment.get_alignment_length())
