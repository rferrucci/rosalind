#!/usr/bin/env python
from rosalind import *

"""Solution to the Transitions and Transversions problem in the Bioinformatics Stronghold section of Rosalind
location: http://rosalind.info/problems/tran/
Problem: find transition/transversion ratio of a given dna sequence
"""

if __name__ == "__main__":
	f = open('data/rosalind_tran.txt','r').readlines()
	fasta, sequences = fastaParse(f)
	seq = sequence(sequences[0])
	trans = seq.transitionsTransversions(sequences[1])
	print(trans)