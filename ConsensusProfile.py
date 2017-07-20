#!/usr/bin/env python
from rosalind import *

"""Solution to the Consensus and Profile problem in the Bioinformatics 
Stronghold section of Rosalind
location: http://rosalind.info/problems/cons/
Problem: given a set of sequences of DNA, return the consensus sequence and profile of nucleotide frequencies at each loci 
"""

if __name__ == "__main__":
	f = open("data/rosalind_cons.txt", 'r').readlines()
	fasta, sequences = fastaParse(f)
	seq = sequence(sequences[0])
	dat = sequences
	cons, profile = seq.getConsensusSequence(dat)
	nuc = ['A','C','G','T']
	print (cons)
	for n in nuc:
		l = "%s: %s" % (n, ' '.join([str(profile[i][n]) for i in range(len(profile))]))
		print(l)
