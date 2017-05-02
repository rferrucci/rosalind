#!/usr/bin/env python
from rosalind import *

# outputs counts of nucleotides from a sequence
if __name__ == "__main__":	
	f = open("data/rosalind_dna.txt",'r').readlines()
	fasta, seqs = fastaParse(f)
	seq = sequence(seqs[0])
	counts = seq.nucleotideCount()
	print(" ".join([str(c) for c in counts.values()]))
	
