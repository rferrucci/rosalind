#!/usr/bin/env python
from rosalind import *

"""Solution to the Counting DNA Nucleotides problem in the Bioinformatics Stronghold section of Rosalind
location: http://rosalind.info/problems/dna/
Problem: return count of nucleotides for given DNA sequence
"""

# outputs counts of nucleotides from a sequence
if __name__ == "__main__":	
	f = open("data/rosalind_dna.txt",'r').readlines()
	fasta, seqs = fastaParse(f)
	seq = sequence(seqs[0])
	counts = seq.nucleotideCount()
	print(" ".join([str(c) for c in counts.values()]))
	
