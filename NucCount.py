#!/usr/bin/env python
from rosalind import *

"""Solution to the Transcribing DNA into RNA problem in the Bioinformatics Stronghold section of Rosalind
location: http://rosalind.info/problems/dna/
Problem: transcribe dna sequence to rna
"""

# outputs counts of nucleotides from a sequence
if __name__ == "__main__":	
	f = open("data/rosalind_rna.txt",'r').readlines()
	fasta, seqs = fastaParse(f)
	seq = sequence(seqs[0])
	rna = seq.convertToRNA()
	print(rna)
