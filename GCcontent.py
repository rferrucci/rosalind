#!/usr/bin/env python
from rosalind import *

"""Solution to the Computing GC Content in the Bioinformatics 
Stronghold section of Rosalind
location: http://rosalind.info/problems/revc/
Problem: return the reverse complement of a dna sequence
"""

# return the reverse complement of a strand of DNA
if __name__ == "__main__":
	f = open('data/rosalind_gc.txt').readlines()
	fasta, sequences = fastaParse(f)
	fas=""
	GC=0
	for i in range(len(fasta)):
		seq = sequence(sequences[i].strip())
		gc= seq.getGCContent()
		if gc > GC: 
			fas = fasta[i]
			GC = gc
	
	print "%s\n%.6f" %(fas,GC)