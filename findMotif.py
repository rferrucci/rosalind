#!/usr/bin/env python
from rosalind import *

"""Solution to the Finding a Motif in DNA problem in the Bioinformatics 
Stronghold section of Rosalind
location: http://rosalind.info/problems/subs/
Problem: given a sequence of DNA, find the locations of a given motif
"""

if __name__ == "__main__":
	f = open("data/rosalind_subs.txt", 'r').readlines()
	dna = sequence(f[0].strip())
	motif = f[1].strip()
	
	print(' '.join(dna.findMotif(motif)))
