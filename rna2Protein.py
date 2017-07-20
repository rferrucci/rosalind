#!/usr/bin/env python
from rosalind import *

"""Solution to the Translating RNA into Protein problem in the Bioinformatics 
Stronghold section of Rosalind
location: http://rosalind.info/problems/prot/
Problem: convert an RNA sequence in corresponding amino acid sequence
"""

if __name__ == "__main__":
	f = open("data/rosalind_prot.txt", 'r').readlines()
	rna = sequence(f[0])
	rna.getProteinSequence()
	print(rna.protein)
