#!/usr/bin/env python
from rosalind import *

"""Solution to the Inferring mRNA from Protein problem in the Bioinformatics Stronghold section of Rosalind
location: http://rosalind.info/problems/mrna/
Problem: find number of possible mrna sequences from a given amino acid.
"""

if __name__ == "__main__":
	f = open('data/rosalind_mrna.txt').readlines()
	prot = f[0].strip()
	prot = sequence(prot)
	n = prot.prot2mrna()

	print(n)
