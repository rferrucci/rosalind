#!/usr/bin/env python
from rosalind import *

"""Solution to the Inferring mRNA from Protein problem in the Bioinformatics 
Stronghold section of Rosalind
location: http://rosalind.info/problems/mrna/
Problem: return the number of differences between a pair of DNA sequences
"""
if __name__ == "__main__":

	seq = open('data/rosalind_mrna.txt').read()
	seq = sequence(seq.strip())
	
	print (seq.prot2mrna())
