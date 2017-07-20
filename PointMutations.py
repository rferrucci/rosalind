#!/usr/bin/env python
from rosalind import *

"""Solution to the Counting Point Mutations problem in the Bioinformatics 
Stronghold section of Rosalind
location: http://rosalind.info/problems/hamm/
Problem: return the number of differences between a pair of DNA sequences
"""

if __name__ == "__main__":
	f = open("data/rosalind_hamm.txt",'r').readlines()
	seq1 = sequence(f[0].strip())
	seq2 = sequence(f[1].strip())
	pointMut = seq1.pointMutations(seq2)
	print(pointMut)
