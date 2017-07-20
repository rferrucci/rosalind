#!/usr/bin/env python
from rosalind import *

"""
given an array of integers, if A[i] = -A[j], return indices i, j. else, return -1
"""

if __name__ == "__main__":	
	f = open('data/rosalind_sseq.txt','r').readlines()
	dat = map(lambda i: i.strip(), f)
	fasta, sequences = fastaParse(dat)
	sub = f[3].strip()
	seq = sequence(sequences[0])
	indices = seq.splicedMotif(sub)
