#!/usr/bin/env python

from rosalind import *

"""Solution to Longest Increasing Subsequence problem in the Bioinformatics Stronghold  
section of Rosalind
location: http://rosalind.info/problems/lgis/
Problem: Given k arrays of size n, output a list of paired indices i, j such that 
A[i] = -A[j] for each array A. Otherwise, output -1
"""

if __name__ == "__main__":
	f = open('data/rosalind_lgis.txt','r').readlines()

	n = int(f[0].strip())
	seq = map(int, f[1].split())
	inc = longestIncreasingSubsequence(len(seq), seq)
	seq.reverse()
	dec = longestIncreasingSubsequence(len(seq), seq)
	dec.reverse()
	print ' '.join(map(str, inc))
	print ' '.join(map(str, dec))
