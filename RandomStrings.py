#!/usr/bin/env python
from rosalind import *
from math import *
"""Solution to the Introduction to Random Strings problem in the Bioinformatics
Stronghold section of Rosalind
location: http://rosalind.info/problems/prob/
Problem: given a DNA sequence and an array of GC content, output the log of the
probability constructing the DNA sequence at random
"""
 
if __name__ == "__main__":	
	f = open("data/rosalind_prob.txt",'r').readlines()
	seq = f[0].strip()
	dat = map(float, f[1].split())
	Prob = randomStrings(seq, dat)
	print(" ".join(map(str, Prob)))
			
