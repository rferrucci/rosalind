#!/usr/bin/env python
from rosalind import *

"""Solution to the Mortal Fibonacci Rabbits problem in the Bioinformatics 
Stronghold section of Rosalind
location: http://rosalind.info/problems/fibd/
Problem: return the total number of pairs of rabbits that will remain after the nth month if all rabbits live for m months.
"""
if __name__ == "__main__":
	f = open('data/rosalind_fibd.txt','r').readlines()
	m, n = map(int, f[0].split())
	print sum(mortalFibRabbits(m, n))
