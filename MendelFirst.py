#!/usr/bin/env python
from rosalind import *

"""Solution to the Mendel's First Law problem in the Bioinformatics 
Stronghold section of Rosalind
location: http://rosalind.info/problems/iprb/
Problem: with a population of k homozygous dominant, m heterozygous, and and n are homozygous recessive individuals in a population, what is the probability that two randomly mating individuals will produce offspring that possesses a dominant allele
"""

if __name__ == "__main__":
	k, m, n = 2, 2, 2
	pops = k,m,n
	"""returns simulated and calculated estimates"""
	print(Mendel(pops,10000))
