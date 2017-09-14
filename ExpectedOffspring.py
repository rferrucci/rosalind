#!/usr/bin/env python
from rosalind import *

"""Solution to the Calculating Expected Offspring problem in the Bioinformatics Stronghold section of Rosalind
location: http://rosalind.info/problems/iev/
Problem: given the number of couples possessing certain genotypes, what is the expected number of offspring 
displaying to dominant phenotype.
"""

if __name__ == "__main__":
	f='data/rosalind_iev.txt'
	f = open(f,'r')
	couplePairs = [float(x) for x in f.readline().split()]
	prob = [1.0,1.0,1.0,0.75,0.5,0]
	
	p = calculatingExpOffspring(couplePairs, prob)
	print(p)
