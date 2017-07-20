#!/usr/bin/env python
from rosalind import *

"""Solution to the Calculating Protein Mass problem in the Bioinformatics 
Stronghold section of Rosalind
location: http://rosalind.info/problems/prtm/
Problem: return the protein mass of a given aa sequence
"""

if __name__ == "__main__":
	f = open("data/rosalind_prtm.txt", 'r').readlines()
	prot = sequence(f[0].strip())
	prot.proteinMass()