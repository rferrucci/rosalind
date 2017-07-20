#!/usr/bin/env python
from rosalind import *

"""Solution to the Locating Restriction Sites problem in the Bioinformatics Stronghold section of Rosalind
location: http://rosalind.info/problems/revp/
Problem: find restriction sites, i.e. reverse palindromes, in a given dna sequnce.
returns list of tuples: position, palindromic sequence, length
"""

if __name__ == "__main__":
	f = open('data/rosalind_revp.txt','r').readlines()

	fasta, seqs = fastaParse(f)
	seq = sequence(seqs[0].strip())

	seq.reverseComplement()
	sites  = seq.findRestrictionSites(12)
	for s in sites:
		print('{} {}'.format(s[0], s[2]))