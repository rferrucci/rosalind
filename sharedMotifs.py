#!/usr/bin/env python
from rosalind import *

"""Solution to the Finding a Shared Motif problem in the Bioinformatics Stronghold section of Rosalind
location: http://rosalind.info/problems/lcsm/
Problem: find longest common motif shared amongst a set of sequences.
"""

if __name__ == "__main__":
	f = open('data/rosalind_lcsm.txt').readlines()
	fasta, sequences = fastaParse(f)
	seq = sequence(sequences[0])
	motifs = seq.getSharedMotifs(sequences[1:])
	print(motifs)
	