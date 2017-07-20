#!/usr/bin/env python
from rosalind import *

"""Solution to the RNA Splicing problem in the Bioinformatics Stronghold section of Rosalind
location: http://rosalind.info/problems/splc/
Problem: find longest common motif shared amongst a set of sequences.
"""

if __name__ == "__main__":
	f=open('data/rosalind_splc.txt','r').readlines()

	fasta, sequences = fastaParse(f)

	dna = sequences[0].strip()

	introns=[sequences[i].strip() for i in range(1, len(sequences))]
	
	seq = sequence(dna)	
	seq.spliceIntrons(introns)

	dna = sequence(seq.mDNA)
	dna.convertToRNA()
	rna = sequence(dna.rna)

	rna.getProteinSequence()

