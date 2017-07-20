#!/usr/bin/env python

import rosalind

"""Solution to the Openm Reading Frames problem in the Bioinformatics Stronghold section of Rosalind
location: http://rosalind.info/problems/lcsm/
Problem: find all possible protein strings along a sequence of dna.
"""

if __name__ == "__main__":
	f = open('data/rosalind_orf.txt','r').readlines()
	dat = map(lambda i: i.strip(), f)
	fasta, sequences = rosalind.fastaParse(dat)
	seq = rosalind.sequence(sequences[0])

	seq.convertToRNA()
	rna = rosalind.sequence(seq.rna)
	seq.reverseComplement()
	comp = rosalind.sequence(seq.complement)

	ORF = rna.getORF()

	comp.convertToRNA()
	rna = rosalind.sequence(comp.rna)
	ORF += rna.getORF()

	for p in set(ORF):
		print p
