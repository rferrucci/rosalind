#!/usr/bin/env python

import rosalind

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
