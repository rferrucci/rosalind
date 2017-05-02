#!/usr/bin/env python
from rosalind import *

# return the reverse complement of a strand of DNA
if __name__ == "__main__":	
	f = open("data/rosalind_cmp.txt",'r').readlines()
	fasta, seqs = fastaParse(f)
	seq = sequence(seqs[0])
	cmp = seq.reverseComplement()
	print(cmp)
