#!/usr/bin/env python

#get number of mRNA sequences possible from protein sequence, modelo 1,000,000

def prot2mrna(seq):
	f = open('./geneticcode.dat','r').readlines()
	revCode = {}
	N = 1
	for i in f:
		x = i.strip().split()
		if x[1] in revCode:
			revCode[x[1]] += 1
		else: revCode[x[1]] = 1
		
   
	#get number of possible sequences
	for i in range(len(seq)):
		N *= revCode[seq[i]]
	
	#don't forget stop codon
	return N % 1000000
	
if __name__ == "__main__":

    seq = open('data/rosalind_mrna.txt').read().strip()
    print (prot2mrna(seq))
