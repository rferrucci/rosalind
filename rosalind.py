#!/usr/bin/env python
import random
from itertools import permutations, product
from math import *
#from string import maketrans

# I wrote this class while working through exercises from the Rosalind project 
# (http://rosalind.info/) of diyBIO. It is a set of functions for analyzing DNA
# sequences. Some background information on the name of the project: Rosalind 
# Franklin was a biologist whose work contributed to the elucidation of the 
# double helical structure of DNA. Dr Franklin did not share in the Noble Prize
# in Chemistry with Watson, Crick, and Maurice Wilkins---Watson and Crick 
# actually deduced the structure, while  Wilkins developed a method for 
# obtaining crystalized structures of nucleic acids that led to Franklin being 
# able to determine that DNA was helical in the first place for her 
# contribution as Noble Prizes are not awarded posthumously

import re
import itertools
import string
#from containers import Counter
	
def fastaParse(dat):
	#takes a set of fasta and sequences and returns two lists, one of the fasta headers and the other of the sequences
	seqs = []
	fasta = []
	count = -1
	for i in dat:
		if (i[0] == '>' and count == -1): 
			fasta.append(i.strip())
			seqs.append("")
			count += 1
			#dna = ""
		elif (i[0] == '>' and count != -1):
			fasta.append(i.strip())
			seqs.append("")
			count += 1
			#dna = ""
		elif (i == '\n'):
			continue
		else:
			seqs[count]+=i.strip()
	return fasta, seqs

class sequence(str):
	#since a DNA sequence is essentially a string, I felt it appropriate to inherit from the string class
	def __init__(self, dnaseq):
		#initiate class with the dnaseq, making sure that our sequence is all uppercase
		self.seq = dnaseq.upper()
		self.n = len(dnaseq)

	def nucleotideCount(self, dna='y'):
		#return dictionary of nucleotide counts from DNA sequence
		seq = self.seq
		if dna=='y': nuc = ['A','T','C','G']
		else: nuc = ['A','U','C','G']
		nucleotides = dict(zip(nuc, [seq.count(n) for n in nuc]))			
		return nucleotides
	
	def convertToRNA(self):
		#return transcribed RNA string based on coding strand
		seq = self.seq
		rna = re.sub('T','U',seq)
		#seq.replace("T", "U") works as well
		self.rna = rna
		return rna
	
	def reverseComplement(self):
		#returns the reverse complement of a strand of DNA
		#compDict= {'A':'T','T':'A','C':'G','G':'C'}
		seq = self.seq
		complement = seq[::-1].translate(str.maketrans('ATCG','TAGC'))
		self.complement = complement
		return complement
	
	def getGCContent(self):
		#returns the GC content of a dna sequence
		seq = self.seq
		n = len(''.join(seq.splitlines()))
		GC = (seq.count('G') + seq.count('C'))/float(n) * 100
		return GC
	
	def pointMutations(self, seq2, freq=0):
		#returns number of point mutations between a pair of sequences. changing option freq to be equal to one will output frequency of point mutations
		seq1 = self.seq
		n = self.n
		m = 0

		mut = sum(a != b for a, b in zip(seq1, seq2))  

		if freq==0: return mut
		if freq==1: return mut/float(n)

	def getProteinSequence(self):
		#translates RNA sequence in corresponding aa sequence. Genetic code is obtained from geneticcode.dat
		f = open('geneticcode.dat','r').readlines()
		geneticcode = {c: a for c, a in [l.split() for l in f]}
			
		seq=self.seq
		protein = ""
		protein = ''.join([geneticcode[seq[i:i+3]] for i in range(0, len(seq), 3)])
		self.protein = protein.rstrip('Stop')
	
	def findMotif(self,motif):
		#returns a given motif in a DNA sequence
		seq = self.seq
		locations = []
		n = len(motif)
		for i in range(0, len(seq)-n):
			if motif == seq[i:i+n]:
				locations.append(str(i+1))
		
		return locations
	
	def spliceIntrons(self, introns):
		# removes introns, from list of introns, from dna sequence outputting a sequence containing just exons
		seq =  self.seq
		for i in introns:
			seq = seq.replace(i,"")
			
		self.mDNA = seq

	def getSharedMotifs(self, seqs):
		#takes a list of sequences and returns a list of motifs found amongst the set of sequences
		seq1 = self.seq
		seq2 = seqs[0]
		#seq3 = seqs[1]
		motifs = []
		for i in range(1, len(seq1)):
			seqLength = len(seq1) - i + 1
			#first we prepare a list of motifs shared between the 
			#first two sequences
			for j in range(0, len(seq2)-seqLength + 1):
				motif=seq1[j:j+seqLength]
				if seq2.find(motif)!=-1: 
					motifs.append(motif)
			
			#compare the list of motifs to each additional sequence, removing
			# motifs not found
			for seq in seqs[1:]:
				for m in motifs:
					if seq.find(m)==-1:
						motifs.remove(m)

		return motifs

	def getConsensusSequence(self, sequences):
		#returns consensus sequences amongst set of sequences
		seq = self.seq
		nuc = 'ATCG'

		#group nucleoties by position
		columns = list(zip(*sequences))

		#create profile of counts for each nucleotide for each position
		n = len(columns)
		profile = [{n: columns[i].count(n) for n in nuc} for i in range(n)]
		
		#generate consensus sequence based on the maximum frequency for each nucleatde
		consensus = ''.join([max(profile[i].keys(), key=(lambda k: profile[i][k])) for i in range(len(profile))])
		
		self.consensus = consensus
		self.profile = profile
		return consensus,profile

	def transitionsTransversions(self, seq2):
		#given two sequences, calculates the ration of transitions to transversions
		seq1 = self.seq
		trans={'A':'G','G':'A','C':'T','T':'C'}
		mut = [0 if trans[a]==b else 1 for a, b in zip(seq1, seq2) if a!= b]
		r = mut.count(0)/float(mut.count(1))
		return r
		
	def getORF(self):
		#takes rna sequence and outputs list of protein sequences based on ORF from start to stop codons
		rna = self.seq
		
		f = open('geneticcode.dat','r').readlines()
		genCode = {i.split()[0]: i.split()[1] for i in f}
		#get locations of start codons
		
		Prot= []
		for i in range(len(rna)-2):
			if rna[i:i+3] == "AUG":
				#start protein sequence once start codon is found
				prot = ""
				for j in range(i,len(rna)-2,3):
					if genCode[rna[j:j+3]] != 'Stop':
						#if stop codon is not encountered, add aa to sequence
						prot+=genCode[rna[j:j+3]]
					else:
						#once stop codon is encountered, add protein to list of sequences and start over.
						Prot.append(prot)
						break	
				
		return Prot
	
	def proteinMass(self):
		prot = self.seq
		f = open("proteinmass.txt", 'r').readlines()
		dat = [l.split() for l in f]
		masstable = {a:float(m) for a, m in dat}
		mass = sum(masstable[aa] for aa in prot)
		print(mass)
		
	def findRestrictionSites(self, R):
		#finds restriction sites in length of DNA. First obtain the reverse complement using the appopriate method from this pacakge
		#returns tuple with location, sequence, and length. R is maximum length of restriction sites
		seq = self.seq
		trans = str.maketrans('ACTG','TGAC')
		#complement = self.complement
		restrictionSites=[]
		print(seq)
		#print(complement)
		for i in range(0, len(seq)):	
			for n in range(4,R+1,2):
				s = seq[i:i+n]
				if len(s) < n: continue
				r = s[::-1].translate(trans)
				if len(seq) -i < n: continue
				if s==r: restrictionSites.append((i+1, r, n))

		return restrictionSites

	def prot2mrna(self):
		f = open('./geneticcode.dat','r').readlines()
		seq = self.seq

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
		N *= 3 

		return N % 1000000

	def splicedMotif(self, mot):
		#works to find first spliced motif, looking for finding all
		#given a DNA motif, find the indices of nuclotides when the motif is not contiguous
		s = 0
		indices = []
		seq = self.seq			 
					
		i = 0
					
		for m in mot:
			i = (i + seq[i:].index(m)) + 1
			indices.append(i)
			
		#for i in range(len(seq)):
		#	if s == len(mot): break
		#	if mot[s] == seq[i]:
		#		indices.append(i+1)
		#		s+=1
		return indices
		
#------------------------------------------------------------------------------#
#accesory functions
#------------------------------------------------------------------------------#

# not a part of the DNA class. May have to make a parse class or something
def fastaParse(dat):
	#takes a set of fasta and sequences and returns two lists, one of the fasta headers and the other of the sequences
	seqs = []
	fasta = []
	count = -1
	for i in dat:
		if (i[0] == '>' and count == -1): 
			fasta.append(i.strip())
			seqs.append("")
			count += 1
			#dna = ""
		elif (i[0] == '>' and count != -1):
			fasta.append(i.strip())
			seqs.append("")
			count += 1
			#dna = ""
		elif (i == '\n'):
			continue
		else:
			seqs[count]+=i.strip()
	return fasta, seqs
	
#------------------------------------------------------------------------------#
#not included in the class. demonstrative
#------------------------------------------------------------------------------#

def Mendel(pops, P):
	#demonstrates Mendel's First Law. P is number of simulations to run.
	p = 0
	A = 0
	a = 1
	k,m,n = pops[0],pops[1],pops[2]
	pop = []
	for i in range(k): pop.append((A,A))
	for i in range(m): pop.append((A,a))
	for i in range(n): pop.append((a,a))

	for i in range(P):
		samp = random.sample(pop,2)

		a1 = random.sample(samp[0],1)[0]
		a2 = random.sample(samp[1],1)[0]
		if a1==0 or a2 == 0: p+=1
		
	simulated = p/float(P)
	#to figure out mathmatically
	N = float(k + m + n)

	r_r = (n /N) * ((n - 1) / float(N - 1))
	h_h = (m /N) * ((m - 1) / float(N - 1))
	h_r = (m /N) * (n / (N - 1)) +  (n / N) * (m / (N - 1))

	recessive = r_r + h_h * 1/4.0 + h_r * 1/2.0
	calculated= 1 - recessive
	return simulated, calculated

def calculatingExpOffspring(couplePairs, prob):
	#takes list of number of couple pairs with a given genotype: AA-AA, AA-Aa, 
	#AA-aa, Aa-Aa, Aa-aa, aa-aa and returns the expected number of those with the 
	#dom phenotype
	p = sum([2 * a * b for a, b in zip(couplePairs, prob)])
	return p

def rabbitPairs(numMonths, numOffspring):
	#code for Rabbits and Recurrence Relations problem.
	F = []
	for n in range(numMonths):
		if n == 0:
		    F.append(1)
		elif n == 1:
			F.append(1)
		else:
			f = F[n-1] + F[n-2]*numOffspring
			F.append(f)
	return F[-1]

def mortalFibRabbits(numMonths, lifespan):
	F = []
	Ages = [0 for i in range(lifespan)]
	Ages[0] = 1

	for n in range(numMonths - 1):
		f = sum(Ages[1:])
		Ages.pop()
		Ages.insert(0, f)
	return Ages

def longestIncreasingSubsequence(n, seq):
	"""find longest increasing or decreasing subsequence given an array of integers
	"""
	P = [None for i in range(n)]
	M = [None for i in range(n + 1)]
	L = 0

	for i in range(n):
		lo = 1
		hi = L

		while (lo <= hi):
			mid = (lo+hi)/2
			if seq[M[mid]] < seq[i]:
				lo = mid + 1
			else:
				hi = mid - 1

		newL = lo
		P[i] = M[newL - 1]
		M[newL] = i

		if newL > L:
			L = newL

	S = []
	k = M[L]
	for i in range(L-1, -1, -1):
		S.append(seq[k])

		k = P[k]

	return(S[::-1])

def kmerLex(letters, sequence, k=4):
	kmers = [''.join(p) for p in product(letter, repeat=4)]
	kMers = {k:0 for k in kmers}

	for i in range(len(sequence) - k):
		kMers[sequence[i:i+k]] +=1

	return kmers, kMers

def overlapGraphs(fasta, seqs, k):
	fasta = [f.lstrip(">") for f in fasta]
	adjList = []
	for i in range(0, len(seqs)-1):
		for j in range(i + 1, len(seqs)):
			if seqs[i][-3:] == seqs[j][:3]:
				adjList.append((fasta[i], fasta[j]))
	return adjList

def randomStrings(seq, GC):
	Prob = []
	for gc in GC:
		GC = {'C': gc/2, 'G':gc/2, 'T': (1-gc)/2, 'A': (1-gc)/2}

		prob = sum(map(log10,[GC[s] for s in seq]))
		Prob.append(prob)
	return Prob

def geneOrders(n):
	N = [i + 1 for i in range(0, n)]
	perm = []
	
	for p in permutations(N):
		perm.append(' '.join(map(str, p)))
	
	return perm
