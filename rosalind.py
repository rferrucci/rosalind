#!/usr/bin/env python

# I wrote this class while working through exercises from the Rosalind project (http://rosalind.info/) of diyBIO. It is a set of functions for analyzing DNA sequences.
#
# Some background information on the name of the project: Rosalind Franklin was a biologist whose work contributed to the elucidation of the double helical structure 
# of DNA. Dr Franklin did not share in the Noble Prize in Chemistry with Watson, Crick, and Maurice Wilkins---Watson and Crick actually deduced the structure, while
# Wilkins developed a method for obtaining crystalized structures of nucleic acids that led to Franklin being able to determine that DNA was helical in the first place
# for her contribution as Noble Prizes are not awarded posthumously

import re
import itertools
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

	def nucleotideCount(self):
		#return dictionary of nucleotide counts from DNA sequence
		seq = self.seq
		nuc = ['A','T','C','G']
		nucleotides = dict(zip(nuc, [seq.count(n) for n in nuc]))			
		return nucleotides
	
	def convertToRNA(self):
		#return transcribed RNA string based on coding strand
		seq = self.seq
		rna = re.sub('T','U',seq)
		#seq.replace("T", "U") works as well
		self.rna = rna
	
	def reverseComplement(self):
		#returns the reverse complement of a strand of DNA
		compDict= {'A':'T','T':'A','C':'G','G':'C'}
		seq = self.seq
		complement = ''.join([compDict[n] for n in seq])[::-1]
		#complement = (seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))) #use this approach for python 3.0
		self.complement = complement
	
	def getGCContent(self):
		#returns the GC content of a dna sequence
		seq = self.seq
		n = len(''.join(seq.splitlines()))
		GC = (seq.count('G') + seq.count('C'))/float(n) * 100
		return GC
	
	def pointMutations(self, seq2, freq=0):
		#returns number of point mutations between a pair of sequences. changing option freq to be equal to one will output frequency of point mutations
		seq = self.seq
		n = self.n
		m = 0

		mut = sum([(1 if seq[i] != seq2[i] else 0) for i in range(len(seq))])
		#mut = sum(a != b for a, b in itertools.izip(s1, s2))   #apparently this is much faster and saves memory
		if freq==0: return mut
		if freq==1: return mut/float(n)

	def getProteinSequence(self):
		#translates RNA sequence in corresponding aa sequence. Genetic code is obtained from geneticcode.dat
		geneticcode = {}
		for i in open('geneticcode.dat','r'):
			x = i.split()
			geneticcode[x[0]] = x[1]
			
		#geneticcode =	{c:a for c,a in [s.strip().split(' ') for s in open('geneticcode.dat','r') ]} more elegant approach than mine

		seq=self.seq
		protein = ""
		for i in range(0,len(seq),3):
			codon = ''.join(seq[i:i+3])
			protein += geneticcode[codon]
			
		#protein = ''.join([geneticcode[''.join(seq[i:i+3])]for i in range(0,len(seq),3) ]) bit more elegant
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
			for j in range(0, len(seq2)-seqLength + 1):
				#print seqLength
				motif=seq1[j:j+seqLength]
				if seq2.find(motif)!=-1: 
					motifs.append(motif)
		#print seq
 		for seq in seqs[1:]:
			for m in motifs:
				if seq.find(m)==-1:
					motifs.remove(m)
		
		return motifs
	
	def findRestrictionSites(self, R):
		#finds restriction sites in length of DNA. First obtain the reverse complement using the appopriate method from this pacakge
		#returns tuple with location, sequence, and length. R is maximum length of restriction sites
		seq = self.seq
		complement = self.complement
		
		restrictionSites=[]
		for i in range(0, len(seq)):
			for n in range(4,R+1,2):
				s = seq[i:i+n]
				r = complement[i:i+n]
				r = r[::-1]
				if len(seq) -i < n: continue
				if s==r: restrictionSites.append((i+1, r, n))
				
		return restrictionSites

	def getConsensusSequence(self, sequences):
		#returns consensus sequences amongst set of sequences
		seq = self.seq
		# use list comprehension for create dictionary of nucleotides for each nucleotide position. Dictionary comprehension makes dictionary of nucleaotides
		# with counts set to zero for each position. 
		nuc = 'ATCG'
		#group nucleoties by position
		columns = zip(*sequences)
		
		#create profile of counts for each nucleotide for each position
		profile = [{n: columns[i].count(n) for n in nuc} for i in range(len(columns))]
		
		#generate consensus sequence based on the maximum frequency for each nucleatde
		consensus = ''.join([max(profile[i].keys(), key=(lambda k: profile[i][k])) for i in range(len(profile))])
		
		#profile = [ {n: 0 for n in nuc} for i in range(len(sequences[0]))]
		
		#tabulate numbers for each nucleotide along each position of the sequence
		#profile = [{base: [sequences[j][i] for j in range(len(sequences))].count(base) for base in 'ACTG'} for i in range(len(sequences[0]))]

		self.consensus = consensus #maybe should set up profile as attribute as well
		#print consensus
		
		#this just creates the profile table for the problem
		#prof = '\n'.join(["%s: %s\n" %(n, ' '.join([str(profile[i][n]) for i in range(len(seq))])) for n in nuc])
		#print prof
		return consensus,profile

	def transitionsTransversions(self, seq2):
		#given two sequences, calculates the ration of transitions to transversions
		seq1 = self.seq
		trans={'A':'G','G':'A','C':'T','T':'C'}
		mut = [0 if trans[a]==b else 1 for a, b in izip(seq1, seq2) if a!= b]
		r = mut.count(0)/float(mut.count(1))
		return r
		
		for i in range(len(seq1)):
			if seq1[i]==seq2[i]: continue
			if seq1[i] in PUR and seq2[i] in PUR: trans+=1
			elif seq1[i] in PYR and seq2[i] in PYR: trans+=1
			else: transv+=1
		
		return trans/float(transv)

	def getORF(self):
		#takes rna sequence and outputs list of protein sequences based on ORF from start to stop codons
		rna = self.seq
		
		f = open('geneticcode.dat','r').readlines()
		genCode = {i.split()[0]: i.split()[1] for i in f}
		#get locations of start codons
		
		PROT= []
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
						PROT.append(prot)
						break	
				
		return PROT

	def prot2mrna(self):
		f = open('./geneticcode.dat','r').readlines()
		revCode = {}
		N = 0
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
		#given a DNA motif, find the indices of nuclotides when the motif is not contiguous
		s = 0
		indices = []
		seq = self.seq
		for i in range(len(seq)):
			if s == len(mot): break
			if mot[s] == seq[i]:
				indices.append(i+1)
				s+=1
		return indices
		
#seq.transitionsTransversions(sequences[1])
f = open('data/rosalind_sseq.txt','r').readlines()
dat = map(lambda i: i.strip(), f)
fasta, sequences = fastaParse(dat)
sub = f[3].strip()
seq = sequence(sequences[0])
indices = seq.splicedMotif(sub)


indices = " ".join([str(i) for i in indices])
print(indices)
"""f = open('data/rosalind_gc.txt').readlines()
fasta, sequences = fastaParse(f)
fas=""
GC=0
for i in range(len(fasta)):
	seq = sequence(sequences[i].strip())
	gc= seq.getGCContent()
	if gc > GC: 
		fas = fasta[i]
		GC = gc
	
print "%s\n%.6f" %(fas,GC)"""
	#gc = 
"""f=open('rosalind_splc.txt','r').readlines()

fasta, sequences = fastaParse(f)

dna = sequences[0].strip()

introns=[]
for i in sequences[1:]:
	INT = i.strip()
	introns.append(INT)
seq = sequence(dna)
seq.spliceIntrons(introns)


dna = sequence(seq.mDNA)
dna.convertToRNA()
rna = sequence(dna.rna)

rna.getProteinSequence()"""

#dna = sequence(seq.mDNA)

#print seq.rna


#print prot

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


"""fasta, sequences = fastaParse(f)
introns=[]
for i in sequences[1:]:
	INT = i.strip()
	introns.append(INT)"""
	

"""dat = map(i.strip(), open('data/rosalind_revp.txt','r').readlines())

fasta, seqs = fastaParse(dat)
seq = sequence(seqs[0].strip())

seq.reverseComplement()
seq.findRestrictionSites(12)"""

#------------------------------------------------------------------------------#
#not included in the class. demonstrative
#------------------------------------------------------------------------------#

def Mendel(pops, P):
	#demonstrates Mendel's First Law. P is number of simulations to run.
	import random
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
	#print p/float(P)
	#to figure out mathmatically
	N = float(k + m + n)

	r_r = (n /N) * ((n - 1) / float(N - 1))
	h_h = (m /N) * ((m - 1) / float(N - 1))
	h_r = (m /N) * (n / (N - 1)) +  (n / N) * (m / (N - 1))

	recessive = r_r + h_h * 1/4.0 + h_r * 1/2.0
	dominant= 1 - recessive

#seq = 'AAAACCCGGT'

#pops = k,m,n
#Mendel(pops,10000)

def calculatingExpOffspring(f='rosalind_iev.txt'):
	#takes list of number of couple pairs with a given genotype: AA-AA, AA-Aa, AA-aa, Aa-Aa, Aa-aa, aa-aa and returns the expected number of those with the 
	#dom phenotype
	f = open(f,'r')
	couplePairs = [float(x) for x in f.readline().split()]

	prob = [1.0,1.0,1.0,0.75,0.5,0]

	p = sum([2 * a * b for a, b in itertools.izip(couplePairs, prob)])
	#p = 0
	#for i in range(len(couplePairs)):
	#    p += prob[i] * couplePairs[i] * 2 
	return p

def rabbitPairs(numMonths, numOffspring):
	#code for Rabbits and Recurrence Relations problem. Cannot seem to get it to solve the problem properly
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
	
#print rabbitPairs(5, 3)
#pairs = rabbitPairs(5, 3)
#print pairs
#use to test the GC code


#print(rna.seq)		
#nucleotideCount(dna)		
