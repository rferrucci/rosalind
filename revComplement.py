#!/usr/bin/env python
from rosalind import *

#!/usr/bin/env python
from rosalind import *

"""Solution to the Complementing a Strand of DNA problem in the Bioinformatics Stronghold section of Rosalind
location: http://rosalind.info/problems/revc/
Problem: return the reverse complement of a dna sequence
"""

# return the reverse complement of a strand of DNA
if __name__ == "__main__":	
	f = open("data/rosalind_cmp.txt",'r').readlines()
	fasta, seqs = fastaParse(f)
	seq = sequence(seqs[0])
	cmp = seq.reverseComplement()
	print(cmp)
