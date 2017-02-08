#!/usr/bin/env python

"""	merges two sorted arrays into a second sorted array
"""

def mergeSort(A, B):
	C = []
	i, j = 0, 0
	while len(A) > 0 and len(B) > 0:
		if A[0] < B[0]:
			C.append(A.pop(0))
		else:
			C.append(B.pop(0))
			
	C.extend(A)
	C.extend(B)

	return C
		
if __name__ == "__main__":
	f = open('data/rosalind_mer.txt','r').readlines()
	A = list(map(int, f[1].split()))
	B = list(map(int, f[3].split()))
	C = mergeSort(A, B)

	print " ".join([str(c) for c in C])
