#!/usr/bin/env python

#sorts an array A using the insertion sort algorithim, returns number of swaps requried

def insertionSort(A):
	
	n = 0
	for i in range(1, len(A)):
		k = i
		while k > 0 and A[k] < A[k-1]:
			tmp = A[k-1]
			A[k-1] = A[k]
			A[k] = tmp
			k -= 1
			n +=1
	return n

if __name__ == "__main__":
	f = open('data/rosalind_ins.txt','r').readlines()
	dat = [int(i) for i in f[1].split()]
	n = insertionSort(dat)
	print n
