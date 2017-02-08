#!/usr/bin/env python

def fib(N):
	#return list of Fibonacci Numbers given n
	F = []
	for n in range(N+1):
		if n == 0:
		    F.append(0)
		elif n == 1:
			F.append(1)
		else:
			f = F[n-1] + F[n-2]
			F.append(f)
	return F
	
def binSearch(alist, item):
    #from list alist, find location of item
    first = 0
    last = len(alist)-1
    found = False

    while first<=last and not found:
        midpoint = (first + last)//2
        if alist[midpoint] == item:
            found = True
        else:
            if item < alist[midpoint]:
                last = midpoint-1
            else:
                first = midpoint+1

    if found == True: return midpoint + 1 	#if item found, return location
    else: return -1							# else, return -1

def degreeArray(n, dat):
	"""	given n vertices, find the degree of each, that is how many neighbors each one has, from 
		from the list of edge pairs in dat
	"""
	deg = [0 for i in range(size)]

	for i in range(0, len(dat)):
		deg[dat[i][0]-1] += 1
		deg[dat[i][1]-1] += 1

	return " ".join([str(i) for i in deg])
		
def insertionSort(A):
	#sorts an array A using the insertion sort algorithim, returns number of swaps requried
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

def majorityElement(n, A):
	#from a list of arrays, determines prescence of majority elements, 
	Maj = []
	for a in A:
		maj = -1
		for i in set(a):
			if a.count(i) > n/2.0:
				maj = i
		Maj.append(maj)
	return Maj

def mergeSort(A, B):
	#merges two sorted arrays into a second sorted array
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
f = open('data/rosalind_mer.txt','r').readlines()
A = list(map(int, f[1].split()))
B = list(map(int, f[3].split()))
C = mergeSort(A, B)

print " ".join([str(c) for c in C]) 
