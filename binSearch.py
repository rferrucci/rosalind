#!/usr/bin/env python

#finds location in sorted list of given item

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

l = ""

f = open('data/rosalind_bins.txt','r').readlines();
A = [int(a) for a in f[2].split()]
B = [int(b) for b in f[3].split()]


for b in range(len(B)):
	l += '%d ' % (binSearch(A, B[b]))

print l
