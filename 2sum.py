#!/usr/bin/env python

"""
given an array of integers, if A[i] = -A[j], return indices i, j. else, return -1
"""

def TwoSUM(dat):
	p = -1
	for i, v in enumerate(dat):
	 	if -v in dat:
			p = i + 1, dat.index(-v) + 1
			return p
	
	return p

if __name__ == "__main__":
	f = open('data/rosalind_2sum.txt','r').readlines()
	o = open("ros_2sum.txt",'w')
	m, n = map(int, f[0].split())
	dat = [list(map(int, f[i].split())) for i in range (1, m+1)]
	for d in dat:
		result =TwoSUM(d)
		if type(result) == int: print "%d" % result	
		else: print "%s" % " ".join([str(j) for j in result])
