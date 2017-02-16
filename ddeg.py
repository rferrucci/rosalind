#!/usr/bin/env python

#returns sum of degrees of each nodes neighbors 

from collections import defaultdict

def doubleDegreeArray(n,dat):
	ad_dict = defaultdict(list)
	
	for i in dat:
		ad_dict[i[0]].append(i[1])
		ad_dict[i[1]].append(i[0])
	
	edges = [0 for i in range(n+1)]

	for i in range(1, len(edges)):
		for j in ad_dict[i]:
			edges[i] += len(ad_dict[j])

	return edges


if __name__ == '__main__':
	f = open('data/rosalind_ddeg.txt','r').readlines()
	n = int (f[0].split()[0])

	dat = [list(map(int, f[i].split())) for i in range (1, len(f))]

	results  = doubleDegreeArray(n,dat)
	results = " ".join([str(results[i]) for i in range(1, len(results))])
	print results
