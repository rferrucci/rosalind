#!/usr/bin/env python

"""	given n vertices, find the degree of each, that is how many neighbors each one has, from 
	from the list of edge pairs in dat
"""

def degreeArray(n, dat):
	deg = [0 for i in range(n)]

	for i in range(0, len(dat)):
		deg[dat[i][0]-1] += 1
		deg[dat[i][1]-1] += 1

	return " ".join([str(i) for i in deg])
		
if __name__ == "__main__":
	f = open("data/rosalind_deg.txt",'r').readlines()
	#dat = [(int(f[i].split()[0].strip()), int(f[i].split()[1].strip())) for i in range(1,len(f))]
	
	n,m = map(int,f[0].split())
	edges = [map(int,f[i].split()) for i in xrange(1,m+1)]
	#size = int(f[0].split()[0])
	
	#print(edges)
	deg = degreeArray(n, dat)
	print (deg)
	#o = open('degree.txt','w')
	#o.write(deg)
	#o.close()
