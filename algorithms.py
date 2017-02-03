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
