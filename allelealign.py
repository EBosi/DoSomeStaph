import os,sys
from BioWrapper import rev_complement

def pairwise_alleles(s1,s2,mm=1):
	N=0
	if len(s1) != len(s2):
		return False,'no same length'
	for i in range(len(s1)):
		if s1[i] != s2[i]: N+=1
		if N>mm: break
	if N<=mm: return True,N
	N=0
	s1=rev_complement(s1)
	for i in range(len(s1)):
		if s1[i] != s2[i]: N+=1
		if N>mm: break
	return False,N
	
