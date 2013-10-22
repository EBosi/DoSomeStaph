import os
from Bio import Entrez
from Bio import SeqIO,SeqRecord
from BioWrapper import *


#######################

def get_amplicon(pr1,pr2,seq,stringency=1):
	""" get putative amplified sequence between primers.
		stringency is a constraint for sequence similarity search
		"""	
	# check if there is aspecific annealing
	if ((seq.count(pr1) > 1) or
		(seq.count(rev_complement(pr1)) > 1) or
		(seq.count(pr2) > 1) or
		(seq.count(rev_complement(pr2)) > 1)):
		return ''
		#
	# if some primers don't anneal, allow some gaps (1)
	if	((seq.count(pr1) == 0) and
		(seq.count(rev_complement(pr1)) == 0)):
			# if we know the other primer orientation, this is easier...
			if seq.count(pr2) == 1:
				try: new_pr1=rev_complement(find_RE(rev_complement(pr1),seq))
				except: pass
			if seq.count(rev_complement(pr2)) == 1:
				new_pr1=find_RE(pr1,seq)
			else:
				# if not, try all orientation
				new_pr1=find_RE(pr1,seq)
				if new_pr1==None:
					new_pr1=find_RE(rev_complement(pr1),seq)
			if new_pr1==None: return ''
			else: pr1=new_pr1

	# same			
	if	((seq.count(pr2) == 0) and
		(seq.count(rev_complement(pr2)) == 0)):
			if seq.count(pr1) == 1:
				try: pr2=rev_complement(find_RE(rev_complement(pr2),seq))
				except: pass
			if seq.count(rev_complement(pr1)) == 1:
				pr2=find_RE(pr2,seq)
			if pr2 == None: return ''
		#
	# if there are not problems, return the amplicon
	if seq.find(pr1) != -1:
		if seq.find(rev_complement(pr2)) != -1:
			from_ = seq.find(pr1) + len(pr1)
			to_ = seq.find(rev_complement(pr2))
			return seq[from_ : to_]
	# same
	if seq.find(rev_complement(pr1)) != -1:
		if seq.find(pr2) != -1:
			from_ = seq.find(pr2) + len(pr2)
			to_ = seq.find(rev_complement(pr1))
			return seq[from_ : to_]
	#
	return ''

def find_RE(pr,seq,mm=1):
	""" just in case the primer don't anneal 100%, allow some mismatches (1) """
	import re
	from itertools import combinations
	seqStr = pr
	if mm == 1: searchSeqREStr = seqStr + '|' + \
	'|'.join(seqStr[:i]+"[ACTGN]".replace(c,'') +seqStr[i+1:] 
			 for i,c in enumerate(seqStr))
	elif mm == 2: searchSeqREStr = seqStr + '|' + \
		'|'.join([seqStr[:i]+"[ACTGN]"+seqStr[i+1:j]+"[ACTGN]"+seqStr[j+1:]
			for i,j in combinations(range(len(seqStr)),2)])
	searchSeqRE=re.compile(searchSeqREStr)
	matches=[match
	for match in searchSeqRE.finditer(seq)]
	if len(matches) > 1: print 'aspecific amplification! (RE)'
	if len(matches) == 0:
		searchSeqREStr = seqStr + '|' + \
		'|'.join([seqStr[:i]+"[ACTGN]"+seqStr[i+1:j]+"[ACTGN]"+seqStr[j+1:]
			for i,j in combinations(range(len(seqStr)),2)])
		searchSeqRE=re.compile(searchSeqREStr)
		matches=[match
			for match in searchSeqRE.finditer(seq)]
	try: return matches[0].group(0)
	except: return None



def try_trimmings(seq,lefts,rights,final_length,border_length):
	good_trimmed=[]
	for l in lefts:
		for r in rights:
			from_= min(l,r)
			to   = max(l,r) + border_length
			trimmed=seq[from_:to]
			if len(trimmed) == final_length:
				good_trimmed.append(trimmed)
	return good_trimmed


def trim_amplicon(seq,borders,marker_length):
	""" given the (conserved) alleles borders,
		the amplicon is trimmed """
	if seq == None: return None
	border_l,border_r=borders[0],borders[1]
	from_,to=None,None
	left_matches=([seq.find(l) for l in border_l if seq.find(l) != -1] + [seq.find(rev_complement(l)) for l in border_l if seq.find(rev_complement(l)) != -1])
	right_matches=([seq.find(r) for r in border_r if seq.find(r) != -1] + [seq.find(rev_complement(r)) for r in border_r if seq.find(rev_complement(r)) != -1])	
	#
	border_length=len(border_r.copy().pop())
	if len(left_matches) == 1 and len(right_matches) == 1:
		from_=min(left_matches[0],right_matches[0])
		to=max(left_matches[0],right_matches[0]) + border_length 
		return seq[from_:to]
	if len(left_matches) > 1 or len(right_matches) > 1:
		possible_trimmings=try_trimmings(seq,left_matches,right_matches,marker_length,border_length)
		if len(possible_trimmings)!= 1:
			print "multiple possible trimming!"
		else:
			return possible_trimmings[0]
	print 'trimming failed!'
	return 

def matching_allele(query,alleles):
	for a in alleles:
		seq=a.seq.tostring()
		if (query == seq) or (rev_complement(query) == seq):
			return a
	#print "new allele!"
	return SeqRecord.SeqRecord(query,id='new allele')

def get_alleles(file_):
	alleles=[s.seq.tostring()
			for s in SeqIO.parse(file_,"fasta") ]
	return alleles

def get_borders(alleles,n=15):
	lefts,rights=set(),set()
	for a in alleles:
		seq=a.seq.tostring()
		lefts.add(seq[:n])
		rights.add(seq[-n:])
	out=[lefts,rights]
	return out

#######################

class PCRexperiment(object):
	""" object for in silico PCR  experiment:
		as input it takes the region to amplify (locus name)
		and the primers (3' to 5'). It gives back the amplicon.
		You may want to trim the amplicon to get the exact marker
		for a MLST analysis."""
	#
	def __init__(self,locus_names,primers,accessions,seqs=['']):
		self.name=locus_names
		self.primer1=primers[0]
		self.primer2=primers[1]
		self.accession=accessions
		self.seqs=seqs
		#
	def amplify(self):
		all_amplicons= [get_amplicon(self.primer1,self.primer2,s) for s in self.seqs]
		amplicons=[a for a in all_amplicons if a != '']
		if len(amplicons) == 1 :
			return amplicons[0]
		if len(amplicons) > 1:
			print 'aspecific amplicons!'
		if len(amplicons) == 0:
			print 'no amplification!'
		return
		#
	def trim(self,amplicon,borders,marker):
		return trim_amplicon(amplicon,borders,marker)
		#
	def identify_allele(self,borders,alleles,marker_length):
		amplicon = self.amplify()
		if amplicon == None:
			return self.name,amplicon
		query_amplicon = self.trim(amplicon,borders,marker_length)
		if query_amplicon == None:
			return self.name,query_amplicon
		return self.name , matching_allele(query_amplicon,alleles)

class Locus_infos(object):
	""" Locus object modelling information for MLST analysis, such as:
		- allele sequences
		- primers
		- compute conserved borders for amplicon trimming
	"""
	def __init__(self,alleles_file,primers=''):
		self.file_=alleles_file
		self.name=self.file_.split('/')[-1].split('_allele')[0]
		self.primers=primers
		self.alleles=[ a for a in SeqIO.parse(open(alleles_file),'fasta')]
		self.marker_length=len(self.alleles[0])
		self.borders=get_borders(self.alleles)
	def add_primers(self,primers):
		self.primers=primers

class Multi_locus(object):
	""" All-in-one multiple locus objects """
	def __init__(self,dir_):
		self.dir_=dir_
		self.files=[dir_ + f
					for f in os.listdir(dir_)
						if f.endswith('_alleles')]
		self.loci=[Locus_infos(l) for l in self.files]
	def add_primers(self,primers_diz):
		for locus in self.loci:
			primers=primers_diz[locus.name]
			locus.add_primers(primers)

# To implement: - circular genome in amplification constraint
# 				- stringency constraint
