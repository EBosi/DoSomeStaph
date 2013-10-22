import os,sys
from SilicoPCR import *
from Bio import Entrez
from Bio import SeqIO,SeqRecord
from Bio.Seq import Seq
from BioWrapper import *

######################

if __name__ == '__main__':
	usage=' python MLST.py ALLELES_DIR GENOMES_DIR OUTPUT '
	primers={
		'arcc':['TTGATTCACCAGCGCGTATTGTC','AGGTATCTGCTTCAATCAGCG']  ,
		'aroe':['ATCGGAAATCCTATTTCACATTC','GGTGTTGTATTAATAACGATATC'],
		'glpf':['CTAGGAACTGCAATCTTAATCC' ,'TGGTAAAATCGCATGTCCAATTC'],
		'gmk':['ATCGTTTTATCGGGACCATC'   ,'TCATTAACTACAACGTAATCGTA'] ,
		'pta':['GTTAAAATCGTATTACCTGAAGG','GACCCTTTTGTTGAAAAGCTTAA'] ,
		'tpi':['TCGTTCATTCTGAACGTCGTGAA','TTTGCACCTTCTAACAATTGTAC'] ,
		'yqil':['CAGCATACAGGACACCTATTGGC','CGTTGAGGAATCGATACTGGAAC']}

	key_order=['arcc',
		'aroe',
		'glpf',
		'gmk',
		'pta',
		'tpi',
		'yqil']

######################

if __name__ == '__main__':
	args = sys.argv
	inp_dir = args[1]
	inp_genomes_dir = args[2]
	out_file = args[3]


######################

def get_accessions(genome,dir_):
	dir_=dir_+genome + '/'
	accessions=[dir_ + i for i in os.listdir(dir_) if i.endswith('.fna')]
	return accessions
	
def get_seqs(files):
	seqs=[seq for f_ in files for seq in SeqIO.parse(open(f_),'fasta')]
	return seqs

######################

class MLST_profile(object):
	def __init__(self,name):
		self.name=name
		self.profile={}
		self.sequences={}
	def add_allele(self,locus_name,allele):
		self.profile[locus_name]=allele
	def add_sequence(self,locus_name,seq):
		self.sequences[locus_name]=seq
	def to_write(self,key_order):
		out=self.name
		for k in key_order:
			out= out + '\t' + self.profile[k]
		out += '\n'
		return out
	def write_fasta(self,key_order,out=None):
		if out == None: out=self.name + '_trimmed_amplicons'
		out=open(out,'w')	
		records=[SeqRecord.SeqRecord(Seq(self.sequences.get(k,'')),id=k,description='') for k in key_order]
		SeqIO.write(records,out,'fasta')

######################

if __name__ == '__main__':
	ML = Multi_locus(inp_dir)
	ML.add_primers(primers)
	genomes=[g for g in os.listdir(inp_genomes_dir)]
	out=open(out_file,'w')
	#
	for genome in genomes:
		mlst_profile=MLST_profile(genome)
		accessions=get_accessions(genome,inp_genomes_dir)
		genome_seqs=[s.seq.tostring() for s in get_seqs(accessions)]
		for locus in ML.loci:
			PCR=PCRexperiment(locus.name,locus.primers,accessions,genome_seqs)
			locus_name,allele=PCR.identify_allele(locus.borders,locus.alleles,locus.marker_length)
			if allele == None:
				allele_name='no_locus_found'
			else:
				allele_name=allele.id
			mlst_profile.add_allele(locus_name,allele_name)
		out.write(mlst_profile.to_write(key_order))
	#
	print 'task done'
