import os,sys
from SilicoPCR import *
from Bio import Entrez
from Bio import SeqIO,SeqRecord
from BioWrapper import *
from MLST import *

######################

if __name__ == '__main__':
	usage=' python MLST.py ALLELES_DIR GENOME '

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

	args = sys.argv
	if len(args) != 3:
		print usage
		sys.exit()
	inp_dir = args[1]
	inp_genome = args[2]
	out_file = inp_genome + 'MLST.fasta'


######################

if __name__=='__main__':
	ML = Multi_locus(inp_dir)
	ML.add_primers(primers)
	#genomes=[g for g in os.listdir(inp_genomes_dir)]
	#
	genome_name=inp_genome.split('/')[-1]
	genome_dir='/'.join(inp_genome.split('/')[:-1]) + '/'
	mlst_profile=MLST_profile(genome_name)
	accessions=get_accessions(genome_name,genome_dir)
	genome_seqs=[s.seq.tostring() for s in get_seqs(accessions)]
	for locus in ML.loci:
		PCR=PCRexperiment(locus.name,locus.primers,accessions,genome_seqs)
		amplicon=PCR.amplify()
		sequence=PCR.trim(amplicon,locus.borders,locus.marker_length)
		if sequence == None: sequence=''
		mlst_profile.add_sequence(locus.name,sequence)
	mlst_profile.write_fasta(key_order,out_file)
	#
	print 'task done'
