from Bio import Entrez
from Bio import SeqIO,SeqRecord
from BioWrapper import *
from SilicoPCR import *
import os,sys

###################################

if __name__ == '__main__':

	args=sys.argv

	alleles_dir=args[1]
	genomes_dir=args[2]
	out_name=args[3]

###################################

if __name__ == '__main__':

	ML=Multi_locus(alleles_dir)
	out=open(out_name,'w')
	#d={l.name:l.alleles for l in ML.loci}
	d={}
	for l in ML.loci: d[l.name]=l.alleles
	for genome in os.listdir(genomes_dir):
		to_write=genome+'\t'
		f=genomes_dir+genome+'/MLST.fasta'
		for s in SeqIO.parse(open(f),'fasta'):
			if s.seq.tostring() =='':
				to_write+='no amplification\t'
				continue
			allele=matching_allele(s.seq.tostring(),d[s.name])
			to_write+=allele.id+'\t'
		out.write(to_write + '\n')
	out.close()
	print 'task done!'
