#ze python

#######################

import os,sys

from Bio import Phylo
from Bio import SeqIO
from Bio import AlignIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.Applications import PhymlCommandline

#######################

args=sys.argv
inp_dir=args[1]

#######################

def fasta_getter(file_):
	for line in open(file_):
		line=line.strip()
		if line.startswith('>'):
			try: yield contig
			except: pass
			contig=''
		contig=contig+line+'\n'
	yield contig

def make_malign(fasta):
	align='%s.aln' %fasta
	dendro='%s.dnd' %fasta
	cline='clustalw -infile=%s' %fasta
	os.system(cline)
	return align,dendro 

def get_genome(id_):
	genome=id_.split('__')[1]
	return genome

#######################

# make multi align

for fasta in os.listdir(inp_dir):
	if not fasta.endswith('.fasta'): continue
	try:make_malign(fasta)
	except: print 'couldnt align %s' %fasta

# wrap it
alignments=[AlignIO.read(f,'clustal')
				for f in os.listdir('.')
					if f.endswith('.aln')]

# sort it
for a in alignments: a.sort(key=lambda x:x.id)
labels=[r.id for r in a]
 
# concatenate it
concatenamers=[SeqRecord(Seq(''),id=l) for l in labels]
for a in alignments:
	for i in range(len(a)):
		concatenamers[i].seq = concatenamers[i].seq + a[i].seq

# get concatenamer aln
align=MultipleSeqAlignment(concatenamers)

# write the concatenamer file
AlignIO.write(align,'concatenamer.phy','phylip-relaxed')

# make a tree
## make tree file
cmdline = PhymlCommandline(input='concatenamer.phy', alpha='e', bootstrap=100)
out_log, err_log = cmdline()
## read it / draw it
egfr_tree = Phylo.read("concatenamer.phy_phyml_tree.txt", "newick")
Phylo.draw_ascii(egfr_tree)


