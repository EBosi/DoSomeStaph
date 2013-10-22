#!/usr/bin/python
from Bio import SeqIO
from Bio import Entrez

def get_from_genbank(accession,mail='emanuele.bosi@unifi.it'):
    Entrez.email=mail
    handle = Entrez.efetch(db="nucleotide",
                        id=accession,
                        rettype="gbwithparts")
    record = SeqIO.read(handle,"genbank")
    handle.close()
    return record
    
