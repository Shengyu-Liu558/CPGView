from Bio import Entrez
from Bio import SeqIO
import os, sys

acc      = sys.argv[1]
file_gb  = sys.argv[2]
#file_gb = outdir + "/" + acc + '.gb'
Entrez.email = "cliu6688@hotmail.com"
hd1 = Entrez.efetch(db="nucleotide",id=[acc],rettype='gb')
seq = SeqIO.read(hd1,'gb')
fw = open(file_gb,'w')
SeqIO.write(seq,fw,'gb')
fw.close()
os.getcwd()
