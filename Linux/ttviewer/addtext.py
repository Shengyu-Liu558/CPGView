import sys
from fitz import fitz
from Bio import SeqIO

arch1 = sys.argv[1]
name = sys.argv[1]

for rec in SeqIO.parse(arch1,'gb'):
    m = (rec.annotations)["organism"]
    
name1 = str(str(name).replace(".gbf",""))

doc = fitz.open(name1 + '_circle0.pdf') 

def add_footer(pdf):
    for i in range(0, pdf.pageCount):
        page = pdf[i]
        page.insertText((15,50),m,fontname="Times-Italic",fontsize=15)


add_footer(doc)
result = name1 + '_circle.pdf'
doc.save(result)