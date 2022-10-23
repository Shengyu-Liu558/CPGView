import sys
import codecs
import os
import PyPDF2 as PyPDF2

arch1 = sys.argv[1]
arch2 = sys.argv[2]
files = list()
for filename in os.listdir(arch1):

    if filename.endswith(".pdf"):
        files.append(filename)

newfiles = sorted(files, key=lambda d: str(d.split(".pdf")[0]))

os.chdir(arch1)

pdfwriter = PyPDF2.PdfFileWriter()
for item in newfiles:
    pdfreader = PyPDF2.PdfFileReader(open(item, "rb"))
    for page in range(pdfreader.numPages):
        pdfwriter.addPage(pdfreader.getPage(page))
result = str(arch2) + "_all_images.pdf"
with codecs.open(result, "wb") as f:
    pdfwriter.write(f)
