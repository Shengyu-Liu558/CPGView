import sys
from fitz import fitz, Rect

def add_footer(pdf):
    img  = open(infile2, "rb").read()
    rect = Rect(298, 247, 444, 392)

    for i in range(0, pdf.pageCount):
        page = pdf[i]
        page.insertImage(rect, stream=img)
 
if __name__ == "__main__":

    if len(sys.argv) != 4:
        print ("Usage: "+sys.argv[0]+" infile1 infile2 outfile")
        exit()
    infile1 = sys.argv[1]
    infile2 = sys.argv[2]
    outfile = sys.argv[3]
    doc = fitz.open(infile1)
    
    add_footer(doc)
    #result = str(outfile) + '/new_chloroplot0.pdf'
    doc.save(outfile)
