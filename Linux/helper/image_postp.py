import sys
from fitz import fitz, Rect
import PyPDF2

def add_footer(pdf):
    img  = open(infile2, "rb").read()
    rect = Rect(298, 247, 444, 392)

    for i in range(0, pdf.pageCount):
        page = pdf[i]
        page.insertImage(rect, stream=img)
 
if __name__ == "__main__":

    if len(sys.argv) != 5:
        print ("Usage: "+sys.argv[0]+" infile1 infile2 outfile1 outfile2")
        exit()
    infile1  = sys.argv[1]
    infile2  = sys.argv[2]
    outfile1 = sys.argv[3]
    outfile2 = sys.argv[4]
    doc      = fitz.open(infile1)
    
    add_footer(doc)
    doc.save(outfile1)

    pdf     = PyPDF2.PdfFileReader(outfile1)
    page0   = pdf.getPage(0)
    num     = 451/741
    page0.scaleBy(num)  # float representing scale factor - this happens in-place
    writer  = PyPDF2.PdfFileWriter()  # create a writer to save the updated results
    writer.addPage(page0)

    #result = str(arch) + "_chloroplot_images.pdf"
    with open(outfile2, "wb+") as f:
        writer.write(f)

