import sys
import PyPDF2

infile  = sys.argv[1]
outfile = sys.argv[2]

pdf     = PyPDF2.PdfFileReader(infile)
page0   = pdf.getPage(0)
num     = 451/741
page0.scaleBy(num)  # float representing scale factor - this happens in-place
writer  = PyPDF2.PdfFileWriter()  # create a writer to save the updated results
writer.addPage(page0)

#result = str(arch) + "_chloroplot_images.pdf"
with open(outfile, "wb+") as f:
    writer.write(f)
