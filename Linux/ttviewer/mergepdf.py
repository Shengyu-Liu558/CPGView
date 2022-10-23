import sys

file1 = sys.argv[1]
file2 = sys.argv[2]
file3 = sys.argv[3]


from PyPDF2 import PdfFileReader, PdfFileMerger
pdf_file1 = PdfFileReader(file1)
pdf_file2 = PdfFileReader(file2)
output = PdfFileMerger()
output.append(pdf_file1)
output.append(pdf_file2)
output.write(file3)
