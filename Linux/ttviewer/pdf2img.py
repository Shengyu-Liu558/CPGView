import sys, fitz, os, datetime

def pyMuPDF_fitz(pdfPath, imagePath):
    startTime_pdf2img = datetime.datetime.now()#开始时间

    print("imagePath="+imagePath)
    pdfDoc = fitz.open(pdfPath)
    for pg in range(pdfDoc.pageCount):
        page   = pdfDoc[pg]
        rotate = int(0)
        
        
        zoom_x = 2 #(1.33333333-->1056x816)   (2-->1584x1224)
        zoom_y = 2
        mat    = fitz.Matrix(zoom_x, zoom_y).preRotate(rotate)
        pix    = page.getPixmap(matrix=mat, alpha=False)

        if not os.path.exists(imagePath):
            os.makedirs(imagePath) 

        pix.writePNG(imagePath+'/'+'images_%s.png' % pg)

    endTime_pdf2img = datetime.datetime.now()
    print('pdf2img time=',(endTime_pdf2img - startTime_pdf2img).seconds)

if __name__ == "__main__":
    pdfPath   = sys.argv[1]; #'../path/demo.pdf'
    imagePath = sys.argv[2]
    #pdfPath = '../path/demo.pdf'
    #imagePath = '../path/image'
    pyMuPDF_fitz(pdfPath, imagePath)
