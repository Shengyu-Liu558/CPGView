import os
import sys
from PIL import Image
###########################################
print ("Usage:\t inputImageFile diameter outputImageFile")

inputImageFile = sys.argv[1] #inputImageFile
diameter = int(sys.argv[2]) #diameter, integer
outputImageFile = sys.argv[3] #outputImageFile

###########################################

im = Image.open(inputImageFile)
region = im.resize((diameter, diameter),Image.ANTIALIAS)
region.save(outputImageFile)

