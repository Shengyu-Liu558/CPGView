import sys
print ("Usage:\toutRingMap innerRingMap InnterRingTmpFile oiRingMapFile")

basemapFile = sys.argv[1] #
upmapFile   = sys.argv[2] #
roundFile   = sys.argv[3] #
outmapFile  = sys.argv[4] #

import os, math 
from PIL import Image

def circle(ImageFile):
    
 ima = Image.open(ImageFile).convert("RGBA")
 #ima = ima.resize((500,500),Image.ANTIALIAS)
 ima = ima.resize((3428,3428),Image.ANTIALIAS)
 
 # ima = ima.resize((600, 600), Image.ANTIALIAS) 
 size = ima.size 
 #print(size) 
  
 r2 = min(size[0], size[1])
 if size[0] != size[1]: 
  ima = ima.resize((r2, r2), Image.ANTIALIAS)
  #ima.show()
  
 r3 = 800
 imb = Image.new('RGBA', (r3*2, r3*2),(255,255,255,0))
 
 pima = ima.load()  
 pimb = imb.load()
 
 r = float(r2/2) 
  
 for i in range(r2): 
  for j in range(r2): 
   lx = abs(i-r) 
   ly = abs(j-r)
   l = (pow(lx,2) + pow(ly,2))** 0.5 #
  
   if l < r3: 
    pimb[i-(r-r3),j-(r-r3)] = pima[i,j]
    
 imb.save(roundFile)
    

####Start of Main Program

circle(upmapFile)

#region = Image.open("D:/test_circle.png")
region = Image.open(roundFile)

print (region.size)

from PIL import Image
base_img = Image.open(basemapFile)

base_img = base_img.convert("RGBA")

print (base_img.size)

#base_img.show()
print ("--------------------------------------------------")


target = Image.new('RGBA', base_img.size, (0, 0, 0, 0))


r_centre_point = int(base_img.size[0]/2)

print (r_centre_point)


radius = 800

left = r_centre_point-int(math.sqrt((radius**2)/2))
upper = r_centre_point-int(math.sqrt((radius**2)/2))
right = r_centre_point+int(math.sqrt((radius**2)/2))
lower = r_centre_point+int(math.sqrt((radius**2)/2))

box = (left, upper, right, lower)

print (box)

region = region.convert("RGBA")
target.paste(base_img,(0,0),base_img)

target.paste(region,(r_centre_point-radius,r_centre_point-radius), mask=region)
target.save(outmapFile)


