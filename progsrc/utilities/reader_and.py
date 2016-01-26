import numpy
import math
import sys
import re

from scipy import stats

####################################################################

def file_len(fname):
  
  with open(fname) as f:
    for i, l in enumerate(f):
      pass

  return i + 1

####################################################################

filename = "noname.txt"

if (len(sys.argv) == 1):
  print >> sys.stderr, "usage: ", sys.argv[0], " filename.txt"
  exit(1)
else:
  filename = sys.argv[1]

numofline = file_len(filename)

fp = open(filename, "r")

# jump first line
fp.readline()

meanval1 = []
meanval2 = []
meanval3 = []
meanval4 = []
meanval5 = []

for i in range(numofline):
  l = fp.readline()

  if not l: 
    break

  p = re.compile(r'\s+')
  line = p.sub(' ', l)
  line = line.lstrip()
  line = line.rstrip()
           
  plist = line.split(" ")

  numof = int(plist[1])

  layersids = ""
  layerold = ""

  xold = 0.0
  yold = 0.0
  zold = 0.0

  val1 = 0.0
  val2 = 0.0
  val3 = 0.0
  val4 = 0.0 
  val5 = 0.0 

  for j in range(numof):
    coordline = fp.readline()

    coordline = p.sub(' ', coordline)
    coordline = coordline.lstrip()
    coordline = coordline.rstrip()
           
    coordlinelist = coordline.split(" ")

    layersids += coordlinelist[3]

    pid = int(coordlinelist[7])

    if (pid > 0):
      charge = 1.0
    else:
      charge = -1.0

    xi = float(coordlinelist[0])
    yi = float(coordlinelist[1])
    zi = float(coordlinelist[2])

    if (j > 0):
      dist = math.sqrt((xi-xold)**2 + (yi-yold)**2 + (zi-zold)**2)

      if j == 1:
        val1 = dist
      elif j == 2:
        val2 = dist
      elif j == 3:
        val3 = dist
      elif j == 4:
        val4 = dist
      elif j == 5:
        val5 = dist

    layerold = coordlinelist[3]

    xold = xi
    yold = yi
    zold = zi

  paramline = fp.readline()
  
  # quick check for layers id
  if (layersids == "5678910"):
    meanval1.append(val1)
    meanval2.append(val2)
    meanval3.append(val3)
    meanval4.append(val4)
    meanval5.append(val5)

print " 5  6 ", numpy.mean(meanval1), " ", numpy.std(meanval1)
print " 6  7 ", numpy.mean(meanval2), " ", numpy.std(meanval2)
print " 7  8 ", numpy.mean(meanval3), " ", numpy.std(meanval3)
print " 8  9 ", numpy.mean(meanval4), " ", numpy.std(meanval4)
print " 9 10 ", numpy.mean(meanval5), " ", numpy.std(meanval5)

fp.close()
