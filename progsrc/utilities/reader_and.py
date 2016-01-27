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

rzmeanval1 = []
rzmeanval2 = []
rzmeanval3 = []
rzmeanval4 = []
rzmeanval5 = []

rphimeanval1 = []
rphimeanval2 = []
rphimeanval3 = []
rphimeanval4 = []
rphimeanval5 = []

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
  rold = 0.0 
  phiold = 0.0

  val1 = 0.0
  val2 = 0.0
  val3 = 0.0
  val4 = 0.0 
  val5 = 0.0 

  rzval1 = 0.0
  rzval2 = 0.0
  rzval3 = 0.0
  rzval4 = 0.0 
  rzval5 = 0.0 

  rphival1 = 0.0
  rphival2 = 0.0
  rphival3 = 0.0
  rphival4 = 0.0 
  rphival5 = 0.0 

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
    ri = math.sqrt(math.pow(xi, 2.0) + math.pow (yi, 2.0))
    phii = math.acos(xi/ri)

    if (j > 0):
      dist = math.sqrt((xi-xold)**2 + (yi-yold)**2 + (zi-zold)**2)
      distrz = math.sqrt((ri-rold)**2 + (zi-zold)**2)
      distrphi = math.sqrt((ri-rold)**2 + (phii-phiold)**2)

      if j == 1:
        val1 = dist
        rzval1 = distrz 
        rphival1 = distrphi
      elif j == 2:
        val2 = dist
        rzval2 = distrz 
        rphival2 = distrphi
      elif j == 3:
        val3 = dist
        rzval3 = distrz 
        rphival3 = distrphi
      elif j == 4:
        val4 = dist
        rzval4 = distrz 
        rphival4 = distrphi
      elif j == 5:
        val5 = dist
        rzval5 = distrz 
        rphival5 = distrphi

    layerold = coordlinelist[3]

    xold = xi
    yold = yi
    zold = zi
    rold = ri
    phiold = phii

  paramline = fp.readline()

  paramline = p.sub(' ', paramline)
  paramline = paramline.lstrip()
  paramline = paramline.rstrip()
  
  paramlinelist = paramline.split(" ")
  
  pt = float(paramlinelist[0])
  phi = float(paramlinelist[1])
  eta = float(paramlinelist[3])
  z0 = float(paramlinelist[4])
  theta = 2.0 * math.atan (math.exp(-eta))
  pz = pt * math.cos(theta)
  
  # quick check for layers id
  if (layersids == "5678910"):
    meanval1.append(val1)
    meanval2.append(val2)
    meanval3.append(val3)
    meanval4.append(val4)
    meanval5.append(val5)

    rzmeanval1.append(val1)
    rzmeanval2.append(val2)
    rzmeanval3.append(val3)
    rzmeanval4.append(val4)
    rzmeanval5.append(val5)

    rphimeanval1.append(val1)
    rphimeanval2.append(val2)
    rphimeanval3.append(val3)
    rphimeanval4.append(val4)
    rphimeanval5.append(val5)

print "Using ", meanval1.size() , " events "

print " XYZ "
print " 5  6 ", numpy.mean(meanval1), " ", numpy.std(meanval1)
print " 6  7 ", numpy.mean(meanval2), " ", numpy.std(meanval2)
print " 7  8 ", numpy.mean(meanval3), " ", numpy.std(meanval3)
print " 8  9 ", numpy.mean(meanval4), " ", numpy.std(meanval4)
print " 9 10 ", numpy.mean(meanval5), " ", numpy.std(meanval5)
print " RZ "
print " 5  6 ", numpy.mean(rzmeanval1), " ", numpy.std(rzmeanval1)
print " 6  7 ", numpy.mean(rzmeanval2), " ", numpy.std(rzmeanval2)
print " 7  8 ", numpy.mean(rzmeanval3), " ", numpy.std(rzmeanval3)
print " 8  9 ", numpy.mean(rzmeanval4), " ", numpy.std(rzmeanval4)
print " 9 10 ", numpy.mean(rzmeanval5), " ", numpy.std(rzmeanval5)
print " RPHI "
print " 5  6 ", numpy.mean(rphimeanval1), " ", numpy.std(rphimeanval1)
print " 6  7 ", numpy.mean(rphimeanval2), " ", numpy.std(rphimeanval2)
print " 7  8 ", numpy.mean(rphimeanval3), " ", numpy.std(rphimeanval3)
print " 8  9 ", numpy.mean(rphimeanval4), " ", numpy.std(rphimeanval4)
print " 9 10 ", numpy.mean(rphimeanval5), " ", numpy.std(rphimeanval5)

fp.close()
