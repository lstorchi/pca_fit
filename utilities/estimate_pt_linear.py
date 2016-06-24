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

ptdiffvalues = []

pt3_10diffvalues = []
pt10_50diffvalues = []
pt50_100diffvalues = []
pt100_200diffvalues = []

x = []
y = []

q = 0.0057544097038
b = 1.64483470352e-06

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

  zval = numpy.zeros(numof)
  rval = numpy.zeros(numof)
  zval13 = numpy.zeros(numof)
  rval13 = numpy.zeros(numof)
  phival = numpy.zeros(numof)
  charge = 0.0

  layersids = ""

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
    phii = math.atan2(yi, xi)

    rval[j] = ri
    phival[j] = phii
    zval[j] = zi

    if (j < 3):
      rval13[j] = ri
      zval13[j] = zi

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
  if (layersids != "5678910"):
    print >> sys.stderr, "Wrong seq: ", layersids
  else:
    if pt >= 3.0:
      slope = (phival[5]-phival[0])/(rval[5]-rval[0])
      est = (slope*q + b)
      print "c/pt: ", charge/pt, "     estimated: ", est
      ptdiffvalues.append((charge/pt) - est)

      if pt <= 10:
        pt3_10diffvalues.append((charge/pt) - est)
      elif pt > 10 and pt <= 50:
        pt10_50diffvalues.append((charge/pt) - est)
      elif pt > 50 and pt <= 100:
        pt50_100diffvalues.append((charge/pt) - est)
      elif pt > 100 and pt <= 200:
        pt100_200diffvalues.append((charge/pt) - est)

print "Num of events: ", len(ptdiffvalues)
print "Mean val: ", numpy.mean(ptdiffvalues)
print "STD  val: ", numpy.std(ptdiffvalues)
print "3 a 10"
print "Mean val: ", numpy.mean(pt3_10diffvalues)
print "STD  val: ", numpy.std(pt3_10diffvalues)
print "10 a 50"
print "Mean val: ", numpy.mean(pt10_50diffvalues)
print "STD  val: ", numpy.std(pt10_50diffvalues)
print "50 a 100"
print "Mean val: ", numpy.mean(pt50_100diffvalues)
print "STD  val: ", numpy.std(pt50_100diffvalues)
print "100 a 200"
print "Mean val: ", numpy.mean(pt100_200diffvalues)
print "STD  val: ", numpy.std(pt100_200diffvalues)


fp.close()
