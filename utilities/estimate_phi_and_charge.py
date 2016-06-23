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

phidiffvalues = []
x = []
y = []

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
    print "RZPhi plane using layers 1 and 6: "
    slope = (phival[5]-phival[0])/(rval[5]-rval[0])
    print "layers 1 6 c/pt: ", charge/pt, "     slope: ", slope
    x.append(charge/pt)
    y.append(slope)
    intercept = phival[0] - slope*rval[0]
    print "layers 1 6  phi: ", phi,  " intercept: ", intercept, " diff: ", phi-intercept

    phidiffvalues.append(phi-intercept)

print "phi layers 1 6: " 
print "Num of events: ", len(phidiffvalues)
print "Mean val: ", numpy.mean(phidiffvalues)
print "STD  val: ", numpy.std(phidiffvalues)

slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

print "Lin Regr slope, intercept, r_value, p_value, std_err"
print slope, intercept, r_value, p_value, std_err

fp.close()
