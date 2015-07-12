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

  rval = numpy.zeros(numof)
  phival = numpy.zeros(numof)
  zval = numpy.zeros(numof)
  charge = 0.0

  for j in range(numof):
    coordline = fp.readline()

    coordline = p.sub(' ', coordline)
    coordline = coordline.lstrip()
    coordline = coordline.rstrip()
           
    coordlinelist = coordline.split(" ")

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


    rval[j] = ri
    phival[j] = phii
    zval[j] = zi

  paramline = fp.readline()

  paramline = p.sub(' ', paramline)
  paramline = paramline.lstrip()
  paramline = paramline.rstrip()

  paramlinelist = paramline.split(" ")

  pt = float(paramlinelist[0])
  phi = float(paramlinelist[1])
  eta = float(paramlinelist[3])
  z0 = float(paramlinelist[4])

  print "RPhi plane: "
  slope, intercept, r_value, p_value, std_err = stats.linregress(rval,phival)
  print slope, intercept, r_value
  print "c/pt: ", charge/pt, " phi: ", paramlinelist[1]

  print "RZ plane: "
  slope, intercept, r_value, p_value, std_err = stats.linregress(rval,zval)
  print slope, intercept, r_value
  print " eta: ", paramlinelist[3], "  z0: ", paramlinelist[4]

  print " "

fp.close()
