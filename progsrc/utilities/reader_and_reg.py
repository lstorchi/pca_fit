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

etadiffvalues = []
#etadiffvalues_pzbin = []
etadiffvalues_etabin1 = []
etadiffvalues_etabin2 = []
etadiffvalues_etabin3 = []
etadiffvalues_etabin4 = []
etadiffvalues_etabin5 = []

z0diffvalues = []
z0diffvalues_etabin1 = []
z0diffvalues_etabin2 = []
z0diffvalues_etabin3 = []
z0diffvalues_etabin4 = []
z0diffvalues_etabin5 = []


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
    phii = math.acos(xi/ri)

    rval[j] = ri
    phival[j] = phii
    zval[j] = zi

    #print ri, " ", zi

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
    #print " "
    #print "RPhi plane: "
    #slope, intercept, r_value, p_value, std_err = stats.linregress(rval,phival)
    #print "r: ", r_value
    #print "c/pt: ", charge/pt, " slope: ", slope
    #print "phi: ", phi , " intercept: ", intercept, " diff: ", intercept-phi
  
    #print "RZ plane: "
    #slope, intercept, r_value, p_value, std_err = stats.linregress(rval,zval)
    #print "r: ", r_value
    #print "eta: ", eta, " slope: ", slope, " diff: ", eta-slope
    #print "z0: ", z0, " intercept: ", intercept, " diff: ", z0-intercept
  
    print "RZ plane using layers 1 and 3: "
    slope = (zval[2]-zval[0])/(rval[2]-rval[0])
    print "layers 1 3 eta: ", eta, "     slope: ", slope,     " diff: ", eta-slope, \
        " theta: ", theta*180.0/math.pi, " pz: ", pz
    intercept = zval[0] - slope*rval[0]
    print "layers 1 3  z0: ", z0,  " intercept: ", intercept, " diff: ", z0-intercept, \
        " eta: " , eta

    #print "RZPhi plane using layers 1 and 6: "
    #slope = (phival[5]-phival[0])/(rval[5]-rval[0])
    #print "layers 1 6 c/pt: ", charge/pt, "     slope: ", slope
    #intercept = phival[0] - slope*rval[0]
    #print "layers 1 6  phi: ", phi,  " intercept: ", intercept, " diff: ", phi-intercept
    #print " "

    #if (pz >= 0.0) and (pz <= 1.0):
    #  etadiffvalues_pzbin.append(eta-slope)

    etadiffvalues.append(eta-slope)
    z0diffvalues.append(z0-intercept)

    if (eta >= -0.6) and (eta < -0.4):
      etadiffvalues_etabin1.append(eta-slope)
      z0diffvalues_etabin1.append(z0-intercept)
    elif (eta >= -0.4) and (eta < -0.2):
      etadiffvalues_etabin2.append(eta-slope)
      z0diffvalues_etabin2.append(z0-intercept)
    elif (eta >= -0.2) and (eta < 0.0):
      etadiffvalues_etabin3.append(eta-slope)
      z0diffvalues_etabin3.append(z0-intercept)
    elif (eta >= 0.0) and (eta < 0.2):
      etadiffvalues_etabin4.append(eta-slope)
      z0diffvalues_etabin4.append(z0-intercept)
    elif (eta >= 0.2) and (eta < 0.4):
      etadiffvalues_etabin5.append(eta-slope)
      z0diffvalues_etabin4.append(z0-intercept)


print "Eta"
print "Num of events: ", len(etadiffvalues)
print "Mean val: ", numpy.mean(etadiffvalues)
print "STD  val: ", numpy.std(etadiffvalues)
#print " "
#print "pz bin 0.0 1.0"
#print "Num of events: ", len(etadiffvalues_pzbin)
#print "Mean val: ", numpy.mean(etadiffvalues_pzbin)
#print "STD  val: ", numpy.std(etadiffvalues_pzbin)
#print " "
print "eta bin -0.6 -0.4"
print "Num of events: ", len(etadiffvalues_etabin1)
print "Mean val: ", numpy.mean(etadiffvalues_etabin1)
print "STD  val: ", numpy.std(etadiffvalues_etabin1)
print "eta bin -0.4 -0.2"
print "Num of events: ", len(etadiffvalues_etabin2)
print "Mean val: ", numpy.mean(etadiffvalues_etabin2)
print "STD  val: ", numpy.std(etadiffvalues_etabin2)
print "eta bin -0.2 0.0"
print "Num of events: ", len(etadiffvalues_etabin3)
print "Mean val: ", numpy.mean(etadiffvalues_etabin3)
print "STD  val: ", numpy.std(etadiffvalues_etabin3)
print "eta bin 0.0 0.2"
print "Num of events: ", len(etadiffvalues_etabin4)
print "Mean val: ", numpy.mean(etadiffvalues_etabin4)
print "STD  val: ", numpy.std(etadiffvalues_etabin4)
print "eta bin 0.2 0.4"
print "Num of events: ", len(etadiffvalues_etabin5)
print "Mean val: ", numpy.mean(etadiffvalues_etabin5)
print "STD  val: ", numpy.std(etadiffvalues_etabin5)

print "z0"
print "Num of events: ", len(z0diffvalues)
print "Mean val: ", numpy.mean(z0diffvalues)
print "STD  val: ", numpy.std(z0diffvalues)
print "eta bin -0.6 -0.4"
print "Num of events: ", len(z0diffvalues_etabin1)
print "Mean val: ", numpy.mean(z0diffvalues_etabin1)
print "STD  val: ", numpy.std(z0diffvalues_etabin1)
print "eta bin -0.4 -0.2"
print "Num of events: ", len(z0diffvalues_etabin2)
print "Mean val: ", numpy.mean(z0diffvalues_etabin2)
print "STD  val: ", numpy.std(z0diffvalues_etabin2)
print "eta bin -0.2 0.0"
print "Num of events: ", len(z0diffvalues_etabin3)
print "Mean val: ", numpy.mean(z0diffvalues_etabin3)
print "STD  val: ", numpy.std(z0diffvalues_etabin3)
print "eta bin 0.0 0.2"
print "Num of events: ", len(z0diffvalues_etabin4)
print "Mean val: ", numpy.mean(z0diffvalues_etabin4)
print "STD  val: ", numpy.std(z0diffvalues_etabin4)
print "eta bin 0.2 0.4"
print "Num of events: ", len(z0diffvalues_etabin5)
print "Mean val: ", numpy.mean(z0diffvalues_etabin5)
print "STD  val: ", numpy.std(z0diffvalues_etabin5)

fp.close()
