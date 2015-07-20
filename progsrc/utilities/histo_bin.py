import numpy
import sys
import re

filename = ""

bintouse = 100

if (len(sys.argv) == 2 ):
  filename = sys.argv[1]
elif (len(sys.argv) == 3 ):
  filename = sys.argv[1]
  bintouse = int(sys.argv[2])

file = open(filename, "r")

values = []

for l in file:
  p = re.compile(r'\s+')
  line = p.sub(' ', l)
  line = line.lstrip()
  line = line.rstrip()
  linelist = line.split(" ")

  firstval = float(linelist[0])

  #if (firstval >= -0.2 and firstval <= 0.0):
  values.append(float(linelist[1]))

y, x = numpy.histogram(values, bins=bintouse)

for i in range(0, len(y)):
  print (x[i]+x[i+1])/2.0, " " , y[i]
