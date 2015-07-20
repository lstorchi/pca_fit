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
  line = p.sub(' ', l)
  line = coordline.lstrip()
  line = coordline.rstrip()
  linelist = coordline.split(" ")

  firstval = float(linelist[0])

  if (firstval >= -0.4 and firstval <= -0.2):
    values.append(float(linelist[1]))

y, x = numpy.histogram(values, bins=bintouse)

for i in range(0, len(y)):
  print (x[i]+x[i+1])/2.0, " " , y[i]
