import numpy
import sys

filename = ""

if (len(sys.argv) == 2 ):
  filename = sys.argv[1]

file = open(filename, "r")

values = []

for l in file:
  values.append(float(l))

y, x = numpy.histogram(values, bins=100)

for i in range(0, len(y)):
  print (x[i]+x[i+1])/2.0, " " , y[i]
