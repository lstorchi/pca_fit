import numpy
import sys

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
  values.append(float(l))

y, x = numpy.histogram(values, bins=bintouse)

for i in range(0, len(y)):
  print (x[i]+x[i+1])/2.0, " " , y[i]
