import numpy
import sys
import re

filename = ""


if (len(sys.argv) == 2 ):
  filename = sys.argv[1]

file = open(filename, "r")

w, h = 20, 16
mat = [[0 for x in range(w)] for y in range(w)]
xval = [0 for x in range(w)]

for itow in range(0,16):
  for j in range(0,20):
    l = file.readline()

    p = re.compile(r'\s+')
    line = p.sub(' ', l)
    line = line.lstrip()
    line = line.rstrip()
    linelist = line.split(",")
    #print linelist[1], linelist[2]
    
    eta = float(linelist[1])
    val = float(linelist[2])

    xval[j] = eta
    mat[j][itow] = val 


for i in range(0,20):
  print xval[i] , " ", numpy.mean(mat[i]), " " ,  numpy.std(mat[i])

