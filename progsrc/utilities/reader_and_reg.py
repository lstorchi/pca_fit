import sys
import re

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

numofline = (numofline-1) 

rint 

fp = open(filename, "r")

# jump first line
fp.readline()

for i in range(numofline):
  l = fp.readline()

  p = re.compile(r'\s+')
  line = p.sub(' ', l)
  line = line.lstrip()
  line = line.rstrip()
           
  plist = line.split(" ")

  numof = int(plist[1])

  for j in range(numof):
    coordline = fp.readline()

  paramline = fp.readline()

  print paramline




fp.close()
