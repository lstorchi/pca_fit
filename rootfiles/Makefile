#!gmake
##################################################
#
# Makefile for gigafitter
#
# $Id$
#
###################################################

CFLAGS = $(shell root-config --cflags)
LIBS = $(shell root-config --libs)

OBJ = readrootfile.o

all: readrootfile

readrootfile.o: readrootfile.cxx
	g++ -c $(CFLAGS) -o $@ $< 

readrootfile: $(OBJ) 
	g++ -o $@ $(OBJ) $(LIBS)

clean:
	rm -f *.o readrootfile
