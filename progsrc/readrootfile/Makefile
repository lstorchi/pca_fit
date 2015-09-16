#!gmake
##################################################
#
# Makefile for gigafitter
#
# $Id$
#
###################################################

include ../../config.mk

CFLAGS = $(shell root-config --cflags) -I./
LIBS = $(shell root-config --libs)
# to use xrootd suggested by I. Tomalin
#LIBS = $(shell root-config --libs) –lNetxNG2

ifeq ($(USEINTBITEWISE),yes)
  CFLAGS += -DINTBITEWISE 
endif

OBJ = readrootfile.o \
      dicts.o

all: readrootfile

readrootfile.o: readrootfile.cxx
	g++ -c $(CFLAGS) -o $@ $< 

dicts.cpp: linkdef.h
	rootcint -v4 -f $@ -c $(CFLAGS) -I./ $<

dicts.o: dicts.cpp
	g++ -c $(CFLAGS) -o $@ $<

readrootfile: $(OBJ) 
	g++ -o $@ $(OBJ) $(LIBS)
	mv -f $@ ../../bin

clean:
	rm -f *.o ../../bin/readrootfile dicts.cpp dicts.h
