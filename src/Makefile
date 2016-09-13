#!gmake
##################################################
#
# Makefile for gigafitter
#
# $Id$
#
###################################################

include ../config.mk

CFLAGS += $(ARMAINCPATH)

CFLAGS +=-I./ -I../include 

CFLAGS += $(shell root-config --cflags) -I./
ROOTCI = rootcint
# custom root 6 
#CFLAGS += $(shell /usr/local/root-6/bin/root-config  --exec-prefix=/usr/local/root-6/ --cflags) -I./
#ROOTCI = usr/local/root-6/bin/rootcint

LIBNAME = libpcafitter.a

OBJ = pcafitter.o \
      utility.o \
      dicts.o \
      rootfilereader.o

all: $(LIBNAME)

dicts.cpp: linkdef.h
	$(ROOTCI) -v4 -f $@ -c $(CFLAGS) -I./ $<

dicts.o: dicts.cpp
	g++ -c $(CFLAGS) -o $@ $<

$(LIBNAME): $(OBJ)
	ar -rcs $@ $(OBJ)
	cp -f $(LIBNAME) ../lib

clean:
	rm -f *.o dicts.cpp dicts.h $(LIBNAME) ../lib/$(LIBNAME)
