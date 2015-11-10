#!gmake
##################################################
#
# Makefile for gigafitter
#
# $Id$
#
###################################################

include ../config.mk

ifeq ($(CERN),yes)
  CFLAGS += -I$(ARMAINCPATH)
endif

CFLAGS +=-I./ -I../include 

CFLAGS += $(shell root-config --cflags) -I./

LIBNAME = libpcafitter.a

OBJ = $(LIBNAME)(pcafitter.o) \
      $(LIBNAME)(utility.o) \
      $(LIBNAME)(dicts.o) \
      $(LIBNAME)(rootfilereader.o)

all: $(LIBNAME)

dicts.cpp: linkdef.h
	rootcint -v4 -f $@ -c $(CFLAGS) -I./ $<

dicts.o: dicts.cpp
	g++ -c $(CFLAGS) -o $@ $<

$(LIBNAME): $(OBJ)
	cp -f $(LIBNAME) ../lib

clean:
	rm -f *.o dicts.cpp dicts.h $(LIBNAME) ../lib/$(LIBNAME)
