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
  CFLAGS += -I/afs/cern.ch/user/l/$(USER)/armadillo-3.930.4/include
endif

CFLAGS +=-I./ -I../include 

LIBNAME = libpcafitter.a

OBJ = $(LIBNAME)(pcafitter.o) \
      $(LIBNAME)(utility.o) 

all: $(LIBNAME)

$(LIBNAME): $(OBJ)
	mv -f $(LIBNAME) ../lib

clean:
	rm -f *.o  ../lib/$(LIBNAME)
