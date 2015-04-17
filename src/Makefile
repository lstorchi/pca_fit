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

LIBNAME = libpcafitter.a

OBJ = $(LIBNAME)(pcafitter.o) \
      $(LIBNAME)(utility.o) 

all: $(LIBNAME)

$(LIBNAME): $(OBJ)
	cp -f $(LIBNAME) ../lib

clean:
	rm -f *.o  $(LIBNAME) ../lib/$(LIBNAME)
