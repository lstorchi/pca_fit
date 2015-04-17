#!gmake
##################################################
#
# Makefile for gigafitter
#
# $Id$
#
###################################################

include ../../config.mk

#CFLAGS += -DUSEALLVALUES 
#CFLAGS += -DDEBUG

ifeq ($(CERN),yes)
  CFLAGS += -I$(ARMAINCPATH)
  LIBS += -L$(ARMALIBPATH)
endif

CFLAGS +=-I./ -I../../include 

LIBS += -larmadillo -lm

OBJ = fitpca_split.o

all: fitpca_split 

fitpca_split: $(OBJ) 
	g++ -o $@ $(OBJ) -L../../lib -lpcafitter $(LIBS) 
	mv -f $@ ../../bin 

clean:
	rm -f *.o  ../../bin/fitpca_split
