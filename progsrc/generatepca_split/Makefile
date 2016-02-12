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
CFLAGS += $(shell root-config --cflags) 

LIBS += -larmadillo -lm
LIBS += $(shell root-config --libs)

OBJ = generatepca_split.o

all: generatepca_split 

generatepca_split: $(OBJ) 
	g++ -o $@ $(OBJ) -L../../lib -lpcafitter $(LIBS) 
	mv -f $@ ../../bin 

clean:
	rm -f *.o  ../../bin/generatepca_split 

generatepca_split.o: ../../include/pcaconst.hpp ../../include/pcaffunctype.hpp
