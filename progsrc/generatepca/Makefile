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
  CFLAGS += -I/afs/cern.ch/user/l/$(USER)/armadillo-3.930.4/include
  LIBS += -L/afs/cern.ch/user/l/$(USER)/armadillo-3.930.4 
endif

CFLAGS +=-I./ -I../../include 

LIBS += -larmadillo -lm

OBJ = generatepca.o

all: generatepca 

generatepca: $(OBJ) 
	g++ -o $@ $(OBJ) -L../../lib -lpcafitter $(LIBS) 
	mv -f $@ ../../bin 

clean:
	rm -f *.o  ../../bin/generatepca 
