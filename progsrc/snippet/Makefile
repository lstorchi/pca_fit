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

CFLAGS += $(ARMAINCPATH)
LIBS += $(ARMALIBPATH)

CFLAGS +=-I./ -I../../include 
CFLAGS += $(shell root-config --cflags)

LIBS += -larmadillo -lm
LIBS += $(shell root-config --libs)

OBJ1 = fitpca_split_bin.o

OBJ2 = fitpca_split_singleconstset.o

fitpca_split_bin: $(OBJ1) 
	g++ -o $@ $(OBJ1) -L../../lib -lpcafitter $(LIBS) 
	mv -f $@ ../../bin 

fitpca_split_singleconstset: $(OBJ2) 
	g++ -o $@ $(OBJ2) -L../../lib -lpcafitter $(LIBS) 
	mv -f $@ ../../bin 

clean:
	rm -f *.o  ../../bin/fitpca_split ../../bin/fitpca_split_bin
