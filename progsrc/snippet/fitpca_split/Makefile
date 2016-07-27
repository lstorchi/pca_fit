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

OBJ = fitpca_split.o

all: fitpca_split 

fitpca_split: $(OBJ) 
	g++ -o $@ $(OBJ) -L../../lib -lpcafitter $(LIBS) 
	mv -f $@ ../../bin 

clean:
	rm -f *.o  ../../bin/fitpca_split 
