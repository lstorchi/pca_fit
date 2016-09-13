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
#CFLAGS += $(shell /usr/local/root-6/bin/root-config  --exec-prefix=/usr/local/root-6/ --cflags) 

LIBS += -larmadillo -lm

LIBS += $(shell root-config --libs)
#LIBS += $(shell /usr/local/root-6/bin/root-config  --exec-prefix=/usr/local/root-6/ --libs)

OBJ = generatepca_split.o

all: generatepca_split 

generatepca_split: $(OBJ) 
	g++ -o $@ $(OBJ) -L../../lib -lpcafitter $(LIBS) 
	mv -f $@ ../../bin 

clean:
	rm -f *.o  ../../bin/generatepca_split 

generatepca_split.o: ../../include/pcaconst.hpp ../../include/pcaffunctype.hpp
