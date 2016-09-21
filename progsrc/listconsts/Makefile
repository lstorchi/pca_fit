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

OBJ = listconsts.o

all: listconsts listconsts_int

listconsts: $(OBJ) 
	g++ -o $@ $(OBJ) -L../../lib -lpcafitter $(LIBS) 
	mv -f $@ ../../bin 

listconsts_int: listconsts_int.o
	g++ -o $@ $< -L../../lib -lpcafitter $(LIBS) 
	mv -f $@ ../../bin 

clean:
	rm -f *.o  ../../bin/listconsts ../../bin/listconsts_int
