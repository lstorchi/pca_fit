#!gmake
##################################################
#
# Makefile for filter
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

OBJ = filter.o \
      jobparams.o

all: filter

filter: $(OBJ) 
	g++ -o $@ $(OBJ) -L../../lib -lpcafitter $(LIBS) 
	mv -f $@ ../../bin 

clean:
	rm -f *.o  ../../bin/filter
