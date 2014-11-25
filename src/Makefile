#!gmake
##################################################
#
# Makefile for gigafitter
#
# $Id$
#
###################################################

include ../config.mk

//CFLAGS += -DUSEALLVALUES 
//CFLAGS += -DDEBUG

ifeq ($(CERN),yes)
  CFLAGS += -I/afs/cern.ch/user/l/$(USER)/armadillo-3.930.4/include
  LIBS = -L/afs/cern.ch/user/l/$(USER)/armadillo-3.930.4 
endif

CFLAGS +=-I./ -I../include 

LIBS += -larmadillo -lm

OBJG = pcafitter.o \
       generatepca.o

OBJF = pcafitter.o \
       fitpca.o 

all: generatepca fitpca

generatepca: $(OBJG) 
	g++ -o $@ $(OBJG) $(LIBS) -larmadillo -lm
	mv -f $@ ../bin 

fitpca: $(OBJF) 
	g++ -o $@ $(OBJF) $(LIBS) -larmadillo -lm
	mv -f $@ ../bin

clean:
	rm -f *.o  ../bin/generatepca ../bin/fitpca
