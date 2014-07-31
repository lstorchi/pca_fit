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
  CFLAGS += -I/afs/cern.ch/user/l/lstorchi/armadillo-3.930.4/include
  LIBS = -L/afs/cern.ch/user/l/lstorchi/armadillo-3.930.4 
endif

LIBS += -larmadillo -lm

OBJG = generatepca.o

OBJF = fitpca.o

all: generatepca fitpca

generatepca: $(OBJG) 
	g++ -o $@ $(OBJG) $(LIBS) -larmadillo -lm

fitpca: $(OBJF) 
	g++ -o $@ $(OBJF) -lm

clean:
	rm -f *.o  generatepca fitpca
