#!gmake
##################################################
#
# Makefile for gigafitter
#
# $Id$
#
###################################################

include ../config.mk

CFLAGS += -O2 -W -Wall

OBJG = generatepca.o

OBJF = fitpca.o

all: generatepca fitpca

generatepca: $(OBJG) 
	g++ -o $@ $(OBJG) -lm

fitpca: $(OBJF) 
	g++ -o $@ $(OBJF) -lm

clean:
	rm -f *.o  generatepca fitpca
