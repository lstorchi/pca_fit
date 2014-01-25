#!gmake
##################################################
#
# Makefile for gffit
#
# $Id$
#
###################################################

libgffit:
	$(MAKE) -C src

all: libgffit

clean:
	$(MAKE) clean -C src
