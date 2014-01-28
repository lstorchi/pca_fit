#!gmake
##################################################
#
# Makefile for gffit
#
# $Id$
#
###################################################

all: libgffit 
	$(MAKE) -C gigafitter

libgffit:
	$(MAKE) -C src

clean:
	$(MAKE) clean -C src
	$(MAKE) clean -C gigafitter
