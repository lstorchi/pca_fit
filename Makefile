#!gmake
##################################################
#
# Makefile for gffit
#
# $Id$
#
###################################################

all: 
	$(MAKE) -C src
	$(MAKE) -C rootfilereader

clean:
	$(MAKE) clean -C src
	$(MAKE) clean -C rootfilereader
