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
	$(MAKE) -C ./progsrc/rootfilereader

clean:
	$(MAKE) clean -C src
	$(MAKE) clean -C ./progsrc/rootfilereader
