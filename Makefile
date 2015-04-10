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
	$(MAKE) -C ./progsrc/generatepca
	$(MAKE) -C ./progsrc/fitpca
	$(MAKE) -C ./progsrc/rootfilereader

clean:
	$(MAKE) clean -C src
	$(MAKE) clean -C ./progsrc/generatepca
	$(MAKE) clean -C ./progsrc/fitpca
	$(MAKE) clean -C ./progsrc/rootfilereader
