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
	$(MAKE) -C ./progsrc/generatepca_split
	$(MAKE) -C ./progsrc/fitpca_split
	$(MAKE) -C ./progsrc/readrootfile

clean:
	$(MAKE) clean -C src
	$(MAKE) clean -C ./progsrc/generatepca_split
	$(MAKE) clean -C ./progsrc/fitpca_split
	$(MAKE) clean -C ./progsrc/readrootfile
