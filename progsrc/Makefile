#!gmake
##################################################
#
# Makefile for pca_fit
#
# $Id$
#
###################################################

all: 
	$(MAKE) -C ./generatepca_split
	$(MAKE) -C ./fitpca_split
	$(MAKE) -C ./listconsts
	$(MAKE) -C ./readrootfile

clean:
	$(MAKE) clean -C ./generatepca_split
	$(MAKE) clean -C ./fitpca_split
	$(MAKE) clean -C ./listconsts
	$(MAKE) clean -C ./readrootfile
