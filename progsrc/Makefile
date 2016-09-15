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
	$(MAKE) -C ./fitpca_split_single
	$(MAKE) -C ./listconsts
	$(MAKE) -C ./coeff_file_gen

clean:
	$(MAKE) clean -C ./generatepca_split
	$(MAKE) clean -C ./fitpca_split_single
	$(MAKE) clean -C ./listconsts
	$(MAKE) clean -C ./coeff_file_gen
