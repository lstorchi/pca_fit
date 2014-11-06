##################################################
#
# config.mk lstorchi to be used for all the basic 
#           configuration related to all sources
#
# $Id$
#
###################################################

# Try to identify machine architecture
PLATFORM = $(shell uname -s)

# DEBUG: 
# yes, 
# no
DEBUG=no

# CERN 
# yes if in lxplus
# no
CERN = yes

CC = gcc
CXX = g++
LIBNAME = gffit

ifeq ($(DEBUG),no)
	CFLAGS = -Wall -W -O2
else
	CFLAGS = -Wall -W -O0 -g -DDEBUG
endif

CXXFLAGS = $(CFLAGS)
