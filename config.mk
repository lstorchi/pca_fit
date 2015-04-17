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

# if needed you sould specify here armadillo library path 
#CFLAGS += -I/afs/cern.ch/user/l/$(USER)/armadillo-3.930.4/include
#LIBS += -L/afs/cern.ch/user/l/$(USER)/armadillo-3.930.4 
CFLAGS =
LIBS =
ARMALIBPATH = /afs/cern.ch/user/l/lstorchi/armadillo-3.930.4
ARMAINCPATH = /afs/cern.ch/user/l/lstorchi/armadillo-3.930.4/include

CC = gcc
CXX = g++
LIBNAME = gffit

ifeq ($(DEBUG),no)
	CFLAGS = -Wall -W -O2
else
	CFLAGS = -Wall -W -O0 -g -DDEBUG
endif

CXXFLAGS = $(CFLAGS)
