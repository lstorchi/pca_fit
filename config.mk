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

# BITEWISE 
# yes to use integer, either generation and fit 
# no
USEINTBITEWISE = no

# if needed you sould specify here armadillo library path 
#CFLAGS += -I/afs/cern.ch/user/l/$(USER)/armadillo-3.930.4/include
#LIBS += -L/afs/cern.ch/user/l/$(USER)/armadillo-3.930.4 
#CFLAGS =
#LIBS =
ARMALIBPATH = /afs/cern.ch/user/l/lstorchi/armadillo-6.600.4
ARMAINCPATH = /afs/cern.ch/user/l/lstorchi/armadillo-6.600.4/include

CC = gcc
CXX = g++
LIBNAME = gffit

ifeq ($(DEBUG),no)
	CFLAGS = -Wall -W -O2
else
	CFLAGS = -Wall -W -O0 -g -DDEBUG
endif

ifeq ($(USEINTBITEWISE),yes)
  CFLAGS += -DINTBITEWISE
endif

CXXFLAGS = $(CFLAGS)
