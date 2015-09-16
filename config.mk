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
CERN = no

# BITEWISE 
# yes to use integer 
# no
USEINTBITEWISE = yes

# if needed you sould specify here armadillo library path 
CFLAGS += -I/home/atanu/Downloads/armadillo-5.200.2/include
LIBS += -L/home/atanu/Downloads/armadillo-5.200.2 
#CFLAGS += -I/afs/cern.ch/user/l/$(USER)/armadillo-3.930.4/include
#LIBS += -L/afs/cern.ch/user/l/$(USER)/armadillo-3.930.4 
#CFLAGS =
#LIBS =
ARMALIBPATH = /afs/cern.ch/user/l/lstorchi/armadillo-3.930.4
ARMAINCPATH = /afs/cern.ch/user/l/lstorchi/armadillo-3.930.4/include
#ARMALIBPATH = /home/atanu/Downloads/armadillo-5.200.2
#ARMAINCPATH = /home/atanu/Downloads/armadillo-5.200.2/include

CC = gcc
CXX = g++
LIBNAME = gffit

ifeq ($(DEBUG),no)
	CFLAGS = -Wall -W -O2
else
	CFLAGS = -Wall -W -O0 -g -DDEBUG
endif

ifeq ($(USEINTBITEWISE),yes)
  CFLAGS += -DINTBITEWISE -std=c++11
endif

CXXFLAGS = $(CFLAGS)
