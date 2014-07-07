##################################################
#
# config.mk 
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

CC = gcc
CXX = g++
LIBNAME = gffit

ifeq ($(DEBUG),no)
	CFLAGS = -Wall -W -O2
else
	CFLAGS = -Wall -W -O0 -g -DDEBUG
endif

CXXFLAGS = $(CFLAGS)
