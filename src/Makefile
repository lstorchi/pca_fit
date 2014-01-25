#!gmake
##################################################
#
# Makefile for libgffit
#
# $Id$
#
###################################################

include ../config.mk

CFLAGS += -I../inc

LIB = lib$(LIBNAME).a

OBJ = $(LIB)(svtsim_functions.o)

all: $(LIB)

$(LIB): $(OBJ) 
	cp -f $@ ../lib/

clean:
	rm -f *.o $(LIB) ../lib/$(LIB) 
