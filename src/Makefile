#!gmake
##################################################
#
# Makefile for libgffit
#
# $Id$
#
###################################################

include ../config.mk

#CFLAGS += -I../inc -D DEBUG_SVT
CFLAGS += -I../inc 

LIB = lib$(LIBNAME).a

OBJ = $(LIB)(svtsim_functions.o) \
      $(LIB)(utility.o) \
      $(LIB)(gf_functions.o) \
      $(LIB)(blob_functions.o)

all: $(LIB)

$(LIB): $(OBJ) 
	cp -f $@ ../lib/

clean:
	rm -f *.o $(LIB) ../lib/$(LIB) 
