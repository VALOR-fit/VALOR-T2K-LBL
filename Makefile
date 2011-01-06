#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include Make.include

PACKAGE           = T2KDisapp2010a
DICTIONARY        = _RootDict$(PACKAGE)
LIBNAME           = lib$(PACKAGE)

all : rootcint lib lib-link

# Include standard package makefile targets
#
include Make.template

FORCE:

