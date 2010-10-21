#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(T2KANA)/src/make/Make.include

PACKAGE           = Yr1NuMuOscAnalysis
DICTIONARY        = _RootDict$(PACKAGE)
LIBNAME           = libT2KAna$(PACKAGE)

all : rootcint lib lib-link

# Include standard package makefile targets
#
include $(T2KANA)/src/make/Make.template

FORCE:

