#
# Makefile for VALOR/T2XK test binary
#
# Steve Dennis <s.r.dennis \at liverpool.ac.uk>
# April 2016
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(T2KOSCANA)/make/Make.include

TGT = initialise_skexp_and_read_only

all: $(TGT)

# 
#
initialise_skexp_and_read_only: FORCE
	$(CXX) $(CXXFLAGS) -c initialise_skexp_and_read_only.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) initialise_skexp_and_read_only.o $(LIBRARIES) $(VALOR_LIBRARIES) -o $(T2KOSCANA)/bin/initialise_skexp_and_read_only

#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(T2KOSCANA)/bin/initialise_skexp_and_read_only

FORCE:

# DO NOT DELETE
