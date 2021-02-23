#
# Makefile for VALOR/T2K apps
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(VALOR_EXPT)/make/Make.include

T2K_LIBRARIES =

TGT = make_table_systematic_effects_t2k

all: $(TGT)

# Throw a user-specified parameter set N times and output a table containing the mean, RMS and RMS/mean for all samples, with T2K-specific calculations.
#
make_table_systematic_effects_t2k: FORCE
	$(CXX) $(CXXFLAGS) -c make_table_systematic_effects_t2k.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) make_table_systematic_effects_t2k.o $(LIBRARIES) -o $(VALOR_EXPT)/bin/valor_t2k-make_table_systematic_effects_t2k
	
#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(VALOR_EXPT)/bin/valor_t2k-make_table_systematic_effects_t2k

distclean: FORCE
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_t2k-make_table_systematic_effects_t2k

FORCE:
	
# DO NOT DELETE
