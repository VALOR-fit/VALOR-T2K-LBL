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

TGT = make_table_systematic_effects_t2k \
      bievent_build \
      bievent_plot

#write_fd_prediction_pdf

all: $(TGT)
# Throw a user-specified parameter set N times and output a table containing the mean, RMS and RMS/mean for all samples, with T2K-specific calculations.
make_table_systematic_effects_t2k: FORCE
	$(CXX) $(CXXFLAGS) -c make_table_systematic_effects_t2k.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) make_table_systematic_effects_t2k.o $(LIBRARIES) -o $(VALOR_EXPT)/bin/valor_t2k-make_table_systematic_effects_t2k

bievent_build: FORCE
	$(CXX) $(CXXFLAGS) -c bievent_build.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) bievent_build.o $(LIBRARIES) -o $(VALOR_EXPT)/bin/valor-bievent_build

bievent_plot: FORCE
	$(CXX) $(CXXFLAGS) -c bievent_plot.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) bievent_plot.o $(LIBRARIES) -o $(VALOR_EXPT)/bin/valor-bievent_plot


#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(VALOR_EXPT)/bin/valor_t2k-make_table_systematic_effects_t2k
	$(RM) $(VALOR_EXPT)/bin/valor-bievent_build
	$(RM) $(VALOR_EXPT)/bin/valor-bievent_plot

distclean: FORCE
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_t2k-make_table_systematic_effects_t2k
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor-bievent_build
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor-bievent_plot

FORCE:
	
# DO NOT DELETE
