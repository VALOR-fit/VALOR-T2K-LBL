#
# Makefile for VALOR/T2K apps
#
# NOTE 2020-06-25: Currently nothing to make, but keeping the instructions of an old DUNE app as an example
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

TGT = bievent_build \
      bievent_plot
#write_fd_prediction_pdf

all: $(TGT)

# Write out fd prediction pdf
#
#write_fd_prediction_pdf: FORCE
#	$(CXX) $(CXXFLAGS) -c write_fd_prediction_pdf.cxx $(INCLUDES)
#	$(LD) $(LDFLAGS) write_fd_prediction_pdf.o $(LIBRARIES) -o $(VALOR_EXPT)/bin/valor_dune-write_fd_prediction_pdf
	
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
	$(RM) $(VALOR_EXPT)/bin/valor-bievent_build
	$(RM) $(VALOR_EXPT)/bin/valor-bievent_plot
	#$(RM) $(VALOR_EXPT)/bin/valor_dune-write_fd_prediction_pdf

distclean: FORCE
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor-bievent_build
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor-bievent_plot
	#$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/write_fd_prediction_pdf

FORCE:
	
# DO NOT DELETE
