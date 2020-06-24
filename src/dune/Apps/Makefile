#
# Makefile for VALOR/DUNE apps
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(VALOR_EXPT)/make/Make.include

DUNE_LIBRARIES = 

TGT = make_caf_minidst \
      make_nd_caf_minidst 

all: $(TGT)

# Generate CAFs based DSTs
#
make_caf_minidst: FORCE
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c make_caf_minidst.cxx $(GINCLUDES)
	$(LD) $(LDFLAGS) make_caf_minidst.o $(GLIBRARIES) $(VALOR_GLIBRARIES) -o $(VALOR_EXPT)/bin/valor_dune-make_caf_minidst
else
	@echo "Skipping make_caf_minidst due to no GENIE"
endif

make_nd_caf_minidst: FORCE
ifdef GENIE
	$(CXX) $(CXXFLAGS) -c make_nd_caf_minidst.cxx $(GINCLUDES)
	$(LD) $(LDFLAGS) make_nd_caf_minidst.o $(GLIBRARIES) $(VALOR_GLIBRARIES) -o $(VALOR_EXPT)/bin/valor_dune-make_nd_caf_minidst
else
	@echo "Skipping make_nd_caf_minidst due to no GENIE"
endif

# Write out fd prediction pdf
#
write_fd_prediction_pdf: FORCE
	$(CXX) $(CXXFLAGS) -c write_fd_prediction_pdf.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) write_fd_prediction_pdf.o $(LIBRARIES) -o $(VALOR_EXPT)/bin/valor_dune-write_fd_prediction_pdf
	
#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(VALOR_EXPT)/bin/valor_dune-make_caf_minidst
	$(RM) $(VALOR_EXPT)/bin/valor_dune-make_nd_caf_minidst
	$(RM) $(VALOR_EXPT)/bin/valor_dune-write_fd_prediction_pdf

distclean: FORCE
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_dune-make_caf_minidst
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/valor_dune-make_nd_caf_minidst
	$(RM) $(VALOR_EXPT_INSTALLATION_PATH)/bin/write_fd_prediction_pdf

FORCE:
	
# DO NOT DELETE
