#
# Makefile for VALOR/T2XK apps
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# Sep 10, 2015
#
# Tom Dealtry <t.dealtry \at lancaster.ac.uk>
# Sep 24, 2015
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(T2KOSCANA)/make/Make.include

TGT = generate_toymc_ensemble run_fit_marginalization

all: $(TGT)

# 
#
generate_toymc_ensemble: FORCE
	$(CXX) $(CXXFLAGS) -c generate_toymc_ensemble.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) generate_toymc_ensemble.o $(LIBRARIES) $(VALOR_LIBRARIES) -o $(T2KOSCANA)/bin/generate_toymc_ensemble

run_fit_marginalization: FORCE
	$(CXX) $(CXXFLAGS) -c run_fit_marginalization.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) run_fit_marginalization.o $(LIBRARIES) $(VALOR_LIBRARIES) -o $(T2KOSCANA)/bin/run_fit_marginalization


#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(T2KOSCANA)/bin/generate_toymc_ensemble
	$(RM) $(T2KOSCANA)/bin/run_fit_marginalization

#distclean: FORCE
#	$(RM) $(VALOR_T2XK_INSTALLATION_PATH)/bin/generate_toymc_ensemble
#	$(RM) $(VALOR_T2XK_INSTALLATION_PATH)/bin/run_fit_marginalization

FORCE:

# DO NOT DELETE
