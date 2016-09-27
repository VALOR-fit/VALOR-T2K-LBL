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

TGT = generate_toymc_ensemble run_fit_marginalization  run_fc_marginalization run_mh_marginalization generate_toymc_ensemble_nuebarapp get_chisq_scale_nuebarapp

all: $(TGT)

# 
#
generate_toymc_ensemble_nuebarapp: FORCE
	$(CXX) $(CXXFLAGS) -c generate_toymc_ensemble_nuebarapp.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) generate_toymc_ensemble_nuebarapp.o $(LIBRARIES) $(VALOR_LIBRARIES) -o $(T2KOSCANA)/bin/generate_toymc_ensemble_nuebarapp

get_chisq_scale_nuebarapp: FORCE
	$(CXX) $(CXXFLAGS) -c get_chisq_scale_nuebarapp.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) get_chisq_scale_nuebarapp.o $(LIBRARIES) $(VALOR_LIBRARIES) -o $(T2KOSCANA)/bin/get_chisq_scale_nuebarapp

generate_toymc_ensemble: FORCE
	$(CXX) $(CXXFLAGS) -c generate_toymc_ensemble.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) generate_toymc_ensemble.o $(LIBRARIES) $(VALOR_LIBRARIES) -o $(T2KOSCANA)/bin/generate_toymc_ensemble

run_fit_marginalization: FORCE
	$(CXX) $(CXXFLAGS) -c run_fit_marginalization.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) run_fit_marginalization.o $(LIBRARIES) $(VALOR_LIBRARIES) -o $(T2KOSCANA)/bin/run_fit_marginalization

run_fc_marginalization: FORCE
	$(CXX) $(CXXFLAGS) -c run_fc_marginalization.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) run_fc_marginalization.o $(LIBRARIES) $(VALOR_LIBRARIES) -o $(T2KOSCANA)/bin/run_fc_marginalization	

run_mh_marginalization: FORCE
	$(CXX) $(CXXFLAGS) -c run_mh_marginalization.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) run_mh_marginalization.o $(LIBRARIES) $(VALOR_LIBRARIES) -o $(T2KOSCANA)/bin/run_mh_marginalization   

#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(T2KOSCANA)/bin/generate_toymc_ensemble_nuebarapp
	$(RM) $(T2KOSCANA)/bin/get_chisq_scale_nuebarapp
	$(RM) $(T2KOSCANA)/bin/generate_toymc_ensemble
	$(RM) $(T2KOSCANA)/bin/run_fit_marginalization
	$(RM) $(T2KOSCANA)/bin/run_fc_marginalization
	$(RM) $(T2KOSCANA)/bin/run_mh_marginalization

#distclean: FORCE
#	$(RM) $(VALOR_T2XK_INSTALLATION_PATH)/bin/generate_toymc_ensemble
#	$(RM) $(VALOR_T2XK_INSTALLATION_PATH)/bin/run_fit_marginalization
#       $(RM) $(VALOR_T2XK_INSTALLATION_PATH)/bin/run_fc_marginalization
#       $(RM) $(VALOR_T2XK_INSTALLATION_PATH)/bin/run_mh_marginalization 

FORCE:

# DO NOT DELETE
