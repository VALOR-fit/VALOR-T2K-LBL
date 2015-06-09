#! /usr/bin/env python

import sys
import optparse
import os
import math # cos, log, pi etc
import string
import logging
import glob

class T2KAnaError(Exception): pass

optionParser = optparse.OptionParser(usage ="""\
usage: %prog [options]
    Run generate_toymc_ensemple_joint_bar.C many times""")

optionParser.add_option("-a","--all", action="store_true", dest="do_all", default=False, help="If set, will run indiviudal systematics, as well as the categories")
(opt, args) = optionParser.parse_args()

#tag = "toymc_nuebarapp_NuOscOnly_Toy_23_th13dcp_Corr_MH_betabar_0_4p0108E20POT_DEV_10k_C"
tag = "toymc_nuebarapp_NuOscOnly_Toy_23_th13dcp_Corr_MH_betabar_1_8E20POT_DEV_Toy_50k"
#tag = "Prova_rootfile_8E20POT_DEV_10_A"

run = 0

first_experiment = 0                
last_experiment = 50100               
#nexperiments = 10 #10100

hypothesis = "3f_osc_normbar"
printsystvariation = 0

input_mixing12_nh = -99    
input_mixing13_nh = -999   
input_mixing23_nh = -999  
input_deltacp_nh  = -999   
input_dm12sq_nh   = -99   
input_dm23sq_nh   = -999 
input_mixing12_ih = -99   
input_mixing13_ih = -999    
input_mixing23_ih = -999  
input_deltacp_ih  = -999  
input_dm12sq_ih   = -99   
input_dm23sq_ih   = -999  
input_beta_bar    = 1

#config_file = "analysis_configuration_run1to6_test_app_4.0108E20POT_fluxJune.xml"
#config_file = "analysis_configuration_run1to6_prefit_test_app_4.0108E20POT_fluxJune.xml"
#config_file = "analysis_configuration_run1to6_test_app_9E20POT_fluxJune.xml"
config_file = "analysis_configuration_run1to6_test_app_8E20POT_fluxJune.xml"
#config_file = "analysis_configuration_run1to6_test_app_6.5E20POT_fluxJune.xml"

rethrow_if_parameters_unphysical = 1
stat_fluctuations = 1
tweak_systematics = 1

doMargMH = 1
PriorProbNH = 0.684

save_in_xml_format = 0
save_in_root_format = 1

points_to_add_1Rmu = 0 # char
points_to_add_1Re = 0 # char
points_to_add_1Rmu_RHC = 0 # char
points_to_add_1Re_RHC = 0 # char
points_to_add_MultiRmu = 0 # char
points_to_add_MultiRe = 0 # char
xml_file_to_add = 0 # char

save_fluct_in_xml = 0 

dmsqscheme = "auto"
hierarchy  = "normal"
mixparamchoice = "sinsqtheta"
mixparamoctant = "two"

file_priors = "../../../inputs/priors/joint_T2Kreact_run14/Priors_T2Kreact_run14_Mach3_22Maj15.root"

command = 'root2k -l -b -q \'$T2KOSCANA/macros/toymcgen/generate_toymc_ensemble_joint_bar.C("'+tag+'",'+str(run)+','+str(first_experiment)+','+str(last_experiment)+',"'+hypothesis+'",'+str(printsystvariation)+','+str(input_mixing12_nh)+','+str(input_mixing13_nh)+','+str(input_mixing23_nh)+','+str(input_deltacp_nh)+','+str(input_dm12sq_nh)+','+str(input_dm23sq_nh)+','+str(input_mixing12_ih)+','+str(input_mixing13_ih)+','+str(input_mixing23_ih)+','+str(input_deltacp_ih)+','+str(input_dm12sq_ih)+','+str(input_dm23sq_ih)+','+str(input_beta_bar)+',"'+config_file+'",'+str(rethrow_if_parameters_unphysical)+','+str(stat_fluctuations)+','+str(tweak_systematics)+','+str(doMargMH)+','+str(PriorProbNH)+','+str(save_in_xml_format)+','+str(save_in_root_format)+','+str(points_to_add_1Rmu)+','+str(points_to_add_1Re)+','+str(points_to_add_1Rmu_RHC)+','+str(points_to_add_1Re_RHC)+','+str(points_to_add_MultiRmu)+','+str(points_to_add_MultiRe)+','+str(xml_file_to_add)+','+str(save_fluct_in_xml)+',"'+dmsqscheme+'","'+hierarchy+'","'+mixparamchoice+'","'+mixparamoctant+'","'+file_priors+'")\''

print command
if os.system(command) == 1:
  sys.exit(1)
