#! /usr/bin/env python

import sys
import optparse
import os
import math # cos, log, pi etc
import string
import logging

"""
This script will run call submit_toyexpt_fitting_job_fitsyst and pass it the relevant options, for a given setting

Run this by /path/to/common_runs.py (you may need to chmod +x)
"""

optionParser = optparse.OptionParser(usage ="""\
usage: %prog [options]
    Run test_bestfit_point_shift_due_to_syst_twk.C""")
optionParser.add_option("--resource", action="store", type="string", dest="resource", default="Oxford", help="RAL,IN2P3, Oxford, Grid, interactive. Default Oxford")
optionParser.add_option("--cputime", action="store", type="string", dest="cputime", default="11:50:00", help="Default: 11:50:00")
optionParser.add_option("--runmode", action="store", type="string", dest="runmode", default="", help="\nPULLS (10000 fits at 2 oscillation points. \nIncludes systematics fluctuations if applicable), STARSTAT (10000 fits at 9 oscillation points. No systematic fluctuations. Only fits for osc params")
optionParser.add_option("--fitparamlist", action="store", type="string", dest="fitparamlist", default="nuosc3f_numudisapp_fix4", help="Default: nuosc3f_numudisapp_fix4 (this is ignored in run mode 'starstat'")
optionParser.add_option("--minuit", action="store_true", dest="minuit", default="false", help="Run with Minuit. Otherwise, use Minuit2")
(opt, args) = optionParser.parse_args()

#will create directories relative to the current working directory, if they don't already exist
basedir = os.getenv("PWD")

#initialise some variables
fitstr = ""
if not opt.minuit:
   fitstr = "--minuit2"
conffile = "analysis_configuration.xml"
oscparams = []
syst = -1
tag = ""
nruns = 1
nexpts = 1
stat = ""
store_chisqsurf = ""
store_chisq_withsyst_surf = ""
store_toy_expts = ""
store_bestfit_spectra = ""

resource = opt.resource
cputime = opt.cputime
runmode = opt.runmode
fitparamlist = opt.fitparamlist

if not runmode:
   print "please specify a run mode!"
   sys.exit()

if runmode == "PULLS":
   #fill with oscillation parameters, in order mixing12, mixing13, mixing23, deltacp, dmsq12, dmsq23
   oscparams.append([-99, -99, 0.9, -99, -99, -99]) #mixing23=0.9 & global best fits
   oscparams.append([-99, -99, -99, -99, -99, -99]) #use global best fits!
   stat = " --stat-fluctuations "
   tag = "pull"
   if fitparamlist == "nuosc3f_numudisapp_fix4":
      syst = 0
   elif fitparamlist == "nuosc3f_numudisapp_fix4+nuisanceselection":
      syst = 6
   else:
      syst = 1
   if fitparamlist == "nuosc3f_numudisapp_fix4" or fitparamlist == "nuosc3f_numudisapp_fix4+nuisanceselection" or fitparamlist == "nuosc3f_numudisapp_fix4+nuisanceunbanffselection" or fitparamlist == "nuosc3f_numudisapp_fix4+nuisanceselectionnobanff":
      nruns = 100
      nexpts = 100
   else:
      nruns = 200
      nexpts = 50

if runmode == "STARSTAT":
   #fill with oscillation parameters, in order mixing12, mixing13, mixing23, deltacp, dmsq12, dmsq23
   oscparams.append([-99, -99, -99, -99, -99, -99]) #use global best fits!
   oscparams.append([-99, -99, 0.9, -99, -99, -99]) #mixing23=0.9 & global best fits
   oscparams.append([-99, -99, 0.95, -99, -99, -99]) #mixing23=0.95 & global best fits
   oscparams.append([-99, -99, -99, -99, -99, 0.0022]) #use global best fits! & dm32sq = 2.2E-3
   oscparams.append([-99, -99, 0.9, -99, -99, 0.0022]) #mixing23=0.9 & global best fits & dm32sq = 2.2E-3
   oscparams.append([-99, -99, 0.95, -99, -99, 0.0022]) #mixing23=0.95 & global best fits & dm32sq = 2.2E-3
   oscparams.append([-99, -99, -99, -99, -99, 0.0026]) #use global best fits! & dm32sq = 2.6E-3
   oscparams.append([-99, -99, 0.9, -99, -99, 0.0026]) #mixing23=0.9 & global best fits & dm32sq = 2.6E-3
   oscparams.append([-99, -99, 0.95, -99, -99, 0.0026]) #mixing23=0.95 & global best fits & dm32sq = 2.6E-3
   syst = 0
   tag = "starstat"
   nruns = 100
   nexpts = 100
   fitparamlist = "nuosc3f_numudisapp_fix4"
   stat = " --stat-fluctuations "


for osc in oscparams:

   if not os.path.isdir(fitparamlist):
      try:
         os.makedirs(fitparamlist)
      except:
         print "Could not make directory: ", fitparamlist
         sys.exit()
   os.chdir(fitparamlist)

   dirname = 'pulls_'+str(osc[0])+'_'+str(osc[1])+'_'+str(osc[2])+'_'+str(osc[3])+'_'+str(osc[4])+'_'+str(osc[5])
   if not os.path.isdir(dirname):
      try:
         os.makedirs(dirname)
      except:
         print "Could not make directory: ", dirname
         sys.exit()
   os.chdir(dirname)

   command = '$T2KOSCANA/production/submit_toyexpt_fitting_job_fitsyst_minuit2 --resource '+resource+' --softw-dir $T2KOSCANA --job-dir $PWD --cpu-time '+cputime+' --tweak-systematics '+str(syst)+' --nruns '+str(nruns)+' --nexpts '+str(nexpts)+' --sinsq2theta12 '+str(osc[0])+' --sinsq2theta13 '+str(osc[1])+' --sinsq2theta23 '+str(osc[2])+' --deltacp '+str(osc[3])+' --dm12sq '+str(osc[4])+' --dm23sq '+str(osc[5])+' --fit-param-list '+fitparamlist+' --config-file '+conffile+' --tag '+tag+' '+store_chisqsurf+' '+store_chisq_withsyst_surf+' '+store_toy_expts+' '+store_bestfit_spectra+' '+stat+' '+fitstr

   print command
   os.system(command)

   os.chdir(basedir)

"""
./common_runs.py --runmode=PULLS --fitparamlist=nuosc3f_numudisapp_fix4
./common_runs.py --runmode=PULLS --fitparamlist=nuosc3f_numudisapp_fix4+nuisanceunbanffselection
./common_runs.py --runmode=PULLS --fitparamlist=nuosc3f_numudisapp_fix4+nuisanceselection
./common_runs.py --runmode=PULLS --fitparamlist=nuosc3f_numudisapp_fix4+nuisancefullbanffselection
"""
