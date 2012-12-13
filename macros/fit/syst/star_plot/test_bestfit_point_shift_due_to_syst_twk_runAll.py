#! /usr/bin/env python

import sys
import optparse
import os
import math # cos, log, pi etc
import string
import logging

"""
This script will run the macro test_bestfit_point_shift_due_to_syst_twk.C for all systematic possiblities, producing
lots of .data files with the outputs, and a variety of .pdf images. Things you can change are:
 - conffile - configuration file
 - samples - the SK samples you want to consider (default: 1Rmu)
 - fitparamlist - the FitParamList you want to use (default: nuosc3f_numudisapp_fix4)
 - fitter - the fitter to use (default: minuit)
 - oscparam - the oscillation parameters you want to calculate shifts relative to (defaults to global best fits - you can specify more than 1 set of parameters to run over)
 - systematics - the systematics you want to run over (defaults to all currently known)

Run this by /path/to/test_bestfit_point_shift_due_to_syst_twk_runAll.py (you may need to chmod +x)
It will then create the .pdf & .data files in new subfolders from the current working directory
"""

optionParser = optparse.OptionParser(usage ="""\
usage: %prog [options]
    Run test_bestfit_point_shift_due_to_syst_twk.C""")
optionParser.add_option("--syst", action="store_true", dest="do_syst", default=False, help="If set, will produce .data files for each systematic")
optionParser.add_option("--stat", action="store_true", dest="do_stat", default=False, help="Stat is not implemented in this macro. Please run the production script, and use the code snippet at the bottom of $T2KOSCANA/macros/fit/syst/star_plot/test_bestfit_point_shift_due_to_syst_twk.C to get stat errors, then replace the numbers instat.data: $T2KOSCANA/production/submit_toyexpt_fitting_job_fitsyst --resource Oxford --softw-dir $T2KOSCANA --job-dir $PWD --cpu-time 01:59:00 --tweak-systematics 1 --nruns 100 --nexpts 100 --dm23sq 0.0024 --sinsq2theta23 1. --stat-fluctuations --fit-param-list nuosc3f_numudisapp_fix4")
optionParser.add_option("--draw", action="store_true", dest="do_draw", default=False, help="If set, will draw all the plots")
(opt, args) = optionParser.parse_args()

if not (opt.do_syst or opt.do_stat or opt.do_draw):
   optionParser.print_usage()
   sys.exit(1)

if opt.do_stat:
   print "Stat is not implemented in this macro. Please run the production script, and use the code snippet at the bottom of $T2KOSCANA/macros/fit/syst/star_plot/test_bestfit_point_shift_due_to_syst_twk.C to get stat errors, then replace the numbers in stat.data: $T2KOSCANA/production/submit_toyexpt_fitting_job_fitsyst --resource Oxford --softw-dir $T2KOSCANA --job-dir $PWD --cpu-time 01:59:00 --tweak-systematics 1 --nruns 100 --nexpts 100 --dm23sq 0.0024 --sinsq2theta23 1. --stat-fluctuations --fit-param-list nuosc3f_numudisapp_fix4"
   sys.exit(1)
   
if (opt.do_syst and opt.do_stat):
   optionParser.print_usage()
   sys.exit(1)

conffile = "analysis_configuration.xml"

#will create directories relative to the current working directory, if they don't already exist
basedir = os.getenv("PWD")

fitparamlist = "nuosc3f_numudisapp_fix4"
fitter = "minuit"

#specify the cuts you want to apply here
systematics = []
#oscillations
systematics.append("osc_mixing12");
systematics.append("osc_mixing13");
systematics.append("osc_dm12sq");
systematics.append("osc_deltacp");
#BANFF
#systematics.append("banff")
#BANFF flux
for i in xrange(16):
   systematics.append("banff_" + '%(number)02d' % {"number": i})
#BANFF xsec
systematics.append("banff_maqe")
systematics.append("banff_mares")
systematics.append("banff_ccqeE0")
systematics.append("banff_ccqeE1")
systematics.append("banff_ccqeE2")
systematics.append("banff_cc1piE0")
systematics.append("banff_cc1piE1")
#SK
#systematics.append("sk_1Rmu_ccqe_ring_counting_eff")
systematics.append("sk_1Rmu_ccqe_ring_counting_eff_0")
systematics.append("sk_1Rmu_ccqe_ring_counting_eff_1")
systematics.append("sk_1Rmu_ccqe_ring_counting_eff_2")
systematics.append("sk_1Rmu_ccqe_other_eff")
systematics.append("sk_1Rmu_ccnonqe_eff")
systematics.append("sk_1Rmu_nc_eff")
systematics.append("sk_1Rmu_nuecc_eff")
systematics.append("sk_energy_scale")
#XSec splines
systematics.append("xsec_pf")
systematics.append("xsec_sf")
systematics.append("xsec_wshape")
systematics.append("xsec_shape_ccoth")
systematics.append("xsec_pilessdelta")
systematics.append("xsec_bindE")
#XSec normalisations
systematics.append("xsec_norm_ccqeE1")
systematics.append("xsec_norm_ccqeE2")
systematics.append("xsec_norm_cc1piE1")
systematics.append("xsec_norm_cccoh")
systematics.append("xsec_norm_nc1piC")
systematics.append("xsec_norm_ncoth")
systematics.append("xsec_norm_nuetonumu")
systematics.append("xsec_norm_antinutonu")
#FSI
for i in xrange(6):
   systematics.append("fsi_" + '%(number)03d' % {"number": i})
"""
#UnBANFF
#systematics.append("unbanff")
for i in xrange(24):
   systematics.append("unbanff_" + '%(number)02d' % {"number": i})

#EDBANFF
#systematics.append("edbanff")
for i in xrange(24):
   systematics.append("edbanff_" + '%(number)02d' % {"number": i})
"""
   
samples = []
samples.append("1Rmu")
#samples.append("1Re")

#fill with oscillation parameters, in order mixing12, mixing13, mixing23, deltacp, dmsq12, dmsq23
oscparams = [] 
oscparams.append([-99, -99, -99, -99, -99, -99]) #use global best fits!
oscparams.append([-99, -99, 0.9, -99, -99, -99]) #mixing23=0.9 & global best fits
oscparams.append([-99, -99, 0.95, -99, -99, -99]) #mixing23=0.9 & global best fits
oscparams.append([-99, -99, -99, -99, -99, 0.0022]) #use global best fits! & dm32sq = 2.2E-3
oscparams.append([-99, -99, 0.9, -99, -99, 0.0022]) #mixing23=0.9 & global best fits & dm32sq = 2.2E-3
oscparams.append([-99, -99, 0.95, -99, -99, 0.0022]) #mixing23=0.95 & global best fits & dm32sq = 2.2E-3
oscparams.append([-99, -99, -99, -99, -99, 0.0026]) #use global best fits! & dm32sq = 2.6E-3
oscparams.append([-99, -99, 0.9, -99, -99, 0.0026]) #mixing23=0.9 & global best fits & dm32sq = 2.6E-3
oscparams.append([-99, -99, 0.95, -99, -99, 0.0026]) #mixing23=0.95 & global best fits & dm32sq = 2.6E-3
#some at the edges of the space
oscparams.append([-99, -99, 0.8, -99, -99, 0.0020])
oscparams.append([-99, -99, 0.8, -99, -99, 0.0030])
oscparams.append([-99, -99, 0.8, -99, -99, 0.0024])
oscparams.append([-99, -99, 1.0, -99, -99, 0.0020])
oscparams.append([-99, -99, 1.0, -99, -99, 0.0030])

for osc in oscparams:

   for sam in samples:
      
      dirname = 'bestfit_point_shift_'+sam+'_'+str(osc[0])+'_'+str(osc[1])+'_'+str(osc[2])+'_'+str(osc[3])+'_'+str(osc[4])+'_'+str(osc[5])
      if not os.path.isdir(dirname):
         try:
            os.makedirs(dirname)
         except:
            sys.exit("Could not make installation directory")

      os.chdir(dirname)

      for syst in systematics:

         if opt.do_syst:
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/test_bestfit_point_shift_due_to_syst_twk.C("'+syst+'","'+sam+'","'+fitparamlist+'","'+fitter+'",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+',"'+conffile+'",false)\''
            print command
            os.system(command)

            if syst.find("banff_") == 0:
               command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/test_bestfit_point_shift_due_to_syst_twk.C("'+syst+'","'+sam+'","'+fitparamlist+'","'+fitter+'",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+',"'+conffile+'",true)\''
               print command
               os.system(command)

         if syst == systematics[len(systematics)-1] and opt.do_draw:

            file = open('stat.data', 'r')
            errmixing23 = file.readline().strip('\n')
	    errdm23sq   = file.readline().strip('\n')
            file.close()

            command = 'rm star_plot_line-all.data star_plot_line-banff[n,x]*.data star_plot_line-edbanff.data star_plot_line-unbanff.data star_plot_line-fsi.data star_plot_line-sk.data star_plot_line-xsec.data star_plot_line-correlated_banff.data'
            print command
            os.system(command)
            command = 'cat star_plot_line-*.data > star_plot_line-all.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-all.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)
            command = 'cat star_plot_line-osc*.data > star_plot_line-osc.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-osc.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)
            command = 'cat star_plot_line-sk_*.data > star_plot_line-sk.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-sk.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-xsec_*.data > star_plot_line-xsec.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-xsec.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'for i in {00..10}; do printf -v ii "%02i" $i; cat star_plot_line-banff_$ii* >> star_plot_line-banffnumu.data; done'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-banffnumu.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'for i in {11..21}; do printf -v ii "%02i" $i; cat star_plot_line-banff_$ii* >> star_plot_line-banffnumubar.data; done'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-banffnumubar.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            """
            command = 'for i in {22..28}; do printf -v ii "%02i" $i; cat star_plot_line-banff_$ii* >> star_plot_line-banffnue.data; done'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-banffnue.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)
            
            command = 'for i in {29..30}; do printf -v ii "%02i" $i; cat star_plot_line-banff_$ii* >> star_plot_line-banffnuebar.data; done'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-banffnuebar.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)
            """

            command = 'cat star_plot_line-banff_ma*.data cat star_plot_line-banff_cc*.data > star_plot_line-banffxsec.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-banffxsec.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-unbanff*.data cat > star_plot_line-unbanff.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-unbanff.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-correlated_banff*.data cat > star_plot_line-correlated_banff.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-correlated_banff.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-edbanff*.data cat > star_plot_line-edbanff.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-edbanff.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-fsi_*.data > star_plot_line-fsi.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-fsi.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-sk_1Rmu_ccnonqe_eff.data star_plot_line-sk_1Rmu_nc_eff.data star_plot_line-sk_energy_scale.data star_plot_line-xsec_pilessdelta.data star_plot_line-xsec_sf.data star_plot_line-xsec_bindE.data star_plot_line-banff_maqe.data star_plot_line-banff_ccqeE0.data > star_plot_line-important.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/fit/syst/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-important.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)
            
         if syst == systematics[len(systematics)-1]:

            os.chdir(basedir)
