#! /usr/bin/env python

import sys
import optparse
import os
import math # cos, log, pi etc
import string
import logging

"""
This script will run the macro test_bestfit_point_shift_due_to_syst_twk.C for all systematic possiblities, producing
lots of .data files with the outputs, and a variety of .png images. Things you can change are:
 - conffile - configuration file
 - samples - the SK samples you want to consider (default: 1Rmu)
 - fitparamlist - the FitParamList you want to use (default: nuosc3f_numudisapp_fix4)
 - fitter - the fitter to use (default: minuit)
 - oscparam - the oscillation parameters you want to calculate shifts relative to (defaults to global best fits - you can specify more than 1 set of parameters to run over)
 - systematics - the systematics you want to run over (defaults to all currently known)

Run this by /path/to/test_bestfit_point_shift_due_to_syst_twk_runAll.py (you may need to chmod +x)
It will then create the .png & .data files in new subfolders from the current working directory
"""

conffile = "analysis_configuration.xml"

#will create directories relative to the current working directory, if they don't already exist
basedir = os.getenv("PWD")

fitparamlist = "nuosc3f_numudisapp_fix4"
fitter = "minuit"

#specify the cuts you want to apply here
systematics = []
#stats
systematics.append("stat");
#oscillations
systematics.append("osc_mixing12");
systematics.append("osc_mixing13");
systematics.append("osc_dm12sq");
systematics.append("osc_deltacp");
#SK
"""
#systematics.append("sk_1Rmu_ccqe_ring_counting_eff")
systematics.append("sk_1Rmu_ccqe_ring_counting_eff_0")
systematics.append("sk_1Rmu_ccqe_ring_counting_eff_1")
systematics.append("sk_1Rmu_ccqe_ring_counting_eff_2")
systematics.append("sk_1Rmu_ccqe_other_eff")
"""
systematics.append("sk_1Rmu_ccnonqe_eff")
systematics.append("sk_1Rmu_nc_eff")
"""
systematics.append("sk_1Rmu_nuecc_eff")
systematics.append("sk_energy_scale")
#BANFF
#systematics.append("banff")
#BANFF flux
for i in xrange(31):
   systematics.append("banff_" + '%(number)02d' % {"number": i})
#BANFF xsec
"""
systematics.append("banff_maqe")
"""
systematics.append("banff_mares")
"""
systematics.append("banff_ccqeE0")
"""
systematics.append("banff_ccqeE1")
systematics.append("banff_ccqeE2")
"""
systematics.append("banff_cc1piE0")
"""
systematics.append("banff_cc1piE1")
#XSec splines
systematics.append("xsec_pf")
"""
systematics.append("xsec_sf")
"""
systematics.append("xsec_wshape")
systematics.append("xsec_shape_ccoth")
"""
systematics.append("xsec_pilessdelta")
"""
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
for i in xrange(15):
   systematics.append("fsi_" + '%(number)03d' % {"number": i})
"""

samples = []
samples.append("1Rmu")
#samples.append("1Re")

#fill with oscillation parameters, in order mixing12, mixing13, mixing23, deltacp, dmsq12, dmsq23
oscparams = [] 
oscparams.append([-99, -99, -99, -99, -99, -99]) #use global best fits!
oscparams.append([-99, -99, 0.9, -99, -99, -99]) #mixing23=0.9 & global best fits

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
          
         command = 'root2k -l -b -q \'$T2KOSCANA/macros/misc/star_plot/test_bestfit_point_shift_due_to_syst_twk.C("'+syst+'","'+sam+'","'+fitparamlist+'","'+fitter+'",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+',"'+conffile+'")\''
         print command
         #os.system(command)

         if syst == systematics[len(systematics)-1]:

            file = open('star_plot_line-stat.data', 'r')
            file.readline()
            errmixing12 = file.readline().strip('\n')
            errmixing13 = file.readline().strip('\n')
            errmixing23 = file.readline().strip('\n')
	    errdeltacp  = file.readline().strip('\n')
	    errdm12sq   = file.readline().strip('\n')
	    errdm23sq   = file.readline().strip('\n')
            file.close()

            command = 'rm star_plot_line-all.data star_plot_line-banff[a-z]*.data star_plot_line-fsi.data star_plot_line-sk.data star_plot_line-xsec.data'
            print command
            os.system(command)

            command = 'cat star_plot_line-osc*.data > star_plot_line-osc.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/misc/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-osc.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)
            command = 'cat star_plot_line-*.data > star_plot_line-all.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/misc/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-all.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-sk_*.data > star_plot_line-sk.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/misc/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-sk.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-xsec_*.data > star_plot_line-xsec.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/misc/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-xsec.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'for i in {00..10}; do printf -v ii "%02i" $i; cat star_plot_line-banff_$ii* >> star_plot_line-banffnumu.data; done'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/misc/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-banffnumu.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'for i in {11..21}; do printf -v ii "%02i" $i; cat star_plot_line-banff_$ii* >> star_plot_line-banffnumubar.data; done'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/misc/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-banffnumubar.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'for i in {22..28}; do printf -v ii "%02i" $i; cat star_plot_line-banff_$ii* >> star_plot_line-banffnue.data; done'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/misc/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-banffnue.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'for i in {29..30}; do printf -v ii "%02i" $i; cat star_plot_line-banff_$ii* >> star_plot_line-banffnuebar.data; done'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/misc/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-banffnuebar.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-banff_ma*.data cat star_plot_line-banff_cc*.data > star_plot_line-banffxsec.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/misc/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-banffxsec.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-fsi_*.data > star_plot_line-fsi.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/misc/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-fsi.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-sk_1Rmu_ccnonqe_eff.data star_plot_line-sk_1Rmu_nc_eff.data star_plot_line-xsec_pilessdelta.data star_plot_line-xsec_sf.data star_plot_line-banff_maqe.data star_plot_line-banff_ccqeE0.data star_plot_line-banff_cc1piE0.data > star_plot_line-important.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/misc/star_plot/make_star_plot_for_disappearance_analysis.C("star_plot_line-important.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)
            
            
            os.chdir(basedir)
