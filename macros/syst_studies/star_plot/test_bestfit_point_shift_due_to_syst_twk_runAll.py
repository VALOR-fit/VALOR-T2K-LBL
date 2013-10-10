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
optionParser.add_option("--stat", action="store_true", dest="do_stat", default=False, help="If set, will submit lots of 2+0 fitting jobs with no systematic fluctuations, to get the stat error at every point")
optionParser.add_option("--staterr", action="store_true", dest="do_staterr", default=False, help="If set, will calculate the stat error at every point")
optionParser.add_option("--draw", action="store_true", dest="do_draw", default=False, help="If set, will draw all the plots")
optionParser.add_option("--table", action="store_true", dest="do_table", default=False, help="If set, will make the max deviation table")
optionParser.add_option("--interactive", action="store_true", dest="interactive", default=False, help="If set, will run +- 1 sigma jobs interactively; otherwise, throw 1000 toys")
optionParser.add_option("--batch-is-oxford", action="store_true", dest="isoxford", default=False, help="If set, will use Oxford setup for BATCH scripts. Otherwise, use RAL")
(opt, args) = optionParser.parse_args()

if not (opt.do_syst or opt.do_stat or opt.do_draw or opt.do_table or opt.do_staterr):
   optionParser.print_usage()
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
#no syst
systematics.append("nominal");
#oscillations
systematics.append("osc_mixing12");
systematics.append("osc_mixing13");
systematics.append("osc_dm12sq");
systematics.append("osc_deltacp");
#BANFF
#systematics.append("banff")
for i in xrange(30):
   systematics.append("banff_" + '%(number)02d' % {"number": i})
#SK energy scale
systematics.append("sk_energy_scale")
#SKDet FSI
for i in xrange(6):
   systematics.append("skdetfsi_" + '%(number)03d' % {"number": i})
"""
#SK
#systematics.append("sk_1Rmu_ccqe_ring_counting_eff")
systematics.append("sk_1Rmu_ccqe_ring_counting_eff_0")
systematics.append("sk_1Rmu_ccqe_ring_counting_eff_1")
systematics.append("sk_1Rmu_ccqe_ring_counting_eff_2")
systematics.append("sk_1Rmu_ccqe_other_eff")
systematics.append("sk_1Rmu_ccnonqe_eff")
systematics.append("sk_1Rmu_nc_eff")
systematics.append("sk_1Rmu_nuecc_eff")
#FSI
for i in xrange(6):
   systematics.append("fsi_" + '%(number)03d' % {"number": i})
"""
#XSec splines
systematics.append("xsec_pf")
systematics.append("xsec_sf")
systematics.append("xsec_wshape")
systematics.append("xsec_shape_ccoth")
systematics.append("xsec_pilessdelta")
systematics.append("xsec_bindE")
#XSec normalisations
systematics.append("xsec_norm_cccoh")
systematics.append("xsec_norm_nc1piC")
systematics.append("xsec_norm_ncoth")
systematics.append("xsec_norm_nuetonumu")
systematics.append("xsec_norm_antinutonu")
#MEC
systematics.append("mec")


samples = []
samples.append("1Rmu")
#samples.append("1Re")

#fill with oscillation parameters, in order mixing12, mixing13, mixing23, deltacp, dmsq12, dmsq23
oscparams = []
"""
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
#some experimental best fits
oscparams.append([0.857, 0.098, 1.000, 0, 7.5E-5, 0.00244]) #t2k run1-3 best fit
oscparams.append([0.857, 0.098, 0.957, 0, 7.5E-5, 0.00239]) #minos neutrino2012 best fit
"""
#using sinsqtheta
"""
oscparams.append([-99, -99, 0.513, -99, -99, 0.0025])   #use global best fits!
oscparams.append([-99, -99, 0.513, -99, -99, 0.0028])   #dm32sq=0.0028 & global best fits!
oscparams.append([-99, -99, 0.513, -99, -99, 0.0022])   #dm32sq=0.0022 & global best fits!

oscparams.append([-99, -99, 0.5,  -99, -99, 0.0022])   #dm32sq=0.0022 & mixing23=0.5 & global best fits
oscparams.append([-99, -99, 0.5,  -99, -99, 0.0025])   #dm32sq=0.0025 & mixing23=0.5 & global best fits
oscparams.append([-99, -99, 0.5,  -99, -99, 0.0028])   #dm32sq=0.0028 & mixing23=0.5 & global best fits
"""
oscparams.append([-99, -99, 0.37,  -99, -99, 0.0025])   #mixing23=0.37 & global best fits
oscparams.append([-99, -99, 0.37,  -99, -99, 0.0028])   #dm32sq=0.0028 & mixing23=0.37 & global best fits
oscparams.append([-99, -99, 0.37,  -99, -99, 0.0022])   #dm32sq=0.0022 & mixing23=0.37 & global best fits

oscparams.append([-99, -99, 0.65,  -99, -99, 0.0025])   #mixing23=0.65 & global best fits
oscparams.append([-99, -99, 0.65,  -99, -99, 0.0028])   #dm32sq=0.0028 & mixing23=0.65 & global best fits
oscparams.append([-99, -99, 0.65,  -99, -99, 0.0022])   #dm32sq=0.0022 & mixing23=0.65 & global best fits


for osc in oscparams:

   for sam in samples:
      
      dirname = 'bestfit_point_shift_'+sam+'_'+str(osc[0])+'_'+str(osc[1])+'_'+str(osc[2])+'_'+str(osc[3])+'_'+str(osc[4])+'_'+str(osc[5])
      if not os.path.isdir(dirname):
         try:
            os.makedirs(dirname)
         except:
            sys.exit("Could not make installation directory")

      if opt.do_staterr:
          command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("'+dirname+'")\''
          print command
          os.system(command)
          continue

      os.chdir(dirname)

      if opt.do_stat:
          farmname = "RAL"
          if opt.isoxford:
              farmname = "Oxford"
          command = '$T2KOSCANA/production/submit_toyexpt_fitting_job --resource '+farmname+' --dm12sq 7.5E-5 --mix12 0.306 --mix13 0.0251 --deltacp -999 --cpu-time 10:00:00 --fit-param-list nuosc3f_numudisapp_fix4 --fit-option minuit --job-dir $PWD --softw-dir $T2KOSCANA --nexpts 100 --nruns 100 --dmsq-scheme abs23 --hierarchy normal --mix-param-choice sinsqtheta --stat-fluctuations --tag stat_mix'+str(osc[2])+'_dm'+str(osc[5])+' --dm23sq '+str(osc[5])+' --mix23 '+str(osc[2])
          print command
          os.system(command)
          os.chdir(basedir)
          continue

      for syst in systematics:

         if opt.do_syst:
            rootcommand = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/test_bestfit_point_shift_due_to_syst_twk.C("'+syst+'","'+sam+'","'+fitparamlist+'","'+fitter+'",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+',"'+conffile
            if opt.interactive:
               command = rootcommand + '",false,0)\''
               print command
               os.system(command)
            else:
               command = rootcommand + '",false,0)\''
               print command
               jobname = syst
               f = open(jobname+".pbs",'w')
               f.write("#!/bin/bash \n")
               f.write("#PBS -N "+jobname+" \n")
               f.write("#PBS -l walltime=24:00:00 \n")
               f.write("#PBS -l cput=24:00:00 \n")
               f.write("#PBS -o "+os.getcwd()+"/"+jobname+".out \n")
               f.write("#PBS -e "+os.getcwd()+"/"+jobname+".err \n")
               f.write("cd "+os.getcwd() +"  \n")
               if opt.isoxford:
                   f.write("source "+os.environ['T2KOSCANA']+"/setup/setup_at_oxford.sh "+os.environ['T2KOSCANA']+" \n")
               else:
                   f.write("source "+os.environ['T2KOSCANA']+"/setup/setup_at_ral.sh "+os.environ['T2KOSCANA']+" \n")
               f.write(command+"\n")
               f.close()
               submit = "qsub "+jobname+".pbs"
               print submit
               os.system(submit)

            if syst.find("banff_") == 0 or syst.find("skdetfsi") == 0 or syst.find("fsi") == 0:
               print command
               if opt.interactive:
                  command = rootcommand + '",true,0)\''
                  print command
                  os.system(command)
               else:
                  command = rootcommand + '",true,0)\''
                  print command
                  jobname = syst+"_correlated"
                  f = open(jobname+".pbs",'w')
                  f.write("#!/bin/bash \n")
                  f.write("#PBS -N "+jobname+" \n")
                  f.write("#PBS -l walltime=24:00:00 \n")
                  f.write("#PBS -l cput=24:00:00 \n")
                  f.write("#PBS -o "+os.getcwd()+"/"+jobname+".out \n")
                  f.write("#PBS -e "+os.getcwd()+"/"+jobname+".err \n")
                  f.write("cd "+os.getcwd() +"  \n")
                  if opt.isoxford:
                      f.write("source "+os.environ['T2KOSCANA']+"/setup/setup_at_oxford.sh "+os.environ['T2KOSCANA']+" \n")
                  else:
                      f.write("source "+os.environ['T2KOSCANA']+"/setup/setup_at_ral.sh "+os.environ['T2KOSCANA']+" \n")
                  f.write(command+"\n")
                  f.close()
                  submit = "qsub "+jobname+".pbs"
                  print submit
                  os.system(submit)

         if syst == systematics[len(systematics)-1] and opt.do_draw:

            file = open('stat.data', 'r')
            errmixing23 = file.readline().strip('\n')
	    errdm23sq   = file.readline().strip('\n')
            file.close()

            command = 'rm combined*.data'
            print command
            os.system(command)
            command = 'cat star_plot_line-*.data > combined-all.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-all.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)
            command = 'cat star_plot_line-[a-z]*.data > combined-alluncorrelated.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-alluncorrelated.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)
            command = 'cat star_plot_line-_correlated*.data star_plot_line-[n,x]*.data star_plot_line-sk_energy_scale.data > combined-allcorrelated.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-allcorrelated.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)
            command = 'cat star_plot_line-osc*.data > combined-osc.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-osc.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)
            command = 'cat star_plot_line-sk_*.data > combined-sk.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-sk.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-xsec_*.data > combined-xsec.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-xsec.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'for i in {00..10}; do printf -v ii "%02i" $i; cat star_plot_line-banff_$ii* >> combined-banffnumu.data; done'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-banffnumu.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'for i in {11..21}; do printf -v ii "%02i" $i; cat star_plot_line-banff_$ii* >> combined-banffnumubar.data; done'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-banffnumubar.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-banff_ma*.data cat star_plot_line-banff_ccqeE*.data star_plot_line-banff_cc1piE* > combined-banffxsec.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-banffxsec.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line*correlated*banff*.data cat > combined-correlated_banff.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-correlated_banff.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)
            
            command = 'cat star_plot_line-fsi_*.data > combined-fsi.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-fsi.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)

            command = 'cat star_plot_line-sk_1Rmu_ccnonqe_eff.data star_plot_line-sk_1Rmu_nc_eff.data star_plot_line-sk_energy_scale.data star_plot_line-xsec_pilessdelta.data star_plot_line-xsec_sf.data star_plot_line-xsec_bindE.data star_plot_line-banff_maqe.data star_plot_line-banff_ccqeE0.data > combined-important.data'
            print command
            os.system(command)
            command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-important.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
            print command
            os.system(command)
            
         if syst == systematics[len(systematics)-1]:

            os.chdir(basedir)

            
      if opt.do_table:
         os.chdir(basedir)
         file = open(dirname+'/stat.data', 'r')
         errmixing23 = file.readline().strip('\n')
         errdm23sq   = file.readline().strip('\n')
         file.close()
         command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/text_to_root.C('+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+',"'+dirname+'/combined-all.data")\''
         print command
         os.system(command)
         command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/text_to_root.C('+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+',"'+dirname+'/combined-alluncorrelated.data","uncorr")\''
         print command
         os.system(command)
         command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/text_to_root.C('+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+',"'+dirname+'/combined-allcorrelated.data","corr")\''
         print command
         os.system(command)


if opt.do_table:
   command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/max_deviation_table.C("readfile*.root")\''
   print command
   os.system(command)
   command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/max_deviation_table.C("corr*.root")\''
   print command
   os.system(command)
   command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/max_deviation_table.C("uncorr*.root")\''
   print command
   os.system(command)
