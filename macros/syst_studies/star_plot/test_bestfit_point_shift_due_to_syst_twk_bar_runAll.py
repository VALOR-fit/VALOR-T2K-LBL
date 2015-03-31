#! /usr/bin/env python

import sys
import optparse
import os
import math # cos, log, pi etc
import string
import logging
import ROOT
import operator

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

env_var = "T2KOSCANA"
try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)
  
valor = ROOT.valor
util    =  valor.util
conv    =  valor.conv
opt     =  valor.Opt.Instance()
t2kopt  =  valor.T2KOpt.Instance()
env     =  valor.Env.Instance()
rnd     =  valor.Rnd.Instance()
arridx  = util.ArrIdx.Instance()

param_of_interest_1 = "mixing23bar"
param_of_interest_2 = "dmsq23bar"

def main(options):

  conffile = options.config
  #will create directories relative to the current working directory, if they don't already exist
  basedir = os.getenv("PWD")

  fitparamlist = options.fpl
  fitter = "minuit"
  
  full_conf = os.environ[env_var]+"/config/"+conffile
  t2kopt.ReadXmlPath(full_conf)
  fparm = util.NamedParameterList(options.fpl)
  
  samples = []
  samples.append(options.sample)
  
  systematics = ["nominal"]
  for p in xrange(conv.knppmax):
    if fparm.Exists(p) and fparm.IsNuisance(p) and not fparm.IsFixed(p):
      systematics.append(util.AsStringSimple(p))
  
  if options.do_syst:
    print "List of systematics: "
    for syst in systematics:
      print syst
     
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
  oscparams.append([-99, -99, 0.513, -99, -99, 0.0022])   #dm32sq=0.0022 & global best fits!
  oscparams.append([-99, -99, 0.513, -99, -99, 0.0025])   #use global best fits!
  oscparams.append([-99, -99, 0.513, -99, -99, 0.0028])   #dm32sq=0.0028 & global best fits!

  oscparams.append([-99, -99, 0.5,  -99, -99, 0.0022])   #dm32sq=0.0022 & mixing23=0.5 & global best fits
  oscparams.append([-99, -99, 0.5,  -99, -99, 0.0025])   #dm32sq=0.0025 & mixing23=0.5 & global best fits
  oscparams.append([-99, -99, 0.5,  -99, -99, 0.0028])   #dm32sq=0.0028 & mixing23=0.5 & global best fits

  oscparams.append([-99, -99, 0.37,  -99, -99, 0.0022])   #dm32sq=0.0022 & mixing23=0.37 & global best fits
  oscparams.append([-99, -99, 0.37,  -99, -99, 0.0025])   #mixing23=0.37 & global best fits
  oscparams.append([-99, -99, 0.37,  -99, -99, 0.0028])   #dm32sq=0.0028 & mixing23=0.37 & global best fits

  oscparams.append([-99, -99, 0.65,  -99, -99, 0.0022])   #dm32sq=0.0022 & mixing23=0.65 & global best fits
  oscparams.append([-99, -99, 0.65,  -99, -99, 0.0025])   #mixing23=0.65 & global best fits
  oscparams.append([-99, -99, 0.65,  -99, -99, 0.0028])   #dm32sq=0.0028 & mixing23=0.65 & global best fits
  """
  
  # Numbers for 2015 numubar disappearance
  oscparams.append([-99, -99, 0.524, -99, -99, 0.0020])   #dm32sq=0.0022 & global best fits!
  oscparams.append([-99, -99, 0.524, -99, -99, 0.0025])   #use global best fits!
  oscparams.append([-99, -99, 0.524, -99, -99, 0.0030])   #dm32sq=0.0028 & global best fits!

  oscparams.append([-99, -99, 0.5  , -99, -99, 0.0020])   #dm32sq=0.0022 & mixing23=0.5 & global best fits
  oscparams.append([-99, -99, 0.5  , -99, -99, 0.0025])   #dm32sq=0.0025 & mixing23=0.5 & global best fits
  oscparams.append([-99, -99, 0.5  , -99, -99, 0.0030])   #dm32sq=0.0028 & mixing23=0.5 & global best fits

  oscparams.append([-99, -99, 0.37 , -99, -99, 0.0020])   #dm32sq=0.0022 & mixing23=0.37 & global best fits
  oscparams.append([-99, -99, 0.37 , -99, -99, 0.0025])   #mixing23=0.37 & global best fits
  oscparams.append([-99, -99, 0.37 , -99, -99, 0.0030])   #dm32sq=0.0028 & mixing23=0.37 & global best fits

  oscparams.append([-99, -99, 0.65 , -99, -99, 0.0020])   #dm32sq=0.0022 & mixing23=0.65 & global best fits
  oscparams.append([-99, -99, 0.65 , -99, -99, 0.0025])   #mixing23=0.65 & global best fits
  oscparams.append([-99, -99, 0.65 , -99, -99, 0.0030])   #dm32sq=0.0028 & mixing23=0.65 & global best fits
  
  nonbar_osc = [-99,-99,-99,-99,-99,-99]
  
  dirnames = []
  
  for osc in oscparams:
    
    for sam in samples:
      
      dirname = 'bestfit_point_shift_'+sam+'_'+str(osc[0])+'_'+str(osc[1])+'_'+str(osc[2])+'_'+str(osc[3])+'_'+str(osc[4])+'_'+str(osc[5])
      if not os.path.isdir(dirname):
        try:
          os.makedirs(dirname)
        except:
          sys.exit("Could not make installation directory")
          
      dirnames.append(dirname)

      if options.do_staterr:
        file_names = dirname+"/star_plot_stat*root"
        stat_err_args = [dirname,param_of_interest_1,param_of_interest_2]
        stat_err_args_string = ",".join(macro_prep(stat_err_args))
        command = ('root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C('+stat_err_args_string+')\'' )
        print command
        os.system(command)
        continue

      os.chdir(dirname)
      
      if options.do_stat:
        stat_only_fitparam_list = fitparamlist.replace("+nuisance","")
        farmname = "RAL"
        if options.isoxford:
          farmname = "Oxford"
        command = ('"${T2KOSCANA}/production/submit_run_fit_xml_toys_joint_bar" --resource RAL '+
           ' --tag star_plot_stat --config-file '+conffile+' --job-dir "${PWD}" --softw-dir "${T2KOSCANA}"'+
           ' --nruns 25 --firstexpt 1 --lastexpt 400 --xml-file "dummy.xml" --use-xml-file false --dotoys true '+
           ' --input-mixing12 '+`nonbar_osc[0]`+' --input-mixing13 '+`nonbar_osc[1]`+' --input-mixing23 '+`nonbar_osc[2]`+
           ' --input-dm23sq '  +`nonbar_osc[5]`+' --input-dm12sq '  +`nonbar_osc[4]`+' --input-deltacp ' +`nonbar_osc[3]`+
           ' --input-mixing12-bar '+`osc[0]`   +' --input-mixing13-bar '+`osc[1]`   +' --input-mixing23-bar '+`osc[2]`+
           ' --input-dm23sq-bar '  +`osc[5]`   +' --input-dm12sq-bar '  +`osc[4]`   +' --input-deltacp-bar ' +`osc[3]`+
           ' --stat-fluctuations true --tweak-systematics false --fit-param-list '+stat_only_fitparam_list +
           ' --fit-option "minuit" --store-toy-expts true --store-bestfit-spectra true '
           ' --dcp-value '+`osc[3]`+' --hierarchy normal --reactor false --rethrow-unphysical false')
        print command
        os.system(command)
        os.chdir(basedir)
        continue
      #end do stat

      for syst in systematics:
        
        if options.do_syst:
          rootcommand = ( 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/test_bestfit_point_shift_due_to_syst_twk.C("'+
            syst+'","'+sam+'","'+fitparamlist+'","'+fitter+'",'+str(osc[0])+','+str(osc[1])+','+str(osc[2])+','+str(osc[3])+','+str(osc[4])+','+str(osc[5])+',"'+conffile )
          if options.interactive:
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
            if options.isoxford:
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
            if options.interactive:
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
              if options.isoxford:
                f.write("source "+os.environ['T2KOSCANA']+"/setup/setup_at_oxford.sh "+os.environ['T2KOSCANA']+" \n")
              else:
                f.write("source "+os.environ['T2KOSCANA']+"/setup/setup_at_ral.sh "+os.environ['T2KOSCANA']+" \n")
              f.write(command+"\n")
              f.close()
              submit = "qsub "+jobname+".pbs"
              print submit
              os.system(submit)
        #end do_syst
        
        
        if syst == systematics[len(systematics)-1] and any([options.do_draw,options.do_table]): # On last iteration, make combined files
          command = 'rm combined*.data'
          print command
          os.system(command)
          
          command = 'cat star_plot_line-*.data > combined-all.data'
          print command
          os.system(command)
          
          command = 'cat star_plot_line-[a-z]*.data > combined-alluncorrelated.data'
          print command
          os.system(command)
          
          command = 'cat star_plot_line-osc*.data > combined-osc.data'
          print command
          os.system(command)
          
          command = 'cat star_plot_line-_correlated*.data star_plot_line-[n,x]*.data star_plot_line-sk_energy_scale.data > combined-allcorrelated.data'
          print command
          os.system(command)
          
          #~ command = 'cat star_plot_line-sk_*.data > combined-sk.data'
          #~ print command
          #~ os.system(command)
          
          command = 'cat star_plot_line-xsec_*.data > combined-xsec.data'
          print command
          os.system(command)
          
          command = 'cat `ls *data | grep \'star_plot_line-\\w*banff\' | grep -v "banff[0-9][0-9]"`  > combined-banffxsec.data'
          print command
          os.system(command)
          
          command = 'cat `ls *data | grep \'star_plot_line-\\w*banff\' | grep "banff[0-9][0-9]"`  > combined-banffflux.data'
          print command
          os.system(command)
          
          command = 'cat star_plot_line*correlated*banff*.data cat > combined-correlated_banff.data'
          print command
          os.system(command)
          
          #~ command = 'cat star_plot_line-fsi_*.data > combined-fsi.data'
          #~ print command
          #~ os.system(command)
          
          command = 'cat star_plot_line-*skdetfsi_*.data > combined-skdetfsi.data'
          print command
          os.system(command)
          
        if syst == systematics[len(systematics)-1] and options.do_draw:

          file = open('stat.data', 'r')
          errmixing23 = file.readline().strip('\n')
          errdm23sq   = file.readline().strip('\n')
          file.close()
          
          command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-all.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
          print command
          os.system(command)
          
          command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-alluncorrelated.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
          print command
          os.system(command)
          
          command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-allcorrelated.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
          print command
          os.system(command)
          
          command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-osc.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
          print command
          os.system(command)
          
          command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-sk.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
          print command
          os.system(command)
          
          command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-xsec.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
          print command
          os.system(command)

          #~ command = 'for i in {00..10}; do printf -v ii "%02i" $i; cat star_plot_line-banff_$ii* >> combined-banffnumu.data; done'
          #~ print command
          #~ os.system(command)
          #~ command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-banffnumu.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
          #~ print command
          #~ os.system(command)

          #~ command = 'for i in {11..21}; do printf -v ii "%02i" $i; cat star_plot_line-banff_$ii* >> combined-banffnumubar.data; done'
          #~ print command
          #~ os.system(command)
          #~ command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-banffnumubar.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
          #~ print command
          #~ os.system(command)
          
          command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-banffxsec.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
          print command
          os.system(command)
          
          command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-correlated_banff.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
          print command
          os.system(command)
          
          command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-fsi.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
          print command
          os.system(command)

          #~ command = 'cat star_plot_line-sk_1Rmu_ccnonqe_eff.data star_plot_line-sk_1Rmu_nc_eff.data star_plot_line-sk_energy_scale.data star_plot_line-xsec_pilessdelta.data star_plot_line-xsec_sf.data star_plot_line-xsec_bindE.data star_plot_line-banff_maqe.data star_plot_line-banff_ccqeE0.data > combined-important.data'
          #~ print command
          #~ os.system(command)
          #~ command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/make_star_plot_for_disappearance_analysis.C("combined-important.data",'+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+')\''
          #~ print command
          #~ os.system(command)
        #end do draw
          
        if syst == systematics[len(systematics)-1]:
          os.chdir(basedir)
          
      #end systematic loop
      
      #~ if options.do_table:
        #~ os.chdir(basedir)
        #~ file = open(dirname+'/stat.data', 'r')
        #~ errmixing23 = file.readline().strip('\n')
        #~ errdm23sq   = file.readline().strip('\n')
        #~ file.close()
        #~ command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/text_to_root.C('+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+',"'+dirname+'/combined-all.data")\''
        #~ print command
        #~ os.system(command)
        #~ command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/text_to_root.C('+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+',"'+dirname+'/combined-alluncorrelated.data","uncorr")\''
        #~ print command
        #~ os.system(command)
        #~ command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/text_to_root.C('+str(osc[2])+','+str(osc[5])+','+errmixing23+','+errdm23sq+',"'+dirname+'/combined-allcorrelated.data","corr")\''
        #~ print command
        #~ os.system(command)

  # end oscparams loop
  
  #~ if options.do_table:
    #~ command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/max_deviation_table.C("readfile*.root")\''
    #~ print command
    #~ os.system(command)
    #~ command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/max_deviation_table.C("corr*.root")\''
    #~ print command
    #~ os.system(command)
    #~ command = 'root2k -l -b -q \'$T2KOSCANA/macros/syst_studies/star_plot/max_deviation_table.C("uncorr*.root")\''
    #~ print command
    #~ os.system(command)
    
  if options.do_table:
    make_star_plot_table(dirnames,"stat.data","combined-all.data")
    
#________________________________________________________________________________________________________
def make_star_plot_table(directories, statfile_name, systfile_name):
  param_order = [conv.kPhypMixing23,conv.kPhypMixing13, conv.kPhypMixing12,
                 conv.kPhypCPV, conv.kPhypDm23sq, conv.kPhypDm12sq,
                 conv.kPhypMixing23Bar,conv.kPhypMixing13Bar, conv.kPhypMixing12Bar,
                 conv.kPhypCPVBar, conv.kPhypDm23sqBar, conv.kPhypDm12sqBar]
  
  p1 = util.PhysParamIdFromString(param_of_interest_1)
  p2 = util.PhysParamIdFromString(param_of_interest_2)
  
  deviations_p1 = dict()
  deviations_p2 = dict()
  for d in directories:
    deviations_p1[d] = dict()
    deviations_p2[d] = dict()
    nominals = dict()
    bf_per_syst_parameter = dict()
    statfile = open(d+'/'+statfile_name, 'r')
    err_p1 = float(statfile.readline().strip('\n'))
    err_p2 = float(statfile.readline().strip('\n'))
    statfile.close()
    
    systfile = open(d+'/'+systfile_name, 'r')
    for l in systfile:
      if not l: continue
      items = l.split("\t")
      if items[0] == "nominal":
        for ip,p in enumerate(param_order):
          nominals[p] = float(items[ip+1])
      else:
        syst_p = util.PhysParamIdFromString(items[0])
        if syst_p == conv.kPhypUndefined:
          print "Bad row:",l
          exit(-1)
          
        bf_per_syst_parameter[syst_p] = dict()
        for ip,p in enumerate(param_order):
          bf_per_syst_parameter[syst_p][p] = float(items[ip+1])
    
    syst_params = deviations_p1[d].keys()
    for syst_p in syst_params:
      deviations_p1[d][syst_p] = fabs((nominals[p1]-bf_per_syst_parameter[syst_p][p1])/err_p1)
      deviations_p2[d][syst_p] = fabs((nominals[p2]-bf_per_syst_parameter[syst_p][p2])/err_p2)
  
  max_deviations_p1 = dict()
  max_deviations_p2 = dict()
  for syst_p in syst_params:
    max_deviations_p1[syst_p] = max([deviations_p1[d][syst_p] for d in directories])
    max_deviations_p2[syst_p] = max([deviations_21[d][syst_p] for d in directories])
    
  rows = []
  for syst_p in syst_params:
    rows.append([util.ParamAsStringRooTex(syst_p),max_deviations_p1[syst_p],max_deviations_p2[syst_p]])
    
  rows_sorted_by_p1_dev = sorted(rows,operator.itemgetter(1),reverse=True)
  rows_sorted_by_p2_dev = sorted(rows,operator.itemgetter(2),reverse=True)
  
  out_p1 = open("star_plot_table_sorted_by_"+util.ParamAsStringSimple(p1)+".txt","w")
  for r in rows_sorted_by_p1_dev:
    out_p1.write(" & ".join([r[0].ljust(30),to3dp.format(r[1]),to3dp.format(r[2])]))
    out_p1.write(" \\\\ \hline \n")
  out_p1.close()
    
  out_p2 = open("star_plot_table_sorted_by_"+util.ParamAsStringSimple(p2)+".txt","w")
  for r in rows_sorted_by_p2_dev:
    out_p2.write(" & ".join([r[0].ljust(30),to3dp.format(r[1]),to3dp.format(r[2])]))
    out_p2.write(" \\\\ \hline \n")
  out_p2.close()
#________________________________________________________________________________________________________
def macro_prep(args):
  formatted_args = []
  for arg in args:
    s = str(arg)
    if s==arg:
      s='"'+arg+'"'
    formatted_args.append(s)
  return formatted_args 
#________________________________________________________________________________________________________
if __name__=="__main__":
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
  optionParser.add_option('-s', '--sample'        , dest='sample' , default="RHC1Rmu"              , help="Sample (1Rmu, 1Re, RHC1Rmu, RHC1Re)")
  optionParser.add_option('-c', '--config'        , dest='config' , default="analysis_configuration.xml", help="Configuration file")
  optionParser.add_option('-f', '--fit-param-list', dest='fpl'    , default="nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance", help="Fit parameter list")
  
  (options, args) = optionParser.parse_args()
  
  print options

  if not (options.do_syst or options.do_stat or options.do_draw or options.do_table or options.do_staterr):
    optionParser.print_usage()
    sys.exit(1)
     
  if (options.do_syst and options.do_stat):
    optionParser.print_usage()
    sys.exit(1)
  
  if args:
    optionParser.print_usage()
    print "Unknown arguments: ",args
     
  main(options)
