#!/usr/bin/env python

# Runs a sensitivity job, using the config specified, on a batch system and tags results

import ROOT
import itertools as it
import os
import numpy as np
import sys
import datetime
import tempfile
import subprocess
import optparse
import glob
import sys
import shutil
import inspect

store_bf_spectrum      = False
use_reactor_constraint = True

fit_par                 = "nuosc3f_12fix+nuisance"
fitter                  = "minuit+hesse"
dmsqscheme              = "abs23"
hierarchy               = "normal"
mixparamchoice          = "sinsqtheta"
tolerance               = 0.01

macroname = "${T2KOSCANA}/macros/sensitivity/run_sensitivity_joint.C"

progname  = "root2k"
flags = "-l -q -b"

time_str = datetime.datetime.today().strftime("%y-%m-%d_%H-%M-%S")

scratch = "/tmp"
user = os.environ['LOGNAME']
init_dir = os.getcwd()

#~ full_macroname = os.path.realpath(macroname)
full_macroname = macroname

for possible_outdir in ["OUTFOLDER","PWD"]:
  try:
    output_dir = os.environ[possible_outdir]
    break
  except:
    continue

def main(setname          ,
         hierarchy        ,
         param1           , param2           ,
         min_param1_bin   , max_param1_bin   ,
         min_param2_bin   , max_param2_bin   ,
         seed             , nexp             ,
         config_file      , dataset_id):
  
  # Get list of function argument values, and build a unique tag
  args = [locals()[arg] for arg in inspect.getargvalues(inspect.currentframe())[0]]
  tag = "_".join(map(str,args))
  
  print "Using tag: "+tag
           
  tmp_dir = tempfile.mkdtemp(suffix="_osc_sensitivity_"+`hash("_".join(sys.argv[1:]))`,prefix=user+"_",dir=scratch)
  print "Made directory "+tmp_dir
  
  os.chdir(tmp_dir)
  print "Moved to directory: "+os.getcwd()
  
  data_file = "toymc-0.xml"
  
  macroargs = [tag, data_file, dataset_id, config_file,
               fit_par, fitter,
               int(store_bf_spectrum), int(use_reactor_constraint),
               dmsqscheme, hierarchy, mixparamchoice, tolerance, 
               param1, param2,
               min_param1_bin, max_param1_bin,
               min_param2_bin, max_param2_bin,
              ]
              
  prepped = macro_prep(macroargs)
               
  macrocall = "'"+full_macroname+"("+(",".join(prepped))+")'"
  
  mycmd = " ".join([progname,flags,macrocall])
  
  print mycmd
  
  result = subprocess.call( mycmd , shell=True,)
  
  #Move output files.
  output_files = glob.glob("*")
  print "ROOT finished, moving files"
  for f in output_files:
    shutil.move(f,output_dir)
    
  # Cleanup
  os.rmdir(tmp_dir)

# Quotes strings in arguments, turns numerals into strings
def macro_prep(args):
  formatted_args = []
  for arg in args:
    s = str(arg)
    if s==arg:
      s='"'+arg+'"'
    formatted_args.append(s)
  return formatted_args 

if __name__=="__main__":
  
  parser = optparse.OptionParser(usage="%prog [options]..",description="Run many VALOR sensitivity jobs on writing to local /tmp")
  parser.add_option('--p1', dest='param1'   , default='mixing23', help='First parameter')
  parser.add_option('--p2', dest='param2'   , default='none'    , help='Second parameter (or "none")')
  parser.add_option('-n'  , dest='nexp'     , default=1         , help='Number of toy experiments to run', type = 'int')
  parser.add_option('-r'  , dest='seed0'    , default=0         , help='Seed for RNG. 0 is truly random' , type = 'int')
  parser.add_option('-t'  , dest='tag'      , default="sens"    , help='Output filename tag')
  parser.add_option('-m'  , dest="hierarchy", default="normal"  , help='Mass hierarchy: "normal" or "inverted"')
  
  parser.add_option('-f', dest='fitpar' ,default="nuosc3f_12fix+nuisance",help='Fit parameter list')
  parser.add_option('-c', '--config', dest='config', default="analysis_configuration.xml", help="Configuration file")
  
  
  parser.add_option('--min-param1-bin', dest='min_param1_bin' , default=-1, help="Bin range for variable, -1 sets full range" , type='int')
  parser.add_option('--max-param1-bin', dest='max_param1_bin' , default=-1, help="Bin range for variable, -1 sets full range" , type='int')
  parser.add_option('--min-param2-bin', dest='min_param2_bin' , default=-1, help="Bin range for variable, -1 sets full range" , type='int')
  parser.add_option('--max-param2-bin', dest='max_param2_bin' , default=-1, help="Bin range for variable, -1 sets full range" , type='int')
  
  parser.add_option('-d', '--dataset', dest='dataset', default=0, help="Which dataset in the XML file to use", type='int')
  
  options,other = parser.parse_args()
  
  param1 = options.param1
  param2 = options.param2
  setname = options.tag
  seed = options.seed0
  nexp = options.nexp
  min_param1_bin   = options.min_param1_bin
  max_param1_bin   = options.max_param1_bin
  min_param2_bin   = options.min_param2_bin
  max_param2_bin   = options.max_param2_bin
  fitpar    = options.fitpar
  hierarchy = options.hierarchy
  config    = options.config
  dataset   = options.dataset
  
  print options
  
  main(setname          , hierarchy        ,
       param1           , param2           ,
       min_param1_bin   , max_param1_bin   ,
       min_param2_bin   , max_param2_bin   ,
       seed             , nexp,
       config           , dataset)
       
