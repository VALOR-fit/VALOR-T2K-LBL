#!/usr/bin/env python

# Script to make table of variations in systematics (tab7 in tn155)
# Generates only body of table
# Adding hlines is left as an exercise for the user
# Takes a while to run (maybe 12 hours)

#Steve Dennis, October 2013

import ROOT
import math
import os
import subprocess
import glob

# Have to finish with osc turned on or the osc params won't be written
columns = ["prefit_no_osc","postfit_no_osc","prefit_with_osc","postfit_with_osc"]

prefit_config  = "2013a_numu_run1+2+3+4_flux11bv3.2_neut13a_banff13v5prefit_sk+fsi.xml"
postfit_config = "2013a_numu_run1+2+3+4_flux11bv3.2_neut13a_banff13v5_sk+fsi.xml"

sample = "1Rmu"
hier   = "normal"
param_choice = "sinsqtheta"

program = "root2k"
prog_args = ["-q","-b"]
ana_dir = os.environ["T2KOSCANA"]

other_modes = [ (200 , "xsec_pf",              ),
                (201 , "xsec_sf",              ),
                (202 , "xsec_wshape",          ),
                (203 , "xsec_shape_ccoth",     ),
                (204 , "xsec_pilessdelta",     ),
                (215 , "xsec_bindE",           ),
                (209 , "xsec_norm_cccoh",      ),
                (210 , "xsec_norm_nc1piC",     ),
                (212 , "xsec_norm_ncoth",      ),
                (213 , "xsec_norm_nuetonumu",  ),
                (214 , "xsec_norm_antinutonu", ),
                (700 , "mixing13",             ),
                (701 , "mixing12",             ),
                (702 , "dm12sq",               ),
                (703 , "deltacp",              ),
              ]

col_width = 15
              
# Load T2K Oscillation analysis libraries
# Load valor libraries if necessary
try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x "+ana_dir+"/root/loadincs.C");    
  ROOT.gROOT.ProcessLine(".x "+ana_dir+"/root/loadlibs.C");
  #~ ROOT.gROOT.ProcessLine(".x "+ana_dir+"/root/root2k.C");    

ROOT.gROOT.SetBatch(1)

conv = ROOT.valor.conv
util = ROOT.valor.util
Opt = ROOT.valor.Opt.Instance()
Rnd = ROOT.valor.Rnd.Instance()

Opt.gOscMixingParamChoice = conv.kNuMixSinsqTheta

expected_oscillation = [-99,-99,0.5,-99,-99, 2.4E-3]

def main():
  for column in columns:
    normalisation_errors(column)
    other_errors(column)
    
  make_table()
  print "Done!"

# Run the macro get normalisation errors
def normalisation_errors(column):
  
  if "no_osc"      in column : osc_params = [0]*6    # No osc
  elif "with_osc"  in column : osc_params = expected_oscillation
  else: raise T2KAnaError("Could not determine whether to use oscillations")
  
  if "prefit"    in column : config = prefit_config
  elif "postfit" in column : config = postfit_config
  else: raise T2KAnaError("Could not determine whether to use prefit or postfit")
  
  macro = ana_dir+"/macros/spectra/test_systematic_error_propagation_all.C"
  macroargs = macro_prep([sample]+osc_params+[config])
  macrocall = "'"+macro+"("+",".join(macroargs)+")'"
  
  fullcall = " ".join([program]+prog_args+[macrocall])
  
  mydir = "syst_"+column
  if not os.path.exists(mydir):
    os.makedirs(mydir)
        
  os.chdir(mydir)
  
  result = subprocess.call( fullcall , shell=True,)

  os.chdir("..")

# Run the macro to get osc/nuxsec errors
def other_errors(column):
  if "no_osc"  in column : osc_params = [0]*6    # No osc
  elif "with_osc"  in column : osc_params = expected_oscillation
  else: raise T2KAnaError("Could not determine whether to use oscillations")
  
  if "prefit"    in column : config = prefit_config
  elif "postfit" in column : config = postfit_config
  else: raise T2KAnaError("Could not determine whether to use prefit or postfit")
  
  macro = ana_dir+"/macros/spectra/test_sk_ereco_spectrum_twk_sys_throw.C"
  
  mydir = "syst_"+column
  if not os.path.exists(mydir):
    os.makedirs(mydir)
        
  os.chdir(mydir)
  
  for mode,name in other_modes: #unconstrained xsec and then osc params
    nexp=20000
    macroargs = macro_prep([mode,1,"3f_osc"]+osc_params+[0,config,nexp,param_choice,hier])
    
    macrocall = "'"+macro+"("+",".join(macroargs)+")'"
    fullcall = " ".join([program]+prog_args+[macrocall])
    
    result = subprocess.call( fullcall , shell=True,)

  os.chdir("..")

# Read the output of the other two macros and make a .tex file
def make_table():
  col_dicts = []
  for column in columns:
    # Reads the file representing the column into a dict
    # with parameter name (in latex) as key and formatted number as value
    # Remembers list of rows from last column file read
    d = dict()
    col_dicts.append(d)
    rows = []
    
    # Read normalisations
    summaries = glob.glob("syst_"+column+"/*.tex.summary")
    if not len(summaries)==1:
      raise T2KAnaError("Found "+`len(summaries)`+ "tex.summary files for "+column+", expected 1")
    f = open(summaries[0])
    
    for l in f:
      k,junk,v = l.partition("&")
      k = k.strip()
      v = v.strip()
      rows.append(k)
      d[k] = format_to2sf_decimal(v)
    
    f.close()
    
    # Read others errors
    # If present, these override the values from the other macro
    thrown = glob.glob("syst_"+column+"/*.txt")
    if not len(thrown)==1:
      raise T2KAnaError("Found "+`len(thrown)`+ "txt files for "+column+", expected 1")
        
    f = open(thrown[0])
    for l in f:
      name,junk,val = l.split("&")
      name = name.strip()
      name = name.replace("prefit","") # the prefit params have 'prefit' appended, which we don't need
      val = val.strip().partition(" ")[0] # Get just the number
      
      physparm_id = util.PhysParamIdFromString(name)
      key = util.AsStringPaperTex(physparm_id)
      if "no_osc" in column and not util.IsSystematic(physparm_id):
        val = "0"
      else:
        val = float(val)*100 # The values from this macro are percentages
        val = format_to2sf_decimal(val)
      #rows.append(key) # All these rows should exist from the other macro.
      d[key] = val
    f.close()
  
  # Actually start writing it out
  f = open("effect_of_systematic_variations.tex","w")
  for param in rows:
    cols = [param.rjust(2*col_width)]
    for column,d in zip(columns,col_dicts):
      cols.append(d[param].rjust(col_width))
    f.write(" & ".join(cols))
    f.write("  \\\\   \n")
  
  f.close()

# Rounds to 2sf, and then writes out as a decimal with up to 5dp
# This is ugly.
def format_to2sf_decimal(value):
  to2sf = float("{0:.1e}".format(float(value)))
  form = "{0:.5f}".format(to2sf).rstrip("0").rstrip(".")
  if "." not in form and form!="0": form+=".0"
  return form
  
# Quotes strings, turns other variables into strings
# Used for preparing arguments to be used as ROOT macro calls
def macro_prep(args):
  formatted_args = []
  for arg in args:
    s = str(arg)
    if s==arg:
      s='"'+arg+'"'
    formatted_args.append(s)
  return formatted_args 

class T2KAnaError(Exception):
  pass

if __name__=="__main__":
  main()
