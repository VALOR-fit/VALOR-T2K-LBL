#!/usr/bin/python

import ROOT
import math
import array
import numpy as np
import sys
import os
import ctypes

env_var = "T2KOSCANA"
try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()
arridx= util.ArrIdx.Instance()


pattern = "spline2015v1_t2ksk.14a.neut5.3.2.13a_tuned_v1r0.[CURRENT].fine.[TYPE]_[SELECTION].root" # The file naming convention

location = os.environ[env_var]+"/inputs/official/2015/xsec/splines_v1/" # Location of the spline files

types = [("nuebar_x_nuebar"  ,"nuebar"    ),  # First entry is new-style names, second entry is VALOR type names.
         ("nue_x_nue"        ,"nue"       ),
         ("numubar_x_nuebar" ,"oscnuebar" ),
         ("numubar_x_numubar","numubar"   ),
         ("numu_x_nue"       ,"oscnue"    ),
         ("numu_x_numu"      ,"numu"      ),
        ]

# Name/mode as first entry, whether or not to reformat the spline to use sigmas is in second
systematics = [("maqe_ccqe"       , True  ),
               ("pfo_ccqe"        , True  ),
               ("ebo_ccqe"        , True  ),
               ("ca5_cc1pi"       , True  ),
               ("ca5_ncpiz"       , True  ),
               ("ca5_ncpipm"      , True  ),
               ("manff_cc1pi"     , True  ),
               ("manff_ncpiz"     , True  ),
               ("manff_ncpipm"    , True  ),
               ("bgscl_cc1pi"     , True  ),
               ("bgscl_ncpiz"     , True  ),
               ("bgscl_ncpipm"    , True  ),
               ("dismpishp_ccoth" , True  ),
               ("sccvec_ccqe"     , True  ),
               ("sccvec_ncoth"    , True  ),
               ("sccaxl_ccqe"     , True  ),
               ("sccaxl_ncoth"    , True  ),
               ("rpa_ccqe"        , False )
]
               
samples   = ["1Rmu","1Re","RHC1Rmu","RHC1Re"]

selections = dict()
selections["1Rmu"]    = "numuselec" # Our samples to Kirsty's selection names
selections["RHC1Rmu"] = "numuselec"
selections["1Re"]     = "nueselec"
selections["RHC1Re" ] = "nueselec"

currents = dict()
currents["1Rmu"]    = "250ka"      # Our samples to current modes
currents["1Re"]     = "250ka"
currents["RHC1Rmu"] = "m250ka"
currents["RHC1Re" ] = "m250ka"


replacements = [("pipm","piC"),
                ("piz" ,"pi0") ] # VALOR modes have pi0 instead of piz

nknots_sigma = 13
sigmas = np.array([0.5*i - 3. for i in xrange(nknots_sigma)])
#______________________________________________________________________________________________________________________________

def main():
  for sample in samples:
    current   = currents  [sample]
    selection = selections[sample]
    
    output_file  = ROOT.TFile("splines_"+sample+".root","RECREATE") # Output for this sample's splines
    
    for orig_type,valor_type in types:
      
      # Find the right input file
      infile = pattern
      for old,new in [("[CURRENT]"  ,current),
                      ("[TYPE]"     ,orig_type),
                      ("[SELECTION]",selection)]:
        infile = infile.replace(old,new)
      
      print "Processing file: ",infile
      path = location+"/"+infile
      tf = ROOT.TFile(path)
      if not tf: print "Failed to open file:\n",path
      
      # Loop over spline types - systematic+mode
      for syst, change_spline_limits in systematics:
        
        par,mode = syst.split("_")
        
        for (old,new) in replacements:
          mode = mode.replace(old,new)
        
        # Get all the splines for this systematic+mode
        spline_dict = dict()
        for key in tf.GetListOfKeys():
          name = key.GetName()
          if not (syst in name and "_sp_" in name): continue
          spline_dict[name] = key.ReadObj()
      
        output_file.cd()
        
        # Loop over the splines - rename them and write them to the file.
        for orig_spl_name, spline in spline_dict.iteritems():
          
          junk,bins = orig_spl_name.rsplit("_sp_")
          
          new_spl_name = "spl_"+sample+"_"+par+"_"+valor_type+mode+"_"+bins
          
          nknots = spline.GetNp()
          knots = np.zeros(nknots,dtype=float)
          
          if change_spline_limits:
            if nknots != nknots_sigma:
              print "Warning: trying to attempt to remake splint with wrong number of knots: ",par,nknots
              
            x = ROOT.Double()
            y = ROOT.Double()
            
            for i in xrange(nknots):
              spline.GetKnot(i,x,y)
              knots[i] = y
          
            new_spline = ROOT.TSpline3(new_spl_name,sigmas,knots,nknots)
          
          else:
            new_spline = ROOT.TSpline3(spline)
            new_spline.SetName(new_spl_name)
            
          new_spline.Write(new_spl_name)
      
      tf.Close()
      # end of input file loop
      
    output_file.Close()
  
if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Reformat splines")
  
  #~ parser.add_option('-f', '--file', dest='output_file', default="", help="Name for output.")
  
  o,other = parser.parse_args()
  
  main()
  
  if other:
    print "Warning - ignored arguments: ",other
