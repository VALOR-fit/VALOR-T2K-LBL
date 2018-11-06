#!/usr/bin/env python

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


pattern = "spline2017fitqun[VERSION]erectheta_t2ksk.14a.neut5.3.2.13a_tuned_v1r1.[CURRENT].fine.[TYPE]_[SELECTION].root"
#pattern = "[TYPE]_[SELECTION]_[CURRENT].root" # The file naming convention

#location = os.environ[env_var]+"/inputs/official/2015/xsec/2d_splines/theta/1bin/float" # Location of the spline files
#location = os.environ[env_var]+"/inputs/official/2017/xsec/splinesandweights2017fitqunv1_namesfixed/splines2017fitqunv1/"
location = "/data/t2k/phrdqd/analysis/irods/splines/splines2017fitqunv4/"

types = [("nuebar_x_nuebar"  ,"nuebar"    ),  # First entry is new-style names, second entry is VALOR type names.
         ("nue_x_nue"        ,"nue"       ),
         ("numubar_x_nuebar" ,"oscnuebar" ),
         ("numubar_x_numubar","numubar"   ),
         ("numu_x_nue"       ,"oscnue"    ),
         ("numu_x_numu"      ,"numu"      ),
        ]

# Name/mode as first entry, whether or not to reformat the spline to use sigmas is in second
systematics = [("maqe_ccqe"       , False ),
               ("pfo_ccqe"        , False ),
#               ("ebo_ccqe"        , False ),
               ("ca5_cc1pi"       , False ),
               ("ca5_ncpiz"       , False ),
               ("ca5_ncpipm"      , False ),
               ("manff_cc1pi"     , False ),
               ("manff_ncpiz"     , False ),
               ("manff_ncpipm"    , False ),
               ("bgscl_cc1pi"     , False ),
               ("bgscl_ncpiz"     , False ),
               ("bgscl_ncpipm"    , False ),
               ("dismpishp_ccoth" , False ),
               ("sccvec_ccqe"     , False ),
               ("sccvec_ncoth"    , False ),
               ("sccaxl_ccqe"     , False ),
               ("sccaxl_ncoth"    , False ),
               ("rpa_ccqe"        , False ),
               ("pfsf_ccqe"       , False ),
               ("maqesf_ccqe"     , False ),
               ("mecpdd_mec"      , False ),
               ("erpa_A_ccqe"     , False ),
               ("erpa_B_ccqe"     , False ),
               ("erpa_D_ccqe"     , False ),
               ("erpa_E_ccqe"     , False )
#              ("erpa_U_ccqe"     , False )               
]

samples   = ["1Re","RHC1Re","MultiRe"]
#samples   = ["1Rmu","1Re","RHC1Rmu","RHC1Re","MultiRe"]
#samples = ["MultiRe"]


selections = dict()
selections["1Re"]     = "nueselec"
selections["RHC1Re" ] = "nueselec"
selections["MultiRe"] = "nueselec"

currents = dict()
currents["1Re"]     = "250ka"
currents["RHC1Re" ] = "m250ka"
currents["MultiRe"] = "250ka"

replacements = [("pipm","1piC"),
                ("piz" ,"1pi0") ] # VALOR modes have pi0 instead of piz

nknots_sigma = 13
sigmas = np.array([0.5*i - 3. for i in xrange(nknots_sigma)])
#______________________________________________________________________________________________________________________________

def main(version):
  for sample in samples:
    current   = currents  [sample]
    selection = selections[sample]
    
    output_file  = ROOT.TFile("splines_fitqun"+version+"erectheta_"+sample+".root","RECREATE") # Output for this sample's splines
    
    for orig_type,valor_type in types:
      
      # Find the right input file
      infile = pattern
      indir  = location
      for old,new in [("[CURRENT]"  ,current),
                      ("[TYPE]"     ,orig_type),
                      ("[SELECTION]",selection),
                      ("[VERSION]"  ,version) ]:
        infile   = infile .replace(old,new)
        indir    = indir.replace(old,new)
      
      print "Processing file: ",infile
      path = indir+"/"+infile
      tf = ROOT.TFile(path)
      if not tf: print "Failed to open file:\n",path
      
      # Loop over spline types - systematic+mode
      for syst, change_spline_limits in systematics:
        par  = ""
        mode = ""
        if "erpa" in syst:
          ppar,coeff,mode = syst.split("_")
          par = ppar +"_"+ coeff
        else:  
          par,mode = syst.split("_")

        #print par, mode  
        
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
          
          true_bin,reco_bin,q2bin = map(int,bins.split("_"))
          
          new_spl_name = "spl_"+sample+"_"+par+"_"+valor_type+mode+"_"+`true_bin`+"_"+`reco_bin`+"_"+`q2bin`
          
          nknots = spline.GetNp()
          knots = np.zeros(nknots,dtype=float)
          
          if change_spline_limits:
            if nknots != nknots_sigma:
              print "Warning: trying to attempt to remake spline with wrong number of knots: ",par,nknots
              
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
  
  parser.add_option('-v', '--version', dest='version', default="v3", help="Current spline version.")
  
  o,other = parser.parse_args()
  
  main(o.version)
  
  if other:
    print "Warning - ignored arguments: ",other
