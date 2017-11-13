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


pattern = "WithFSI_spline2017[VERSION]_t2ksk.14a.neut5.3.2.13a_tuned_v1r1.[CURRENT].fine.[TYPE]_[SELECTION].root" # The file naming convention
#pattern = "[TYPE]_[SELECTION]_[CURRENT].root" # The file naming convention

location = os.environ[env_var]+"/inputs/official/2018/xsec/splines2018fitqunv9_withFSI"

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
               ("erpa_E_ccqe"     , False ),
               ("NC1gam_EtrueEreco_FEFQE"   , False ),
               ("2p2h_EtrueEreco_FEFQE"     , False ),
               ("NCoth_EtrueEreco_FEFQE"    , False ),
               ("NCcoh_EtrueEreco_FEFQE"    , False ),
               ("NC1pipm_EtrueEreco_FEFQE"  , False ),
               ("NC1pi0_EtrueEreco_FEFQE"   , False ),
               ("CCDIS_EtrueEreco_FEFQE"    , False ),
               ("CCcoh_EtrueEreco_FEFQE"    , False ),
               ("CC1pi_EtrueEreco_FEFQE"    , False ),
               ("CCQE_EtrueEreco_FEFQE"     , False ),
               ("NC1gam_EtrueEreco_FEFCX"   , False ),
               ("2p2h_EtrueEreco_FEFCX"     , False ),
               ("NCoth_EtrueEreco_FEFCX"    , False ),
               ("NCcoh_EtrueEreco_FEFCX"    , False ),
               ("NC1pipm_EtrueEreco_FEFCX"  , False ),
               ("NC1pi0_EtrueEreco_FEFCX"   , False ),
               ("CCDIS_EtrueEreco_FEFCX"    , False ),
               ("CCcoh_EtrueEreco_FEFCX"    , False ),
               ("CC1pi_EtrueEreco_FEFCX"    , False ),
               ("CCQE_EtrueEreco_FEFCX"     , False ),
               ("NC1gam_EtrueEreco_FEFQEH"  , False ),
               ("2p2h_EtrueEreco_FEFQEH"    , False ),
               ("NCoth_EtrueEreco_FEFQEH"   , False ),
               ("NCcoh_EtrueEreco_FEFQEH"   , False ),
               ("NC1pipm_EtrueEreco_FEFQEH" , False ),
               ("NC1pi0_EtrueEreco_FEFQEH"  , False ),
               ("CCDIS_EtrueEreco_FEFQEH"   , False ),
               ("CCcoh_EtrueEreco_FEFQEH"   , False ),
               ("CC1pi_EtrueEreco_FEFQEH"   , False ),
               ("CCQE_EtrueEreco_FEFQEH"    , False ),
               ("NC1gam_EtrueEreco_FEFABS"  , False ),
               ("2p2h_EtrueEreco_FEFABS"    , False ),
               ("NCoth_EtrueEreco_FEFABS"   , False ),
               ("NCcoh_EtrueEreco_FEFABS"   , False ),
               ("NC1pipm_EtrueEreco_FEFABS" , False ),
               ("NC1pi0_EtrueEreco_FEFABS"  , False ),
               ("CCDIS_EtrueEreco_FEFABS"   , False ),
               ("CCcoh_EtrueEreco_FEFABS"   , False ),
               ("CC1pi_EtrueEreco_FEFABS"   , False ),
               ("CCQE_EtrueEreco_FEFABS"    , False ),
               ("NC1gam_EtrueEreco_FEFINEL" , False ),
               ("2p2h_EtrueEreco_FEFINEL"   , False ),
               ("NCoth_EtrueEreco_FEFINEL"  , False ),
               ("NCcoh_EtrueEreco_FEFINEL"  , False ),
               ("NC1pipm_EtrueEreco_FEFINEL", False ),
               ("NC1pi0_EtrueEreco_FEFINEL" , False ),
               ("CCDIS_EtrueEreco_FEFINEL"  , False ),
               ("CCcoh_EtrueEreco_FEFINEL"  , False ),
               ("CC1pi_EtrueEreco_FEFINEL"  , False ),
               ("CCQE_EtrueEreco_FEFINEL"   , False )
#              ("erpa_U_ccqe"     , False )               
]
               
#samples   = ["1Rmu","1Re","RHC1Rmu","RHC1Re","MultiRe"]
samples   = ["1Re","RHC1Re","MultiRe"]
#samples = ["MultiRe"]


selections = dict()
selections["1Rmu"]    = "numuselec" # Our samples to Kirsty's selection names
selections["RHC1Rmu"] = "numuselec"
selections["1Re"]     = "nueselec"
selections["RHC1Re" ] = "nueselec"
selections["MultiRe"] = "cc1piselec"

currents = dict()
currents["1Rmu"]    = "250ka"      # Our samples to current modes
currents["1Re"]     = "250ka"
currents["RHC1Rmu"] = "m250ka"
currents["RHC1Re" ] = "m250ka"
currents["MultiRe"] = "250ka"

ereco_bin_offsets = dict()
ereco_bin_offsets["1Rmu"]     = 0
ereco_bin_offsets["RHC1Rmu"]  = 0
ereco_bin_offsets["1Re"]      = 0 #-1 # bin 0 of the 1Re samples from Kirsty are negative energy, remove it.
ereco_bin_offsets["RHC1Re" ]  = 0 #-1 # bin 0 of the 1Re samples from Kirsty are negative energy, remove it.
ereco_bin_offsets["MultiRe" ] = 0 #-1 # bin 0 of the 1Re samples from Kirsty are negative energy, remove it.

# Mapping between VALOR modes and input modes
replacements = [("NC1pipm","nc1piC"),
                ("NC1gam" ,"nc1gamma"),
                ("2p2h" ,"mec"),
                ("NCoth" ,"ncoth"),
                ("NCcoh" ,"nccoh"),
                ("NC1" ,"nc1"),
                ("CCDIS" ,"ccoth"),
                ("CCcoh" ,"cccoh"),
                ("CC1pi" ,"cc1pi"),
                ("CCQE" ,"ccqe"),
                ("pipm","1piC"),
                ("piz" ,"1pi0")]

nknots_sigma = 13
sigmas = np.array([0.5*i - 3. for i in xrange(nknots_sigma)])
#______________________________________________________________________________________________________________________________

def main(version, loc):
  location = os.environ[env_var]+"/inputs/official/2018/xsec/" + loc
  for sample in samples:
    current   = currents  [sample]
    selection = selections[sample]
    
    output_file  = ROOT.TFile(version+"_"+sample+".root","RECREATE") # Output for this sample's splines
    
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
        elif "FEF" in syst:
          mode, _, par = syst.split("_")
        else:  
          par,mode = syst.split("_")

        #print par, mode  
        
        for (old,new) in replacements:
          mode = mode.replace(old,new)
        
        # Get all the splines for this systematic+mode
        spline_dict = dict()
        for key in tf.GetListOfKeys():
          name = key.GetName()
          if not (syst in name and ("_sp_" in name or ("FEF" in syst and "sp_" in name))): continue
          spline_dict[name] = key.ReadObj()
      
        output_file.cd()
        
        # Loop over the splines - rename them and write them to the file.
        for orig_spl_name, spline in spline_dict.iteritems():
          true_bin,reco_bin,q2bin = map(int, orig_spl_name.split("_")[-3:])
          
          reco_bin += ereco_bin_offsets[sample]
          if reco_bin < 0: continue
          
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
  
  parser.add_option('-v', '--version', dest='version', default="fitqunv9", help="Current spline version.")
  parser.add_option('-l', '--location', dest='location', default="fitqunv9", help="Input spline location.")
  
  o,other = parser.parse_args()
  
  main(o.version, o.location)
  
  if other:
    print "Warning - ignored arguments: ",other
