#!/usr/bin/python

import ROOT
import math
import sys
import os
import glob
import numpy as np

types = [("nuebar_x_nuebar"  ,"nuebar"    ),  # First entry is new-style names, second entry is VALOR type names.
         ("nue_x_nue"        ,"nue"       ),
         ("numubar_x_nuebar" ,"signuebar" ),
         ("numubar_x_numubar","numubar"   ),
         ("numu_x_nue"       ,"signue"    ),
         ("numu_x_numu"      ,"numu"      ),
        ]

# Scales for 14a MC, 13a nominal flux from http://www.t2k.org/t2ksk/documentation-page/MCWeighting/MCReweighting14a
scales_FHC = dict()
scales_FHC["numu_x_numu"        ] = ( 1396.245 , 634881. )
scales_FHC["nue_x_nue"          ] = (   26.663 , 126899. )
scales_FHC["numubar_x_numubar"  ] = (   50.964 , 634600. )
scales_FHC["nuebar_x_nuebar"    ] = (    2.884 , 126728. )
scales_FHC["numu_x_nue"         ] = ( 1483.994 , 126916. )
scales_FHC["numubar_x_nuebar"   ] = (   52.939 , 126984. )

scales_RHC = dict()
scales_RHC["numu_x_numu"        ] = ( 208.369 , 634659. )
scales_RHC["nue_x_nue"          ] = (   8.941 , 127262. )
scales_RHC["numubar_x_numubar"  ] = ( 353.789 , 633988. )
scales_RHC["nuebar_x_nuebar"    ] = (   6.811 , 127013. )
scales_RHC["numu_x_nue"         ] = ( 215.099 , 127278. )
scales_RHC["numubar_x_nuebar"   ] = ( 376.420 , 126812. )

orig_treename = "mtuple"
#______________________________________________________________________________________________________________________________

def main(input_dir):
  
  files = glob.glob(input_dir+"/*root")
  
  f = files[0].rsplit("/")[-1]
  
  if not "t2ksk.14a" in f:
    print "Unknown MC scaling"
    exit(-1)
  
  if "m250ka" in f:
    print "Using RHC scaling\n"
    this_scale = scales_RHC
    current = "RHC"
  elif "250ka" in f:
    print "Using FHC scaling\n"
    this_scale = scales_FHC
    current = "FHC"
  else:
    print "Unknown horn current"
  
  output_head = files[0].split("fine")[0]
  
  output_file_numu = ROOT.TFile(output_head+"numuselec.root" ,"RECREATE")
  output_file_nue  = ROOT.TFile(output_head+"nueselec.root"  ,"RECREATE")
  output_file_pi0  = ROOT.TFile(output_head+"2pi0selec.root" ,"RECREATE")
  
  output_files = [output_file_numu, output_file_nue, output_file_pi0]
  #~ output_files = [output_file_numu, output_file_nue]
  
  print "Created output files: "
  for f in output_files:
    print f.GetName()
  print
  
  for typename,ourtype in types:
    
    appropriate_files = [f for f in files if typename in f]
    print appropriate_files
    #~ assert(len(appropriate_files)==1)
    
    ROOT.gROOT.cd()
    in_tree = ROOT.TChain(orig_treename)
    assert(in_tree)
    for fn in appropriate_files:
      in_tree.AddFile(fn)
    
    print "Reading "+ourtype+" from tree "+orig_treename+" in "+`appropriate_files`
    
    n_entries = in_tree.GetEntries()
    
    output_file_numu.cd()
    trmu  = in_tree.CloneTree(0)
    output_file_nue.cd()
    tre   = in_tree.CloneTree(0)
    output_file_pi0.cd()
    trpi0 = in_tree.CloneTree(0)
    
    #~ out_trees = [trmu,tre]
    out_trees = [trmu,tre,trpi0]
    
    #~ enu = np.array([1.])
    #~ for ot in out_trees:
      #~ ot.Branch("enu"      , enu ,  "enu/D"     )
    
    for event in in_tree:
      #~ enu[0] = event.pnu[0]
      if event.iclass == 3:
        trmu.Fill()
      if event.iclass == 2:
        trpi0.Fill()
        #~ tre.Fill()
      if event.iclass == 1:
        tre.Fill()
    
    out_tree_name = "t"+ourtype
    out_hist_name = "h_stat_snt_"+out_tree_name
    
    for t in out_trees:
      t.SetName(out_tree_name+"_"+current)
    
    scale_events,scale_per_mc = this_scale[typename]
        
    norm_hist = ROOT.TH1D(out_hist_name,out_hist_name,10,0,10)
    norm_hist.Fill(0.5, scale_events )
    norm_hist.Fill(1.5, scale_per_mc * 1E21)
    
    output_file_numu.cd()
    trmu.Write(out_tree_name)
    norm_hist.Write(out_hist_name)
    
    output_file_nue.cd()
    tre.Write(out_tree_name)
    norm_hist.Write(out_hist_name)
    
    output_file_pi0.cd()
    trpi0.Write(out_tree_name)
    norm_hist.Write(out_hist_name)
    
    print
    
  
if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Reformat splines")
  
  #~ parser.add_option('-f', '--file', dest='output_file', default="", help="Name for output.")
  
  o,other = parser.parse_args()
  
  assert(len(other)) == 1
  main(other[0])
  
