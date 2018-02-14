#!/usr/bin/env python

#
# Split the minituple, mtuple from the sk group into useful format for VALOR group
#

import ROOT
import math
import sys
import os
import glob
import numpy as np

samples = dict() # Map of selection name and value of iclass variable
samples["numu"] = 13
samples["nue" ] = 11
samples["cc1pi"] = 14 
# samples["2pi0"] = 2


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
#orig_treename = "h1"
#orig_treename = "nom"

branches_we_care_about = ["iclass",
                          "erec",
                          "pnu",
                          "mode",
                          "wgtflx",
                          "wgtRPA_nr",
                          "wgtRPA_r",
                          "wgtBeRPA_nom",
                          "wgtsfrfg",
                          "wght_CCcoh",
                          "wght_NCgm",
                          "amomm",
                          "amome",
                          "dir",
                          "dirnu",
                          "fqmomm",
                          "fqmome",
                          "fq1rdir"]
#_______________________________________________________________________

def main(input_dir,simplify_tree):
  
  files = glob.glob(input_dir+"/*root")

  f = files[0].rsplit("/")[-1]
  print files
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
  
  output_head = f.split("fine")[0]

  print "Output_file head"
  print output_head
  
  
  output_files = dict()
  for s in samples:
    output_files[s] = ROOT.TFile(output_head+s+"selec.root" ,"RECREATE")
  
  print "Created output files: "
  for f in output_files.itervalues():
    print f.GetName()
  print
  
  if (simplify_tree):
    print "Simplifying tree down to ",len(branches_we_care_about)," branches:"
    print ", ".join(branches_we_care_about)
  
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
    
    if (simplify_tree):
      in_tree.SetBranchStatus("*",0)
      for branch in branches_we_care_about:
        in_tree.SetBranchStatus(branch,1)
    
    n_entries = in_tree.GetEntries()

    output_trees = dict()
    for s,f in output_files.iteritems():
      f.cd()
      output_trees[s] = in_tree.CloneTree(0)
      friends = output_trees[s].GetListOfFriends()
      if friends:
        friends.RemoveAll()
      
    #~ enu = np.array([1.])
    #~ for ot in output_trees:
      #~ ot.Branch("enu"      , enu ,  "enu/D"     )
    
    for event in in_tree:
      #~ enu[0] = event.pnu[0]
      for s,selected_iclass in samples.iteritems():
        if event.iclass == selected_iclass:
          output_trees[s].Fill()

     
    out_tree_name = "t"+ourtype
    out_hist_name = "h_stat_snt_"+out_tree_name
    
    for t in output_trees.itervalues():
      t.SetName(out_tree_name+"_"+current)
    
    scale_events,scale_per_mc = this_scale[typename]
        
    norm_hist = ROOT.TH1D(out_hist_name,out_hist_name,10,0,10)
    norm_hist.Fill(0.5, scale_events )
    norm_hist.Fill(1.5, scale_per_mc * 1E21)
    
    for s in samples.iterkeys():
      output_files[s].cd()
      output_trees[s].Write(out_tree_name)
      norm_hist.Write(out_hist_name)
    
    print
    
  
if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Reformat and optionally reduce input SK MC trees")
  
  #~parser.add_option('-f', '--file', dest='output_file', default="", help="Name for output.")
  parser.add_option ('-s',"--small",dest="make_small",default=False,action="store_true",
                     help="Save only the branches we really need")

  o,other = parser.parse_args()
  
  assert(len(other)) == 1
  main(other[0],o.make_small)
  
