#!/usr/bin/env python

import os
import ROOT as R
import array

R.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C");

valor = R.valor
util = R.valor.util
conv = R.valor.conv

base = os.environ["VALOR_EXPT"]

files_fhc = (
 "fhc",
 base+"data/dune/ndfit/2016a/lartpc/lartpc_nu_minidst_pt2.root",
 base+"data/dune/ndfit/2016a/flux/flux_histos/histos_g4lbne_v3r4p2_QGSP_BERT_Nominal_80GeV_neutrino_LBNEND.root",
 base+"data/dune/ndfit/2016a/flux/flux_histos/histos_g4lbne_v3r4p2_QGSP_BERT_Nominal_80GeV_neutrino_LBNEFD.root"
)

files_rhc = (
 "rhc",
 base+"data/dune/ndfit/2016a/lartpc/lartpc_nubar_minidst_pt2.root",
 base+"data/dune/ndfit/2016a/flux/flux_histos/histos_g4lbne_v3r4p2_QGSP_BERT_Nominal_80GeV_antineutrino_LBNEND.root",
 base+"data/dune/ndfit/2016a/flux/flux_histos/histos_g4lbne_v3r4p2_QGSP_BERT_Nominal_80GeV_antineutrino_LBNEFD.root"
)

def convert_selection(selection, bar):
  if (not bar) :
    if util.IsMuLike(selection):
      if util.IsCCQELike(selection):
        return conv.ksNuMuCChigh
      else:
        return conv.ksNuMuCClow
        
    if util.IsELike(selection):
      if util.IsCCQELike(selection):
        return conv.ksNuECChigh
      else:
        return conv.ksNuECClow
  else:
    if util.IsMuLike(selection):
      if util.IsCCQELike(selection):
        return conv.ksNuMuBarCChigh
      else:
        return conv.ksNuMuBarCClow
        
    if util.IsELike(selection):
      if util.IsCCQELike(selection):
        return conv.ksNuEBarCChigh
      else:
        return conv.ksNuEBarCClow
  
  return conv.ksUndefined
  
def do_reweight(current,input_minidst_file, input_nd_hist_file, input_fd_hist_file):
  
  tf_dst = R.TFile(input_minidst_file)
  assert(tf_dst)
  tree = tf_dst.Get("valor_dune_nd_mini_dst")
  assert(tree)
  
  tf_nd_hist = R.TFile(input_nd_hist_file)
  tf_fd_hist = R.TFile(input_fd_hist_file)
  
  nd_hists = read_all_objects_of_class(tf_nd_hist,"TH1D")
  fd_hists = read_all_objects_of_class(tf_fd_hist,"TH1D")
  
  tokens = input_minidst_file.split(".root")
  output_file = tokens[0]+"_fd.root"
  tf_output = R.TFile(output_file,"RECREATE")
  
  tree.SetBranchStatus("wght",0)
  out_tree = tree.CloneTree(0)
  tree.SetBranchStatus("wght",1)
  
  out_weight = array.array("d",[1.])
  out_mode   = array.array("i",[-1])
  out_sample = array.array("i",[-1])
  
  out_tree.Branch("wght",out_weight,"wght/D")
  out_tree.SetBranchAddress("reac_mode",out_mode)
  out_tree.SetBranchAddress("sample",out_sample)
  
  for iev,event in enumerate(tree):
    try:
      old_weight = event.wght
    except AttributeError:
      old_weight = 1.
      
    if (iev%10000 == 0 and iev):
      print iev
    
    final_flav_name = util.FlavourAsString(util.ReacModeToFinalFlavor(event.reac_mode)).lower()
    
    # Fill the version without oscillation
    nd_hist = nd_hists[final_flav_name+"_flux_forplots"]
    nd_bin  = nd_hist.FindBin(event.enu_true)
    nd_flux = nd_hist.GetBinContent(nd_bin)
    
    fd_hist = fd_hists[final_flav_name+"_flux_forplots"]
    fd_bin  = fd_hist.FindBin(event.enu_true)
    fd_flux = fd_hist.GetBinContent(fd_bin)
    
    flux_reweight = fd_flux/nd_flux
    out_weight[0] = old_weight * flux_reweight
    out_mode  [0] = event.reac_mode
    out_sample[0] = convert_selection(event.sample, current=="rhc")
    out_tree.Fill()
    
    # Fill the version with oscillation, if needed
    osc_mode = util.NuEModeToOscNuEMode(event.reac_mode)
    if osc_mode >=0 and osc_mode < 1E6:
      init_flav_name = util.FlavourAsString(util.ReacModeToInitFlavor(osc_mode)).lower()
      
      fd_hist = fd_hists[init_flav_name+"_flux_forplots"]
      fd_bin  = fd_hist.FindBin(event.enu_true)
      fd_flux = fd_hist.GetBinContent(fd_bin)
      
      flux_reweight = fd_flux/nd_flux
      out_weight[0] = old_weight * flux_reweight
      out_mode[0] = osc_mode
      out_tree.Fill()
      
  out_tree.Write()
  print "Written! Input tree has "+`iev+1`+" events. Output tree has "+`out_tree.GetEntries()`+" events."
  tf_output.Close()
    
  return
  
def read_all_objects_of_class(tfile,oftype="TH1D"):
  objects = dict()
  for k in tfile.GetListOfKeys():
    if k.GetClassName()==oftype or oftype is None:
      objects[k.GetName()] = k.ReadObj()
  return objects
  
def main():
  do_reweight(*files_fhc)
  do_reweight(*files_rhc)

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Use flux histograms to reweight near detector MC into far detector MC")
  
  #~ parser.add_option('-c', '--config', dest='config', default="fit/dune_ndfit.xml", help="Configuration file")
  #~ parser.add_option('-o', '--outfile', dest='outfile', default="initial_errors.root", help="Output ROOT file")
  
  o,args = parser.parse_args()
  
  main()
