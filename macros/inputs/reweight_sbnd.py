#!/usr/bin/env python

import os
import ROOT as R
import array

R.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

valor = R.valor
util = R.valor.util
conv = R.valor.conv

base = os.environ["VALOR_EXPT"]

def wrong_sign_sample_to_right_sign(s):
  if (s == conv.ksNuMuCCWS)     : return conv.ksNuMuCCoth
  if (s == conv.ksNuMuCCWS0pi ) : return conv.ksNuMuCCWS0pi
  if (s == conv.ksNuMuCCWS1piC) : return conv.ksNuMuCCWS1piC
  if (s == conv.ksNuMuCCWS1pi0) : return conv.ksNuMuCCWS1pi0
  if (s == conv.ksNuECCWS     ) : return conv.ksNuECCWS
  
  return s
  
def do_reweight(input_minidst_file, ratio_hists_file, output_file):
  
  tf_dst = R.TFile(input_minidst_file)
  assert(tf_dst)
  tree = tf_dst.Get("valor_dune_nd_mini_dst")
  assert(tree)
  
  tf_ratio_hist = R.TFile(ratio_hists_file)
  
  ratio_hists = read_all_objects_of_class(tf_ratio_hist,"TH1D")
  
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
    ratio_hist = ratio_hists[final_flav_name+"_110m_ratio"]
    ratio_bin  = ratio_hist.FindBin(event.enu_true)
    ratio_flux = ratio_hist.GetBinContent(ratio_bin)
    
    out_weight[0] = old_weight * ratio_flux
    out_mode  [0] = event.reac_mode
    out_sample[0] = wrong_sign_sample_to_right_sign(event.sample)
    out_tree.Fill()
    
    #~ # Fill the version with oscillation, if needed
    #~ osc_mode = util.NuEModeToOscNuEMode(event.reac_mode)
    #~ if osc_mode >=0 and osc_mode < 1E6:
      #~ init_flav_name = util.FlavourAsString(util.ReacModeToInitFlavor(osc_mode)).lower()
      
      #~ fd_hist = fd_hists[init_flav_name+"_flux_forplots"]
      #~ fd_bin  = fd_hist.FindBin(event.enu_true)
      #~ fd_flux = fd_hist.GetBinContent(fd_bin)
      
      #~ flux_reweight = fd_flux/nd_flux
      #~ out_weight[0] = old_weight * ratio_flux
      #~ out_mode[0] = osc_mode
      #~ out_tree.Fill()
      
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
  do_reweight(
    os.environ['VALOR_EXPT']+'data/dune/ndfit/2017a/fgt_nu_minidst_pt4.root',
    os.environ['VALOR_EXPT']+'data/sbn/2017a/flux/sbn_flux_ratio_hist_fhc.root',
    os.environ['VALOR_EXPT']+'data/sbn/2017a/sbnd_fhc_from_dune_fgt.root'
    )

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Use flux histograms to reweight near detector MC into far detector MC")
  
  #~ parser.add_option('-c', '--config', dest='config', default="fit/dune_ndfit.xml", help="Configuration file")
  #~ parser.add_option('-o', '--outfile', dest='outfile', default="initial_errors.root", help="Output ROOT file")
  
  o,args = parser.parse_args()
  
  main()
