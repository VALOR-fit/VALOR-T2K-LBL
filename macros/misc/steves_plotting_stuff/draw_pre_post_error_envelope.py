#!/usr/bin/python

# Script to make plots of relative error envelope sizes
# eg plots 9-12 in TN-183
# Takes output of macros/spectrum/test_sk_ereco_spectrum_twk_sys_throw_runAll.py

#Steve Dennis, October 2013

import ROOT
import math
import os
import glob

linewidth = 3
title_offset = 0.8
bottom_margin = 0.15
use_log_plots = False

default_oscillated_input_file = (
  "sk_ereco_spectrum_RHC1Rmu-twk_sys_throw_sinsqtheta[12_0.304,13_0.042,23_0.524],dcp_1.910,dmsq_Normal[12_0.00007530,23_0.002510]eV2.root")
default_unoscillated_input_file = (
  "sk_ereco_spectrum_RHC1Rmu-twk_sys_throw_sinsqtheta[12_0.000,13_0.000,23_0.000],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002510]eV2.root")

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kMagenta+1,  ROOT.kGreen+2, ROOT.kOrange+7, ROOT.kCyan+2,]

param_sets = ["SKDet_FSI"   ,
              "BANFF_xsec"  ,
              "BANFF_flux"  ,
              "all"         ,]


readable = dict()
readable["SKDet_FSI"         ] = "SK Detector + FSI"
readable["BANFF_xsec"        ] = "Cross-section"
readable["BANFF_flux"        ] = "Flux"
readable["all"               ] = "All systematics"

plot_opt = ""

try:
  ROOT.t2kana
except:
  #~ ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/loadincs.C")
  #~ ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/loadlibs.C")
  ROOT.gROOT.ProcessLine(".L $T2KOSCANA/root/root2k.C");
  pass

  

def main(sample,unoscillated_input_filename, oscillated_input_filename):
  f_osc = ROOT.TFile(oscillated_input_filename)
  f_unosc = ROOT.TFile(unoscillated_input_filename)
  if not f_osc and f_unosc:
    raise T2KAnaError("Could not open input files")
    
  pre_and_postfit_envelope(sample,f_unosc,"unoscillated")
  pre_and_postfit_envelope(sample,f_osc  ,"oscillated" )
  

def pre_and_postfit_envelope(sample,rootfile,name):
  
  nominal_hist = ROOT.TH1D(rootfile.Get(sample+"_nominal"))
  nominal_hist.SetTitle("")
  nominal_hist.GetXaxis().SetRangeUser(0.,7.)
  nominal_hist.GetXaxis().SetRangeUser(0.,7.)
  nominal_hist.GetXaxis().SetTitle("E_{reco}")
  nominal_hist.GetYaxis().SetTitle("Events per bin")
  nominal_hist.GetYaxis().SetTitleOffset(title_offset)
  nominal_hist.SetLineColor(ROOT.kBlack)
  nominal_hist.SetLineStyle(1)
  nominal_hist.SetLineWidth(linewidth)
  
  all_hists = [nominal_hist]
  
  for ii,par_set in enumerate(param_sets):
    
    par_hists = [nominal_hist]
    
    hname = sample+"_"+par_set+"prefit"+"_envpluserr"
    h = rootfile.Get(hname)
    if not h:
      print "Failed to get histogram: ",hname ; exit()
    hplus_prefit     = ROOT.TH1D( h )
    
    hname = sample+"_"+par_set+"prefit"+"_envminuserr"
    h = rootfile.Get(hname)
    if not h:
      print "Failed to get histogram: ",hname ; exit()
    hminus_prefit     = ROOT.TH1D( h )
    
    hname = sample+"_"+par_set+"_envpluserr"
    h = rootfile.Get(hname)
    if not h:
      print "Failed to get histogram: ",hname ; exit()
    hplus_postfit     = ROOT.TH1D( h )
    
    hname = sample+"_"+par_set+"_envminuserr"
    h = rootfile.Get(hname)
    if not h:
      print "Failed to get histogram: ",hname ; exit()
    hminus_postfit     = ROOT.TH1D( h )
    
      
    hists = [hplus_prefit,hminus_prefit,hplus_postfit,hminus_postfit]
    for h in hists:
      par_hists.append(h)
      all_hists.append(h)
    
    for h in [hplus_prefit,hminus_prefit]:
      h.SetLineColor(ROOT.kRed)
      h.SetLineStyle(2)
      h.SetLineWidth(linewidth)
    for h in [hplus_postfit,hminus_postfit]:
      h.SetLineColor(ROOT.kBlue)
      h.SetLineStyle(3)
      h.SetLineWidth(linewidth)
      
    legend = ROOT.TLegend(0.6,0.6,0.9,0.9);
    legend.SetFillStyle(0);
    legend.SetBorderSize(0);
    legend.SetTextSize(0.04);
      
    legend.AddEntry(hplus_prefit ,"Prefit" ,"L")
    legend.AddEntry(hplus_postfit,"Postfit","L")
  
    limit = max([abs(h.GetMinimum()) for h in par_hists] + [abs(h.GetMaximum()) for h in par_hists])
    
    canv = ROOT.TCanvas()
    canv.SetBottomMargin(bottom_margin)
    if use_log_plots: canv.SetLogx()
    
    for h in par_hists:
      h.SetMinimum(0.)
      h.SetMaximum(1.05*limit)
      h.Draw("same"+plot_opt)
    
    legend.Draw()
    canv.Print(name+"_"+par_set+"_envelope.eps")
    
  return

class T2KAnaError(Exception):
  pass

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make plots of size of systematic/statistical effects")
  parser.add_option('-s', dest='sample',default="RHC1Rmu"  ,help='The sample')
  parser.add_option('-o', dest='osc'  ,default=default_oscillated_input_file  ,help='The unoscillated input ROOT file')
  parser.add_option('-u', dest='unosc',default=default_unoscillated_input_file,help='The oscillated input ROOT file')
  
  options,other = parser.parse_args()
  
  main(options.sample,options.unosc,options.osc)
