#!/usr/bin/env python

# Script to make plots of relative error envelope sizes
# eg plots 9-12 in TN-183
# Takes output of macros/spectrum/test_sk_ereco_spectrum_twk_sys_throw_runAll.py

#Steve Dennis, October 2013

import ROOT
import math
import os
import glob

ROOT.PyConfig.IgnoreCommandLineOptions = True

linewidth = 1
title_offset = 0.8
bottom_margin = 0.15
use_log_plots = False

default_oscillated_input_file = (
  "tab6/sk_ereco_spectrum-twk_sys_throw_sinsqtheta[12_0.306,13_0.025,23_0.500],dcp_0.000,dmsq_Normal[12_0.00007500,23_0.002400]eV2.root")
default_unoscillated_input_file = (
  "tab6/sk_ereco_spectrum-twk_sys_throw_sinsqtheta[12_0.000,13_0.000,23_0.000],dcp_0.000,dmsq_Normal[12_0.00000000,23_0.000000]eV2.root")

ROOT.gROOT.SetBatch(1)

colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kMagenta+1,  ROOT.kGreen+2, ROOT.kOrange+7, ROOT.kCyan+2,]

param_sets = ["all",
              "SKDet_FSI",
              "BANFF",
              #~ "unconstrained_xsec",
              #~ "Osc_param"
              ]
readable = dict()
readable["SKDet_FSI"         ] = "SK Detector + FSI"
readable["BANFF"             ] = "BANFF"
readable["unconstrained_xsec"] = "Uncorrelated cross-section"
readable["all"               ] = "All systematics"
readable["Osc_param"         ] = "Non-23 sector oscillations"

plot_opt = ""

try:
  ROOT.t2kana
except:
  #~ ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/loadincs.C")
  #~ ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/loadlibs.C")
  ROOT.gROOT.ProcessLine(".L $T2KOSCANA/root/root2k.C");
  pass

  

def main(unoscillated_input_filename, oscillated_input_filename):
  f_osc = ROOT.TFile(oscillated_input_filename)
  f_unosc = ROOT.TFile(unoscillated_input_filename)
  if not f_osc and f_unosc:
    raise T2KAnaError("Could not open input files")
    
  ratio_to_nominal(f_unosc,"unoscillated")
  ratio_to_nominal(f_osc  ,"oscillated" )
  
  #~ ratio_syst_to_stat_err(f_unosc,"unoscillated") # Not currently used
  #~ ratio_syst_to_stat_err(f_osc  ,"oscillated" )  # Not currently used
  
  syst_err_ratio_to_total_syst(f_unosc,"unoscillated")
  syst_err_ratio_to_total_syst(f_osc  ,"oscillated" )
  

def ratio_to_nominal(rootfile,name):
  # Make the plot with ratio of stat/syst to nominal
  canv = ROOT.TCanvas()
  canv.SetBottomMargin(bottom_margin)
  if use_log_plots: canv.SetLogx()
  legend = ROOT.TLegend(0.5,0.8,0.9,0.88);
  legend.SetFillStyle(0);
  legend.SetBorderSize(0);
  legend.SetTextSize(0.03);
  
  all_hists = []
  
  for ii,par_set in enumerate(param_sets):
    if par_set != "all": continue
    try:
      hplus  = ROOT.TH1D(rootfile.Get(par_set+"_pluserr") )
    except TypeError:
      continue  # Osc params won't load for unoscillated
    hists = [hplus] # Superfluous loop as an artifact of when this made plus/minus plots
    for h in hists:
      #~ h.SetTitle("Effects of stat + syst ("+name+")")
      h.SetTitle("")
      h.GetXaxis().SetRangeUser(0.,10.)
      h.GetYaxis().SetTitle("Ratio of events at 1 #sigma to nominal")
      h.GetYaxis().SetTitleOffset(title_offset)
      h.SetLineColor(colors[ii+1])
      h.SetLineWidth(linewidth)
      
      all_hists.append(h)
      
    legend.AddEntry(hplus,readable[par_set],"L")
  
  nominal = rootfile.Get("nominal")
  stat_err = ROOT.TH1D(hplus)
  stat_err.SetLineColor(colors[0])
  for ibin in xrange(1,nominal.GetNbinsX()):
    try: err = 1./math.sqrt(nominal.GetBinContent(ibin))
    except ZeroDivisionError: err=0.
    stat_err.SetBinContent(ibin,err)
    
  stat_err_m = ROOT.TH1D(stat_err)
  stat_err_m.Scale(-1.)
  
  legend.AddEntry(stat_err,"Statistical error","L")
  
  all_hists = [stat_err,stat_err_m]+all_hists
  
  limit = max([abs(h.GetMinimum()) for h in all_hists] + [abs(h.GetMaximum()) for h in all_hists])
  
  for h in all_hists:
    h.SetMinimum(-0.05)
    h.SetMaximum(1.1*limit)
    h.Draw("same"+plot_opt)
  
  legend.Draw()
  canv.Print(name+"_onesigma_ratio_to_nominal.eps")
  
def ratio_syst_to_stat_err(rootfile,name):
  # Make the plot with ratio of syst to stat_err
  canv = ROOT.TCanvas()
  canv.SetBottomMargin(bottom_margin)
  if use_log_plots: canv.SetLogx()
  legend = ROOT.TLegend(0.55,0.75,0.9,0.88);
  legend.SetFillStyle(0);
  legend.SetBorderSize(0);
  legend.SetTextSize(0.03);
  
  nominal = rootfile.Get("nominal")
  stat_err = ROOT.TH1D(ROOT.TH1D(rootfile.Get(param_sets[0]+"_pluserr") ))
  stat_err.SetLineColor(colors[0])
  for ibin in xrange(1,nominal.GetNbinsX()):
    try: err = 1./math.sqrt(nominal.GetBinContent(ibin))
    except ZeroDivisionError: err=0.
    stat_err.SetBinContent(ibin,err)
    
  all_hists = []
  
  for ii,par_set in enumerate(param_sets):
    if par_set != "all": continue
    try:
      hplus  = ROOT.TH1D(rootfile.Get(par_set+"_pluserr") )
    except TypeError:
      continue  # Osc params won't load for unoscillated
    hists = [hplus]
    for h in hists:
      #~ h.SetTitle("Ratio of effect of systematic to statistical error ("+name+")")
      h.SetTitle("")
      h.GetXaxis().SetRangeUser(0.,10.)
      h.GetYaxis().SetTitle("Systematic effect over statistical effect")
      h.GetYaxis().SetTitleOffset(title_offset)
      h.SetLineColor(colors[ii+1])
      h.SetLineWidth(linewidth)
      
      h.Divide(stat_err)
      
      all_hists.append(h)
      
    legend.AddEntry(hplus,readable[par_set],"L")
  
  limit = max([abs(h.GetMinimum()) for h in all_hists] + [abs(h.GetMaximum()) for h in all_hists])
  
  for h in all_hists:
    h.SetMinimum(-0.05)
    h.SetMaximum(1.1*limit)
    h.Draw("same"+plot_opt)
  
  legend.Draw()
  canv.Print(name+"_onesigma_ratio_to_statistical.eps")
  
def syst_err_ratio_to_total_syst(rootfile,name):
  # Make the plot with ratio of each syst group to all systematics
  canv = ROOT.TCanvas()
  canv.SetBottomMargin(bottom_margin)
  if use_log_plots: canv.SetLogx()
  legend = ROOT.TLegend(0.55,0.75,0.9,0.88);
  legend.SetFillStyle(0);
  legend.SetBorderSize(0);
  legend.SetTextSize(0.03);
    
  all_hists = []
  
  hall = ROOT.TH1D(rootfile.Get("all_pluserr"))
  
  for ii,par_set in enumerate(param_sets):
    if par_set == "all": continue
    try:
      hplus  = ROOT.TH1D(rootfile.Get(par_set+"_pluserr") )
    except TypeError:
      continue  # Osc params won't load for unoscillated
    hists = [hplus]
    for h in hists:
      #~ h.SetTitle("Ratio of effect of each systematic source to total systematic ("+name+")")
      h.SetTitle("")
      h.GetXaxis().SetRangeUser(0.,10.)
      h.GetYaxis().SetTitle("This syst / All syst")
      h.GetYaxis().SetTitleOffset(title_offset)
      h.SetLineColor(colors[ii+1])
      h.SetLineWidth(linewidth)
      
      h.Divide(hall)
      
      all_hists.append(h)
      
    legend.AddEntry(hplus,readable[par_set],"L")
  
  limit = max([abs(h.GetMinimum()) for h in all_hists] + [abs(h.GetMaximum()) for h in all_hists])
  
  for h in all_hists:
    h.SetMinimum(-0.05)
    h.SetMaximum(1.3)
    h.Draw("same"+plot_opt)
  
  legend.Draw()
  canv.Print(name+"_onesigma_ratio_to_syst.eps")

class T2KAnaError(Exception):
  pass

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make plots of size of systematic/statistical effects")
  parser.add_option('-o', dest='osc'  ,default=default_oscillated_input_file  ,help='The unoscillated input ROOT file')
  parser.add_option('-u', dest='unosc',default=default_unoscillated_input_file,help='The oscillated input ROOT file')
  
  options,other = parser.parse_args()
  
  main(options.unosc,options.osc)
