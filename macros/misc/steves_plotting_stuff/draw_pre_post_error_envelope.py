#!/usr/bin/python

# Script to make plots of relative error envelope sizes
# eg plots 9-12 in TN-183
# Takes output of macros/spectrum/test_sk_ereco_spectrum_twk_sys_throw_runAll.py

#Steve Dennis, October 2013

import ROOT
import math
import os
import glob

max_energy = 4.5
linewidth = 3
title_offset = 0.8
bottom_margin = 0.15
use_log_plots = False

to3dp = "{0:.3f}"

default_oscillated_input_file = (
  "sk_ereco_spectrum_RHC1Rmu-twk_sys_throw_sinsqtheta[12_0.304,13_0.042,23_0.524],dcp_1.910,dmsq_Normal[12_0.00007530,23_0.002510]eV2.root")
default_unoscillated_input_file = (
  "sk_ereco_spectrum_RHC1Rmu-twk_sys_throw_sinsqtheta[12_0.000,13_0.000,23_0.000],dcp_0.000,dmsq_Normal[12_0.00007530,23_0.002510]eV2.root")

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

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

is_constrained  = dict()
is_constrained ["SKDet_FSI"         ] = True
is_constrained ["BANFF_xsec"        ] = True
is_constrained ["BANFF_flux"        ] = True
is_constrained ["all"               ] = True

plot_opt = ""

try:
  ROOT.t2kana
except:
  #~ ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/loadincs.C")
  #~ ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/loadlibs.C")
  ROOT.gROOT.ProcessLine(".L $T2KOSCANA/root/root2k.C");
  pass

  

def main(sample,unoscillated_input_filename, oscillated_input_filename):
  f_unosc = ROOT.TFile(unoscillated_input_filename)
  if not f_unosc:
    raise T2KAnaError("Could not open unoscillated input file")
    
  pre_and_postfit_envelope(sample,f_unosc,"unoscillated")
  
  f_osc   = ROOT.TFile(oscillated_input_filename)
  if not f_unosc:
    raise T2KAnaError("Could not open oscillated input file")
    
  pre_and_postfit_envelope(sample,f_osc  ,"oscillated" )
  

def pre_and_postfit_envelope(sample,rootfile,name):
  
  assert(rootfile)
  
  nh = rootfile.Get(sample+"_nominal")
  nominal_hist = ROOT.TH1D(nh)
  nominal_hist.SetTitle("")
  nominal_hist.GetXaxis().SetRangeUser(0.,max_energy)
  nominal_hist.GetXaxis().SetTitle("E_{reco}")
  nominal_hist.GetYaxis().SetTitle("Events per bin")
  nominal_hist.GetYaxis().SetTitleOffset(title_offset)
  nominal_hist.SetLineColor(ROOT.kBlack)
  nominal_hist.SetLineStyle(1)
  nominal_hist.SetLineWidth(linewidth)
  
  nom_events = nominal_hist.Integral()
  
  prefit_nominal = ROOT.TH1D(rootfile.Get(sample+"_nominalprefit"))
  prefit_nominal.SetLineStyle(1)
  prefit_nominal.SetLineWidth(linewidth)
  prefit_nominal.SetLineColor(ROOT.kBlue)
  
  all_hists = [nominal_hist,prefit_nominal]
  
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
      #~ h.SetLineColor(ROOT.kRed)
      h.SetLineStyle(2)
      h.SetLineWidth(linewidth)
      h.GetXaxis().SetRangeUser(0.,max_energy)
      
      if not is_constrained[par_set]:
        h.SetLineStyle(0) # Don't draw prefit lines for unconstrained parameters
        h.SetLineWidth(0) # Don't draw prefit lines for unconstrained parameters
        
    for h in [hplus_postfit,hminus_postfit]:
      h.SetLineStyle(1)
      h.SetLineWidth(0)
      h.GetXaxis().SetRangeUser(0.,max_energy)
      print h.GetTitle(),h.GetLineWidth()
      
    legend = ROOT.TLegend(0.6,0.6,0.9,0.9);
    legend.SetFillStyle(0);
    legend.SetBorderSize(0);
    legend.SetTextSize(0.04);
  
    limit = max([abs(h.GetMinimum()) for h in par_hists] + [abs(h.GetMaximum()) for h in par_hists])
    
    canv = ROOT.TCanvas()
    canv.SetBottomMargin(bottom_margin)
    if use_log_plots: canv.SetLogx()
    
    for h in par_hists:
      h.SetMinimum(0.)
      h.SetMaximum(1.05*limit)
      
    hplus_prefit.GetYaxis().SetTitleOffset(0.8)
    hplus_prefit.SetLineColor(ROOT.kBlue)
    
    hplus_prefit.Draw()
    hplus_postfit.SetFillColor(ROOT.kRed)
    hplus_postfit.SetFillStyle(1001)
    hplus_postfit.Draw("same")
    hminus_postfit.SetFillColor(10)
    hminus_postfit.SetFillStyle(1001)
    hminus_postfit.Draw("same")
    hminus_prefit.SetFillColor(10)
    hminus_prefit.SetFillStyle(1001)
    hminus_prefit.SetLineColor(ROOT.kBlue)
    hminus_prefit.Draw("same")
    hplus_prefit.Draw("same")
    nominal_hist.Draw("same")
    prefit_nominal.Draw("same")
    hplus_prefit.Draw("same")
    
    legend_prefit = ROOT.TH1D(hminus_prefit)
    legend_prefit.SetLineStyle(prefit_nominal.GetLineStyle())
    legend_prefit.SetLineColor(prefit_nominal.GetLineColor())
    legend_prefit.SetLineWidth(prefit_nominal.GetLineWidth())
    
    legend_postfit = ROOT.TH1D(hplus_postfit)
    legend_postfit.SetLineStyle(nominal_hist.GetLineStyle())
    legend_postfit.SetLineColor(nominal_hist.GetLineColor())
    legend_postfit.SetLineWidth(nominal_hist.GetLineWidth())
    
    all_hists.extend([legend_prefit,legend_postfit])
    
    legend.AddEntry(legend_prefit   ,"Prefit" ,"FL")
    legend.AddEntry(legend_postfit  ,"Postfit","FL")
    
    if is_constrained[par_set]: legend.Draw()
    
    canv.Print(name+"_"+par_set+"_envelope.eps")
    
    nom_events     = nominal_hist  .Integral()
    p1_events_pre  = hplus_prefit  .Integral()
    m1_events_pre  = hminus_prefit .Integral()
    p1_events_post = hplus_postfit .Integral()
    m1_events_post = hminus_postfit.Integral()
    
    p1_percent_pre  = 100.*( (p1_events_pre -nom_events)/nom_events )
    m1_percent_pre  = 100.*( (m1_events_pre -nom_events)/nom_events )
    p1_percent_post = 100.*( (p1_events_post-nom_events)/nom_events )
    m1_percent_post = 100.*( (m1_events_post-nom_events)/nom_events )
    
    items = [name+"_"+par_set,to3dp.format(nom_events)]
    items.append(to3dp.format(p1_events_pre ) +" ("+to3dp.format(p1_percent_pre )+"%)")
    items.append(to3dp.format(m1_events_pre ) +" ("+to3dp.format(m1_percent_pre )+"%)")
    items.append(to3dp.format(p1_events_post) +" ("+to3dp.format(p1_percent_post)+"%)")
    items.append(to3dp.format(m1_events_post) +" ("+to3dp.format(m1_percent_post)+"%)")
    print "\t".join(items)
    
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
  
  if other:
    print "Bad arguments: ",other
    exit(-1)
  
  main(options.sample,options.unosc,options.osc)
