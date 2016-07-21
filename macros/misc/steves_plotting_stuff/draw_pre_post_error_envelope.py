#!/usr/bin/env python

# Script to make plots of relative error envelope sizes
# eg plots 9-12 in TN-183
# Takes output of macros/spectrum/test_sk_ereco_spectrum_twk_sys_throw_runAll.py

#Steve Dennis, October 2013

import ROOT
import math
import os
import glob

R = ROOT

max_energy = 4.5
linewidth = 3
title_offset = 0.8
left_margin   = 0.16
bottom_margin = 0.18
right_margin  = 0.05
use_log_plots = False
title_size = 0.08
label_size = 0.07

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
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

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
      
      h.GetXaxis().SetTitleSize(title_size)
      h.GetYaxis().SetTitleSize(title_size)
      h.GetXaxis().SetLabelSize(label_size)
      h.GetYaxis().SetLabelSize(label_size)
    
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
      
    legend = ROOT.TLegend(0.7,0.6,0.93,0.92);
    legend.SetFillStyle(0);
    legend.SetBorderSize(0);
    #~ legend.SetTextSize(0.04);
  
    limit = max([abs(h.GetMinimum()) for h in par_hists] + [abs(h.GetMaximum()) for h in par_hists])
    
    canv = ROOT.TCanvas()
    canv.SetLeftMargin  (left_margin)
    canv.SetBottomMargin(bottom_margin)
    canv.SetRightMargin (right_margin )
    if use_log_plots: canv.SetLogx()
    
    for h in par_hists:
      h.SetMinimum(0.)
      h.SetMaximum(1.05*limit)
      
    hplus_prefit.GetYaxis().SetTitleOffset(0.8)
    hplus_prefit.SetLineColor(ROOT.kBlue)
    
    hplus_prefit.Draw("AXIS")
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
    
    #~ if is_constrained[par_set]: legend.Draw()
    
    canv.Print(name+"_"+par_set+"_envelope.eps")
    
    # Now do it with ratios
    ratio_plus_prefit   = R.TH1D(hplus_prefit)
    ratio_minus_prefit  = R.TH1D(hminus_prefit)
    ratio_plus_postfit  = R.TH1D(hplus_postfit)
    ratio_minus_postfit = R.TH1D(hminus_postfit)
    
    ratio_plus_prefit  .Divide(prefit_nominal)
    ratio_minus_prefit .Divide(prefit_nominal)
    ratio_plus_postfit .Divide(nominal_hist)
    ratio_minus_postfit.Divide(nominal_hist)
    for i in xrange(nominal_hist.GetXaxis().GetNbins()):
      if abs(nominal_hist.GetBinContent(i+1))<conv.kASmallNum:
        ratio_plus_postfit .SetBinContent(i+1,1)
        ratio_minus_postfit.SetBinContent(i+1,1)
      if abs(prefit_nominal.GetBinContent(i+1))<conv.kASmallNum:
        ratio_plus_prefit .SetBinContent(i+1,1)
        ratio_minus_prefit.SetBinContent(i+1,1)
    
    ratio_hists = [(ratio_plus_prefit  , ratio_minus_prefit ) ,
                   (ratio_plus_postfit , ratio_minus_postfit)]
    
    ratio_graph_prefit  = ROOT.TGraphAsymmErrors(ratio_plus_prefit  .GetXaxis().GetNbins())
    ratio_graph_postfit = ROOT.TGraphAsymmErrors(ratio_plus_postfit .GetXaxis().GetNbins())
    
    ratio_graphs = [ratio_graph_prefit  ,
                    ratio_graph_postfit , ]
    
    for (hp,hm),g in zip(ratio_hists,ratio_graphs):
      for i in xrange(h.GetXaxis().GetNbins()):
        low_edge  = ratio_plus_prefit.GetBinLowEdge(i+1)
        high_edge = low_edge + ratio_plus_prefit.GetBinWidth(i+1)
        g.SetPoint      (i,ratio_plus_prefit.GetBinCenter(i+1),1)
        g.SetPointEYlow (i,abs(hm.GetBinContent(i+1)-1))
        g.SetPointEYhigh(i,abs(hp.GetBinContent(i+1)-1))
        g.SetPointEXlow (i,ratio_plus_prefit.GetBinWidth(i+1)*0.5 - 0.005)
        g.SetPointEXhigh(i,ratio_plus_prefit.GetBinWidth(i+1)*0.5 - 0.005)
        
      g.SetFillStyle(1001)
      g.SetLineStyle(0)
      g.SetLineWidth(0)
      g.GetXaxis().SetRangeUser(0,max_energy)
      CopyHistAxisFormatting(hplus_prefit,g)
      g.GetYaxis().SetTitle("Ratio of # Total Events")
      
    ratio_graph_prefit .SetFillColor(R.kBlue)
    ratio_graph_postfit.SetFillColor(R.kRed)
    ratio_graph_prefit .SetLineColor(R.kBlue)
    ratio_graph_postfit.SetLineColor(R.kRed)
    
    ratio_graph_prefit .SetLineWidth(3)
    ratio_graph_postfit.SetLineWidth(3)
    ratio_graph_prefit .SetFillStyle(1001)
    ratio_graph_postfit.SetFillStyle(1001)
    
    canv2 = ROOT.TCanvas()
    canv2.SetLeftMargin  (left_margin)
    canv2.SetBottomMargin(bottom_margin)
    canv2.SetRightMargin (right_margin )
    if use_log_plots: canv2.SetLogx()
    canv2.cd()
    
    ratio_max = 1.3*ratio_plus_prefit .GetMaximum()
    
    ratio_graph_postfit.SetMaximum(ratio_max)
    ratio_graph_postfit.SetMinimum(2-ratio_max)
    
    ratio_graph_prefit .Draw("a2")
    ratio_graph_postfit.Draw("2")
    
    legend2 = ROOT.TLegend(0.5,0.85,0.93,0.94);
    legend2.SetFillStyle(0);
    legend2.SetBorderSize(0);
    legend2.SetNColumns(2)
    legend2.AddEntry(ratio_graph_prefit  ,"Prefit" ,"F")
    legend2.AddEntry(ratio_graph_postfit ,"Postfit","F")
    
    legend2.Draw("same")
    
    line_at_1 = ROOT.TLine(0,1,max_energy,1)
    line_at_1.SetLineWidth(2)
    line_at_1.Draw("same")
    
    canv2.Update()
    
    canv2.Print(name+"_"+par_set+"_envelope_ratio.eps")
    
    
    # Print out some event rates
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
  
def CopyHistAxisFormatting(hfrom,hto):
  for (afrom,ato) in [(hfrom.GetXaxis(),hto.GetXaxis()),
                      (hfrom.GetYaxis(),hto.GetYaxis())]:
    ato.SetTitle      (afrom.GetTitle      ())
    ato.SetTitleOffset(afrom.GetTitleOffset())
    ato.SetTitleSize  (afrom.GetTitleSize  ())
    ato.SetLabelSize  (afrom.GetLabelSize  ())
    ato.SetLabelOffset(afrom.GetLabelOffset())
    
  return

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
