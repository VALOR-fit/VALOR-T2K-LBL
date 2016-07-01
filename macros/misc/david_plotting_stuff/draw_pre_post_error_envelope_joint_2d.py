#!/usr/bin/python

# Script to make plots of relative error envelope sizes
# eg plots 9-12 in TN-183
# Takes output of macros/spectrum/test_sk_ereco_spectrum_twk_sys_throw_runAll_nuebar.py

#Steve Dennis, October 2013

import ROOT
import math
import os
import glob

R = ROOT # Steve

# max_energy = 4.5
# linewidth = 3
# title_offset = 0.8
# bottom_margin = 0.15
# use_log_plots = False

#max_energy = 1.25
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
  "RHC1Re/sk_ereco_spectrum_RHC1Re-twk_sys_throw_dcp_-99.000.root"
)
default_unoscillated_input_file = (
  "RHC1Re/sk_ereco_spectrum_RHC1Re-twk_sys_throw_dcp_-99.000.root"
)

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kMagenta+1,  ROOT.kGreen+2, ROOT.kOrange+7, ROOT.kCyan+2,]

param_sets = [
  "SKDet_FSI"     ,
  "BANFF_xsec"    ,
  "BANFF_flux"    ,
  "BANFF"         ,
  #"skdetfsi_banff",
  "Osc_param"     ,
  "all"           ,
]

readable = dict()
readable["SKDet_FSI"         ] = "SK Detector + FSI + SI"
readable["BANFF_xsec"        ] = "Cross-section"
readable["BANFF_flux"        ] = "Flux"
readable["BANFF"             ] = "BANFF all"
#readable["skdetfsi_banff"    ] = "Flux + XSec + SKDetector + FSI + SI"
readable["Osc_param"         ] = "Oscillations"
readable["all"               ] = "All systematics"

is_constrained  = dict()
is_constrained ["SKDet_FSI"         ] = True
is_constrained ["BANFF_xsec"        ] = True
is_constrained ["BANFF_flux"        ] = True
is_constrained ["BANFF"             ] = True
#is_constrained ["skdetfsi_banff"    ] = True
is_constrained ["Osc_param"         ] = True
is_constrained ["all"               ] = True

plot_opt = ""

try:
  ROOT.t2kana
  #ROOT.valor # Steve
except:
  ROOT.gROOT.ProcessLine(".L $T2KOSCANA/root/root2k.C");
  pass

#valor = ROOT.valor # Steve
#util =  valor.util # Steve
#conv =  valor.conv # Steve
#opt  =  valor.Opt.Instance() # Steve
#env  =  valor.Env.Instance() # Steve

  

def main(sample,unoscillated_input_filename, oscillated_input_filename,usemedian,doprefit,tag):
  
  f_unosc = ROOT.TFile(unoscillated_input_filename)
  if not f_unosc:
    raise T2KAnaError("Could not open unoscillated input file")
  pre_and_postfit_envelope(sample,f_unosc,"unoscillated",usemedian,doprefit,tag)
  
  f_osc   = ROOT.TFile(oscillated_input_filename)
  if not f_unosc:
    raise T2KAnaError("Could not open oscillated input file")
  pre_and_postfit_envelope(sample,f_osc  ,"oscillated", usemedian, doprefit,tag)
  

def pre_and_postfit_envelope(sample,rootfile,name,usemedian,doprefit,tag):
  
  assert(rootfile)

  if sample=="RHC1Re" or sample=="1Re" :
    max_energy = 1.25
  elif sample=="RHC1Rmu" or sample=="1Rmu" :
    max_energy = 2.9
    
  nh = rootfile.Get(sample+"_nominal")
  nominal_hist = ROOT.TH1D(nh)
  nominal_hist.SetTitle("")
  nominal_hist.GetXaxis().SetRangeUser(0.,max_energy)
  nominal_hist.GetXaxis().SetTitle("E_{reco} (GeV)")
  nominal_hist.GetYaxis().SetTitle("Events per bin")
  nominal_hist.GetYaxis().SetTitleOffset(title_offset)
  nominal_hist.SetLineColor(ROOT.kBlack)
  nominal_hist.SetLineStyle(1)
  nominal_hist.SetLineWidth(linewidth)
  
  nom_events = nominal_hist.Integral()
  
  if doprefit==True :
    prefit_nominal = ROOT.TH1D(rootfile.Get(sample+"_nominalprefit"))
    prefit_nominal.SetLineStyle(1)
    prefit_nominal.SetLineWidth(linewidth)
    prefit_nominal.SetLineColor(ROOT.kBlue)
  

  canv_nsk_prefit = ROOT.TCanvas()
  canv_nsk_postfit = ROOT.TCanvas()
  legend_nsk_prefit = ROOT.TLegend(0.7,0.6,0.9,0.9);
  legend_nsk_prefit.SetFillStyle(0);
  legend_nsk_prefit.SetBorderSize(0);
  legend_nsk_prefit.SetTextSize(0.04);
  legend_nsk_postfit = ROOT.TLegend(0.7,0.6,0.9,0.9);
  legend_nsk_postfit.SetFillStyle(0);
  legend_nsk_postfit.SetBorderSize(0);
  legend_nsk_postfit.SetTextSize(0.04);

  if doprefit==True :
    all_hists = [nominal_hist,prefit_nominal]
  else :
    all_hists = [nominal_hist]

  all_hists_nsk_prefit = [] # declare it as global only for legend
  all_hists_nsk_postfit = [] # declare it as global only for legend

  for ii,par_set in enumerate(param_sets):
      
    par_hists = [nominal_hist]

    ##### Draw histo of total # of SK events for each group of systematics

    if doprefit==True :
      hname = sample+"_"+par_set+"prefit"+"_nsk"
      h = rootfile.Get(hname)
      if not h:
        print "Failed to get histogram: ",hname ; exit()
      prefit_nsk_hist = ROOT.TH1D( h )
      prefit_nsk_hist.SetLineColor(ii+1)
      prefit_nsk_hist.SetLineStyle(1)
      prefit_nsk_hist.SetLineWidth(linewidth)
      canv_nsk_prefit.cd()
      if ii==0 :
        prefit_nsk_hist.DrawCopy()
      else :
        prefit_nsk_hist.DrawCopy("same")
      all_hists_nsk_prefit.append(prefit_nsk_hist)
      legend_nsk_prefit.AddEntry(all_hists_nsk_prefit[ii],readable[par_set],"L")
      legend_nsk_prefit.Draw()

    hname = sample+"_"+par_set+"_nsk"
    h = rootfile.Get(hname)
    if not h:
      print "Failed to get histogram: ",hname ; exit()
    postfit_nsk_hist = ROOT.TH1D( h )
    postfit_nsk_hist.SetLineColor(ii+1)
    postfit_nsk_hist.SetLineStyle(1)
    postfit_nsk_hist.SetLineWidth(linewidth)
    canv_nsk_postfit.cd()
    if ii==0 :
      postfit_nsk_hist.DrawCopy()
    else :
      postfit_nsk_hist.DrawCopy("same")
    all_hists_nsk_postfit.append(postfit_nsk_hist)
    legend_nsk_postfit.AddEntry(all_hists_nsk_postfit[ii],readable[par_set],"L")
    legend_nsk_postfit.Draw()

    #######

    if doprefit==True :
      #hname = sample+"_"+par_set+"prefit"+"_average"
      hname = sample+"_"+par_set+"prefit"+"_median"
      h = rootfile.Get(hname)
      if not h:
        print "Failed to get histogram: ",hname ; exit()
      prefit_average_hist = ROOT.TH1D( h )
      prefit_average_hist.SetLineColor(ROOT.kCyan)
      prefit_average_hist.SetLineStyle(1)
      prefit_average_hist.SetLineWidth(linewidth)

    #hname = sample+"_"+par_set+"_average"
    hname = sample+"_"+par_set+"_median"
    h = rootfile.Get(hname)
    if not h:
      print "Failed to get histogram: ",hname ; exit()
    postfit_average_hist = ROOT.TH1D( h )
    postfit_average_hist.SetLineColor(ROOT.kYellow)
    postfit_average_hist.SetLineStyle(1)
    postfit_average_hist.SetLineWidth(linewidth)

    if doprefit==True :

      if usemedian==False :
        hname = sample+"_"+par_set+"prefit"+"_envpluserr"
      elif usemedian==True :
        #hname = sample+"_"+par_set+"prefit"+"_up_average"
        hname = sample+"_"+par_set+"prefit"+"_up_median"
        h = rootfile.Get(hname)
        if not h:
          print "Failed to get histogram: ",hname ; exit()
        hplus_prefit     = ROOT.TH1D( h )
    
        if usemedian==False :
          hname = sample+"_"+par_set+"prefit"+"_envminuserr"
        elif usemedian==True :
          #hname = sample+"_"+par_set+"prefit"+"_down_average"
          hname = sample+"_"+par_set+"prefit"+"_down_median"
          h = rootfile.Get(hname)
          if not h:
            print "Failed to get histogram: ",hname ; exit()
          hminus_prefit     = ROOT.TH1D( h )
    
    if usemedian==False :
      hname = sample+"_"+par_set+"_envpluserr"
    elif usemedian==True :
      #hname = sample+"_"+par_set+"_up_average"
      hname = sample+"_"+par_set+"_up_median"
    h = rootfile.Get(hname)
    if not h:
      print "Failed to get histogram: ",hname ; exit()
    hplus_postfit     = ROOT.TH1D( h )
    
    if usemedian==False :
      hname = sample+"_"+par_set+"_envminuserr"
    elif usemedian==True :
      #hname = sample+"_"+par_set+"_down_average"
      hname = sample+"_"+par_set+"_down_median"
    h = rootfile.Get(hname)
    if not h:
      print "Failed to get histogram: ",hname ; exit()
    hminus_postfit     = ROOT.TH1D( h )
    
    if doprefit==True :
      hists = [hplus_prefit,hminus_prefit,hplus_postfit,hminus_postfit]
    else :
      hists = [hplus_postfit,hminus_postfit]

    for h in hists:
      par_hists.append(h)
      all_hists.append(h)
    
    if doprefit==True :
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
      
    legend = ROOT.TLegend(0.65,0.6,0.8,0.9);
    legend.SetFillStyle(0);
    legend.SetBorderSize(0);
    legend.SetTextSize(0.05);
  
    limit = max([abs(h.GetMinimum()) for h in par_hists] + [abs(h.GetMaximum()) for h in par_hists])
    
    canv = ROOT.TCanvas()
    canv.SetLeftMargin  (left_margin) # Steve
    canv.SetBottomMargin(bottom_margin)
    canv.SetRightMargin (right_margin ) # Steve
    if use_log_plots: canv.SetLogx()
    
    for h in par_hists:
      h.SetMinimum(0.)
      h.SetMaximum(1.05*limit)
      
    if doprefit==True :
      hplus_prefit.SetLineColor(ROOT.kBlue)
      hplus_prefit.GetYaxis().SetTitleOffset(0.8)
      hplus_prefit.SetTitle("")
      hplus_prefit.GetXaxis().SetRangeUser(0.,max_energy)
      hplus_prefit.GetXaxis().SetTitle("E_{reco} (GeV)")
      hplus_prefit.GetYaxis().SetTitle("Events per bin")
      hplus_prefit.GetYaxis().SetTitleOffset(title_offset)
      #hplus_prefit.SetLineStyle(1)
      hplus_prefit.SetLineWidth(linewidth)

    hplus_postfit.GetYaxis().SetTitleOffset(0.8)
    hplus_postfit.SetTitle("")
    hplus_postfit.GetXaxis().SetRangeUser(0.,max_energy)
    hplus_postfit.GetXaxis().SetTitle("E_{reco} (GeV)")
    hplus_postfit.GetYaxis().SetTitle("Events per bin")
    hplus_postfit.GetYaxis().SetTitleOffset(title_offset)

    hplus_postfit.SetFillColor(ROOT.kRed)
    hplus_postfit.SetFillStyle(1001)
    
    hminus_postfit.SetFillColor(10)
    hminus_postfit.SetFillStyle(1001)

    if doprefit==True :
      hminus_prefit.SetFillColor(10)
      hminus_prefit.SetFillStyle(1001)
      hminus_prefit.SetLineColor(ROOT.kBlue)

    if doprefit==True : 
      hplus_prefit.DrawCopy()
      hplus_postfit.DrawCopy("same")
    else :
      hplus_postfit.DrawCopy("")
    hminus_postfit.DrawCopy("same")
    if doprefit==True : 
      hminus_prefit.DrawCopy("same")
    if doprefit==True : 
      hplus_prefit.DrawCopy("same")
    nominal_hist.DrawCopy("same")
    #postfit_average_hist.DrawCopy("same")
    if doprefit==True : 
      prefit_nominal.DrawCopy("same")
      ##prefit_average_hist.DrawCopy("same") 
    if doprefit==True : 
      hplus_prefit.DrawCopy("same")
    
    ROOT.gPad.RedrawAxis();

    if doprefit==True :
      legend_prefit = ROOT.TH1D(hminus_prefit)
      legend_prefit.SetLineStyle(prefit_nominal.GetLineStyle())
      legend_prefit.SetLineColor(prefit_nominal.GetLineColor())
      legend_prefit.SetLineWidth(prefit_nominal.GetLineWidth())
    
    legend_postfit = ROOT.TH1D(hplus_postfit)
    legend_postfit.SetLineStyle(nominal_hist.GetLineStyle())
    legend_postfit.SetLineColor(nominal_hist.GetLineColor())
    legend_postfit.SetLineWidth(nominal_hist.GetLineWidth())
    
    if doprefit==True :
      all_hists.extend([legend_prefit,legend_postfit])
    else :
      all_hists.extend([legend_postfit])

    
    if doprefit==True : legend.AddEntry(legend_prefit   ,"Nominal error" ,"FL")
    legend.AddEntry(legend_postfit  ,"#pm 1#sigma error","FL")
    
    if is_constrained[par_set]: legend.Draw()
    
    canv.Update()

    canv.Print(tag+"_"+sample+name+"_"+par_set+"_envelope.eps")
    







    # Now do it with ratios

    # ratio_plus_prefit   = R.TH1D(hplus_prefit)
    # ratio_minus_prefit  = R.TH1D(hminus_prefit)
    # ratio_plus_postfit  = R.TH1D(hplus_postfit)
    # ratio_minus_postfit = R.TH1D(hminus_postfit)
 
    # if usemedian==True:
    #   ratio_plus_prefit  .Divide(prefit_average_hist)
    #   ratio_minus_prefit .Divide(prefit_average_hist)
    #   ratio_plus_postfit .Divide(postfit_average_hist)
    #   ratio_minus_postfit.Divide(postfit_average_hist)
    # else :
    #   ratio_plus_prefit  .Divide(prefit_nominal)
    #   ratio_minus_prefit .Divide(prefit_nominal)
    #   ratio_plus_postfit .Divide(nominal_hist)
    #   ratio_minus_postfit.Divide(nominal_hist)

    # for i in xrange(nominal_hist.GetXaxis().GetNbins()):
    #   if abs(nominal_hist.GetBinContent(i+1))<0.001:
    #     ratio_plus_postfit .SetBinContent(i+1,1)
    #     ratio_minus_postfit.SetBinContent(i+1,1)
    #   if abs(prefit_nominal.GetBinContent(i+1))<0.001:
    #     ratio_plus_prefit .SetBinContent(i+1,1)
    #     ratio_minus_prefit.SetBinContent(i+1,1)
    
    # ratio_hists = [(ratio_plus_prefit  , ratio_minus_prefit ) ,
    #                (ratio_plus_postfit , ratio_minus_postfit)]
    
    # ratio_graph_prefit  = ROOT.TGraphAsymmErrors(ratio_plus_prefit  .GetXaxis().GetNbins())
    # ratio_graph_postfit = ROOT.TGraphAsymmErrors(ratio_plus_postfit .GetXaxis().GetNbins())
    
    # ratio_graphs = [ratio_graph_prefit  ,
    #                 ratio_graph_postfit , ]
    
    # for (hp,hm),g in zip(ratio_hists,ratio_graphs):
    #   for i in xrange(h.GetXaxis().GetNbins()):
    #     low_edge  = ratio_plus_prefit.GetBinLowEdge(i+1)
    #     high_edge = low_edge + ratio_plus_prefit.GetBinWidth(i+1)
    #     g.SetPoint      (i,ratio_plus_prefit.GetBinCenter(i+1),1)
    #     g.SetPointEYlow (i,abs(hm.GetBinContent(i+1)-1))
    #     g.SetPointEYhigh(i,abs(hp.GetBinContent(i+1)-1))
    #     g.SetPointEXlow (i,ratio_plus_prefit.GetBinWidth(i+1)*0.5 - 0.005)
    #     g.SetPointEXhigh(i,ratio_plus_prefit.GetBinWidth(i+1)*0.5 - 0.005)
        
    #   g.SetFillStyle(1001)
    #   g.SetLineStyle(0)
    #   g.SetLineWidth(0)
    #   g.GetXaxis().SetRangeUser(0,max_energy)
    #   CopyHistAxisFormatting(hplus_prefit,g)
    #   g.GetYaxis().SetTitle("Ratio of # Total Events")
      
    # ratio_graph_prefit .SetFillColor(R.kBlue)
    # ratio_graph_postfit.SetFillColor(R.kYellow)
    # ratio_graph_prefit .SetLineColor(R.kBlue)
    # ratio_graph_postfit.SetLineColor(R.kYellow)    
    # #ratio_graph_prefit .SetLineWidth(3)
    # #ratio_graph_postfit.SetLineWidth(3)
    
    # #ratio_graph_prefit.SetFillStyle(1001)
    # ratio_graph_prefit.SetFillStyle(3001)
    # ratio_graph_prefit.SetLineColor(R.kBlue)
    # ratio_graph_prefit.SetLineWidth(2)
        
    # #ratio_graph_postfit.SetFillStyle(1001)
    # ratio_graph_postfit.SetFillStyle(3001)
    # ratio_graph_postfit.SetLineWidth(2)

    # canv2 = ROOT.TCanvas()
    # canv2.SetLeftMargin  (left_margin)
    # canv2.SetBottomMargin(bottom_margin)
    # canv2.SetRightMargin (right_margin )
    # if use_log_plots: canv2.SetLogx()
    # canv2.cd()
    
    # ratio_max = 1.3*ratio_plus_prefit .GetMaximum()
    # ratio_graph_postfit.SetMaximum(ratio_max)
    # ratio_graph_postfit.SetMinimum(2-ratio_max)
    # ratio_graph_prefit .Draw("a2")
    # ratio_graph_postfit.Draw("2")

    # # ratio_max = 2.;
    # # ratio_min = 0.;
    # # ratio_plus_postfit.SetMaximum(ratio_max)
    # # ratio_minus_postfit.SetMinimum(ratio_min)
    # # ratio_plus_postfit.DrawCopy("")
    # # ratio_minus_postfit.DrawCopy("same")
    
    # legend2 = ROOT.TLegend(0.6,0.75,0.93,0.92);
    # legend2.SetFillStyle(0);
    # legend2.SetBorderSize(0);
    # legend2.AddEntry(ratio_graph_prefit  ,"Prefit" ,"F")
    # legend2.AddEntry(ratio_graph_postfit ,"Postfit","F")
    # legend2.Draw("same")
    
    # line_at_1 = ROOT.TLine(0,1,max_energy,1)
    # line_at_1.SetLineWidth(2)
    # line_at_1.Draw("same")
    
    # canv2.Update()
    # canv2.Print(name+"_"+par_set+"_envelope_ratio.eps")






    

    if doprefit==True :
      ratio_plus_prefit   = R.TH1D(hplus_prefit)
      ratio_minus_prefit  = R.TH1D(hminus_prefit)
    ratio_plus_postfit  = R.TH1D(hplus_postfit)
    ratio_minus_postfit = R.TH1D(hminus_postfit)
 
    if usemedian==True:
      if doprefit==True :
        ratio_plus_prefit  .Divide(prefit_average_hist)
        ratio_minus_prefit .Divide(prefit_average_hist)
      ratio_plus_postfit .Divide(postfit_average_hist)
      ratio_minus_postfit.Divide(postfit_average_hist)
    else :
      if doprefit==True :
        ratio_plus_prefit  .Divide(prefit_nominal)
        ratio_minus_prefit .Divide(prefit_nominal)
      ratio_plus_postfit .Divide(nominal_hist)
      ratio_minus_postfit.Divide(nominal_hist)

    for i in xrange(nominal_hist.GetXaxis().GetNbins()):
      if abs(nominal_hist.GetBinContent(i+1))<0.001:
        ratio_plus_postfit .SetBinContent(i+1,1)
        ratio_minus_postfit.SetBinContent(i+1,1)
      if doprefit==True :
        if abs(prefit_nominal.GetBinContent(i+1))<0.001:
          ratio_plus_prefit .SetBinContent(i+1,1)
          ratio_minus_prefit.SetBinContent(i+1,1)
    
    if doprefit==True :
      ratio_hists = [(ratio_plus_prefit  , ratio_minus_prefit ) ,
                     (ratio_plus_postfit , ratio_minus_postfit)]
    else :
      ratio_hists = [(ratio_plus_postfit , ratio_minus_postfit)]

    # ratio_graph_prefit  = ROOT.TGraphAsymmErrors(ratio_plus_prefit  .GetXaxis().GetNbins())
    # ratio_graph_postfit = ROOT.TGraphAsymmErrors(ratio_plus_postfit .GetXaxis().GetNbins())
    
    # ratio_graphs = [ratio_graph_prefit  ,
    #                 ratio_graph_postfit , ]
    
    # for (hp,hm),g in zip(ratio_hists,ratio_graphs):
    #   for i in xrange(h.GetXaxis().GetNbins()):
    #     low_edge  = ratio_plus_prefit.GetBinLowEdge(i+1)
    #     high_edge = low_edge + ratio_plus_prefit.GetBinWidth(i+1)
    #     g.SetPoint      (i,ratio_plus_prefit.GetBinCenter(i+1),1)
    #     g.SetPointEYlow (i,abs(hm.GetBinContent(i+1)-1))
    #     g.SetPointEYhigh(i,abs(hp.GetBinContent(i+1)-1))
    #     g.SetPointEXlow (i,ratio_plus_prefit.GetBinWidth(i+1)*0.5 - 0.005)
    #     g.SetPointEXhigh(i,ratio_plus_prefit.GetBinWidth(i+1)*0.5 - 0.005)
        
    #   g.SetFillStyle(1001)
    #   g.SetLineStyle(0)
    #   g.SetLineWidth(0)
    #   g.GetXaxis().SetRangeUser(0,max_energy)
    #   CopyHistAxisFormatting(hplus_prefit,g)
    #   g.GetYaxis().SetTitle("Ratio of # Total Events")
      
    # ratio_graph_prefit .SetFillColor(R.kBlue)
    # ratio_graph_postfit.SetFillColor(R.kYellow)
    # ratio_graph_prefit .SetLineColor(R.kBlue)
    # ratio_graph_postfit.SetLineColor(R.kYellow)    
    # #ratio_graph_prefit .SetLineWidth(3)
    # #ratio_graph_postfit.SetLineWidth(3)
    
    # #ratio_graph_prefit.SetFillStyle(1001)
    # ratio_graph_prefit.SetFillStyle(3001)
    # ratio_graph_prefit.SetLineColor(R.kBlue)
    # ratio_graph_prefit.SetLineWidth(2)
        
    # #ratio_graph_postfit.SetFillStyle(1001)
    # ratio_graph_postfit.SetFillStyle(3001)
    # ratio_graph_postfit.SetLineWidth(2)

    canv2 = ROOT.TCanvas()
    canv2.SetLeftMargin  (left_margin)
    canv2.SetBottomMargin(bottom_margin)
    canv2.SetRightMargin (right_margin )
    if use_log_plots: canv2.SetLogx()
    canv2.cd()
    
    # ratio_max = 1.3*ratio_plus_prefit .GetMaximum()
    # ratio_graph_postfit.SetMaximum(ratio_max)
    # ratio_graph_postfit.SetMinimum(2-ratio_max)
    # ratio_graph_prefit .Draw("a2")
    # ratio_graph_postfit.Draw("2")

    ratio_plus_postfit.SetFillStyle(0)
    ratio_minus_postfit.SetFillStyle(0)
 
    ratio_max = 1.5
    ratio_min = 0.5
    
    ratio_plus_postfit.SetLineWidth(2)
    ratio_minus_postfit.SetLineWidth(2)
    if doprefit==True :
      ratio_plus_prefit.SetLineWidth(2)
      ratio_minus_prefit.SetLineWidth(2)

    ratio_plus_postfit.SetLineStyle(1)
    ratio_minus_postfit.SetLineStyle(1)
    if doprefit==True :
      ratio_plus_prefit.SetLineStyle(1)
      ratio_minus_prefit.SetLineStyle(1)

    ratio_plus_postfit.SetLineColor(R.kRed)
    ratio_minus_postfit.SetLineColor(R.kRed)
    if doprefit==True :
      ratio_plus_prefit.SetLineColor(R.kBlue)
      ratio_minus_prefit.SetLineColor(R.kBlue)
    
    ratio_plus_postfit.SetMaximum(ratio_max)
    ratio_plus_postfit.SetMinimum(ratio_min)
    ratio_minus_postfit.SetMaximum(ratio_max)
    ratio_minus_postfit.SetMinimum(ratio_min)
    if doprefit==True :
      ratio_plus_prefit.SetMaximum(ratio_max)
      ratio_plus_prefit.SetMinimum(ratio_min)
      ratio_minus_prefit.SetMaximum(ratio_max)
      ratio_minus_prefit.SetMinimum(ratio_min)
 
      
    if doprefit==True :
      ratio_plus_prefit.GetYaxis().SetTitle("Fractional error")
      ratio_plus_prefit.DrawCopy("")
      ratio_minus_prefit.DrawCopy("same")
      ratio_plus_postfit.DrawCopy("same")
      ratio_minus_postfit.DrawCopy("same")
    else :
      ratio_plus_postfit.GetYaxis().SetTitle("Fractional error")
      ratio_plus_postfit.DrawCopy("")
      ratio_minus_postfit.DrawCopy("same")

    ROOT.gPad.RedrawAxis();

    legend2 = ROOT.TLegend(0.4,0.70,0.73,0.87);
    legend2.SetTextSize(0.05);
    legend2.SetFillStyle(0);
    legend2.SetBorderSize(0);
    if doprefit==True :
      legend2.AddEntry(ratio_plus_prefit  ,"Nominal error" ,"L")
    legend2.AddEntry(ratio_plus_postfit ,"#pm 1#sigma error","L")
    legend2.Draw("same")
    
    line_at_1 = ROOT.TLine(0,1,max_energy,1)
    line_at_1.SetLineWidth(2)
    line_at_1.Draw("same")
    
    canv2.Update()
    canv2.Print(tag+"_"+sample+name+"_"+par_set+"_envelope_ratio.eps")

    # Print out some event rates
    nom_events     = nominal_hist  .Integral()
    if doprefit==True :
      p1_events_pre  = hplus_prefit  .Integral()
      m1_events_pre  = hminus_prefit .Integral()
    p1_events_post = hplus_postfit .Integral()
    m1_events_post = hminus_postfit.Integral()
    
    if doprefit==True :    
      p1_percent_pre  = 100.*( (p1_events_pre -nom_events)/nom_events )
      m1_percent_pre  = 100.*( (m1_events_pre -nom_events)/nom_events )
    p1_percent_post = 100.*( (p1_events_post-nom_events)/nom_events )
    m1_percent_post = 100.*( (m1_events_post-nom_events)/nom_events )
    
    items = [name+"_"+par_set,to3dp.format(nom_events)]
    if doprefit==True :
      items.append(to3dp.format(p1_events_pre ) +" ("+to3dp.format(p1_percent_pre )+"%)")
      items.append(to3dp.format(m1_events_pre ) +" ("+to3dp.format(m1_percent_pre )+"%)")
    items.append(to3dp.format(p1_events_post) +" ("+to3dp.format(p1_percent_post)+"%)")
    items.append(to3dp.format(m1_events_post) +" ("+to3dp.format(m1_percent_post)+"%)")
    print "\t".join(items)
    

  canv_nsk_prefit.Update()
  canv_nsk_prefit.Print(tag+"_"+sample+"_prefit_nsk.eps")
  canv_nsk_postfit.Update()
  canv_nsk_postfit.Print(tag+"_"+sample+"_nsk.eps")

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
  parser.add_option('-s', dest='sample',default="RHC1Re"  ,help='The sample')
  parser.add_option('-o', dest='osc'  ,default=default_oscillated_input_file  ,help='The unoscillated input ROOT file')
  parser.add_option('-u', dest='unosc',default=default_unoscillated_input_file,help='The oscillated input ROOT file')
  parser.add_option('-a', dest='usemedian', action='store_true' , default=False, help='Use envelopes with respect to the median or average histogram')
  parser.add_option('-p', dest='doprefit'  , action='store_true' , default=False, help='Draw pre-banff histograms')
  parser.add_option('-t', dest='tag',default="test"  ,help='Tag of the output file')

  options,other = parser.parse_args()
  
  if other:
    print "Bad arguments: ",other
    exit(-1)
  
  main(options.sample,options.unosc,options.osc,options.usemedian,options.doprefit,options.tag)
