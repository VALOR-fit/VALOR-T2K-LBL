#!/usr/bin/env python

# Plots contours
# Steve Dennis, March 2015

import ROOT
import math
import array
import numpy as np
import sys
import os

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
R=ROOT
  
ROOT.gROOT.SetBatch(1)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

opt.gOscMixingParamChoice = conv.kNuMixSinsqTheta
opt.SetDefaultParameterValues()

image_format = "eps"
conf90 = 4.61
conf68 = 2.30

to3dp = "{0:.3f}"
to6dp = "{0:.6f}"

p1 = conv.kPhypMixing23Bar
p2 = conv.kPhypDm23sqBar

# First file is stat-only, second is stat+syst
def main(stat_file, syst_file, stat_bf_normal = None, syst_bf_normal =None, stat_bf_inverted = None, syst_bf_inverted =None, tag=None):
  
  ROOT.gStyle.SetPalette(56);
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  
  # m = stat only
  # p = stat+syst
  
  assert(syst_file)
  
  print "Reading stat from ",stat_file
  print "Reading stat+SYST from ",syst_file
  
  for hierarchy in ["nh","ih"]:
    f_stat = ROOT.TFile(stat_file)
    h_stat = f_stat.Get("h_chisq_dmsq_mixing23_Ave_"+hierarchy)
    assert(h_stat)
    h_stat.Print()
    min_stat = h_stat.GetMinimum()
    
    f_syst = ROOT.TFile(syst_file)
      
    h_syst = f_syst.Get("h_chisq_dmsq_mixing23_Ave_"+hierarchy)
    assert(h_syst)
    h_syst.Print()
    min_syst = h_syst.GetMinimum()
    
    print min_stat
    print min_syst
    
    # Precalculate the best-fit points from the histogram
    bf_bin_syst = h_syst.GetMinimumBin()
    bf_bin_stat = h_stat.GetMinimumBin()
    
    bf_chisq_syst = h_syst.GetBinContent(bf_bin_syst)
    bf_chisq_stat = h_stat.GetBinContent(bf_bin_stat)
    
    x_bin_syst = ROOT.Long() ; y_bin_syst = ROOT.Long()
    x_bin_stat = ROOT.Long() ; y_bin_stat = ROOT.Long()
    dummy = ROOT.Long()
    
    h_syst.GetBinXYZ(bf_bin_syst,x_bin_syst,y_bin_syst,dummy)
    h_stat.GetBinXYZ(bf_bin_stat,x_bin_stat,y_bin_stat,dummy)
    
    bf_x_syst = h_syst.GetXaxis().GetBinCenter(x_bin_syst)
    bf_x_stat = h_stat.GetXaxis().GetBinCenter(x_bin_stat)
    bf_y_syst = h_syst.GetYaxis().GetBinCenter(y_bin_syst)
    bf_y_stat = h_stat.GetYaxis().GetBinCenter(y_bin_stat)
    
    
    # If we have files, go and get the bestfits from those instead.
    if (hierarchy=="nh"):
      if stat_bf_normal is not None:
        tf_stat_normal = ROOT.TFile(stat_bf_normal)
        tr = tf_stat_normal.Get("fitnt")
        tr.GetEntry(0)
        bf_x_stat = tr.best_fit_params.BestFitValue(p1)
        bf_y_stat = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_stat = tr.best_fit_params.GetBestFitStatistic()
      if syst_bf_normal is not None:
        tf_syst_normal = ROOT.TFile(syst_bf_normal)
        tr = tf_syst_normal.Get("fitnt")
        tr.GetEntry(0)
        bf_x_syst = tr.best_fit_params.BestFitValue(p1)
        bf_y_syst = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_syst = tr.best_fit_params.GetBestFitStatistic()
    else:
      if stat_bf_inverted is not None:
        tf_stat_inverted = ROOT.TFile(stat_bf_inverted)
        tr = tf_stat_inverted.Get("fitnt")
        tr.GetEntry(0)
        bf_x_stat = tr.best_fit_params.BestFitValue(p1)
        bf_y_stat = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_stat = tr.best_fit_params.GetBestFitStatistic()
      if syst_bf_inverted is not None:
        tf_syst_inverted = ROOT.TFile(syst_bf_inverted)
        tr = tf_syst_inverted.Get("fitnt")
        tr.GetEntry(0)
        bf_x_syst = tr.best_fit_params.BestFitValue(p1)
        bf_y_syst = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_syst = tr.best_fit_params.GetBestFitStatistic()
        
    
    h_delta_stat = GetDeltaHist(h_stat,bf_chisq_stat)
    
    h_contour_stat = ROOT.TH2D(h_delta_stat)
    h_contour_stat.SetContour(2)
    h_contour_stat.SetContourLevel(1,conf90)
    h_contour_stat.SetContourLevel(0,conf68)
    h_contour_stat.SetLineWidth(1)
    h_contour_stat.SetLineColor(R.kBlack)
    
    h_delta_syst = GetDeltaHist(h_syst,bf_chisq_syst)
    
    h_contour_syst = ROOT.TH2D(h_delta_syst)
    h_contour_syst.SetContour(2)
    h_contour_syst.SetContourLevel(1,conf90)
    h_contour_syst.SetContourLevel(0,conf68)
    h_contour_syst.SetLineWidth(1)
    h_contour_syst.SetLineStyle(2)
    h_contour_syst.SetLineColor(R.kBlack)
    
    contours_syst = GetContours(h_contour_syst)
    contours_stat = GetContours(h_contour_stat)
    
    marker_syst = ROOT.TMarker(bf_x_syst,bf_y_syst,ROOT.kFullStar)
    marker_stat = ROOT.TMarker(bf_x_stat,bf_y_stat,ROOT.kFullStar)
    
    if (bf_chisq_syst > bf_chisq_stat):
      bf_is_normal = True
      marker_stat.SetMarkerStyle(ROOT.kOpenStar)
    else:
      bf_is_normal = False
      marker_syst.SetMarkerStyle(ROOT.kOpenStar)
      
    marker_syst.SetMarkerColor(ROOT.kViolet);
    marker_stat.SetMarkerColor(ROOT.kBlue+1);
    marker_syst.SetMarkerSize(1.5);
    marker_stat.SetMarkerSize(1.5);
    
    CopyMarkerToHist(marker_syst,h_contour_syst)
    CopyMarkerToHist(marker_stat,h_contour_stat)
    
    print "Bestfit point (stat)  : "   ,bf_x_stat,bf_y_stat,"\twith chisq: ",bf_chisq_stat
    print "Bestfit point (stat+syst): ",bf_x_syst,bf_y_syst,"\twith chisq: ",bf_chisq_syst
      
    print "\n\n"
    
    # Start drawing
    
    #~ for zoom in [False,True]:
    for zoom in [True]:
      if zoom: h_contour_stat.GetYaxis().SetRangeUser(1E-3,4.5E-3)
      
      for g in contours_syst[0]: g.SetLineColor(R.kBlack)
      for g in contours_syst[1]: g.SetLineColor(R.kRed)
      for g in contours_syst[0]: g.SetLineStyle(h_contour_syst.GetLineStyle())
      for g in contours_syst[1]: g.SetLineStyle(h_contour_syst.GetLineStyle())
      
      for g in contours_stat[0]: g.SetLineColor(R.kBlack)
      for g in contours_stat[1]: g.SetLineColor(R.kRed)
      for g in contours_stat[0]: g.SetLineStyle(h_contour_stat.GetLineStyle())
      for g in contours_stat[1]: g.SetLineStyle(h_contour_stat.GetLineStyle())
    
      
      canv = ROOT.TCanvas("canv1","canv1",1280,720)
      canv.cd()
      
      h_contour_stat.Draw("AXIS")
      for cont_list in contours_syst + contours_stat:
        for cont in cont_list:
          cont.Draw("same")
          
      # Draw markers
      marker_syst.Draw()
      marker_stat.Draw()
      
      # Make legend
      
      legend = ROOT.TLegend(0.7,0.7,0.9,0.9)
      legend.AddEntry(h_contour_stat,"Stat only","lp")
      legend.AddEntry(h_contour_syst,"Stat+syst","lp")
      legend.Draw()
      
      if tag is None:
        plot_name = "profiled_stat_syst"
      else:
        plot_name = tag
        
      if zoom: plot_name += "_zoomed"
      
      canv.Print(plot_name+"_"+hierarchy+"."+image_format)
      
      canv_nh = ROOT.TCanvas("canv_nh","canv_nh",1280,720)
      canv_nh.cd()
      h_delta_stat.SetContour(100)
      h_delta_stat.Draw("colz")
      marker_stat.Draw()
      for g in contours_stat[1]:
        g.SetLineStyle(1)
        g.SetLineColor(R.kBlue)
        g.Draw("same")
      canv_nh.Print(plot_name+"_"+hierarchy+"_statonly."+image_format)
      
      canv_ih = ROOT.TCanvas("canv_ih","canv_ih",1280,720)
      canv_ih.cd()
      h_delta_syst.SetContour(100)
      h_delta_syst.Draw("colz")
      marker_syst.Draw()
      for g in contours_syst[1]:
        g.SetLineStyle(1)
        g.SetLineColor(R.kBlue)
        g.Draw("same")
      canv_ih.Print(plot_name+"_"+hierarchy+"_stat+syst."+image_format)

def GetContours(h):
  dummy_canv = ROOT.TCanvas()
  dummy_canv.cd()
  h.Draw("cont list")
  dummy_canv.Update()
  obj_array = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
  del dummy_canv
  tlists         = [obj_array[i]         for i in xrange(obj_array.GetEntries())]
  contour_levels = [h.GetContourLevel(i) for i in xrange(obj_array.GetEntries())]
  to_return = []
  for tl,level in zip(tlists,contour_levels):
    cont_list = [ROOT.TGraph(tl[i]) for i in xrange(tl.GetEntries())]
    for c in cont_list:
      c.SetName ("contour_"+to3dp.format(level))
      c.SetTitle("Contour "+to3dp.format(level))
    to_return.append(cont_list)
  return to_return
  
def CopyMarkerToHist(marker,histogram):
  histogram.SetMarkerSize (marker.GetMarkerSize ())
  histogram.SetMarkerStyle(marker.GetMarkerStyle())
  histogram.SetMarkerColor(marker.GetMarkerColor())
  
def GetDeltaHist(h,minval=None):
  h_delta = ROOT.TH2D(h)
  h_delta.Reset()
  if minval is None: minval = h.GetMinimum()
  for ix in xrange(h_delta.GetXaxis().GetNbins()):
    for iy in xrange(h_delta.GetYaxis().GetNbins()):
      h_delta.SetBinContent(ix+1,iy+1,h.GetBinContent(ix+1,iy+1)-minval)
      
  return h_delta
  
if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make contours")
  
  parser.add_option('--stat-bf-nh' , dest='stat_bf_normal' ,  default=None, help="File containing true best-fit (stat, NH)")
  parser.add_option('--stat-bf-ih' , dest='stat_bf_inverted', default=None, help="File containing true best-fit (stat, IH)")
  parser.add_option('--syst-bf-nh' , dest='syst_bf_normal' ,  default=None, help="File containing true best-fit (syst, NH)")
  parser.add_option('--syst-bf-ih' , dest='syst_bf_inverted', default=None, help="File containing true best-fit (syst, IH)")
  parser.add_option('-t','--tag' , dest='tag', default=None, help="Tag for output")
  
  o,args = parser.parse_args()
  
  main(args[0],args[1],o.stat_bf_normal,o.syst_bf_normal,o.stat_bf_inverted,o.syst_bf_inverted,o.tag)

