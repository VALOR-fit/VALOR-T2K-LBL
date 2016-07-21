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

# First file is asim-only, second is asim+data
def main(asim_file, data_file, asim_bf_normal = None, data_bf_normal =None, asim_bf_inverted = None, data_bf_inverted =None, tag=None):
  
  ROOT.gStyle.SetPalette(56);
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  
  # m = asim only
  # p = asim+data
  
  assert(data_file)
  
  print "Reading asim from ",asim_file
  print "Reading data from ",data_file
  
  for hierarchy in ["nh","ih"]:
    m_f = ROOT.TFile(asim_file)
    h_m = m_f.Get("h_chisq_dmsq_mixing23_Ave_"+hierarchy)
    assert(h_m)
    h_m.Print()
    min_m = h_m.GetMinimum()
    
    m_p = ROOT.TFile(data_file)
      
    h_p = m_p.Get("h_chisq_dmsq_mixing23_Ave_"+hierarchy)
    assert(h_p)
    h_p.Print()
    min_p = h_p.GetMinimum()
    
    print min_m
    print min_p
    
    # Precalculate the best-fit points from the histogram
    bf_bin_p = h_p.GetMinimumBin()
    bf_bin_m = h_m.GetMinimumBin()
    
    bf_chisq_p = h_p.GetBinContent(bf_bin_p)
    bf_chisq_m = h_m.GetBinContent(bf_bin_m)
    
    x_bin_p = ROOT.Long() ; y_bin_p = ROOT.Long()
    x_bin_m = ROOT.Long() ; y_bin_m = ROOT.Long()
    dummy = ROOT.Long()
    
    h_p.GetBinXYZ(bf_bin_p,x_bin_p,y_bin_p,dummy)
    h_m.GetBinXYZ(bf_bin_m,x_bin_m,y_bin_m,dummy)
    
    bf_x_p = h_p.GetXaxis().GetBinCenter(x_bin_p)
    bf_x_m = h_m.GetXaxis().GetBinCenter(x_bin_m)
    bf_y_p = h_p.GetYaxis().GetBinCenter(y_bin_p)
    bf_y_m = h_m.GetYaxis().GetBinCenter(y_bin_m)
    
    
    # If we have files, go and get the bestfits from those instead.
    if (hierarchy=="nh"):
      if asim_bf_normal is not None:
        tf_asim_normal = ROOT.TFile(asim_bf_normal)
        tr = tf_asim_normal.Get("fitnt")
        tr.GetEntry(0)
        bf_x_m = tr.best_fit_params.BestFitValue(p1)
        bf_y_m = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_m = tr.best_fit_params.GetBestFitStatistic()
      if data_bf_normal is not None:
        tf_data_normal = ROOT.TFile(data_bf_normal)
        tr = tf_data_normal.Get("fitnt")
        tr.GetEntry(0)
        bf_x_p = tr.best_fit_params.BestFitValue(p1)
        bf_y_p = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_p = tr.best_fit_params.GetBestFitStatistic()
    else:
      if asim_bf_inverted is not None:
        tf_asim_inverted = ROOT.TFile(asim_bf_inverted)
        tr = tf_asim_inverted.Get("fitnt")
        tr.GetEntry(0)
        bf_x_m = tr.best_fit_params.BestFitValue(p1)
        bf_y_m = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_m = tr.best_fit_params.GetBestFitStatistic()
      if data_bf_inverted is not None:
        tf_data_inverted = ROOT.TFile(data_bf_inverted)
        tr = tf_data_inverted.Get("fitnt")
        tr.GetEntry(0)
        bf_x_p = tr.best_fit_params.BestFitValue(p1)
        bf_y_p = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_p = tr.best_fit_params.GetBestFitStatistic()
        
    
    h_delta_m = GetDeltaHist(h_m,bf_chisq_m)
    
    h_contour_m = ROOT.TH2D(h_delta_m)
    h_contour_m.SetContour(2)
    h_contour_m.SetContourLevel(1,conf90)
    h_contour_m.SetContourLevel(0,conf68)
    h_contour_m.SetLineWidth(1)
    h_contour_m.SetLineColor(R.kBlack)
    
    h_delta_p = GetDeltaHist(h_p,bf_chisq_p)
    
    h_contour_p = ROOT.TH2D(h_delta_p)
    h_contour_p.SetContour(2)
    h_contour_p.SetContourLevel(1,conf90)
    h_contour_p.SetContourLevel(0,conf68)
    h_contour_p.SetLineWidth(1)
    h_contour_p.SetLineStyle(2)
    h_contour_p.SetLineColor(R.kBlack)
    
    contours_p = GetContours(h_contour_p)
    contours_m = GetContours(h_contour_m)
    
    marker_p = ROOT.TMarker(bf_x_p,bf_y_p,ROOT.kFullStar)
    marker_m = ROOT.TMarker(bf_x_m,bf_y_m,ROOT.kFullStar)
    
    if (bf_chisq_p > bf_chisq_m):
      bf_is_normal = True
      marker_m.SetMarkerStyle(ROOT.kOpenStar)
    else:
      bf_is_normal = False
      marker_p.SetMarkerStyle(ROOT.kOpenStar)
      
    marker_p.SetMarkerColor(ROOT.kViolet);
    marker_m.SetMarkerColor(ROOT.kBlue+1);
    marker_p.SetMarkerSize(1.5);
    marker_m.SetMarkerSize(1.5);
    
    CopyMarkerToHist(marker_p,h_contour_p)
    CopyMarkerToHist(marker_m,h_contour_m)
    
    print "Bestfit point (asim)  : ",bf_x_m,bf_y_m,"\twith chisq: ",bf_chisq_m
    print "Bestfit point (data): ",bf_x_p,bf_y_p,"\twith chisq: ",bf_chisq_p
      
    print "\n\n"
    
    # Start drawing
    for zoom in [False,True]:
      if zoom: h_contour_m.GetYaxis().SetRangeUser(1E-3,4.5E-3)
      
      for g in contours_p[0]: g.SetLineColor(R.kBlack)
      for g in contours_p[1]: g.SetLineColor(R.kRed)
      for g in contours_p[0]: g.SetLineStyle(2)
      for g in contours_p[1]: g.SetLineStyle(2)
      
      for g in contours_m[0]: g.SetLineColor(R.kBlack)
      for g in contours_m[1]: g.SetLineColor(R.kRed)
      for g in contours_m[0]: g.SetLineStyle(1)
      for g in contours_m[1]: g.SetLineStyle(1)
    
      
      canv = ROOT.TCanvas("canv1","canv1",1280,720)
      canv.cd()
      
      h_contour_m.Draw("AXIS")
      for cont_list in contours_p + contours_m:
        for cont in cont_list:
          cont.Draw("same")
          
      # Draw markers
      marker_p.Draw()
      marker_m.Draw()
      
      # Make legend
      
      legend = ROOT.TLegend(0.7,0.7,0.9,0.9)
      legend.AddEntry(h_contour_m,"Asimov dataset","lp")
      legend.AddEntry(h_contour_p,"Real data","lp")
      legend.Draw()
      
      if tag is None:
        plot_name = "profiled_asim_data"
      else:
        plot_name = tag
        
      if zoom: plot_name += "_zoomed"
      
      canv.Print(plot_name+"_"+hierarchy+"."+image_format)
      
      canv_nh = ROOT.TCanvas("canv_nh","canv_nh",1280,720)
      canv_nh.cd()
      h_delta_m.SetContour(100)
      h_delta_m.Draw("colz")
      marker_m.Draw()
      for g in contours_m[1]:
        g.SetLineStyle(1)
        g.SetLineColor(R.kBlue)
        g.Draw("same")
      canv_nh.Print(plot_name+"_"+hierarchy+"_asimonly."+image_format)
      
      canv_ih = ROOT.TCanvas("canv_ih","canv_ih",1280,720)
      canv_ih.cd()
      h_delta_p.SetContour(100)
      h_delta_p.Draw("colz")
      marker_p.Draw()
      for g in contours_p[1]:
        g.SetLineStyle(1)
        g.SetLineColor(R.kBlue)
        g.Draw("same")
      canv_ih.Print(plot_name+"_"+hierarchy+"_data."+image_format)

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
  
  parser.add_option('--asim-bf-nh' , dest='asim_bf_normal' ,  default=None, help="File containing true best-fit (asim, NH)")
  parser.add_option('--asim-bf-ih' , dest='asim_bf_inverted', default=None, help="File containing true best-fit (asim, IH)")
  parser.add_option('--data-bf-nh' , dest='data_bf_normal' ,  default=None, help="File containing true best-fit (data, NH)")
  parser.add_option('--data-bf-ih' , dest='data_bf_inverted', default=None, help="File containing true best-fit (data, IH)")
  parser.add_option('-t','--tag' , dest='tag', default=None, help="Tag for output")
  
  o,args = parser.parse_args()
  
  main(args[0],args[1],o.asim_bf_normal,o.data_bf_normal,o.asim_bf_inverted,o.data_bf_inverted,o.tag)

