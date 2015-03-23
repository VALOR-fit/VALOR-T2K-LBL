#!/usr/bin/python

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
conf68 = 2.71

to3dp = "{0:.3f}"
  
def main(marg_file, tag=None):
  
  ROOT.gStyle.SetPalette(56);
  
  mf = ROOT.TFile(marg_file)
  h_m = mf.Get("h_chisq_dmsq_mixing23_Ave_nh")
  assert(h_m)
  h_m.Print()
  min_m = h_m.GetMinimum()
  print min_m
  
  h_delta_m = ROOT.TH2D(h_m)
  h_delta_m.Clear()
  for ix in xrange(h_delta_m.GetXaxis().GetNbins()):
    for iy in xrange(h_delta_m.GetYaxis().GetNbins()):
      h_delta_m.SetBinContent(ix,iy,h_m.GetBinContent(ix,iy)-min_m)
  
  h_contour_m = ROOT.TH2D(h_delta_m)
  h_contour_m.SetContour(2)
  h_contour_m.SetContourLevel(1,conf90)
  h_contour_m.SetContourLevel(0,conf68)
  h_contour_m.SetLineWidth(1)
  h_contour_m.SetLineColor(R.kBlack)
  
  
  h_p = mf.Get("h_chisq_dmsq_mixing23_Ave_ih")
  assert(h_p)
  h_p.Print()
  min_p = h_p.GetMinimum()
  print min_p
  
  h_delta_p = ROOT.TH2D(h_p)
  h_delta_p.Clear()
  for ix in xrange(h_delta_p.GetXaxis().GetNbins()):
    for iy in xrange(h_delta_p.GetYaxis().GetNbins()):
      h_delta_p.SetBinContent(ix,iy,h_p.GetBinContent(ix,iy)-min_p)
  
  h_contour_p = ROOT.TH2D(h_delta_p)
  h_contour_p.SetContour(2)
  h_contour_p.SetContourLevel(1,conf90)
  h_contour_p.SetContourLevel(0,conf68)
  h_contour_p.SetLineWidth(1)
  h_contour_p.SetLineStyle(2)
  h_contour_p.SetLineColor(R.kBlack)
  
  contours_p = GetContours(h_contour_p)
  contours_m = GetContours(h_contour_m)
  
  for g in contours_p[0]: g.SetLineColor(R.kBlack)
  for g in contours_p[1]: g.SetLineColor(R.kRed)
  for g in contours_p[0]: g.SetLineStyle(1)
  for g in contours_p[1]: g.SetLineStyle(1)
  
  for g in contours_m[0]: g.SetLineColor(R.kBlack)
  for g in contours_m[1]: g.SetLineColor(R.kRed)
  for g in contours_m[0]: g.SetLineStyle(2)
  for g in contours_m[1]: g.SetLineStyle(2)
  
  # Pick bestfit points
  bf_bin_p = h_contour_p.GetMinimumBin()
  bf_bin_m = h_contour_m.GetMinimumBin()
  
  bf_chisq_p = h_contour_p.GetBinContent(bf_bin_p)
  bf_chisq_m = h_contour_p.GetBinContent(bf_bin_m)
  
  x_bin_p = ROOT.Long() ; y_bin_p = ROOT.Long()
  x_bin_m = ROOT.Long() ; y_bin_m = ROOT.Long()
  dummy = ROOT.Long()
  
  h_contour_p.GetBinXYZ(bf_bin_p,x_bin_p,y_bin_p,dummy)
  h_contour_m.GetBinXYZ(bf_bin_m,x_bin_m,y_bin_m,dummy)
  
  bf_x_p = h_contour_p.GetXaxis().GetBinCenter(x_bin_p)
  bf_x_m = h_contour_m.GetXaxis().GetBinCenter(x_bin_m)
  bf_y_p = h_contour_p.GetYaxis().GetBinCenter(y_bin_p)
  bf_y_m = h_contour_m.GetYaxis().GetBinCenter(y_bin_m)
  
  marker_p = ROOT.TMarker(bf_x_p,bf_y_p,ROOT.kFullStar)
  marker_m = ROOT.TMarker(bf_x_m,bf_y_m,ROOT.kFullStar)
  
  if (bf_chisq_p > bf_chisq_m):
    marker_m.SetMarkerStyle(ROOT.kOpenStar)
  else:
    marker_p.SetMarkerStyle(ROOT.kOpenStar)
    
  marker_p.SetMarkerColor(ROOT.kViolet);
  marker_m.SetMarkerColor(ROOT.kBlue+1);
  marker_p.SetMarkerSize(1.5);
  marker_m.SetMarkerSize(1.5);
  
  CopyMarkerToHist(marker_p,h_contour_p)
  CopyMarkerToHist(marker_m,h_contour_m)
  
  # Start drawing
  canv = ROOT.TCanvas("canv1","canv1",1280,720)
  canv.cd()
  
  h_contour_m.Draw("AXIS")
  for cont_list in contours_p + contours_m:
    for cont in cont_list:
      print cont.GetName()
      cont.Draw("same")
      
  # Draw markers
  marker_p.Draw()
  marker_m.Draw()
  
  # Make legend
  
  legend = ROOT.TLegend(0.7,0.7,0.9,0.9)
  legend.AddEntry(h_contour_m,"NH","lp")
  legend.AddEntry(h_contour_p,"IH","lp")
  legend.Draw()
  
  if tag is None:
    plot_name = "profiled_nh_vs_ih"
  else:
    plot_name = tag
  
  canv.Print(plot_name+"."+image_format)
  
  canv_nh = ROOT.TCanvas("canv_nh","canv_nh",1280,720)
  canv_nh.cd()
  h_delta_m.Draw("colz")
  marker_m.Draw()
  canv_nh.Print(plot_name+"_nh."+image_format)
  
  canv_ih = ROOT.TCanvas("canv_ih","canv_ih",1280,720)
  canv_ih.cd()
  h_delta_p.Draw("colz")
  marker_p.Draw()
  canv_ih.Print(plot_name+"_ih."+image_format)

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
  
if __name__=="__main__":
  main(*sys.argv[1:])

