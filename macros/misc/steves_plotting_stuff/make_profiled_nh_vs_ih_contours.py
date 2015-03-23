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
  h_contour_m.SetLineColor(R.kBlack)
  
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
  
  # Start drawing
  canv = ROOT.TCanvas("canv1","canv1",1280,720)
  canv.cd()
  
  h_contour_m.Draw("AXIS")
  for cont_list in contours_p + contours_m:
    for cont in cont_list:
      cont.Draw("same")
  
  # Make legend
  
  
  legend = ROOT.TLegend(0.7,0.7,0.9,0.9)
  legend.AddEntry(h_contour_m,"NH")
  legend.AddEntry(h_contour_p,"IH")
  legend.Draw()
  
  if tag is None:
    plot_name = "profiled_nh_vs_ih"
  else:
    plot_name = tag
  
  canv.Print(plot_name+"."+image_format)

def GetContours(h):
  dummy_canv = ROOT.TCanvas()
  dummy_canv.cd()
  h.Draw("cont list")
  dummy_canv.Update()
  obj_array = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
  del dummy_canv
  tlists = [obj_array[i] for i in xrange(obj_array.GetEntries())]
  l = []
  for tl in tlists:
    l.append([ROOT.TGraph(tl[i]) for i in xrange(tl.GetEntries())])
  return l
  
if __name__=="__main__":
  main(*sys.argv[1:])

