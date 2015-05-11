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
conf68 = 2.30

to3dp = "{0:.3f}"
to6dp = "{0:.6f}"

p1 = conv.kPhypMixing23Bar
p2 = conv.kPhypDm23sqBar

# First file is asim-only, second is asim+data
def main(asim_file, data_file, asim_bf_normal = None, data_bf_normal =None, asim_bf_inverted = None, data_bf_inverted =None,
         tag=None, draw_68=True):
  
  ROOT.gStyle.SetPalette(56);
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  
  title_size = 0.08
  
  line_style_68 = ROOT.kDashed
  line_style_90 = ROOT.kSolid
  
  # m = asim only
  # p = asim+data
  
  assert(data_file)
  
  print "Reading asim from ",asim_file
  print "Reading data from ",data_file
  
  for hierarchy in ["nh","ih"]:
    f_asim = ROOT.TFile(asim_file)
    h_asim = f_asim.Get("h_chisq_dmsq_mixing23_Ave_"+hierarchy)
    assert(h_asim)
    h_asim.Print()
    min_asim = h_asim.GetMinimum()
    
    f_data = ROOT.TFile(data_file)
      
    h_data = f_data.Get("h_chisq_dmsq_mixing23_Ave_"+hierarchy)
    assert(h_data)
    h_data.Print()
    min_data = h_data.GetMinimum()
    
    print min_asim
    print min_data
    
    # Set title offset
    h_data.GetXaxis().SetTitleOffset(0.8)
    h_asim.GetXaxis().SetTitleOffset(0.8)
    h_asim.GetYaxis().SetTitleOffset(0.9)
    h_data.GetYaxis().SetTitleOffset(0.9)
    
    h_asim.GetXaxis().SetTitleSize(title_size)
    h_data.GetXaxis().SetTitleSize(title_size)
    h_asim.GetYaxis().SetTitleSize(title_size)
    h_data.GetYaxis().SetTitleSize(title_size)
    
    # Precalculate the best-fit points from the histogram
    bf_bin_data = h_data.GetMinimumBin()
    bf_bin_asim = h_asim.GetMinimumBin()
    
    bf_chisq_data = h_data.GetBinContent(bf_bin_data)
    bf_chisq_asim = h_asim.GetBinContent(bf_bin_asim)
    
    x_bin_data = ROOT.Long() ; y_bin_data = ROOT.Long()
    x_bin_asim = ROOT.Long() ; y_bin_asim = ROOT.Long()
    dummy = ROOT.Long()
    
    h_data.GetBinXYZ(bf_bin_data,x_bin_data,y_bin_data,dummy)
    h_asim.GetBinXYZ(bf_bin_asim,x_bin_asim,y_bin_asim,dummy)
    
    bf_x_data = h_data.GetXaxis().GetBinCenter(x_bin_data)
    bf_x_asim = h_asim.GetXaxis().GetBinCenter(x_bin_asim)
    bf_y_data = h_data.GetYaxis().GetBinCenter(y_bin_data)
    bf_y_asim = h_asim.GetYaxis().GetBinCenter(y_bin_asim)
    
    
    # If we have files, go and get the bestfits from those instead.
    if (hierarchy=="nh"):
      if asim_bf_normal is not None:
        tf_asim_normal = ROOT.TFile(asim_bf_normal)
        tr = tf_asim_normal.Get("fitnt")
        tr.GetEntry(0)
        bf_x_asim = tr.best_fit_params.BestFitValue(p1)
        bf_y_asim = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_asim = tr.best_fit_params.GetBestFitStatistic()
      if data_bf_normal is not None:
        tf_data_normal = ROOT.TFile(data_bf_normal)
        tr = tf_data_normal.Get("fitnt")
        tr.GetEntry(0)
        bf_x_data = tr.best_fit_params.BestFitValue(p1)
        bf_y_data = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_data = tr.best_fit_params.GetBestFitStatistic()
    else:
      if asim_bf_inverted is not None:
        tf_asim_inverted = ROOT.TFile(asim_bf_inverted)
        tr = tf_asim_inverted.Get("fitnt")
        tr.GetEntry(0)
        bf_x_asim = tr.best_fit_params.BestFitValue(p1)
        bf_y_asim = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_asim = tr.best_fit_params.GetBestFitStatistic()
      if data_bf_inverted is not None:
        tf_data_inverted = ROOT.TFile(data_bf_inverted)
        tr = tf_data_inverted.Get("fitnt")
        tr.GetEntry(0)
        bf_x_data = tr.best_fit_params.BestFitValue(p1)
        bf_y_data = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_data = tr.best_fit_params.GetBestFitStatistic()
        
    
    h_delta_asim = GetDeltaHist(h_asim,bf_chisq_asim)
    h_delta_asim.SetMinimum(0)
    
    h_contour_asim = ROOT.TH2D(h_delta_asim)
    h_contour_asim.SetContour(2)
    h_contour_asim.SetContourLevel(1,conf90)
    h_contour_asim.SetContourLevel(0,conf68)
    h_contour_asim.SetLineWidth(1)
    h_contour_asim.SetLineColor(R.kRed)
    h_contour_asim.SetMinimum(0)
    
    h_delta_data = GetDeltaHist(h_data,bf_chisq_data)
    h_delta_data.SetMinimum(0)
    
    h_contour_data = ROOT.TH2D(h_delta_data)
    h_contour_data.SetContour(2)
    h_contour_data.SetContourLevel(1,conf90)
    h_contour_data.SetContourLevel(0,conf68)
    h_contour_data.SetLineWidth(1)
    h_contour_data.SetLineStyle(2)
    h_contour_data.SetLineColor(R.kBlack)
    h_contour_data.SetMinimum(0)
    
    contours_data = GetContours(h_contour_data)
    contours_asim = GetContours(h_contour_asim)
    
    marker_data = ROOT.TMarker(bf_x_data,bf_y_data,ROOT.kFullStar)
    marker_asim = ROOT.TMarker(bf_x_asim,bf_y_asim,ROOT.kFullStar)
    
    if (bf_chisq_data > bf_chisq_asim):
      bf_is_asim = True
      marker_asim.SetMarkerStyle(ROOT.kOpenStar)
    else:
      bf_is_asim = False
      marker_data.SetMarkerStyle(ROOT.kOpenStar)
      
    marker_data.SetMarkerColor(ROOT.kViolet);
    marker_asim.SetMarkerColor(ROOT.kBlue+1);
    marker_data.SetMarkerSize(1.5);
    marker_asim.SetMarkerSize(1.5);
    
    CopyMarkerToHist(marker_data,h_contour_data)
    CopyMarkerToHist(marker_asim,h_contour_asim)
    
    print "Bestfit point (asim): ",bf_x_asim,bf_y_asim,"\twith chisq: ",bf_chisq_asim
    print "Bestfit point (data): ",bf_x_data,bf_y_data,"\twith chisq: ",bf_chisq_data
      
    print "\n\n"
    
    # Start drawing
    for zoom in [False,True]:
    #~ for zoom in [True]:
      if zoom: h_contour_asim.GetYaxis().SetRangeUser(1E-3,4.5E-3)
      
      for g in contours_data[0]: g.SetLineColor(h_contour_data.GetLineColor())
      for g in contours_data[1]: g.SetLineColor(h_contour_data.GetLineColor())
      for g in contours_data[0]: g.SetLineStyle(line_style_68)
      for g in contours_data[1]: g.SetLineStyle(line_style_90)
      
      for g in contours_asim[0]: g.SetLineColor(h_contour_asim.GetLineColor())
      for g in contours_asim[1]: g.SetLineColor(h_contour_asim.GetLineColor())
      for g in contours_asim[0]: g.SetLineStyle(line_style_68)
      for g in contours_asim[1]: g.SetLineStyle(line_style_90)
    
      
      canv = ROOT.TCanvas("canv1","canv1",1280,720)
      canv.cd()
      SetBLMargins(canv)
      SetRMargins (canv,False)
      
      h_contour_asim.Draw("AXIS")
      for icont,cont_list in enumerate(contours_data + contours_asim):
        if icont%2==0 and not draw_68: continue
        for cont in cont_list:
          cont.Draw("same")
          
      # Draw markers
      marker_data.Draw()
      marker_asim.Draw()
      
      # Make legend
      
      legend = ROOT.TLegend(0.6,0.7,0.9,0.9)
      legend.AddEntry(h_contour_asim,"Sensitivity","lp")
      legend.AddEntry(h_contour_data,"Real data","lp")
      legend.SetFillStyle(0)
      legend.SetBorderSize(0)
      legend.Draw()
      
      if tag is None:
        plot_name = "profiled_asim_data"
      else:
        plot_name = tag
        
      if zoom: plot_name += "_zoomed"
      
      canv.Print(plot_name+"_"+hierarchy+"."+image_format)
      
      canv_nh = ROOT.TCanvas("canv_nh","canv_nh",1280,720)
      canv_nh.cd()
      SetBLMargins(canv_nh)
      SetRMargins (canv_nh,True)
      h_delta_asim.SetContour(100)
      h_delta_asim.Draw("colz")
      marker_asim.Draw()
      for g in contours_asim[1]:
        g.SetLineStyle(1)
        g.SetLineColor(R.kBlue)
        g.Draw("same")
      canv_nh.Print(plot_name+"_"+hierarchy+"_asimonly."+image_format)
      
      canv_ih = ROOT.TCanvas("canv_ih","canv_ih",1280,720)
      canv_ih.cd()
      SetBLMargins(canv_ih)
      SetRMargins (canv_ih,True)
      h_delta_data.SetContour(100)
      h_delta_data.Draw("colz")
      marker_data.Draw()
      for g in contours_data[1]:
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
  
def SetBLMargins(canv):
  canv.SetBottomMargin(0.16)
  canv.SetLeftMargin  (0.14)
  return
  
def SetRMargins(canv,needZ = False):
  if needZ: 
    canv.SetRightMargin(0.1)
  else: 
    canv.SetRightMargin(0.05)
  return
  
if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make contours")
  
  parser.add_option('--asim-bf-nh' , dest='asim_bf_normal' ,  default=None, help="File containing true best-fit (asim, NH)")
  parser.add_option('--asim-bf-ih' , dest='asim_bf_inverted', default=None, help="File containing true best-fit (asim, IH)")
  parser.add_option('--data-bf-nh' , dest='data_bf_normal' ,  default=None, help="File containing true best-fit (data, NH)")
  parser.add_option('--data-bf-ih' , dest='data_bf_inverted', default=None, help="File containing true best-fit (data, IH)")
  parser.add_option('-t','--tag' , dest='tag', default=None, help="Tag for output")
  parser.add_option('--no-68' , dest='no68', default=False,action="store_true", help="Should we remove the 68% confidence region")
  
  o,args = parser.parse_args()
  
  main(args[0],args[1],o.asim_bf_normal,o.data_bf_normal,o.asim_bf_inverted,o.data_bf_inverted,o.tag,not o.no68)

