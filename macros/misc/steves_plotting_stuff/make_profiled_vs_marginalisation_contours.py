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
to4sf = "{0:.3e}"

p1 = conv.kPhypMixing23Bar
p2 = conv.kPhypDm23sqBar

# First file is profiled, second is marginalised
def main(prof_file, marg_file, prof_bf_normal = None, marg_bf_normal =None, prof_bf_inverted = None, marg_bf_inverted =None,
         tag=None,draw_68 = True):
  
  ROOT.gStyle.SetPalette(56);
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  
  title_size = 0.08
  label_size = 0.07
  
  line_style_68 = ROOT.kDashed
  line_style_90 = ROOT.kSolid
  
  assert(marg_file)
  
  print "Reading profiled from ",prof_file
  print "Reading marginalised from ",marg_file
  
  for hierarchy in ["nh"]:
    f_prof = ROOT.TFile(prof_file)
    h_prof = f_prof.Get("h_chisq_dmsq_mixing23_Ave_"+hierarchy)
    assert(h_prof)
    h_prof.Print()
    min_prof = h_prof.GetMinimum()
    
    f_marg = ROOT.TFile(marg_file)
      
    h_marg = f_marg.Get("h_chisq_marg")
    assert(h_marg)
    h_marg.Print()
    min_marg = h_marg.GetMinimum()
    
    print min_prof
    print min_marg
    
    # Set title offset
    h_prof.GetXaxis().SetTitleOffset(0.8)
    h_marg.GetXaxis().SetTitleOffset(0.8)
    h_prof.GetYaxis().SetTitleOffset(0.9)
    h_marg.GetYaxis().SetTitleOffset(0.9)
      
    h_prof.GetXaxis().SetTitleSize(title_size)
    h_marg.GetXaxis().SetTitleSize(title_size)
    h_prof.GetYaxis().SetTitleSize(title_size)
    h_marg.GetYaxis().SetTitleSize(title_size)
    h_prof.GetXaxis().SetLabelSize(label_size)
    h_marg.GetXaxis().SetLabelSize(label_size)
    #~ h_prof.GetYaxis().SetLabelSize(label_size)
    #~ h_marg.GetYaxis().SetLabelSize(label_size)
    
    # Precalculate the best-fit points from the histogram
    bf_bin_marg = h_marg.GetMinimumBin()
    bf_bin_prof = h_prof.GetMinimumBin()
    
    bf_chisq_marg = h_marg.GetBinContent(bf_bin_marg)
    bf_chisq_prof = h_prof.GetBinContent(bf_bin_prof)
    
    x_bin_marg = ROOT.Long() ; y_bin_marg = ROOT.Long()
    x_bin_prof = ROOT.Long() ; y_bin_prof = ROOT.Long()
    dummy = ROOT.Long()
    
    h_marg.GetBinXYZ(bf_bin_marg,x_bin_marg,y_bin_marg,dummy)
    h_prof.GetBinXYZ(bf_bin_prof,x_bin_prof,y_bin_prof,dummy)
    
    bf_x_marg = h_marg.GetXaxis().GetBinCenter(x_bin_marg)
    bf_x_prof = h_prof.GetXaxis().GetBinCenter(x_bin_prof)
    bf_y_marg = h_marg.GetYaxis().GetBinCenter(y_bin_marg)
    bf_y_prof = h_prof.GetYaxis().GetBinCenter(y_bin_prof)
    
    
    # If we have files, go and get the bestfits from those instead.
    if (hierarchy=="nh"):
      if prof_bf_normal is not None:
        tf_prof_normal = ROOT.TFile(prof_bf_normal)
        tr = tf_prof_normal.Get("fitnt")
        tr.GetEntry(0)
        bf_x_prof = tr.best_fit_params.BestFitValue(p1)
        bf_y_prof = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_prof = tr.best_fit_params.GetBestFitStatistic()
      if marg_bf_normal is not None:
        tf_marg_normal = ROOT.TFile(marg_bf_normal)
        tr = tf_marg_normal.Get("fitnt")
        tr.GetEntry(0)
        bf_x_marg = tr.best_fit_params.BestFitValue(p1)
        bf_y_marg = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_marg = tr.best_fit_params.GetBestFitStatistic()
    else:
      if prof_bf_inverted is not None:
        tf_prof_inverted = ROOT.TFile(prof_bf_inverted)
        tr = tf_prof_inverted.Get("fitnt")
        tr.GetEntry(0)
        bf_x_prof = tr.best_fit_params.BestFitValue(p1)
        bf_y_prof = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_prof = tr.best_fit_params.GetBestFitStatistic()
      if marg_bf_inverted is not None:
        tf_marg_inverted = ROOT.TFile(marg_bf_inverted)
        tr = tf_marg_inverted.Get("fitnt")
        tr.GetEntry(0)
        bf_x_marg = tr.best_fit_params.BestFitValue(p1)
        bf_y_marg = tr.best_fit_params.BestFitValue(p2)
        bf_chisq_marg = tr.best_fit_params.GetBestFitStatistic()
        
    
    h_delta_prof = GetDeltaHist(h_prof,bf_chisq_prof)
    
    h_contour_prof = ROOT.TH2D(h_delta_prof)
    h_contour_prof.SetContour(2)
    h_contour_prof.SetContourLevel(1,conf90)
    h_contour_prof.SetContourLevel(0,conf68)
    h_contour_prof.SetLineWidth(1)
    h_contour_prof.SetLineColor(R.kBlack)
    
    h_delta_marg = GetDeltaHist(h_marg,bf_chisq_marg)
    
    h_contour_marg = ROOT.TH2F(h_delta_marg)
    h_contour_marg.SetContour(2)
    h_contour_marg.SetContourLevel(1,conf90)
    h_contour_marg.SetContourLevel(0,conf68)
    h_contour_marg.SetLineWidth(1)
    h_contour_marg.SetLineStyle(2)
    h_contour_marg.SetLineColor(R.kBlack)
    
    contours_marg = GetContours(h_contour_marg)
    contours_prof = GetContours(h_contour_prof)
    
    marker_marg = ROOT.TMarker(bf_x_marg,bf_y_marg,ROOT.kFullStar)
    marker_prof = ROOT.TMarker(bf_x_prof,bf_y_prof,ROOT.kFullStar)
    
    if (bf_chisq_marg > bf_chisq_prof):
      bf_is_profiled = True
      marker_prof.SetMarkerStyle(ROOT.kOpenStar)
    else:
      bf_is_profiled = False
      marker_marg.SetMarkerStyle(ROOT.kOpenStar)
      
    marker_marg.SetMarkerColor(ROOT.kViolet);
    marker_prof.SetMarkerColor(ROOT.kBlue+1);
    marker_marg.SetMarkerSize(1.5);
    marker_prof.SetMarkerSize(1.5);
    
    CopyMarkerToHist(marker_marg,h_contour_marg)
    CopyMarkerToHist(marker_prof,h_contour_prof)
    
    print "Bestfit point (profiled)  :   ",util.ParamAsStringSimple(p1)+" = ",to6dp.format(bf_x_prof),", ",
    print util.ParamAsStringSimple(p2)+" = "+to6dp.format(bf_y_prof),"\twith chisq: ",to4sf.format(bf_chisq_prof)
    print "Bestfit point (marginalised): ",util.ParamAsStringSimple(p1)+" = ",to6dp.format(bf_x_marg),", ",
    print util.ParamAsStringSimple(p2)+" = "+to6dp.format(bf_y_marg),"\twith chisq: ",to4sf.format(bf_chisq_marg)
      
    print "\n\n"
    
    
    # Dump graphs for Kirsty
    outfile = ROOT.TFile(tag+"_valor_marg_contours.root","RECREATE")
    outfile.cd()
    
    h_delta_prof.Write("valor_prof_"+hierarchy+"_dchisq")
    h_delta_marg.Write("valor_marg_"+hierarchy+"_dchisq")
  
    contours_prof[0][0].Write("valor_prof_"+hierarchy+"_68_contour")
    contours_prof[1][0].Write("valor_prof_"+hierarchy+"_90_contour")
    contours_marg[0][0].Write("valor_marg_"+hierarchy+"_68_contour")
    contours_marg[1][0].Write("valor_marg_"+hierarchy+"_90_contour")
    
    bf_prof_graph = ROOT.TGraph(1)
    bf_marg_graph = ROOT.TGraph(1)
    bf_prof_graph.SetTitle("Bestfit: Normal")
    bf_marg_graph.SetTitle("Bestfit: Inverted")
    bf_prof_graph.SetPoint(0,bf_x_prof,bf_y_prof)
    bf_marg_graph.SetPoint(0,bf_x_marg,bf_y_marg)
    bf_prof_graph.Write("valor_prof_"+hierarchy+"_bestfit")
    bf_marg_graph.Write("valor_marg_"+hierarchy+"_bestfit")
    outfile.Close()
  
    
    # Start drawing
    
    
    #~ for zoom in [False,True]:
    for zoom in [True]:
      if zoom: h_contour_prof.GetYaxis().SetRangeUser(1E-3,4.5E-3)
      
      for g in contours_marg[0]: g.SetLineColor(R.kBlack)
      for g in contours_marg[1]: g.SetLineColor(R.kRed)
      for g in contours_marg[0]: g.SetLineStyle(h_contour_marg.GetLineStyle())
      for g in contours_marg[1]: g.SetLineStyle(h_contour_marg.GetLineStyle())
      
      for g in contours_prof[0]: g.SetLineColor(R.kBlack)
      for g in contours_prof[1]: g.SetLineColor(R.kRed)
      for g in contours_prof[0]: g.SetLineStyle(h_contour_prof.GetLineStyle())
      for g in contours_prof[1]: g.SetLineStyle(h_contour_prof.GetLineStyle())
    
      
      canv = ROOT.TCanvas("canv1","canv1",1280,720)
      canv.cd()
      SetBLMargins(canv)
      SetRMargins (canv,False)
      
      h_contour_prof.Draw("AXIS")
      for cont_list in contours_marg + contours_prof:
        for cont in cont_list:
          cont.Draw("same")
          
      # Draw markers
      marker_marg.Draw()
      marker_prof.Draw()
      
      # Make legend
      
      legend = ROOT.TLegend(0.7,0.7,0.93,0.93)
      legend.AddEntry(h_contour_prof,"Profiling","lp")
      legend.AddEntry(h_contour_marg,"Marginalisation","lp")
      legend.SetFillStyle(0)
      legend.SetBorderSize(0)
      legend.Draw()
      
      if tag is None:
        plot_name = "profiled_vs_marginalised"
      else:
        plot_name = tag
        
      if zoom: plot_name += "_zoomed"
      
      canv.Print(plot_name+"_"+hierarchy+"."+image_format)
      
      canv_nh = ROOT.TCanvas("canv_nh","canv_nh",1280,720)
      canv_nh.cd()
      SetBLMargins(canv_nh)
      SetRMargins (canv_nh,True)
      h_delta_prof.SetContour(100)
      h_delta_prof.Draw("colz")
      marker_prof.Draw()
      for g in contours_prof[1]:
        g.SetLineStyle(1)
        g.SetLineColor(R.kBlue)
        g.Draw("same")
      canv_nh.Print(plot_name+"_"+hierarchy+"_profiled."+image_format)
      
      canv_ih = ROOT.TCanvas("canv_ih","canv_ih",1280,720)
      canv_ih.cd()
      SetBLMargins(canv_ih)
      SetRMargins (canv_ih,True)
      h_delta_marg.SetContour(100)
      h_delta_marg.Draw("colz")
      marker_marg.Draw()
      for g in contours_marg[1]:
        g.SetLineStyle(1)
        g.SetLineColor(R.kBlue)
        g.Draw("same")
      canv_ih.Print(plot_name+"_"+hierarchy+"_marginalised."+image_format)

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
  h_delta = h.Clone()
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
  
  parser.add_option('--prof-bf-nh' , dest='prof_bf_normal' ,  default=None, help="File containing true best-fit (prof, NH)")
  parser.add_option('--prof-bf-ih' , dest='prof_bf_inverted', default=None, help="File containing true best-fit (prof, IH)")
  parser.add_option('--marg-bf-nh' , dest='marg_bf_normal' ,  default=None, help="File containing true best-fit (marg, NH)")
  parser.add_option('--marg-bf-ih' , dest='marg_bf_inverted', default=None, help="File containing true best-fit (marg, IH)")
  parser.add_option('-t','--tag' , dest='tag', default=None, help="Tag for output")
  parser.add_option('--no-68' , dest='no68', default=False,action="store_true", help="Should we remove the 68% confidence region")

  
  o,args = parser.parse_args()
  
  main(args[0],args[1],o.prof_bf_normal,o.marg_bf_normal,o.prof_bf_inverted,o.marg_bf_inverted,o.tag, not o.no68)

