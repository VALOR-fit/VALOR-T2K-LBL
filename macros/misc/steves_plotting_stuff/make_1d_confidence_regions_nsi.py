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
conf90 = 2.71
conf68 = 1

to3dp = "{0:.3f}"
to3sf = "{0:.3g}"
to6dp = "{0:.6f}"

params = dict()
params["epsmutau"] = conv.kPhypNSIepsMuTau
precision = dict()
precision["epsmutau"] = to3dp

flat_style  = R.kDotted
minos_style = R.kSolid

def main(data_file_epsmutau_flat, data_file_epsmutau_minos, 
         data_bf_normal_epsmutau_flat  = None, data_bf_inverted_epsmutau_flat  = None,
         data_bf_normal_epsmutau_minos = None, data_bf_inverted_epsmutau_minos = None,
         tag=None,draw_hierarchies_individually = False):
  
  ROOT.gStyle.SetPalette(56);
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  
  # m = nh
  # p = ih
  
  param = "epsmutau"
  
  bf_nh_flat  = dict([("epsmutau",data_bf_normal_epsmutau_flat  )])
  bf_ih_flat  = dict([("epsmutau",data_bf_inverted_epsmutau_flat)])
  bf_nh_minos = dict([("epsmutau",data_bf_normal_epsmutau_minos  )])
  bf_ih_minos = dict([("epsmutau",data_bf_inverted_epsmutau_minos)])
            
  f_nh_flat = ROOT.TFile(data_file_epsmutau_flat)
  h_nh_flat = f_nh_flat.Get("h_chisq_epsmutau_Ave_nh")
  assert(h_nh_flat)
  h_nh_flat.SetName(h_nh_flat.GetName()+"_flat")
  h_nh_flat.Print()
  min_nh_flat = h_nh_flat.GetMinimum()
  
  f_ih_flat = ROOT.TFile(data_file_epsmutau_flat)
  h_ih_flat = f_ih_flat.Get("h_chisq_epsmutau_Ave_ih")
  assert(h_ih_flat)
  h_ih_flat.SetName(h_ih_flat.GetName()+"_flat")
  h_ih_flat.Print()
  min_ih_flat = h_ih_flat.GetMinimum()
  
  f_nh_minos = ROOT.TFile(data_file_epsmutau_minos)
  h_nh_minos = f_nh_minos.Get("h_chisq_epsmutau_Ave_nh")
  assert(h_nh_minos)
  h_nh_minos.SetName(h_nh_minos.GetName()+"_minos")
  h_nh_minos.Print()
  min_nh_minos = h_nh_minos.GetMinimum()
  
  f_ih_minos = ROOT.TFile(data_file_epsmutau_minos)
  h_ih_minos = f_ih_minos.Get("h_chisq_epsmutau_Ave_ih")
  assert(h_ih_minos)
  h_ih_minos.SetName(h_ih_minos.GetName()+"_minos")
  h_ih_minos.Print()
  min_ih_minos = h_ih_minos.GetMinimum()
  
  # Precalculate the best-fit points from the histogram
  bf_bin_nh_flat  = h_nh_flat .GetMinimumBin()
  bf_bin_ih_flat  = h_ih_flat .GetMinimumBin()
  bf_bin_nh_minos = h_nh_minos.GetMinimumBin()
  bf_bin_ih_minos = h_ih_minos.GetMinimumBin()
  
  bf_chisq_nh_flat  = h_nh_flat .GetBinContent(bf_bin_nh_flat )
  bf_chisq_ih_flat  = h_ih_flat .GetBinContent(bf_bin_ih_flat )
  bf_chisq_nh_minos = h_nh_minos.GetBinContent(bf_bin_nh_minos)
  bf_chisq_ih_minos = h_ih_minos.GetBinContent(bf_bin_ih_minos)
  
  dummy = ROOT.Long()
  
  x_bin_ih_flat  = ROOT.Long()
  x_bin_nh_flat  = ROOT.Long()
  
  h_ih_flat .GetBinXYZ(bf_bin_ih_flat ,x_bin_ih_flat ,dummy,dummy)
  h_nh_flat .GetBinXYZ(bf_bin_nh_flat ,x_bin_nh_flat ,dummy,dummy)
  
  bf_x_ih_flat  = h_ih_flat .GetXaxis().GetBinCenter(x_bin_ih_flat )
  bf_x_nh_flat  = h_nh_flat .GetXaxis().GetBinCenter(x_bin_nh_flat )
  
  x_bin_ih_minos = ROOT.Long()
  x_bin_nh_minos = ROOT.Long()
  
  h_ih_minos.GetBinXYZ(bf_bin_ih_minos,x_bin_ih_minos,dummy,dummy)
  h_nh_minos.GetBinXYZ(bf_bin_nh_minos,x_bin_nh_minos,dummy,dummy)
  
  bf_x_ih_minos = h_ih_minos.GetXaxis().GetBinCenter(x_bin_ih_minos)
  bf_x_nh_minos = h_nh_minos.GetXaxis().GetBinCenter(x_bin_nh_minos)
  
  
  # Big block of flat prior stuff
  
  # If we have files, go and get the bestfits from those instead.
  if bf_nh_flat [param] is not None:
    print "Flat NH BF: ",bf_nh_flat [param]
    tf_data_normal = ROOT.TFile(bf_nh_flat [param])
    tr = tf_data_normal.Get("fitnt")
    tr.GetEntry(0)
    bf_x_nh_flat  = tr.best_fit_params.BestFitValue(params[param])
    
    bf_chisq_nh_flat  = tr.best_fit_params.GetBestFitStatistic()
    
  if bf_ih_flat[param] is not None:
    tf_data_inverted = ROOT.TFile(bf_ih_flat[param])
    print "Flat IH BF: ",bf_ih_flat[param]
    tr = tf_data_inverted.Get("fitnt")
    tr.GetEntry(0)
    bf_x_ih_flat = tr.best_fit_params.BestFitValue(params[param])
    bf_chisq_ih_flat = tr.best_fit_params.GetBestFitStatistic()
      
  h_delta_nh_flat = GetDeltaHist(h_nh_flat,bf_chisq_nh_flat)
  h_delta_ih_flat = GetDeltaHist(h_ih_flat,bf_chisq_ih_flat)
  
  h_delta_nh_flat.SetLineColor(R.kBlack)
  h_delta_ih_flat.SetLineColor(R.kBlack)
  h_delta_nh_flat.SetLineStyle(flat_style)
  h_delta_ih_flat.SetLineStyle(flat_style)
  
  print "Bestfit point (nh_flat)  : ",bf_x_nh_flat,"\twith chisq: ",bf_chisq_nh_flat 
  print "Bestfit point (ih_flat)  : ",bf_x_ih_flat,"\twith chisq: ",bf_chisq_ih_flat
  
  neg_90_nh_flat_x  = None
  ibin_nh_flat  = bf_bin_nh_flat 
  while ibin_nh_flat  > 0 and neg_90_nh_flat_x  is None:
    if h_delta_nh_flat.GetBinContent(ibin_nh_flat) > conf90:
      neg_90_nh_flat_x  = h_delta_nh_flat.GetXaxis().GetBinUpEdge(ibin_nh_flat)
    ibin_nh_flat  -= 1
    
  pos_90_nh_flat_x  = None
  ibin_nh_flat  = bf_bin_nh_flat 
  while ibin_nh_flat  < h_delta_nh_flat.GetXaxis().GetNbins()+1 and pos_90_nh_flat_x  is None:
    if h_delta_nh_flat.GetBinContent(ibin_nh_flat) > conf90:
      pos_90_nh_flat_x  = h_delta_nh_flat.GetXaxis().GetBinLowEdge(ibin_nh_flat)
    ibin_nh_flat  += 1
    
  neg_90_ih_flat_x = None
  ibin_ih_flat = bf_bin_ih_flat
  while ibin_ih_flat > 0 and neg_90_ih_flat_x is None:
    if h_delta_ih_flat.GetBinContent(ibin_ih_flat) > conf90:
      neg_90_ih_flat_x = h_delta_ih_flat.GetXaxis().GetBinUpEdge(ibin_ih_flat)
    ibin_ih_flat -= 1
    
  pos_90_ih_flat_x = None
  ibin_ih_flat = bf_bin_ih_flat
  while ibin_ih_flat < h_delta_ih_flat.GetXaxis().GetNbins()+1 and pos_90_ih_flat_x is None:
    if h_delta_ih_flat.GetBinContent(ibin_ih_flat) > conf90:
      pos_90_ih_flat_x = h_delta_ih_flat.GetXaxis().GetBinLowEdge(ibin_ih_flat)
    ibin_ih_flat += 1
   
  print "\n\n"
  
  h_delta_nh_flat.GetYaxis().SetTitleOffset(0.6)
  # End big block of flat prior stuff
  
  
  # Big block of minos prior stuff
  
  # If we have files, go and get the bestfits from those instead.
  if bf_nh_minos [param] is not None:
    tf_data_normal = ROOT.TFile(bf_nh_minos [param])
    tr = tf_data_normal.Get("fitnt")
    tr.GetEntry(0)
    bf_x_nh_minos  = tr.best_fit_params.BestFitValue(params[param])
    bf_chisq_nh_minos  = tr.best_fit_params.GetBestFitStatistic()
  
  
  if bf_ih_minos[param] is not None:
    tf_data_inverted = ROOT.TFile(bf_ih_minos[param])
    tr = tf_data_inverted.Get("fitnt")
    tr.GetEntry(0)
    bf_x_ih_minos = tr.best_fit_params.BestFitValue(params[param])
    bf_chisq_ih_minos = tr.best_fit_params.GetBestFitStatistic()
    
  h_delta_nh_minos = GetDeltaHist(h_nh_minos,bf_chisq_nh_minos)
  h_delta_ih_minos = GetDeltaHist(h_ih_minos,bf_chisq_ih_minos)

  
  h_delta_nh_minos.SetLineColor(R.kBlack)
  h_delta_ih_minos.SetLineColor(R.kBlack)
  h_delta_nh_minos.SetLineStyle(minos_style)
  h_delta_ih_minos.SetLineStyle(minos_style)
  
  print "Bestfit point (nh_minos)  : ",bf_x_nh_minos,"\twith chisq: ",bf_chisq_nh_minos 
  print "Bestfit point (ih_minos)  : ",bf_x_ih_minos,"\twith chisq: ",bf_chisq_ih_minos
  
  neg_90_nh_minos_x = None
  ibin_nh_minos  = bf_bin_nh_minos 
  while ibin_nh_minos  > 0 and neg_90_nh_minos_x is None:
    if h_delta_nh_minos.GetBinContent(ibin_nh_minos) > conf90:
      neg_90_nh_minos_x = h_delta_nh_minos.GetXaxis().GetBinUpEdge(ibin_nh_minos)
    ibin_nh_minos  -= 1
    
  pos_90_nh_minos_x = None
  ibin_nh_minos  = bf_bin_nh_minos 
  while ibin_nh_minos  < h_delta_nh_minos.GetXaxis().GetNbins()+1 and pos_90_nh_minos_x is None:
    if h_delta_nh_minos.GetBinContent(ibin_nh_minos) > conf90:
      pos_90_nh_minos_x = h_delta_nh_minos.GetXaxis().GetBinLowEdge(ibin_nh_minos)
    ibin_nh_minos  += 1
    
  neg_90_ih_minos_x = None
  ibin_ih_minos = bf_bin_ih_minos
  while ibin_ih_minos > 0 and neg_90_ih_minos_x is None:
    if h_delta_ih_minos.GetBinContent(ibin_ih_minos) > conf90:
      neg_90_ih_minos_x = h_delta_ih_minos.GetXaxis().GetBinUpEdge(ibin_ih_minos)
    ibin_ih_minos -= 1
    
  pos_90_ih_minos_x = None
  ibin_ih_minos = bf_bin_ih_minos
  while ibin_ih_minos < h_delta_ih_minos.GetXaxis().GetNbins()+1 and pos_90_ih_minos_x is None:
    if h_delta_ih_minos.GetBinContent(ibin_ih_minos) > conf90:
      pos_90_ih_minos_x = h_delta_ih_minos.GetXaxis().GetBinLowEdge(ibin_ih_minos)
    ibin_ih_minos += 1
   
  print "\n\n"
  
  h_delta_nh_minos.GetYaxis().SetTitleOffset(0.6)
  
  # End big block of minos prior stuff
  
  h_delta_nh_flat .SetMinimum(0)
  h_delta_ih_flat .SetMinimum(0)
  h_delta_nh_minos.SetMinimum(0)
  h_delta_ih_minos.SetMinimum(0)
  
  param_title_nh_flat  = h_delta_nh_flat .GetXaxis().GetTitle()
  param_title_ih_flat  = h_delta_ih_flat .GetXaxis().GetTitle()
  if (param_title_nh_flat  != param_title_ih_flat ):
    param_title = param_title_nh_flat +" or "+param_title_ih_flat 
    h_delta_nh_flat .GetXaxis().SetTitle(param_title)
    h_delta_ih_flat .GetXaxis().SetTitle(param_title)
  else:
    param_title = param_title_nh_flat 
    
  param_title_nh_minos = h_delta_nh_minos.GetXaxis().GetTitle()
  param_title_ih_minos = h_delta_ih_minos.GetXaxis().GetTitle()
  if (param_title_nh_minos != param_title_ih_minos):
    param_title = param_title_nh_minos+" or "+param_title_ih_minos
    h_delta_nh_minos.GetXaxis().SetTitle(param_title)
    h_delta_ih_minos.GetXaxis().SetTitle(param_title)
  else:
    param_title = param_title_nh_minos
  
  ymin = 0.
  ymax_nh = 1.1*max([h_delta_nh_flat .GetMaximum(),h_delta_nh_minos.GetMaximum()])
  ymax_ih = 1.1*max([h_delta_ih_flat .GetMaximum(),h_delta_ih_minos.GetMaximum()])
  
  h_delta_nh_flat .SetMaximum(ymax_nh)
  h_delta_ih_flat .SetMaximum(ymax_ih)
  
  neg_line_nh_flat  = ROOT.TLine(neg_90_nh_flat_x ,ymin,neg_90_nh_flat_x ,ymax_nh)
  pos_line_nh_flat  = ROOT.TLine(pos_90_nh_flat_x ,ymin,pos_90_nh_flat_x ,ymax_nh)
  neg_line_ih_flat  = ROOT.TLine(neg_90_ih_flat_x ,ymin,neg_90_ih_flat_x ,ymax_ih)
  pos_line_ih_flat  = ROOT.TLine(pos_90_ih_flat_x ,ymin,pos_90_ih_flat_x ,ymax_ih)
  neg_line_nh_minos = ROOT.TLine(neg_90_nh_minos_x,ymin,neg_90_nh_minos_x,ymax_nh)
  pos_line_nh_minos = ROOT.TLine(pos_90_nh_minos_x,ymin,pos_90_nh_minos_x,ymax_nh)
  neg_line_ih_minos = ROOT.TLine(neg_90_ih_minos_x,ymin,neg_90_ih_minos_x,ymax_ih)
  pos_line_ih_minos = ROOT.TLine(pos_90_ih_minos_x,ymin,pos_90_ih_minos_x,ymax_ih)
  
  vertical_lines_nh = [neg_line_nh_flat ,pos_line_nh_flat ,neg_line_nh_minos,pos_line_nh_minos]
  vertical_lines_ih = [neg_line_ih_flat ,pos_line_ih_flat ,neg_line_ih_minos,pos_line_ih_minos]
  
  vertical_lines = vertical_lines_nh + vertical_lines_ih
  
  neg_line_nh_flat .SetLineStyle(flat_style)
  pos_line_nh_flat .SetLineStyle(flat_style)
  neg_line_nh_minos.SetLineStyle(minos_style)
  pos_line_nh_minos.SetLineStyle(minos_style)
  
  neg_line_ih_flat .SetLineStyle(flat_style)
  pos_line_ih_flat .SetLineStyle(flat_style)
  neg_line_ih_minos.SetLineStyle(minos_style)
  pos_line_ih_minos.SetLineStyle(minos_style)
  
  line_h90 = ROOT.TLine(h_delta_nh_flat .GetXaxis().GetXmin(),conf90,h_delta_nh_flat .GetXaxis().GetXmax(),conf90)
  line_h90.SetLineColor(R.kRed)
  
  
  # Make NH plot
  canv_nh = ROOT.TCanvas("canv_nh ","canv_nh ",1280,720)
  canv_nh.cd()
  SetBLMargins(canv_nh)
  SetRMargins (canv_nh ,False)
  
  h_delta_nh_flat .Draw()
  h_delta_nh_minos.Draw("same")
  line_h90.Draw("same")
  
  for l in vertical_lines_nh:
    l.SetLineColor(R.kBlue)
    l.Draw("same")
  
  # Make legend
  
  legend = ROOT.TLegend(0.45,0.25,0.6,0.4)
  legend.AddEntry(h_delta_nh_flat ,"Flat prior","l")
  legend.AddEntry(h_delta_nh_minos,"MINOS prior","l")
  legend.Draw()
  
  prec = precision[param]
  
  pave_nh = ROOT.TPaveText(0.365,0.7,0.68,0.89,"ndc")
  pave_nh.SetTextFont(22)
  pave_nh.SetShadowColor(0)
  pave_nh.SetFillColor(0)
  pave_nh.SetFillStyle(0)
  pave_nh.SetLineColor(0)
  pave_nh.AddText("Flat prior:  "+prec.format(neg_90_nh_flat_x )+" < "+param_title_nh_flat +" < "+prec.format(pos_90_nh_flat_x )+" (at 90%)")
  pave_nh.AddText("MINOS prior: "+prec.format(neg_90_nh_minos_x)+" < "+param_title_nh_minos+" < "+prec.format(pos_90_nh_minos_x)+" (at 90%)")
  pave_nh.Draw("same")
  
  if tag is None:
    plot_name = "confidence_1D_nh"
  else:
    plot_name = tag+"_nh"
    
  canv_nh.Print(plot_name+"_epsmutau."+image_format)
  
  
  # Make IH plot
  canv_ih = ROOT.TCanvas("canv_ih ","canv_ih ",1280,720)
  canv_ih.cd()
  SetBLMargins(canv_ih)
  SetRMargins (canv_ih ,False)
  
  h_delta_ih_flat .Draw()
  h_delta_ih_minos.Draw("same")
  line_h90.Draw("same")
  
  for l in vertical_lines_ih:
    l.SetLineColor(R.kBlue)
    l.Draw("same")
  
  # Make legend
  
  legend = ROOT.TLegend(0.45,0.25,0.6,0.4)
  legend.AddEntry(h_delta_ih_flat ,"Flat prior","l")
  legend.AddEntry(h_delta_ih_minos,"MINOS prior","l")
  legend.Draw()
  
  prec = precision[param]
  
  pave_ih = ROOT.TPaveText(0.365,0.7,0.68,0.89,"ndc")
  pave_ih.SetTextFont(22)
  pave_ih.SetShadowColor(0)
  pave_ih.SetFillColor(0)
  pave_ih.SetFillStyle(0)
  pave_ih.SetLineColor(0)
  pave_ih.AddText("Flat prior:  "+prec.format(neg_90_ih_flat_x )+" < "+param_title_ih_flat +" < "+prec.format(pos_90_ih_flat_x )+" (at 90%)")
  pave_ih.AddText("MINOS prior: "+prec.format(neg_90_ih_minos_x)+" < "+param_title_ih_minos+" < "+prec.format(pos_90_ih_minos_x)+" (at 90%)")
  pave_ih.Draw("same")
  
  if tag is None:
    plot_name = "confidence_1D_ih"
  else:
    plot_name = tag+"_ih"
    
  canv_ih.Print(plot_name+"_epsmutau."+image_format)
  
def CopyMarkerToHist(marker,histogram):
  histogram.SetMarkerSize (marker.GetMarkerSize ())
  histogram.SetMarkerStyle(marker.GetMarkerStyle())
  histogram.SetMarkerColor(marker.GetMarkerColor())
  
def GetDeltaHist(h,minval=None):
  h_delta = ROOT.TH1D(h)
  h_delta.Reset()
  if minval is None: minval = h.GetMinimum()
  for ix in xrange(h_delta.GetXaxis().GetNbins()):
    h_delta.SetBinContent(ix+1,h.GetBinContent(ix+1)-minval)
  h_delta.GetYaxis().SetTitle("\Delta \chi^{2}");
  h_delta.SetDirectory(0)
      
  return h_delta
  
def SetBLMargins(canv):
  canv.SetBottomMargin(0.15)
  canv.SetLeftMargin  (0.08)
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
  
  parser.add_option('--data-bf-nh-flat'  , dest='data_bf_normal_flat'   , default=None, help="File containing flat prior best-fit (data, NH)")
  parser.add_option('--data-bf-ih-flat'  , dest='data_bf_inverted_flat' , default=None, help="File containing flat prior best-fit (data, IH)")
  parser.add_option('--data-bf-nh-minos' , dest='data_bf_normal_minos'  , default=None, help="File containing MINOS prior best-fit (data, NH)")
  parser.add_option('--data-bf-ih-minos' , dest='data_bf_inverted_minos', default=None, help="File containing MINOS prior best-fit (data, IH)")
  parser.add_option('-t','--tag' , dest='tag', default=None, help="Tag for output")
  parser.add_option('-i','--individual-priors' , dest='individual_priors', default=False,action="store_true", help="Should we draw the independent prior plots")
  
  o,args = parser.parse_args()
  
  main(args[0],args[1],
    o.data_bf_normal_flat , o.data_bf_inverted_flat ,
    o.data_bf_normal_minos, o.data_bf_inverted_minos,
    o.tag,o.individual_priors)

