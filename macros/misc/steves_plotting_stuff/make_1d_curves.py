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
conf90 = 2.71
conf68 = 1

to3dp = "{0:.3f}"
to3sf = "{0:.3g}"
to6dp = "{0:.6f}"

params = dict()
params["mixing23"] = conv.kPhypMixing23Bar
params["dmsq"    ] = conv.kPhypDm23sqBar
precision = dict()
precision["mixing23"] = to3dp
precision["dmsq"    ] = to3sf

def main(data_file_mixing, data_file_dmsq, data_bf_normal_mixing =None, data_bf_inverted_mixing =None, data_bf_normal_dmsq =None, data_bf_inverted_dmsq =None, tag=None):
  
  ROOT.gStyle.SetPalette(56);
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  
  # m = nh
  # p = ih
  
  bf_nh = dict([("mixing23",data_bf_normal_mixing),
                ("dmsq"    ,data_bf_normal_dmsq)])
  bf_ih = dict([("mixing23",data_bf_inverted_mixing),
                ("dmsq"    ,data_bf_inverted_dmsq) ])
            
  data_files = dict([("mixing23",data_file_mixing),
                     ("dmsq",data_file_dmsq) ])
            
  for param in ["mixing23","dmsq"]:
    
    m_f = ROOT.TFile(data_files[param])
    h_m = m_f.Get("h_chisq_"+param+"_Ave_nh")
    assert(h_m)
    h_m.Print()
    min_m = h_m.GetMinimum()
    
    m_p = ROOT.TFile(data_files[param])
      
    h_p = m_p.Get("h_chisq_"+param+"_Ave_ih")
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
    
    x_bin_p = ROOT.Long()
    x_bin_m = ROOT.Long()
    dummy = ROOT.Long()
    
    h_p.GetBinXYZ(bf_bin_p,x_bin_p,dummy,dummy)
    h_m.GetBinXYZ(bf_bin_m,x_bin_m,dummy,dummy)
    
    bf_x_p = h_p.GetXaxis().GetBinCenter(x_bin_p)
    bf_x_m = h_m.GetXaxis().GetBinCenter(x_bin_m)
    
    
    # If we have files, go and get the bestfits from those instead.
    if bf_nh[param] is not None:
      tf_data_normal = ROOT.TFile(bf_nh[param])
      tr = tf_data_normal.Get("fitnt")
      tr.GetEntry(0)
      bf_x_p = tr.best_fit_params.BestFitValue(params[param])
      
      bf_chisq_p = tr.best_fit_params.GetBestFitStatistic()
      
    if bf_ih[param] is not None:
      tf_data_inverted = ROOT.TFile(bf_ih[param])
      tr = tf_data_inverted.Get("fitnt")
      tr.GetEntry(0)
      bf_x_p = tr.best_fit_params.BestFitValue(params[param])
      bf_chisq_p = tr.best_fit_params.GetBestFitStatistic()
        
    
    h_delta_m = GetDeltaHist(h_m,bf_chisq_m)
    h_delta_p = GetDeltaHist(h_p,bf_chisq_p)
    
    h_delta_m.SetLineColor(R.kBlack)
    h_delta_p.SetLineColor(R.kBlack)
    h_delta_m.SetLineStyle(1)
    h_delta_p.SetLineStyle(2)
    
    if (bf_chisq_p > bf_chisq_m):
      bf_is_normal = True
    else:
      bf_is_normal = False
    
    print "Bestfit point (nh)  : ",bf_x_m,"\twith chisq: ",bf_chisq_m
    print "Bestfit point (ih)  : ",bf_x_p,"\twith chisq: ",bf_chisq_p
    
    neg_90_nh_x = None
    ibin_nh = bf_bin_m
    while ibin_nh > 0 and neg_90_nh_x is None:
      if h_delta_m.GetBinContent(ibin_nh) > conf90:
        neg_90_nh_x = h_delta_m.GetXaxis().GetBinUpEdge(ibin_nh)
      ibin_nh -= 1
      
    pos_90_nh_x = None
    ibin_nh = bf_bin_m
    while ibin_nh < h_delta_m.GetXaxis().GetNbins()+1 and pos_90_nh_x is None:
      if h_delta_m.GetBinContent(ibin_nh) > conf90:
        pos_90_nh_x = h_delta_m.GetXaxis().GetBinLowEdge(ibin_nh)
      ibin_nh += 1
      
    neg_90_ih_x = None
    ibin_ih = bf_bin_p
    while ibin_ih > 0 and neg_90_ih_x is None:
      if h_delta_p.GetBinContent(ibin_ih) > conf90:
        neg_90_ih_x = h_delta_p.GetXaxis().GetBinUpEdge(ibin_ih)
      ibin_ih -= 1
      
    pos_90_ih_x = None
    ibin_ih = bf_bin_p
    while ibin_ih < h_delta_p.GetXaxis().GetNbins()+1 and pos_90_ih_x is None:
      if h_delta_p.GetBinContent(ibin_ih) > conf90:
        pos_90_ih_x = h_delta_p.GetXaxis().GetBinLowEdge(ibin_ih)
      ibin_ih += 1
     
    print "\n\n"
  
    canv = ROOT.TCanvas("canv1","canv1",1280,720)
    canv.cd()
    
    h_delta_m.SetMinimum(0)
    h_delta_p.SetMinimum(0)
    
    param_title = util.ParamAsStringPaperRooTex(params[param])
    h_delta_m.GetXaxis().SetTitle(param_title)
    h_delta_p.GetXaxis().SetTitle(param_title)
    
    ymin = 0
    ymax = h_delta_m.GetMaximum()*1.1
    
    h_delta_m.Draw();
    h_delta_p.Draw("same")
    
    h_delta_m.SetMaximum(ymax)
    h_delta_p.SetMaximum(ymax)
    
    neg_line_nh = ROOT.TLine(neg_90_nh_x,ymin,neg_90_nh_x,ymax)
    pos_line_nh = ROOT.TLine(pos_90_nh_x,ymin,pos_90_nh_x,ymax)
    neg_line_ih = ROOT.TLine(neg_90_ih_x,ymin,neg_90_ih_x,ymax)
    pos_line_ih = ROOT.TLine(pos_90_ih_x,ymin,pos_90_ih_x,ymax)
    
    vertical_lines = [neg_line_nh,pos_line_nh,neg_line_ih,pos_line_ih]
    neg_line_ih.SetLineStyle(2)
    pos_line_ih.SetLineStyle(2)
    
    line_h90 = ROOT.TLine(h_delta_m.GetXaxis().GetXmin(),conf90,h_delta_m.GetXaxis().GetXmax(),conf90)
    line_h90.SetLineColor(R.kRed)
    line_h90.Draw("same")
    
    for l in vertical_lines:
      l.SetLineColor(R.kBlue)
      l.Draw("same")
    
    # Make legend
    
    legend = ROOT.TLegend(0.7,0.7,0.9,0.9)
    legend.AddEntry(h_delta_m,"NH","l")
    legend.AddEntry(h_delta_p,"IH","l")
    legend.Draw()
    
    prec = precision[param]
    
    pave = ROOT.TPaveText(0.4,0.6,0.6,0.89,"ndc")
    pave.SetShadowColor(0)
    pave.SetFillColor(0)
    pave.SetFillStyle(0)
    pave.SetLineColor(0)
    pave.AddText("NH: "+prec.format(neg_90_nh_x)+" < "+param_title+" "+prec.format(pos_90_nh_x)+" (at 90%)")
    pave.AddText("IH: "+prec.format(neg_90_ih_x)+" < "+param_title+" "+prec.format(pos_90_ih_x)+" (at 90%)")
    pave.Draw("same")
    
    if tag is None:
      plot_name = "confidence_1D_"+param
    else:
      plot_name = tag
      
    #~ if zoom: plot_name += "_zoomed"
    
    canv.Print(plot_name+"_"+param+"."+image_format)
    
    canv_nh = ROOT.TCanvas("canv_nh","canv_nh",1280,720)
    canv_nh.cd()
    h_delta_m.SetContour(100)
    h_delta_m.Draw("")
    canv_nh.Print(plot_name+"_"+param+"_nh."+image_format)
    
    canv_ih = ROOT.TCanvas("canv_ih","canv_ih",1280,720)
    canv_ih.cd()
    h_delta_p.SetContour(100)
    h_delta_p.Draw("")
    canv_ih.Print(plot_name+"_"+param+"_ih."+image_format)
  
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
      
  return h_delta
  
if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make contours")
  
  parser.add_option('--data-bf-nh' , dest='data_bf_normal' ,  default=None, help="File containing true best-fit (data, NH)")
  parser.add_option('--data-bf-ih' , dest='data_bf_inverted', default=None, help="File containing true best-fit (data, IH)")
  parser.add_option('-t','--tag' , dest='tag', default=None, help="Tag for output")
  
  o,args = parser.parse_args()
  
  main(args[0],args[1],o.data_bf_normal,o.data_bf_normal,o.data_bf_inverted,o.data_bf_inverted,o.tag)

