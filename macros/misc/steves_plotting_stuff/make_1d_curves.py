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
params["mixing23"] = conv.kPhypMixing23Bar
params["dmsq"    ] = conv.kPhypDm23sqBar
precision = dict()
precision["mixing23"] = to3dp
precision["dmsq"    ] = to3sf

def main(data_file_mixing, data_file_dmsq,
         data_bf_normal_mixing = None, data_bf_inverted_mixing = None,
         data_bf_normal_dmsq   = None, data_bf_inverted_dmsq   = None,
         tag=None,draw_hierarchies_individually = False):
  
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
    
    f_nh = ROOT.TFile(data_files[param])
    h_nh = f_nh.Get("h_chisq_"+param+"_Ave_nh")
    assert(h_nh)
    h_nh.Print()
    min_nh = h_nh.GetMinimum()
    
    f_ih = ROOT.TFile(data_files[param])
      
    h_ih = f_ih.Get("h_chisq_"+param+"_Ave_ih")
    assert(h_ih)
    h_ih.Print()
    min_ih = h_ih.GetMinimum()
    
    print min_nh
    print min_ih
    
    # Precalculate the best-fit points from the histogram
    bf_bin_ih = h_ih.GetMinimumBin()
    bf_bin_nh = h_nh.GetMinimumBin()
    
    bf_chisq_ih = h_ih.GetBinContent(bf_bin_ih)
    bf_chisq_nh = h_nh.GetBinContent(bf_bin_nh)
    
    x_bin_ih = ROOT.Long()
    x_bin_nh = ROOT.Long()
    dummy = ROOT.Long()
    
    h_ih.GetBinXYZ(bf_bin_ih,x_bin_ih,dummy,dummy)
    h_nh.GetBinXYZ(bf_bin_nh,x_bin_nh,dummy,dummy)
    
    bf_x_ih = h_ih.GetXaxis().GetBinCenter(x_bin_ih)
    bf_x_nh = h_nh.GetXaxis().GetBinCenter(x_bin_nh)
    
    
    # If we have files, go and get the bestfits from those instead.
    if bf_nh[param] is not None:
      print "NH BF: ",bf_nh[param]
      tf_data_normal = ROOT.TFile(bf_nh[param])
      tr = tf_data_normal.Get("fitnt")
      tr.GetEntry(0)
      bf_x_nh = tr.best_fit_params.BestFitValue(params[param])
      
      bf_chisq_nh = tr.best_fit_params.GetBestFitStatistic()
      
    if bf_ih[param] is not None:
      tf_data_inverted = ROOT.TFile(bf_ih[param])
      print "IH BF: ",bf_ih[param]
      tr = tf_data_inverted.Get("fitnt")
      tr.GetEntry(0)
      bf_x_ih = tr.best_fit_params.BestFitValue(params[param])
      bf_chisq_ih = tr.best_fit_params.GetBestFitStatistic()
        
    
    h_delta_nh = GetDeltaHist(h_nh,bf_chisq_nh)
    h_delta_ih = GetDeltaHist(h_ih,bf_chisq_ih)
    
    h_delta_nh.SetLineColor(R.kBlack)
    h_delta_ih.SetLineColor(R.kBlack)
    h_delta_nh.SetLineStyle(1)
    h_delta_ih.SetLineStyle(2)
    
    if (bf_chisq_ih > bf_chisq_nh):
      bf_is_normal = True
    else:
      bf_is_normal = False
    
    print "Bestfit point (nh)  : ",bf_x_nh,"\twith chisq: ",bf_chisq_nh
    print "Bestfit point (ih)  : ",bf_x_ih,"\twith chisq: ",bf_chisq_ih
    
    neg_90_nh_x = None
    ibin_nh = bf_bin_nh
    while ibin_nh > 0 and neg_90_nh_x is None:
      if h_delta_nh.GetBinContent(ibin_nh) > conf90:
        neg_90_nh_x = h_delta_nh.GetXaxis().GetBinUpEdge(ibin_nh)
      ibin_nh -= 1
      
    pos_90_nh_x = None
    ibin_nh = bf_bin_nh
    while ibin_nh < h_delta_nh.GetXaxis().GetNbins()+1 and pos_90_nh_x is None:
      if h_delta_nh.GetBinContent(ibin_nh) > conf90:
        pos_90_nh_x = h_delta_nh.GetXaxis().GetBinLowEdge(ibin_nh)
      ibin_nh += 1
      
    neg_90_ih_x = None
    ibin_ih = bf_bin_ih
    while ibin_ih > 0 and neg_90_ih_x is None:
      if h_delta_ih.GetBinContent(ibin_ih) > conf90:
        neg_90_ih_x = h_delta_ih.GetXaxis().GetBinUpEdge(ibin_ih)
      ibin_ih -= 1
      
    pos_90_ih_x = None
    ibin_ih = bf_bin_ih
    while ibin_ih < h_delta_ih.GetXaxis().GetNbins()+1 and pos_90_ih_x is None:
      if h_delta_ih.GetBinContent(ibin_ih) > conf90:
        pos_90_ih_x = h_delta_ih.GetXaxis().GetBinLowEdge(ibin_ih)
      ibin_ih += 1
     
    print "\n\n"
    
    h_delta_nh.GetYaxis().SetTitleOffset(0.6)
  
    canv = ROOT.TCanvas("canv1","canv1",1280,720)
    canv.cd()
    SetBLMargins(canv)
    SetRMargins (canv,False)
    
    h_delta_nh.SetMinimum(0)
    h_delta_ih.SetMinimum(0)
    
    param_title_nh = h_delta_nh.GetXaxis().GetTitle()
    param_title_ih = h_delta_ih.GetXaxis().GetTitle()
    if (param_title_nh != param_title_ih):
      param_title = param_title_nh+" or "+param_title_ih
      h_delta_nh.GetXaxis().SetTitle(param_title)
      h_delta_ih.GetXaxis().SetTitle(param_title)
    else:
      param_title = param_title_nh
    
    ymin = 0
    ymax = h_delta_nh.GetMaximum()*1.1
    
    h_delta_nh.Draw();
    h_delta_ih.Draw("same")
    
    h_delta_nh.SetMaximum(ymax)
    h_delta_ih.SetMaximum(ymax)
    
    neg_line_nh = ROOT.TLine(neg_90_nh_x,ymin,neg_90_nh_x,ymax)
    pos_line_nh = ROOT.TLine(pos_90_nh_x,ymin,pos_90_nh_x,ymax)
    neg_line_ih = ROOT.TLine(neg_90_ih_x,ymin,neg_90_ih_x,ymax)
    pos_line_ih = ROOT.TLine(pos_90_ih_x,ymin,pos_90_ih_x,ymax)
    
    vertical_lines = [neg_line_nh,pos_line_nh,neg_line_ih,pos_line_ih]
    neg_line_ih.SetLineStyle(2)
    pos_line_ih.SetLineStyle(2)
    
    line_h90 = ROOT.TLine(h_delta_nh.GetXaxis().GetXmin(),conf90,h_delta_nh.GetXaxis().GetXmax(),conf90)
    line_h90.SetLineColor(R.kRed)
    line_h90.Draw("same")
    
    for l in vertical_lines:
      l.SetLineColor(R.kBlue)
      l.Draw("same")
    
    # Make legend
    
    legend = ROOT.TLegend(0.8,0.7,0.9,0.9)
    legend.AddEntry(h_delta_nh,"NH","l")
    legend.AddEntry(h_delta_ih,"IH","l")
    legend.Draw()
    
    prec = precision[param]
    
    pave = ROOT.TPaveText(0.365,0.7,0.68,0.89,"ndc")
    pave.SetTextFont(22)
    pave.SetShadowColor(0)
    pave.SetFillColor(0)
    pave.SetFillStyle(0)
    pave.SetLineColor(0)
    pave.AddText("NH: "+prec.format(neg_90_nh_x)+" < "+param_title_nh+" < "+prec.format(pos_90_nh_x)+" (at 90%)")
    pave.AddText("IH: "+prec.format(neg_90_ih_x)+" < "+param_title_ih+" < "+prec.format(pos_90_ih_x)+" (at 90%)")
    pave.Draw("same")
    
    if tag is None:
      plot_name = "confidence_1D_"+param
    else:
      plot_name = tag
      
    #~ if zoom: plot_name += "_zoomed"
    
    canv.Print(plot_name+"_"+param+"."+image_format)
    
    if not draw_hierarchies_individually: continue
    
    canv_nh = ROOT.TCanvas("canv_nh","canv_nh",1280,720)
    canv_nh.cd()
    SetBLMargins(canv_nh)
    SetRMargins (canv_nh,False)
    h_delta_nh.SetContour(100)
    h_delta_nh.Draw("")
    canv_nh.Print(plot_name+"_"+param+"_nh."+image_format)
    
    canv_ih = ROOT.TCanvas("canv_ih","canv_ih",1280,720)
    canv_ih.cd()
    SetBLMargins(canv_ih)
    SetRMargins (canv_ih,False)
    h_delta_ih.SetContour(100)
    h_delta_ih.Draw("")
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
  
  parser.add_option('--data-bf-nh' , dest='data_bf_normal' ,  default=None, help="File containing true best-fit (data, NH)")
  parser.add_option('--data-bf-ih' , dest='data_bf_inverted', default=None, help="File containing true best-fit (data, IH)")
  parser.add_option('-t','--tag' , dest='tag', default=None, help="Tag for output")
  parser.add_option('-i','--individual-hierarchies' , dest='individual_hierarchy', default=False,action="store_true", help="Should we draw the independent hierarchy plots")
  
  o,args = parser.parse_args()
  
  main(args[0],args[1],o.data_bf_normal,o.data_bf_inverted,o.data_bf_normal,o.data_bf_inverted,o.tag,o.individual_hierarchy)

