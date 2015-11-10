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

to2dp = "{0:.2f}"
to3dp = "{0:.3f}"
to3sf = "{0:.3g}"

chisq_max = 30
nbins = 30

linewidth = 2

chisq_pdf = R.TF1("chisq","1/(TMath::Power(2,[0]/2) * TMath::Gamma([0]/2)) * TMath::Power(x,[0]/2 - 1) * TMath::Exp(-x/2)")
chisq_pdf.SetParName(0,"NDoF")
chisq_pdf.SetParLimits(0,0.,1E4)

def main(toy_pattern_nh,toy_pattern_ih,data_fit_nh,data_fit_ih,tag,fit_degrees):
  
  ROOT.gStyle.SetPalette(56);
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  
  assert(data_fit_nh)
  assert(data_fit_ih)
  
  hierarchies = ["nh","ih"]
  
  patterns = dict([ ("nh",toy_pattern_nh),
                    ("ih",toy_pattern_ih)])
  data_bf  = dict([ ("nh",data_fit_nh),
                    ("ih",data_fit_ih)])
                    
  hists = dict()
  canvs = dict()
  paves = dict()
  lines = dict()
                    
  for hier in hierarchies:
    bf_file = R.TFile(data_bf[hier])
    bf_tree = bf_file.Get("fitnt")
    bf_tree.GetEntry(0)
    bf_chisq = bf_tree.best_fit_params.GetBestFitStatistic()
    
    toy_chain = R.TChain("fitnt")
    toy_chain.Add(patterns[hier])
    
    n_fits = 0
    n_fits_below = 0
    
    h = R.TH1D("chisq_"+hier,"Goodness of Fit #Chi^{2}",nbins,0,chisq_max)
    
    for entry in toy_chain:
      chisq = toy_chain.best_fit_params.GetBestFitStatistic()
      h.Fill(chisq)
      n_fits+=1
      if chisq < bf_chisq:
        n_fits_below+=1
    
    pvalue = 1.-float(n_fits_below)/n_fits
    
    h.Scale(1./n_fits)
    
    if fit_degrees:
      chisq_pdf.SetParameter(0,10)
      hfit = h.Clone()
      hfit.Fit(chisq_pdf)
      fitted_ndof = chisq_pdf.GetParameter(0)
      chisq_pdf.SetRange(hfit.GetXaxis().GetXmin(),hfit.GetXaxis().GetXmax())
      fitted_hist = R.TH1D(chisq_pdf.GetHistogram())
      fitted_hist.SetLineColor(R.kBlue)
    
    ymin = 0
    ymax = h.GetMaximum()*1.1
    h.SetMaximum(ymax)
    
    line = ROOT.TLine(bf_chisq,ymin,bf_chisq,ymax)
    line.SetLineColor(R.kRed)
    line.SetLineWidth(linewidth)
    
    h.SetLineWidth(linewidth)
    h.GetXaxis().SetTitle("#chi^{2}_{GOF}")
    h.GetYaxis().SetTitle("Fraction of fits")
    
    pave = R.TPaveText(0.7,0.7,0.9,0.9,"ndc")
    pave.SetShadowColor(0)
    pave.SetFillColor(0)
    pave.SetFillStyle(0)
    pave.SetLineColor(0)
    pave.AddText("#chi^{2}_{data} = "+to2dp.format(bf_chisq))
    pave.AddText("p-value = "+to2dp.format(pvalue))
    
    if fit_degrees:
      pave2 = R.TPaveText(0.7,0.6,0.9,0.7,"ndc")
      pave2.SetShadowColor(0)
      pave2.SetFillColor(0)
      pave2.SetFillStyle(0)
      pave2.SetLineColor(0)
      pave2.AddText("# DoF = "+to2dp.format(fitted_ndof))
    
    canv = ROOT.TCanvas(hier,hier,1280,1024)
    canv.cd()
    SetBLMargins(canv)
    SetRMargins(canv)
    
    h.Draw()
    if fit_degrees:fitted_hist.Draw("C same")
    line.Draw()
    pave.Draw()
    if fit_degrees: pave2.Draw()
    canv.Print(tag+"_"+hier+"."+image_format)
    
    hists[hier] = h
    canvs[hier] = canv
    paves[hier] = pave
    lines[hier] = line
    
    
  
  
def SetBLMargins(canv):
  canv.SetBottomMargin(0.15)
  canv.SetLeftMargin  (0.12)
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
  
  parser.add_option('--data-bf-nh'  , dest='data_bf_nh' , default=None, help="File containing data best-fit (NH)")
  parser.add_option('--data-bf-ih'  , dest='data_bf_ih' , default=None, help="File containing data best-fit (IH)")
  parser.add_option('-t','--tag'    , dest='tag', default="gof", help="Tag for output")
  parser.add_option('-d','--degrees', dest='degrees',default=False, action="store_true",help='Draw fitted NDOF on plots')

  
  o,args = parser.parse_args()
  
  main(args[0],args[1],
    o.data_bf_nh , o.data_bf_ih ,
    o.tag,o.degrees)

