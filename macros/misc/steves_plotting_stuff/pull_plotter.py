#!/usr/bin/env python

# Plots pulls
# Steve Dennis, 2nd June, 2014 

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
  
ROOT.gROOT.SetBatch(1)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

opt.gOscMixingParamChoice = conv.kNuMixSinsqTheta
opt.SetDefaultParameterValues()

image_format = "png"
n_x = 2
n_y = 5
  
def main(files, treename="fitnt",make_combined=True):
  chain = ROOT.TChain(treename)
  for filename in files:
    chain.Add(filename)
  
  chain.GetEntry(0)
  
  fpl = chain.best_fit_params
  
  params = []
  pull_hists = dict()
  
  nbins = int( 2 * math.sqrt(chain.GetEntries()))
  
  for iparam in xrange(fpl.Size()):
    p = fpl.WhichParam(iparam)
    if fpl.IsFixed(p): continue
    name = util.ParamAsStringCompact(p)
    hist = ROOT.TH1D("pull_"+name,"Pull of "+name,nbins,-5,5)
    hist.SetLineColor(ROOT.kBlack)
    hist.SetLineWidth(5)
    
    pull_hists[p] = hist
    params.append(p)
  
  for event in chain:
    for p,hist in pull_hists.iteritems():
      starting      = chain.input_params.Value(p)
      bestfit       = chain.best_fit_params.BestFitValue(p)
      postfit_error = chain.best_fit_params.HesseError(p)
      if postfit_error == 0.: continue
      
      hist.Fill( (bestfit-starting)/ postfit_error)
  
  ROOT.gStyle.SetOptStat(111110)
  canv = ROOT.TCanvas("c","c",1280,720)
  for hist in pull_hists.itervalues():
    canv.cd()
    hist.Draw()
    canv.Print(hist.GetName()+"."+image_format)
    
  if make_combined:
    big_canv = ROOT.TCanvas("big_canv","big_canv",2000,2400)
    
    n_per_page = n_x * n_y
    
    for ipage in xrange( int(math.ceil(float(len(params))/n_per_page)) ):
      big_canv.Divide(n_x,n_y)
      for iplot in xrange(n_per_page):
        i = ipage * n_per_page + iplot
        if i > len(params)-1: continue
        p = params[i]
        hist = pull_hists[p]
        big_canv.cd(iplot+1)
        hist.Draw()
      big_canv.Print("combined_"+str(ipage)+"."+image_format)
      big_canv.Clear()


if __name__=="__main__":
  main(sys.argv[1:])

