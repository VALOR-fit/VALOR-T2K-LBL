#!/usr/bin/env python

import ROOT as R
import numpy as np
import math
import os
import sys

env_var = "VALOR_EXPT";
try:
  R.valor
except:
  R.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C");
  
R.gROOT.SetBatch(1)

R.gStyle.SetOptStat(0)
R.gStyle.SetOptTitle(0)
  
valor = R.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

colors = [R.kBlack, R.kRed, R.kBlue, R.kMagenta+1,  R.kGreen+2, R.kOrange+7, R.kCyan+2,]

param_x = conv.kPhypMixingMuMu
param_y = conv.kPhypDm41sq

x_min = 1E-3
x_max = 1.
x_npoints = 50

y_min = 1E-2
y_max = 1E2
y_npoints = 50

ndim = 2

def GetHist(filename):
  chain = R.TChain("valor_dune_nd_ntp")
  chain.Add(filename)
  print "Reading ",chain.GetEntries()," fits"
  assert(chain.GetEntries())
  chain.GetEntry(0)
  
  # ~ step_x = (math.log10(x_max) - math.log10(x_min))/(x_npoints-1)
  # ~ bin_edge_x = [10**(step_x*i + math.log10(x_min)) for i in xrange(x_npoints)] # Obviously a bit wrong, puts points at lower bin edges
  # ~ step_y = (math.log10(y_max) - math.log10(y_min))/(y_npoints-1)
  # ~ bin_edge_y = [10**(step_y*i + math.log10(y_min)) for i in xrange(y_npoints)] # Obviously a bit wrong, puts points at lower bin edges
  # ~ hist = R.TH2D("hist","",
    # ~ x_npoints,np.array(bin_edge_x),
    # ~ y_npoints,np.array(bin_edge_y))
    
  hist = R.TH2D("hist","",
    x_npoints,math.log10(x_min),math.log10(x_max),
    y_npoints,math.log10(y_min),math.log10(y_max))
    
  values = dict()
  for ev in chain:
    fp = ev.best_fit_params
    true = ev.true_params
    x = true.Value(param_x)
    y = true.Value(param_y)
    stat = math.sqrt(max(0.,fp.GetBestFitStatistic()))
    
    ibin = hist.FindBin(math.log10(x),math.log10(y))
    hist.SetBinContent(ibin,stat)

  return hist,
  # ~ return (hist,hist2)

def main(files,labels):
  
  hists = [GetHist(f) for f in files]
  canv = R.TCanvas("canv","canv",1200,800)
  canv.SetTopMargin(0.03)
  canv.SetRightMargin(0.03)
  # ~ canv.SetLogx()
  # ~ canv.SetLogy()
  
  assert(len(hists)==1)
  hist  = hists[0][0]
  
  hist.GetXaxis().SetTitle("log10(sin^{2} 2 #theta_{#mu #mu})")
  hist.GetYaxis().SetTitle("log10(#Delta m^{2})")
  
  contours = np.array(
    [util.DeltaChisquare(conv.kCL90,ndim),
     util.DeltaChisquare(conv.kCL99,ndim)])
     
  for h in hists:
    h.SetContour(len(contours),contours)
  
  # ~ frame = R.TH1D("f","f",1,0.001,1)

  # ~ frame.Draw()
  
  print
  
  hist.Draw("cont3")
  
  # ~ hist2.SetLineStyle(R.kDashed)
  # ~ hist2.Draw("cont3 same")
  
  # ~ legend = R.TLegend(0.15,0.40,0.45,0.97)
  # ~ for i,g in enumerate(graphs):
    # ~ g.SetLineColor(colors[i%2])
    # ~ g.Draw("C")
    # ~ legend.AddEntry(g,labels[i],"L")
    
    # ~ print labels[i],"\t:","Max =",R.TMath.MaxElement(g.GetN(),g.GetY())
  
  # ~ for l in lines:
    # ~ l.Draw()
    
  # ~ legend.Draw()
    
  canv.Print("sterile_sens.png")
  canv.Print("sterile_sens.eps")
  
if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options] patterns",description="Make plots of deltaCP sensitivity")
  parser.add_option('-l','--labels', dest='labels',default="Postfit",help='The labels for postfit for the legend')
  
  options,patterns = parser.parse_args()
  labels = options.labels.split(",")
  assert(len(labels)==len(patterns))
  
  main(patterns,labels)
