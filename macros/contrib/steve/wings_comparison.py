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
  R.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")
  
R.gROOT.SetBatch(1)

R.gStyle.SetOptStat(0)
R.gStyle.SetOptTitle(0)
  
valor = R.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

colors = [R.kBlack, R.kRed, R.kBlue, R.kMagenta+1,  R.kGreen+2, R.kOrange+7, R.kCyan+2,]
markers = [R.kFullSquare, R.kFullCircle, R.kFullTriangleUp, R.kFullCross, R.kFullDiamond, R.kFullStar, R.kFullTriangleDown]

ymax = 8.

def GetGraph(filename):
  chain = R.TChain("valor_dune_nd_ntp")
  chain.Add(filename)
  print "Reading ",chain.GetEntries()," fits"
  assert(chain.GetEntries())
  chain.GetEntry(0)

  values = dict()
  for ev in chain:
    fp = ev.best_fit_params
    true = ev.true_params
    dcp = true.Value(conv.kPhypCPV)/math.pi
    stat = math.sqrt(max(0.,fp.GetBestFitStatistic()))
    try:
      values[dcp] = min(values[dcp],stat)
    except KeyError:
      values[dcp] = stat
      
  # Mirror the highest dcp point (for dcp=pi)
  dcp_max = max(values.iterkeys())
  values[dcp_max-2.] = values[dcp_max]

  xl = sorted(values.iterkeys())
  yl = [values[dcp] for dcp in xl]
  
  mygraph = R.TGraph(len(xl),np.array(xl),np.array(yl))
  mygraph.SetLineWidth(3)
  
  return mygraph

def main(files,labels):
  
  graphs = [GetGraph(f) for f in files]
  canv = R.TCanvas("canv","canv",1200,800)
  canv.SetTopMargin(0.03)
  canv.SetRightMargin(0.03)
  canv.SetGridy(1)

  frame = R.TH1D("f","f",1,-1,1)
  frame.GetXaxis().SetTitle("#delta_{CP} / #pi")
  frame.GetYaxis().SetTitle("#sigma")
  frame.GetYaxis().SetRangeUser(0,ymax)
  frame.Draw()

  #~ cdr_file = R.TFile(os.environ[env_var]+"/data/test/CDR_sinsq23_05.root")
  #~ cdr_graph = cdr_file.Get("g")
  #~ cdr_graph.SetLineWidth(3)
  #~ cdr_graph.SetLineColor(R.kBlack)
  #~ cdr_graph.Draw("C")
  
  print
  
  legend = R.TLegend(0.35,0.70,0.78,0.97)
  for i,g in enumerate(graphs):
    g.SetLineColor(colors[i])
    g.Draw("C")
    legend.AddEntry(g,labels[i],"L")
    
    print labels[i],"\t:","Max =",R.TMath.MaxElement(g.GetN(),g.GetY())
    
  legend.Draw()
    
  canv.Print("dcp_sens.png")
  canv.Print("dcp_sens.eps")
  canv.Print("dcp_sens.root")
  
if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options] patterns",description="Make plots of deltaCP sensitivity")
  parser.add_option('-l','--labels', dest='labels',default="Postfit",help='The labels for postfit for the legend')
  
  options,patterns = parser.parse_args()
  labels = options.labels.split(",")
  assert(len(labels)==len(patterns))
  
  main(patterns,labels)
