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
nr   =  util.NameRegistry.Instance()

colors = [R.kBlack, R.kRed, R.kBlue, R.kMagenta+1,  R.kGreen+2, R.kOrange+7, R.kCyan+2,]
markers = [R.kFullSquare, R.kFullCircle, R.kFullTriangleUp, R.kFullCross, R.kFullDiamond, R.kFullStar, R.kFullTriangleDown]

param = "mixing23"

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
    xval = fp.BestFitValue(nr.PhysParamFromString(param))
    stat = math.sqrt(max(0.,fp.GetBestFitStatistic()))
    try:
      values[xval] = min(values[xval],stat)
    except KeyError:
      values[xval] = stat

  xl = sorted(values.iterkeys())
  yl = [values[xval] for xval in xl]
  
  mygraph = R.TGraph(len(xl),np.array(xl),np.array(yl))
  mygraph.SetLineWidth(3)
  
  return mygraph

def main(files,labels,min_val,max_val,ymax=None):
  
  graphs = [GetGraph(f) for f in files]
  canv = R.TCanvas("canv","canv",1200,800)
  canv.SetTopMargin(0.03)
  canv.SetRightMargin(0.03)
  canv.SetGridy(1)
  
  
  if min_val is None:
    min_val = min( (g.GetX()[0] for g in graphs) )
    
  if max_val is None:
    max_val = max( (g.GetX()[g.GetN()-1] for g in graphs) )
  
  print min_val,max_val
  
  frame = R.TH1D("f","f",1,min_val,max_val)
  frame.GetXaxis().SetTitle(param)
  frame.GetYaxis().SetTitle("#sigma")
  
  if ymax is None:
    for g in graphs:
      ycurr = max( (g.GetY()[i] for i in xrange(g.GetN()) if g.GetX()[i] >= min_val and g.GetX()[i] <= max_val) )
      ymax = max(ymax,ycurr)
  
  draw_legend = len(graphs) > 1
  if draw_legend:
    ymax *= 1.4
  else:
    ymax *= 1.1
    
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
  
  if draw_legend:
    legend.Draw()
    
  canv.Print("param_scan_%s.png" %param)
  canv.Print("param_scan_%s.eps" %param)
  canv.Print("param_scan_%s.root"%param)
  
if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options] patterns",description="Make plots of deltaCP sensitivity")
  parser.add_option('-l','--labels', dest='labels',default="Postfit",help='The labels for postfit for the legend')
  parser.add_option('--min', type=float,dest='min_val',default=None,help='The minimum value on the X axis')
  parser.add_option('--max', type=float,dest='max_val',default=None,help='The maximum value on the X axis')
  
  options,patterns = parser.parse_args()
  labels = options.labels.split(",")
  assert(len(labels)==len(patterns))
  
  main(patterns,labels,options.min_val,options.max_val)
