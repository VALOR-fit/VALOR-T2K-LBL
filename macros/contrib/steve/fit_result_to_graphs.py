#!/usr/bin/env python

import ROOT
import numpy as np
import math
import os
import sys

R=ROOT

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

# Add these as command line arguments
param_x = conv.kPhypMixingMuMu # (disappearance), conv.kPhypMixingMuE for mu-e (appearance)
param_y = conv.kPhypDm41sq

x_min = 1E-4
x_max = 1.
x_npoints = 50

y_min = 1E-2
y_max = 1E2
y_npoints = 50

ndim = 2


eps = 1E-10 # a small number to guarantee our point goes in its bin

def GetLoggedHist(filename):
  chain = R.TChain("valor_dune_nd_ntp")
  chain.Add(filename)
  print "Reading ",chain.GetEntries()," fits"
  assert(chain.GetEntries())
  chain.GetEntry(0)

  hist = R.TH2D("hist","",
    x_npoints,math.log10(x_min),math.log10(x_max),
    y_npoints,math.log10(y_min),math.log10(y_max))

  values = dict()
  for ev in chain:
    fp = ev.best_fit_params
    true = ev.true_params
    x = true.Value(param_x)
    y = true.Value(param_y)
    stat = max(0.,fp.GetBestFitStatistic())

    logx = math.log10(x)
    logy = math.log10(y)
    if logx >= hist.GetXaxis().GetXmax(): logx -= eps
    if logy >= hist.GetYaxis().GetXmax(): logy -= eps
    ibin = hist.FindBin(logx,logy)

    hist.SetBinContent(ibin,stat)

  return hist


def MakeContour(h,dchisq):
  tmp_hist = h.Clone()
  tmp_hist.SetContour(1,np.array([dchisq]))

  canvas = ROOT.TCanvas("c", "c", 0, 0, 700, 500)
  tmp_hist.Draw("cont3same")
  tmp_hist.Draw("CONT LIST same")
  canvas.Update()
  contours = R.gROOT.GetListOfSpecials().FindObject("contours")[0]
  graphs = [g.Clone() for g in contours]

  return graphs

def UnlogGraph(graph):
  n = graph.GetN()
  x = np.array([graph.GetX()[i] for i in xrange(n)])
  y = np.array([graph.GetY()[i] for i in xrange(n)])

  new_graph = R.TGraph(n,10**x,10**y)
  return new_graph

def main(files,tag):

  logged_hists = [GetLoggedHist(f) for f in files]
  canv = R.TCanvas("canv","canv",1200,800)
  canv.SetTopMargin(0.03)
  canv.SetRightMargin(0.03)
  # ~ canv.SetLogx()
  # ~ canv.SetLogy()

  f = R.TFile("contours_as_graphs.root","UPDATE")

  for h in logged_hists:

    h.SetName("%s_log_hist"%(tag))
    h.Write()

    for cl in xrange(conv.kNConfLevels):
      clname = util.ConfLevelAsString(cl)
      dchisq = util.DeltaChisquare(cl,2)

      if (dchisq < 0): continue
      logged_graphs = MakeContour(h,dchisq)
      unlogged_graphs = map(UnlogGraph,logged_graphs)

      for ig in xrange(len(logged_graphs)):
        logged_graphs  [ig].SetName("%s_%s_log_%u"%(tag,clname,ig))
        unlogged_graphs[ig].SetName("%s_%s_lin_%u"%(tag,clname,ig))

        logged_graphs  [ig].Write()
        unlogged_graphs[ig].Write()

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options] patterns",description="Make plots of deltaCP sensitivity")
  parser.add_option('-t','--tag', dest='tag',default="contour",help='The names to give the graphs')

  options,patterns = parser.parse_args()

  main(patterns,options.tag)
