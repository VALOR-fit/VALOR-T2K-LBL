#!/usr/bin/env python

import ROOT as R
import numpy as np
import math
import os

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


chain = R.TChain("valor_dune_nd_ntp")
chain.Add("test*root")
chain.GetEntry(0)

print "Reading ",chain.GetEntries()," fits"

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

xl = sorted(values.iterkeys())
yl = [values[dcp] for dcp in xl]

canv = R.TCanvas("canv","canv",1200,800)

frame = R.TH1D("f","f",1,-1,1)
frame.GetXaxis().SetTitle("#delta_{CP} / #pi")
frame.GetYaxis().SetTitle("#sigma")
frame.GetYaxis().SetRangeUser(0,10)
frame.Draw()

cdr_file = R.TFile(os.environ[env_var]+"/data/test/CDR_sinsq23_05.root")
cdr_graph = cdr_file.Get("g")
cdr_graph.SetLineWidth(3)
cdr_graph.Draw("C")

mygraph = R.TGraph(len(xl),np.array(xl),np.array(yl))
mygraph.SetLineWidth(3)
mygraph.SetLineColor(R.kRed)
mygraph.Draw("C")

canv.Print("dcp_sens.png")
canv.Print("dcp_sens.eps")
