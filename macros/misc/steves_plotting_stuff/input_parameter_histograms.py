#!/usr/bin/env python
import sys
import ROOT 
import os

env_var = "T2KOSCANA"

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
img_format = "png"

valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

ROOT.gROOT.SetBatch(1)

def main(files):
  ch = ROOT.TChain("fitnt")
  for f in files:
    ch.Add(f)
  
  ch.GetEntry(0)
  
  hists = dict()
  params = []
  for param in xrange(conv.knppmax):
    if not ch.input_params.Contains(param): continue
    params.append(param)
  
  min_val = dict([(p, 1E100) for p in params])
  max_val = dict([(p,-1E100) for p in params])
  for entry in ch:
    for param in params:
      val = entry.input_params.Value(param)
      if val < min_val[param]: min_val[param] = val
      if val > max_val[param]: max_val[param] = val
  
  for param in params:
    pname = util.ParamAsStringSimple(param)
    hists[param] = ROOT.TH1D(pname,pname,25,min_val[param],max_val[param])
    
  for entry in ch:
    for param in params:
      hists[param].Fill(entry.input_params.Value(param))
  
  for param in params:
    pname = util.ParamAsStringSimple(param)
    plotname = pname+"_values."+img_format 
    c = ROOT.TCanvas()
    hists[param].Draw()
    c.Update()
    c.Print(plotname)
    
  
if __name__=="__main__":
  main(sys.argv[1:])

