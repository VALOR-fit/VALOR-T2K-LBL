#!/usr/bin/python

import sys
import math
import ROOT
import re
import glob

ROOT.PyConfig.IgnoreCommandLineOptions = True

R = ROOT

R.gROOT.SetBatch(1)

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

files = glob.glob("LBNE_smear_plot_*.root")

output_file = R.TFile("smearing_histograms.root","RECREATE")

for f in files:
  pattern = r"LBNE_smear_plot_(\w*)_(\w*)_ccQE_Q1_(\w*)\.root"
  
  d_flav,d_beam,d_ana = re.match(pattern, f).groups()
  
  if (d_flav=="NC"):
    output_flav = "NC"
  else:
    if d_flav=="numu" : nuflav = conv.kNuMu
    if d_flav=="anumu": nuflav = conv.kNuMuBar
    if d_flav=="nue"  : nuflav = conv.kNuE
    if d_flav=="anue" : nuflav = conv.kNuEBar
    if d_flav=="nut"  : nuflav = conv.kNuTau
    if d_flav=="anut" : nuflav = conv.kNuTauBar
    output_flav = util.FlavourAsString(nuflav)
  
  if d_beam == "nu" :
    beam = conv.kbLBNF_FHC
  elif d_beam == "anu":
    beam = conv.kbLBNF_RHC
    
  if d_ana=="app": sample = "nueCChigh"
  if d_ana=="dis": sample = "numuCChigh"
    
  tf = R.TFile(f)
  canv = tf.GetListOfKeys()[0].ReadObj()
  canv.Draw()
  hists = [h for h in canv.GetListOfPrimitives() if h.ClassName()=="TH2D"]
  assert(len(hists)==1)
  hist = R.TH2D(hists[0])
  
  hist.SetName("smearing_%s_%s_%s"%(util.BeamConfAsString(beam),output_flav,sample))
  
  output_file.cd()
  hist.Write()
  print "Written: ",hist.GetName()
  
  del canv
  
