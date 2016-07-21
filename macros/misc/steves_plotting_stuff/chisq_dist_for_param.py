#!/usr/bin/env python

# Plots chi-square distribution for physics parameters at T2K
# Steve Dennis, June 2014

import ROOT
import math
import array
import numpy as np
import sys
import os

env_var = "T2KOSCANA"
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

config_file = "analysis_configuration.xml"
fit_param_list = "nuosc3f_numudisapp_nuis3fixdcp+nuisance"
stat_fluctuations = True
syst_tweaks = True
nbins = 200

def main(params):
  
  # Set up superKexpectation and read config file
  config =  os.environ[env_var]+"/config/"+config_file
  if not os.path.exists(config):
    raise T2KAnaError("chisq_dist_for_param >> XML file does not exist: "+config)
    
  print "Reading analysis configuration from: "+config;
  skexp = valor.SuperkExpectation()
  skexp.Read(config)
  
  if syst_tweaks: skexp.RandomizeSystematics();
  
  fparm = valor.util.NamedParameterList(fit_param_list)
      
  skobs = valor.SuperkObservation()
  applied_negevents_condition = skobs.Generate(skexp, fparm, stat_fluctuations)
  
  skexp.ResetSystematics();
  
  outfile = ROOT.TFile("chisq_output.root","RECREATE")
  
  hists = []
  
  for param_name in params:
    #~ param = util.PhysParamIdFromString(param_name)
    param = conv.kPhypUndefined
    print "param: ",param
    if (param == conv.kPhypUndefined):
      try:
        param = getattr(conv,param_name)
        print param
      except:
        print "Failed to get param ",param_name
        continue
      
    if not fparm.Exists(param):
      print "Param ",param,"=",util.ParamAsString(param)," not found in FitParamList ",fit_param_list
      continue
      
    fparm.ResetCurrToInitValues()
    
    lower_limit = fparm.MinValue(param)
    upper_limit = fparm.MaxValue(param)
    
    h = ROOT.TH1D(util.ParamAsStringCompact(param)+"_chisq","Chisq of "+util.ParamAsString(param),nbins,lower_limit,upper_limit)
    
    for ibin in xrange(nbins):
      par_value = h.GetBinCenter(ibin+1)
      fparm.Update(param, par_value)
      h.Fill(par_value , skobs.Statistic(skexp,fparm))
      
    hists.append(h)
    h.Write()
  
  canv = ROOT.TCanvas()
  canv.cd()
  for h in hists:
    h.Draw()
    canv.Update()
    canv.Print(h.GetName()+".png")
    
class T2KAnaError(Exception):
  pass
  
if __name__=="__main__":
  main(sys.argv[1:])

