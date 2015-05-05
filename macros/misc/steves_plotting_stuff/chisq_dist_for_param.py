#!/usr/bin/python

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

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()

opt.gOscMixingParamChoice = conv.kNuMixSinsqTheta
opt.SetDefaultParameterValues()

#~ fit_param_list = "nuosc3f_numudisapp_nuis3fixdcp+nuisance"
#~ fit_param_list = "nuosc3fbar_numudisapp_constrainall+nuisance"
fit_param_list = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance"
stat_fluctuations = True
syst_tweaks       = False
nbins = 1000

write_bin_spectra = False

python_undefined = 2L*sys.maxint+1

def main(params,config_file,sample):
  
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
  
  canv = ROOT.TCanvas("canv","canv",1280,1024)
  canv.cd()
  
  for param_name in params:
    param = util.PhysParamIdFromString(param_name)
    
    if (param == conv.kPhypUndefined or param == python_undefined):
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
    
    h = ROOT.TH1D(util.ParamAsStringSimple(param)+"_chisq","Chisq of "+util.ParamAsStringSimple(param),nbins,lower_limit,upper_limit)
    
    for ibin in xrange(nbins):
      par_value = h.GetBinCenter(ibin+1)
      fparm.Update(param, par_value)
      chisq = skobs.Statistic(skexp,fparm)
      h.Fill(par_value , chisq )
      
      if (write_bin_spectra):
        tmp_hist = skexp.CurrPredictedErecoSpectrum(s)
        tmp_hist.Draw()
        canv.Print(util.ParamAsStringSimple(param)+"_chisq_bin_"+(str(ibin).zfill(4))+"_chisq="+str(chisq)+".png")
        del tmp_hist;
      
    hists.append(h)
    h.Write()
  
  canv = ROOT.TCanvas()
  canv.cd()
  for h in hists:
    h.Draw()
    canv.Update()
    canv.Print(h.GetName()+".eps")
    
class T2KAnaError(Exception):
  pass
  
if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of applied systematics")
  
  #~ parser.add_option('-f', '--file', dest='output_file', default="effects_of_tuning.tex", help="Name for output TeX file.")
  
  parser.add_option('-s', '--sample' , dest='sample'      , default="1Rmu"                    , help="Sample (1Rmu, 1Re, RHC1Rmu, RHC1Re)")
  
  parser.add_option('-c', '--config'        , dest='config' , default="analysis_configuration.xml", help="Configuration file")
  
  opt,params = parser.parse_args()
  
  main(params,opt.config,opt.sample)
