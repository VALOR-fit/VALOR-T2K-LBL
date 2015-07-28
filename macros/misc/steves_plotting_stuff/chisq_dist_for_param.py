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

#~ default_fit_param_list = "nuosc3f_numudisapp_nuis3fixdcp+nuisance"
#~ default_fit_param_list = "nuosc3fbar_numudisapp_constrainall+nuisance"
#~ default_fit_param_list = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance"
default_fit_param_list = "nsi3f_fix4+nuisance"

stat_fluctuations = True
syst_tweaks       = False
nbins = 201

write_bin_spectra = False

to3sf = "{0:.3g}"
python_undefined = 2L*sys.maxint+1

input_mixing12     = 0.304
input_mixing13     = 0.0248
input_mixing23     = 0.527
input_deltacp      = -1.55
input_dm12sq       = 7.53E-5
input_dm23sq       = 2.51E-3
input_mixing12Bar  = 0.304
input_mixing13Bar  = 0.0248
input_mixing23Bar  = 0.527
input_deltacpBar   = -1.55
input_dm12sqBar    = 7.53E-5
input_dm23sqBar    = 2.51E-3
input_epsmutau     = 0

color = ROOT.kBlack

img_format = "png"

def main(params,config_file,fit_param_list,sample):
  
  # Set up superKexpectation and read config file
  config =  os.environ[env_var]+"/config/"+config_file
  if not os.path.exists(config):
    raise T2KAnaError("chisq_dist_for_param >> XML file does not exist: "+config)
    
  print "Reading analysis configuration from: "+config;
  skexp = valor.SuperkExpectation()
  skexp.Read(config)
  
  if syst_tweaks: skexp.RandomizeSystematics();
  
  if fit_param_list is None: fit_param_list = default_fit_param_list
  
  fparm = valor.util.NamedParameterList(fit_param_list)
  
  if fparm.Exists(conv.kPhypMixing12    ) : fparm.SetInitValue(conv.kPhypMixing12   , input_mixing12);
  if fparm.Exists(conv.kPhypMixing13    ) : fparm.SetInitValue(conv.kPhypMixing13   , input_mixing13);
  if fparm.Exists(conv.kPhypMixing23    ) : fparm.SetInitValue(conv.kPhypMixing23   , input_mixing23);
  if fparm.Exists(conv.kPhypCPV         ) : fparm.SetInitValue(conv.kPhypCPV        , input_deltacp);
  if fparm.Exists(conv.kPhypDm12sq      ) : fparm.SetInitValue(conv.kPhypDm12sq     , input_dm12sq);
  if fparm.Exists(conv.kPhypDm23sq      ) : fparm.SetInitValue(conv.kPhypDm23sq     , input_dm23sq);
  if fparm.Exists(conv.kPhypMixing12Bar ) : fparm.SetInitValue(conv.kPhypMixing12Bar, input_mixing12Bar);
  if fparm.Exists(conv.kPhypMixing13Bar ) : fparm.SetInitValue(conv.kPhypMixing13Bar, input_mixing13Bar);
  if fparm.Exists(conv.kPhypMixing23Bar ) : fparm.SetInitValue(conv.kPhypMixing23Bar, input_mixing23Bar);
  if fparm.Exists(conv.kPhypCPVBar      ) : fparm.SetInitValue(conv.kPhypCPVBar     , input_deltacpBar);
  if fparm.Exists(conv.kPhypDm12sqBar   ) : fparm.SetInitValue(conv.kPhypDm12sqBar  , input_dm12sqBar);
  if fparm.Exists(conv.kPhypDm23sqBar   ) : fparm.SetInitValue(conv.kPhypDm23sqBar  , input_dm23sqBar);
  if fparm.Exists(conv.kPhypNSIepsMuTau ) : fparm.SetInitValue(conv.kPhypNSIepsMuTau, input_epsmutau);

      
  skobs = valor.SuperkObservation()
  applied_negevents_condition = skobs.Generate(skexp, fparm, stat_fluctuations)
  
  skexp.ResetSystematics();
  
  outfile = ROOT.TFile("chisq_output.root","RECREATE")
  
  hists = []
  
  canv = ROOT.TCanvas("canv","canv",1280,1024)
  canv.cd()
  
  if len(params)==0:
    params = [util.ParamAsStringSimple(p) for p in xrange(conv.knppmax) if fparm.Exists(p)]
    
  print "Making chisq distribution for "+`len(params)`+" parameters."
  
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
      
    print "Doing param: ",util.ParamAsString(param)
      
    fparm.ResetCurrToInitValues()
    
    lower_limit = fparm.MinValue(param)
    upper_limit = fparm.MaxValue(param)
    
    hname = util.ParamAsStringSimple(param)+"_chisq"
    htitle = "Chisq of "+util.ParamAsStringSimple(param)
    
    h = ROOT.TH1D(hname,htitle,nbins,lower_limit,upper_limit)
    h.SetLineColor(color)
    
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
    h.Write(hname)
  
  canv = ROOT.TCanvas()
  canv.cd()
  for h in hists:
    h.Draw()
    canv.Update()
    canv.Print(h.GetName()+"."+img_format)
    
class T2KAnaError(Exception):
  pass
  
if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of applied systematics")
  
  parser.add_option('-f', '--fit-param-list', dest='fpl', default=None, help="Name of FPL.")
  
  parser.add_option('-s', '--sample' , dest='sample'      , default="1Rmu"                    , help="Sample (1Rmu, 1Re, RHC1Rmu, RHC1Re)")
  
  parser.add_option('-c', '--config'        , dest='config' , default="analysis_configuration.xml", help="Configuration file")
  
  opt,params = parser.parse_args()
  
  main(params,opt.config,opt.fpl,opt.sample)
