#!/usr/bin/env python

# Script to make plots of true values vs fitted values for DUNE

# Steve Dennis, May 2016

import ROOT
import math
import os
import glob
import numpy as np

ROOT.PyConfig.IgnoreCommandLineOptions = True

R = ROOT

#~ title_offset = 0.2
left_margin   = 0.12
bottom_margin = 0.30
right_margin  = 0.05
top_margin  = 0.05
use_log_plots = False
title_size = 0.08
label_size = 0.07

to3dp = "{0:.3f}"

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

plot_opt = ""

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

def main(config, postfit_fn,output_tag,img_format,entry):
  config_path = os.environ["VALOR_EXPT"]+"config/"+config
  
  elem = valor.CompoundFitElements()
  elem.Read(config_path,"nosplines nomc")
  model = elem.FitModel()
  prefit_cov = model.GetSystematics().CovMatrix()
  prefit_err = util.CovarianceMatrixToErrorVector(prefit_cov)
  
  postfit_chain = ROOT.TChain("valor_dune_nd_ntp")
  postfit_chain.Add(postfit_fn)
  assert(postfit_chain.GetEntries()>entry)
  postfit_chain.GetEntry(entry)
  postfit_cov = R.TMatrixD(postfit_chain.best_fit_params.GetCovMatrix())
  postfit_err = util.CovarianceMatrixToErrorVector(postfit_cov)
  
  fp = postfit_chain.best_fit_params
  nfp = fp.FreeParameterSubset()
  truepar = postfit_chain.true_params
  
  npar = nfp.Size()
  
  h_prefit  = R.TH1D("prefit" ,"Prefit",npar,0,npar)
  h_postfit = R.TH1D("postfit","Postfit",npar,0,npar)
  
  #~ h_prefit.GetXaxis().SetTitle("Parameter")
  h_prefit.GetXaxis().SetTitle("")
  h_prefit.GetYaxis().SetTitle("Value")
  
  for ii in xrange(npar):
    p = nfp.WhichParam(ii)
    trueval    = truepar.Value(p)
    prefit_err = model.GetSystematics().OneSigmaError(p)
    
    fitval = fp.BestFitValue(p)
    fiterr = fp.HesseError(p)
    
    h_prefit.Fill(ii+0.5,trueval)
    h_prefit.SetBinError(ii+1,prefit_err)
    
    h_postfit.Fill(ii+0.5,fitval)
    h_postfit.SetBinError(ii+1,postfit_err[ii])
  
  h_prefit.SetMarkerStyle(conv.kOpenSquare)
  h_prefit.SetMarkerColor(conv.kBlack)
  h_prefit.SetLineColor  (conv.kBlack)
  
  h_postfit.SetMarkerStyle(conv.kFullCircle)
  h_postfit.SetMarkerColor(conv.kRed)
  h_postfit.SetLineColor  (conv.kRed)
  
  h_prefit.GetYaxis().SetRangeUser(0,2)
  label_bins(h_prefit,fp,4) # label every fourth bin
  
  div_lines = get_dividing_lines(fp,0,2)
  for l in div_lines:
    l.SetLineColor(R.kBlue)
  
  legend = R.TLegend(left_margin,0.8-top_margin,left_margin+0.2,1.-top_margin)
  legend.AddEntry(h_prefit ,"Prefit" ,"PE")
  legend.AddEntry(h_postfit,"Postfit","PE")
  
  canv = R.TCanvas("c","c",1280,800)
  canv.SetMargin(left_margin,right_margin,bottom_margin,top_margin)
  h_prefit .Draw("EP")
  h_postfit.Draw("EP same")
  for l in div_lines:
    l.Draw()
  legend.Draw()
  canv.Update()
  canv.Print(output_tag+"."+img_format)
  
  return

class DUNEAnaError(Exception):
  pass

def label_bins(h,fp,div = 1):
  for ii in xrange(fp.Size()):
    if (ii%div==0):
      h.GetXaxis().SetBinLabel(ii+1,fp.NameByIndex(ii))
  return

def get_dividing_lines(fit_param_list, range_y_min, range_y_max):
  div_lines = []
  for ip in xrange(1,fit_param_list.Size()):
    p_last = fit_param_list.WhichParam(ip-1);
    p_this = fit_param_list.WhichParam(ip);
    g_last = util.GetParameterGroup(p_last);
    g_this = util.GetParameterGroup(p_this);
    
    if (g_last!=g_this):
      div_lines.append(R.TLine(ip,range_y_min,ip,range_y_max));
  
  return div_lines
  
  
if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make plots of reduction of fractional error")
  parser.add_option('-c','--config', dest='config'  ,default="fit/dunend_fit.xml",help='The config file we get prefit errors from.')
  parser.add_option('-p','--postfit', dest='postfit',default="fit_ntp.root",help='The postfit input ROOT file')
  parser.add_option('-f','--format', dest='imgformat',default='eps',help='The image format')
  parser.add_option('-o','--output', dest='output',default='true_par_vs_fit',help='Output name')
  parser.add_option('-n','--entry', dest='entry',default=0,type=int,help='The entry in the tree to get')
  
  options,other = parser.parse_args()
  
  if other:
    print "Bad arguments: ",other
    exit(-1)
  
  main(options.config,options.postfit,options.output,options.imgformat,options.entry)
