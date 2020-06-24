#!/usr/bin/env python

# Script to make plot of relative error magnitude, per parameter for prefit and N postfit results
# Eg fig 39 in DUNE-doc-1291-v2 (docs/internal_notes/dune/v2016a)

# Steve Dennis, May 2016

import ROOT
import math
import os
import glob
import numpy as np

ROOT.PyConfig.IgnoreCommandLineOptions = True

R = ROOT

max_energy = 10.
linewidth = 2
title_offset = 1.1
left_margin   = 0.12
bottom_margin = 0.30
right_margin  = 0.05
top_margin  = 0.05
use_log_plots = False
title_size = 0.08
label_size = 0.07

prefit_color  = ROOT.kGray+1
postfit_color = ROOT.kRed

#~ prefit_style = 3244 # hatched
prefit_style = 1001 # hatched
postfit_style = 1001 # Solid

to3dp = "{0:.3f}"

default_prefit_input_file = (
  "prefit.root")
default_postfit_input_file = (
  "postfit.root")

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kMagenta+1,  ROOT.kGreen+2, ROOT.kOrange+7, ROOT.kCyan+2,]
markers = [R.kFullSquare, R.kFullCircle, R.kFullTriangleUp, R.kFullCross, R.kFullDiamond, R.kFullStar, R.kFullTriangleDown]

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

def main(config, postfit_fns, postfit_labels,output,img_format,entry, spacing):
  config_path = os.environ["VALOR_EXPT"]+"config/"+config
  
  elem = valor.CompoundFitElements()
  elem.Read(config_path,"nosplines nomc")
  model = elem.FitModel()
  syst = model.GetSystematics()
  prefit_cov = syst.CovMatrix()
  prefit_err = util.CovarianceMatrixToErrorVector(prefit_cov)
  
  min_ = min(prefit_err)
  max_ = max(prefit_err)
  
  all_postfit_err = dict()
  
  postfit_files  = postfit_fns.split(",")
  postfit_labels = postfit_labels.split(",")
  
  assert(len(postfit_files)==len(postfit_labels))
  
  for label, postfit_fn in zip(postfit_labels,postfit_files):
    postfit_chain = ROOT.TChain("valor_dune_nd_ntp")
    postfit_chain.Add(postfit_fn)
    assert(postfit_chain.GetEntries())
    postfit_chain.GetEntry(entry)
    postfit_cov = R.TMatrixD(postfit_chain.best_fit_params.GetCovMatrix())
    postfit_err = util.CovarianceMatrixToErrorVector(postfit_cov)
    all_postfit_err[label] = postfit_err
    
    min_ = min(min_,min(postfit_err))
    max_ = max(max_,max(postfit_err))
  
  n_params = prefit_err.GetNrows()
  x = np.array([0.5 + i for i in xrange(n_params)])
  prefit_g  = R.TGraph(n_params,x,prefit_err .GetMatrixArray())
  all_postfit_graphs = dict()
  all_graph_list = [prefit_g]
  for label,postfit_err in all_postfit_err.iteritems():
    g = R.TGraph(n_params,x,postfit_err.GetMatrixArray())
    all_postfit_graphs[label] = g
    all_graph_list.append(g)
  
  for i,g in enumerate(all_graph_list):
    g.SetMarkerColor(colors[i])
    g.SetMarkerStyle(markers[i])
  
  canv = R.TCanvas("c","c",1280,800)
  canv.SetLogy()
  canv.SetMargin(left_margin,right_margin,bottom_margin,top_margin)
                
  miny = min(map(log10_floor,[min_,max_]))
  maxy = max(map(log10_ceil ,[min_,max_]))
  axis_hist = R.TH1D("axis","",n_params,0,n_params)
  axis_hist.GetXaxis().SetRangeUser(0,n_params)
  axis_hist.GetYaxis().SetRangeUser(miny,maxy)
  axis_hist.GetYaxis().SetTitle("Fractional Error")
  axis_hist.GetXaxis().SetTitle("")
  label_bins(axis_hist,postfit_chain.best_fit_params,spacing)
  dividing_lines = get_dividing_lines(postfit_chain.best_fit_params, miny, maxy)
  
  legend = R.TLegend(left_margin,0.95-top_margin,1.-right_margin,1.-top_margin)
  legend.SetNColumns(1+len(postfit_labels))
  legend.AddEntry(prefit_g, "Prefit", "P")
  axis_hist.Draw("axis")
  prefit_g.Draw("P")
  
  for label in postfit_labels:
    postfit_g = all_postfit_graphs[label]
    postfit_g.Draw("P")
    legend.AddEntry(postfit_g,label,"P")
  
  for l in dividing_lines:
    l.Draw()
  legend.Draw()
  canv.Update()
  canv.Print(output+"."+img_format)
  
  return

class DUNEAnaError(Exception):
  pass
  
def log10_floor(x):
  try:
    n = math.log10(x)
  except ValueError:
    n = -5
  return 10**(math.floor(n))
  
def log10_ceil(x):
  try:
    n = math.log10(x)
  except ValueError:
    n = 0
  return 10**(math.ceil(n))
  
def get_hist(rootfile,hname):
  h = rootfile.Get(hname)
  if not h:
    print "Failed to get histogram: ",hname ; exit()
  return ROOT.TH1D( h )
  
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
  
def label_bins(h,fp,div = 1):
  for ii in xrange(fp.Size()):
    if (ii%div==0):
      h.GetXaxis().SetBinLabel(ii+1,fp.NameByIndex(ii))
  return

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make plots of reduction of fractional error")
  parser.add_option('-c','--config', dest='config'  ,default=default_prefit_input_file  ,help='The config file we get prefit errors from.')
  parser.add_option('-p','--postfit', dest='postfit',default=default_postfit_input_file,help='The postfit input ROOT files, comma separated')
  parser.add_option('-l','--postfit-labels', dest='labels',default="Postfit",help='The labels for postfit for the legend')
  parser.add_option('-f','--format', dest='imgformat',default='eps',help='The image format')
  parser.add_option('-o','--output', dest='output',default='prefit_postfit_errors',help='Output name')
  parser.add_option('-n','--entry', dest='entry',default=0,type=int,help='The entry in the tree to get')
  #~ parser.add_option('--skip-fixed',dest='skip_fixed',default=False,action='store_true',help='Do not include fixed parameters in the plot')
  #~ parser.add_option('--skip-det',dest='skip_det',default=False,action='store_true',help='Do not include detector parameters in the plot')
  parser.add_option('--spacing',dest='spacing',type=int,default=4,help='Spacing between labels on the x-axis')
  
  options,other = parser.parse_args()
  
  if other:
    print "Bad arguments: ",other
    exit(-1)
    
  main(options.config,options.postfit,options.labels,options.output,options.imgformat,options.entry,options.spacing)
