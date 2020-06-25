#!/usr/bin/python

# Make parameter matrices using spectra made with varied universes
# These come from eg write_sbn_spectra_using_universes
# Steve Dennis
# November 2019

import sys
import math
import ROOT
import re
import glob
import numpy as np

ROOT.PyConfig.IgnoreCommandLineOptions = True

R = ROOT

R.gROOT.SetBatch(1)

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C");
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

def main(files,config,output,input_tag,nominal,kinematic,grouping):
  
  elements = valor.CompoundFitElements()
  elements.Read(util.CompleteConfigPath(config),"nomc nosplines")
  
  kv = util.KineVarFromString(kinematic)
  is_true = util.KineVarIsTrue(kv)
  is_reco = util.KineVarIsReco(kv)
  
  # Calculate the number of entries that we'll put in the matrix
  row_size = 0
  model_nbins = []
  for model in elements.FitModel().GetSubModels():
    full_binning = model.GetBinning()
    binning = None
    if (is_true):
      binning = full_binning.TrueSmpBinning(model.GetReactionModes()[0]) # only take the binning from the first mode
    if (is_reco):
      binning = full_binning.RecoOptBinning()
      
    n_bins = binning.NBinsByKV(kv)
    model_nbins.append(n_bins)
    # ~ print row_size,n_bins
    
    row_size += n_bins
  
  print "Our covariance will have",row_size,"elements."
  
  # Make the TPrincipal we'll use to actually get covariances
  tp = R.TPrincipal(row_size)
  
  # Loop over each file (containing a universe each)
  for ifile,f in enumerate(files):
    myfile = R.TFile(f)
    
    row = np.array([0.]*row_size)
    ientry = 0
    for imodel,model in enumerate(elements.FitModel().GetSubModels()):
      first_entry_for_sample = ientry
      for mode in model.GetReactionModes():
        if mode < 0 or mode > 1E5: continue
        ientry = first_entry_for_sample
        histname = "%s_%s_%s_%s_%s"%(input_tag,"pdfcs_ob_rf",model.Name(),kinematic,util.ModeAsString(mode))
        hist = myfile.Get(histname)
        if not hist:
          print "Failed to get",histname,"from",f
          exit(-1)
        assert(model_nbins[imodel] == hist.GetXaxis().GetNbins())
        for ibin in xrange(hist.GetXaxis().GetNbins()):
          # We're summing over all reaction modes here
          row[ientry] += hist.GetBinContent(ibin+1)
          ientry+=1
        #end loop ibin
      #end loop mode
    #end loop model
    assert(ientry == row_size)
    print "Adding row to covariance",ifile
    tp.AddRow(row)
  #end loop files
  
  cov = R.TMatrixD(util.LowerTriangularMatrixToSymMatrix(tp.GetCovarianceMatrix()))
  mean = tp.GetMeanValues()
  frac_cov = util.CovarianceMatrixToFractionalCovarianceMatrix(cov,mean)
  corr = util.CovarianceMatrixToCorrelationMatrix(frac_cov)
  err  = util.CovarianceMatrixToErrorVector(frac_cov)
  
  # Do central values here if a nominal file is specified
  # ~ central_values = 
  
  outfile = R.TFile(output,"RECREATE")
  frac_cov.Write("frac_cov")
  corr.Write("corr")
  
  print "Done"
  return
  

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options] patterns",description="Make covariance matrix using universes")
  parser.add_option('-c','--config', dest='config',default="analysis/default.xml",help='The config file')
  parser.add_option('-o','--output', dest='output',default="covariance.root",help='The output file to write')
  parser.add_option('-k','--kinematic',dest='kinematic',default="etrue",help="The kinematic variable to parameterise in")
  parser.add_option('-n','--nominal',dest='nominal',default=None,help="If a nominal spectrum is given, we can write the central values too [NOT YET IMPLEMENTED")
  parser.add_option('-g','--grouping',dest='grouping',default="all",help="How to group modes. Options all,flav")
  parser.add_option('--input-tag',dest='input_tag',default="valor",help="The first part of the histogram name")
  parser.add_option('--write-xml',dest='write_xml',action="store_true",default=False,help="Write parameter XML definitions [NOT YET IMPLEMENTED]")
  
  options,files = parser.parse_args()
  
  main(files,options.config,options.output,options.input_tag,options.kinematic,options.nominal,options.grouping)
