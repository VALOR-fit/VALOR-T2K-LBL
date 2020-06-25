#!/usr/bin/env python

import os
import ROOT as R

R.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

valor = R.valor

def main(config,outfile):
  
  config_path = os.environ["VALOR_EXPT"]+"config/"+config
  
  elements = valor.CompoundFitElements()
  elements.Read(config_path)
  
  syst = elements.FitModel().GetSystematics()
  
  ppl = valor.PhysParamList()
  ppl.SetName("prefit_errors")
  syst.CopyOneSigmaErrorsTo(ppl)
  
  f = R.TFile(outfile,"RECREATE")
  ppl.Write()
  f.Close()
  
  print "Done!"

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Submit many VALOR sensitivity jobs")
  
  parser.add_option('-c', '--config', dest='config', default="fit/dune_ndfit.xml", help="Configuration file")
  parser.add_option('-o', '--outfile', dest='outfile', default="initial_errors.root", help="Output ROOT file")
  
  o,args = parser.parse_args()
  
  main(o.config,o.outfile)
