#!/usr/bin/env python

import ROOT
import math
import array
import sys
import os

ROOT.PyConfig.IgnoreCommandLineOptions = True

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C");
  
ROOT.gROOT.SetBatch(1)
  
valor   = ROOT.valor
util    = valor.util
conv    = valor.conv
opt     = valor.Opt.Instance()
env     = valor.Env.Instance()
rnd     = valor.Rnd.Instance()
arridx  = util.ArrIdx.Instance()

def main(files):
  splines = dict()
  for f in files:
    tf = ROOT.TFile(f)
    for k in tf.GetListOfKeys():
      spl = k.ReadObj()
      if "pline" not in spl.ClassName(): continue
      splines[spl.GetName()] = spl
  
  print "Read",len(splines)," splines"
  
  n_written = 0
  n_files = 0
  for ib in xrange(conv.kNBeamConf):
    b = util.BeamConfFromInt(ib)
    bname = util.BeamConfAsString(b)
    for ismp in xrange(conv.kMaxNSamples):
      #~ s = util.SampleFromInt(ismp)
      s = ismp
      sname = util.SampleAsString(s)
      
      relevant_splines = [splines[k] for k in splines.iterkeys() if bname in k and sname in k]
      
      if not relevant_splines: continue
      
      output_file = ROOT.TFile("splines_"+bname+"_"+sname+".root","RECREATE")
      print "Writing file ",output_file.GetName()
      n_files+=1
      for spl in relevant_splines:
        spl.Write()
        n_written+=1
      output_file.Close()
  
  print "Wrote",n_written,"splines to",n_files," files."
  return

if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Merge DUNE spline files")

  o,files = parser.parse_args()
  
  main(files)
