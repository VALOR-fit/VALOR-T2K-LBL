#!/usr/bin/env python

#
# For each FC point (for fixed FC hierarchy):
#   Collect the FC Delta chisq value for each FC toy.
#   Use this calculate critical values
#
# Francis Bench <F.Bench \at liverpool.ac.uk>
# August 2019

import ROOT
import sys
import os

sys.path.append(os.path.abspath(os.environ["VALOR_FMWK"] + "/src/valor/Utilities/"))

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

python_undefined = 2L*sys.maxint+1
python_undefined_short = 2**32 -1

_vec_name  = "Dchi2_toy_vec"
_hist_name = "Dchi2_hist"

#_______________________________________________________________________

def main(in_filenames, out_filename):

  # First iterate over the input file to get the limits of the combined hist
  lim_low  =  999
  lim_high = -999
  nbins    = -1
  for filename in in_filenames:
    in_file = ROOT.TFile(filename, "read")
    hist = in_file.Get(_hist_name)
    lim_low  = min(lim_low , hist.GetXaxis().GetBinLowEdge(1))
    lim_high = max(lim_high, hist.GetXaxis().GetBinLowEdge(hist.GetNbinsX()) + hist.GetXaxis().GetBinWidth(1))
    nbins    = max(nbins, hist.GetNbinsX())
    in_file.Close()
    del in_file

  vec_combined = ROOT.vector('double')()
  hist_combined = ROOT.TH1D("hist_combined", "hist_combined", nbins, lim_low, lim_high)

  for filename in in_filenames:
    print "Reading from file", filename
    in_file = ROOT.TFile(filename, "read")

    dchisq_vals = list(in_file.Get(_vec_name))
    for dchisq_val in dchisq_vals:
      vec_combined.push_back(dchisq_val)
      hist_combined.Fill(dchisq_val)

    in_file.Close()
    del in_file

  print "Saving to file", out_filename
  outfile = ROOT.TFile(out_filename, "recreate")
  outfile.WriteObject(vec_combined , _vec_name)
  outfile.WriteObject(hist_combined, _hist_name)
  outfile.Close()
  del outfile

  print "Done!"

#_______________________________________________________________________
if __name__=="__main__":
  import argparse

  parser = argparse.ArgumentParser(description="Index splines")
  parser.add_argument("-f", "--filenames", dest="in_filenames", nargs="+", help="Input ROOT filenames containing dchisq vectors to combine.")
  parser.add_argument("-o", "--output"   , dest="out_filename"           , help="Output ROOT filename")

  args = parser.parse_args()

  main(args.in_filenames, args.out_filename)
