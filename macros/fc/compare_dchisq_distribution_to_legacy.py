#!/usr/bin/env python

#
# For a single FC grid point, compare the distribution of FC Delta chisq values to that of the legacy analysis
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
arridx = util.ArrIdx.Instance()

python_undefined = 2L*sys.maxint+1
python_undefined_short = 2**32 -1

_dchisq_vec_name_sdk = "Dchi2_toy_vec"
_legacy_tree_name    = "fitnt"

_hist_n_bins     = 100
_hist_dchisq_max = 10


#_______________________________________________________________________

def main(output_tag, in_filename_sdk, in_filename_legacy):

  hist_legacy = ROOT.TH1D("legacy", "legacy", _hist_n_bins, 0, _hist_dchisq_max)
  hist_sdk    = ROOT.TH1D("sdk"   , "sdk"   , _hist_n_bins, 0, _hist_dchisq_max)

  # Read in the legacy Delta chisq values from a tree
  file_legacy = ROOT.TFile(in_filename_legacy, "read")
  if not file_legacy:
    print "ERROR: Could not find file", in_filename_legacy
    sys.exit()
  tree_legacy = file_legacy.Get(_legacy_tree_name)
  n_toys_legacy = tree_legacy.GetEntries()
  for toy in tree_legacy:
    hist_legacy.Fill(toy.deltachisq0BF_glob)
  file_legacy.Close()

  # Read in the SDK Delta chisq values from a vector
  file_sdk = ROOT.TFile(in_filename_sdk, "read")
  if not file_sdk:
    print "ERROR: Could not find file", in_filename_sdk
    sys.exit()
  dchisq_vals_sdk = list(file_sdk.Get(_dchisq_vec_name_sdk))
  n_toys_sdk = len(dchisq_vals_sdk)
  for dchisq in dchisq_vals_sdk:
    hist_sdk.Fill(dchisq)
  file_sdk.Close()

  if n_toys_legacy != n_toys_sdk:
    print "WARNING: # legacy toys =", n_toys_legacy, ", but # SDK toys =", n_toys_sdk

  hist_sdk.GetXaxis().SetTitle("#Delta #chi^{2} FC")
  hist_sdk.GetYaxis().SetTitle("Entries")
  hist_legacy.SetLineColor(2)
  hist_legacy.SetLineStyle(2)

  maximum = max(hist_sdk.GetMaximum(), hist_legacy.GetMaximum())
  hist_sdk.SetMaximum(1.05 * maximum)

  leg = ROOT.TLegend(0.70, 0.80, 0.90, 0.95)
  leg.AddEntry(hist_sdk   , "SDK"   , "l")
  leg.AddEntry(hist_legacy, "Legacy", "l")

  canv = ROOT.TCanvas("Dchisq_dist", "Dchisq_dist", 700, 700)
  hist_sdk.Draw("hist")
  hist_legacy.Draw("hist,same")
  leg.Draw("same")
  canv.SaveAs(output_tag+".pdf")

  print "Done!"

#_______________________________________________________________________
if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Calculate FC Delta chisq values")

  # Input/Output
  parser.add_option('-o', '--out-tag', dest='output_tag'     , default="FC_dchisq_dist_legacy_comparison", help="Tag used to name the output image")
  parser.add_option('-s', '--in-sdk' , dest='in_filename_sdk', default="", help="Input SDK file containing a vector of Delta chisq values")
  parser.add_option('-l', '--in-leg' , dest='in_filename_leg', default="", help="Input T2K legacy file containing a tree of Delta chisq values")

  o,other = parser.parse_args()
  main(o.output_tag, o.in_filename_sdk, o.in_filename_leg)

  if other:
    print "Warning - ignored arguments: ", other
