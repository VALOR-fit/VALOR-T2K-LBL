#!/usr/bin/env python

# Plot the distribution of FC fake data toy throws for 2D mixing23-Dm23sq, and 1D deltacp, mixing12, mixing13 and Dm12sq
#
# Francis Bench
# Nov. 2020

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

_tree_name = "marg_parameter_values"
_hist_nbins_1d = 100
_hist_nbins_2d = 50
_params = [
  "mixing13",
  "mixing23",
  "mixing12",
  "dm23sq",
  "dm12sq",
  "cpv"
]
_min_mixing23 = 0.3
_max_mixing23 = 0.7
_min_Dmsq = 2.25e-3
_max_Dmsq = 2.70e-3
_margin_l = 0.12
_margin_r = 0.14
_margin_t = 0.1
_margin_b = 0.1
_title_size = 0.05
_label_size = 0.04
_title_offset = 0.9

#_______________________________________________________________________

def main(output_tag, in_tag_nh):

  if not util.plotting.SetupCustomPalette("valor_inverted_dark_body_radiator"):
    print "WARNING: the custom color palette could not be set up. Reverting to default"


  # Get both normal and inverted ordering trees
  fnh = ROOT.TFile(in_tag_nh, "read")
  if not fnh or not fnh.IsOpen():
    sys.exit("ERROR: Could not find file \'{}\'".format(in_tag_nh))
  tnh = fnh.Get(_tree_name)
  if not tnh or type(tnh) == ROOT.TObject:
    sys.exit("ERROR: could not get tree named \'{}\' from NH file".format(_tree_name))

  fih_name = in_tag_nh.replace("normal", "inverted")
  fih = ROOT.TFile(fih_name, "read")
  if not fih or not fih.IsOpen():
    sys.exit("ERROR: Could not find file \'{}\'".format(fih_name))
  tih = fih.Get(_tree_name)
  if not tih or type(tih) == ROOT.TObject:
    sys.exit("ERROR: could not get tree named \'{}\' from IH file".format(_tree_name))

  # Loop over all entries to get the max and min param values
  max_osc = {}
  for param in _params:
    max_osc[param] = -999.0
  min_osc = {}
  for param in _params:
    min_osc[param] = 999.0
  assert(tnh.GetEntries() == tih.GetEntries())
  nev = tnh.GetEntries()
  for iev in range(0, nev):
    tnh.GetEntry(iev)
    tih.GetEntry(iev)
    for param in _params:
      val_nh = getattr(tnh, param)
      max_osc[param] = max(max_osc[param], val_nh)
      min_osc[param] = min(min_osc[param], val_nh)
      val_ih = getattr(tih, param)
      max_osc[param] = max(max_osc[param], val_ih)
      min_osc[param] = min(min_osc[param], val_ih)

  # Make hists
  hists_osc_1d = {}
  for param in _params:
    if param == "mixing23" or param == "dm23sq":
      continue
    hists_osc_1d[param] = ROOT.TH1D(param, param, _hist_nbins_1d, 0.9*min_osc[param], 1.1*max_osc[param])
    hists_osc_1d[param].GetXaxis().SetTitle(util.RawParamAsStringRooTex(util.PhysParamIdFromRawString(param)))
    hists_osc_1d[param].GetYaxis().SetTitle("Fraction of Pseudo-Experiments")
    format_hist(hists_osc_1d[param])
  hist_atmos_nh = ROOT.TH2D("atmos_nh", "atmos_nh", _hist_nbins_2d, _min_mixing23, _max_mixing23, _hist_nbins_2d, _min_Dmsq, _max_Dmsq)
  hist_atmos_ih = ROOT.TH2D("atmos_ih", "atmos_ih", _hist_nbins_2d, _min_mixing23, _max_mixing23, _hist_nbins_2d, _min_Dmsq, _max_Dmsq)
  hist_atmos_nh.GetXaxis().SetTitle("sin^{2} #theta_{23}")
  hist_atmos_nh.GetYaxis().SetTitle("#Delta m^{2}_{23} (eV^{2} / c^{4})")
  hist_atmos_nh.GetZaxis().SetTitle("Fraction of Pseudo-Experiments")
  format_hist(hist_atmos_nh)
  hist_atmos_ih.GetXaxis().SetTitle("sin^{2} #theta_{23}")
  hist_atmos_ih.GetYaxis().SetTitle("|#Delta m^{2}_{31}| (eV^{2} / c^{4})")
  hist_atmos_ih.GetZaxis().SetTitle("Fraction of Pseudo-Experiments")
  format_hist(hist_atmos_ih)

  # Fill hists
  for iev in range(0, nev):
    tnh.GetEntry(iev)
    tih.GetEntry(iev)
    for param in _params:
      if param == "mixing23" or param == "dm23sq":
        continue
      val_nh = getattr(tnh, param)
      hists_osc_1d[param].Fill(val_nh)
    val_mixing23 = getattr(tnh, "mixing23")
    val_Dmsq_nh  = getattr(tnh, "dm23sq")
    val_Dmsq_ih  = getattr(tih, "dm23sq")
    hist_atmos_nh.Fill(val_mixing23, val_Dmsq_nh)
    hist_atmos_ih.Fill(val_mixing23, val_Dmsq_ih)

  # Normalise
  for param in _params:
    if param == "mixing23" or param == "dm23sq":
      continue
    hists_osc_1d[param].Scale(1.0 / hists_osc_1d[param].Integral())
  hist_atmos_nh.Scale(1.0 / hist_atmos_nh.Integral())
  hist_atmos_ih.Scale(1.0 / hist_atmos_ih.Integral())

  # Draw
  canv = ROOT.TCanvas("canv", "canv", 800, 600)
  ROOT.TGaxis.SetMaxDigits(3)
  canv.SetMargin(_margin_l, _margin_r, _margin_t, _margin_b)
  hist_atmos_nh.Draw("colz")
  canv.SaveAs(output_tag + ".pdf(")
  hist_atmos_ih.Draw("colz")
  canv.SaveAs(output_tag + ".pdf")
  canv.SetMargin(_margin_l, _margin_l, _margin_t, _margin_b)
  for ip, param in enumerate(_params):
    if param == "mixing23" or param == "dm23sq":
      continue
    hists_osc_1d[param].Draw("hist")
    if ip == len(_params)-1:
      canv.SaveAs(output_tag + ".pdf)")
    else:
      canv.SaveAs(output_tag + ".pdf")

#_______________________________________________________________________
def format_hist(h):
  h.GetZaxis().SetLabelSize(_label_size)
  h.GetYaxis().SetLabelSize(_label_size)
  h.GetXaxis().SetLabelSize(_label_size)
  h.GetZaxis().SetTitleSize(_title_size)
  h.GetYaxis().SetTitleSize(_title_size)
  h.GetXaxis().SetTitleSize(_title_size)
  h.GetZaxis().SetTitleOffset(_title_offset)
  h.GetYaxis().SetTitleOffset(_title_offset)
  h.GetXaxis().SetTitleOffset(_title_offset)

#_______________________________________________________________________
if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Calculate FC critical values for a single mass hierarchy")

  # Input/Output
  parser.add_option('-o', '--out-tag', dest='output_tag', default="FC_toy_throws"                                                    , help="Tag used to name the output file(s)")
  parser.add_option('-i', '--in-tag' , dest='in_tag'    , default="run1-10_FC_fake_data_normal_cpv_reactor_0_1_40000_toy_throws.root", help="Normal ordering input file name")

  o,other = parser.parse_args()

  main(o.output_tag, o.in_tag)

  if other:
    print "Warning - ignored arguments: ", other
