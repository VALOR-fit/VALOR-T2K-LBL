#!/usr/bin/env python

# This macro produces a plot showing how the true neutrino energy depends on the energy of the pions which decay to
# produce them, and how this also depends on the off-axis angle.
#
# Francis Bench <F.Bench \at liverpool.ac.uk>
# March 2021

import ROOT
import sys
import os

sys.path.append(os.path.abspath(os.environ["VALOR_FMWK"] + "/src/valor/Utilities/"))

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

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


#_______________________________________________________________________
def main(output_filename):
  # Define conventions
  OA = [0.0, 1.0, 2.5, 4.0] # Degrees
  OA_labels = {
    0.0 : "0#circ",
    1.0 : "1#circ",
    2.5 : "2.5#circ",
    4.0 : "4#circ"
  }
  OA_cols = {
    0.0 : 1,
    1.0 : 856,
    2.5 : 628,
    4.0 : 419,
  }
  nOA = len(OA)
  Epi_min =  0.0  # All in GeV
  Epi_max = 10.0
  Enu_min =  0.0
  Enu_max =  4.5

  # Setup canvas
  canv = ROOT.TCanvas("canv", "canv", 800, 800)
  leg = ROOT.TLegend(0.20, 0.73, 0.30, 0.90)
  leg.SetFillStyle(1001)
  leg.SetBorderSize(0)
  leg.SetTextFont(132)
  leg.SetTextSize(0.04)
  draw_opt = "C"

  # Make hists
  hists = []
  for iOA in range(0, nOA):
    hname = "h_" +str(iOA)
    h = ROOT.TH1D(hname, hname, 500, Epi_min, Epi_max)
    for ib in range(1, h.GetNbinsX()+1):
      Epi = h.GetBinCenter(ib)
      thetaOA = OA[iOA]
      Enu = calc_Enu(Epi, thetaOA)
      if Enu == -999:
        # unphysical as E < m
        continue
      h.SetBinContent(ib, Enu)
    h.SetMinimum(Enu_min)
    h.SetMaximum(Enu_max)
    h.SetLineColor(OA_cols[OA[iOA]])
    h.GetXaxis().SetTitle("E_{#pi} (GeV)")
    h.GetYaxis().SetTitle("E_{#nu} (GeV)")
    h.GetYaxis().SetTitleOffset(0.92)
    leg.AddEntry(h, OA_labels[OA[iOA]], "l")
    h.Draw(draw_opt)
    if iOA == 0:
      draw_opt = draw_opt + ",SAME"
    hists.append(h)

  leg.Draw()
  canv.Update()
  canv.SaveAs(output_filename)


#_______________________________________________________________________
# Calculate Enu in GeV
def calc_Enu(Epi, thetaOA_deg):
  mpi = 0.13957039    # pi^+ mass in GeV, PDG2020
  mmu = 0.1056583745  # muon mass in GeV, PDG2020

  if Epi < mpi:
    # unphysical
    return -999

  costheta = ROOT.TMath.Cos( thetaOA_deg * (ROOT.TMath.Pi() / 180.0) )
  ppi = ROOT.TMath.Sqrt(Epi*Epi - mpi*mpi)
  Enu = ( mpi*mpi - mmu*mmu ) / ( 2*( Epi - ppi*costheta ) )
  return Enu


#_______________________________________________________________________
if __name__ == "__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="")

  # Input/Output
  parser.add_option('-o', '--out-name'       , dest='output_name'    , default="", help="Name of the output file")

  o,other = parser.parse_args()

  main(o.output_name)

  if other:
    print "Warning - ignored arguments: ", other
