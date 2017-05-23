#!/usr/bin/env python

import ROOT
import math
import sys
import os

ROOT.PyConfig.IgnoreCommandLineOptions = True

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)

components = ["NumuToNue","NumuBarToNueBar", "NC", "NueBeam","NumuBeam","NumuBarBeam"]
labels = dict()
labels["NueBeam"] = r"#nu_{e}/#bar{#nu}_{e} intrinsic"
labels["NumuBeam"] = r"#nu_{#mu} intrinsic"
labels["NumuBarBeam"] = r"#bar{#nu}_{#mu} intrinsic"
labels["NumuToNue"] = r"#nu_{#mu}#rightarrow#nu_{e}"
labels["NumuBarToNueBar"] = r"#bar{#nu}_{#mu}#rightarrow#bar{#nu}_{e}"
labels["NC"] = r"NC"

colors = dict()
colors["NueBeam"] = ROOT.kYellow
colors["NumuBeam"] = ROOT.kCyan
colors["NumuBarBeam"] = ROOT.kBlue
colors["NumuToNue"] = ROOT.kWhite
colors["NumuBarToNueBar"] = ROOT.kGreen
colors["NC"] = ROOT.kRed

def make_ereco_spectrum(histos, outfile):
    out = ROOT.TFile(outfile, "RECREATE")
    out.cd()
    total = histos[components[0]].Clone()
    for component in components[1:]:
        total.Add(histos[component])
    total.SetName("totalEreco")
    total.Write("totalEreco")
    out.Close()

def make_ereco_theta_spectrum(histos, totHistos, outfile):
    out = ROOT.TFile(outfile, "RECREATE")
    out.cd()
    total = histos[components[0]].Clone()
    for component in components[1:]:
        total.Add(histos[component])
    total.SetName("totalErecoTheta")
    total.Write("totalErecoTheta")
    for histo in totHistos:
        histo.Write()
    out.Close()

def main(infile, outfile, sample):
    file = ROOT.TFile(infile)
    histos = dict()
    if sample in ["1Rmu", "RHC1Rmu"]:
        for component in components:
            histos[component] = file.Get("hErecoMode_" + component)
        make_ereco_spectrum(histos, outfile)
    else:
        for component in components:
            histos[component] = file.Get("hErecoVsQ2Mode_" + component)
        erecoTot = file.Get("hErecoTot")
        erecoTot.SetName("totalEreco")
        thetaTot = file.Get("hQ2Tot")
        thetaTot.SetName("totalTheta")
        totHistos = [ erecoTot, thetaTot ]
        make_ereco_theta_spectrum(histos, totHistos, outfile)

if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Plot spectra")
    parser.add_option("-f", "--input-file", dest="infile", help='Input file containing the spectra histograms')
    parser.add_option("-o", "--output-file", dest="outfile", help='Output file for spectrum plot')
    parser.add_option("-s", "--sample", dest="sample", help='The sample to plot. 1Re, 1Rmu, MultiRe, RHC1Re, RHC1Rmu.')
    
    v = optparse.Values()
    parser.parse_args(values = v)   # v contains only the values that were specified, not defaults
    o, other = parser.parse_args()
  
    main(o.infile, o.outfile, o.sample)
  
    if other:
        print "Warning - ignored arguments: ", other
