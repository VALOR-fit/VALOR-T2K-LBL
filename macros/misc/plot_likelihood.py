#!/usr/bin/env python

import ROOT
import argparse
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

def plot_likelihood(in_filename, out_filename):
    in_file = ROOT.TFile(in_filename, "READ")
    canvas = ROOT.TCanvas("likelihood_canvas", "likelihood_canvas", 1024, 768)
    canvas.SetLeftMargin(0.125)
    canvas.SetBottomMargin(0.125)
    canvas.SetRightMargin(0.125)

    likelihood = in_file.Get("likelihood_2D")
    print likelihood.GetMaximum()
    likelihood.GetXaxis().SetTitle("sin^{2}#theta_{23}")
    if "nh" in in_filename or "normal" in in_filename:
        likelihood.GetYaxis().SetTitle("|#Deltam^{2}_{32}| (eV^{2}c^{-4})")
    else:
        likelihood.GetYaxis().SetTitle("|#Deltam^{2}_{13}| (eV^{2}c^{-4})")
    likelihood.GetXaxis().SetTitleOffset(1.2)
    likelihood.GetYaxis().SetTitleOffset(1.8)
    likelihood.Draw("COLZ")
    canvas.Print(out_filename)
    in_file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the likelihood surface")
    parser.add_argument("--input", type=str, help="The root file containing the surface", required=True)
    parser.add_argument("--output", type=str, help="Name of the output PDF file", required=True)
    args = parser.parse_args()
    plot_likelihood(args.input, args.output)
