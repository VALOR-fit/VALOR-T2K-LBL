#!/usr/bin/env python

import ROOT
import argparse
import re
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

def plot_fc_inputs(filename, outputname):
    file = ROOT.TFile(filename, "READ")
    th12 = file.Get("hToyHisto_Mix12_nh")
    th13 = file.Get("hToyHisto_Mix13_nh")
    dm12 = file.Get("hToyHisto_Dm12_nh")

    canvas = ROOT.TCanvas("canvas", "canvas", 10, 10, 1034, 778)
    canvas.SetLeftMargin(0.1)
    canvas.SetBottomMargin(0.15)
    canvas.SetRightMargin(0.1)
    canvas.SetTopMargin(0.1)
    th12.GetXaxis().SetTitle("sin^{2}#theta_{12}")
    th13.GetXaxis().SetTitle("sin^{2}#theta_{13}")
    dm12.GetXaxis().SetTitle("#Deltam^{2}_{12}")
    for i, histo in enumerate([th12, dm12, th13]):
        histo.Draw("")
        histo.GetXaxis().SetTitleOffset(1.5)
        histo.GetXaxis().CenterTitle()
        histo.GetYaxis().SetTitle("#Throws")
        histo.GetYaxis().SetTitleOffset(1.2)
        histo.GetYaxis().CenterTitle()
        if i == 0:
            canvas.Print(outputname + "(")
        elif i == 2:
            canvas.Print(outputname + ")")
        else:
            canvas.Print(outputname)
    file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the FC inputs")
    parser.add_argument("filename", type=str, help="The name of the file containing the FC inputs")
    parser.add_argument("output", type=str, help="Name of the PDF output")
    args = parser.parse_args()
    plot_fc_inputs(args.filename, args.output)
