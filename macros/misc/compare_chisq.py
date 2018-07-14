#!/usr/bin/env python

import ROOT
import argparse
import re
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

def make_palette():
    r = np.array([0.0, 1.0, 1.0])
    g = np.array([0.0, 1.0, 0.0])
    b = np.array([1.0, 1.0, 0.0])
    s = np.array([0.0, 0.5, 1.0])
    npoints = len(s)
    ncontours = 255
    ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)

def get_chisq(filename, obj):
    file = ROOT.TFile(filename, "READ")
    histo = file.Get(obj)
    histo.SetDirectory(0)
    file.Close()
    return histo

def plot_difference(histo1, histo2, outputname, fractional=False):
    canvas = ROOT.TCanvas("canvas", "canvas", 10, 10, 1000, 1000)
    canvas.SetLeftMargin(0.2)
    canvas.SetBottomMargin(0.2)
    canvas.SetRightMargin(0.2)
    diff_histo = histo1.Clone("diff_histo")
    for i in range(1, diff_histo.GetXaxis().GetNbins() + 1):
        for j in range(1, diff_histo.GetYaxis().GetNbins() + 1):
            val1 = histo1.GetBinContent(i, j)
            val2 = histo2.GetBinContent(i, j)
            if not fractional:
                diff = val1 - val2
            else:
                if val1 != 0:
                    diff = (val1 - val2) / val1
                else:
                    diff = -val2
            diff_histo.SetBinContent(i, j, diff)
    if not fractional:
        diff_histo.GetZaxis().SetRangeUser(-20, 20)
    else:
        diff_histo.GetZaxis().SetRangeUser(-1, 1)
    diff_histo.GetZaxis().SetTitle("#sqrt{#left|V_{ij}#right|} x Sign#left(V_{ij}#right)")
    diff_histo.GetZaxis().SetTitleOffset(1.5)
    diff_histo.GetZaxis().CenterTitle()
    diff_histo.GetXaxis().SetLabelSize(0.02)
    diff_histo.GetYaxis().SetLabelSize(0.02)
    diff_histo.GetXaxis().SetTitleOffset(1.0)
    diff_histo.GetYaxis().SetTitleOffset(1.5)
    diff_histo.Draw("COLZ")
    canvas.SaveAs(outputname + ".pdf")
    canvas.SaveAs(outputname + ".eps")
    canvas.SaveAs(outputname + ".png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the BANFF or xsec matrix")
    parser.add_argument("-v", "--valor", type=str, help="The VALOR contour file")
    parser.add_argument("-p", "--ptheta", type=str, help="The p-theta contour file")
    parser.add_argument("-o", "--output", type=str, help="Name of the PDF output")
    parser.add_argument("--fractional", action="store_true", help="Compute fractional difference")
    args = parser.parse_args()
    make_palette()
    valor_histo = get_chisq(args.valor, "dchisq_map_2D")
    ptheta_histo = get_chisq(args.ptheta, "Asimov_2D")
    plot_difference(valor_histo, ptheta_histo, args.output, args.fractional)
