#!/usr/bin/env python

import ROOT
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

def make_palette():
    r = np.array([0.0, 1.0, 1.0])
    g = np.array([0.0, 1.0, 0.0])
    b = np.array([1.0, 1.0, 0.0])
    s = np.array([0.0, 0.5, 1.0])
    npoints = len(s)
    ncontours = 255
    ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)

def plot_banff(filename, matrixname, outputname):
    file = ROOT.TFile(filename, "READ")

    xsec_names = file.Get("param_list_banff_xsec")
    banff_names = file.Get("param_list")
    nxsec = xsec_names.GetEntries()
    matrix = file.Get(matrixname)
    xsec_matrix = ROOT.TMatrixD(nxsec, nxsec)

    for i in range(xsec_names.GetEntries()):
        x = banff_names.IndexOf(xsec_names[i])
        for j in range(xsec_names.GetEntries()):
            y = banff_names.IndexOf(xsec_names[j])
            xsec_matrix[i][j] = matrix(x, y)
    histo = ROOT.TH2D(xsec_matrix)
    for i in range(1,histo.GetNbinsX() + 1):
        for j in range(1, histo.GetNbinsY() + 1):
            val = histo.GetBinContent(i, j)
            val = np.sqrt(np.abs(val)) * np.sign(val)
            histo.SetBinContent(i, j, val)
    for i in range(xsec_names.GetEntries()):
        histo.GetYaxis().SetBinLabel(i + 1, xsec_names[i].GetName())
        histo.GetXaxis().SetBinLabel(i + 1, xsec_names[i].GetName())
    histo.LabelsOption("v")        
    canvas = ROOT.TCanvas("canvas", "canvas", 10, 10, 1000, 1000)
    canvas.SetLeftMargin(0.2)
    canvas.SetBottomMargin(0.2)
    canvas.SetRightMargin(0.2)
    histo.Draw("COLZ")
    histo.GetZaxis().SetRangeUser(-1,1)
    histo.GetZaxis().SetTitle("#sqrt{#left|V_{ij}#right|} x Sign#left(V_{ij}#right)")
    histo.GetZaxis().SetTitleOffset(1.8)
    histo.GetZaxis().CenterTitle()
    canvas.SaveAs(outputname)
    file.Close()

if __name__ == "__main__":
    make_palette()
    plot_banff(sys.argv[1], sys.argv[2], sys.argv[3])
