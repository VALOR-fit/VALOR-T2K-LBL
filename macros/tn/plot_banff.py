#!/usr/bin/env python

import ROOT
import argparse
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

    banff_names = file.Get("param_list")
    if "prefit" in matrixname:
        nominal = file.Get("prefit_params")
    else:
        nominal = file.Get("postfit_params")
    banff_matrix = file.Get(matrixname)
    last_index = -1
    for i,name in enumerate(banff_names):
        if name.GetString().Contains("NuMode"):
            last_index = i
    start_index = last_index + 1
    Nremain = banff_names.GetEntries() - start_index
    xsec_names = ROOT.TObjArray(Nremain)
    for i in range(Nremain):
        xsec_names.Add(banff_names[i + start_index])

    xsec_matrix = ROOT.TMatrixD(Nremain, Nremain)
    for x in range(Nremain):
        for y in range(Nremain):
            xsec_matrix[y][x] = banff_matrix[y + start_index][x + start_index]
    histo = ROOT.TH2D(xsec_matrix)
    for i in range(1,histo.GetNbinsX() + 1):
        for j in range(1, histo.GetNbinsY() + 1):
            val = histo.GetBinContent(i, j)
            val = np.sqrt(np.abs(val)) * np.sign(val)
            histo.SetBinContent(i, j, val)
    for i in range(xsec_names.GetEntries()):
        histo.GetYaxis().SetBinLabel(i + 1, xsec_names[i].GetName())
        histo.GetXaxis().SetBinLabel(i + 1, xsec_names[i].GetName())
        print(banff_names[i + start_index], nominal[i + start_index])

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

def plot_xsec(filename, matrixname, outputname):
    file = ROOT.TFile(filename, "READ")

    xsec_names = file.Get("xsec_param_names")
    nxsec = xsec_names.GetEntries()
    xsec_matrix = file.Get(matrixname)
    xsec_nominal = file.Get("xsec_param_prior")

    histo = ROOT.TH2D(xsec_matrix)
    for i in range(1,histo.GetNbinsX() + 1):
        for j in range(1, histo.GetNbinsY() + 1):
            val = histo.GetBinContent(i, j)
            val = np.sqrt(np.abs(val)) * np.sign(val)
            histo.SetBinContent(i, j, val)
    for i in range(xsec_names.GetEntries()):
        histo.GetYaxis().SetBinLabel(i + 1, xsec_names[i].GetName())
        histo.GetXaxis().SetBinLabel(i + 1, xsec_names[i].GetName())
        print(xsec_names[i].GetName(), xsec_nominal[i])
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
    parser = argparse.ArgumentParser(description="Plot the BANFF or xsec matrix")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-x", "--xsec", action="store_true", help="Option to specify reading of a standalone xsec matrix")
    group.add_argument("-b", "--banff", action="store_true", help="Option to specify reading of the full BANFF matrix")
    parser.add_argument("filename", type=str, help="The name of the file containing the desired matrix")
    parser.add_argument("matrixname", type=str, help="The ROOT object containing the desired matrix")
    parser.add_argument("output", type=str, help="Name of the PDF output")
    args = parser.parse_args()
    make_palette()
    if args.banff:
        plot_banff(args.filename, args.matrixname, args.output)
    else:
        plot_xsec(args.filename, args.matrixname, args.output)
