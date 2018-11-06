#!/usr/bin/env python

import ROOT
import sys
import numpy as np
import argparse

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

def sk_names():
    params = []
    params.append("osc #nu_{e} CC (1Re)")
    params.append("")
    params.append("")
    params.append("#nu_{#mu}/#bar{#nu}_{#mu} CC (1Re)")
    params.append("")
    params.append("")
    params.append("#nu_{e}/#bar{#nu}_{e} CC (1Re)")
    params.append("")
    params.append("")
    params.append("all NC (1Re)")
    params.append("")
    params.append("")

    params.append("#nu_{#mu}/#bar{#nu}_{#mu} CCQE (1Rmu)")
    params.append("")
    params.append("")
    params.append("#nu_{#mu}/#bar{#nu}_{#mu} CCnQE (1Rmu)")
    params.append("#nu_{e}/#bar{#nu}_{e}/sig#nu_{e} CC (1Rmu)")
    params.append("all NC (1Rmu)")

    params.append("osc #nu_{e} CC (1Re); RHC")
    params.append("")
    params.append("")
    params.append("#nu_{#mu}/#bar{#nu}_{#mu} CC (1Re); RHC")
    params.append("")
    params.append("")
    params.append("#nu_{e}/#bar{#nu}_{e} CC (1Re); RHC")
    params.append("")
    params.append("")
    params.append("all NC (1Re); RHC")
    params.append("")
    params.append("")

    params.append("#nu_{#mu}/#bar{#nu}_{#mu} CCQE (1Rmu); RHC")
    params.append("")
    params.append("")
    params.append("#nu_{#mu}/#bar{#nu}_{#mu} CCnQE (1Rmu); RHC")
    params.append("#nu_{e}/#bar{#nu}_{e}/sig#nu_{e} CC (1Rmu); RHC")
    params.append("all NC (1Rmu); RHC")

    params.append("osc nue CC (CC1#pi)")
    params.append("")
    params.append("#nu_{#mu}/#bar{#nu}_{#mu} CC (CC1#pi)")
    params.append("")
    params.append("#nu_{e}/#bar{#nu}_{e} CC (CC1#pi)")
    params.append("")
    params.append("all NC (CC1#pi)")
    params.append("")

    params.append("SK energy scale")

    return params

def plot_sk_matrix(filename, matrixname, outputname):
    file = ROOT.TFile(filename, "READ")

    param_names = sk_names()
    matrix = file.Get(matrixname)
    histo = ROOT.TH2D(matrix)
    histo.GetXaxis().SetLabelSize(0.02)
    histo.GetYaxis().SetLabelSize(0.02)
    for i, name in enumerate(param_names):
        histo.GetYaxis().SetBinLabel(i + 1, param_names[i])
        histo.GetXaxis().SetBinLabel(i + 1, param_names[i])
    histo.LabelsOption("v")
    canvas = ROOT.TCanvas("canvas", "canvas", 10, 10, 1000, 1000)
    canvas.SetLeftMargin(0.2)
    canvas.SetBottomMargin(0.2)
    canvas.SetRightMargin(0.2)
    histo.Draw("COLZ")
    minimum = histo.GetMinimum()
    maximum = histo.GetMaximum()
    norm = max(abs(minimum), abs(maximum))
    histo.GetZaxis().SetRangeUser(-norm, norm)
    histo.GetZaxis().SetTitle("#sqrt{#left|V_{ij}#right|} x Sign#left(V_{ij}#right)")
    histo.GetZaxis().SetTitleOffset(1.8)
    histo.GetZaxis().CenterTitle()
    canvas.SaveAs(outputname)
    file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots the SK matrix')
    parser.add_argument("-f", "--in-filename", dest="in_filename", type=str, help="Input SK matrix filename")
    parser.add_argument("-m", "--object-name", dest="objname", type=str, help="SK matrix object name (skdetfsi, skdet, fsi)")
    parser.add_argument("-o", "--out-filename", dest="out_filename", type=str, help="Output PDF file")
    args = parser.parse_args()

    make_palette()
    plot_sk_matrix(args.in_filename, args.objname, args.out_filename)
