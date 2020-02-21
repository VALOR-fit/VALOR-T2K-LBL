#!/usr/bin/env python

import ROOT
import argparse
import re
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

def flux_name_map(input, plot_raw):
    m = re.search(r"\d+", input)     # Searches for the last digits in the name
    output = ""
    # Only apply a label to the first parameter of a group
    if plot_raw:
        if m.group() == str(0):
            if "SKNuMode" in input:
                output = "FHC "
            else:
                output = "RHC "
            if "Numub" in input:
                output += "#bar{#nu}_{#mu}" + " Flux " + m.group()
            elif "Numu" in input:
                output += "#nu_{#mu}" + " Flux " + m.group()
            elif "Nueb" in input:
                output += "#bar{#nu}_{e}" + " Flux " + m.group()
            else:
                output += "#nu_{e}" + " Flux " + m.group()
    else:
        if m.group() == "00":
            if "bar" in input:
                output = "RHC "
            else:
                output = "FHC "
            if "numubar" in input:
                output += "#bar{#nu}_{#mu}" + " Flux " + m.group()[1]
            elif "numu" in input:
                output += "#nu_{#mu}" + " Flux " + m.group()[1]
            elif "nuebar" in input:
                output += "#bar{#nu}_{e}" + " Flux " + m.group()[1]
            else:
                output += "#nu_{e}" + " Flux " + m.group()[1]
    return output

def make_palette():
    r = np.array([0.0, 1.0, 1.0])
    g = np.array([0.0, 1.0, 0.0])
    b = np.array([1.0, 1.0, 0.0])
    s = np.array([0.0, 0.5, 1.0])
    npoints = len(s)
    ncontours = 255
    ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)

def plot_flux(filename, matrixname, outputname, plot_raw):
    label_size = 0.03
    lmargin = 0.15
    rmargin = 0.25
    bmargin = lmargin
    tmargin = 0.10
    canv_size_x = 1000
    canv_size_y = canv_size_x * (1.0 - bmargin - tmargin)   # Adjust the canvas size to make the drawn matrix appear approximately square

    file = ROOT.TFile(filename, "READ")

    flux_name_tag = "flux"
    if plot_raw:
        flux_name_tag = "NuMode"

    banff_names = file.Get("param_list")
    if "prefit" in matrixname:
        nominal = file.Get("prefit_params")
    else:
        nominal = file.Get("postfit_params")
    banff_matrix = file.Get(matrixname)
    fhc_first_index = -1
    fhc_last_index = -1
    for i,name in enumerate(banff_names):
        if name.GetString().Contains(flux_name_tag):
            if fhc_first_index == -1:
                fhc_first_index = i
            else:
                fhc_last_index = i

    # FHC matrix
    Nfhc = 1 + fhc_last_index - fhc_first_index
    flux_fhc_names = ROOT.TObjArray()
    for i in range(fhc_first_index, fhc_last_index + 1):
        flux_fhc_names.Add(banff_names[i])

    flux_fhc_matrix = ROOT.TMatrixD(Nfhc, Nfhc)
    for x in range(Nfhc):
        for y in range(Nfhc):
            flux_fhc_matrix[y][x] = banff_matrix[y + fhc_first_index][x + fhc_first_index]
    histo = ROOT.TH2D(flux_fhc_matrix)
    for i in range(1,histo.GetNbinsX() + 1):
        for j in range(1, histo.GetNbinsY() + 1):
            val = histo.GetBinContent(i, j)
            val = np.sqrt(np.abs(val)) * np.sign(val)
            histo.SetBinContent(i, j, val)
    for i in range(flux_fhc_names.GetEntries()):
        histo.GetYaxis().SetBinLabel(i + 1, flux_name_map(flux_fhc_names[i].GetName(), plot_raw))
        histo.GetXaxis().SetBinLabel(i + 1, flux_name_map(flux_fhc_names[i].GetName(), plot_raw))
        print flux_fhc_names[i], "Nom:", nominal[i + fhc_first_index], "Err:", np.sqrt(flux_fhc_matrix[i][i])

    histo.LabelsOption("v")
    canvas = ROOT.TCanvas("canvas", "canvas", 10, 10, 1000, int(canv_size_y))
    canvas.SetLeftMargin  (lmargin)
    canvas.SetBottomMargin(bmargin)
    canvas.SetRightMargin (rmargin)
    canvas.SetTopMargin   (tmargin)
    histo.Draw("COLZ")
    minimum = histo.GetMinimum()
    maximum = histo.GetMaximum()
    norm = max(abs(minimum), abs(maximum))
    histo.GetZaxis().SetRangeUser(-norm, norm)
    histo.GetZaxis().SetTitle("#sqrt{#left|V_{ij}#right|} x Sign#left(V_{ij}#right)")
    histo.GetZaxis().SetTitleOffset(1.7)
    histo.GetZaxis().CenterTitle()
    histo.GetXaxis().SetLabelSize(label_size)
    histo.GetYaxis().SetLabelSize(label_size)
    canvas.SaveAs(outputname + ".pdf")
    canvas.SaveAs(outputname + ".eps")
    canvas.SaveAs(outputname + ".png")
    file.Close()

def plot_xsec(filename, matrixname, outputname, plot_raw):
    label_size = 0.03
    lmargin = 0.25
    rmargin = 0.25
    bmargin = lmargin
    tmargin = 0.10
    canv_size_x = 1200
    canv_size_y = canv_size_x * (1.0 - bmargin - tmargin)   # Adjust the canvas size to make the drawn matrix appear approximately square

    file = ROOT.TFile(filename, "READ")

    flux_name_tag = "flux"
    if plot_raw:
        flux_name_tag = "NuMode"

    banff_names = file.Get("param_list")
    if "prefit" in matrixname:
        nominal = file.Get("prefit_params")
    else:
        nominal = file.Get("postfit_params")
    banff_matrix = file.Get(matrixname)
    last_index = -1
    for i,name in enumerate(banff_names):
        if name.GetString().Contains(flux_name_tag):
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
        print banff_names[i + start_index], "Nom:", nominal[i + start_index], "Err:", np.sqrt(banff_matrix[i + start_index][i + start_index])
        if banff_names[i + start_index] == "nue_numu":
            print "nue_numu/nueb_numub correlation:", banff_matrix[i + start_index][i + 1 + start_index]
        elif banff_names[i + start_index] == "nuebar_numubar":
            print "nueb_numub/nue_numu correlation:", banff_matrix[i + start_index][i - 1 + start_index]

    histo.LabelsOption("v")
    canvas = ROOT.TCanvas("canvas", "canvas", 10, 10, 1000, int(canv_size_y))
    canvas.SetLeftMargin  (lmargin)
    canvas.SetBottomMargin(bmargin)
    canvas.SetRightMargin (rmargin)
    canvas.SetTopMargin   (tmargin)
    histo.Draw("COLZ")
    histo.GetZaxis().SetRangeUser(-1,1)
    histo.GetZaxis().SetTitle("#sqrt{#left|V_{ij}#right|} x Sign#left(V_{ij}#right)")
    histo.GetZaxis().SetTitleOffset(1.7)
    histo.GetZaxis().CenterTitle()
    histo.GetXaxis().SetLabelSize(label_size)
    histo.GetYaxis().SetLabelSize(label_size)
    canvas.SaveAs(outputname + ".pdf")
    canvas.SaveAs(outputname + ".eps")
    canvas.SaveAs(outputname + ".png")
    file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the BANFF or xsec matrix")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-x", "--xsec", action="store_true", help="Option to plot the xsec component of a BANFF matrix")
    group.add_argument("-f", "--flux", action="store_true", help="Option to plot the flux component of a BANFF matrix")
    parser.add_argument("filename", type=str, help="The name of the file containing the desired matrix")
    parser.add_argument("matrixname", type=str, help="The ROOT object containing the desired matrix")
    parser.add_argument("output", type=str, help="Name of the PDF output")
    parser.add_argument("-r", "--raw", action="store_true", help="Option to plot the raw BANFF matrix rather than the VALOR format BANFF matrix.")
    args = parser.parse_args()
    make_palette()
    if args.xsec:
        plot_xsec(args.filename, args.matrixname, args.output, args.raw)
    elif args.flux:
        plot_flux(args.filename, args.matrixname, args.output, args.raw)
