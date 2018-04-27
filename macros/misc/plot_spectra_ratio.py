#!/usr/bin/env python

import ROOT
import argparse
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

def get_erec(histo):
    x_axis = histo.GetXaxis()
    Nx = x_axis.GetNbins()
    y_axis = histo.GetYaxis()
    Ny = y_axis.GetNbins()
    y_low = y_axis.GetBinLowEdge(1)
    y_high =  y_axis.GetBinUpEdge(Ny)

    name = histo.GetName() + "_erec"
    bin_edges = np.zeros(Ny + 1)
    for r in range(1, Ny + 1):
        bin_edges[r - 1] = y_axis.GetBinLowEdge(r)
    bin_edges[-1] = y_axis.GetBinUpEdge(Ny)
    output_histo = ROOT.TH1D(name, name, Ny, bin_edges)
    output_histo.SetDirectory(0)
    for r in range(1, Ny + 1):
        total = 0
        for c in range(1, Nx + 1):
            total += histo.GetBinContent(c, r)
        output_histo.SetBinContent(r, total)
    return output_histo

def make_coarse(histo, min=0.099):
    x_axis = histo.GetXaxis()
    Nx = x_axis.GetNbins()
    name = histo.GetName() + "_coarse" + str(min)
    bin_edges = [x_axis.GetBinLowEdge(1)]
    for i in range(1, Nx + 1):
        up_edge = x_axis.GetBinUpEdge(i)
        if up_edge - bin_edges[-1] < min and i != Nx:
            continue
        else:
            bin_edges.append(up_edge)
    bin_edges = np.array(bin_edges)
    output_histo = ROOT.TH1D(name, name, len(bin_edges) - 1, bin_edges)
    output_histo.SetDirectory(0)
    for i in range(1, Nx + 1):
        j = output_histo.FindBin(histo.GetBinCenter(i))
        output_histo.SetBinContent(j, output_histo.GetBinContent(j) + histo.GetBinContent(i))
    print histo.Integral(), output_histo.Integral()
    return output_histo

def details(histo):
    x_axis = histo.GetXaxis()
    y_axis = histo.GetYaxis()
    for r in range(1, histo.GetNbinsY() + 1):
        print "Row", r, y_axis.GetBinLowEdge(r), y_axis.GetBinUpEdge(r)

    for c in range(1, histo.GetNbinsX() + 1):
        print "Col", c, x_axis.GetBinLowEdge(c), x_axis.GetBinUpEdge(c)

def get_ratio(data_filename, asimov_filename):
    samples = [ "1Rmu", "RHC1Rmu", "1Re", "RHC1Re", "MultiRe" ]
    f = ROOT.TFile(data_filename, "READ")
    data_histograms = dict()
    for s in samples:
#        data_histograms[s] = make_coarse(get_erec(f.Get("mock_data_" + s + "_0")))
        data_histograms[s] = get_erec(f.Get("mock_data_" + s + "_0"))
        data_histograms[s].SetDirectory(0)
    f.Close()

    f = ROOT.TFile(asimov_filename, "READ")
    asimov_histograms = dict()
    for s in samples:
#        asimov_histograms[s] = make_coarse(get_erec(f.Get("mock_data_" + s + "_0")))
        asimov_histograms[s] = get_erec(f.Get("mock_data_" + s + "_0"))
        asimov_histograms[s].SetDirectory(0)
    f.Close()

    for s in samples:
        data_histograms[s].Divide(asimov_histograms[s])
    return data_histograms

def plot_ratios(ratio_histograms1, ratio_histograms2, output_filename):
    samples = [ "1Rmu", "RHC1Rmu", "1Re", "RHC1Re", "MultiRe" ]
    titles = [ "FHC #mu-like", "RHC #mu-like", "FHC e-like", "RHC e-like","FHC CC1#pi" ]
    canvas = ROOT.TCanvas("canvas", "canvas", 10, 10, 1000, 1000)
    canvas.SetLeftMargin(0.1)
    canvas.SetBottomMargin(0.1)
    canvas.SetRightMargin(0.2)
    for i, s in enumerate(samples):
        legend = ROOT.TLegend(0.85, 0.7, 0.95, 0.76)
        legend.SetFillColor(0);
        legend.SetLineColor(0);
        legend.SetTextFont(132);
        legend.SetTextSize(0.03);
        xaxis = ratio_histograms1[s].GetXaxis()
        ratio_histograms1[s].SetTitle(titles[i])
        if xaxis.GetBinUpEdge(xaxis.GetNbins()) > 1.2:
            xaxis.SetRangeUser(0, 1.2)
        xaxis.SetTitle("Erec (GeV)")
        xaxis.SetLabelSize(0.02)
        yaxis = ratio_histograms1[s].GetYaxis()
        yaxis.SetTitle("Data / Asimov")
        yaxis.SetLabelSize(0.02)
        ratio_histograms1[s].Draw("axis")
        line = ROOT.TLine(xaxis.GetBinLowEdge(1), 1, xaxis.GetBinUpEdge(xaxis.GetNbins()), 1)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineWidth(1)
        ratio_histograms1[s].SetLineColor(ROOT.kRed)
        ratio_histograms1[s].SetLineWidth(2)
        ratio_histograms2[s].SetLineWidth(2)
        ratio_histograms1[s].Draw("hist")
        ratio_histograms2[s].Draw("histsame")
        line.Draw()
        legend.AddEntry(ratio_histograms1[s], "Run 1-8", "l")
        legend.AddEntry(ratio_histograms2[s], "Run 1-9c", "l")
        legend.Draw()
        if s == "1Rmu":
            canvas.Print(output_filename + "(")
        elif s == "MultiRe":
            canvas.Print(output_filename + ")")
        else:
            canvas.Print(output_filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the ratio of data and Asimov spectra")
    parser.add_argument("-d1", "--data1", type=str, help="The file containing the first data spectra")
    parser.add_argument("-a1", "--asimov1", type=str, help="The file containing the first Asimov spectra")
    parser.add_argument("-d2", "--data2", type=str, help="The file containing the second data spectra")
    parser.add_argument("-a2", "--asimov2", type=str, help="The file containing the second Asimov spectra")
    parser.add_argument("-o", "--output", type=str, help="Output PDF file")
    args = parser.parse_args()
    ratio1 = get_ratio(args.data1, args.asimov1)
    ratio2 = get_ratio(args.data2, args.asimov2)
    plot_ratios(ratio1, ratio2, args.output)

