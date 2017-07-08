#!/usr/bin/env python

import ROOT
import math
import sys
import os

from array import array
import numpy as np

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

def plot_ereco_spectrum(histos, data, outfile):
    canvas = ROOT.TCanvas("canvas","Spectra",10, 10, 1034, 778)
    hs = ROOT.THStack("hs", "Spectrum")
    legend = ROOT.TLegend(0.7, 0.7, 0.95, 0.95);
    for component in components:
        histos[component].SetFillColor(colors[component])
        histos[component].SetMarkerStyle(21)
        histos[component].SetMarkerColor(colors[component])
        hs.Add(histos[component])
        legend.AddEntry(histos[component], labels[component], "f")
    hs.Draw()
    hs.GetXaxis().SetTitle("#nu Reconstructed Energy (GeV)")
    hs.GetYaxis().SetTitle("Number of Events")
    hs.GetXaxis().SetRangeUser(0,3)
    hs.Draw()
    if not data is None:
        data1D = data.ProjectionY()
        eAxis = data1D.GetXaxis()
        nEBins = eAxis.GetNbins()
        eEdges = []
        eCentre = []
        for i in range(1, nEBins + 1):
            eEdges.append(eAxis.GetBinLowEdge(i))
            eCentre.append(eAxis.GetBinCenter(i))
        eEdges.append(eAxis.GetBinUpEdge(i))
        pointsX = []
        pointsY = []
        for e in range(1, nEBins + 1):
            if data1D.GetBinContent(e) > 0:
                pointsX.append(eCentre[e - 1])
                pointsY.append(data1D.GetBinContent(e))
        x = np.array(pointsX)
        y = np.array(pointsY)
        ex = np.zeros(len(pointsX))
        ey = np.sqrt(y)
        graph = ROOT.TGraphErrors(len(pointsX), x, y, ex, ey)
        graph.SetMarkerStyle(21)
        graph.Draw("P")
    legend.Draw()
    canvas.Print(outfile)

def plot_ereco_theta_spectrum(histos, data, outfile):
    canvas = ROOT.TCanvas("canvas","Spectra", 0, 0, 1000, 750)
    total = histos[components[0]].Clone()
    for component in components[1:]:
        total.Add(histos[component])

    thetaAxis = total.GetXaxis()
    nThetaBins = thetaAxis.GetNbins()
    thetaEdges = []
    thetaCentre = []
    for i in range(1, nThetaBins + 1):
        thetaEdges.append(thetaAxis.GetBinLowEdge(i))
        thetaCentre.append(thetaAxis.GetBinCenter(i))
    thetaEdges.append(thetaAxis.GetBinUpEdge(i))

    eAxis = total.GetYaxis()
    nEBins = eAxis.GetNbins()
    eEdges = []
    eCentre = []
    for i in range(1, nEBins + 1):
        eEdges.append(eAxis.GetBinLowEdge(i))
        eCentre.append(eAxis.GetBinCenter(i))
    eEdges.append(eAxis.GetBinUpEdge(i))

    swappedTotal = ROOT.TH2D("hErecoVsQ2Tot","hErecoVsQ2Tot", nEBins, array('d', eEdges), nThetaBins, array('d', thetaEdges))
    for e in range(1, nEBins + 1):
        for theta in range(1, nThetaBins + 1):
            swappedTotal.SetBinContent(e, theta, total.GetBinContent(theta, e))
    swappedTotal.Draw("COLZ")
    swappedTotal.GetYaxis().SetTitle("#theta (degrees)")
    swappedTotal.GetXaxis().SetTitle("#nu Reconstructed Energy (GeV)")
    swappedTotal.GetZaxis().SetTitle("Number of Events")
    swappedTotal.GetYaxis().SetRangeUser(0,180)
    swappedTotal.GetXaxis().SetRangeUser(0,1.25)
    swappedTotal.Draw("COLZ")
    canvas.SetRightMargin(0.2)
    if not data is None:
#        swappedData = ROOT.TH2D("hEvents","hEvents", nEBins, array('d', eEdges), nThetaBins, array('d', thetaEdges))
#        pointsX = []
#        pointsY = []
#        for e in range(1, nEBins + 1):
#            for theta in range(1, nThetaBins + 1):
#                swappedData.SetBinContent(e, theta, data.GetBinContent(theta, e))
#                if data.GetBinContent(theta, e) > 0:
#                    pointsX.append(eCentre[e - 1])
#                    pointsY.append(thetaCentre[theta - 1])
#        graph = ROOT.TGraph(len(pointsX), np.array(pointsX), np.array(pointsY))
#        graph.Draw("P0")
#        overlay = ROOT.TPad("overlay","",0,0,1,1)
#        overlay.SetFillColor(0)
#        overlay.Draw()
#        overlay.cd()
        data.Draw("P0")
#        overlay.Draw()

    canvas.Print(outfile)

def plot_data_ereco_theta_spectrum(histo, outfile):
    canvas = ROOT.TCanvas("canvas","Spectra", 0, 0, 1000, 750)
    total = histo.Clone()

    thetaAxis = total.GetXaxis()
    nThetaBins = thetaAxis.GetNbins()
    thetaEdges = []
    thetaCentre = []
    for i in range(1, nThetaBins + 1):
        thetaEdges.append(thetaAxis.GetBinLowEdge(i))
        thetaCentre.append(thetaAxis.GetBinCenter(i))
    thetaEdges.append(thetaAxis.GetBinUpEdge(i))

    eAxis = total.GetYaxis()
    nEBins = eAxis.GetNbins()
    eEdges = []
    eCentre = []
    for i in range(1, nEBins + 1):
        eEdges.append(eAxis.GetBinLowEdge(i))
        eCentre.append(eAxis.GetBinCenter(i))
    eEdges.append(eAxis.GetBinUpEdge(i))

    swappedTotal = ROOT.TH2D("hErecoVsQ2Tot","hErecoVsQ2Tot", nEBins, array('d', eEdges), nThetaBins, array('d', thetaEdges))
    for e in range(1, nEBins + 1):
        for theta in range(1, nThetaBins + 1):
            swappedTotal.SetBinContent(e, theta, total.GetBinContent(theta, e))
    swappedTotal.Draw("COLZ")
    swappedTotal.GetYaxis().SetTitle("#theta (degrees)")
    swappedTotal.GetXaxis().SetTitle("#nu Reconstructed Energy (GeV)")
    swappedTotal.GetZaxis().SetTitle("Number of Events")
    swappedTotal.GetYaxis().SetRangeUser(0,180)
    swappedTotal.GetXaxis().SetRangeUser(0,1.25)
    swappedTotal.Draw("COLZ")
    canvas.SetRightMargin(0.2)
    canvas.Print(outfile)

def asimov(infile, outfile, datafile, sample):
    file = ROOT.TFile(infile)
    if not datafile == "":
        dfile = ROOT.TFile(datafile)
        data = dfile.Get("mock_data_" + sample + "_0")
        graph = dfile.Get("g" + sample)
    else:
        data = None
        graph = None
    histos = dict()
    if sample in ["1Rmu", "RHC1Rmu"]:
        for component in components:
            histos[component] = file.Get("hErecoMode_" + component)
        plot_ereco_spectrum(histos, data, outfile)
    else:
        for component in components:
            histos[component] = file.Get("hErecoVsQ2Mode_" + component)
        plot_ereco_theta_spectrum(histos, graph, outfile)

def data(infile, outfile, sample):
    file = ROOT.TFile(infile)
    if sample in ["1Rmu", "RHC1Rmu"]:
        histo = file.Get("mock_data_" + sample + "_0")
        plot_data_ereco_spectrum(histo, outfile)
    else:
        histo = file.Get("mock_data_" + sample + "_0")
        plot_data_ereco_theta_spectrum(histo, outfile)

if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Plot spectra")
    parser.add_option("-f", "--input-file", dest="infile", help='Input file containing the spectra histograms')
    parser.add_option("-o", "--output-file", dest="outfile", help='Output file for spectrum plot')
    parser.add_option("-d", "--data-file", dest="datafile", help='Data file for spectrum plot')
    parser.add_option("-s", "--sample", dest="sample", help='The sample to plot. 1Re, 1Rmu, MultiRe, RHC1Re, RHC1Rmu.')
    
    v = optparse.Values()
    parser.parse_args(values = v)   # v contains only the values that were specified, not defaults
    o, other = parser.parse_args()
  
    asimov(o.infile, o.outfile, o.datafile, o.sample)
    #data(o.infile, o.outfile, o.sample)
  
    if other:
        print "Warning - ignored arguments: ", other
