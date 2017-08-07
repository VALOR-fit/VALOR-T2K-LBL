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
colors["NumuToNue"] = ROOT.kMagenta
colors["NumuBarToNueBar"] = ROOT.kGreen
colors["NC"] = ROOT.kRed

def combine_elike(fhc_1re, fhc_cc1pi):
    sum_histo = fhc_1re.Clone()
    sum_histo.Add(fhc_cc1pi)
    return sum_histo

def negate_histo(histo):
    new_histo = histo.Clone()
    axis = new_histo.GetXaxis()
    nbins = axis.GetNbins()
    for i in range(1, nbins + 1):
        new_histo.SetBinContent(i, -new_histo.GetBinContent(i))
    return new_histo

def process_data(data, negate=False):
    data1D = data.ProjectionY()
    data1D = data1D.Rebin(2)
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
    if negate:
        y = -y
    graph = ROOT.TGraphErrors(len(pointsX), x, y, ex, ey)
    graph.SetMarkerStyle(21)
    return graph

def plot_bestfit_dcp0_spectrum(fhc_mc, rhc_mc, fhc_data, rhc_data, outfile):
    canvas = ROOT.TCanvas("canvas","Spectra",10, 10, 1034, 778)
    hs = ROOT.THStack("hs", "Spectrum")
    legend = ROOT.TLegend(0.7, 0.73, 0.8, 0.8);
    histos = [fhc_mc, rhc_mc, fhc_data, rhc_data]
    min = histos[0].GetMinimum()
    max = histos[0].GetMaximum()
    for i, histo in enumerate(histos):
        if i < 2:
            histo.SetLineWidth(2)
            histo.SetLineColor(ROOT.kRed)
            histo.SetMarkerColor(ROOT.kRed)
            if histo.GetMinimum() < min:
                min = histo.GetMinimum()
            if histo.GetMaximum() > max:
                max = histo.GetMaximum()
        else:
            histo.SetMarkerColor(ROOT.kBlack)
            histo.SetLineColor(ROOT.kBlack)
            high = ROOT.TMath.MaxElement(histo.GetN(), histo.GetY())
            low = ROOT.TMath.MinElement(histo.GetN(), histo.GetY())
            low -= np.sqrt(np.abs(low))
            high += np.sqrt(np.abs(high))
            if low < min:
                min = low
            if high > max:
                max = high
    min *= 1.05
    max *= 1.05
    histos[0].Draw()
    axis = histos[0].GetXaxis()
    title = ROOT.TLatex(axis.GetBinLowEdge(7), 1.01 * max, "T2K Run1-8 Preliminary")
    title.SetTextSize(0.05)
    title.Draw()
    line = ROOT.TLine(0, 0, axis.GetBinUpEdge(axis.GetNbins()), 0)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(2)
    line.Draw("same")
    histos[0].SetMinimum(min)
    histos[0].SetMaximum(max)
    legend.AddEntry(fhc_mc, "#delta_{CP} = 0", "l")
    legend.AddEntry(fhc_data, "Data", "lp")
    histos[0].GetXaxis().SetTitle("#nu Reconstructed Energy (GeV)")
    histos[0].GetYaxis().SetTitle("Number of Events")
    for i, histo in enumerate(histos):
        if i < 2:
            histo.Draw("same")
        else:
            histo.Draw("P")
    legend.Draw()
    canvas.SetTopMargin(0.1)
    canvas.Print(outfile)

def make_asimovb_variation_spectrum(outfile):
    files = []
    histos = []
    rhc_histos = []
    legend = ROOT.TLegend(0.7, 0.7, 0.95, 0.95)
    rhc_legend = ROOT.TLegend(0.7, 0.7, 0.95, 0.95)

    file_1Re_nh_dcp_0 = ROOT.TFile("spectra_b_nh_dcp_0_1Re.root", "READ")
    files.append(file_1Re_nh_dcp_0)
    h1Re_nh_dcp_0 = file_1Re_nh_dcp_0.Get("hErecoTot")
    h1Re_nh_dcp_0.SetTitle("MC 1Re NH #delta_{cp} = 0")
    h1Re_nh_dcp_0.SetLineWidth(2)
    h1Re_nh_dcp_0.SetLineColor(ROOT.kBlack)
    histos.append(h1Re_nh_dcp_0)

    file_1Re_ih_dcp_pi = ROOT.TFile("spectra_b_ih_dcp_pi_1Re.root", "READ")
    files.append(file_1Re_ih_dcp_pi)
    h1Re_ih_dcp_pi = file_1Re_ih_dcp_pi.Get("hErecoTot")
    h1Re_ih_dcp_pi.SetTitle("MC 1Re IH #delta_{cp} = #pi")
    h1Re_ih_dcp_pi.SetLineWidth(2)
    h1Re_ih_dcp_pi.SetLineColor(ROOT.kRed)
    histos.append(h1Re_ih_dcp_pi)

    file_1Re_ih_dcp_0 = ROOT.TFile("spectra_b_ih_dcp_0_1Re.root", "READ")
    files.append(file_1Re_ih_dcp_0)
    h1Re_ih_dcp_0 = file_1Re_ih_dcp_0.Get("hErecoTot")
    h1Re_ih_dcp_0.SetTitle("MC 1Re IH #delta_{cp} = 0")
    h1Re_ih_dcp_0.SetLineWidth(2)
    h1Re_ih_dcp_0.SetLineColor(ROOT.kBlue)
    histos.append(h1Re_ih_dcp_0)

    file_1Re_nh_dcp_pi = ROOT.TFile("spectra_b_nh_dcp_pi_1Re.root", "READ")
    files.append(file_1Re_nh_dcp_pi)
    h1Re_nh_dcp_pi = file_1Re_nh_dcp_pi.Get("hErecoTot")
    h1Re_nh_dcp_pi.SetTitle("MC 1Re NH #delta_{cp} = #pi")
    h1Re_nh_dcp_pi.SetLineWidth(2)
    h1Re_nh_dcp_pi.SetLineColor(ROOT.kGreen)
    histos.append(h1Re_nh_dcp_pi)

    file_RHC1Re_nh_dcp_0 = ROOT.TFile("spectra_b_nh_dcp_0_RHC1Re.root", "READ")
    files.append(file_RHC1Re_nh_dcp_0)
    hRHC1Re_nh_dcp_0 = file_RHC1Re_nh_dcp_0.Get("hErecoTot")
    hRHC1Re_nh_dcp_0.SetTitle("MC RHC1Re NH #delta_{cp} = 0")
    hRHC1Re_nh_dcp_0.SetLineWidth(2)
    hRHC1Re_nh_dcp_0.SetLineStyle(2)
    hRHC1Re_nh_dcp_0.SetLineColor(ROOT.kBlack)
    rhc_histos.append(hRHC1Re_nh_dcp_0)

    file_RHC1Re_ih_dcp_pi = ROOT.TFile("spectra_b_ih_dcp_pi_RHC1Re.root", "READ")
    files.append(file_RHC1Re_ih_dcp_pi)
    hRHC1Re_ih_dcp_pi = file_RHC1Re_ih_dcp_pi.Get("hErecoTot")
    hRHC1Re_ih_dcp_pi.SetTitle("MC RHC1Re IH #delta_{cp} = #pi")
    hRHC1Re_ih_dcp_pi.SetLineWidth(2)
    hRHC1Re_ih_dcp_pi.SetLineStyle(2)
    hRHC1Re_ih_dcp_pi.SetLineColor(ROOT.kRed)
    rhc_histos.append(hRHC1Re_ih_dcp_pi)

    file_RHC1Re_ih_dcp_0 = ROOT.TFile("spectra_b_ih_dcp_0_RHC1Re.root", "READ")
    files.append(file_RHC1Re_ih_dcp_0)
    hRHC1Re_ih_dcp_0 = file_RHC1Re_ih_dcp_0.Get("hErecoTot")
    hRHC1Re_ih_dcp_0.SetTitle("MC RHC1Re IH #delta_{cp} = 0")
    hRHC1Re_ih_dcp_0.SetLineWidth(2)
    hRHC1Re_ih_dcp_0.SetLineStyle(2)
    hRHC1Re_ih_dcp_0.SetLineColor(ROOT.kBlue)
    rhc_histos.append(hRHC1Re_ih_dcp_0)

    file_RHC1Re_nh_dcp_pi = ROOT.TFile("spectra_b_nh_dcp_pi_RHC1Re.root", "READ")
    files.append(file_RHC1Re_nh_dcp_pi)
    hRHC1Re_nh_dcp_pi = file_RHC1Re_nh_dcp_pi.Get("hErecoTot")
    hRHC1Re_nh_dcp_pi.SetTitle("MC RHC1Re NH #delta_{cp} = #pi")
    hRHC1Re_nh_dcp_pi.SetLineWidth(2)
    hRHC1Re_nh_dcp_pi.SetLineStyle(2)
    hRHC1Re_nh_dcp_pi.SetLineColor(ROOT.kGreen)
    rhc_histos.append(hRHC1Re_nh_dcp_pi)

    for histo in histos:
        histo.GetXaxis().SetTitle("#nu Reconstructed Energy (GeV)")
        histo.GetYaxis().SetTitle("Number of Events")
        histo.Scale(1.0 / histo.Integral())
    canvas = ROOT.TCanvas("canvas", "canvas", 10, 10, 1034, 778)
    for i, histo in enumerate(histos):
        if i > 0:
            histo.Draw("same")
            legend.AddEntry(histo, histo.GetTitle(), "l")
        else:
            histo.Draw()
            legend.AddEntry(histo, histo.GetTitle(), "l")
    legend.Draw()
    canvas.SaveAs(outfile)

    for histo in rhc_histos:
        histo.GetXaxis().SetTitle("#nu Reconstructed Energy (GeV)")
        histo.GetYaxis().SetTitle("Number of Events")
        histo.Scale(1.0 / histo.Integral())
    rhc_canvas = ROOT.TCanvas("canvas", "canvas", 10, 10, 1034, 778)
    for i, histo in enumerate(rhc_histos):
        if i > 0:
            histo.Draw("same")
            rhc_legend.AddEntry(histo, histo.GetTitle(), "l")
        else:
            histo.Draw()
            rhc_legend.AddEntry(histo, histo.GetTitle(), "l")
    rhc_legend.Draw()
    rhc_canvas.SaveAs("rhc_" + outfile)


    for file in files:
        file.Close()

def drop_bin(histo):
    # Drop the highest energy bin from muon sample for rebinning purposes,
    # we don't need it for the plot
    name = histo.GetName()
    title = histo.GetTitle()
    eAxis = histo.GetXaxis()
    nBins = eAxis.GetNbins() - 1
    eEdges = [0.0]
    for i in range(1, nBins + 1):
        eEdges.append(eAxis.GetBinUpEdge(i))
    eEdges = np.array(eEdges)
    newHisto = ROOT.TH1D(name, title, nBins, eEdges)
    for i in range(1, nBins + 1):
        newHisto.SetBinContent(i, histos[component].GetBinContent(i))
    return newHisto

def plot_ereco_spectrum(histos, data, outfile):
    canvas = ROOT.TCanvas("canvas","Spectra",10, 10, 1034, 778)
    hs = ROOT.THStack("hs", "Spectrum")
    legend = ROOT.TLegend(0.7, 0.7, 0.95, 0.95);
    for component in components:
        if "1Rmu" in outfile:
            newHisto = drop_bin(histos[component])
            if "RHC1Rmu" in outfile:
                histos[component] = newHisto.Rebin(3)
            else:
                histos[component] = newHisto.Rebin(2)
        histos[component].SetFillColor(colors[component])
        histos[component].SetMarkerStyle(21)
        histos[component].SetMarkerColor(colors[component])
        hs.Add(histos[component])
        legend.AddEntry(histos[component], labels[component], "f")
    hs.Draw("A")
    hs.GetXaxis().SetTitle("#nu Reconstructed Energy (GeV)")
    hs.GetYaxis().SetTitle("Number of Events")
    hs.GetXaxis().SetRangeUser(0,3)
    max = hs.GetMaximum()
    if not data is None:
        data1D = data.ProjectionY()
        if "1Rmu" in outfile:
            data1D = drop_bin(data1D)
            if "RHC1Rmu" in outfile:
                data1D = data1D.Rebin(3)
            else:
                data1D = data1D.Rebin(2)
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
        if np.max(y + ey) > max:
            max = np.max(y + ey)
        hs.SetMaximum(np.ceil(max))
        hs.Draw()
        graph = ROOT.TGraphErrors(len(pointsX), x, y, ex, ey)
        graph.SetMarkerStyle(21)
        graph.Draw("P")
    else:
        hs.Draw()
    legend.Draw()
    canvas.Print(outfile)

def plot_theta_spectrum(histos, data, outfile):
    canvas = ROOT.TCanvas("canvas","Spectra",10, 10, 1034, 778)
    hs = ROOT.THStack("hs", "Spectrum")
    legend = ROOT.TLegend(0.7, 0.7, 0.95, 0.95);
    for component in components:
        histos[component].SetFillColor(colors[component])
        histos[component].SetMarkerStyle(21)
        histos[component].SetMarkerColor(colors[component])
        hs.Add(histos[component])
        legend.AddEntry(histos[component], labels[component], "f")
    hs.Draw("A")
    hs.GetXaxis().SetTitle("#theta (degrees)")
    hs.GetYaxis().SetTitle("Number of Events")
#    hs.GetXaxis().SetRangeUser(0,3)
    max = hs.GetMaximum()
    if not data is None:
        data1D = data.ProjectionX()
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
        if np.max(y + ey) > max:
            max = np.max(y + ey)
        hs.SetMaximum(np.ceil(max))
        hs.Draw()
        graph = ROOT.TGraphErrors(len(pointsX), x, y, ex, ey)
        graph.SetMarkerStyle(21)
        graph.Draw("P")
    else:
        hs.Draw()
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
        data.Draw("P0")
    canvas.Print(outfile)

def plot_data_ereco_theta_spectrum(histo, outfile):
    canvas = ROOT.TCanvas("canvas","Spectra", 0, 0, 1000, 750)
    total = histo.Clone()

    thetaAxis = total.GetXaxis()
    nThetaBins = thetaAxis.GetNbins()
    thetaEdSKdata_Run1to8wGraph.rootges = []
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
        erecHistos = dict()
        thetaHistos = dict()
        for component in components:
            histos[component] = file.Get("hErecoVsQ2Mode_" + component)
            erecHistos[component] = file.Get("hErecoMode_" + component)
            thetaHistos[component] = file.Get("hQ2Mode_" + component)

        plot_ereco_theta_spectrum(histos, graph, outfile)
        plot_ereco_spectrum(erecHistos, data, "erec_" + outfile)
        plot_theta_spectrum(thetaHistos, data, "theta_" + outfile)

def data(infile, outfile, sample):
    file = ROOT.TFile(infile)
    if sample in ["1Rmu", "RHC1Rmu"]:
        histo = file.Get("mock_data_" + sample + "_0")
        plot_data_ereco_spectrum(histo, outfile)
    else:
        histo = file.Get("mock_data_" + sample + "_0")
        plot_data_ereco_theta_spectrum(histo, outfile)

def dcp_0(outfile):
    fhc_1re_file = ROOT.TFile("spectra_osc_ichikawa_1Re.root")
    rhc_1re_file = ROOT.TFile("spectra_osc_ichikawa_RHC1Re.root")
    fhc_cc1pi_file = ROOT.TFile("spectra_osc_ichikawa_MultiRe.root")
    data_file = ROOT.TFile("SKdata_Run1to8wGraph.root")

    fhc_mc = combine_elike(fhc_1re_file.Get("hErecoTot"), fhc_cc1pi_file.Get("hErecoTot"))
    rhc_mc = rhc_1re_file.Get("hErecoTot")
    fhc_mc = fhc_mc.Rebin(2)
    rhc_mc = rhc_mc.Rebin(2)
    fhc_data = combine_elike(data_file.Get("mock_data_1Re_0"), data_file.Get("mock_data_MultiRe_0"))
    rhc_data = data_file.Get("mock_data_RHC1Re_0")
    rhc_mc = negate_histo(rhc_mc)
    fhc_data = process_data(fhc_data)
    rhc_data = process_data(rhc_data, True)
    plot_bestfit_dcp0_spectrum(fhc_mc, rhc_mc, fhc_data, rhc_data, outfile)

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
  
    #asimov(o.infile, o.outfile, o.datafile, o.sample)
    #data(o.infile, o.outfile, o.sample)

    #make_asimovb_variation_spectrum(o.outfile)
  
    dcp_0("ichikawa2.pdf")

    if other:
        print "Warning - ignored arguments: ", other
