#!/usr/bin/env python

import ROOT
import numpy as np
import sys
  
ROOT.gROOT.SetBatch(1)

samples = ["1Rmu", "RHC1Rmu", "1Re", "RHC1Re", "MultiRe"]
groups = ["skdet", "skfsi", "banff_constrained", "xsec_nue", "xsec_ncgamma", "xsec_ncoth", "osc", "allsyst", "all"]
banff_fits = ["prefit", "postfit"]
names = [
    "SK Detector",
    "SK FSI+SI+PN",
    "Flux + Xsec constrained",
    r"$\sigma(\nu_e)/\sigma(\bar{\nu}_e)$",
    r"NC1$\gamma$",
    "NC Other",
    "Osc",
    "All Systematics",
    "All with osc"
    ]

def CreateEnvelope(N, outfile, filename1, filename2=None):
    ROOT.gStyle.SetOptTitle(0)
    binEvents = np.zeros(N)
    file = ROOT.TFile(filename1)
    histos = []
    for i in range(0, N):
        histos.append(file.Get("events_%d" % i))
    nominal = histos[0].Clone()
    lower = histos[0].Clone()
    upper = histos[0].Clone()
    nbins = nominal.GetNbinsX()
    for i in range(1, nbins + 1):
        for j in range(0, N):
            binEvents[j] = histos[j].GetBinContent(i)
        nom = ROOT.TMath.Mean(N, binEvents)
        err = ROOT.TMath.RMS(N, binEvents)
        nominal.SetBinContent(i, nom)
        upper.SetBinContent(i, nom + err)
        lower.SetBinContent(i, nom - err)
    c = ROOT.TCanvas("c1","c1", 10,10, 1034,778)
    nominal.SetMinimum(0.0)
    nominal.SetMaximum(upper.GetMaximum() + 0.5)
    max = nominal.GetXaxis().GetBinCenter(nbins)
    if(max > 2.9):
        nominal.GetXaxis().SetRangeUser(0.0, 2.9)
    nominal.Draw("")
    nominal.GetXaxis().SetTitleOffset(0.8)
    nominal.GetYaxis().SetTitleOffset(0.8)
    upper.SetFillColor(ROOT.kRed)
    upper.SetFillStyle(1001)
    upper.SetLineWidth(1)
    upper.Draw("same")
    lower.SetFillColor(10)
    lower.SetFillStyle(1001)
    lower.SetLineWidth(1)
    lower.Draw("same")
    nominal.SetLineColor(ROOT.kBlack)
    nominal.SetLineStyle(1)
    nominal.SetLineWidth(2)
    nominal.Draw("same")

    if filename2 is not None:
        binEvents2 = np.zeros(N)
        file2 = ROOT.TFile(filename2)
        histos2 = []
        for i in range(0, N):
            histos2.append(file2.Get("events_%d" % i))
        nominal2 = histos2[0].Clone()
        lower2 = histos2[0].Clone()
        upper2 = histos2[0].Clone()
        nbins2 = nominal2.GetNbinsX()
        for i in range(1, nbins2 + 1):
            for j in range(0, N):
                binEvents2[j] = histos2[j].GetBinContent(i)
            nom = ROOT.TMath.Mean(N, binEvents2)
            err = ROOT.TMath.RMS(N, binEvents2)
            nominal2.SetBinContent(i, nom)
            upper2.SetBinContent(i, nom + err)
            lower2.SetBinContent(i, nom - err)
        nominal2.SetMinimum(0.0)
        nominal2.SetMaximum(upper2.GetMaximum() + 0.5)
        if(lower2.GetMinimum() < nominal.GetMinimum()):
            nominal.SetMinimum(lower2.GetMinimum() * 1.1)
        if(upper2.GetMaximum() > nominal.GetMaximum()):
            nominal.SetMaximum(upper2.GetMaximum() * 1.1)
        upper2.SetLineColor(ROOT.kBlue)
        upper2.SetLineStyle(2)
        upper2.SetLineWidth(2)
        upper2.Draw("same")
        lower2.SetLineColor(ROOT.kBlue)
        lower2.SetLineStyle(2)
        lower2.SetLineWidth(2)
        lower2.Draw("same")
        nominal2.SetLineColor(ROOT.kBlue)
        nominal2.SetLineStyle(1)
        nominal2.SetLineWidth(1)
        nominal2.Draw("same")

    tx = nominal.GetXaxis().GetBinLowEdge(nominal.GetXaxis().GetLast()) * 0.5
    ty = nominal.GetMaximum() * 1.01
    title = ROOT.TLatex(tx, ty,"T2K Run1-8 Preliminary");
    title.SetTextSize(0.05);
    title.Draw("same");

    tx = nominal.GetXaxis().GetBinLowEdge(nominal.GetXaxis().GetLast()) * 0.1 
    ty = nominal.GetMaximum() * 1.01
    caveat = ROOT.TLatex(tx, ty,"Final systematics pending");
    caveat.SetTextSize(0.025);
    caveat.Draw("same");

    c.SaveAs(outfile + ".pdf")
    c.SaveAs(outfile + ".png")
    c.SaveAs(outfile + ".eps")

def CreateFractionalError(N, outfile, filename1, filename2=None):
    ROOT.TGaxis.SetMaxDigits(2)
    ROOT.gStyle.SetOptTitle(0)
    binEvents = np.zeros(N)
    file = ROOT.TFile(filename1)
    histos = []
    for i in range(0, N):
        histos.append(file.Get("events_%d" % i))
    lower = histos[0].Clone()
    upper = histos[0].Clone()
    nbins = histos[0].GetNbinsX()

    max = histos[0].GetXaxis().GetBinCenter(nbins)
    xaxisLimit = histos[0].GetXaxis().GetBinUpEdge(nbins)
    if(max > 2.9):
#        upper.GetXaxis().SetRange(0, 70)
#        lower.GetXaxis().SetRange(0, 70)
        upper.GetXaxis().SetRangeUser(0.0, 2.9)
        lower.GetXaxis().SetRangeUser(0.0, 2.9)
        xaxisLimit = 2.9
    for i in range(1, nbins + 1):
        for j in range(0, N):
            binEvents[j] = histos[j].GetBinContent(i)
        nom = ROOT.TMath.Mean(N, binEvents)
        err = ROOT.TMath.RMS(N, binEvents)
        if nom == 0:
            frac = 0
        else:
            frac = err / nom
        upper.SetBinContent(i, frac)
        lower.SetBinContent(i, -frac)
    c = ROOT.TCanvas("c1","c1", 10,10, 1034,778)

    upper.SetMinimum(lower.GetMinimum() * 1.1)
    upper.SetMaximum(upper.GetMaximum() * 1.1)
    upper.Draw("")
    upper.GetYaxis().SetTitle("Fractional error")
    line = ROOT.TLine(0, 0, xaxisLimit, 0)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(2)
    line.Draw("same")
    upper.GetXaxis().SetTitleOffset(0.8)
    upper.GetYaxis().SetTitleOffset(0.8)
    upper.SetLineColor(ROOT.kRed)
    upper.SetLineStyle(1)
    upper.SetLineWidth(2)
    upper.Draw("same")
    lower.SetLineColor(ROOT.kRed)
    lower.SetLineStyle(1)
    lower.SetLineWidth(2)
    lower.Draw("same")

    if filename2 is not None:
        binEvents2 = np.zeros(N)
        file2 = ROOT.TFile(filename2)
        histos2 = []
        for i in range(0, N):
            histos2.append(file2.Get("events_%d" % i))
        lower2 = histos2[0].Clone()
        upper2 = histos2[0].Clone()
        nbins2 = histos2[0].GetNbinsX()
        for i in range(1, nbins2 + 1):
            for j in range(0, N):
                binEvents2[j] = histos2[j].GetBinContent(i)
            nom = ROOT.TMath.Mean(N, binEvents2)
            err = ROOT.TMath.RMS(N, binEvents2)
            if nom == 0:
                frac = 0
            else:
                frac = err / nom
            upper2.SetBinContent(i, frac)
            lower2.SetBinContent(i, -frac)
        if(lower2.GetMinimum() < upper.GetMinimum()):
            upper.SetMinimum(lower2.GetMinimum() * 1.1)
        if(upper2.GetMaximum() > upper.GetMaximum()):
            upper.SetMaximum(upper2.GetMaximum() * 1.1)
        upper2.SetLineColor(ROOT.kBlue)
        upper2.SetLineStyle(2)
        upper2.SetLineWidth(2)
        upper2.Draw("same")
        lower2.SetLineColor(ROOT.kBlue)
        lower2.SetLineStyle(2)
        lower2.SetLineWidth(2)
        lower2.Draw("same")

    tx = upper.GetXaxis().GetBinLowEdge(upper.GetXaxis().GetLast()) * 0.5
    ty = upper.GetMaximum() * 1.01
    title = ROOT.TLatex(tx, ty,"T2K Run1-8 Preliminary");
    title.SetTextSize(0.05);
    title.Draw("same");

    tx = upper.GetXaxis().GetBinLowEdge(upper.GetXaxis().GetLast()) * 0.1 
    ty = upper.GetMaximum() * 1.01
    caveat = ROOT.TLatex(tx, ty,"Final systematics pending");
    caveat.SetTextSize(0.025);
    caveat.Draw("same");

    c.SaveAs(outfile + "_fractional.pdf")
    c.SaveAs(outfile + "_fractional.png")
    c.SaveAs(outfile + "_fractional.eps")

def CreateSummaryTable():
    for b, banff in enumerate(banff_fits):
        for g, group in enumerate(groups):
            line = names[g]
            # Skip prefit for now, histos aren't created. In the future skip conditionally
            if banff == "prefit":
                continue
            if banff == "prefit" and not group in ["banff_constrained", "allsyst", "all"]:
                continue
            fhcrhc = []
            for s, sample in enumerate(samples):
                events = []
                histos = []
                histofilename = "sk_error_table_%s_%s_%s_histos.root" % (sample, group, banff)
                f = ROOT.TFile(histofilename)
                if not f:
                    continue
                for i in range(10000):
                    histo = f.Get("events_%d" % i)
                    histos.append(histo)
                    events.append(histo.GetEntries())
                    if sample == "1Re":
                        fhcrhc.append(events[i])
                    elif sample == "RHC1Re":
                        fhcrhc[i] /= events[i]
                events = np.array(events)
                mean = ROOT.TMath.Mean(len(events), events)
                rms = ROOT.TMath.RMS(len(events), events)
                line += " & %.2f" % (100 * rms / mean)
                if sample == "MultiRe":
                    fhcrhc = np.array(fhcrhc)
                    mean = ROOT.TMath.Mean(len(fhcrhc), fhcrhc)
                    rms = ROOT.TMath.RMS(len(fhcrhc), fhcrhc)
                    line += " & %.2f" % (100 * rms / mean)
            line += r" \\ \hline"
            filename = "sk_error_table_%s.tex" % group
            target = open(filename, 'w')
            target.truncate()
            target.write(line)
            target.close()

if __name__=="__main__":
#    CreateSummaryTable()
    if len(sys.argv) == 3:
        CreateEnvelope(10000, sys.argv[1], sys.argv[2])
        CreateFractionalError(10000, sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        CreateEnvelope(10000, sys.argv[1], sys.argv[2], sys.argv[3])
        CreateFractionalError(10000, sys.argv[1], sys.argv[2], sys.argv[3])
        
