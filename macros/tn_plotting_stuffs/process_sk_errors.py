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

def CreateEnvelope(N, filename):
    binEvents = np.zeros(N)
    file = ROOT.TFile(filename)
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
    c.SaveAs(filename + ".pdf")

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
    CreateEnvelope(10000, sys.argv[1])
