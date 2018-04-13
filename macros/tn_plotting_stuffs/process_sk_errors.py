#!/usr/bin/env python

import ROOT
import numpy as np
import sys
  
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

samples = ["1Rmu", "RHC1Rmu", "1Re", "RHC1Re", "MultiRe"]
groups = ["skdet", "skfsi", "banff_constrained", "xsec_nue", "xsec_ncgamma", "xsec_ncoth", "xsec_eb", "osc", "allsyst", "all"]
banff_fits = ["prefit", "postfit"]
names = [
    "SK Detector",
    "SK FSI+SI+PN",
    "Flux + Xsec constrained",
    r"$\sigma(\nu_e)/\sigma(\bar{\nu}_e)$",
    r"NC1$\gamma$",
    "NC Other",
    r"E$_\text{b}$",
    "Osc",
    "All Systematics",
    "All with osc"
    ]
modes = [
   "CCQE",
   "CC2p2h",
   "CCcoh",
   "CC1pi",
   "CCoth",
   "NC1pi0",
   "NC1piC",
   "NCcoh",
   "NC1gamma",
   "NCoth"
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
#    caveat.Draw("same");

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
    title = ROOT.TLatex(tx, ty,"T2K Run1-9c Preliminary");
    title.SetTextSize(0.05);
    title.Draw("same");

    tx = upper.GetXaxis().GetBinLowEdge(upper.GetXaxis().GetLast()) * 0.1 
    ty = upper.GetMaximum() * 1.01
    caveat = ROOT.TLatex(tx, ty,"Final systematics pending");
    caveat.SetTextSize(0.025);
#    caveat.Draw("same");

    c.SaveAs(outfile + "_fractional.pdf")
    c.SaveAs(outfile + "_fractional.png")
    c.SaveAs(outfile + "_fractional.eps")

def CreateSummaryTable(N):
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
                for i in range(N):
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

def CreateByModeErecPlots(N):
    banff = "postfit"
    group = "allsyst"

    tot_events_by_erec_mean = dict()
    for s, sample in enumerate(samples):
        histofilename = "sk_error_table_%s_%s_%s_histos.root" % (sample, group, banff)
        f = ROOT.TFile(histofilename)
        if not f:
            continue
        dummy = f.Get("events_0")
        Nbins = dummy.GetNbinsX()
        tot_events_by_erec = np.zeros((N,Nbins))
        for i in range(N):
            histo = f.Get("events_%d" % i)
            for j in range(1, Nbins + 1):
                tot_events_by_erec[i, j - 1] = histo.GetBinContent(j)
        tot_events_by_erec_mean[s] = np.mean(tot_events_by_erec, axis=0)
        print "%s sum bin means %f" % (samples[s], np.sum(tot_events_by_erec_mean[s]))
        f.Close()

    sample_erec = dict()
    for s, sample in enumerate(samples):
        histofilename = "sk_error_table_%s_%s_%s_histos.root" % (sample, group, banff)
        f = ROOT.TFile(histofilename)
        if not f:
            continue
        dummy = f.Get("events_0")
        Nbins = dummy.GetNbinsX()
        for p in range(len(modes)):
            mode_events_by_erec = np.zeros((N,Nbins))
            for i in range(N):
                histo = f.Get("events_%d_%d" % (i,p))
                for j in range(1, Nbins + 1):
                    mode_events_by_erec[i, j - 1] = histo.GetBinContent(j)
            mode_events_by_erec_rms = np.zeros(Nbins)
            for i in range(Nbins):
                events = np.array(mode_events_by_erec[:,i])
                #mode_events_by_erec_rms[i] = ROOT.TMath.Mean(len(events), events)
                mode_events_by_erec_rms[i] = np.std(events)
            err_histo = ROOT.TH1D("err_histo_" + str(s) + "_" + str(p), samples[s] + " " + modes[p], Nbins, histo.GetXaxis().GetXmin(), histo.GetXaxis().GetXmax())
            if p == len(modes) - 1:
                err_histo.SetLineColor(p + 2)
            else:
                err_histo.SetLineColor(p + 1)
            for j in range(Nbins):
                if tot_events_by_erec_mean[s][j] != 0:
                    err_histo.SetBinContent(j + 1, mode_events_by_erec_rms[j] / tot_events_by_erec_mean[s][j])
                else:
                    err_histo.SetBinContent(j + 1, 0)
            err_histo.SetDirectory(0)
            sample_erec[str(s) + "_" + str(p)] = err_histo
        f.Close()

    canvas = ROOT.TCanvas("err_canvas", "err_canvas", 1024, 768)
    canvas.SetTopMargin(0.1)
    canvas.SetLeftMargin(0.15)
    canvas.SetBottomMargin(0.15)
    canvas.SetRightMargin(0.2)

    sample_max = np.zeros(len(samples))

    for s, sample in enumerate(samples):
        for p in range(len(modes)):
            histo = sample_erec[str(s) + "_" + str(p)]
            if histo.GetMaximum() > sample_max[s]:
                sample_max[s] = histo.GetMaximum()

    for s, sample in enumerate(samples):
        opt = ""
        dummy = sample_erec["0_0"]
        thesum = np.zeros(dummy.GetNbinsX())
        legend = ROOT.TLegend(0.85, 0.5, 0.95, 0.8)
        legend.SetFillColor(0);
        legend.SetLineColor(0);
        legend.SetTextFont(132);
        legend.SetTextSize(0.03);

        for p in range(len(modes)):
            histo = sample_erec[str(s) + "_" + str(p)]
            histo.SetTitle(samples[s])
            histo.GetYaxis().SetRangeUser(0, sample_max[s])
            histo.Draw(opt)
            legend.AddEntry(histo, modes[p], "l")
            opt = "same"
            for i in range(1, histo.GetNbinsX()):
                thesum[i - 1] += histo.GetBinContent(i)
        legend.Draw()
        if s == 0:
            canvas.Print("mode_errors.pdf(")
        elif s == (len(samples) - 1):
            canvas.Print("mode_errors.pdf)")
        else:
            canvas.Print("mode_errors.pdf")

def CreateByModeSummaryTable(N):
    banff = "postfit"
    group = "allsyst"
    line = r"\begin{tabular}{|l|l|l|l|l|l|}" + "\n"
    line += r"\hline" + "\n"
    line += r"& \multicolumn{2}{c||}{\textbf{1-Ring $\mu$}} & \multicolumn{3}{c|}{\textbf{1-Ring $e$}} \\ \hline" + "\n"
    line += r"Interaction mode & FHC & RHC & FHC & RHC & FHC 1 d.e. \\ \hline\hline" + "\n"

    for p in range(len(modes)):
        line += modes[p]
        for s, sample in enumerate(samples):
            events = []
            histofilename = "sk_error_table_%s_%s_%s_histos.root" % (sample, group, banff)
            f = ROOT.TFile(histofilename)
            if not f:
                continue
            for i in range(N):
                histo = f.Get("events_%d_%d" % (i,p))
                events.append(histo.GetEntries())
            events = np.array(events)
            mean = ROOT.TMath.Mean(len(events), events)
#            rms = ROOT.TMath.RMS(len(events), events)
            rms = np.std(events)

            totevents = []
            for i in range(N):
                tothisto = f.Get("events_%d" % i)
                totevents.append(tothisto.GetEntries())
            totevents = np.array(totevents)
            totmean = ROOT.TMath.Mean(len(totevents), totevents)

            if mean != 0:
                line += " & %.4f" % (rms / totmean)
            else:
                line += " & -"
            f.Close()
        line += r" \\ \hline" + "\n"
    line += r"\end{tabular}"
    filename = "sk_by_mode_error_table.tex"
    target = open(filename, 'w')
    target.truncate()
    target.write(line)
    target.close()

if __name__=="__main__":
    N = 10000
    if len(sys.argv) == 1:
#        CreateSummaryTable(N)
#        CreateByModeSummaryTable(N)
        CreateByModeErecPlots(N)
    elif len(sys.argv) == 3:
        CreateEnvelope(N, sys.argv[1], sys.argv[2])
        CreateFractionalError(N, sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        CreateEnvelope(N, sys.argv[1], sys.argv[2], sys.argv[3])
        CreateFractionalError(N, sys.argv[1], sys.argv[2], sys.argv[3])
        
