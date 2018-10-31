#!/usr/bin/env python

import ROOT
import argparse
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

def plot(chisqs, chisq_bf, out_filename):
    canvas = ROOT.TCanvas("fc_nh", "fc_nh", 1024, 768)
    canvas.SetLeftMargin(0.1)
    canvas.SetBottomMargin(0.1)
    canvas.SetRightMargin(0.1)
    min_chisq = int(np.floor(min(chisqs)))
    max_chisq = int(np.ceil(max(chisqs)))
    num_bins = max_chisq - min_chisq

    h_chisq_sample_min = ROOT.TH1D("toy_chisq", "toy_chisq", num_bins, min_chisq, max_chisq)
    h_chisq_sample_min.SetDirectory(0)
    for chisq in chisqs:
        h_chisq_sample_min.Fill(chisq)
    h_chisq_sample_min.GetXaxis().SetTitle("#chi^{2}")
    h_chisq_sample_min.GetYaxis().SetTitle("Toys")
    h_chisq_sample_min.SetLineColor(ROOT.kBlack)
    h_chisq_sample_min.SetLineWidth(1)
    h_chisq_sample_min.Draw()

    chisqs.sort()
    pvalue = 1
    num_toys = float(len(chisqs))
    for i, chisq in enumerate(chisqs):
#        if i < 2000:
#            print i, chisq, chisq_bf
        if chisq > chisq_bf:
            pvalue = i / num_toys
            break
    # One-sided p-value
    pvalue = 1 - pvalue
    print "p-value: ", pvalue
    ptext = ROOT.TLatex(chisq_bf + 1, h_chisq_sample_min.GetMaximum(), "p = %0.3f" % pvalue)
    ptext.SetTextColor(ROOT.kBlue)
    ptext.SetTextSize(0.05)
    ptext.Draw()

    line = ROOT.TLine(chisq_bf, 0, chisq_bf, h_chisq_sample_min.GetMaximum())
    line.SetLineColor(ROOT.kBlue)
    line.SetLineWidth(2)
    line.Draw("SAME")

    canvas.Print(out_filename)

def calculate_pvalue(toys_filename, data_filename, sample):
    f_toys = ROOT.TFile(toys_filename, "READ")
    toys_tree = f_toys.Get("fitnt")
    chisq_min_list = []
    for event in toys_tree:
        h_chisq_5sample_nh = event.hchisq1D_nh
        h_chisq_5sample_ih = event.hchisq1D_ih
        chisq_5sample_min = min(h_chisq_5sample_nh.GetMinimum(), h_chisq_5sample_ih.GetMinimum())
        if sample.lower() == "all":
            chisq_min_list.append(chisq_5sample_min)
        else:
            if sample.lower() == "1rmu":
                h_chisq_sample_nh = event.hchisq1D_1Rmu_nh
                h_chisq_sample_ih = event.hchisq1D_1Rmu_ih
            elif sample.lower() == "rhc1rmu":
                h_chisq_sample_nh = event.hchisq1D_RHC1Rmu_nh
                h_chisq_sample_ih = event.hchisq1D_RHC1Rmu_ih
            elif sample.lower() == "1re":
                h_chisq_sample_nh = event.hchisq1D_1Re_nh
                h_chisq_sample_ih = event.hchisq1D_1Re_ih
            elif sample.lower() == "rhc1re":
                h_chisq_sample_nh = event.hchisq1D_RHC1Re_nh
                h_chisq_sample_ih = event.hchisq1D_RHC1Re_ih
            else:
                h_chisq_sample_nh = event.hchisq1D_MultiRe_nh
                h_chisq_sample_ih = event.hchisq1D_MultiRe_ih
            if h_chisq_5sample_nh.GetMinimum() <= h_chisq_5sample_ih.GetMinimum():
                minbin = h_chisq_5sample_nh.GetMinimumBin()
                chisq = h_chisq_sample_nh.GetBinContent(minbin)
                if chisq > 0:
                    chisq_min_list.append(chisq)
            else:
                minbin = h_chisq_5sample_ih.GetMinimumBin()
                chisq = h_chisq_sample_ih.GetBinContent(minbin)
                if chisq > 0:
                    chisq_min_list.append(chisq)
    f_toys.Close()

    f_data = ROOT.TFile(data_filename, "READ")
    data_tree = f_data.Get("fitnt")
    pvalue = 0
    for event in data_tree:
        h_chisq_5sample_fit = event.hchisq1D
        if sample.lower() == "all":
            chisq_bf = h_chisq_5sample_fit.GetMinimum()
        else:
            if sample.lower() == "1rmu":
                chisq_bf = event.hchisq1D_1Rmu.GetBinContent(h_chisq_5sample_fit.GetMinimumBin())
            elif sample.lower() == "rhc1rmu":
                chisq_bf = event.hchisq1D_RHC1Rmu.GetBinContent(h_chisq_5sample_fit.GetMinimumBin())
            elif sample.lower() == "1re":
                chisq_bf = event.hchisq1D_1Re.GetBinContent(h_chisq_5sample_fit.GetMinimumBin())
            elif sample.lower() == "rhc1re":
                chisq_bf = event.hchisq1D_RHC1Re.GetBinContent(h_chisq_5sample_fit.GetMinimumBin())
            else:
                chisq_bf = event.hchisq1D_MultiRe.GetBinContent(h_chisq_5sample_fit.GetMinimumBin())
    f_data.Close()

    return (chisq_min_list, chisq_bf)

def check_cc1pi(toys_filename):
    f_toys = ROOT.TFile(toys_filename, "READ")
    toys_tree = f_toys.Get("fitnt")
    chisq_min_list = []
    for t, event in enumerate(toys_tree):
        h_chisq_sample_nh = event.hchisq1D_MultiRe_nh
        h_chisq_sample_ih = event.hchisq1D_MultiRe_ih
        for i in range(1, h_chisq_sample_nh.GetNbinsX() + 1):
            if h_chisq_sample_nh.GetBinContent(i) == 0:
                print "Toy", t, "NH", i, "found 0"
            if h_chisq_sample_ih.GetBinContent(i) == 0:
                print "Toy", t, "IH", i, "found 0"

    f_toys.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate per sample p-values")
    parser.add_argument("--toys", type=str, help="FC toys file", required=True)
    parser.add_argument("--data", type=str, help="Data fit file", required=True)
    parser.add_argument("--sample", type=str, help="SK Sample", required=True)
    parser.add_argument("--output", type=str, help="Name of the output PDF file", required=True)
    args = parser.parse_args()
#    chisq, pvalue = calculate_pvalue(args.toys, args.data, args.sample)
#    plot(chisq, pvalue, args.output)
    check_cc1pi(args.toys)
