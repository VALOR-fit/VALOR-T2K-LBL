#!/usr/bin/env python

import ROOT
import argparse
import sys
import numpy as np

from os.path import exists
from array import array

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

hierarchy = dict(zip(["nh", "ih"], ["normal", "inverted"]))
points = dict(\
    zip(\
    ["mpi", "m3pi4", "mpi2", "mpi4", "0", "ppi4", "ppi2", "p3pi4"],\
    [-3.14, -2.36, -1.57, -0.79, 0.0, 0.79, 1.57, 2.36])\
    )

def get_tree(f, treename):
    keys = [key.GetName() for key in f.GetListOfKeys()]
    if not treename in keys:
        print "Error: Could not find key \'%s\'" % treename
        sys.exit()
    tree = f.Get(treename)
    branches = [branch.GetName() for branch in tree.GetListOfBranches()]
    req_branches = ["chisqmin_nh", "chisqmin_ih", "trueMH"]
    for branch in req_branches:
        if not branch in branches:
            print "Error: Could not find branch \'%s\'" % branch
            sys.exit()
    return tree

def get_chisq_from_histos(hchisq1D_nh, hchisq1D_ih, method):
    if method.lower() == "average":
        N = hchisq1D_nh.GetNbinsX()
        nh_vals = np.zeros(N)
        ih_vals = np.zeros(N)
        sum_likelihood_nh = 0.0
        sum_likelihood_ih = 0.0
        for i in range(N):
            sum_likelihood_nh += np.exp(-hchisq1D_nh.GetBinContent(i + 1) / 2)
            sum_likelihood_ih += np.exp(-hchisq1D_ih.GetBinContent(i + 1) / 2)
        chisq_nh = -2 * np.log(sum_likelihood_nh / N)
        chisq_ih = -2 * np.log(sum_likelihood_ih / N)
        return chisq_nh, chisq_ih
    elif method.lower() == "bestfit":
        chisq_nh = hchisq1D_nh.GetBinContent(hchisq1D_nh.GetMinimumBin())
        chisq_ih = hchisq1D_ih.GetBinContent(hchisq1D_ih.GetMinimumBin())
        return chisq_nh, chisq_ih

def get_toy_mh_chisq(fc_dir, true_dcp, true_mh):
    filename = fc_dir + "/" + true_mh + "/fc_deltacp_" + true_dcp + "_" + hierarchy[true_mh] + "_merge.root"
    if not exists(filename):
        print "Error: No such file \'%s\'" % filename
        sys.exit()
    f = ROOT.TFile(filename, "READ")
    tree = get_tree(f, "fitnt")
    print "%s: %d" % (filename, tree.GetEntries())
    chisq_nh = np.zeros(tree.GetEntries())
    chisq_ih = np.zeros(tree.GetEntries())
    for t, toy in enumerate(tree):
        if not toy.trueMH == hierarchy[true_mh]:
            print "Error: Expected true hierarchy %s in %s, found %s" % (true_mh, filename, true_mh[0])
        chisq_nh[t], chisq_ih[t] = get_chisq_from_histos(toy.hchisq1D_nh, toy.hchisq1D_ih, "bestfit")
    f.Close()
    return chisq_nh, chisq_ih

def get_data_mh_chisq(filename_nh, filename_ih):
    likelihoods = np.zeros(2)
    f = ROOT.TFile(filename_nh, "READ")
    hchisq_nh = f.Get("h_chisq")
    hchisq_nh.SetDirectory(0)
    f.Close()
    f = ROOT.TFile(filename_ih, "READ")
    hchisq_ih = f.Get("h_chisq")
    hchisq_ih.SetDirectory(0)
    f.Close()
    chisq_nh, chisq_ih = get_chisq_from_histos(hchisq_nh, hchisq_ih, "bestfit")
    return chisq_nh, chisq_ih

def get_dchisq_mh(chisq_nh, chisq_ih):
    return chisq_nh - chisq_ih

def get_histograms_mh(dchisq_toys_tnh, dchisq_toys_tih, dchisq_data):
    min_dchisq = np.floor(min(np.min(dchisq_toys_tnh), np.min(dchisq_toys_tih), np.min(dchisq_data)))
    max_dchisq = np.floor(max(np.max(dchisq_toys_tnh), np.max(dchisq_toys_tih), np.max(dchisq_data)))
    h_toys_tnh = ROOT.TH1D("toys_tnh", "toys_tnh", 10 * int(max_dchisq - min_dchisq), min_dchisq, max_dchisq)
    h_toys_tnh.SetDirectory(0)
    h_toys_tnh.GetXaxis().SetTitle("#Delta#chi^{2}")
    h_toys_tnh.GetYaxis().SetTitle("Number of toys")
    h_toys_tnh.SetLineColor(ROOT.kBlue)
    for dchisq in dchisq_toys_tnh:
        h_toys_tnh.Fill(dchisq)
    h_toys_tih = h_toys_tnh.Clone("toys_tih")
    h_toys_tih.SetTitle("toys_tih")
    h_toys_tih.Reset()
    h_toys_tih.SetDirectory(0)
    h_toys_tih.SetLineColor(ROOT.kRed)
    for dchisq in dchisq_toys_tih:
        h_toys_tih.Fill(dchisq)
    maxval = max(h_toys_tnh.GetMaximum(), h_toys_tih.GetMaximum())
    h_toys_tnh.GetYaxis().SetRangeUser(0, maxval)
    h_toys_tih.GetYaxis().SetRangeUser(0, maxval)
    return h_toys_tnh, h_toys_tih  

def get_pvalue_two_sided(toys, dchisq):
    Ntoys = len(toys)
    mid_point = 0
    if (Ntoys % 2 == 0):
        mid_point = Ntoys // 2
        mid_dchisq = (toys[mid_point] + toys[mid_point - 1]) / 2.0
    else:
        mid_point = Ntoys // 2
        mid_dchisq = toys[mid_point]
    pval = 0.0
    if dchisq > mid_dchisq:
        for i in range(Ntoys - 1, mid_point, -1):
            if dchisq > toys[i]:
                pval = 2 * pval / Ntoys
                break
            pval += 1.0
    else:
        for i in range(0, mid_point):
            if dchisq < toys[i]:
                pval = 2 * float(i) / Ntoys
                break
            pval += 1.0
    return pval

def plot_dchisq_single(h_toys_tnh, h_toys_tih, dchisq_data, output_pdf):
    canvas = ROOT.TCanvas("canvas", "canvas", 1024, 768)
    h_toys_tnh.Draw("")
    h_toys_tih.Draw("same")
    maxval = h_toys_tnh.GetMaximum()
    data_line = ROOT.TLine(dchisq_data, 0.0, dchisq_data, maxval)
    data_line.SetLineColor(ROOT.kBlack)
    data_line.SetLineWidth(2)
    data_line.Draw("SAME")
    ptext = ROOT.TLatex(dchisq_data + 0.5, maxval * 0.75, "Data")
    ptext.SetTextColor(ROOT.kBlack)
    ptext.SetTextSize(0.04)
    ptext.Draw()

    canvas.SaveAs(output_pdf)

def plot_dchisq_all(dchisq_toys_tnh_dict, dchisq_toys_tih_dict, dchisq_data, pvalue):
    ROOT.gStyle.SetEndErrorSize(6)
    canvas = ROOT.TCanvas("canvas", "canvas", 1024, 768)
    canvas.SetTopMargin(0.1)
    canvas.SetLeftMargin(0.1)
    canvas.SetBottomMargin(0.1)
    canvas.SetRightMargin(0.15)

    N = len(dchisq_toys_tnh_dict)
    dcp = np.zeros(N)
    data_nh = np.zeros(N)
    data_ih = np.zeros(N)
    toys_95_low_nh = np.zeros(N)
    toys_95_high_nh = np.zeros(N)
    toys_95_low_ih = np.zeros(N)
    toys_95_high_ih = np.zeros(N)
    toys_100_low_nh = np.zeros(N)
    toys_100_high_nh = np.zeros(N)
    toys_100_low_ih = np.zeros(N)
    toys_100_high_ih = np.zeros(N)

    data_95_nh = np.zeros(N)
    data_95_ih = np.zeros(N)
    data_100_nh = np.zeros(N)
    data_100_ih = np.zeros(N)

    toys_95_nh = np.zeros(N)
    toys_95_ih = np.zeros(N)
    toys_100_nh = np.zeros(N)
    toys_100_ih = np.zeros(N)

    keys = list(points.keys())
    for i in range(0, N):
        dcp[i] = points[keys[i]]
        dchisq_toys_tnh_sorted = np.sort(dchisq_toys_tnh_dict[keys[i]])
        Ntnh = len(dchisq_toys_tnh_sorted)
        dchisq_toys_tih_sorted = np.sort(dchisq_toys_tih_dict[keys[i]])
        Ntih = len(dchisq_toys_tih_sorted)

        data_nh[i] = data_ih[i] = dchisq_data
        plow = (1.0 - pvalue) / 2
        phigh = 1.0 - plow
        toys_95_low_nh[i] = dchisq_toys_tnh_sorted[int(np.floor(plow * Ntnh) - 1)]
        toys_95_high_nh[i] = dchisq_toys_tnh_sorted[int(np.floor(phigh * Ntnh) - 1)]
        toys_95_low_ih[i] = dchisq_toys_tih_sorted[int(np.floor(plow * Ntih) - 1)]
        toys_95_high_ih[i] = dchisq_toys_tih_sorted[int(np.floor(phigh * Ntih) - 1)]
        toys_100_low_nh[i] = dchisq_toys_tnh_sorted[0]
        toys_100_high_nh[i] = dchisq_toys_tnh_sorted[-1]
        toys_100_low_ih[i] = dchisq_toys_tih_sorted[0]
        toys_100_high_ih[i] = dchisq_toys_tih_sorted[-1]

        # A bit of a hack to symmetrize the "errors" to allow use of TGraphErrors to plot
        # both the 95% interval and full interval
        toys_95_nh[i] = toys_95_high_nh[i] - (toys_95_low_nh[i] + toys_95_high_nh[i]) / 2
        toys_95_ih[i] = toys_95_high_ih[i] - (toys_95_low_ih[i] + toys_95_high_ih[i]) / 2
        toys_100_nh[i] = toys_100_high_nh[i] - (toys_100_low_nh[i] + toys_100_high_nh[i]) / 2
        toys_100_ih[i] = toys_100_high_ih[i] - (toys_100_low_ih[i] + toys_100_high_ih[i]) / 2

        data_95_nh[i] = (toys_95_low_nh[i] + toys_95_high_nh[i]) / 2
        data_95_ih[i] = (toys_95_low_ih[i] + toys_95_high_ih[i]) / 2
        data_100_nh[i] = (toys_100_low_nh[i] + toys_100_high_nh[i]) / 2
        data_100_ih[i] = (toys_100_low_ih[i] + toys_100_high_ih[i]) / 2

        # Compute two-sided p-value
        p_nh = get_pvalue_two_sided(dchisq_toys_tnh_sorted, dchisq_data)
        p_ih = get_pvalue_two_sided(dchisq_toys_tih_sorted, dchisq_data)
        print "p-value (NH, %.2f) = %.4f" % (dcp[i], p_nh)
        print "p-value (IH, %.2f) = %.4f" % (dcp[i], p_ih)

    min_dchisq = min(toys_100_low_ih[0], toys_100_low_nh[0])
    max_dchisq = max(toys_100_high_ih[0], toys_100_high_nh[0])
    for i in range(0, N):
        min_dchisq = min(toys_100_low_nh[i], min_dchisq)
        min_dchisq = min(toys_100_low_ih[i], min_dchisq)
        max_dchisq = max(toys_100_high_nh[i], max_dchisq)
        max_dchisq = max(toys_100_high_ih[i], max_dchisq)
    axis = ROOT.TH2D("axis", "axis", 10, -ROOT.TMath.Pi() - ROOT.TMath.Pi() / 8, ROOT.TMath.Pi() + ROOT.TMath.Pi() / 8, \
        10, np.floor(min_dchisq), np.ceil(max_dchisq))
    axis.GetXaxis().SetTitle("#delta_{CP}")
    axis.GetYaxis().SetTitle("#Delta#chi^{2}")
    axis.GetXaxis().SetTitleOffset(1.2)
    axis.GetYaxis().SetTitleOffset(1.2)
    axis.Draw("axis")

    g95_nh = []
    g100_nh = []
    g95_ih = []
    g100_ih = []
    for i in range(0, N):
        # Plot NH full interval
        g100_nh.append(ROOT.TGraphErrors(N, dcp, data_100_nh, np.zeros(N), toys_100_nh))
        g100_nh[i].SetLineColor(ROOT.kBlue)
        g100_nh[i].SetLineStyle(3)
        g100_nh[i].SetMarkerSize(0)
        g100_nh[i].Draw("P SAME")

        # Plot IH full interval
        g100_ih.append(ROOT.TGraphErrors(N, dcp, data_100_ih, np.zeros(N), toys_100_ih))
        g100_ih[i].SetLineColor(ROOT.kRed)
        g100_ih[i].SetLineStyle(3)
        g100_ih[i].Draw("P SAME")
        g100_ih[i].SetMarkerSize(0)

        # Plot NH 95% interval
        g95_nh.append(ROOT.TGraphErrors(N, dcp, data_95_nh, np.zeros(N), toys_95_nh))
        g95_nh[i].SetLineColor(ROOT.kBlue)
        g95_nh[i].SetMarkerColor(ROOT.kBlack)
        g95_nh[i].SetMarkerStyle(0)
        g95_nh[i].SetMarkerSize(0)
        g95_nh[i].Draw("P SAME")

        # Plot IH 95% interval
        g95_ih.append(ROOT.TGraphErrors(N, dcp, data_95_ih, np.zeros(N), toys_95_ih))
        g95_ih[i].SetLineColor(ROOT.kRed)
        g95_ih[i].SetMarkerColor(ROOT.kRed)
        g95_ih[i].SetMarkerStyle(0)
        g95_ih[i].SetMarkerSize(0)
        g95_ih[i].Draw("P SAME")

    # Data line - for a single point
    data_line = ROOT.TLine(- ROOT.TMath.Pi(), dchisq_data, +ROOT.TMath.Pi(), dchisq_data)
    data_line.SetLineColor(ROOT.kBlack)
    data_line.SetLineWidth(2)
    data_line.Draw("SAME")

    textsize = 0.04
    lymin = 0.7
    lymax = lymin + 2 * textsize
    legend = ROOT.TLegend(0.86, lymin, 0.98, lymax)
    legend.SetFillStyle(0)
    legend.SetLineColor(0)
    legend.SetTextSize(textsize)
    legend.AddEntry(g95_nh[0], "True NO", "l")
    legend.AddEntry(g95_ih[0], "True IO", "l")
    legend.Draw()

    canvas.SaveAs("fc_gof_all_%f.pdf" % pvalue)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute GoF")
    parser.add_argument("-r", "--root-output", type=str, dest="root_out", help="ROOT histograms of posterior probabilities", required=True)
    parser.add_argument("-o", "--pdf-output", type=str, dest="pdf_out", help="PDF plots of posterior probabilities", required=True)
    parser.add_argument("-d", "--fc-dir", type=str, dest="fc_dir", help="Top-level directory containing FC files", required=True)
    parser.add_argument("-v", "--delta", type=str, dest="delta", help="The true value of delta", required=True)
    parser.add_argument("-n", "--data-nh", type=str, dest="filename_nh", help="ROOT file containing the normal hierarchy data contours", required=True)
    parser.add_argument("-i", "--data-ih", type=str, dest="filename_ih", help="ROOT file containing the inverted hierarchy data contours", required=True)
    parser.add_argument("--print", dest="print_only", action="store_true", help="Print to PDF from existing root file")

    args = parser.parse_args()

    chisq_data_nh, chisq_data_ih = get_data_mh_chisq(args.filename_nh, args.filename_ih)
    dchisq_data = get_dchisq_mh(chisq_data_nh, chisq_data_ih)
    print args.print_only

    if not args.print_only:
        true_points = [args.delta]
        dchisq_toys_tnh_dict = dict()
        dchisq_toys_tih_dict = dict()

        if args.delta.lower() == "all":
            true_points = list(points.keys())
        
        for dcp in true_points:
            chisq_tnh_nh, chisq_tnh_ih = get_toy_mh_chisq(args.fc_dir, dcp, "nh")
            chisq_tih_nh, chisq_tih_ih = get_toy_mh_chisq(args.fc_dir, dcp, "ih")

            dchisq_toys_tnh = get_dchisq_mh(chisq_tnh_nh, chisq_tnh_ih)
            dchisq_toys_tih = get_dchisq_mh(chisq_tih_nh, chisq_tih_ih)

            h_toys_tnh, h_toys_tih = get_histograms_mh(dchisq_toys_tnh, dchisq_toys_tih, dchisq_data)
            f = ROOT.TFile("fc_gof_" + dcp + ".root", "RECREATE")
            h_toys_tnh.Write()
            h_toys_tih.Write()
            f.Close()

            plot_dchisq_single(h_toys_tnh, h_toys_tih, dchisq_data, "fc_gof_" + dcp + ".pdf")

            dchisq_toys_tnh_dict[dcp] = dchisq_toys_tnh
            dchisq_toys_tih_dict[dcp] = dchisq_toys_tih

        if args.delta.lower() == "all":
            plot_dchisq_all(dchisq_toys_tnh_dict, dchisq_toys_tih_dict, dchisq_data, 0.9545)
            plot_dchisq_all(dchisq_toys_tnh_dict, dchisq_toys_tih_dict, dchisq_data, 0.6827)
            plot_dchisq_all(dchisq_toys_tnh_dict, dchisq_toys_tih_dict, dchisq_data, 0.9973)
    else:
        f = ROOT.TFile(args.root_out, "READ")
        h_toys_tnh = f.Get("toys_tnh")
        h_toys_tnh.SetDirectory(0)
        h_toys_tih = f.Get("toys_tih")
        h_toys_tih.SetDirectory(0)
        f.Close()
        plot_dchisq_single(h_toys_tnh, h_toys_tih, dchisq_data, args.pdf_out)

