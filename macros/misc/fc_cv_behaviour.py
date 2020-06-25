#!/usr/bin/env python

import ROOT
import argparse
import sys
import numpy as np

from os.path import exists
from array import array

ROOT.gROOT.SetBatch(1)
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

def plot_bestfit_distributions(bestfits, critvals, output):
    canvas = ROOT.TCanvas("canvas_bf", "canvas_bf", 1024, 768)
    canvas.SetTopMargin(0.1)
    canvas.SetLeftMargin(0.1)
    canvas.SetBottomMargin(0.1)
    canvas.SetRightMargin(0.2)
    canvas.SetLogy()
    opt = ""
   
    legend = ROOT.TLegend(0.81, 0.89 - 0.03 * len(bestfits), 0.99, 0.89)
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.SetTextSize(0.03);

    for i, histo in enumerate(bestfits):
        histo.GetXaxis().SetTitleOffset(1.2)
        histo.GetYaxis().SetTitleOffset(1.2)
        histo.SetLineColor(i + 1)
        histo.Draw(opt)
        if i == (len(critvals) - 1):
            legend.AddEntry(histo, "#Delta#chi^{2} #geq %.2f" % critvals[i], "l")
        else:
            legend.AddEntry(histo, "#Delta#chi^{2} %.2f - %.2f" % (critvals[i], critvals[i + 1]), "l")
        opt = "SAME"
    legend.Draw()
    canvas.Print(output)


def get_bestfit_distributions(tree, critical_vals):
    tree.GetEntry(0)
    print "Toys:", tree.GetEntries()
    h_bestfits_nh = []
    h_bestfits_ih = []
    h_bestfits_glob = []
    for i, critical in enumerate(critical_vals):
        h_bestfits_nh.append(tree.hchisq1D_nh.Clone("bestfits_nh_%.2f" % critical))
        h_bestfits_nh[i].SetDirectory(0)
        h_bestfits_nh[i].Reset()
        h_bestfits_nh[i].GetYaxis().SetTitle("Number of toys")
        h_bestfits_ih.append(tree.hchisq1D_ih.Clone("bestfits_ih_%.2f" % critical))
        h_bestfits_ih[i].SetDirectory(0)
        h_bestfits_ih[i].Reset()
        h_bestfits_ih[i].GetYaxis().SetTitle("Number of toys")
        h_bestfits_glob.append(tree.hchisq1D_nh.Clone("bestfits_glob_%.2f" % critical))
        h_bestfits_glob[i].SetDirectory(0)
        h_bestfits_glob[i].Reset()
        h_bestfits_glob[i].GetYaxis().SetTitle("Number of toys")
    for toy in tree:
        dchisq = toy.deltachisq0BF_glob
        for i, critical in enumerate(critical_vals):
            cv_min = critical
            if i == (len(critical_vals) - 1):
                cv_max = np.Infinity
            else:
                cv_max = critical_vals[i + 1]
            if dchisq >= cv_min and dchisq < cv_max:
                bin_bf_nh = toy.hchisq1D_nh.GetMinimumBin()
                bin_bf_ih = toy.hchisq1D_ih.GetMinimumBin()
                chisq_true = toy.chisqtrue0
                chisq_bf_nh = toy.hchisq1D_nh.GetBinContent(bin_bf_nh)
                chisq_bf_ih = toy.hchisq1D_ih.GetBinContent(bin_bf_ih)
                h_bestfits_nh[i].Fill(toy.hchisq1D_nh.GetXaxis().GetBinCenter(bin_bf_nh))
                h_bestfits_ih[i].Fill(toy.hchisq1D_ih.GetXaxis().GetBinCenter(bin_bf_ih))
                if chisq_bf_nh < chisq_bf_ih:
                    h_bestfits_glob[i].Fill(toy.hchisq1D_nh.GetXaxis().GetBinCenter(bin_bf_nh))
                else:
                    h_bestfits_glob[i].Fill(toy.hchisq1D_ih.GetXaxis().GetBinCenter(bin_bf_ih))

    for i, critical in enumerate(critical_vals):
        print i, h_bestfits_nh[i].Integral()
    return h_bestfits_nh, h_bestfits_ih, h_bestfits_glob
                
def get_toy_info(fc_dir, true_dcp, true_mh):
    filename = fc_dir + "/" + true_mh + "/fc_deltacp_" + true_dcp + "_" + hierarchy[true_mh] + "_merge.root"
    if not exists(filename):
        print "Error: No such file \'%s\'" % filename
        sys.exit()
    f = ROOT.TFile(filename, "READ")
    tree = get_tree(f, "fitnt")
    dchisq_bf = np.zeros(tree.GetEntries())
    rate_fhc = np.zeros(tree.GetEntries())
    rate_rhc = np.zeros(tree.GetEntries())
    for t, toy in enumerate(tree):
        if not toy.trueMH == hierarchy[true_mh]:
            print "Error: Expected true hierarchy %s in %s, found %s" % (true_mh, filename, true_mh[0])
        dchisq_bf[t] = toy.deltachisq0BF_glob
        rate_fhc[t] = toy.hErecoToy_1Re.Integral() + toy.hErecoToy_MultiRe.Integral()
        rate_rhc[t] = toy.hErecoToy_RHC1Re.Integral()
    f.Close()
    return dchisq_bf, rate_fhc, rate_rhc

def plot_rate_dchisq(dchisq_bf, rate_fhc, rate_rhc):
    fhc_min = np.floor(np.min(rate_fhc))
    fhc_max = np.ceil(np.max(rate_fhc))
    rhc_min = np.floor(np.min(rate_rhc))
    rhc_max = np.ceil(np.max(rate_rhc))
    dchisq_min = np.floor(np.min(dchisq_bf))
    dchisq_max = np.ceil(np.max(dchisq_bf))

    h_fhc = ROOT.TH2D("h_fhc", "h_fhc", int(min(35, fhc_max - fhc_min)), fhc_min, fhc_max, 100, 0, 9)
    h_fhc.GetXaxis().SetTitle("#nu_{e} events")
    h_fhc.GetYaxis().SetTitle("#Delta#chi^{2}")
    h_fhc.SetDirectory(0)
    h_rhc = ROOT.TH2D("h_rhc", "h_rhc", int(min(35, rhc_max - rhc_min)), rhc_min, rhc_max, 100, 0, 9)
    h_rhc.GetXaxis().SetTitle("#bar{#nu}_{e} events")
    h_rhc.GetYaxis().SetTitle("#Delta#chi^{2}")
    h_rhc.SetDirectory(0)
    for t in range(0, len(dchisq_bf)):
        h_fhc.Fill(rate_fhc[t], dchisq_bf[t])
        h_rhc.Fill(rate_rhc[t], dchisq_bf[t])
    return h_fhc, h_rhc

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute GoF")
    parser.add_argument("-r", "--root-output", type=str, dest="root_out", help="ROOT histograms of posterior probabilities", required=True)
    parser.add_argument("-o", "--pdf-output", type=str, dest="pdf_out", help="PDF plots of posterior probabilities", required=True)
    parser.add_argument("-d", "--fc-dir", type=str, dest="fc_dir", help="Top-level directory containing FC files", required=True)
    parser.add_argument("-v", "--delta", type=str, dest="delta", help="The true value of delta", required=True)

    args = parser.parse_args()

    canvas = ROOT.TCanvas("canvas", "canvas", 1024, 768)
    canvas.SetTopMargin(0.1)
    canvas.SetLeftMargin(0.1)
    canvas.SetBottomMargin(0.1)
    canvas.SetRightMargin(0.1)

    dcp = args.delta

    critvals = [4.0, 6.0, 9.0]
#    critvals = [0.0, 1.0, 4.0]
    fc_filename = args.fc_dir + "/nh/fc_deltacp_" + dcp + "_" + hierarchy["nh"] + "_merge.root"
    f = ROOT.TFile(fc_filename, "READ")
    tree = get_tree(f, "fitnt")
    h_bestfit_nh, h_bestfit_ih, h_bestfit_glob = get_bestfit_distributions(tree, critvals)
    f.Close()

    plot_bestfit_distributions(h_bestfit_nh, critvals, "fc_bestfit_distributions_%s_nh.pdf" % dcp)
    plot_bestfit_distributions(h_bestfit_ih, critvals, "fc_bestfit_distributions_%s_ih.pdf" % dcp)
    plot_bestfit_distributions(h_bestfit_glob, critvals, "fc_bestfit_distributions_%s_global.pdf" % dcp)

    f = ROOT.TFile("fc_bestfits_" + dcp + ".root", "RECREATE")
    for h in h_bestfit_nh:
        h.Write()
    for h in h_bestfit_ih:
        h.Write()
    for h in h_bestfit_glob:
        h.Write()
    f.Close()

    dchisq_bf, rate_fhc, rate_rhc = get_toy_info(args.fc_dir, dcp, "nh")
    h_fhc, h_rhc = plot_rate_dchisq(dchisq_bf, rate_fhc, rate_rhc)

    f = ROOT.TFile("fc_cv_nh_" + dcp + ".root", "RECREATE")
    h_fhc.Write()
    h_rhc.Write()
    f.Close()

    h_fhc.SetTitle("True NO, #delta_{CP} = %.2f - FHC" % points[dcp])
    h_fhc.Draw("COLZ")
    canvas.Print("fc_cv_" + dcp + ".pdf(")
    h_rhc.SetTitle("True NO, #delta_{CP} = %.2f - RHC" % points[dcp])
    h_rhc.Draw("COLZ")
    canvas.Print("fc_cv_" + dcp + ".pdf")

    dchisq_bf, rate_fhc, rate_rhc = get_toy_info(args.fc_dir, dcp, "ih")
    h_fhc, h_rhc = plot_rate_dchisq(dchisq_bf, rate_fhc, rate_rhc)

    f = ROOT.TFile("fc_cv_ih_" + dcp + ".root", "RECREATE")
    h_fhc.Write()
    h_rhc.Write()
    f.Close()

    h_fhc.SetTitle("True IO, #delta_{CP} = %.2f - FHC" % points[dcp])
    h_fhc.Draw("COLZ")
    canvas.Print("fc_cv_" + dcp + ".pdf")
    h_rhc.SetTitle("True IO, #delta_{CP} = %.2f - RHC" % points[dcp])
    h_rhc.Draw("COLZ")
    canvas.Print("fc_cv_" + dcp + ".pdf)")

