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
#points = dict(zip(range(9), ["mpi", "m3pi4", "mpi2", "mpi4", "0", "ppi4", "ppi2", "p3pi4", "ppi"]))
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

def read_fc(dir, point):
    chisq = dict(zip(["nh", "ih"], [ dict(), dict() ] ))
    chisq["nh"]["mnh"] = []
    chisq["nh"]["mih"] = []
    chisq["ih"]["mnh"] = []
    chisq["ih"]["mih"] = []
    for mh in hierarchy:
        filename = dir + "/" + mh + "/fc_deltacp_" + point + "_" + hierarchy[mh] + "_merged.root"
        if not exists(filename):
            print "Error: No such file \'%s\'" % filename
            sys.exit()
        f = ROOT.TFile(filename, "READ")
        tree = get_tree(f, "fitnt")
        for entry in tree:
            if not tree.trueMH == hierarchy[mh]:
                print "Error: Expected true hierarchy %s in %s, found %s" % (mh, filename, true_mh[0])
            N = tree.hchisq1D_nh.GetNbinsX()
            nh_vals = np.zeros(N)
            ih_vals = np.zeros(N)
            for i in range(N):
                nh_vals[i] = tree.hchisq1D_nh.GetBinContent(i + 1)
                ih_vals[i] = tree.hchisq1D_ih.GetBinContent(i + 1)
            chisq[mh]["mnh"].append(nh_vals)
            chisq[mh]["mih"].append(ih_vals)
        f.Close()
    return chisq

def get_likelihoods(chisqs):
    likelihoods = dict(zip(["nh", "ih"], [ dict(), dict() ]))
    likelihoods["nh"]["mnh"] = []
    likelihoods["nh"]["mih"] = []
    likelihoods["ih"]["mnh"] = []
    likelihoods["ih"]["mih"] = []

    for mh in hierarchy:
        N = len(chisqs[mh]["mnh"])
        for x in range(N):
            npoints = len(chisqs[mh]["mnh"][x])
            likelihoods[mh]["mnh"].append([])
            likelihoods[mh]["mih"].append([])
            for i in range(npoints):
                chisq = chisqs[mh]["mnh"][x][i]
                likelihoods[mh]["mnh"][x].append(np.exp(-chisq / 2))
                chisq = chisqs[mh]["mih"][x][i]
                likelihoods[mh]["mih"][x].append(np.exp(-chisq / 2))

    return likelihoods

def get_average_likelihood(likelihoods):
    avg_likelihoods = dict( zip(["nh", "ih"], [ dict(), dict() ] ))
    avg_likelihoods["nh"]["mnh"] = []
    avg_likelihoods["nh"]["mih"] = []
    avg_likelihoods["ih"]["mnh"] = []
    avg_likelihoods["ih"]["mih"] = []

    for mh in hierarchy:
        N = len(likelihoods[mh]["mnh"])
        for x in range(N):
            avg_likelihoods[mh]["mnh"].append(0.0)
            avg_likelihoods[mh]["mih"].append(0.0)
            npoints = len(likelihoods[mh]["mnh"][x])
            for i in range(npoints):
                avg_likelihoods[mh]["mnh"][x] += likelihoods[mh]["mnh"][x][i] / npoints
                avg_likelihoods[mh]["mih"][x] += likelihoods[mh]["mih"][x][i] / npoints
    return avg_likelihoods

def get_posterior_probability(likelihoods):
    true_nh = dict()
    true_nh["nh"] = ROOT.TH1D("nh_posterior", "nh_posterior", 100, 0.0, 1.0)
    true_nh["ih"] = ROOT.TH1D("ih_posterior", "ih_posterior", 100, 0.0, 1.0)
    true_nh["nh"].SetDirectory(0)
    true_nh["ih"].SetDirectory(0)

    N = len(likelihoods["nh"]["mnh"])
    for i in range(N):
        p = likelihoods["nh"]["mnh"][i] / (likelihoods["nh"]["mnh"][i] + likelihoods["nh"]["mih"][i])
        true_nh["nh"].Fill(p)
        p = likelihoods["nh"]["mih"][i] / (likelihoods["nh"]["mnh"][i] + likelihoods["nh"]["mih"][i])
        true_nh["ih"].Fill(p)

    true_ih = dict()
    true_ih["nh"] = ROOT.TH1D("nh_posterior", "nh_posterior", 100, 0.0, 1.0)
    true_ih["ih"] = ROOT.TH1D("ih_posterior", "ih_posterior", 100, 0.0, 1.0)
    true_ih["nh"].SetDirectory(0)
    true_ih["ih"].SetDirectory(0)

    N = len(likelihoods["ih"]["mnh"])
    for i in range(N):
        p = likelihoods["ih"]["mnh"][i] / (likelihoods["ih"]["mnh"][i] + likelihoods["ih"]["mih"][i])
        true_ih["nh"].Fill(p)
        p = likelihoods["ih"]["mih"][i] / (likelihoods["ih"]["mnh"][i] + likelihoods["ih"]["mih"][i])
        true_ih["ih"].Fill(p)

    return true_nh, true_ih

def plot_posterior_probabilities(nh_histograms, ih_histograms, dcp, output):
    canvas = ROOT.TCanvas("canvas", "canvas", 1024, 768)
    nh_histograms["nh"].SetTitle("Posterior probability - True NH, #delta=" + str(dcp))
    nh_histograms["nh"].GetXaxis().SetTitle("Posterior probability")
    nh_histograms["nh"].GetYaxis().SetTitle("Number of toys")
    nh_histograms["nh"].SetLineColor(ROOT.kRed)
    nh_histograms["ih"].SetLineColor(ROOT.kBlue)
    nh_histograms["nh"].Draw()
    nh_histograms["ih"].Draw("SAME")
    legend = ROOT.TLegend(0.79, 0.79, 0.89, 0.89)
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.SetTextFont(132);
    legend.SetTextSize(0.05);
    legend.AddEntry(nh_histograms["nh"], "NH", "l")
    legend.AddEntry(nh_histograms["ih"], "IH", "l")
    legend.Draw()
    canvas.Print(output + "(")

    ih_histograms["nh"].SetTitle("Posterior probability - True IH, #delta=" + str(dcp))
    ih_histograms["nh"].GetXaxis().SetTitle("Posterior probability")
    ih_histograms["nh"].GetYaxis().SetTitle("Number of toys")
    ih_histograms["nh"].SetLineColor(ROOT.kRed)
    ih_histograms["ih"].SetLineColor(ROOT.kBlue)
    ih_histograms["nh"].Draw()
    ih_histograms["ih"].Draw("SAME")
    legend = ROOT.TLegend(0.79, 0.79, 0.89, 0.89)
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.SetTextFont(132);
    legend.SetTextSize(0.05);
    legend.AddEntry(nh_histograms["nh"], "NH", "l")
    legend.AddEntry(nh_histograms["ih"], "IH", "l")
    legend.Draw()

    canvas.Print(output + ")")

    Ntoys_nh = nh_histograms["nh"].Integral()
    Nbins = nh_histograms["nh"].GetNbinsX()
    low = nh_histograms["nh"].FindBin(0.95)
    reject_false95_nh = nh_histograms["nh"].Integral(low, Nbins) / float(Ntoys_nh)
    low = nh_histograms["nh"].FindBin(0.95)
    reject_true95_nh = nh_histograms["ih"].Integral(low, Nbins) / float(Ntoys_nh)
    
    low = nh_histograms["nh"].FindBin(0.99)
    reject_false99_nh = nh_histograms["nh"].Integral(low, Nbins) / float(Ntoys_nh)
    low = nh_histograms["nh"].FindBin(0.99)
    reject_true99_nh = nh_histograms["ih"].Integral(low, Nbins) / float(Ntoys_nh)

    Ntoys_ih = ih_histograms["nh"].Integral()
    low = ih_histograms["nh"].FindBin(0.95)
    reject_false95_ih = ih_histograms["ih"].Integral(low, Nbins) / float(Ntoys_ih)
    low = ih_histograms["nh"].FindBin(0.95)
    reject_true95_ih = ih_histograms["nh"].Integral(low, Nbins) / float(Ntoys_ih)
    
    low = ih_histograms["nh"].FindBin(0.99)
    reject_false99_ih = ih_histograms["ih"].Integral(low, Nbins) / float(Ntoys_ih)
    low = ih_histograms["nh"].FindBin(0.99)
    reject_true99_ih = ih_histograms["nh"].Integral(low, Nbins) / float(Ntoys_ih)

    print "True normal hierarchy"
    print "Reject false at 95%%: %.4f" % reject_false95_nh
    print "Reject true at 95%%: %.4f" % reject_true95_nh
    print "Reject false at 99%%: %.4f" % reject_false99_nh
    print "Reject true at 99%%: %.4f" % reject_true99_nh
    
    print "True inverted hierarchy"
    print "Reject false at 95%%: %.4f" % reject_false95_ih
    print "Reject true at 95%%: %.4f" % reject_true95_ih
    print "Reject false at 99%%: %.4f" % reject_false99_ih
    print "Reject true at 99%%: %.4f" % reject_true99_ih

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the BANFF or xsec matrix")
    parser.add_argument("-r", "--root-output", type=str, dest="root_out", help="ROOT histograms of posterior probabilities")
    parser.add_argument("-o", "--pdf-output", type=str, dest="pdf_out", help="PDF plots of posterior probabilities")
    parser.add_argument("-d", "--fc-dir", type=str, dest="fc_dir", help="Top-level directory containing FC files")
    parser.add_argument("-v", "--delta", type=str, dest="delta", help="The true value of delta")
    parser.add_argument("--print", dest="print_only", action="store_true", help="Print to PDF from existing root file")

    args = parser.parse_args()

    if not args.print_only:
        chisqs = read_fc(args.fc_dir, args.delta)
        likelihoods = get_likelihoods(chisqs)
        likelihoods = get_average_likelihood(likelihoods)
        nh_histograms, ih_histograms = get_posterior_probability(likelihoods)

        f = ROOT.TFile(args.root_out, "RECREATE")
        nh_histograms["nh"].Write("nh_true_nh_marg_histogram")
        nh_histograms["ih"].Write("nh_true_ih_marg_histogram")
        ih_histograms["nh"].Write("ih_true_nh_marg_histogram")
        ih_histograms["ih"].Write("ih_true_ih_marg_histogram")
        plot_posterior_probabilities(nh_histograms, ih_histograms, points[args.delta], args.pdf_out)
        f.Close()
    else:
        f = ROOT.TFile(args.root_out, "READ")
        nh_histograms = dict()
        nh_histograms["nh"] = f.Get("nh_true_nh_marg_histogram")
        nh_histograms["ih"] = f.Get("nh_true_ih_marg_histogram")
        ih_histograms = dict()
        ih_histograms["nh"] = f.Get("ih_true_nh_marg_histogram")
        ih_histograms["ih"] = f.Get("ih_true_ih_marg_histogram")
        plot_posterior_probabilities(nh_histograms, ih_histograms, points[args.delta], args.pdf_out)
        f.Close()
