#!/usr/bin/env python

import argparse
import ROOT
import array
import glob
import numpy as np
import sys

ROOT.gROOT.SetBatch(1)
#ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

def plot_event_distributions(dist1, dist2, out_filename):
    canvas = ROOT.TCanvas("c", "c", 700, 500)
    canvas.SetLeftMargin(0.125)
    canvas.SetBottomMargin(0.125)
    canvas.SetRightMargin(0.125)
    canvas.Print(out_filename + "[")

    for sample in ["1Rmu", "RHC1Rmu", "1Re", "RHC1Re", "MultiRe"]:
        dist1[sample].GetXaxis().SetTitle("Events")
        dist1[sample].GetYaxis().SetTitle("Toys")
        dist1[sample].SetLineColor(ROOT.kBlack)
        dist2[sample].SetLineColor(ROOT.kRed)
        dist1[sample].SetTitle(sample)
    
    for sample in ["1Rmu", "RHC1Rmu", "1Re", "RHC1Re", "MultiRe"]:
        if dist2[sample].GetMaximum() > dist1[sample].GetMaximum():
            dist1[sample].GetYaxis().SetRangeUser(0, dist2[sample].GetMaximum())
        dist1[sample].Draw()
        dist2[sample].Draw("SAME")
        canvas.Print(out_filename)
    canvas.Print(out_filename + "]")

def get_event_distribution(file_pattern):
    file_set = glob.glob(file_pattern)

    h_events = dict()
    h_events["1Rmu"] = ROOT.TH1D("fhc1rmu", "fhc1rmu", 130, 220, 350)
    h_events["RHC1Rmu"] = ROOT.TH1D("rhc1rmu", "rhc1rmu", 100, 90, 190)
    h_events["1Re"] = ROOT.TH1D("fhc1re", "fhc1re", 100, 0, 100)
    h_events["RHC1Re"] = ROOT.TH1D("rhc1re", "rhc1re", 40, 0, 40)
    h_events["MultiRe"] = ROOT.TH1D("multire", "multire", 20, 0, 20)

    for filename in file_set:
        f = ROOT.TFile(filename, "READ")
        for sample in ["1Rmu", "RHC1Rmu", "1Re", "RHC1Re", "MultiRe"]:
            keys = ROOT.TIter(f.GetListOfKeys())
            sample_keys = []
            while True:
                key = keys()
                if keys == ROOT.TIter.End():
                    break
                keyval = key.GetName()
                if "_" + sample in keyval:
                    sample_keys.append(keyval)
            for k in sample_keys:
                h_events[sample].Fill(f.Get(k).Integral())
            h_events[sample].SetDirectory(0)

        f.Close()

    outfile = ROOT.TFile("temp.root", "RECREATE")
    for key, dist in h_events.iteritems():
        dist.Write()
    outfile.Close()

    return h_events

def get_chisq_info(pattern):
    file_set = glob.glob(pattern)

    h_chisq_glob = ROOT.TH1D("chisq_glob", "chisq_glob", 200, 400, 600)
    h_chisq_nh = ROOT.TH1D("chisq_nh", "chisq_glob", 200, 400, 600)
    h_chisq_ih = ROOT.TH1D("chisq_ih", "chisq_glob", 200, 400, 600)

    for filename in file_set:
        f = ROOT.TFile(filename, "READ")
        tree = f.Get("fitnt")
        for event in tree:
            h_chisq_glob.Fill(event.chisqmin_glob)
            h_chisq_nh.Fill(event.chisqmin_nh)
            h_chisq_ih.Fill(event.chisqmin_ih)
        h_chisq_glob.SetDirectory(0)
        h_chisq_nh.SetDirectory(0)
        h_chisq_ih.SetDirectory(0)
        f.Close()
    outfile = ROOT.TFile("temp.root", "RECREATE")
    h_chisq_glob.Write()
    h_chisq_nh.Write()
    h_chisq_ih.Write()
    outfile.Close()
    return (h_chisq_glob, h_chisq_nh, h_chisq_ih)

def plot_chisq_info(chisqs1, chisqs2, out_filename):
    chisq1_glob, chisq1_nh, chisq1_ih = chisqs1
    chisq2_glob, chisq2_nh, chisq2_ih = chisqs2
    chisq1_glob.SetLineColor(ROOT.kBlack)
    chisq2_glob.SetLineColor(ROOT.kBlack)
    chisq1_nh.SetLineColor(ROOT.kBlue)
    chisq2_nh.SetLineColor(ROOT.kBlue)
    chisq1_ih.SetLineColor(ROOT.kRed)
    chisq2_ih.SetLineColor(ROOT.kRed)
    chisq2_glob.SetLineStyle(2)
    chisq2_nh.SetLineStyle(2)
    chisq2_ih.SetLineStyle(2)
    
    canvas = ROOT.TCanvas("c", "c", 700, 500)
    canvas.SetLeftMargin(0.125)
    canvas.SetBottomMargin(0.125)
    canvas.SetRightMargin(0.125)
    chisq1_glob.GetXaxis().SetTitle("#chi^{2}")
    chisq1_glob.GetYaxis().SetTitle("Toys")

    chisq1_glob.Draw()
    chisq1_nh.Draw("SAME")
    chisq1_ih.Draw("SAME")
    chisq2_glob.Draw("SAME")
    chisq2_nh.Draw("SAME")
    chisq2_ih.Draw("SAME")
    
    textsize = 0.04
    lymin = 0.65
    legend = ROOT.TLegend(0.755, lymin, 0.95, lymin + textsize * 6)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(textsize)
    legend.AddEntry(chisq1_glob, "-#pi/2 NH Glob", "l")
    legend.AddEntry(chisq1_nh, "-#pi/2 NH NH", "l")
    legend.AddEntry(chisq1_ih, "-#pi/2 NH IH", "l")
    legend.AddEntry(chisq2_glob, "+#pi/2 IH Glob", "l")
    legend.AddEntry(chisq2_nh, "+#pi/2 IH NH", "l")
    legend.AddEntry(chisq2_ih, "+#pi/2 IH IH", "l")
    legend.Draw()

    canvas.Print(out_filename)

def toys(pattern1, pattern2, out_filename):
#    dists_1 = get_event_distribution(args.dir1 + "/batch1_*/*toys.root")
#    dists_2 = get_event_distribution(args.dir2 + "/batch2_*/*toys.root")
    dists_1 = get_event_distribution(pattern1)
    dists_2 = get_event_distribution(pattern2)

    plot_event_distributions(dists_1, dists_2, out_filename)

def results(pattern1, pattern2, out_filename):
    chisqs_1 = get_chisq_info(pattern1)
    chisqs_2 = get_chisq_info(pattern2)
    plot_chisq_info(chisqs_1, chisqs_2, out_filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Perform checks on FC")
    parser.add_argument("--dir1", dest="dir1", help="A root directory of FC toys and critical values.")
    parser.add_argument("--dir2", dest="dir2", help="Another root directory of FC toys and critical values.")
    parser.add_argument("-o", "--output-filename", dest="out_filename", help="Output PDF file.")
    args = parser.parse_args()

#    toys(args.dir1 + "/batch*/*toys.root", args.dir2 + "/batch*/*toys.root", args.out_filename)
    results(args.dir1 + "/fc_deltacp_mpi2_normal_merge.root", args.dir2 + "/fc_deltacp_ppi2_inverted_merge.root", args.out_filename)

