#!/usr/bin/env python

import ROOT
import argparse
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

def get_histogram(f, objname):
    histogram = f.Get(objname)
    histogram.SetDirectory(0)
    return histogram

def subtract_spectra(left, right, output):
    samples = [ "1Rmu", "1Re", "MultiRe", "RHC1Rmu", "RHC1Re" ]
    l_histograms = []
    l_file = ROOT.TFile(left, "READ")
    for s in samples:
        l_histograms.append(get_histogram(l_file, "mock_data_" + s + "_0"))
    l_file.Close()

    r_histograms = []
    r_file = ROOT.TFile(right, "READ")
    for s in samples:
        r_histograms.append(get_histogram(r_file, "mock_data_" + s + "_0"))
    r_file.Close()

    o_histograms = []
    for i in range(len(l_histograms)):
        o_histograms.append(l_histograms[i].Clone("o_" + l_histograms[i].GetName()))
        o_histograms[i].Reset()
        for x in range(1, l_histograms[i].GetNbinsX() + 1):
            for y in range(1, l_histograms[i].GetNbinsY() + 1):
                lval = l_histograms[i].GetBinContent(x, y)
                rval = r_histograms[i].GetBinContent(x, y)
                o_histograms[i].SetBinContent(x, y, lval - rval)

    o_file = ROOT.TFile(output, "RECREATE")
    for i, s in enumerate(samples):
        name = "mock_data_" + s + "_0"
        o_histograms[i].Write(name)
        if "Re" in s:
            N = int(o_histograms[i].Integral())
            graph = ROOT.TGraph(N)
            graph.SetName("g" + s)
            graph.SetTitle("g" + s)
            n = 0
            for theta in range(1, o_histograms[i].GetNbinsX() + 1):
                for erec in range(1, o_histograms[i].GetNbinsY() + 1):
                    if o_histograms[i].GetBinContent(theta, erec):
                        graph.SetPoint(n, erec, theta)
                        n += 1
            graph.Write("g" + s)
        print "%s: %.2f" % (s, o_histograms[i].Integral())
    o_file.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute the difference of two spectra (left - right)")
    parser.add_argument("-l", "--left", type=str, help="The name of the spectra file for the left operand ")
    parser.add_argument("-r", "--right", type=str, help="The name of the spectra file for the right operand")
    parser.add_argument("-o", "--output", type=str, help="Name of the file for the resulatant spectra")
    args = parser.parse_args()
    subtract_spectra(args.left, args.right, args.output)
