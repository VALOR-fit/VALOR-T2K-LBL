#!/usr/bin/env python

import ROOT
import argparse
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

def officialise(filename_in, filename_out):
    file_in = ROOT.TFile(filename_in, "READ")
    g68 = file_in.Get("g68_valor_0")
    g90 = file_in.Get("g90_valor_0")
    chisq = file_in.Get("chisq_map_2D")
    dchisq = file_in.Get("dchisq_map_2D")
    likelihood = file_in.Get("likelihood_2D")
    bestfit = file_in.Get("best_fit_point_marker")

    file_out = ROOT.TFile(filename_out, "RECREATE")
    g68.SetName("g68_0")
    g90.SetName("g90_0")
    chisq.SetNameTitle("chisq_map_2D", "chisq_map_2D")
    dchisq.SetNameTitle("dchisq_map_2D", "dchisq_map_2D")
    likelihood.SetNameTitle("likelihood_2D", "likelihood_2D")
    g68.Write()
    g90.Write()
    chisq.Write()
    dchisq.Write()
    likelihood.Write()
    bestfit.Write("best_fit_point_marker")
    file_out.Close()

    file_in.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the BANFF or xsec matrix")
    parser.add_argument("-f", "--filename-in", type=str, dest="filename_in", help="File containing the VALOR data")
    parser.add_argument("-o", "--filename-out", type=str, dest="filename_out", help="Official output filename")
    args = parser.parse_args()

    officialise(args.filename_in, args.filename_out)

