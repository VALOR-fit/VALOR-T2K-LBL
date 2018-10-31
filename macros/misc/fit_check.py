#!/usr/bin/env python

import argparse
import ROOT
import array
import glob
import numpy as np
import sys

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

def parse_tree(filename):
    f = ROOT.TFile(filename, "READ")
    tree = f.Get("fitnt")
    for event in tree:
        print "chisq_min", event.chisqmin
        h_chisq = event.hchisq1D
        thebin = h_chisq.GetMinimumBin() 
        print thebin, h_chisq.GetMinimum(), "(", h_chisq.GetBinContent(thebin), ")"
        chi_1rmu = event.hchisq1D_1Rmu.GetBinContent(thebin)
        chi_1re = event.hchisq1D_1Re.GetBinContent(thebin)
        chi_rhc1rmu = event.hchisq1D_RHC1Rmu.GetBinContent(thebin)
        chi_rhc1re = event.hchisq1D_RHC1Re.GetBinContent(thebin)
        chi_multire = event.hchisq1D_MultiRe.GetBinContent(thebin)
        print "1Rmu:", chi_1rmu
        print "1Re:", chi_1re
        print "RHC1Rmu:", chi_rhc1rmu
        print "RHC1Re:", chi_rhc1re
        print "MultiRe:", chi_multire
        print "Sum =", (chi_1rmu + chi_1re + chi_rhc1rmu + chi_rhc1re + chi_multire)
    f.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Perform checks on the fit")
    parser.add_argument("--filename", dest="filename", help="The fit file (non-contour).")
    args = parser.parse_args()

    parse_tree(args.filename)
