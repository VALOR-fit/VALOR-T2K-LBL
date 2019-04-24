#!/usr/bin/env python

import ROOT
import numpy as np
import argparse
import copy
from copy import deepcopy
import math
import array
import sys
import os
import subprocess

ROOT.PyConfig.IgnoreCommandLineOptions = True

env_var = "T2KOSCANA"
try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)
  
valor   = ROOT.valor
util    = valor.util
conv    = valor.conv
opt     = valor.Opt.Instance()
t2kopt  = valor.T2KOpt.Instance()
env     = valor.Env.Instance()
rnd     = valor.Rnd.Instance()
arridx  = util.ArrIdx.Instance()

_max_plots_per_row = 3
_samples = ["1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRe"]

#______________________________________________________________________________________________________________________________

def main(in_filenames, tag):
    nrows = np.ceil(len(in_filenames)/float(_max_plots_per_row))
    print "Splitting canvas into", int(nrows), "rows and", int(_max_plots_per_row), "columns."

    for sample in _samples:
        max = get_max_sample(in_filenames, sample)
        hname = "mock_data_" + sample + "_0"
        canv = ROOT.TCanvas("canv_"+sample, "canv_"+sample, 1200,1000)
        canv.Divide(int(_max_plots_per_row), int(nrows))

        hists = []
        for i in range(0, len(in_filenames)):
            f = ROOT.TFile(in_filenames[i], "read")
            h2D = f.Get(hname)
            if "e" in sample:
                # e-like samples, keep 2D
                hists.append(copy.deepcopy(h2D))
            else:
                # mu-like samples, only want 1D Erec
                h = h2D.ProjectionY(hname+"_1D")
                hists.append(copy.deepcopy(h))


        for i in range(0, len(in_filenames)):
            canv.cd(i+1)
            draw_opt = ""
            if "e" in sample:
                if draw_opt != "":
                    draw_opt += ","
                draw_opt += "COLZ"
            else:
                draw_opt += "HIST"
            hists[i].SetMaximum(max)
            hists[i].Draw(draw_opt)
        outname = tag+".pdf"
        if sample == "1Rmu":
            outname += "("
        elif sample == "MultiRe":
            outname += ")"
        canv.Update()
        canv.cd(0)
        canv.SaveAs(outname, "recreate")

    print "Done!"



def get_max_sample(in_filenames, sample):
    maximum = -999
    hname = "mock_data_" + sample + "_0"
    for name in in_filenames:
        f = ROOT.TFile(name, "read")
        maximum = max(maximum, f.Get(hname).GetMaximum())
    return maximum



if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--filenames", dest="in_filenames", nargs="+", help="Input filename(s) for contour(s).")
    parser.add_argument("-t", "--tag", dest="tag", default="test", help="Output file tag.")

    args = parser.parse_args()
    main(args.in_filenames, args.tag)
