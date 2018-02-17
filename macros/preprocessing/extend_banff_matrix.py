#!/usr/bin/env python

import ROOT
import argparse
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

def extend_banff(N, banff_matrix, banff_nominal, banff_names):
    size = len(banff_nominal)
    # Copy existing BANFF matrix to new matrix with space for new variables
    new_matrix = ROOT.TMatrixD(N, N)
    for x in range(size):
        for y in range(size):
            new_matrix[y][x] = banff_matrix[y][x]

    # Copy existing nominal values to new vector with space for new variables
    new_nominal = ROOT.TVectorD(N)
    for i in range(size):
        new_nominal[i] = banff_nominal[i]

    # Copy existing parameter names to new TObjArray with space for new variables
    new_names = ROOT.TObjArray(N)
    for name in banff_names:
        new_names.Add(ROOT.TObjString(name))

    return new_matrix, new_nominal, new_names

def copy_eb_to_banff(banff_matrix, banff_nominal, banff_names):
    N = len(banff_nominal) + 1
    new_matrix, new_nominal, new_names = extend_banff(N, banff_matrix, banff_nominal, banff_names)
    eb_cov = 2.0 ** 2
    new_matrix[N - 1][N - 1] = eb_cov
    eb_nominal = 0.0
    new_nominal[N - 1] = eb_nominal
    eb_name = "ebo"
    new_names.Add(ROOT.TObjString(eb_name))
    return new_matrix, new_nominal, new_names

def add_eb_to_banff(banff_filename, outputname):
    banff_file = ROOT.TFile(banff_filename, "READ")

    prefit_matrix = banff_file.Get("prefit_cov")
    prefit_nominal = banff_file.Get("prefit_params")
    postfit_matrix = banff_file.Get("postfit_cov")
    postfit_nominal = banff_file.Get("postfit_params")
    banff_names = banff_file.Get("param_list")

    new_prefit_matrix, new_prefit_nominal, new_names = copy_eb_to_banff(prefit_matrix, prefit_nominal, banff_names)
    new_postfit_matrix, new_postfit_nominal, _ = copy_eb_to_banff(postfit_matrix, postfit_nominal, banff_names)

    new_banff_file = ROOT.TFile(outputname, "RECREATE")
    new_names.Write("param_list", 1)
    new_prefit_nominal.Write("prefit_params")
    new_prefit_matrix.Write("prefit_cov");
    new_postfit_nominal.Write("postfit_params")
    new_postfit_matrix.Write("postfit_cov");

    # Copy objects that don't change
    for mode in [ "sk_numode", "sk_anumode" ]:
        for selection in [ "numu", "numub", "nue", "nueb" ]:
            name = mode + "_" + selection + "_bins"
            banff_file.cd()
            obj = banff_file.Get(name)
            new_banff_file.cd()
            obj.Write(name)

    new_banff_file.Close()
    new_banff_file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Insert the new FSI parameters into the unprocessed BANFF matrix")
    parser.add_argument("-b", "--banff", type=str, help="The filename for the original BANFF matrix", required=True)
    parser.add_argument("-v", "--var", type=str, help="The name of the variable to add", required=True)
    parser.add_argument("--output", type=str, help="Name of the output file", required=True)
    args = parser.parse_args()
    if args.var == "eb":
        add_eb_to_banff(args.banff, args.output)
