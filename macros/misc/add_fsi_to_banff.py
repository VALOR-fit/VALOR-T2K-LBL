#!/usr/bin/env python

import ROOT
import argparse
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

fsi_params = [ "FEFQE", "FEFQEH", "FEFINEL", "FEFABS", "FEFCX" ]
cc_params = [ "CCnu", "CCnubar" ]

def make_banff_fsi_names(N, banff_names):
    banff_names_fsi = ROOT.TObjArray(N)
    # Copy existing BANFF params
    for name in banff_names:
        banff_names_fsi.Add(ROOT.TObjString(name))
    # Add CC params
    for name in cc_params:
        banff_names_fsi.Add(ROOT.TObjString(name))
    # Add FSI params
    for name in fsi_params:
        banff_names_fsi.Add(ROOT.TObjString(name))
    return banff_names_fsi

def copy_banff(N, Nbanff, banff_matrix, banff_nominal):
    # Copy existing (prefit) BANFF matrix to new matrix with space for FSI
    banff_matrix_fsi = ROOT.TMatrixD(N, N)
    for x in range(Nbanff):
        for y in range(Nbanff):
            banff_matrix_fsi[y][x] = banff_matrix[y][x]

    # Copy existing (prefit) nominal values to new vector with space for FSI
    banff_nominal_fsi = ROOT.TVectorD(N)
    for i in range(Nbanff):
        banff_nominal_fsi[i] = banff_nominal[i]
    return banff_matrix_fsi, banff_nominal_fsi 

def copy_cc_to_banff(Nbanff, banff_matrix_fsi, banff_nominal_fsi):
    cc_cov = [ [0.02 ** 2, 0.02 * 0.01 * -0.999999], [ 0.02 * 0.01 * -0.999999, 0.01 ** 2 ] ]
#    cc_cov = [ [0.02**2, 0.00], [ 0.0, 0.01**2 ] ]
    for x in range(len(cc_params)):
        for y in range(len(cc_params)):
            banff_matrix_fsi[Nbanff + y][Nbanff + x] = cc_cov[y][x]
    cc_nominal = [ 1.0, 1.0 ]
    for i in range(len(cc_params)):
        banff_nominal_fsi[Nbanff + i] = cc_nominal[i]

def copy_fsi_to_banff(Nbanff, banff_matrix_fsi, banff_nominal_fsi, fsi_filename):
    # Get FSI matrix parameters
    Ncc = Nbanff + len(cc_params)
    fsi_file = ROOT.TFile(fsi_filename, "READ")
    xsec_names = fsi_file.Get("xsec_param_names")
    xsec_matrix = fsi_file.Get("xsec_cov")
    fsi_offset = xsec_names.IndexOf(ROOT.TObjString(fsi_params[0]))
    for x in range(len(fsi_params)):
        for y in range(len(fsi_params)):
            banff_matrix_fsi[Ncc + y][Ncc + x] = xsec_matrix[fsi_offset + y][fsi_offset + x]
    # Get FSI matrix nominal
    fsi_nominal = fsi_file.Get("xsec_param_prior")
    for i in range(len(fsi_params)):
        banff_nominal_fsi[Ncc + i] = fsi_nominal[fsi_offset + i]
    fsi_file.Close()

def add_fsi_to_banff(banff_filename, fsi_filename, outputname):
    banff_file = ROOT.TFile(banff_filename, "READ")
    # Get existing BANFF parameters
    banff_names = banff_file.Get("param_list")
    Nbanff = banff_names.GetEntries()
    # Create new BANFF parameters
    N = Nbanff +len(cc_params) + len(fsi_params)
    banff_names_fsi = make_banff_fsi_names(N, banff_names)

    # Copy existing (prefit) BANFF matrix to new matrix with space for FSI
    banff_matrix = banff_file.Get("prefit_cov")
    banff_nominal = banff_file.Get("prefit_params")
    banff_matrix_fsi, banff_nominal_fsi = copy_banff(N, Nbanff, banff_matrix, banff_nominal)

    copy_cc_to_banff(Nbanff, banff_matrix_fsi, banff_nominal_fsi)
    copy_fsi_to_banff(Nbanff, banff_matrix_fsi, banff_nominal_fsi, fsi_filename)

    banff_fsi_file = ROOT.TFile(outputname, "RECREATE")
    banff_names_fsi.Write("param_list", 1)
    banff_nominal_fsi.Write("prefit_params")
    banff_matrix_fsi.Write("prefit_cov");

    # Nothing provided for postfit matrix yet, so just copy the prefit matrix
    # and nominal to ensure that reformatting to valor doesn't crash due to matrix
    # size mismatch
    banff_nominal_fsi.Write("postfit_params")
    banff_matrix_fsi.Write("postfit_cov");

    # Get objects that don't change
#    for postfit in [ "postfit_params", "postfit_cov" ]:
#        banff_file.cd()
#        obj = banff_file.Get(postfit)
#        banff_fsi_file.cd()
#        obj.Write(postfit)
    for mode in [ "sk_numode", "sk_anumode" ]:
        for selection in [ "numu", "numub", "nue", "nueb" ]:
            name = mode + "_" + selection + "_bins"
            banff_file.cd()
            obj = banff_file.Get(name)
            banff_fsi_file.cd()
            obj.Write(name)

    banff_fsi_file.Close()
    banff_file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Insert the new FSI parameters into the unprocessed BANFF matrix")
    parser.add_argument("-b", "--banff", type=str, help="The filename for the original BANFF matrix", required=True)
    parser.add_argument("-f", "--fsi", type=str, help="The filename for the xsec matrix with the new FSI parameters", required=True)
    parser.add_argument("--output", type=str, help="Name of the output file", required=True)
    args = parser.parse_args()
    add_fsi_to_banff(args.banff, args.fsi, args.output)
