#!/usr/bin/env python

import ROOT
import argparse
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

fsi_params = [ "FEFQE", "FEFQEH", "FEFINEL", "FEFABS", "FEFCX" ]
cc_params = [ "CCnu", "CCnubar" ]
xsec_params = [ "CC_norm_nu", "CC_norm_nubar", "FEFQE", "FEFQEH", "FEFINEL", "FEFABS", "FEFCX" ]

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

def make_banff_xsec_names(N, banff_names):
    banff_names_xsec = ROOT.TObjArray(N)
    # Copy existing BANFF params
    for name in banff_names:
        banff_names_xsec.Add(ROOT.TObjString(name))
    # Add xsec params
    for name in xsec_params:
        banff_names_xsec.Add(ROOT.TObjString(name))
    return banff_names_xsec

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

def copy_xsec_to_banff(Nbanff, banff_matrix, banff_nominal, xsec_filename):
    # Get FSI matrix parameters
    xsec_file = ROOT.TFile(xsec_filename, "READ")
    xsec_names = get_names_from_xsec_matrix(xsec_filename)
    xsec_matrix = xsec_file.Get("xsec_cov")
    for x in range(len(xsec_names)):
        for y in range(len(xsec_names)):
            banff_matrix[Nbanff + y][Nbanff + x] = xsec_matrix[y][x]
    # Get xsec matrix nominal
    xsec_nominal = xsec_file.Get("xsec_param_prior")
    for i in range(len(xsec_names)):
        banff_nominal[Nbanff + i] = xsec_nominal[i]
    xsec_file.Close()

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

def get_flux_from_banff(banff_filename):
    banff_file = ROOT.TFile(banff_filename, "READ")

    # Copy flux parameter names
    banff_names = banff_file.Get("param_list")
    Nflux = banff_names.IndexOf(ROOT.TObjString("SKANuModeNueb6"))
    N = Nflux + 1
    flux_names = ROOT.TObjArray()
    for i in range(N):
        flux_names.Add(banff_names[i])

    # Copy flux parameters from BANFF matrix to new matrix
    banff_matrix = banff_file.Get("prefit_cov")
    flux_matrix = ROOT.TMatrixD(N, N)
    for x in range(N):
        for y in range(N):
            flux_matrix[y][x] = banff_matrix[y][x]
    
    banff_nominal = banff_file.Get("prefit_params")
    flux_nominal = ROOT.TVectorD(N)
    for i in range(N):
        flux_nominal[i] = banff_nominal[i]

    return (flux_matrix, flux_nominal, flux_names)

def streamline_xsec(xsec_filename):
    unused_params = [ "pF_C", "2p2h_shape_C", "NC_other_near", "CC_Coh_C" ]

    xsec_file = ROOT.TFile(xsec_filename, "READ")
    xsec_names = ROOT.TObjArray(xsec_file.Get("xsec_param_names"))

    Nxsec = xsec_names.GetEntries()
    unused_indices = []
    for name in unused_params:
        unused_indices.append(xsec_names.IndexOf(ROOT.TObjString(name)))
    unused_indices.sort()

    N = Nxsec - len(unused_indices)
    xsec_matrix = xsec_file.Get("xsec_cov")
    x_offset = 0
    y_offset = 0
    new_matrix = ROOT.TMatrixD(N, N)
    for x in range(Nxsec):
        if x in unused_indices:
            x_offset -= 1
            continue
        y_offset = 0
        for y in range(Nxsec):
            if y in unused_indices:
                y_offset -= 1
                continue
            new_matrix[y + y_offset][x + x_offset] = xsec_matrix[y][x]
    # Get xsec matrix nominal
    xsec_nominal = xsec_file.Get("xsec_param_prior")
    new_nominal = ROOT.TVectorD(N)
    new_names = ROOT.TObjArray(N)
    i_offset = 0
    for i in range(Nxsec):
        if i in unused_indices:
            i_offset -= 1
            continue
        new_nominal[i + i_offset] = xsec_nominal[i]
        new_names.Add(xsec_names[i])
    xsec_file.Close()
    return (new_matrix, new_nominal, new_names)

def build_banff_from_flux_and_xsec(flux_cov, flux_nom, flux_name, xsec_cov, xsec_nom, xsec_name):
    Nflux = flux_name.GetEntries()
    Nxsec = xsec_name.GetEntries()
    N = Nflux + Nxsec

    # Combine matrices
    banff_cov = ROOT.TMatrixD(N, N)
    for x in range(Nflux):
        for y in range(Nflux):
            banff_cov[y][x] = flux_cov[y][x]
    for x in range(Nxsec):
        for y in range(Nxsec):
            banff_cov[y + Nflux][x + Nflux] = xsec_cov[y][x]

    # Combine nominal
    banff_nom = ROOT.TVectorD(N)
    for i in range(Nflux):
        banff_nom[i] = flux_nom[i]
    for i in range(Nxsec):
        banff_nom[i + Nflux] = xsec_nom[i]

    # Combine names
    banff_name = ROOT.TObjArray(N)
    for i in range(Nflux):
        banff_name.Add(flux_name[i])
    for i in range(Nxsec):
        banff_name.Add(xsec_name[i])

    return (banff_cov, banff_nom, banff_name)

def replace_xsec_in_banff_2017(banff_filename, xsec_filename, outputname):
    flux_cov, flux_nom, flux_name = get_flux_from_banff(banff_filename)
    xsec_cov, xsec_nom, xsec_name = streamline_xsec(xsec_filename)
    banff_cov, banff_nom, banff_name = build_banff_from_flux_and_xsec(flux_cov, flux_nom, flux_name, xsec_cov, xsec_nom, xsec_name)

    output_file = ROOT.TFile(outputname, "RECREATE")
    banff_name.Write("param_list", 1)
    banff_nom.Write("prefit_params")
    banff_cov.Write("prefit_cov");

    # Nothing provided for postfit matrix yet, so just copy the prefit matrix
    # and nominal to ensure that reformatting to valor doesn't crash due to matrix
    # size mismatch
    banff_nom.Write("postfit_params")
    banff_cov.Write("postfit_cov");

    input_file = ROOT.TFile(banff_filename, "READ")
    # Get objects that don't change
    for mode in [ "sk_numode", "sk_anumode" ]:
        for selection in [ "numu", "numub", "nue", "nueb" ]:
            name = mode + "_" + selection + "_bins"
            input_file.cd()
            obj = input_file.Get(name)
            output_file.cd()
            obj.Write(name)

    output_file.Close()
    input_file.Close()

def get_names_from_xsec_matrix(xsec_filename):
    xsec_file = ROOT.TFile(xsec_filename, "READ")
    xsec_names = ROOT.TObjArray(xsec_file.Get("xsec_param_names"))
    xsec_file.Close()
    return xsec_names

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Insert the new FSI parameters into the unprocessed BANFF matrix")
    parser.add_argument("-b", "--banff", type=str, help="The filename for the original BANFF matrix", required=True)
    parser.add_argument("-f", "--fsi", type=str, help="The filename for the xsec matrix with the new FSI parameters", required=True)
    parser.add_argument("--output", type=str, help="Name of the output file", required=True)
    args = parser.parse_args()
#    add_fsi_to_banff(args.banff, args.fsi, args.output)
    replace_xsec_in_banff_2017(args.banff, args.fsi, args.output)
