#!/usr/bin/env python

import ROOT
import sys
import argparse

spline_version = "ebspline2017fitqunv9"
t2k_version = "_t2ksk.14a.neut5.3.2.13a_tuned_v1r1."

max_etrue_bins = 0

modes = {
    0: list(["ccqe"]),
    1: list(["mec"]),
    2: list(["cc1pi"]),
    3: list(["cccoh"]),
    4: list(["ccoth"]),
    5: list(["ncpipm", "ncpiz"]),
    6: list(["nccoh"]),
    7: list(["nc1gamma"]),
    8: list(["ncoth"])
}

def flavour_mapping(in_flavour):
    if "nue_sig" in in_flavour:
        return "numu_x_nue"
    elif "nuebar_sig" in in_flavour:
        return "numubar_x_nuebar"
    elif "nue_" in in_flavour:
        return "nue_x_nue"
    elif "nuebar_" in in_flavour:
        return "nuebar_x_nuebar"
    elif "numu_" in in_flavour:
        return "numu_x_numu"
    elif "numubar_" in in_flavour:
        return "numubar_x_numubar"
    else:
        print "Error: Unrecognised flavour:", in_flavour
        sys.exit(0)
    return ""

def standardise_filename(in_filename):
    tuple = in_filename.split("/")
    if len(tuple) < 3:
        print "Error: Expected path ending rfnu*/nu*_beam/*.root. Found ", in_filename
        sys.exit(0)
    root_dir, sub_dir, filename = tuple[-3], tuple[-2], tuple[-1]  # need only the lowest two directory levels and the filename
    if "bar" in sub_dir:
        current = "m250ka.fine."
    else:
        current = "250ka.fine."
    flavour = flavour_mapping(filename)
    if "numu" in root_dir:  # erec
        selection = "_numuselec"
        out_filename = spline_version + t2k_version + current + flavour + selection + ".root"
    elif "nue" in root_dir: # erec/theta
        if "cc1pi" in sub_dir:
            selection = "_cc1piselec"
        else:
            selection = "_nueselec"
        out_filename = spline_version + "erectheta" + t2k_version + current + flavour + selection + ".root"
    else:
        print "Error: Could not determine binning from filename", in_filename
        sys.exit(0)
    print out_filename
    return out_filename

def get_etrue_bins(in_filename):
    tuple = in_filename.split("/")
    if len(tuple) < 3:
        print "Error: Expected path ending rfnu*/nu*_beam/*.root. Found ", in_filename
        sys.exit(0)
    root_dir, sub_dir, filename = tuple[-3], tuple[-2], tuple[-1]  # need only the lowest two directory levels and the filename
    if "bar" in sub_dir:
        fhc = False
    else:
        fhc = True
    flavour = flavour_mapping(filename)

    if fhc:
        if flavour == "numu_x_numu":
            return 11
        if flavour == "numubar_x_numubar":
            return 5
        if flavour == "nue_x_nue":
            return 7
        if flavour == "nuebar_x_nuebar":
            return 2
        if flavour == "numu_x_nue":
            return 11
        else: # "numubar_x_nuebar"
            return 5
    else:
        if flavour == "numu_x_numu":
            return 5
        if flavour == "numubar_x_numubar":
            return 11
        if flavour == "nue_x_nue":
            return 2
        if flavour == "nuebar_x_nuebar":
            return 7
        if flavour == "numu_x_nue":
            return 5
        else: # "numubar_x_nuebar"
            return 11

def graph_to_spline(in_filename, out_filename):
    fin = ROOT.TFile(in_filename)
    fout = ROOT.TFile(out_filename, "RECREATE")
    names = [ key.GetName() for key in fin.GetListOfKeys() ]
    names.sort()
    for name in names:
        tuple = name.split("_")
        spline_names = []
        if len(tuple) == 7: # erec/theta (prefixes im, enu, binx, biny)
            mode, etrue, erec, theta = tuple[3], tuple[4], tuple[5], tuple[6]
            etrue = str(int(etrue[3:]))
            if int(etrue) >= max_etrue_bins:
                continue
            erec = str(int(erec[4:]) - 1)
            theta = str(int(theta[4:]) -1)
            mapped_modes = modes[int(mode[2:])]
            for mapped_mode in mapped_modes:
                spline_names.append("dev_ebo_" + mapped_mode + "_sp_" + etrue + "_" + erec + "_" + theta)
        elif len(tuple) == 6: # erec
            mode, etrue, erec = tuple[3], tuple[4], tuple[5]
            etrue = str(int(etrue[3:]))
            if int(etrue) >= max_etrue_bins:
                continue
            erec = str(int(erec[4:]) - 1)
            mapped_modes = modes[int(mode[2:])]
            for mapped_mode in mapped_modes:
                spline_names.append("dev_ebo_" + mapped_mode + "_sp_" + etrue + "_" + erec)
        else: # error
            print "Error: Unrecognised graph name:", name
            sys.exit(0)
        graph = fin.Get(name)
        for spline_name in spline_names:
            spline = ROOT.TSpline3(spline_name, graph)
            spline.Write(spline_name)
    fout.Close()
    fin.Close()

import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--file", dest="filename", default="", help="Input file containing splines.")
    args = parser.parse_args()
    max_etrue_bins = get_etrue_bins(args.filename)
    out_filename = standardise_filename(args.filename)
    graph_to_spline(args.filename, out_filename)
