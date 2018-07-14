#!/usr/bin/env python

import ROOT
import argparse
import sys

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

def make_data_file(samples, output_filename):
    output = ROOT.TFile(output_filename, "RECREATE")
    for key, value in samples.items():
        value.SetName("mock_data_" + key + "_0")
        value.SetTitle("mock_data_" + key + "_0")
        print key, value.Integral()
        value.Write()
    output.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create fake data file from individual spectra files")
    parser.add_argument("-f", "--filenames", dest="in_filenames", nargs="+", help="The name of the files containing the sample spectra")
    parser.add_argument("-o", "--output", help="Name of the output root file")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--full", action="store_true", help="Option to build projected data with FHC and RHC projection")
    group.add_argument("--fhc", action="store_true", help="Option to build projected data with FHC projection")
    group.add_argument("--rhc", action="store_true", help="Option to build projected data with RHC projection")
    args = parser.parse_args()

    samples = dict()
    for filename in args.in_filenames:
        if "RHC1Re" in filename:
            sample = "RHC1Re"
        elif "RHC1Rmu" in filename:
            sample = "RHC1Rmu"
        elif "1Rmu" in filename:
            sample = "1Rmu"
        elif "1Re" in filename:
            sample = "1Re"
        elif "MultiRe" in filename:
            sample = "MultiRe"
        else:
            print "Error: Expected sample name in filename; none found."
            sys.exit()
        if args.fhc and sample in ["RHC1Re", "RHC1Rmu"]:
            continue
        elif args.rhc and sample in ["1Re", "MultiRe", "1Rmu"]:
            continue
        key = "curr_ereco_Vs_q2_" + sample
        in_file = ROOT.TFile(filename, "READ")
        if in_file.GetListOfKeys().Contains(key):
            samples[sample] = in_file.Get(key)
            samples[sample].SetDirectory(0)
        else:
            print "Error: No spectrum found"
            sys.exit()
        in_file.Close()

    make_data_file(samples, args.output)
