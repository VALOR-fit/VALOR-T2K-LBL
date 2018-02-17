#!/usr/bin/env python

import ROOT
import sys
import argparse
from os import listdir
from os.path import isfile, join

def file_listing(dir):
    return [ f for f in listdir(dir) if isfile(join(dir, f)) and f.endswith(".root") ]

def combine_files(main_dir, eb_dir, filename):
    eb_listing = file_listing(eb_dir)
    check = filename
    check = check.replace("Eb_FSI_spline2017", "ebsplines2018")
    out_filename = filename.replace("v9", "v1")
    out_filename = out_filename.replace("2017", "2018")
    fmain = ROOT.TFile(join(main_dir, filename), "READ")
    fout = ROOT.TFile(out_filename, "RECREATE")
    names = [ key.GetName() for key in fmain.GetListOfKeys() ]
    names.sort()
    for name in names:
        if "sp_" in name:
            spline = fmain.Get(name)
            fout.cd()
            spline.Write(spline.GetTitle())
        elif "dev_tmp_0_0" in name:
            spline = fmain.Get(name)
            fout.cd()
            spline.Write(spline.GetTitle())

    fmain.Close()
    if check in eb_listing:
        feb = ROOT.TFile(join(eb_dir, check), "READ")
        names = [ key.GetName() for key in feb.GetListOfKeys() ]
        names.sort()
        for name in names:
            if "_sp_" in name:
                spline = feb.Get(name)
                fout.cd()
                spline.Write(spline.GetTitle())
    fout.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add Eb splines to existing splines")
    parser.add_argument("--maindir", dest="main_dir", default="", help="Path to non-Eb splines.")
    parser.add_argument("--ebdir", dest="eb_dir", default="", help="Path to Eb splines.")
    parser.add_argument("--filename", dest="filename", default="", help="The filename to process.")
    args = parser.parse_args()
    combine_files(args.main_dir, args.eb_dir, args.filename)
