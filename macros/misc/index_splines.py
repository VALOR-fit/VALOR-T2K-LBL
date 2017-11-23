#!/usr/bin/env python

import ROOT
import sys
import argparse

def index(filename):
    f = ROOT.TFile(filename)
    names = [ key.GetName() for key in f.GetListOfKeys() ]
    names.sort()
    spline_set = set()
    for name in names:
        tuple = name.split("_")
        sample, systematic, mode = tuple[1], tuple[2], tuple[3]
        spline_set.add(sample + "_" + systematic + "_" + mode)
    splines = list(spline_set)
    splines.sort()
    for spline in splines:
        print(spline)

def index_single(filename, sample, systematic, mode):
    f = ROOT.TFile(filename)
    names = [ key.GetName() for key in f.GetListOfKeys() ]
    names.sort()
    key = "spl_" + sample + "_" + systematic + "_" + mode
    splines = []
    for name in names:
        if name.startswith(key):
            splines.append(name)
    splines.sort()
    for spline in splines:
        print(spline)

def temp(filename):
    f = ROOT.TFile(filename)
    names = [ key.GetName() for key in f.GetListOfKeys() ]
    names.sort()
    for name in names:
        if "maqe" in name and "_sp_" in name:
            print(name)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--file", dest="filename", default="", help="Input file containing splines.")
    args = parser.parse_args()
#    index(args.filename)
#    index_single(args.filename, "1Rmu", "maqe", "numubarccqe")
    temp(args.filename)
