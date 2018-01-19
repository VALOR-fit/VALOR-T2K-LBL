#!/usr/bin/env python

import ROOT
import sys
import argparse
import os

def explore(filename, systematic=None):
    tuple = filename.split("/")
    print tuple[-3] + "/" + tuple[-2] + "/" + tuple[-1]
    f = ROOT.TFile(filename)
    names = [ key.GetName() for key in f.GetListOfKeys() ]
    if not systematic == None:
        names = [ name for name in names if systematic in name ]
    names.sort()
    for name in names:
        spline = ROOT.TSpline3(name, f.Get(name))
        val = spline.Eval(2)
        if abs(val) > 2E3:
            print "  ", name, spline.Eval(0), val

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--file", dest="filename", default="", help="Input file containing splines.")
    parser.add_argument("-d", "--dir", dest="dirname", default="", help="Top-level directory containing splines.")
    args = parser.parse_args()
    if not args.filename == "":
        explore(args.filename)
    elif not args.dirname == "":
        rootfiles = []
        for root, _, filenames in os.walk(args.dirname):
            for filename in filenames: 
                if filename.endswith(".root"):
                    rootfiles.append(os.path.join(root, filename))
        for filename in rootfiles:
            explore(filename)
    else:
        print "Error: Need a filename or directory"
