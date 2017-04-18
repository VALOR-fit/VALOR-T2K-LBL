#!/usr/bin/env python

import ROOT
import sys
import argparse

def plot(filename, var, point):
    f = ROOT.TFile(filename)
    drawn = False
    names = [ key.GetName() for key in f.GetListOfKeys() if var in key.GetName() ]
    for name in names:
        spline = f.Get(name)
        spline.Eval(point)
        if not drawn:
            spline.Draw()
            drawn = True
        else:
            spline.Draw("SAME")
        print name

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot splines")
    parser.add_argument("-f", "--file", dest="filename", default="", help="Input file containing splines.")
    parser.add_argument("-s", "--systematic", dest="systematic", default="", help="Systematic to generate spline for.")
    parser.add_argument("-p", "--point", dest="point", default="0.0", help="Spline evaluation point.")
    args = parser.parse_args()
    plot(args.filename, args.systematic, float(args.point))
    raw_input("Press Enter to continue...")
