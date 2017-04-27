#!/usr/bin/env python

import ROOT
import sys
import argparse

def plot(filename, var, point):
    f = ROOT.TFile(filename)
    drawn = False
    names = [ key.GetName() for key in f.GetListOfKeys() if var in key.GetName() ]
    names.sort()
    canvas = ROOT.TCanvas("c1","c1", 700, 500)
    canvas.Print("spline_" + str(point) + ".pdf[")
    for name in names:
        spline = f.Get(name)
        spline.Eval(point)
        spline.Draw("LP")
        canvas.Modified()
        canvas.Print("spline_" + str(point) + ".pdf")
    canvas.Print("spline_" + str(point) + ".pdf]")
    for name in names:
        print name

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
def state(filename, var, bin):
    f = ROOT.TFile(filename)
    names = [ key.GetName() for key in f.GetListOfKeys() if var in key.GetName() ]
#    names.sort()
#    for name in names:
#        print name
#    exit(0)
    names = [ name for name in names if "_" + str(bin) + "_0" in name ]
    names.sort()
    pdf = matplotlib.backends.backend_pdf.PdfPages(var + ".pdf")
    for name in names:
        print name
        spline = f.Get(name)
        xmin = spline.GetXmin()
        xmax = spline.GetXmax()
        xspan = xmax - xmin
        axismin = xmin - xspan / 2
        axismax = xmax + xspan / 2
        dx = (axismax - axismin) / 20
        xx = []
        yy = []
        for x in np.arange(axismin, axismax + 1e-6, dx):
            val = spline.Eval(x)
            xx.append(x)
            yy.append(val)
        fig = plt.figure()
        plt.plot(xx, yy)
        plt.axvline(x=xmin, color="red")
        plt.axvline(x=xmax, color="red")
        plt.xlabel("$x$")
        plt.ylabel("wght")
        plt.title(name)
        pdf.savefig(fig)
        plt.close()
    pdf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot splines")
    parser.add_argument("-f", "--file", dest="filename", default="", help="Input file containing splines.")
    parser.add_argument("-s", "--systematic", dest="systematic", default="", help="Systematic to generate spline for.")
    parser.add_argument("-p", "--point", dest="point", default="0.0", help="Spline evaluation point.")
    parser.add_argument("-b", "--bin", dest="bin", default="0", help="Ereco bin.")
    args = parser.parse_args()
    state(args.filename, args.systematic, args.bin)
#    plot(args.filename, args.systematic, float(args.point))
