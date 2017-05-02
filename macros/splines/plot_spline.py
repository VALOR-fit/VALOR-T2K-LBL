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

def sigma(name):
    if "ca5" in name:
        return 0.057437
    elif "bgscl" in name:
        return 0.189887
    elif "maqe" in name:
        return 0.068549
    elif "manff" in name:
        return 0.04588
    elif "pfo" in name:
        return 0.088927
    elif "dismpishp" in name:
        return 0.184702
    elif "erpa_A" in name:
        return 0.105128
    elif "erpa_B" in name:
        return 0.139964
    elif "erpa_D" in name:
        return 0.13877
    elif "erpa_E" in name:
        return 0.357223
    elif "mecpdd" in name:
        return 0

def nominal(name):
    if "ca5" in name:
        return -0.00564
    elif "bgscl" in name:
        return -0.07577
    elif "maqe" in name:
        return 0.00022
    elif "manff" in name:
        return -0.2009
    elif "pfo" in name:
        return -0.06443
    elif "dismpishp" in name:
        return 0.446271
    elif "erpa_A" in name:
        return -0.20329
    elif "erpa_B" in name:
        return 0.77547
    elif "erpa_D" in name:
        return -0.08457
    elif "erpa_E" in name:
        return -0.14269
    elif "mecpdd" in name:
        return 1

def niwg(name):
    if "ca5" in name:
        return (-1,9999)
    elif "bgscl" in name:
        return (-1,9999)
    elif "maqe" in name:
        return (-1,9999)
    elif "manff" in name:
        return (-1,9999)
    elif "pfo" in name:
        return (-0.11,0.27)
    elif "dismpishp" in name:
        return (-9999,9999)
    elif "erpa_A" in name:
        return (-1,9999)
    elif "erpa_B" in name:
        return (-1,9999)
    elif "erpa_D" in name:
        return (-1,9999)
    elif "erpa_E" in name:
        return (-1,9999)
    elif "mecpdd" in name:
        return (-9999,9999)


def state(filename, var, mode, flav, etbin, erbin, thetabin):
    f = ROOT.TFile(filename)
    if "t2ksk" in filename:
        names = [ key.GetName() for key in f.GetListOfKeys() if var in key.GetName() ]
        names = [ name for name in names if "_sp_" in name ]
        names = [ name for name in names if mode in name ]
    else:
        names = [ key.GetName() for key in f.GetListOfKeys() if var in key.GetName() ]
        names = [ name for name in names if (flav + mode) in name ]
    newnames = []
    for name in names:
        bins = name.split("_")[-3:]
        inc_etbin = etbin == "*" or etbin == bins[0]
        inc_erbin = erbin == "*" or erbin == bins[1]
        inc_thetabin = thetabin == "*" or thetabin == bins[2]
        inc = inc_etbin and inc_erbin and inc_thetabin
        if inc:
            newnames.append(name)
    names = newnames
    names.sort()
    pdf = matplotlib.backends.backend_pdf.PdfPages(var + ".pdf")
    for name in names:
        print name
        spline = f.Get(name)
        xmin = spline.GetXmin()
        xmax = spline.GetXmax()
        xspan = xmax - xmin
        nom = nominal(name)
        err = sigma(name)
        low, high = niwg(name)
        if low > -9999:
            axismin = low
        else:
            axismin = min(nom - 5 * err, xmin - xspan / 2)
        axismax = max(nom + 5 * err, xmax + xspan / 2)
        dx = (axismax - axismin) / 50
        xx = []
        yy = []
        for x in np.arange(axismin, axismax + 1e-6, dx):
            val = spline.Eval(x)
            xx.append(x)
            yy.append(val)
        fig = plt.figure()
        plt.plot(xx, yy, "black")
        plt.axvline(x=nom, color="black", linestyle="dashed")
        plt.axvline(x=xmin, color="green")
        plt.axvline(x=xmax, color="green")
        plt.axvline(x=nom - 3*err, color="blue")
        plt.axvline(x=nom + 3*err, color="blue")
        if low > -9999:
            plt.axvline(x=low, color="red")
        if high < 9999:
            plt.axvline(x=high, color="red")
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
    parser.add_argument("-m", "--mode", dest="mode", default="ccqe", help="Mode to generate spline for.")
    parser.add_argument("-t", "--etrue-bin", dest="etbin", default="*", help="Etrue bin.")
    parser.add_argument("-e", "--ereco-bin", dest="erbin", default="*", help="Ereco bin.")
    parser.add_argument("-p", "--theta-bin", dest="thetabin", default="*", help="theta bin.")
    parser.add_argument("-v", "--flavour", dest="flavour", default="", help="Neutrino flavour.")
    args = parser.parse_args()
    state(args.filename, args.systematic, args.mode, args.flavour, args.etbin, args.erbin, args.thetabin)
#    plot(args.filename, args.systematic, float(args.point))
