#!/usr/bin/env python

import ROOT
import sys
import argparse

ROOT.gROOT.SetBatch(1)

def plot(in_filename, out_filename):
    f = ROOT.TFile(in_filename)
    drawn = False
    names = [ key.GetName() for key in f.GetListOfKeys() ]
    names.sort()
    canvas = ROOT.TCanvas("c1","c1", 700, 500)
    canvas.Print(out_filename + "[")
    for name in names:
        objarray = f.Get(name)
        for spline in objarray:
            N = spline.GetNp()
            is_flat = True
            for i in range(N):
                x = ROOT.Double(0)
                y = ROOT.Double(0)
                spline.GetKnot(i, x, y)
                if y != 1.0:
                    is_flat = False
                    break
            if not is_flat:
                spline.SetTitle(name)
                spline.Draw("LP")
                canvas.Modified()
                canvas.Print(out_filename)
    canvas.Print(out_filename + "]")
    f.Close()

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf

def sigma(name):
    if "ca5" in name:
        return 0.0591673
    elif "bgscl" in name:
        return 0.189454
    elif "maqe" in name:
        return 0.0551331
    elif "manff" in name:
        return 0.0451161
    elif "pfo" in name:
        return 0.0727042
    elif "dismpishp" in name:
        return 0.186182
    elif "erpa_A" in name:
        return 0.0558764
    elif "erpa_B" in name:
        return 0.117005
    elif "erpa_D" in name:
        return 0.126062
    elif "erpa_E" in name:
        return 0.354058
    elif "mecpdd" in name:
        return 0.30569

def nominal(name):
    if "ca5" in name:
        return 0.973796 - 1
    elif "bgscl" in name:
        return 0.977381 - 1
    elif "maqe" in name:
        return 0.898567 - 1
    elif "manff" in name:
        return 0.824069 - 1
    elif "pfo" in name:
        return 0.917547 - 1
    elif "dismpishp" in name:
        return 0.509591 - 0
    elif "erpa_A" in name:
        return 0.661699 - 1
    elif "erpa_B" in name:
        return 1.64227 - 1
    elif "erpa_D" in name:
        return 0.994518 - 1
    elif "erpa_E" in name:
        return 0.878244 - 1
    elif "mecpdd" in name:
        return 0.973912 - 1

def niwg(name):
    if "ca5" in name:
        return (-1.0,9999.0)
    elif "bgscl" in name:
        return (-1.0,9999.0)
    elif "maqe" in name:
        return (-1.0,9999.0)
    elif "manff" in name:
        return (-1.0,9999.0)
    elif "pfo" in name:
        return (-0.11,0.27)
    elif "dismpishp" in name:
        return (-9999.0,9999.0)
    elif "erpa_A" in name:
        return (-1.0,9999.0)
    elif "erpa_B" in name:
        return (-1.0,9999.0)
    elif "erpa_D" in name:
        return (-1.0,9999.0)
    elif "erpa_E" in name:
        return (-1.0,9999.0)
    elif "mecpdd" in name:
        return (-1.0,1.0)


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
        if high < 9999:
            axismax = high
        else:
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
    parser.add_argument("-f", "--filename", dest="in_filename", default="", help="Input file containing splines.")
    parser.add_argument("-s", "--systematic", dest="systematic", default="", help="Systematic to generate spline for.")
    parser.add_argument("-m", "--mode", dest="mode", default="ccqe", help="Mode to generate spline for.")
    parser.add_argument("-t", "--etrue-bin", dest="etbin", default="*", help="Etrue bin.")
    parser.add_argument("-e", "--ereco-bin", dest="erbin", default="*", help="Ereco bin.")
    parser.add_argument("-p", "--theta-bin", dest="thetabin", default="*", help="theta bin.")
    parser.add_argument("-v", "--flavour", dest="flavour", default="", help="Neutrino flavour.")
    parser.add_argument("-o", "--out-filename", dest="out_filename", default="", help="Output filename.")
    args = parser.parse_args()
#    state(args.in_filename, args.systematic, args.mode, args.flavour, args.etbin, args.erbin, args.thetabin)
    plot(args.in_filename, args.out_filename)
