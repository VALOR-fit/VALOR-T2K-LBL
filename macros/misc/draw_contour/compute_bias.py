#!/usr/bin/env python
import argparse
import sys
import ROOT

# Don't let ROOT ruin --help
ROOT.PyConfig.IgnoreCommandLineOptions = True

# load valor
env_var = "T2KOSCANA"
try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()

def get_contour_and_bestfit(filename, objname):
    f = ROOT.TFile(filename, "READ")
    histo = f.Get(objname)
    histo.SetDirectory(0)
    bestfit = histo.GetXaxis().GetBinCenter(histo.GetMinimumBin())
    f.Close()
    return histo, bestfit

def get_interval(contour, bestfit):
    root_finder = valor.RootFinder()
    root_finder.SetGraph(contour)
    root_finder.SetDchiSq(1)
    roots = root_finder.GetRoots() 
    intervals = dict()
    nvals = roots.size()
    nintervals = nvals / 2
    bestfit_interval = 0
    for i in range(nintervals):
        idx = 2 * i
        intervals[i] = (roots[idx], roots[idx + 1])
        if intervals[i][0] <= bestfit <= intervals[i][1]:
            bestfit_interval = i
    return intervals[bestfit_interval]

def bias_bestfit(fakedata_bestfit, asimov_bestfit, asimov_interval):
    asimov_size = 0.5 * (asimov_interval[1] - asimov_interval[0])
    return (fakedata_bestfit - asimov_bestfit) / asimov_size

def bias_interval(fakedata_interval, asimov_interval):
    fakedata_mid = 0.5 * (fakedata_interval[0] + fakedata_interval[1])
    asimov_mid = 0.5 * (asimov_interval[0] + asimov_interval[1])
    asimov_size = 0.5 * (asimov_interval[1] - asimov_interval[0])
    return (fakedata_mid - asimov_mid) / asimov_size

def setup_plot():
    canvas = ROOT.TCanvas("canvas", "canvas", 1024, 768)
    canvas.SetTopMargin(0.1)
    canvas.SetLeftMargin(0.15)
    canvas.SetBottomMargin(0.15)
    canvas.SetRightMargin(0.2)
    lymin = 0.7
    lymax = lymin + 2 * 0.03
    legend = ROOT.TLegend(0.85, lymin, 0.95, lymax)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetTextFont(132)
    legend.SetTextSize(0.03)
    return canvas, legend

def compute_bias(fakedata_filename, asimov_filename, label, output_filename):
    fakedata_contour, fakedata_bestfit = get_contour_and_bestfit(fakedata_filename, "h_deltachisq")
    asimov_contour, asimov_bestfit = get_contour_and_bestfit(asimov_filename, "h_deltachisq")

    canvas, legend = setup_plot()
    asimov_contour.SetTitle("")
    asimov_contour.Draw("axis")
    line = ROOT.TLine(asimov_contour.GetXaxis().GetXmin(), 1.0, asimov_contour.GetXaxis().GetXmax(), 1.0)
    line.SetLineColor(ROOT.kBlue)
    line.Draw("same")
    asimov_contour.SetLineColor(ROOT.kRed)
    if asimov_contour.GetMaximum() > 4:
        asimov_contour.GetYaxis().SetRangeUser(0, 4.0)
    asimov_contour.Draw("lsame")
    fakedata_contour.SetLineColor(ROOT.kBlack)
    fakedata_contour.Draw("lsame")
    legend.AddEntry(fakedata_contour, label, "l")
    legend.AddEntry(asimov_contour, "Asimov", "l")
    legend.Draw()
    canvas.Print(output_filename)

    fd_interval = get_interval(fakedata_contour, fakedata_bestfit)
    av_interval = get_interval(asimov_contour, asimov_bestfit)
    bias_bf = bias_bestfit(fakedata_bestfit, asimov_bestfit, av_interval)
    bias_int = bias_interval(fd_interval, av_interval)
    fd_int_size = fd_interval[1] - fd_interval[0]
    av_int_size = av_interval[1] - av_interval[0]
    with open("bias.info", "a") as f:
        tag = label
        if "th23" in output_filename:
            tag += "_th23"
        else:
            tag += "_dmsq"
        if "_a_" in output_filename:
            tag += "_a"
        else:
            tag += "_b"
        if "smeared" in output_filename:
            tag += "_smeared"
        if "nosyst" in output_filename:
            tag += "_nosyst"
        f.write("%s,%.5f,(%.5f,%.5f),%.5f,(%.5f,%.5f),%.3f,%.3f\n" % (tag, fakedata_bestfit, fd_interval[0], fd_interval[1], asimov_bestfit, av_interval[0], av_interval[1], bias_bf, bias_int))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot fake data and asimov contours and compute the bias")
    parser.add_argument("-f", "--fakedata", type=str, help="The file containing the fake data contour")
    parser.add_argument("-a", "--asimov", type=str, help="The file containing the Asimov contour")
    parser.add_argument("-l", "--label", type=str, help="Label for the fake data contour")
    parser.add_argument("-o", "--output", type=str, help="The name of the output PDF for plots")
    args = parser.parse_args()
    compute_bias(args.fakedata, args.asimov, args.label, args.output)
