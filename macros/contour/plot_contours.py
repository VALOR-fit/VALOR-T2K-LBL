#!/usr/bin/env python

import argparse
import ROOT
import array
import numpy as np
import sys

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

def plot_contours_1d(dchisq_objects, labels, pdf_filename):
    max_val = 0.0
    for dchisq_obj in dchisq_objects:
        max_val = max(max_val, dchisq_obj.GetMaximum())
    max_val = np.ceil(max_val)
    for dchisq_obj in dchisq_objects:
        dchisq_obj.GetYaxis().SetRangeUser(0, max_val)

    canvas = ROOT.TCanvas("canvas1d", "canvas1d", 700, 500)
   
    options = ""
    for i, dchisq_obj in enumerate(dchisq_objects):
        dchisq_obj.SetLineColor(i + 1)
        dchisq_obj.Draw("l" + options)
        options = "same"
    lymax = 0.89
    lymin = lymax - len(dchisq_objects) * 0.03
    legend = ROOT.TLegend(0.15, lymin, 0.3, lymax);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.SetTextFont(132);
    legend.SetTextSize(0.03);
    for i, dchisq_obj in enumerate(dchisq_objects):
        legend.AddEntry(dchisq_obj, labels[i]);
    legend.Draw();

    canvas.Print(pdf_filename)

def plot_contours_2d(dchisq_objs, labels, pdf_filename):
    canvas = ROOT.TCanvas("canvas2d", "canvas2d", 700, 500)
    lymax = 0.89
    lymin = lymax - len(dchisq_objects) * 0.03
    legend = ROOT.TLegend(0.7, lymin, 0.85, lymax);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.SetTextFont(132);
    legend.SetTextSize(0.03);
    dchisq_objs[-1].GetYaxis().SetTitleOffset(1.25)
    dchisq_objs[-1].SetTitle("")
    dchisq_objs[-1].Draw("axis")

    multigraph = ROOT.TMultiGraph()
    for i, dchisq_obj in enumerate(dchisq_objs):
        dchisq_obj.SetContour(2, np.array([2.3, 4.61]))
        contours = get_contour(dchisq_obj)
        legend.AddEntry(contours[0], labels[i], "l");
        for contour in contours:
            contour.SetLineColor(i + 1)
            multigraph.Add(contour)
    multigraph.Draw()
    legend.Draw();
    canvas.Print(pdf_filename)

def get_contour(dchisq_obj):
    canvas = ROOT.TCanvas("c", "c", 0, 0, 700, 500)
    dchisq_obj.Draw("cont3same")
    dchisq_obj.Draw("CONT LIST same")
    canvas.Update()
    contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    graphs = []
    lst = contours.At(0)
    next = ROOT.TIter(lst);
    for g in next:
        g.SetLineStyle(2)
        graphs.append(g.Clone())
    lst = contours.At(1)
    next = ROOT.TIter(lst);
    for g in next:
        g.SetLineStyle(1)
        graphs.append(g.Clone())
    return graphs

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--filename", dest="in_filenames", nargs="+", help="Input filename(s) for contour(s).")
    parser.add_argument("-l", "--label", dest="labels", nargs="+", help="Legend label(s) for contours(s). This should have the same number of arguments as provided to --filename.")
    parser.add_argument("-o", "--output", dest="out_filename", help="Output PDF filename")
    args = parser.parse_args()
    if len(args.in_filenames) == len(args.labels):
        None
    else:
        print "Error: The number of filenames and labels should match"
        sys.exit()

    dchisq_objects = []
    num_var = 0
    for filename in args.in_filenames:
        in_file = ROOT.TFile(filename, "READ")
        if in_file.GetListOfKeys().Contains("h_deltachisq"):
            if num_var == 0:
                num_var = 1
            elif num_var != 1:
                num_var = -1
            dchisq_objects.append(in_file.Get("h_deltachisq"))
            dchisq_objects[-1].SetDirectory(0)
        elif in_file.GetListOfKeys().Contains("dchisq_map_2D"):
            if num_var == 0:
                num_var = 2
            elif num_var != 2:
                num_var == -1
            dchisq_objects.append(in_file.Get("dchisq_map_2D"))
            dchisq_objects[-1].SetDirectory(0)
        else:
            print "Error: No dchisq surface found"
            sys.exit()
        in_file.Close()

    if num_var == 1:
        plot_contours_1d(dchisq_objects, args.labels, args.out_filename)
    elif num_var == 2:
        plot_contours_2d(dchisq_objects, args.labels, args.out_filename)
    else:
        print "Error: Not all contours have the same dimensionality"
        sys.exit()


