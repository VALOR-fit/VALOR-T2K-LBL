#!/usr/bin/env python

import ROOT
import sys
 
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

colours = range(1, 10)

def GetFilenames():
    filenames = []
    for arg in sys.argv[1:]:
        if arg.endswith(".root"):
            filenames.append(arg)
    return filenames

def GetLabels(filenames):
    labels = []
    for arg in sys.argv[1:]:
        if not arg.endswith(".root"):
            labels.append(arg)
    if len(labels) < len(filenames):
        for i in range(len(labels) + 1, len(filenames) + 1):
            labels.append("Set %d" % i)
    elif len(labels) > len(filenames):
        print "More labels than files. Exiting."
        exit(1)
    return labels

def DrawContours(filenames, labels):
    ff = [None] * len(filenames)
    hh = [None] * len(filenames)
    max = 0.0
    canvas = ROOT.TCanvas("canvas","canvas")
    legend = ROOT.TLegend(0.3, 0.6, 0.45, 0.82);
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetTextFont(132)
    legend.SetTextSize(0.05)
    for i, filename in enumerate(filenames):
        ff[i] = ROOT.TFile(filename, "read")
        hh[i] = ff[i].Get("h_deltachisq")
        hh[i].SetTitle("Marginalisation study")
        hh[i].SetLineWidth(2)
        hh[i].SetLineColor(colours[i])
        if hh[i].GetMaximum() > max:
            max = hh[i].GetMaximum()
        hh[i].GetYaxis().SetTitle("-2#Deltaln(L)")
        legend.AddEntry(hh[i], labels[i])
    for i, filename in enumerate(filenames):
        hh[i].GetYaxis().SetRangeUser(0, max)
        if i == 0:
            hh[i].Draw("l")
        else:
            hh[i].Draw("lsame")
    legend.Draw()
    canvas.SaveAs("output.pdf")
    for f in ff:
        f.Close()

def DrawContours2D(filenames, labels):
    ff = [None] * len(filenames)
    hh = [None] * len(filenames)
    hhp = [None] * len(filenames)
    canvas = ROOT.TCanvas("canvas","canvas")
    legend = ROOT.TLegend(0.2, 0.7, 0.35, 0.82);
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.SetTextFont(132)
    legend.SetTextSize(0.05)
    for i, filename in enumerate(filenames):
        ff[i] = ROOT.TFile(filename, "read")
        hh[i] = ff[i].Get("dchisq_map_2D")
        hh[i].SetTitle("Data v Asimov")
        hh[i].SetLineWidth(2)
        hh[i].SetLineColor(colours[i])
        hhp[i] = hh[i].Clone()
        hh[i].SetContour(1)
        hh[i].SetContourLevel(0, 4.61)
        hhp[i].SetLineStyle(2)
        hhp[i].SetContour(1)
        hhp[i].SetContourLevel(0, 2.3)
        legend.AddEntry(hh[i], labels[i])
        if i == 0:
            hh[i].Draw("axis")
            hh[i].Draw("cont3same")
            hhp[i].Draw("cont3same")
        else:
            hh[i].Draw("cont3same")
            hhp[i].Draw("cont3same")
    legend.Draw()
    canvas.SaveAs("output.pdf")
    for f in ff:
        f.Close()

if __name__ == "__main__":
    filenames = GetFilenames()
    labels = GetLabels(filenames)
#    DrawContours(filenames, labels)
    DrawContours2D(filenames, labels)

