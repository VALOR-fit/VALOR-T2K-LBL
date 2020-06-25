#!/usr/bin/env python

import argparse
import ROOT
import array
import numpy as np
import sys

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

def get_contour(filename):
    f = ROOT.TFile(filename)
    dchisq = f.Get("h_deltachisq")
    dchisq.SetDirectory(0)
    f.Close()
    return dchisq

def plot(unsmeared, smeared, ul, uh, sl, sh, output):
    maxy = 4
    canvas = ROOT.TCanvas("c", "c", 1024, 768)
    canvas.SetLeftMargin(0.125)
    canvas.SetBottomMargin(0.125)
    canvas.SetRightMargin(0.125)
    unsmeared.SetLineColor(ROOT.kBlue)
    unsmeared.SetLineWidth(1)
    unsmeared.SetLineStyle(1)
    unsmeared.GetXaxis().SetRangeUser(2.3e-3, 2.7e-3)
    unsmeared.GetYaxis().SetRangeUser(0, maxy)
    smeared.SetLineColor(ROOT.kRed)
    smeared.SetLineWidth(1)
    smeared.SetLineStyle(1)
    unsmeared.Draw()
    smeared.Draw("same")
    lines = []
    for i, pt in enumerate(list([ul, uh, sl, sh])):
        lines.append(ROOT.TLine(pt, 0, pt, maxy))
        lines[i].SetLineWidth(2)
        if i < 2:
            lines[i].SetLineColor(ROOT.kBlue)
        else:
            lines[i].SetLineColor(ROOT.kRed)
        lines[i].Draw("same")
    uy = 2
    uline = ROOT.TLine(ul, uy, uh, uy)
    uline.SetLineColor(ROOT.kBlue)
    uline.SetLineStyle(2)
    uline.SetLineWidth(2)
    uline.Draw("same")
    utext = ROOT.TLatex()
    utext.SetTextSize(0.04)
    utext.SetTextAlign(11)
    utext.DrawLatex(2.48e-3, uy + 0.1, "1.191#times10^{-4}");
    sy = 3
    sline = ROOT.TLine(sl, sy, sh, sy)
    sline.SetLineColor(ROOT.kRed)
    sline.SetLineStyle(2)
    sline.SetLineWidth(2)
    sline.Draw("same")
    stext = ROOT.TLatex()
    stext.SetTextSize(0.04)
    stext.SetTextAlign(11)
    stext.DrawLatex(2.48e-3, sy + 0.1, "1.427#times10^{-4}");

    canvas.Print(output)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Overlay smeared and unsmeared contours")
    parser.add_argument("-u", "--unsmeared-file", dest="unsmeared_filename", help="Unsmeared contour file.")
    parser.add_argument("-s", "--smeared-file", dest="smeared_filename", help="Smeared contour file.")
    parser.add_argument("-o", "--output-filename", dest="out_filename", help="Output contour file.")
    parser.add_argument("--unsmeared-lower", dest="ul", type=float, help="Unsmeared 1-sigma lower bound.")
    parser.add_argument("--unsmeared-upper", dest="uh", type=float, help="Unsmeared 1-sigma upper bound.")
    parser.add_argument("--smeared-lower", dest="sl", type=float, help="Smeared 1-sigma lower bound.")
    parser.add_argument("--smeared-upper", dest="sh", type=float, help="Smeared 1-sigma upper bound.")

    args = parser.parse_args()
    unsmeared = get_contour(args.unsmeared_filename)
    smeared = get_contour(args.smeared_filename)
    plot(unsmeared, smeared, args.ul, args.uh, args.sl, args.sh, args.out_filename)

