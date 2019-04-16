#!/usr/bin/env python

import argparse
import ROOT
import array
import numpy as np
import sys

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

def plot_contours_1d(dchisq_objects, labels, pdf_filename, is_t2k_only, is_inverted=False):
    max_val = 0.0
    for dchisq_obj in dchisq_objects:
        max_val = max(max_val, dchisq_obj.GetMaximum())
    max_val = np.ceil(max_val)
    for dchisq_obj in dchisq_objects:
        dchisq_obj.GetYaxis().SetRangeUser(0, max_val)

    canvas = ROOT.TCanvas("canvas1d", "canvas1d", 700, 500)
   
    options = ""
    if len(dchisq_objects) == 3 and not "nova" in pdf_filename:
        line_color = 2
    else:
        line_color = 1
    for i, dchisq_obj in enumerate(dchisq_objects):
        dchisq_obj.SetLineColor(line_color)
        dchisq_obj.Draw("l" + options)
        options = "same"
        line_color += 1
    lymax = 0.89
    lymin = lymax - (1 + len(dchisq_objects)) * 0.03
    legend = ROOT.TLegend(0.15, lymin, 0.3, lymax);
    legend.SetFillStyle(0);
    legend.SetBorderSize(0);
    legend.SetTextFont(132);
    legend.SetTextSize(0.03);
    for i, dchisq_obj in enumerate(dchisq_objects):
        legend.AddEntry(dchisq_obj, labels[i]);

    if "mix" in pdf_filename and not is_t2k_only:
        fc_file = ROOT.TFile("FCgraphs.root", "READ")
        if not is_inverted:
            g_3sigma_p1 = fc_file.Get("gFC3sigma_nh_p1")
            g_3sigma_m1 = fc_file.Get("gFC3sigma_nh_m1")
        else:
            g_3sigma_p1 = fc_file.Get("gFC3sigma_ih_p1")
            g_3sigma_m1 = fc_file.Get("gFC3sigma_ih_m1")
        g_3sigma_p1.SetLineColor(ROOT.kBlack)
        g_3sigma_p1.SetLineStyle(2)
        g_3sigma_m1.SetLineColor(ROOT.kBlack)
        g_3sigma_m1.SetLineStyle(2)
        fc_file.Close()
        g_3sigma_p1.Draw("SAME")
        g_3sigma_m1.Draw("SAME")
        legend.AddEntry(g_3sigma_p1, "~FC 3#sigma region", "l")
    else:
        critval2 = 4 
        critical2 = ROOT.TLine(-np.pi, critval2, np.pi, critval2)
        critical2.SetLineColor(ROOT.kBlack)
        critical2.SetLineStyle(2)
        critical2.Draw()
        critval3 = 9
        critical3 = ROOT.TLine(-np.pi, critval3, np.pi, critval3)
        critical3.SetLineColor(ROOT.kBlack)
        critical3.SetLineStyle(2)
        critical3.Draw();


    legend.Draw()

    canvas.Print(pdf_filename)

def plot_contours_2d(dchisq_objs, labels, theta23, dm23, pdf_filename, is_inverted=False):
    canvas = ROOT.TCanvas("canvas2d", "canvas2d", 700, 500)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.05)
    canvas.SetBottomMargin(0.1)
    lymax = 0.89
    lymin = lymax - (len(dchisq_objs) + 1) * 0.03
    legend = ROOT.TLegend(0.7, lymin, 0.85, lymax);
    legend.SetFillStyle(0);
    legend.SetBorderSize(0);
    legend.SetTextFont(132);
    legend.SetTextSize(0.03);
    if "_{32}" in dchisq_objs[-1].GetYaxis().GetTitle():
        if is_inverted:
            dchisq_objs[-1].GetYaxis().SetTitle("#Deltam^{2}_{13} (eV^{2}c^{-4})")
        else:
            dchisq_objs[-1].GetYaxis().SetTitle("#Deltam^{2}_{32} (eV^{2}c^{-4})")
        
    dchisq_objs[-1].GetYaxis().SetTitleOffset(1.7)
    dchisq_objs[-1].SetTitle("")
    dchisq_objs[-1].Draw("axis")

    multigraph = ROOT.TMultiGraph()
    if len(dchisq_objs) == 3 and not "nova" in pdf_filename:
        line_color = 2
    else:
        line_color = 1
    for i, dchisq_obj in enumerate(dchisq_objs):
        dchisq_obj.SetContour(2, np.array([2.3, 4.61]))
        contours = get_contour(dchisq_obj)
        if len(contours) > 1:
            legend.AddEntry(contours[1], labels[i], "l")
        elif len(contours) > 0:
            legend.AddEntry(contours[0], labels[i], "l")
        for contour in contours:
            contour.SetLineColor(line_color)
            if len(contours) < 2:
                # No 68% contour, set 90% CL style to solid
                contour.SetLineStyle(1)
            multigraph.Add(contour)
        line_color += 1
    multigraph.Draw()

    if "atmos" in pdf_filename:
        input = ROOT.TMarker(theta23, dm23, 23)
        input.SetMarkerColor(12)
        input.Draw("SAME")
        legend.AddEntry(input, "Truth", "p")

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

def get_local_histo(filename):
    in_file = ROOT.TFile(filename, "READ")
    if in_file.GetListOfKeys().Contains("h_deltachisq"):
        obj = in_file.Get("h_deltachisq")
        obj.SetDirectory(0)
        in_file.Close()
        return (obj, 1)
    elif in_file.GetListOfKeys().Contains("dchisq_map_2D"):
        obj = in_file.Get("dchisq_map_2D")
        obj.SetDirectory(0)
        in_file.Close()
        return obj, 2
    else:
        print "Error: No dchisq surface found"
        sys.exit()

def get_global_histos(filename):
    if "normal" in filename:
        filename_nh = filename
        filename_ih = filename.replace("normal", "inverted")
    else:
        filename_ih = filename
        filename_nh = filename.replace("inverted", "normal")
    nh_file = ROOT.TFile(filename_nh, "READ")
    ih_file = ROOT.TFile(filename_ih, "READ")
    if nh_file.GetListOfKeys().Contains("h_deltachisq"):
        obj_nh = nh_file.Get("h_chisq")
        obj_nh.SetDirectory(0)
        obj_ih = ih_file.Get("h_chisq")
        obj_ih.SetDirectory(0)
        nh_file.Close()
        ih_file.Close()
        minval = min(obj_nh.GetMinimum(), obj_ih.GetMinimum())
        for i in range(1, obj_nh.GetNbinsX() + 1):
            obj_nh.SetBinContent(i, obj_nh.GetBinContent(i) - minval)
            obj_ih.SetBinContent(i, obj_ih.GetBinContent(i) - minval)
        return (obj_nh, obj_ih, 1)
    elif nh_file.GetListOfKeys().Contains("dchisq_map_2D"):
        obj_nh = nh_file.Get("chisq_map_2D")
        obj_nh.SetDirectory(0)
        obj_ih = ih_file.Get("chisq_map_2D")
        obj_ih.SetDirectory(0)
        nh_file.Close()
        ih_file.Close()
        minval = min(obj_nh.GetMinimum(), obj_ih.GetMinimum())
        for i in range(1, obj_nh.GetNbinsX() + 1):
            for j in range(1, obj_nh.GetNbinsY() + 1):
                obj_nh.SetBinContent(i, j, obj_nh.GetBinContent(i, j) - minval)
                obj_ih.SetBinContent(i, j, obj_ih.GetBinContent(i, j) - minval)
        return obj_nh, obj_ih, 2
    else:
        print "Error: No dchisq surface found"
        sys.exit()


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--filename", dest="in_filenames", nargs="+", help="Input filename(s) for contour(s).")
    parser.add_argument("-l", "--label", dest="labels", nargs="+", help="Legend label(s) for contours(s). This should have the same number of arguments as provided to --filename.")
    parser.add_argument("-o", "--output", dest="out_filename", help="Output PDF filename")
    parser.add_argument("--global", dest="do_global", action="store_true")
    parser.add_argument("--t2k", dest="is_t2k_only", action="store_true")
    parser.add_argument("--theta23", dest="theta23", type=float, help="Input value of theta23.")
    parser.add_argument("--dm23", dest="dm23", type=float, help="Input value of dm23.")
    args = parser.parse_args()
    if len(args.in_filenames) == len(args.labels):
        None
    else:
        print args.in_filenames
        print args.labels
        print "Error: The number of filenames and labels should match"
        sys.exit()

    dchisq_objects = []
    dchisq_nh_objects = []
    dchisq_ih_objects = []
    num_var = 0
    for filename in args.in_filenames:
        if not args.do_global:            
            obj, dim = get_local_histo(filename)
            dchisq_objects.append(obj)
        else:
            obj_nh, obj_ih, dim = get_global_histos(filename)
            dchisq_nh_objects.append(obj_nh)
            dchisq_ih_objects.append(obj_ih)
        if num_var == 0:
            num_var = dim
        elif num_var != dim:
            num_var = -1
    print len(dchisq_objects)

    if num_var == 1:
        if not args.do_global:
            plot_contours_1d(dchisq_objects, args.labels, args.out_filename, args.is_t2k_only)
        else:
            nh_filename = args.out_filename
            ih_filename = args.out_filename.replace("normal", "inverted")
            plot_contours_1d(dchisq_nh_objects, args.labels, nh_filename, args.is_t2k_only)
            plot_contours_1d(dchisq_ih_objects, args.labels, ih_filename, args.is_t2k_only, True)
    elif num_var == 2:
        if not args.do_global:
            plot_contours_2d(dchisq_objects, args.labels, args,theta23, args.dm23, args.out_filename)
        else:
            nh_filename = args.out_filename
            ih_filename = args.out_filename.replace("normal", "inverted")
            plot_contours_2d(dchisq_nh_objects, args.labels, args.theta23, args.dm23, nh_filename)
            plot_contours_2d(dchisq_ih_objects, args.labels, args.theta23, args.dm23, ih_filename, True)
    else:
        print "Error: Not all contours have the same dimensionality"
        sys.exit()


