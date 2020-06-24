#!/usr/bin/env python

import argparse
import ROOT
import array
import numpy as np
import sys
import copy
from copy import deepcopy

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

def plot_contours_1d(dchisq_objects, labels, pdf_filename, is_t2k_only, colours, styles, prelim, title_x, title_y, line_width, text_size, FC=False, is_inverted=False):
    max_val = 0.0
    min_val = 999
    for dchisq_obj in dchisq_objects:
        max_val = max(max_val, dchisq_obj.GetMaximum())
        min_val = min(min_val, dchisq_obj.GetMinimum())
    max_val = np.ceil(max_val)
    min_val = np.floor(min_val)
    for dchisq_obj in dchisq_objects:
        dchisq_obj.GetYaxis().SetRangeUser(0, max_val)

    canvas = ROOT.TCanvas("canvas1d", "canvas1d", 700, 500)

    if "mixing13" in pdf_filename:
        sinsqth13_nom = 0.0212
        sinsqth13_err = 0.0008
        X1 = sinsqth13_nom - sinsqth13_err
        X2 = sinsqth13_nom + sinsqth13_err
        print "*** test:", X1,min_val,X2,max_val
        breactors = ROOT.TBox(X1,min_val,X2,max_val)
        breactors.SetFillColor(ROOT.kYellow+1)
        dchisq_objects[0].Draw("axis")
        breactors.Draw("same")

    options = ""
    if "mixing13" in pdf_filename:
        options = "same"
    if len(dchisq_objects) == 3 and not "nova" in pdf_filename:
        line_color = 2
    else:
        line_color = 1
    for i, dchisq_obj in enumerate(dchisq_objects):
        if not colours or not styles:
            dchisq_obj.SetLineColor(line_color)
        else:
            # Custom line colours and styles have been specified
            dchisq_obj.SetLineColor(colours[i])
            dchisq_obj.SetLineStyle(styles[i])
        if line_width:
            dchisq_obj.SetLineWidth(line_width)
        if title_x:
            dchisq_obj.GetXaxis().SetTitle(title_x)
        if title_y:
            dchisq_obj.GetYaxis().SetTitle(title_y)
        dchisq_obj.GetYaxis().SetTitleOffset(0.8)
        dchisq_obj.Draw("l" + options)
        options = "same"
        line_color += 1

    if not text_size:
        text_size = 0.03
    lymax = 0.89
    lymin = lymax - (1 + len(dchisq_objects)) * text_size
    lxmin = 0.15
    lxmax = 0.3
    if "mixing13" in pdf_filename:
        lxmin = 0.65
        lxmax = 0.8

    legend = ROOT.TLegend(lxmin, lymin, lxmax, lymax)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(132)
    legend.SetTextSize(text_size)
    for i, dchisq_obj in enumerate(dchisq_objects):
        legend.AddEntry(dchisq_obj, labels[i], "l")

    if "mix" in pdf_filename and not is_t2k_only and FC:
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
        # TODO: Take any number of critvals as a command line input and plot all
        # TODO: Not just for deltacp
        xlow = dchisq_objects[0].GetXaxis().GetBinLowEdge(1)
        xhigh = dchisq_objects[0].GetXaxis().GetBinUpEdge(dchisq_objects[0].GetNbinsX())
        xpos_frac = 0.8
        if "mixing13" in pdf_filename:
            xpos_frac = 0.1
        xpos = dchisq_objects[0].GetXaxis().GetBinLowEdge( int(xpos_frac * dchisq_objects[0].GetNbinsX()) )

        critval1 = 1
        critical1 = ROOT.TLine(xlow, critval1, xhigh, critval1)
        critical1.SetLineColor(ROOT.kBlue)
        critical1.SetLineStyle(1)
        critical1.SetLineWidth(2)
        critical1.Draw()

        latLine1 = ROOT.TLatex(xpos, critval1 + 0.15, str(int(ROOT.TMath.Sqrt(critval1))) + " #sigma")
        latLine1.SetTextColor(ROOT.kBlue)
        latLine1.SetTextSize(0.05)
        latLine1.Draw("same")

        critval2 = 4
        critical2 = ROOT.TLine(xlow, critval2, xhigh, critval2)
        critical2.SetLineColor(ROOT.kBlue)
        critical2.SetLineStyle(1)
        critical2.SetLineWidth(2)
        critical2.Draw()

        latLine2 = ROOT.TLatex(xpos, critval2 + 0.15, str(int(ROOT.TMath.Sqrt(critval2))) + " #sigma")
        latLine2.SetTextColor(ROOT.kBlue)
        latLine2.SetTextSize(0.05)
        latLine2.Draw("same")

    if prelim:
        title = ROOT.TLatex(0.5, 0.91, "T2K Run 1-9 Preliminary")
        title.SetTextSize(0.05)
        title.SetNDC(True)
        title.Draw("same")

    if "mixing13" in pdf_filename:
        legend.AddEntry(breactors, "PDG 2018", "f")

    legend.Draw()

    canvas.Print(pdf_filename)

def plot_contours_2d(dchisq_objs, labels, theta23, dm23, pdf_filename, colours, styles, best_fit_styles, best_fit_labels, prelim, title_x, title_y, line_width, marker_size, text_size, is_inverted=False):
    canvas = ROOT.TCanvas("canvas2d", "canvas2d", 700, 500)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.05)
    canvas.SetBottomMargin(0.1)
    if not text_size:
        text_size = 0.03
    lymax = 0.89
    lymin = lymax - (len(dchisq_objs) + 1) * text_size
    lxmax = 0.85
    lxmin = 0.7
    lx_tsize_scale = 0.3*(1 - text_size/0.03)
    legend = ROOT.TLegend(lxmin + lx_tsize_scale, lymin, lxmax, lymax)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(132)
    legend.SetTextSize(text_size)
    if "_{32}" in dchisq_objs[-1].GetYaxis().GetTitle():
        if is_inverted:
            dchisq_objs[-1].GetYaxis().SetTitle("#Deltam^{2}_{13} (eV^{2}c^{-4})")
        else:
            dchisq_objs[-1].GetYaxis().SetTitle("#Deltam^{2}_{32} (eV^{2}c^{-4})")
    if title_x:
            dchisq_objs[-1].GetXaxis().SetTitle(title_x)
    if title_y:
        dchisq_objs[-1].GetYaxis().SetTitle(title_y)
        
    dchisq_objs[-1].GetYaxis().SetTitleOffset(1.3) #1.7
    dchisq_objs[-1].GetXaxis().SetTitleOffset(0.9)
    dchisq_objs[-1].SetTitle("")
    dchisq_objs[-1].Draw("axis")

    multigraph = ROOT.TMultiGraph()
    if len(dchisq_objs) == 3 and not "nova" in pdf_filename:
        line_color = 2
    else:
        line_color = 1
    for i, dchisq_obj in enumerate(dchisq_objs):
        dchisq_obj.SetContour(1, np.array([4.61]))  #2, [2.3, 4.61]
        contours = get_contour(dchisq_obj)
        if len(contours) > 1:
            legend.AddEntry(contours[1], labels[i], "l")
        elif len(contours) > 0:
            legend.AddEntry(contours[0], labels[i], "l")
        for contour in contours:
            if not colours or not styles:
                contour.SetLineColor(line_color)
                if len(contours) < 2:
                    # No 68% contour, set 90% CL style to solid
                    contour.SetLineStyle(1)
            else:
                # Custom line colours and styles have been specified
                contour.SetLineColor(colours[i])
                contour.SetLineStyle(styles[i])
            if line_width:
                contour.SetLineWidth(line_width)
            multigraph.Add(contour)
        line_color += 1
    multigraph.Draw()

    # Plot best-fit point marker
    if best_fit_labels and best_fit_styles:
        bf_points = []
        bf_lxmin = lxmin - 0.25 + lx_tsize_scale
        bf_lxmax = lxmin + lx_tsize_scale
        bf_lymin = lymax - (len(filter(None, best_fit_labels)) + 1) * text_size
        bf_lymax = lymax
        bf_legend = ROOT.TLegend(bf_lxmin, bf_lymin, bf_lxmax , bf_lymax)
        bf_legend.SetFillStyle(0)
        bf_legend.SetBorderSize(0)
        bf_legend.SetTextFont(132)
        bf_legend.SetTextSize(text_size)
        for i, dchisq_obj in enumerate(dchisq_objs):
            if best_fit_styles[i] != 0:
                minx = ROOT.Long(0)
                miny = ROOT.Long(0)
                minz = ROOT.Long(0)
                dchisq_obj.GetMinimumBin(minx, miny, minz)
                bf_point = ROOT.TMarker(dchisq_obj.GetXaxis().GetBinCenter(minx), dchisq_obj.GetYaxis().GetBinCenter(miny), best_fit_styles[i])
                bf_point.SetMarkerColor(colours[i])
                if marker_size:
                    bf_point.SetMarkerSize(marker_size)
                bf_points.append(copy.deepcopy(bf_point))   # Each iteration of the loop would overwrite the last if we didn't copy
                bf_points[len(bf_points)-1].Draw("same")
                if best_fit_labels[i].replace(" ", "") != "":
                    bf_legend.AddEntry(bf_points[len(bf_points)-1], best_fit_labels[i], "p")
        bf_legend.Draw()

    if "atmos" in pdf_filename:
        input = ROOT.TMarker(theta23, dm23, 23)
        input.SetMarkerColor(12)
        input.Draw("SAME")
        legend.AddEntry(input, "Truth", "p")
    legend.Draw()

    if prelim:
        title = ROOT.TLatex(0.6, 0.96, "T2K Run 1-9 Preliminary")
        title.SetTextSize(0.05)
        title.SetNDC(True)
        title.Draw("same")

    canvas.Update()
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
    if contours.GetEntriesFast() > 1:    
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
    parser.add_argument("-l", "--label", dest="labels", nargs="+", help="Legend label(s) for contours(s). This should have the same number of arguments as provided to --filename, or double if --global is used.")
    parser.add_argument("-c", "--colour", dest="colours", nargs="+", help="Custom line colour(s) for contours(s) (optional). If used, this should have the same number of arguments as provided to --filename, or double if --global is used.")
    parser.add_argument("-s", "--style", dest="styles", nargs="+", help="Custom line style(s) for contours(s) (optional). If used, this should have the same number of arguments as provided to --filename, or double if --global is used.")
    parser.add_argument("--best-fit-label", dest="best_fit_labels", nargs="+", help="Labels for best-fit point(s) for 2D contours(s) (optional). If used, this should have the same number of arguments as provided to --filename, or double if --global is used.")
    parser.add_argument("--best-fit-style", dest="best_fit_styles", nargs="+", help="Marker styles for best fit point(s) for 2D contours(s) (optional). If used, this should have the same number of arguments as provided to --filename, or double if --global is used.")
    parser.add_argument("--title-x", dest="title_x", help="Override default x-axis title")
    parser.add_argument("--title-y", dest="title_y", help="Override default y-axis title")
    parser.add_argument("--line-width", dest="lwidth", type=int, help="Contour line width")
    parser.add_argument("--marker-size", dest="msize", type=float, help="Best-fit marker size")
    parser.add_argument("--text-size", dest="tsize", type=float, help="Legend text size")
    parser.add_argument("-o", "--output", dest="out_filename", help="Output PDF filename")
    parser.add_argument("--global", dest="do_global", action="store_true")
    parser.add_argument("--t2k", dest="is_t2k_only", action="store_true")
    parser.add_argument("--fc", dest="FC", action="store_true")
    parser.add_argument("--prelim", dest="prelim", action="store_true")
    parser.add_argument("--theta23", dest="theta23", type=float, help="Input value of theta23.")
    parser.add_argument("--dm23", dest="dm23", type=float, help="Input value of dm23.")
    args = parser.parse_args()
    if not len(args.in_filenames) == len(args.labels) and not args.do_global:
        print "Error: The number of filenames and labels should match"
        sys.exit()
    elif not 2 * len(args.in_filenames) == len(args.labels) and args.do_global:
        print "Error: The number of labels should be twice the number of filenames"
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
            dchisq_objects.append(obj_nh)
            dchisq_objects.append(obj_ih)
        if num_var == 0:
            num_var = dim
        elif num_var != dim:
            num_var = -1
    print len(dchisq_objects)

    colours = []
    if args.colours:
        colours = [int(i) for i in args.colours]
    styles = []
    if args.styles:
        styles = [int(i) for i in args.styles]
    best_fit_styles = []
    if args.best_fit_styles:
        best_fit_styles = [int(i) for i in args.best_fit_styles]

    if num_var == 1:
        plot_contours_1d(dchisq_objects, args.labels, args.out_filename, args.is_t2k_only, colours, styles, args.prelim, args.title_x, "-2#Deltaln(L)", args.lwidth, args.tsize, args.FC)
    elif num_var == 2:
        plot_contours_2d(dchisq_objects, args.labels, args.theta23, args.dm23, args.out_filename, colours, styles, best_fit_styles, args.best_fit_labels, args.prelim, args.title_x, args.title_y, args.lwidth, args.msize, args.tsize)
    else:
        print "Error: Not all contours have the same dimensionality"
        sys.exit()


