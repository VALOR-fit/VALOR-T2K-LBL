#!/usr/bin/env python

import ROOT
import argparse
import re
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

def is_fit_1d(fit):
    if fit == "dcp" or fit == "th23" or fit == "th13" or fit == "dmsq":
        return True
    else:
        return False

def get_axis_histogram(valor_filename):
    f = ROOT.TFile(valor_filename, "READ")
    axis_histo = f.Get("dchisq_map_2D")
    axis_histo.SetDirectory(0)
    f.Close()
    return axis_histo

def get_axis_labels(fit):
    if fit == "th23":
        return ("sin^{2}#theta_{23}", "#Delta#chi^{2}")
    elif fit == "dmsq":
        return ("#Deltam^{2}_{32} (eV^{2} / c^{4})", "#Delta#chi^{2}")

def get_common_axis_range(p_dchisq, v_dchisq):
    x_min = -sys.float_info.max
    x_max = sys.float_info.max
    if not p_dchisq is None:
        p_min = ROOT.TMath.MinElement(p_dchisq.GetN(), p_dchisq.GetX())
        p_max = ROOT.TMath.MaxElement(p_dchisq.GetN(), p_dchisq.GetX())
        if p_min > x_min:
            x_min = p_min
        if p_max < x_max:
            x_max = p_max
    if not v_dchisq is None:
        xaxis = v_dchisq.GetXaxis()
        if xaxis.GetXmin() > x_min:
            x_min = xaxis.GetXmin()
        if xaxis.GetXmax() < x_max:
            x_max = xaxis.GetXmax()
    return (x_min, x_max)

def make_graph_local_minimum(graph):
    N = graph.GetN()
    xx = np.zeros(N)
    yy = np.zeros(N)
    for j in range(N):
        x = np.zeros(1)
        y = np.zeros(1)
        graph.GetPoint(j, x, y)
        xx[j] = x[0]
        yy[j] = y[0]
    minval = np.min(yy)
    yy -= minval
    graph_out = ROOT.TGraph(N, xx, yy)
    graph_out.SetLineColor(graph.GetLineColor())
    graph_out.SetLineStyle(graph.GetLineStyle())

    return graph_out

def make_histo_local_minimum(histo):
    N = histo.GetNbinsX()
    minval = histo.GetMinimum()
    histo_out = histo.Clone(histo.GetName() + "_clone")
    histo_out.SetDirectory(0)
    histo_out.Reset()
    for i in range(1, N + 1):
        histo_out.SetBinContent(i, histo.GetBinContent(i) - minval)

    return histo_out

def get_ptheta_contour(filename, fit, is_fake):
    """Retrieve the p-theta contour from file.

    Keyword arguments:
        filename -- The name of the file containing the contour
    """
    f = ROOT.TFile(filename, "READ")
    keys = [key.GetName() for key in f.GetListOfKeys()]

    if is_fit_1d(fit):
        if not is_fake:
            mh = "_0"
        else:
            mh = "_1"
        if fit == "th23":
            dchisq = f.Get("T23_1D" + mh)
        else:
            dchisq = f.Get("Dm32_1D" + mh)
        dchisq.SetLineColor(ROOT.kBlack)
        dchisq.SetLineWidth(1)
        ret_val = (dchisq, None, None)

    f.Close()
    return ret_val

def get_valor_contour(filename, fit):
    """Retrieve the VALOR contour from file.

    Keyword arguments:
        filename -- The name of the file containing the contour
    """
    f = ROOT.TFile(filename, "READ")
    keys = [key.GetName() for key in f.GetListOfKeys()]
    if is_fit_1d(fit):
        # 1D contour
        dchisq = f.Get("h_deltachisq")
        dchisq.SetDirectory(0)
        dchisq.SetLineColor(ROOT.kRed)
        ret_val = (dchisq, None, None)
    f.Close()
    return ret_val

def plot(dchisq, x_axis_label, y_axis_label, opt):
    dchisq.SetTitle("")
    dchisq.GetXaxis().SetTitle(x_axis_label)
    dchisq.GetYaxis().SetTitle(y_axis_label)
    dchisq.Draw(opt)

def compare(fit, is_fake, ptheta_filename, valor_filename, output_filename):
    """Produce comparisons of the specified contours.
    
    Keyword arguments:
        fit -- The type of fit
        is_normal -- Normal hierarchy if True, inverted otherwise
        mach3_filename -- The name of the MaCh3 contour file (or None)
        ptheta_filename -- The name of the p-theta contour file (or None)
        valor_filename -- The name of the valor contour file (or None)
        output_filename -- The name of the PDF file to produce
    """
    max_chisq = -sys.float_info.max
    num_fitters = 0
    p_dchisq = v_dchisq = None

    if not ptheta_filename is None:
        num_fitters += 1
        p_contour1, p_contour2, p_bestfit = get_ptheta_contour(ptheta_filename, fit, is_fake)
        if p_bestfit is None:
            p_dchisq = make_graph_local_minimum(p_contour1)
            p_max = ROOT.TMath.MaxElement(p_dchisq.GetN(), p_dchisq.GetY())
            if p_max > max_chisq:
                max_chisq = p_max

    if not valor_filename is None:
        num_fitters += 1
        v_contour1, v_contour2, v_bestfit = get_valor_contour(valor_filename, fit)
        if v_bestfit is None:
            v_dchisq = make_histo_local_minimum(v_contour1)
            if v_dchisq.GetMaximum() > max_chisq:
                max_chisq = v_dchisq.GetMaximum()

    x_axis_label, y_axis_label = get_axis_labels(fit)
    opt = "l"

    canvas = ROOT.TCanvas("canvas", "canvas", 1024, 768)
    canvas.SetTopMargin(0.1)
    canvas.SetLeftMargin(0.15)
    canvas.SetBottomMargin(0.15)
    canvas.SetRightMargin(0.2)
    lymin = 0.7
    lymax = lymin + num_fitters * 0.03
    legend = ROOT.TLegend(0.85, lymin, 0.95, lymax)
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.SetTextFont(132);
    legend.SetTextSize(0.03);

    if is_fit_1d(fit):
        x_min, x_max = get_common_axis_range(p_dchisq, v_dchisq)
        if max_chisq > 30:
             max_chisq = 30.0
        if not valor_filename is None:
            v_dchisq.GetYaxis().SetRangeUser(0, np.ceil(max_chisq))
            v_dchisq.GetXaxis().SetRangeUser(x_min, x_max)
            if fit == "dmsq":
                v_dchisq.GetXaxis().SetNdivisions(7);
            plot(v_dchisq, x_axis_label, y_axis_label, opt)
            legend.AddEntry(v_dchisq, "VALOR", "l")
            opt = "lsame"
        if not ptheta_filename is None:
            p_dchisq.GetYaxis().SetLimits(0, np.ceil(max_chisq))
            p_dchisq.GetXaxis().SetLimits(x_min, x_max)
            plot(p_dchisq, x_axis_label, y_axis_label, opt)
            legend.AddEntry(p_dchisq, "p-theta", "l")
            opt = "lsame"
    legend.Draw()

    canvas.Print(output_filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the BANFF or xsec matrix")
    parser.add_argument("-p", "--p-theta", type=str, dest="ptheta", help="p-theta file specified")
    parser.add_argument("-v", "--valor", type=str, dest="valor", help="VALOR file specified")
    parser.add_argument("-o", "--output", type=str, required=True, help="Name of the PDF output")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--theta23", action="store_true", help="1D theta23 fit")
    group.add_argument("--dmsq", action="store_true", help="1D dmsq fit")
    parser.add_argument("--fake", action="store_true", help="Indicate fake data contour rather than Asimov")
    args = parser.parse_args()
    if not args.ptheta is None or not args.valor is None:
        # At least one group specified
        if args.theta23:
            fit = "th23"
        elif args.dmsq:
            fit = "dmsq"
        compare(fit, args.fake, args.ptheta, args.valor, args.output)
    else:
        print "Error: No contours specified. Countour files should be specified via"
        print "   -p ptheta_file"
        print "   -v valor_file"
        sys.exit()

