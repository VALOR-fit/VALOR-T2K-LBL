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

def get_axis_labels(fit, is_normal, mach3_present):
    if fit == "dcp":
        return ("#delta_{CP}", "#Delta#chi^{2}")
    elif fit == "th13":
        return ("sin^{2}#theta_{13}", "#Delta#chi^{2}")
    elif fit == "th23":
        return ("sin^{2}#theta_{23}", "#Delta#chi^{2}")
    elif fit == "dmsq":
        if mach3_present or is_normal:
            return ("#Deltam^{2}_{32} (eV^{2} / c^{4})", "#Delta#chi^{2}")
        else:
            return ("#Deltam^{2}_{13} (eV^{2} / c^{4})", "#Delta#chi^{2}")
    elif fit == "app":
        return ("sin^{2}#theta_{13}", "#delta_{CP}")
    elif fit == "disapp":
        if mach3_present or is_normal:
            return ("sin^{2}#theta_{23}", "#Deltam^{2}_{32} (eV^{2} / c^{4})")
        else:
            return ("sin^{2}#theta_{23}", "#Deltam^{2}_{13} (eV^{2}  / c^{4})")

def get_common_axis_range(m_dchisq, p_dchisq, v_dchisq):
    x_min = -sys.float_info.max
    x_max = sys.float_info.max
    if not m_dchisq is None:
        xaxis = m_dchisq.GetXaxis()
        if xaxis.GetXmin() > x_min:
            x_min = xaxis.GetXmin()
        if xaxis.GetXmax() < x_max:
            x_max = xaxis.GetXmax()
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

def shift_dmsq_graph_1d(g):
    N = g.GetN()
    xx = np.zeros(N)
    yy = np.zeros(N)
    for j in range(N):
        x = np.zeros(1)
        y = np.zeros(1)
        g.GetPoint(j, x, y)
        xx[j] = -(x + 7.53E-5)
        yy[j] = y + 7.53E-5
    g_out = ROOT.TGraph(N, xx, yy)
    g_out.SetLineColor(g.GetLineColor())
    g_out.SetLineStyle(g.GetLineStyle())

    return g_out

def shift_dmsq_histo_1d(histo):
    N = histo.GetNbinsX()
    min_x = -(histo.GetXaxis().GetXmax() + 7.53E-5)
    max_x = -(histo.GetXaxis().GetXmin() + 7.53E-5)
    histo_out = ROOT.TH1D(histo.GetName() + "_shift", histo.GetTitle(), N, min_x, max_x)
    histo_out.GetXaxis().SetTitle(histo.GetXaxis().GetTitle())
    histo_out.GetYaxis().SetTitle(histo.GetYaxis().GetTitle())
    for i in range(0, N):
        histo_out.SetBinContent(i + 1, histo.GetBinContent(N - i))
        histo_out.SetLineColor(histo.GetLineColor())
        histo_out.SetLineStyle(histo.GetLineStyle())

    return histo_out

def shift_dmsq_graph_2d(g68, g90, bestfit):
    g68_out = []
    g90_out = []
    bestfit_out = ROOT.TMarker(bestfit.GetX(), -(bestfit.GetY() + 7.53E-5), 34)
    bestfit_out.SetMarkerColor(bestfit.GetMarkerColor())
    for i, g in enumerate(g68):
        N = g.GetN()
        xx = np.zeros(N)
        yy = np.zeros(N)
        for j in range(N):
            x = np.zeros(1)
            y = np.zeros(1)
            g.GetPoint(j, x, y)
            xx[j] = x
            yy[j] = -(y + 7.53E-5)
        g_out = ROOT.TGraph(N, xx, yy)
        g_out.SetLineColor(g.GetLineColor())
        g_out.SetLineStyle(g.GetLineStyle())
        g68_out.append(g_out)
    for i, g in enumerate(g90):
        N = g.GetN()
        xx = np.zeros(N)
        yy = np.zeros(N)
        for j in range(N):
            x = np.zeros(1)
            y = np.zeros(1)
            g.GetPoint(j, x, y)
            xx[j] = x
            yy[j] = -(y + 7.53E-5)
        g_out = ROOT.TGraph(N, xx, yy)
        g_out.SetLineColor(g.GetLineColor())
        g_out.SetLineStyle(g.GetLineStyle())
        g90_out.append(g_out)

    return (g68_out, g90_out, bestfit_out)

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
    minval = histo_out.GetMinimum()
    return histo_out

def get_mach3_1d_contour(f, object_name):
    histo_in = f.Get(object_name)
    dchisq = histo_in.Clone()
    dchisq.SetDirectory(0)
    dchisq.Reset()
    maxval = histo_in.GetBinContent(histo_in.GetMaximumBin())
    for i in range(1, histo_in.GetNbinsX() + 1):
        bincont = histo_in.GetBinContent(i)
        if bincont == 0:
            bincont = 1
        likelihood = -2 * ROOT.TMath.Log(bincont / maxval)
        dchisq.SetBinContent(i, likelihood)
    dchisq.SetLineColor(ROOT.kBlue)
    dchisq.SetFillStyle(0)
    return dchisq

def get_mach3_contour(filename, fit):
    """Retrieve the MaCh3 contour from file.

    Keyword arguments:
        filename -- The name of the file containing the contour
    """
    f = ROOT.TFile(filename, "READ")
    keys = [key.GetName() for key in f.GetListOfKeys()]
    if is_fit_1d(fit):
        # 1D deltacp contour, implies th13 present as well
        if fit == "dcp":
            dchisq = get_mach3_1d_contour(f, "h_dcp")
        elif fit == "th13":
            dchisq = get_mach3_1d_contour(f, "h_th13")
        elif fit == "th23":
            dchisq = get_mach3_1d_contour(f, "h_th23")
        else:
            dchisq = get_mach3_1d_contour(f, "h_dm32")
        ret_val = (dchisq, None, None)
    else:
        # 2D contour
        i = 0
        g68_contours = []
        g90_contours = []
        done68 = False
        done90 = False
        while True:
            key = "g68_conf_" + str(i)
            if key in keys:
                g = f.Get(key)
                g.SetLineColor(ROOT.kBlue)
                g.SetLineWidth(1)
                g.SetLineStyle(2)
                g68_contours.append(g)
            else:
                done68 = True
            key = "g90_conf_" + str(i)
            if key in keys:
                g = f.Get(key)
                g.SetLineColor(ROOT.kBlue)
                g.SetLineWidth(1)
                g90_contours.append(g)
            else:
                done90 = True
            if done68 and done90:
                break
            i += 1
        bestfit = f.Get("Bestfit")
        x = np.zeros(1)
        y = np.zeros(1)
        bestfit.GetPoint(0, x, y)
        bestfit = ROOT.TMarker(x[0], y[0], 34)
        bestfit.SetMarkerColor(ROOT.kBlue)
        ret_val = (g68_contours, g90_contours, bestfit)

    f.Close()
    return ret_val

def get_ptheta_contour(filename, fit, is_normal):
    """Retrieve the p-theta contour from file.

    Keyword arguments:
        filename -- The name of the file containing the contour
    """
    f = ROOT.TFile(filename, "READ")
    keys = [key.GetName() for key in f.GetListOfKeys()]

    if is_fit_1d(fit):
        if is_normal:
            mh = "_0"
        else:
            mh = "_1"
        if fit == "dcp":
            dchisq = f.Get("cont_like_all_conv" + mh)
        elif fit == "th13":
            dchisq = f.Get("T13_1D" + mh)
        elif fit == "th23":
            dchisq = f.Get("T23_1D" + mh)
        else:
            dchisq = f.Get("Dm32_1D" + mh)
        dchisq.SetLineColor(ROOT.kBlack)
        dchisq.SetLineWidth(1)
        ret_val = (dchisq, None, None)
    else:
        # 2D contour
        if is_normal:
            bestfit = f.Get("BestFit_Point_All")
            mh = "_0"
        else:
            bestfit = f.Get("BestFit_Point_All_FDS")
            mh = "_1"
        x = np.zeros(1)
        y = np.zeros(1)
        bestfit.GetPoint(0, x, y)
        bestfit = ROOT.TMarker(x[0], y[0], 34)
        bestfit.SetMarkerColor(ROOT.kBlack)
        i = 0
        g68_contours = []
        done = False
        while not done:
            if "cont_like_all_conv_68_0_0" in keys:
                key = "cont_like_all_conv_68" + mh + "_" + str(i)
            else:
                key = "cont_like_all_conv_68" + mh
                i = -1
            if key in keys:
                g = f.Get(key)
                g.SetLineColor(ROOT.kBlack)
                g.SetLineWidth(1)
                g.SetLineStyle(2)
                g68_contours.append(g)
                if i == -1:
                    done = True
            else:
                done = True
            i += 1
        i = 0
        g90_contours = []
        done = False
        while not done:
            if "cont_like_all_conv_90_0_0" in keys:
                key = "cont_like_all_conv_90" + mh + "_" + str(i)
            else:
                key = "cont_like_all_conv_90" + mh
                i = -1
            if key in keys:
                g = f.Get(key)
                g.SetLineColor(ROOT.kBlack)
                g.SetLineWidth(1)
                g90_contours.append(g)
                if i == -1:
                    done = True
            else:
                done = True
            i += 1
        ret_val = (g68_contours, g90_contours, bestfit)

    f.Close()
    return ret_val

def get_valor_2d_contours(dchisq):
    canvas = ROOT.TCanvas("temp_canvas", "temp_canvas", 1024, 768)
    dchisq.SetContour(2, np.array([2.3, 4.61]))
    dchisq.Draw("cont3same")
    dchisq.Draw("CONT LIST same")
    canvas.Update()
    contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    g68 = []
    lst = contours.At(0)
    next = ROOT.TIter(lst);
    for g in next:
        g.SetLineColor(ROOT.kRed)
        g.SetLineStyle(2)
        g68.append(g.Clone())
    lst = contours.At(1)
    next = ROOT.TIter(lst);
    g90 = []
    for g in next:
        g.SetLineColor(ROOT.kRed)
        g.SetLineStyle(1)
        g90.append(g.Clone())
    return (g68, g90)

def get_valor_contour(filename, fit):
    """Retrieve the VALOR contour from file.

    Keyword arguments:
        filename -- The name of the file containing the contour
    """
    f = ROOT.TFile(filename, "READ")
    keys = [key.GetName() for key in f.GetListOfKeys()]
    if is_fit_1d(fit):
        # 1D contour
        if "_normal.root" in filename:
            f2 = ROOT.TFile(filename.replace("_normal.root", "_inverted.root"), "READ")
            chisq_nh = f.Get("h_chisq")
            chisq_ih = f2.Get("h_chisq")
            minval = np.min(np.array([chisq_nh.GetMinimum(), chisq_ih.GetMinimum()]))
            dchisq = chisq_nh.Clone("h_dchisq")
            dchisq.SetDirectory(0)
            for i in range(1, dchisq.GetNbinsX() + 1):
                dchisq.SetBinContent(i, dchisq.GetBinContent(i) - minval)
            f2.Close()
        else:
            f2 = ROOT.TFile(filename.replace("_inverted.root", "_normal.root"), "READ")
            chisq_ih = f.Get("h_chisq")
            chisq_nh = f2.Get("h_chisq")
            minval = np.min(np.array([chisq_nh.GetMinimum(), chisq_ih.GetMinimum()]))
            dchisq = chisq_ih.Clone("h_dchisq")
            dchisq.SetDirectory(0)
            for i in range(1, dchisq.GetNbinsX() + 1):
                dchisq.SetBinContent(i, dchisq.GetBinContent(i) - minval)
            f2.Close()
        dchisq.SetLineColor(ROOT.kRed)
        ret_val = (dchisq, None, None)
    else:
        # 2D contour
        dchisq = f.Get("dchisq_map_2D")
        dchisq.SetDirectory(0)
        g68_contours, g90_contours = get_valor_2d_contours(dchisq)
        bestfit = f.Get("best_fit_point_marker")
        bestfit.SetMarkerColor(ROOT.kRed)
        bestfit.SetMarkerStyle(34)
        bestfit.SetMarkerSize(1)
        ret_val = (g68_contours, g90_contours, bestfit)
    f.Close()
    return ret_val

def plot(dchisq, x_axis_label, y_axis_label, opt):
    dchisq.SetTitle("")
    dchisq.GetXaxis().SetTitle(x_axis_label)
    dchisq.GetYaxis().SetTitle(y_axis_label)
    dchisq.Draw(opt)
#    dchisq.GetXaxis().SetLabelSize(0.02)
#    dchisq.GetXaxis().SetTitleOffset(0.5)
#    dchisq.GetYaxis().SetLabelSize(0.02)
#    dchisq.GetYaxis().SetTitleOffset(0.5)

def plot2d(g68, g90, bestfit, x_axis_label, y_axis_label, opt):
    for g in g68:
        g.Draw(opt)
        opt = "same"
    for g in g90:
        g.Draw("same")
    bestfit.Draw("same")

def compare(fit, is_normal, mach3_filename, ptheta_filename, valor_filename, plot_fc, use_global_min, output_filename):
    """Produce comparisons of the specified contours.
    
    Keyword arguments:
        fit -- The type of fit
        is_normal -- Normal hierarchy if True, inverted otherwise
        mach3_filename -- The name of the MaCh3 contour file (or None)
        ptheta_filename -- The name of the p-theta contour file (or None)
        valor_filename -- The name of the valor contour file (or None)
        output_filename -- The name of the PDF file to produce
    """
    m_is_1d = True
    p_is_1d = True
    v_is_1d = True
    max_chisq = -sys.float_info.max
    num_fitters = 0
    m_dchisq = p_dchisq = v_dchisq = None
    if not mach3_filename is None:
        num_fitters += 1
        m_contour1, m_contour2, m_bestfit = get_mach3_contour(mach3_filename, fit)
        if m_bestfit is None:
            if not use_global_min:
                m_dchisq = make_histo_local_minimum(m_contour1)
            else:
                m_dchisq = m_contour1
            if m_dchisq.GetMaximum() > max_chisq:
                max_chisq = m_dchisq.GetMaximum()
        else:
            m_is_1d = False
            m_g68 = m_contour1
            m_g90 = m_contour2

    if not ptheta_filename is None:
        num_fitters += 1
        p_contour1, p_contour2, p_bestfit = get_ptheta_contour(ptheta_filename, fit, is_normal)
        if p_bestfit is None:
            if not use_global_min:
                p_dchisq = make_graph_local_minimum(p_contour1)
            else:
                p_dchisq = p_contour1
            p_max = ROOT.TMath.MaxElement(p_dchisq.GetN(), p_dchisq.GetY())
            if p_max > max_chisq:
                max_chisq = p_max
            if fit == "dmsq" and not is_normal and not mach3_filename is None:
                p_dchisq = shift_dmsq_graph_1d(p_dchisq)
        else:
            p_is_1d = False
            p_g68 = p_contour1
            p_g90 = p_contour2
            if fit == "disapp" and not is_normal and not mach3_filename is None:
                p_g68, p_g90, p_bestfit = shift_dmsq_graph_2d(p_g68, p_g90, p_bestfit)

    if not valor_filename is None:
        num_fitters += 1
        v_contour1, v_contour2, v_bestfit = get_valor_contour(valor_filename, fit)
        if v_bestfit is None:
            if not use_global_min:
                v_dchisq = make_histo_local_minimum(v_contour1)
            else:
                v_dchisq = v_contour1
            if v_dchisq.GetMaximum() > max_chisq:
                max_chisq = v_dchisq.GetMaximum()
            if fit == "dmsq" and not is_normal and not mach3_filename is None:
                v_dchisq = shift_dmsq_histo_1d(v_dchisq)
        else:
            v_is_1d = False
            v_g68 = v_contour1
            v_g90 = v_contour2
            if fit == "disapp" and not is_normal and not mach3_filename is None:
                v_g68, v_g90, v_bestfit = shift_dmsq_graph_2d(v_g68, v_g90, v_bestfit)

    x_axis_label, y_axis_label = get_axis_labels(fit, is_normal, not mach3_filename is None)
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
        x_min, x_max = get_common_axis_range(m_dchisq, p_dchisq, v_dchisq)
        if max_chisq > 30:
             max_chisq = 30.0
        if not mach3_filename is None:
            opt = "hist"
            m_dchisq.GetYaxis().SetRangeUser(0, np.ceil(max_chisq))
            m_dchisq.GetXaxis().SetRangeUser(x_min, x_max)
            if fit == "dmsq":
                m_dchisq.GetXaxis().SetNdivisions(7);
            plot(m_dchisq, x_axis_label, y_axis_label, opt)
            legend.AddEntry(m_dchisq, "MaCh3", "l")
            opt = "lsame"
        if not valor_filename is None:
            v_dchisq.GetYaxis().SetRangeUser(0, np.ceil(max_chisq))
            v_dchisq.GetXaxis().SetRangeUser(x_min, x_max)
            if fit == "dmsq":
                v_dchisq.GetXaxis().SetNdivisions(7);
            plot(v_dchisq, x_axis_label, y_axis_label, opt)
            if fit == "dcp" and plot_fc:
                fc_file = ROOT.TFile("FCgraphs.root")
                if is_normal:
                    fc_graph = fc_file.Get("gFC2sigma_nh")
                else:
                    fc_graph = fc_file.Get("gFC2sigma_ih")
                fc_graph.SetLineColor(ROOT.kRed)
                fc_graph.SetLineStyle(2)
                fc_graph.Draw("same")
                fc_file.Close()
            legend.AddEntry(v_dchisq, "VALOR", "l")
            opt = "lsame"
        if not ptheta_filename is None:
            p_dchisq.GetYaxis().SetLimits(0, np.ceil(max_chisq))
            p_dchisq.GetXaxis().SetLimits(x_min, x_max)
            plot(p_dchisq, x_axis_label, y_axis_label, opt)
            if fit == "dcp" and plot_fc:
                fc_file = ROOT.TFile("ptheta_fc.root")
                if is_normal:
                    fc_graph = fc_file.Get("g2sigma_nh")
                else:
                    fc_graph = fc_file.Get("g2sigma_ih")
                fc_graph.SetLineColor(ROOT.kBlack)
                fc_graph.SetLineStyle(2)
                fc_graph.Draw("same")
                fc_file.Close()

            legend.AddEntry(p_dchisq, "p-theta", "l")
            opt = "lsame"
    else:
        axis_histo = get_axis_histogram(valor_filename)
        if fit == "disapp" and not is_normal and not mach3_filename is None:
            nx = axis_histo.GetNbinsX()
            ny = axis_histo.GetNbinsY()
            xmin = axis_histo.GetXaxis().GetXmin()
            xmax = axis_histo.GetXaxis().GetXmax()
            ymin = -(axis_histo.GetYaxis().GetXmax() + 7.53E-5)
            ymax = -(axis_histo.GetYaxis().GetXmin() + 7.53E-5)
            axis_histo = ROOT.TH2D("dummy", "", nx, xmin, xmax, ny, ymin, ymax)
        axis_histo.SetTitle("")
        axis_histo.GetXaxis().SetTitleOffset(2.0)
        axis_histo.GetYaxis().SetTitleOffset(2.0)
        axis_histo.GetXaxis().SetTitle(x_axis_label)
        axis_histo.GetYaxis().SetTitle(y_axis_label)
        axis_histo.Draw("axis")
        if not mach3_filename is None:
            plot2d(m_g68, m_g90, m_bestfit, x_axis_label, y_axis_label, opt)
            legend.AddEntry(m_g90[0], "MaCh3", "l")
            opt = "same"
        if not ptheta_filename is None:
            plot2d(p_g68, p_g90, p_bestfit, x_axis_label, y_axis_label, opt)
            legend.AddEntry(p_g90[0], "p-theta", "l")
            opt = "same"
        if not valor_filename is None:
            plot2d(v_g68, v_g90, v_bestfit, x_axis_label, y_axis_label, opt)
            legend.AddEntry(v_g90[0], "VALOR", "l")
            opt = "same"
    legend.Draw()

    canvas.Print(output_filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the BANFF or xsec matrix")
    parser.add_argument("-m", "--mach3", type=str, dest="mach3", help="MaCh3 file specified")
    parser.add_argument("-p", "--p-theta", type=str, dest="ptheta", help="p-theta file specified")
    parser.add_argument("-v", "--valor", type=str, dest="valor", help="VALOR file specified")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-n", "--normal", action="store_true", help="Normal hierarchy")
    group.add_argument("-i", "--inverted", action="store_true", help="Inverted hierarchy")
    parser.add_argument("-o", "--output", type=str, required=True, help="Name of the PDF output")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--deltacp", action="store_true", help="1D deltacp fit")
    group.add_argument("--theta23", action="store_true", help="1D theta23 fit")
    group.add_argument("--theta13", action="store_true", help="1D theta13 fit")
    group.add_argument("--dmsq", action="store_true", help="1D dmsq fit")
    group.add_argument("--appearance", action="store_true", help="Appearance fit")
    group.add_argument("--disappearance", action="store_true", help="Disappearance fit")
    parser.add_argument("--fc", action="store_true", help="Display FC critical values")
    parser.add_argument("--global_min", action="store_true", help="Use global hierarchy minima")
    args = parser.parse_args()
    if not args.mach3 is None or not args.ptheta is None or not args.valor is None:
        # At least one group specified
        if args.normal:
            is_normal = True
        else:
            is_normal = False
        if args.deltacp:
            fit = "dcp"
        elif args.theta23:
            fit = "th23"
        elif args.theta13:
            fit = "th13"
        elif args.dmsq:
            fit = "dmsq"
        elif args.appearance:
            fit = "app"
        else:
            fit = "disapp"
        compare(fit, is_normal, args.mach3, args.ptheta, args.valor, args.fc, args.global_min, args.output)
    else:
        print "Error: No contours specified. Countour files should be specified via"
        print "   -m mach3_file"
        print "   -p ptheta_file"
        print "   -v valor_file"
        sys.exit()

