#!/usr/bin/env python

import ROOT
import argparse
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
                                                                                                            
def get_1sigma_ptheta(hierarchy):
    dcp = np.array([
        -3.14159265359,-2.35619449019,-1.57079632679,
        -0.78539816339,0.0,0.78539816339,
        1.57079632679,2.35619449019,3.14159265359
    ])
    if hierarchy == "normal":
        central = np.array([
            1.33685,   0.918013,   0.489131,   1.1639,     1.63422,    1.64574,    1.47577,    1.45326,    1.33685
        ])
        err_p = np.array([
            1.35816,   0.933285,   0.506248,   1.17774,    1.65707,    1.66869,    1.49817,    1.47874,    1.35816
        ])
        err_m = np.array([
            1.31051,   0.901412,   0.470563,   1.14611,    1.60675,    1.62627,    1.45899,    1.42888,    1.31051
        ])
    else:
        central = np.array([
            1.72801,   1.68239,    1.49607,    1.45637,    1.39545,    1.07431,    0.75144,    1.32113,    1.72801
        ])
        err_p = np.array([
            1.76231,   1.70355,    1.5157,     1.47835,    1.4169,     1.09118,    0.769979,   1.33888,    1.76231
        ])
        err_m = np.array([
            1.70604,   1.66413,    1.47351,    1.43674,    1.36697,    1.05863,    0.727023,   1.30319,    1.70604
        ])
    gcentral = ROOT.TGraph(len(central), dcp, central)
    gerr_p = ROOT.TGraph(len(err_p), dcp, err_p)
    gerr_m = ROOT.TGraph(len(err_m), dcp, err_m)
    return (gcentral, gerr_m, gerr_p)
                                                                                            
def get_2sigma_ptheta(hierarchy):
    dcp = np.array([
        -3.14159265359,-2.35619449019,-1.57079632679,
        -0.78539816339,0.0,0.78539816339,
        1.57079632679,2.35619449019,3.14159265359
    ])
    if hierarchy == "normal":
        central = np.array([
            3.7301,    2.86577,    2.94686,    3.08528,    4.10718,    4.19564,    3.89917,    3.94545,    3.7301
        ])
        err_p = np.array([
            3.78355,   2.92726,    3.00886,    3.13782,    4.17493,    4.24121,    3.95521,    4.0093,     3.78355
        ])
        err_m = np.array([
            3.67525,   2.81413,    2.90014,    3.03054,    4.04596,    4.14042,    3.82915,    3.88069,    3.67525
        ])
    else:
        central = np.array([
            4.36348,   4.13774,    3.73848,    3.83816,    3.92575,    3.48139,    3.59933,    3.75345,    4.36348
        ])
        err_p = np.array([
            4.42767,   4.18631,    3.79261,    3.90991,    3.96735,    3.54312,    3.66838,    3.80902,    4.42767
        ])
        err_m = np.array([
            4.32913,   4.06989,    3.67521,    3.77397,    3.86069,    3.4345,     3.53091,    3.66713,    4.32913
        ])
    gcentral = ROOT.TGraph(len(central), dcp, central)
    gerr_p = ROOT.TGraph(len(err_p), dcp, err_p)
    gerr_m = ROOT.TGraph(len(err_m), dcp, err_m)
    return (gcentral, gerr_m, gerr_p)

def write_ptheta_graph(g_nh, g_ih):
    f = ROOT.TFile("ptheta_fc.root", "RECREATE")
    g_nh.Write("g2sigma_nh")
    g_ih.Write("g2sigma_ih")
    f.Close()

def get_90cl_ptheta(hierarchy):
    dcp = np.array([
        -3.14159265359,-2.35619449019,-1.57079632679,
        -0.78539816339,0.0,0.78539816339,
        1.57079632679,2.35619449019,3.14159265359
    ])
    if hierarchy == "normal":
        central = np.array([
            2.79803,   1.96558,    1.83203,    2.21087,    3.13907,    3.14253 ,   2.83913,    2.91232,    2.79803
        ])
        err_p = np.array([
            2.82822,   1.98972,    1.87762,    2.23599,    3.18088,    3.18416,    2.87986,    2.95255,    2.82822
        ])
        err_m = np.array([
            2.76579,   1.93043,    1.78571,    2.18437,    3.10556,    3.09595,    2.79109,    2.87349,    2.76579
        ])
    else:
        central = np.array([
            3.3036,    3.14508,    2.83882,    2.8689  ,   2.88333,    2.36627,    2.37207,    2.62994,    3.3036
        ])
        err_p = np.array([
            3.34912,   3.18417,    2.86743,    2.89726,    2.9393,     2.40165,    2.42656,    2.67093,    3.34912
        ])
        err_m = np.array([
            3.26105,   3.10444,    2.79883,    2.83752,    2.84131,    2.30612,    2.32571,    2.59192,    3.26105
        ])
    gcentral = ROOT.TGraph(len(central), dcp, central)
    gerr_p = ROOT.TGraph(len(err_p), dcp, err_p)
    gerr_m = ROOT.TGraph(len(err_m), dcp, err_m)
    return (gcentral, gerr_m, gerr_p)

def get_valor_graph(valor_file, objname):
    return valor_file.Get(objname)

def get_1sigma_valor(valor_file, hierarchy):
    tag = "gFC1sigma"
    if hierarchy == "normal":
        mh = "_nh"
    elif hierarchy == "inverted":
        mh = "_ih"
    else:
        print "Error: Unknown hierarchy \"" + hierarchy + "\""
        sys.exit()

    central = get_valor_graph(valor_file, tag + mh)
    err_m = get_valor_graph(valor_file, tag + mh + "_m1")
    err_p = get_valor_graph(valor_file, tag + mh + "_p1")
    return (central, err_m, err_p)

def get_2sigma_valor(valor_file, hierarchy):
    tag = "gFC2sigma"
    if hierarchy == "normal":
        mh = "_nh"
    elif hierarchy == "inverted":
        mh = "_ih"
    else:
        print "Error: Unknown hierarchy \"" + hierarchy + "\""
        sys.exit()

    central = get_valor_graph(valor_file, tag + mh)
    err_m = get_valor_graph(valor_file, tag + mh + "_m1")
    err_p = get_valor_graph(valor_file, tag + mh + "_p1")
    return (central, err_m, err_p)

def get_90cl_valor(valor_file, hierarchy):
    tag = "gFC90CL"
    if hierarchy == "normal":
        mh = "_nh"
    elif hierarchy == "inverted":
        mh = "_ih"
    else:
        print "Error: Unknown hierarchy \"" + hierarchy + "\""
        sys.exit()

    central = get_valor_graph(valor_file, tag + mh)
    err_m = get_valor_graph(valor_file, tag + mh + "_m1")
    err_p = get_valor_graph(valor_file, tag + mh + "_p1")
    return (central, err_m, err_p)

def compare_fc(valor_file, out_filename):
    canvas = ROOT.TCanvas("fc_nh", "fc_nh", 1024, 768)
    canvas.SetLeftMargin(0.125)
    canvas.SetBottomMargin(0.125)
    legend = ROOT.TLegend(0.7, 0.2, 0.85, 0.3)

    v1sigmas = get_1sigma_valor(valor_file, "normal")
    p1sigmas = get_1sigma_ptheta("normal")
    v90s = get_90cl_valor(valor_file, "normal")
    p90s = get_90cl_ptheta("normal")
    v2sigmas = get_2sigma_valor(valor_file, "normal")
    p2sigmas = get_2sigma_ptheta("normal")

    p2sigma_nh = p2sigmas[0]

    vmax = ROOT.TMath.MaxElement(v2sigmas[2].GetN(), v2sigmas[2].GetY())
    pmax = ROOT.TMath.MaxElement(p2sigmas[2].GetN(), p2sigmas[2].GetY())
    if vmax > pmax:
        gmax = vmax
    else:
        gmax = pmax
    gmax = ROOT.TMath.Ceil(gmax)
    opt = ""
    # 1 sigma
    for i, g in enumerate(v1sigmas):
        g.GetYaxis().SetRangeUser(0.0, gmax)
        g.SetLineColor(ROOT.kBlack)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.GetXaxis().SetTitle("#delta_{CP}")
        g.GetYaxis().SetTitle("#Delta#chi^2")
        g.SetTitle("FC NH critical values")
        g.Draw(opt)
        opt = "same"
    for i, g in enumerate(p1sigmas):
        g.SetLineColor(ROOT.kRed)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
    # 90CL
    for i, g in enumerate(v90s):
        g.GetYaxis().SetRangeUser(0.0, gmax)
        g.SetLineColor(ROOT.kMagenta)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
        opt = "same"
    for i, g in enumerate(p90s):
        g.SetLineColor(ROOT.kCyan)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
    # 2 sigma
    for i, g in enumerate(v2sigmas):
        g.GetYaxis().SetRangeUser(0.0, gmax)
        g.SetLineColor(ROOT.kBlue)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
        opt = "same"
    for i, g in enumerate(p2sigmas):
        g.SetLineColor(ROOT.kGreen)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)

    legend.AddEntry(v1sigmas[0], "VALOR 1#sigma", "l")
    legend.AddEntry(p1sigmas[0], "p-theta 1#sigma", "l")
    legend.AddEntry(v90s[0], "VALOR 90CL", "l")
    legend.AddEntry(p90s[0], "p-theta 90CL", "l")
    legend.AddEntry(v2sigmas[0], "VALOR 2#sigma", "l")
    legend.AddEntry(p2sigmas[0], "p-theta 2#sigma", "l")

    legend.Draw()

    canvas.Print(out_filename + "(")

    canvas = ROOT.TCanvas("fc_ih", "fc_ih", 1024, 768)
    canvas.SetLeftMargin(0.125)
    canvas.SetBottomMargin(0.125)
    legend = ROOT.TLegend(0.7, 0.2, 0.85, 0.3)

    v1sigmas = get_1sigma_valor(valor_file, "inverted")
    p1sigmas = get_1sigma_ptheta("inverted")
    v90s = get_90cl_valor(valor_file, "inverted")
    p90s = get_90cl_ptheta("inverted")
    v2sigmas = get_2sigma_valor(valor_file, "inverted")
    p2sigmas = get_2sigma_ptheta("inverted")

    p2sigma_ih = p2sigmas[0]

    write_ptheta_graph(p2sigma_nh, p2sigma_ih)

    vmax = ROOT.TMath.MaxElement(v2sigmas[2].GetN(), v2sigmas[2].GetY())
    pmax = ROOT.TMath.MaxElement(p2sigmas[2].GetN(), p2sigmas[2].GetY())
    if vmax > pmax:
        gmax = vmax
    else:
        gmax = pmax
    gmax = ROOT.TMath.Ceil(gmax)
    opt = ""
    # 1 sigma
    for i, g in enumerate(v1sigmas):
        g.GetYaxis().SetRangeUser(0.0, gmax)
        g.SetLineColor(ROOT.kBlack)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.GetXaxis().SetTitle("#delta_{CP}")
        g.GetYaxis().SetTitle("#Delta#chi^2")
        g.SetTitle("FC IH critical values")
        g.Draw(opt)
        opt = "same"
    for i, g in enumerate(p1sigmas):
        g.SetLineColor(ROOT.kRed)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
    # 90CL
    for i, g in enumerate(v90s):
        g.GetYaxis().SetRangeUser(0.0, gmax)
        g.SetLineColor(ROOT.kMagenta)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
        opt = "same"
    for i, g in enumerate(p90s):
        g.SetLineColor(ROOT.kCyan)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
    # 2 sigma
    for i, g in enumerate(v2sigmas):
        g.GetYaxis().SetRangeUser(0.0, gmax)
        g.SetLineColor(ROOT.kBlue)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
        opt = "same"
    for i, g in enumerate(p2sigmas):
        g.SetLineColor(ROOT.kGreen)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)

    legend.AddEntry(v1sigmas[0], "VALOR 1#sigma", "l")
    legend.AddEntry(p1sigmas[0], "p-theta 1#sigma", "l")
    legend.AddEntry(v90s[0], "VALOR 90CL", "l")
    legend.AddEntry(p90s[0], "p-theta 90CL", "l")
    legend.AddEntry(v2sigmas[0], "VALOR 2#sigma", "l")
    legend.AddEntry(p2sigmas[0], "p-theta 2#sigma", "l")

    legend.Draw()

    canvas.Print(out_filename + ")")


def main(in_filename, out_filename):
    valor_file = ROOT.TFile(in_filename, "READ")
    compare_fc(valor_file, out_filename)
    valor_file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Insert the new FSI parameters into the unprocessed BANFF matrix")
    parser.add_argument("--input", type=str, help="The VALOR FC graph", required=True)
    parser.add_argument("--output", type=str, help="Name of the output PDF file", required=True)
    args = parser.parse_args()
    main(args.input, args.output)
