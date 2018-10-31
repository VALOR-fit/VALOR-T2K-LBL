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
            1.42841,   0.971327,     0.517528,     1.21209,     1.63962,     1.69637,     1.49823,     1.4967,     1.42841
        ])
        err_p = np.array([
            1.44494,    0.983721,   0.528202,   1.22372,    1.65732,    1.71246,    1.51303,    1.51073,    1.44494
        ])
        err_m = np.array([
            1.41031,    0.95818,    0.506787,   1.1973,     1.6214,     1.6795,     1.48509,    1.47994,    1.41031
        ])
    else:
        central = np.array([
            1.83198,    1.78041,    1.59909,    1.56084,    1.49057,    1.14834,    0.801538,   1.38609,    1.83198
        ])
        err_p = np.array([
            1.85178,    1.79441,    1.61624,    1.5761,     1.51044,    1.1589,     0.816973,   1.39869,    1.85178
        ])
        err_m = np.array([
            1.81327,    1.76345,    1.58534,    1.54397,    1.47345,    1.13291,    0.787904,   1.37097,    1.81327
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
            3.94117,     2.98058,     3.13813,     3.23689,    4.15363,     4.29318,     3.78187,     4.06287,     3.94117
        ])
        err_p = np.array([
            3.9734,     3.01493,    3.19812,    3.27908,    4.18265,    4.34552,    3.82385,    4.11013,    3.9734
        ])
        err_m = np.array([
            3.90175,    2.94359,    3.06977,    3.19851,    4.11212,    4.24951,    3.73411,    4.01657,    3.90175
        ])
    else:
        central = np.array([
            4.53456,    4.37234,    3.89413,    4.03318,    4.15353,    3.62808,    3.72282,    3.86741,    4.53456            
        ])
        err_p = np.array([
            4.56928,    4.4298,     3.93771,    4.08019,    4.21451,    3.67862,    3.77717,    3.92508,    4.56928
        ])
        err_m = np.array([
            4.49626,    4.32694,    3.85826,    3.99947,    4.11508,    3.58497,    3.64838,    3.82819,    4.49626
        ])
    gcentral = ROOT.TGraph(len(central), dcp, central)
    gerr_p = ROOT.TGraph(len(err_p), dcp, err_p)
    gerr_m = ROOT.TGraph(len(err_m), dcp, err_m)
    return (gcentral, gerr_m, gerr_p)

def get_3sigma_ptheta(hierarchy):
    dcp = np.array([
        -3.14159265359,-2.35619449019,-1.57079632679,
        -0.78539816339,0.0,0.78539816339,
        1.57079632679,2.35619449019,3.14159265359
    ])
    if hierarchy == "normal":
        central = np.array([
            7.14189,    6.9917,     7.25759,    7.35769,    7.49438,    8.51605,    8.1012,     7.94791,    7.14189
        ])
        err_p = np.array([
            7.28868,    7.10152,    7.47141 ,   7.62189,    7.57498,    8.71947,    8.39234,    8.16214,    7.28868
        ])
        err_m = np.array([
            7.04791,    6.88437,    7.04536,    7.21949,    7.3996,     8.24896,    7.89309,    7.76813,    7.04791
        ])
    else:
        central = np.array([
            8.22741,    7.59662,    7.96901,    7.6689,     8.19489,    7.85326,    8.45545,    8.12065,    8.22741
        ])
        err_p = np.array([
            8.34986,    7.71256,    8.10219,    7.84004,    8.31989,    7.98715,    8.72958,    8.39883,    8.34986
        ])
        err_m = np.array([
            8.06399,    7.46685,    7.80227,    7.38685,    7.99173,    7.6282  ,   8.30157,    7.90948,    8.06399
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
            2.92791,     2.06695,     1.94758,     2.33131,     3.15869,     3.21636 ,     2.81765,     2.9884,     2.92791
        ])
        err_p = np.array([
            2.96108,    2.0882,     1.98061,    2.35553,    3.18286,    3.2467,     2.84582,    3.01756,    2.96108
        ])
        err_m = np.array([
            2.89644,    2.03169,    1.91192,    2.30207,    3.13022,    3.18489,    2.7854,     2.95169,    2.89644
        ])
    else:
        central = np.array([
            3.43813,    3.3644,     2.94837,    3.07676 ,   3.07233,    2.49237,    2.47119,    2.76296,    3.43813
        ])
        err_p = np.array([
            3.464,  3.39075,    2.98027,    3.10756,    3.10165,    2.51808,    2.49906,    2.78725,    3.464
        ])
        err_m = np.array([
            3.40271,    3.32448,    2.91927,    3.03507,    3.03662,    2.45925,    2.4254,     2.72632,    3.40271
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

def get_3sigma_valor(valor_file, hierarchy):
    tag = "gFC3sigma"
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
    canvas.SetRightMargin(0.2)

    v1sigmas = get_1sigma_valor(valor_file, "normal")
    p1sigmas = get_1sigma_ptheta("normal")
    v90s = get_90cl_valor(valor_file, "normal")
    p90s = get_90cl_ptheta("normal")
    v2sigmas = get_2sigma_valor(valor_file, "normal")
    p2sigmas = get_2sigma_ptheta("normal")
    v3sigmas = get_3sigma_valor(valor_file, "normal")
    p3sigmas = get_3sigma_ptheta("normal")

    p2sigma_nh = p2sigmas[0]

    vmax = ROOT.TMath.MaxElement(v3sigmas[2].GetN(), v3sigmas[2].GetY())
    pmax = ROOT.TMath.MaxElement(p3sigmas[2].GetN(), p3sigmas[2].GetY())
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
        g.GetYaxis().SetTitle("#Delta#chi^{2}")
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
#    for i, g in enumerate(v90s):
#        g.GetYaxis().SetRangeUser(0.0, gmax)
#        g.SetLineColor(ROOT.kMagenta)
#        if i > 0:
#            g.SetLineStyle(2)
#        else:
#            g.SetLineStyle(1)
#        g.Draw(opt)
#        opt = "same"
#    for i, g in enumerate(p90s):
#        g.SetLineColor(ROOT.kCyan)
#        if i > 0:
#            g.SetLineStyle(2)
#        else:
#            g.SetLineStyle(1)
#        g.Draw(opt)
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
    # 3 sigma
    for i, g in enumerate(v3sigmas):
        g.GetYaxis().SetRangeUser(0.0, gmax)
        g.SetLineColor(ROOT.kMagenta)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
        opt = "same"
    for i, g in enumerate(p3sigmas):
        g.SetLineColor(ROOT.kCyan)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)


    textsize = 0.04
    lymin = 0.65
    legend = ROOT.TLegend(0.755, lymin, 0.95, lymin + textsize * 6)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(textsize)
    legend.AddEntry(v1sigmas[0], "VALOR 1#sigma", "l")
    legend.AddEntry(p1sigmas[0], "p-theta 1#sigma", "l")
#    legend.AddEntry(v90s[0], "VALOR 90% CL", "l")
#    legend.AddEntry(p90s[0], "p-theta 90% CL", "l")
    legend.AddEntry(v2sigmas[0], "VALOR 2#sigma", "l")
    legend.AddEntry(p2sigmas[0], "p-theta 2#sigma", "l")
    legend.AddEntry(v3sigmas[0], "VALOR 3#sigma", "l")
    legend.AddEntry(p3sigmas[0], "p-theta 3#sigma", "l")

    legend.Draw()

    canvas.Print(out_filename + "(")

    canvas = ROOT.TCanvas("fc_ih", "fc_ih", 1024, 768)
    canvas.SetLeftMargin(0.125)
    canvas.SetBottomMargin(0.125)
    canvas.SetRightMargin(0.2)

    v1sigmas = get_1sigma_valor(valor_file, "inverted")
    p1sigmas = get_1sigma_ptheta("inverted")
    v90s = get_90cl_valor(valor_file, "inverted")
    p90s = get_90cl_ptheta("inverted")
    v2sigmas = get_2sigma_valor(valor_file, "inverted")
    p2sigmas = get_2sigma_ptheta("inverted")
    v3sigmas = get_3sigma_valor(valor_file, "inverted")
    p3sigmas = get_3sigma_ptheta("inverted")

    p2sigma_ih = p2sigmas[0]

    write_ptheta_graph(p2sigma_nh, p2sigma_ih)

    vmax = ROOT.TMath.MaxElement(v3sigmas[2].GetN(), v3sigmas[2].GetY())
    pmax = ROOT.TMath.MaxElement(p3sigmas[2].GetN(), p3sigmas[2].GetY())
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
#    for i, g in enumerate(v90s):
#        g.GetYaxis().SetRangeUser(0.0, gmax)
#        g.SetLineColor(ROOT.kMagenta)
#        if i > 0:
#            g.SetLineStyle(2)
#        else:
#            g.SetLineStyle(1)
#        g.Draw(opt)
#        opt = "same"
#    for i, g in enumerate(p90s):
#        g.SetLineColor(ROOT.kCyan)
#        if i > 0:
#            g.SetLineStyle(2)
#        else:
#            g.SetLineStyle(1)
#        g.Draw(opt)
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
    # 3 sigma
    for i, g in enumerate(v3sigmas):
        g.GetYaxis().SetRangeUser(0.0, gmax)
        g.SetLineColor(ROOT.kMagenta)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
        opt = "same"
    for i, g in enumerate(p3sigmas):
        g.SetLineColor(ROOT.kCyan)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)

    textsize = 0.04
    lymin = 0.65
    legend = ROOT.TLegend(0.755, lymin, 0.95, lymin + textsize * 6)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(textsize)
    legend.AddEntry(v1sigmas[0], "VALOR 1#sigma", "l")
    legend.AddEntry(p1sigmas[0], "p-theta 1#sigma", "l")
#    legend.AddEntry(v90s[0], "VALOR 90% CL", "l")
#    legend.AddEntry(p90s[0], "p-theta 90% CL", "l")
    legend.AddEntry(v2sigmas[0], "VALOR 2#sigma", "l")
    legend.AddEntry(p2sigmas[0], "p-theta 2#sigma", "l")
    legend.AddEntry(v2sigmas[0], "VALOR 3#sigma", "l")
    legend.AddEntry(p2sigmas[0], "p-theta 3#sigma", "l")

    legend.Draw()

    canvas.Print(out_filename + ")")

def compare_valor_fc(valor_file1, valor_file2, out_filename):
    canvas = ROOT.TCanvas("fc_nh", "fc_nh", 1024, 768)
    canvas.SetLeftMargin(0.125)
    canvas.SetBottomMargin(0.125)
    canvas.SetRightMargin(0.2)

    v1sigmas = get_1sigma_valor(valor_file1, "normal")
    v1sigmas2 = get_1sigma_valor(valor_file2, "normal")
    v90s = get_90cl_valor(valor_file1, "normal")
    v90s2 = get_90cl_valor(valor_file2, "normal")
    v2sigmas = get_2sigma_valor(valor_file1, "normal")
    v2sigmas2 = get_2sigma_valor(valor_file2, "normal")

    v3sigmas = get_3sigma_valor(valor_file1, "normal")
    v3sigmas2 = get_3sigma_valor(valor_file2, "normal")

    vmax = ROOT.TMath.MaxElement(v3sigmas2[2].GetN(), v3sigmas2[2].GetY())
    gmax = ROOT.TMath.Ceil(vmax)
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
    for i, g in enumerate(v1sigmas2):
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
    for i, g in enumerate(v90s2):
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
    for i, g in enumerate(v2sigmas2):
        g.SetLineColor(ROOT.kGreen)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
    # 3 sigma
    for i, g in enumerate(v3sigmas):
        g.GetYaxis().SetRangeUser(0.0, gmax)
        g.SetLineColor(8)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
        opt = "same"
    for i, g in enumerate(v3sigmas2):
        g.SetLineColor(9)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)

    textsize = 0.04
    lymin = 0.55
    legend = ROOT.TLegend(0.755, lymin, 0.95, lymin + textsize * 8)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(textsize)
#    legend.AddEntry(v1sigmas[0], "1-9d 1#sigma", "l")
#    legend.AddEntry(v1sigmas2[0], "1-9c 1#sigma", "l")
#    legend.AddEntry(v90s[0], "1-9d 90%", "l")
#    legend.AddEntry(v90s2[0], "1-9c 90%", "l")
#    legend.AddEntry(v2sigmas[0], "1-9d 2#sigma", "l")
#    legend.AddEntry(v2sigmas2[0], "1-9c 2#sigma", "l")
#    legend.AddEntry(v3sigmas[0], "1-9d 3#sigma", "l")
#    legend.AddEntry(v3sigmas2[0], "1-9c 3#sigma", "l")
    legend.AddEntry(v1sigmas[0], "Set 1 1#sigma", "l")
    legend.AddEntry(v1sigmas2[0], "Set 2 1#sigma", "l")
    legend.AddEntry(v90s[0], "Set 1 90%", "l")
    legend.AddEntry(v90s2[0], "Set 2 90%", "l")
    legend.AddEntry(v2sigmas[0], "Set 1 2#sigma", "l")
    legend.AddEntry(v2sigmas2[0], "Set 2 2#sigma", "l")
    legend.AddEntry(v3sigmas[0], "Set 1 3#sigma", "l")
    legend.AddEntry(v3sigmas2[0], "Set 2 3#sigma", "l")


    legend.Draw()

    canvas.Print(out_filename + "(")

    # Inverted
    canvas = ROOT.TCanvas("fc_ih", "fc_ih", 1024, 768)
    canvas.SetLeftMargin(0.125)
    canvas.SetBottomMargin(0.125)
    canvas.SetRightMargin(0.2)

    v1sigmas = get_1sigma_valor(valor_file1, "inverted")
    v1sigmas2 = get_1sigma_valor(valor_file2, "inverted")
    v90s = get_90cl_valor(valor_file1, "inverted")
    v90s2 = get_90cl_valor(valor_file2, "inverted")
    v2sigmas = get_2sigma_valor(valor_file1, "inverted")
    v2sigmas2 = get_2sigma_valor(valor_file2, "inverted")

    v3sigmas = get_3sigma_valor(valor_file1, "inverted")
    v3sigmas2 = get_3sigma_valor(valor_file2, "inverted")

    vmax = ROOT.TMath.MaxElement(v3sigmas2[2].GetN(), v3sigmas2[2].GetY())
    gmax = ROOT.TMath.Ceil(vmax)
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
    for i, g in enumerate(v1sigmas2):
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
    for i, g in enumerate(v90s2):
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
    for i, g in enumerate(v2sigmas2):
        g.SetLineColor(ROOT.kGreen)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
    # 3 sigma
    for i, g in enumerate(v3sigmas):
        g.GetYaxis().SetRangeUser(0.0, gmax)
        g.SetLineColor(8)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)
        opt = "same"
    for i, g in enumerate(v3sigmas2):
        g.SetLineColor(9)
        if i > 0:
            g.SetLineStyle(2)
        else:
            g.SetLineStyle(1)
        g.Draw(opt)

    textsize = 0.04
    lymin = 0.55
    legend = ROOT.TLegend(0.755, lymin, 0.95, lymin + textsize * 8)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(textsize)
    legend.AddEntry(v1sigmas[0], "1-9d 1#sigma", "l")
    legend.AddEntry(v1sigmas2[0], "1-9c 1#sigma", "l")
    legend.AddEntry(v90s[0], "1-9d 90%", "l")
    legend.AddEntry(v90s2[0], "1-9c 90%", "l")
    legend.AddEntry(v2sigmas[0], "1-9d 2#sigma", "l")
    legend.AddEntry(v2sigmas2[0], "1-9c 2#sigma", "l")
    legend.AddEntry(v3sigmas[0], "1-9d 3#sigma", "l")
    legend.AddEntry(v3sigmas2[0], "1-9c 3#sigma", "l")

    legend.Draw()

    canvas.Print(out_filename + ")")


def main(in_filename, out_filename):
    valor_file = ROOT.TFile(in_filename, "READ")
    compare_fc(valor_file, out_filename)
    valor_file.Close()

def valor_only(in_filename, in_filename2, out_filename):
    valor_file1 = ROOT.TFile(in_filename, "READ")
    valor_file2 = ROOT.TFile(in_filename2, "READ")
    compare_valor_fc(valor_file1, valor_file2, out_filename)
    valor_file2.Close()
    valor_file1.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Insert the new FSI parameters into the unprocessed BANFF matrix")
    parser.add_argument("--input", type=str, help="The VALOR FC graph", required=True)
    parser.add_argument("--input2", type=str, help="A previous VALOR FC graph", required=False)
    parser.add_argument("--output", type=str, help="Name of the output PDF file", required=True)
    args = parser.parse_args()
    if args.input2:
        valor_only(args.input, args.input2, args.output)
    else:
        main(args.input, args.output)
