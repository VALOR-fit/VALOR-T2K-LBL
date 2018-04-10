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
            1.38459, 0.894935, 0.452939, 1.10762, 1.55957, 1.62541, 1.42472, 1.46937, 1.38459
        ])
        err_p = np.array([
            1.4127, 0.908336, 0.466189, 1.122, 1.59256, 1.64983, 1.45046, 1.49598, 1.4127
        ])
        err_m = np.array([
            1.35677, 0.88047, 0.43791, 1.09131, 1.53699, 1.60075, 1.40786, 1.45178, 1.35677
        ])
    else:
        central = np.array([
            1.73789, 1.60937, 1.46242, 1.49188, 1.42776, 1.16335, 0.925116, 1.33708, 1.73789
        ])
        err_p = np.array([
            1.75929, 1.62937, 1.47615, 1.51508, 1.44882, 1.18413, 0.948731, 1.35603, 1.75929
        ])
        err_m = np.array([
            1.71822, 1.58937, 1.44768, 1.47042, 1.40686, 1.1425, 0.907833, 1.31507, 1.71822
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
            3.71947, 2.84271, 2.73064, 2.92548, 3.87466, 4.02331, 3.7634, 3.94602, 3.71947
        ])
        err_p = np.array([
            3.76586, 2.91378, 2.79506, 2.9797, 3.91693, 4.08164, 3.838, 4.00413, 3.76586
        ])
        err_m = np.array([
            3.66347, 2.77395, 2.65064, 2.86766, 3.82874, 3.98035, 3.70246, 3.88164, 3.66347
        ])
    else:
        central = np.array([
            4.52787, 4.00032, 3.51173, 3.89359, 3.98404, 3.81199, 4.0194, 4.13833, 4.52787
        ])
        err_p = np.array([
            4.60943, 4.09032, 3.55508, 3.95256, 4.02719, 3.87523, 4.12684, 4.22923, 4.60943
        ])
        err_m = np.array([
            4.47151, 3.92032, 3.47272, 3.86361, 3.93045, 3.71958, 3.93889, 4.10371, 4.47151
        ])
    gcentral = ROOT.TGraph(len(central), dcp, central)
    gerr_p = ROOT.TGraph(len(err_p), dcp, err_p)
    gerr_m = ROOT.TGraph(len(err_m), dcp, err_m)
    return (gcentral, gerr_m, gerr_p)

def get_90cl_ptheta(hierarchy):
    dcp = np.array([
        -3.14159265359,-2.35619449019,-1.57079632679,
        -0.78539816339,0.0,0.78539816339,
        1.57079632679,2.35619449019,3.14159265359
    ])
    if hierarchy == "normal":
        central = np.array([
            2.85555, 1.86254, 1.68463, 2.12269, 3.05463, 3.03986, 2.74489, 2.89147, 2.85555
        ])
        err_p = np.array([
            2.88614, 1.89501, 1.72503, 2.1706, 3.08942, 3.07271, 2.79802, 2.93397, 2.88614
        ])
        err_m = np.array([
            2.82137, 1.83111, 1.63454, 2.09169, 3.03082, 3.00367, 2.6958, 2.85256, 2.82137
        ])
    else:
        central = np.array([
            3.37474, 3.01003, 2.68589, 2.93446, 2.91553, 2.56831, 2.6889, 2.87143, 3.37474
        ])
        err_p = np.array([
            3.42149, 3.06003, 2.71874, 2.98288, 2.96902, 2.60559, 2.72538, 2.92337, 3.42149
        ])
        err_m = np.array([
            3.33278, 2.961003, 2.65874, 2.89585, 2.87833, 2.52079, 2.64378, 2.8307, 3.33278
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
