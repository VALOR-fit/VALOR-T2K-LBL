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
            1.34175,   0.927735,   0.486905,   1.15815, 1.79437, 1.75725, 1.62461, 1.63782, 1.58749
        ])
        err_p = np.array([
            1.36518,   0.944802,   0.505869,   1.17436, 1.81659, 1.78123, 1.64911, 1.66089, 1.60635
        ])
        err_m = np.array([
            1.31573,   0.908691,   0.467739,   1.14374, 1.76941, 1.73522, 1.60616, 1.61279, 1.5645
        ])
    else:
        central = np.array([
            1.68051, 1.70032, 1.56235, 1.61821, 1.59378, 1.14895, 0.619072, 1.25884, 1.68051
        ])
        err_p = np.array([
            1.70071,   1.72988, 1.5804, 1.63959, 1.61505, 1.1667, 0.637063, 1.28128, 1.70071
        ])
        err_m = np.array([
            1.66106, 1.67162, 1.54314, 1.59656, 1.57004, 1.13326, 0.603608, 1.2438, 1.66106
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
            3.74145,   2.83597,    2.94668,    3.09374, 4.12707, 4.22363, 3.78154, 4.23962, 4.15285
        ])
        err_p = np.array([
            3.79817,   2.89504,    3.01124,    3.14407, 4.19106, 4.28653, 3.86172, 4.29448, 4.19106
        ])
        err_m = np.array([
            3.70392,   2.78986,    2.89486,    3.03267, 4.10345, 4.17289, 3.72031, 4.17817, 4.10345
        ])
    else:
        central = np.array([
            4.24359, 4.38851, 4.06957, 4.24001, 4.18436, 3.34522, 3.20918, 3.61763, 4.24359
        ])
        err_p = np.array([
            4.29792, 4.48313, 4.14951, 4.29634, 4.24408, 3.42091, 3.30203, 3.70047, 4.29792
        ])
        err_m = np.array([
            4.18865, 4.34497, 4.00264, 4.17686, 4.13504, 3.28231, 3.1549, 3.55688, 4.18865
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
            2.81915,   1.95831,    1.83327,    2.20945, 3.36646, 3.22613, 2.89796, 3.15321, 3.13212
        ])
        err_p = np.array([
            2.85101,   1.99931,    1.87852,    2.23228, 3.42261, 3.25947, 2.92689, 3.19583, 3.17473
        ])
        err_m = np.array([
            2.78842,   1.92658,    1.78674,    2.18037, 3.32349, 3.19282, 2.86628, 3.11515, 3.09777
        ])
    else:
        central = np.array([
            3.24388, 3.42679, 3.01119, 3.14361, 3.14789, 2.32226, 2.17296, 2.53488, 3.24388
        ])
        err_p = np.array([
            3.28679, 3.46938, 3.0458, 3.18512, 3.18903, 2.35401, 2.22108, 2.57791, 3.28679
        ])
        err_m = np.array([
            3.20153, 3.37543, 2.96553, 3.09291, 3.10989, 2.28926, 2.12717, 2.4898, 3.20153
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
