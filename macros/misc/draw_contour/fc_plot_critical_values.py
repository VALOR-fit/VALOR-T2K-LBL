#!/usr/bin/env python

import ROOT
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

def plot_critical_values(root_file, output_file):
    f = ROOT.TFile(root_file, "READ")
    names = ["gFC1sigma", "gFC90CL", "gFC2sigma", "gFC3sigma"]
    orderings = ["_nh", "_ih"]
    pt = ["_p1", "_m1"]
    colours = [ROOT.kOrange, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]
    labels = ["1#sigma", "90%", "2#sigma", "3#sigma"]
    line_styles = [1, 2]
   
    canvas = ROOT.TCanvas("canvas","canvas")
    canvas.DrawFrame(-3.6, 0.0, 3.6, 5.0)
    legend = ROOT.TLegend(0.87,0.5,1.00,0.85)
    legend.SetTextFont(132)
    legend.SetTextSize(0.05)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    g_shade = dict()

    for n_idx, name in enumerate(names):
        for m_idx, mo in enumerate(orderings):
            graph = f.Get(name + mo)
            graph_p1 = f.Get(name + mo + pt[0])
            graph_m1 = f.Get(name + mo + pt[1])

            n = 9
            xx = np.zeros(n)
            yy = np.zeros(n)
            yy_min = np.zeros(n)
            yy_max = np.zeros(n)
            for i in range(0, n):
                x = ROOT.Double(0)
                y = ROOT.Double(0)
                graph.GetPoint(i, x, y)
                xx[i] = x
                yy[i] = y
                graph_p1.GetPoint(i, x, y)
                yy_max[i] = y
                graph_m1.GetPoint(i, x, y)
                yy_min[i] = y
            full_name = name + mo
            g_shade[full_name] = ROOT.TGraph(2 * n)
            for i in range(0, n):
                g_shade[full_name].SetPoint(i, xx[i], yy_min[i])
                g_shade[full_name].SetPoint(n + i, xx[n - i - 1], yy_max[n - i - 1])

            g_shade[full_name].SetFillColorAlpha(colours[n_idx], 0.5)
            g_shade[full_name].SetFillStyle(1000)
            g_shade[full_name].SetTitle("")
            g_shade[full_name].GetXaxis().SetTitle("#delta_{CP} (Radians)")
            g_shade[full_name].GetYaxis().SetTitle("-2#Deltaln(L)_{crit}")
            g_shade[full_name].GetXaxis().SetTitleSize(0.055)
            g_shade[full_name].GetYaxis().SetTitleSize(0.055)
            g_shade[full_name].GetXaxis().SetLabelSize(0.04)
            g_shade[full_name].GetYaxis().SetLabelSize(0.04)
            g_shade[full_name].GetXaxis().CenterTitle()
            g_shade[full_name].GetYaxis().CenterTitle()
            if n_idx == 0 and m_idx == 0:
                g_shade[full_name].GetYaxis().SetRangeUser(0, 10)
                g_shade[full_name].Draw("af")
            else:
                g_shade[full_name].Draw("f same")
            graph.SetLineColor(colours[n_idx])
            graph.SetLineStyle(line_styles[m_idx])
            graph.Draw("lp same")
            if m_idx == 0:
                legend.AddEntry(graph, labels[n_idx], "l")
            if n_idx == len(names)-1 and m_idx == 0:
                legend.AddEntry(graph, "NO", "l")
            elif n_idx == len(names)-1 and m_idx == 1:
                legend.AddEntry(graph, "IO", "l")

    legend.Draw()

    title = ROOT.TLatex(0.0, 1.01 * 5.0, "T2K Run1-9 Preliminary")
    title.SetTextSize(0.05)
    #title.Draw()

    canvas.SaveAs(output_file + ".pdf")
    canvas.SaveAs(output_file + ".root")

    f.Close()

if __name__ == "__main__":
    plot_critical_values("FCgraphs.root", "FC_critical_shaded")

