#!/usr/bin/env python

import ROOT
import collections
from collections import namedtuple
import copy
from copy import deepcopy

ROOT.PyConfig.IgnoreCommandLineOptions = True

env_var = "T2KOSCANA"
try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");

ROOT.gROOT.SetBatch(1)

valor   = ROOT.valor
util    = valor.util
conv    = valor.conv
opt     = valor.Opt.Instance()
t2kopt  = valor.T2KOpt.Instance()
env     = valor.Env.Instance()
rnd     = valor.Rnd.Instance()
arridx  = util.ArrIdx.Instance()

title_size_x     = 0.07
title_size_y     = 0.08
title_offset_y   = 1.05
leg_text_size    = 0.06
label_text_size  = 0.055
label_offset_x   = 0.005
line_thickness   = 2
ymax = 0.13

data_nev    = 15.
data_dchisq = 3.81


def main():
    # TODO: make these input arguments
    outfile_tag       = "nuebarapp_distributions_paper"
    infile_rate_name  = "nuebarapp_graphs_nmargtoys=200k_nexp=20k_data_rate-only_MCsampling_VALOR_paper_rate-only_hists.root"
    infile_shape_name = "nuebarapp_graphs_nmargtoys=200k_nexp=20k_data_rate+shape_MCsampling_VALOR_paper_rate+shape_hists.root"
    hname_beta0       = "beta0"
    hname_beta1       = "beta1"

    # Get the input hists
    infile_rate  = ROOT.TFile(infile_rate_name, "read")
    infile_shape = ROOT.TFile(infile_shape_name, "read")
    hist_beta0_rate_in = infile_rate .Get(hname_beta0)
    hist_beta1_rate_in = infile_rate .Get(hname_beta1)
    hist_beta0_shape   = infile_shape.Get(hname_beta0)
    hist_beta1_shape   = infile_shape.Get(hname_beta1)

    # The rate-only hists start at -1 (so that y-axis tick marks didn't overlap with the hist for log-scale plots
    # in another macro), so rebin it to start at 0. The 1st bin is always empty, so this is fine.
    hist_beta0_rate = ROOT.TH1D(hist_beta0_rate_in.GetName(), hist_beta0_rate_in.GetTitle(),
                                hist_beta0_rate_in.GetNbinsX()-1, 0,
                                hist_beta0_rate_in.GetBinLowEdge(hist_beta0_rate_in.GetNbinsX()) + hist_beta0_rate_in.GetBinWidth(1))
    for ibin in range(2, hist_beta0_rate_in.GetNbinsX()+1):
        hist_beta0_rate.Fill(hist_beta0_rate_in.GetBinCenter(ibin), hist_beta0_rate_in.GetBinContent(ibin))

    hist_beta1_rate = ROOT.TH1D(hist_beta1_rate_in.GetName(), hist_beta1_rate_in.GetTitle(),
                                hist_beta1_rate_in.GetNbinsX()-1, 0,
                                hist_beta1_rate_in.GetBinLowEdge(hist_beta1_rate_in.GetNbinsX()) + hist_beta1_rate_in.GetBinWidth(1))
    for ibin in range(2, hist_beta1_rate_in.GetNbinsX()+1):
        hist_beta1_rate.Fill(hist_beta1_rate_in.GetBinCenter(ibin), hist_beta1_rate_in.GetBinContent(ibin))

    print "Hist integrals:\nrate beta0:", hist_beta0_rate.Integral(), "\nrate beta1:", hist_beta1_rate.Integral(), \
        "\nshape beta0:", hist_beta0_shape.Integral(), "\nshape beta1:", hist_beta1_shape.Integral()

    # Hist formatting
    hist_beta0_rate.SetMaximum(ymax)
    hist_beta0_rate.SetLineColor(3)
    hist_beta1_rate.SetLineColor(2)
    hist_beta0_rate.SetLineWidth(line_thickness)
    hist_beta1_rate.SetLineWidth(line_thickness)
    hist_beta0_rate.GetXaxis().SetTitle("#splitline{Number of observed #bar{#nu} mode}{single-ring e-like events}")
    hist_beta0_rate.GetYaxis().SetTitle("Fraction of pseudo-experiments")
    hist_beta0_rate.GetYaxis().SetTitleSize(title_size_y)
    hist_beta0_rate.GetXaxis().SetTitleSize(title_size_x)
    hist_beta0_rate.GetXaxis().CenterTitle(True)
    hist_beta0_rate.GetYaxis().CenterTitle(True)
    hist_beta0_rate.GetYaxis().SetTitleOffset(title_offset_y)
    hist_beta0_rate.GetXaxis().SetLabelSize(label_text_size)
    hist_beta0_rate.GetYaxis().SetLabelSize(label_text_size)
    hist_beta0_rate.GetXaxis().SetLabelOffset(label_offset_x)

    hist_beta0_shape.SetMaximum(ymax)
    hist_beta0_shape.SetLineColor(3)
    hist_beta1_shape.SetLineColor(2)
    hist_beta0_shape.SetLineWidth(line_thickness)
    hist_beta1_shape.SetLineWidth(line_thickness)
    hist_beta0_shape.GetXaxis().SetTitleSize(title_size_x)
    hist_beta0_shape.GetXaxis().CenterTitle(True)
    hist_beta0_shape.GetXaxis().SetLabelSize(label_text_size*1.1)
    hist_beta0_shape.GetXaxis().SetLabelOffset(label_offset_x)

    # Make data lines
    data_line_rate = ROOT.TLine(data_nev, 0, data_nev, ymax)
    data_line_rate.SetLineStyle(5)
    data_line_rate.SetLineWidth(line_thickness)

    data_line_shape = ROOT.TLine(data_dchisq, 0, data_dchisq, ymax)
    data_line_shape.SetLineStyle(5)
    data_line_shape.SetLineWidth(line_thickness)

    # Setup the canvas and draw everything
    canv = ROOT.TCanvas("nuebarapp_paper_dist", "nuebarapp_paper_dist", 800, 600)
    ROOT.gStyle.SetOptStat(0)
    canv.Divide(2,1,0,0)
    pad_left = canv.cd(1)
    hist_beta0_rate.Draw()
    hist_beta1_rate.Draw("same")
    data_line_rate.Draw("same")

    pad_right = canv.cd(2)
    hist_beta0_shape.Draw()
    hist_beta1_shape.Draw("same")
    data_line_shape.Draw("same")

    # Remove the right-hand plot y-axis
    hist_beta0_shape.GetYaxis().SetLabelOffset(999)
    hist_beta0_shape.GetYaxis().SetLabelSize(0)
    hist_beta0_shape.GetYaxis().SetTitleOffset(999)
    hist_beta0_shape.GetYaxis().SetTitleSize(0)

    # Set the margins so that the two plots are joined together
    pad_left .SetLeftMargin (0.17)
    pad_right.SetRightMargin(0.02)

    # Draw two blank pads over the overlapping axis labels (two as the labels are on both of the main pads)
    erase_xlow  = 0.95
    erase_xhigh = 1.05
    erase_ylow  = 0.10
    erase_yhigh = 0.17
    #print "erase box coords = (",erase_xlow,",",erase_ylow,"), (",erase_xhigh,",",erase_yhigh,")"
    pad_left.cd()
    pe_left = ROOT.TPad("pe_left", "pe_left", erase_xlow, erase_ylow, 1.0, erase_yhigh)
    pe_left.Draw()
    pe_left.SetFillColor(pad_left.GetFillColor())
    pe_left.SetBorderMode(0)
    pad_right.cd()
    pe_right = ROOT.TPad("pe_right", "pe_right", 0.0, erase_ylow, erase_xhigh - 1.0, erase_yhigh)
    pe_right.Draw()
    pe_right.SetFillColor(pad_left.GetFillColor())
    pe_right.SetBorderMode(0)

    # Make and draw the legend
    leg = ROOT.TLegend(0.50, 0.75, 1.0 - pad_right.GetRightMargin(), 1.0)
    leg.SetTextSize(leg_text_size)
    leg.AddEntry(hist_beta0_shape, "#splitline{Pseudo-data}{#beta = 0}", "l")
    leg.AddEntry(hist_beta1_shape, "#beta = 1"                         , "l")
    leg.AddEntry(data_line_shape , "T2K Data"                          , "l")
    leg.Draw("same")

    canv.Update()
    print "saving image to file", outfile_tag + ".pdf"
    canv.SaveAs(outfile_tag+".pdf" ,"recreate")
    canv.SaveAs(outfile_tag+".eps" ,"recreate")
    canv.SaveAs(outfile_tag+".png" ,"recreate")
    canv.SaveAs(outfile_tag+".root","recreate")


if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Plot the Delta chisq and event rate graphs and calculate p-values for the nuebar appearance analysis")
    #parser.add_option('--infilename_beta0', dest='infilename_beta0', default="", help='Root file output by marginalisation, with fake data toys (distribution toys) for beta = 0')


    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o, other = parser.parse_args()

    main()

    if other:
        print "Warning - ignored arguments: ",other



