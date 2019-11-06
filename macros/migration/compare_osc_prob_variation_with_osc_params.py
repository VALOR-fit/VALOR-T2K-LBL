#!/usr/bin/python

# Plots numu -> numu, nue, nutau oscillation probabilities as a function of osc param values, at fixed energy.
# Each standard 3-flavour osc param is varied in turn, with all others held to nominal values
# Francis Bench, November 2019

import ROOT
import sys

R=ROOT

env_var = "T2KOSCANA"
try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C")

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()

hist_prefix = "osc_prob"
params = ["sinsq_theta23", "sinsq_theta12", "sinsq_theta13",  "dm23sq", "dm12sq", "deltacp"]
flavours = ["mu", "e", "tau"]
colour_leg = dict()
colour_leg["mu"]  = 40
colour_leg["e"]   = 41
colour_leg["tau"] = 42

filename_sdk = "oscillation_probability_vs_mixing_params.root"
filename_leg = "oscillation_probability_vs_mixing_params_legacy.root"
save_name    = "oscillation_probability_vs_mixing_params_sdk_vs_legacy.pdf"


def main():

    file_sdk = ROOT.TFile(filename_sdk, "read")
    file_leg = ROOT.TFile(filename_leg, "read")

    for iparam, param in enumerate(params):
        canv = ROOT.TCanvas(param, param, 800, 800)
        hists_sdk = dict()
        hists_leg = dict()
        for f in flavours:
            hname = hist_prefix + "_" + f + "_vs_" + param
            hists_sdk[f] = file_sdk.Get(hname)
            hists_leg[f] = file_leg.Get(hname)

            if not hists_sdk[f]:
                print "ERROR: Could not get the SDK hist named", hname
            if not hists_leg[f]:
                print "ERROR: Could not get the legacy hist named", hname

            hists_leg[f].SetLineStyle(2)
            hists_leg[f].SetLineColor(colour_leg[f])

            hists_sdk[f].Draw("hist,same")
            hists_leg[f].Draw("hist,same")

        legend = R.TLegend(0.90,0.6,0.95,0.9)
        legend.AddEntry(hists_sdk["mu"] , "#nu_{#mu}" , "l")
        legend.AddEntry(hists_sdk["e"]  , "#nu_{e}"   , "l")
        legend.AddEntry(hists_sdk["tau"], "#nu_{#tau}", "l")
        legend.AddEntry(hists_leg["mu"] , "Legacy"    , "l")
        legend.SetFillStyle(0)
        legend.Draw()

        tag_canv = save_name
        if iparam == 0:
            tag_canv += "("
        elif iparam == len(params)-1:
            tag_canv += ")"
        canv.SaveAs(tag_canv)


if __name__=="__main__":
    main()

