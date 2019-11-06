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

param_ranges         = dict()
param_nominal_values = dict()

param_ranges["sinsq_theta23"] = [0.3, 0.7]
param_ranges["sinsq_theta12"] = [0.25, 0.35]
param_ranges["sinsq_theta13"] = [0.018, 0.025]
param_ranges["dm23sq"]        = [2.2E-3, 2.8E-3]
param_ranges["dm12sq"]        = [6.8E-5, 8.2E-5]
param_ranges["deltacp"]       = [-3.141592653589793, 3.141592653589793]

param_nominal_values["sinsq_theta23"] = 0.528
param_nominal_values["sinsq_theta12"] = 0.304
param_nominal_values["sinsq_theta13"] = 0.0212
param_nominal_values["dm23sq"]        = 2.509E-3
param_nominal_values["dm12sq"]        = 7.53E-5
param_nominal_values["deltacp"]       = -1.601

opt.gMassHierarchyAssumed = conv.kNmhNormal
opt.gOscMixingParamChoice = conv.kNuMixSinsqTheta
opt.gOscOctantChoice      = conv.kNuOctantTwo
opt.gOscDmsqScheme        = conv.kNuDmsqAuto
opt.SetDefaultParameterValues()

nbins = 1000

baseline = 295. # In km
energy   = 0.6  # In Gev

save_name = "oscillation_probability_vs_mixing_params_legacy"
file_format = "pdf"

title = "P_{#nu_{#mu} #rightarrow #nu_{#alpha}} at L = " + str(baseline) + " km, E = " + str(energy) + " Gev"


def main():

    all_hists = dict()

    osc = valor.NuOsc3Flavour()
    osc.SetBaseline(baseline)

    tag_root = save_name + "." + "root"
    fout = ROOT.TFile(tag_root, "recreate")

    for iparam, param_name in enumerate(param_ranges.keys()):
        print "Varying param", param_name

        canv = ROOT.TCanvas(param_name, param_name, 800, 800)

        low  = param_ranges[param_name][0]
        high = param_ranges[param_name][1]
        step = (high - low) / float(nbins)

        h_mu  = R.TH1D("osc_prob_"+"mu" +"_vs_"+param_name, "osc_prob_"+"mu" +"_vs_"+param_name, nbins, low, high)
        h_e   = R.TH1D("osc_prob_"+"e"  +"_vs_"+param_name, "osc_prob_"+"e"  +"_vs_"+param_name, nbins, low, high)
        h_tau = R.TH1D("osc_prob_"+"tau"+"_vs_"+param_name, "osc_prob_"+"tau"+"_vs_"+param_name, nbins, low, high)

        hists = [h_mu,h_e,h_tau]

        for h in hists:
            all_hists[h.GetName()] = h
            h.GetXaxis().SetTitle(param_name)
            h.GetYaxis().SetTitle(title)
            h.SetMaximum(1)
            h.SetMinimum(0)
            h.GetXaxis().SetTitleSize(0.045)
            h.GetYaxis().SetTitleSize(0.045)
            h.GetXaxis().SetLabelSize(0.03)
            h.GetYaxis().SetLabelSize(0.03)
            h.SetLineWidth(2)

        h_mu .SetLineColor(ROOT.kRed+1)
        h_e  .SetLineColor(ROOT.kBlue+2)
        h_tau.SetLineColor(ROOT.kGreen+2)

        for ii in xrange(0,nbins):
            val = low + step * (ii + 0.5)

            vals = dict()
            for pn in param_ranges.keys():
                if pn == param_name:
                    vals[pn] = val
                else:
                    vals[pn] = param_nominal_values[pn]

            osc.SetOscParametersSinSqTAll(vals["sinsq_theta12"], vals["sinsq_theta13"], vals["sinsq_theta23"],
                                          vals["dm12sq"], vals["dm23sq"], vals["deltacp"])

            h_mu .Fill(val, osc.CalcMatterOscProb(conv.kPdgNuMu, conv.kPdgNuMu , energy))
            h_e  .Fill(val, osc.CalcMatterOscProb(conv.kPdgNuMu, conv.kPdgNuE  , energy))
            h_tau.Fill(val, osc.CalcMatterOscProb(conv.kPdgNuMu, conv.kPdgNuTau, energy))
            #print "For osc params:", vals, " P(mu->mu) =", osc.CalcMatterOscProb(conv.kPdgNuMu, conv.kPdgNuMu , energy)

        h_mu .Draw("hist")
        h_e  .Draw("hist,same")
        h_tau.Draw("hist,same")

        legend = R.TLegend(0.90,0.6,0.95,0.9)
        legend.AddEntry(h_mu , "#nu_{#mu}" , "l")
        legend.AddEntry(h_e  , "#nu_{e}"   , "l")
        legend.AddEntry(h_tau, "#nu_{#tau}", "l")
        legend.SetFillStyle(0)
        legend.Draw()

        tag_canv = save_name + "." + file_format
        if iparam == 0:
            tag_canv += "("
        elif iparam == len(param_ranges.keys())-1:
            tag_canv += ")"
        canv.SaveAs(tag_canv)

        fout.cd()
        h_mu .Write()
        h_e  .Write()
        h_tau.Write()

    fout.Close()


if __name__=="__main__":
    main()

