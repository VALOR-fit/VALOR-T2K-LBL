#!/usr/bin/python2

# Plots numu -> numu, nue, nutau oscillation probabilities as a function of L/E, at fixed oscillation parameter values
# Steve Dennis, September 2014
# Francis Bench, November 2019

import ROOT
import sys

R=ROOT

env_var = "VALOR_FMWK"
try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $VALOR_FMWK/root/valor_interpreter.C")

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()

# NuFit 2018 best-fit values
sinsq_theta23 = 0.528
sinsq_theta12 = 0.304
sinsq_theta13 = 0.0212
dm23sq = 2.509E-3
dm12sq = 7.42E-5
deltacp = -2.50

l_over_e_t2k = 490.  # km/GeV

opt.gMassHierarchyAssumed = conv.kNmhNormal
opt.gOscMixingParamChoice = conv.kNuMixSinsqTheta
opt.gOscOctantChoice      = conv.kNuOctantTwo
opt.gOscDmsqScheme        = conv.kNuDmsqAuto
opt.SetDefaultParameterValues()

nbins = 1000

baseline = 295.


def main():

    all_hists = dict()

    osc = valor.NuOsc3Flavour()
    osc.SetBaseline(baseline)

    canv = ROOT.TCanvas("canv","canv",800,600)

    for ilimit,max_ratio in enumerate([4000]):     #([4000,40000]):

        h_mu  = R.TH1D("osc_prob_"+"mu" +"_"+`max_ratio`,"osc_prob_"+"mu" +"_"+`max_ratio`,nbins,0,max_ratio)
        h_e   = R.TH1D("osc_prob_"+"e"  +"_"+`max_ratio`,"osc_prob_"+"e"  +"_"+`max_ratio`,nbins,0,max_ratio)
        h_tau = R.TH1D("osc_prob_"+"tau"+"_"+`max_ratio`,"osc_prob_"+"tau"+"_"+`max_ratio`,nbins,0,max_ratio)
        l_exp = R.TLine(l_over_e_t2k, 0.0, l_over_e_t2k, 1.0)
        l_exp.SetLineColor(1)
        l_exp.SetLineStyle(2)
        l_exp.SetLineWidth(2)

        hists = [h_mu,h_e,h_tau]

        for h in hists:
            all_hists[h.GetName()] = h
            h.GetXaxis().SetTitle("L / E (km GeV^{-1})")
            h.GetYaxis().SetTitle("P_{#nu_{#mu} #rightarrow #nu_{#alpha}}")
            h.SetMaximum(1)
            h.SetMinimum(0)
            h.GetXaxis().SetTitleSize(0.05)
            h.GetYaxis().SetTitleSize(0.05)
            h.SetLineWidth(2)

        h_mu .SetLineColor(ROOT.kRed+1)
        h_e  .SetLineColor(ROOT.kBlue+2)
        h_tau.SetLineColor(ROOT.kGreen+2)

        for ii in xrange(1,nbins+1):
            lOverE = h_mu.GetBinCenter(ii)
            energy = baseline / lOverE
            osc.SetOscParametersSinSqTAll(sinsq_theta12,sinsq_theta13,sinsq_theta23,
                                          dm12sq, dm23sq, deltacp)

            h_mu .Fill(lOverE,osc.CalcMatterOscProb(conv.kPdgNuMu, conv.kPdgNuMu , energy))
            h_e  .Fill(lOverE,osc.CalcMatterOscProb(conv.kPdgNuMu, conv.kPdgNuE  , energy))
            h_tau.Fill(lOverE,osc.CalcMatterOscProb(conv.kPdgNuMu, conv.kPdgNuTau, energy))

        h_mu .Draw("hist")
        h_e  .Draw("hist,same")
        h_tau.Draw("hist,same")
        l_exp.Draw("same")


    legend = R.TLegend(0.85,0.65,0.95,0.95)
    legend.AddEntry(h_mu , "#nu_{#mu}"           , "l")
    legend.AddEntry(h_e  , "#nu_{e}"             , "l")
    legend.AddEntry(h_tau, "#nu_{#tau}"          , "l")
    legend.AddEntry(l_exp, "#splitline{T2K}{L/E}", "l")
    legend.AddEntry(None , "", "")
    legend.SetFillStyle(1001)
    legend.SetTextFont(132)
    legend.SetTextSize(0.04)

    legend.Draw()

    canv.Print("oscillation_patterns.pdf")

if __name__=="__main__":
    main()

