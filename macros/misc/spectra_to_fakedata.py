#!/usr/bin/env python

import argparse
import ROOT

ROOT.gROOT.SetBatch(1)

def make_erec_spectrum(erec_theta):
    Nerec = erec_theta.GetNbinsY()
    erec_axis = erec_theta.GetYaxis()
    erec_low = erec_axis.GetBinLowEdge(1)
    erec_high = erec_axis.GetBinUpEdge(Nerec)
    Ntheta = erec_theta.GetNbinsX()
    theta_axis = erec_theta.GetXaxis()
    theta_low = theta_axis.GetBinLowEdge(1)
    theta_high = theta_axis.GetBinUpEdge(Ntheta)
    erec_name = erec_theta.GetName() + "_erec"
    erec = ROOT.TH1D(erec_name, erec_name, Nerec, erec_low, erec_high)
    for e in range(1, Nerec + 1):
        for t in range(1, Ntheta + 1):
            erec.SetBinContent(e, erec.GetBinContent(e) + erec_theta.GetBinContent(t, e))
    erec.SetDirectory(0)
    return erec

def get_sample(in_file, in_obj_name, out_obj_name):
    in_histo = in_file.Get(in_obj_name)
    out_histo = in_histo.Clone(out_obj_name)
    out_histo.SetDirectory(0)

    return out_histo

def combine(tag, out_filename):
    histograms = dict()
    histograms_erec = dict()
    for s in [ "1Rmu", "1Re", "MultiRe", "RHC1Rmu", "RHC1Re" ]:
        name = tag + "_" + s + ".root"
        in_file = ROOT.TFile(name)
        histograms[s] = get_sample(in_file, "curr_ereco_Vs_q2_" + s, "mock_data_" + s + "_0")
        histograms_erec[s] = make_erec_spectrum(histograms[s])
        in_file.Close()
        print histograms[s].GetName(), histograms[s].Integral(), histograms_erec[s].Integral()

    out_file = ROOT.TFile(out_filename, "RECREATE")
    for s in [ "1Rmu", "1Re", "MultiRe", "RHC1Rmu", "RHC1Re" ]:
        histograms[s].Write()
        histograms_erec[s].Write()
    out_file.Close()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-t", "--tag", dest="tag", default="", help="Prefix for sample root file.")
    parser.add_argument("-o", "--output", dest="output", default="", help="Output Valor fake data file.")
    args = parser.parse_args()
    combine(args.tag, args.output)
  
