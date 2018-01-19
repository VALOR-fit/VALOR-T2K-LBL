#!/usr/bin/env python

import argparse
import ROOT

ROOT.gROOT.SetBatch(1)

def convert_e_sample(in_file, in_obj_name, out_obj_name):
    histo = in_file.Get(in_obj_name)
    Nerec = histo.GetNbinsX()
    erec_axis = histo.GetXaxis()
    erec_low = erec_axis.GetBinLowEdge(1)
    erec_high = erec_axis.GetBinUpEdge(Nerec)
    Ntheta = histo.GetNbinsY()
    theta_axis = histo.GetYaxis()
    theta_low = theta_axis.GetBinLowEdge(1)
    theta_high = theta_axis.GetBinUpEdge(Ntheta)
    new_histo = ROOT.TH2D(out_obj_name, out_obj_name, Ntheta, theta_low, theta_high, Nerec, erec_low, erec_high)
    for e in range(1, Nerec + 1):
        for t in range(1, Ntheta + 1):
            new_histo.SetBinContent(t, e, histo.GetBinContent(e, t))
    erec_axis = new_histo.GetYaxis()
    theta_axis = new_histo.GetXaxis()

    return new_histo

def convert_mu_sample(in_file, in_obj_name, out_obj_name):
    histo = in_file.Get(in_obj_name)
    N = histo.GetNbinsX()
    x_axis = histo.GetXaxis()
    low = x_axis.GetBinLowEdge(1)
    high = x_axis.GetBinUpEdge(N)
    new_histo = ROOT.TH2D(out_obj_name, out_obj_name, 1, 0.0, 180.0, N, low, high)
    for e in range(1, N + 1):
        new_histo.SetBinContent(1, e, histo.GetBinContent(e))
    erec_axis = new_histo.GetYaxis()
    theta_axis = new_histo.GetXaxis()

    return new_histo

def convert(in_filename, out_filename):
    in_file = ROOT.TFile(in_filename)
    fhc1re = convert_e_sample(in_file, "hDistribNueEtheta", "mock_data_1Re_0")
    rhc1re = convert_e_sample(in_file, "hDistribNuebarEtheta", "mock_data_RHC1Re_0")
    fhc1re_cc1pi = convert_e_sample(in_file, "hDistribNueDEtheta", "mock_data_MultiRe_0")
    fhc1rmu = convert_mu_sample(in_file, "hDistribNumu", "mock_data_1Rmu_0")
    rhc1rmu = convert_mu_sample(in_file, "hDistribNumubar", "mock_data_RHC1Rmu_0")

    print fhc1re.GetName(), fhc1re.Integral()
    print rhc1re.GetName(), rhc1re.Integral()
    print fhc1re_cc1pi.GetName(), fhc1re_cc1pi.Integral()
    print fhc1rmu.GetName(), fhc1rmu.Integral()
    print rhc1rmu.GetName(), rhc1rmu.Integral()

    out_file = ROOT.TFile(out_filename, "RECREATE")
    out_file.cd()
    fhc1re.Write()
    rhc1re.Write()
    fhc1rmu.Write()
    rhc1rmu.Write()
    fhc1re_cc1pi.Write()
    out_file.Close()
    in_file.Close()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--filename", dest="in_filename", default="", help="Input p-theta fake data file.")
    parser.add_argument("-o", "--output-filename", dest="out_filename", default="", help="Output Valor fake data file.")
    args = parser.parse_args()
    convert(args.in_filename, args.out_filename)
  
