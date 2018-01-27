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
    return erec

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

    return new_histo

def convert(in_filename, out_filename):
    in_file = ROOT.TFile(in_filename)
    fhc1re = convert_e_sample(in_file, "hNueEthetaPDF", "mock_data_1Re_0")
    rhc1re = convert_e_sample(in_file, "hNuebarEthetaPDF", "mock_data_RHC1Re_0")
    fhc1re_cc1pi = convert_e_sample(in_file, "hNueDEthetaPDF", "mock_data_MultiRe_0")
    fhc1rmu = convert_mu_sample(in_file, "hNumuPDF", "mock_data_1Rmu_0")
    rhc1rmu = convert_mu_sample(in_file, "hNumubarPDF", "mock_data_RHC1Rmu_0")

    fhc1re_erec = make_erec_spectrum(fhc1re)
    rhc1re_erec = make_erec_spectrum(rhc1re)
    fhc1re_cc1pi_erec = make_erec_spectrum(fhc1re_cc1pi)
    fhc1rmu_erec = make_erec_spectrum(fhc1rmu)
    rhc1rmu_erec = make_erec_spectrum(rhc1rmu)

    print fhc1re.GetName(), fhc1re.Integral(), fhc1re_erec.Integral()
    print rhc1re.GetName(), rhc1re.Integral(), rhc1re_erec.Integral()
    print fhc1re_cc1pi.GetName(), fhc1re_cc1pi.Integral(), fhc1re_cc1pi_erec.Integral()
    print fhc1rmu.GetName(), fhc1rmu.Integral(), fhc1rmu_erec.Integral()
    print rhc1rmu.GetName(), rhc1rmu.Integral(), rhc1rmu_erec.Integral()

    out_file = ROOT.TFile(out_filename, "RECREATE")
    out_file.cd()
    fhc1re.Write()
    rhc1re.Write()
    fhc1rmu.Write()
    rhc1rmu.Write()
    fhc1re_cc1pi.Write()
    fhc1re_erec.Write()
    rhc1re_erec.Write()
    fhc1rmu_erec.Write()
    rhc1rmu_erec.Write()
    fhc1re_cc1pi_erec.Write()

    out_file.Close()
    in_file.Close()

def get_rates(filename):
    in_file = ROOT.TFile(filename)
    fhc1re = in_file.Get("hNuePDF")
    rhc1re = in_file.Get("hNuebarPDF")
    fhc1rmu = in_file.Get("hNumuPDF")
    rhc1rmu = in_file.Get("hNumubarPDF")
    fhc1re_cc1pi = in_file.Get("hNueDPDF")
    print fhc1re.GetName(), fhc1re.Integral()
    print rhc1re.GetName(), rhc1re.Integral()
    print fhc1re_cc1pi.GetName(), fhc1re_cc1pi.Integral()
    print fhc1rmu.GetName(), fhc1rmu.Integral()
    print rhc1rmu.GetName(), rhc1rmu.Integral()
    in_file.Close()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--filename", dest="in_filename", default="", help="Input p-theta fake data file.")
    parser.add_argument("-o", "--output-filename", dest="out_filename", default="", help="Output Valor fake data file.")
    args = parser.parse_args()
    convert(args.in_filename, args.out_filename)
    get_rates(args.in_filename)
  
