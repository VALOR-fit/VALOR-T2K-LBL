#!/usr/bin/env python

import ROOT
import array
import argparse

fhc = ["run1", "run2", "run3b", "run3c", "run4", "run5a", "run5b", "run7a", "run7c", "run8", "run9a"]
rhc = ["run5c", "run6b", "run6c", "run6d", "run6e", "run7b", "run9b", "run9c"]
fhc_all = ["run1-9a", "run1-9a_v1"]
rhc_all = ["run5c-9c", "run5c-9d_v1"]
special = ["run5c-9c", "run9d"]

ROOT.gStyle.SetOptStat(0)

def get_flux_weight_histograms(flux_filename):
    flux_file = ROOT.TFile(flux_filename, "READ")
    nue = flux_file.Get("enu_sk_tuned13a_nue_ratio")
    nue.GetXaxis().SetRangeUser(0, 2)
    nue.SetDirectory(0)
    nueb = flux_file.Get("enu_sk_tuned13a_nueb_ratio")
    nueb.GetXaxis().SetRangeUser(0, 2)
    nueb.SetDirectory(0)
    numu = flux_file.Get("enu_sk_tuned13a_numu_ratio")
    numu.GetXaxis().SetRangeUser(0, 2)
    numu.SetDirectory(0)
    numub = flux_file.Get("enu_sk_tuned13a_numub_ratio")
    numub.GetXaxis().SetRangeUser(0, 2)
    numub.SetDirectory(0)
    flux_file.Close()
    return (nue, nueb, numu, numub)

def get_flux_ratios(histo_num, histo_denom, title):
    N = histo_num.GetNbinsX()
    ratio = histo_num.Clone(histo_num.GetName() + "_ratio")
    ratio.Reset()
    ratio.SetTitle(title)
    ratio.GetYaxis().SetRangeUser(0.9, 1.1)
    ratio.GetYaxis().SetTitle("Ratio of flux ratios")
    ratio.GetXaxis().SetTitle("E_{#nu} (GeV)")
    ratio.GetXaxis().SetTitleOffset(1.2)
    ratio.GetXaxis().SetTitleSize(0.04)
    ratio.GetYaxis().SetTitleOffset(1.3)
    ratio.GetYaxis().SetTitleSize(0.04)

    for i in range(1, N + 1):
        if histo_denom.GetBinContent(i) == 0:
            ratio.SetBinContent(i, histo_num.GetBinContent(i))
        else:
            ratio.SetBinContent(i, histo_num.GetBinContent(i) / histo_denom.GetBinContent(i))
    ratio.SetDirectory(0)
    ratio.SetMarkerStyle(22)
    ratio.SetMarkerColor(ROOT.kBlack)
    return ratio

def plot_flux_ratios(dir1, dir2, output):
    canvas = ROOT.TCanvas("canvas", "canvas", 1024, 768)
    for i, run in enumerate(fhc):
        filename1 = dir1 + "/" + run + "/sk_tuned13av2_13anom_" + run + "_numode_fine.root"
        nue1, nueb1, numu1, numub1 = get_flux_weight_histograms(filename1)
        filename2 = dir2 + "/" + run + "/sk_tuned13av3_13anom_" + run + "_numode_fine.root"
        nue2, nueb2, numu2, numub2 = get_flux_weight_histograms(filename2)
        nue_ratio = get_flux_ratios(nue2, nue1, run + " (FHC) #nu_{e}")
        nue_ratio.Draw("p")
        if i == 0:
            canvas.Print(output + "(")
        else:
            canvas.Print(output)
        nueb_ratio = get_flux_ratios(nueb2, nueb1, run + " (FHC) #bar{#nu}_{e}")
        nueb_ratio.Draw("p")
        canvas.Print(output)
        numu_ratio = get_flux_ratios(numu2, numu1, run + " (FHC) #nu_{#mu}")
        numu_ratio.Draw("p")
        canvas.Print(output)
        numub_ratio = get_flux_ratios(numub2, numub1, run + " (FHC) #bar{#nu}_{#mu}")
        numub_ratio.Draw("p")
    canvas.Print(output)
    for i, run in enumerate(rhc):
        filename1 = dir1 + "/" + run + "/sk_tuned13av2_13anom_" + run + "_antinumode_fine.root"
        nue1, nueb1, numu1, numub1 = get_flux_weight_histograms(filename1)
        filename2 = dir2 + "/" + run + "/sk_tuned13av3_13anom_" + run + "_antinumode_fine.root"
        nue2, nueb2, numu2, numub2 = get_flux_weight_histograms(filename2)
        nue_ratio = get_flux_ratios(nue2, nue1, run + " (RHC) #nu_{e}")
        nue_ratio.Draw("p")
        if i == 0:
            canvas.Print(output + "(")
        else:
            canvas.Print(output)
        nueb_ratio = get_flux_ratios(nueb2, nueb1, run + " (RHC) #bar{#nu}_{e}")
        nueb_ratio.Draw("p")
        canvas.Print(output)
        numu_ratio = get_flux_ratios(numu2, numu1, run + " (RHC) #nu_{#mu}")
        numu_ratio.Draw("p")
        canvas.Print(output)
        numub_ratio = get_flux_ratios(numub2, numub1, run + " (RHC) #bar{#nu}_{#mu}")
        numub_ratio.Draw("p")
    canvas.Print(output + ")")

def plot_special_flux_ratios(dir1, dir2, output):
    canvas = ROOT.TCanvas("canvas", "canvas", 1024, 768)
    filename1 = dir1 + "/" + special[0] + "/sk_tuned13av2_13anom_" + special[0] + "_antinumode.root"
    nue1, nueb1, numu1, numub1 = get_flux_weight_histograms(filename1)
    filename2 = dir2 + "/" + special[1] + "/sk_tuned13av3_13anom_" + special[1] + "_antinumode.root"
    nue2, nueb2, numu2, numub2 = get_flux_weight_histograms(filename2)
    nue1.Draw("p")
    nue2.Draw("p same")
    canvas.Print(output + "(")
    nueb1.Draw("p")
    nueb2.Draw("p same")
    canvas.Print(output)
    numu1.Draw("p")
    numu2.Draw("p same")
    canvas.Print(output)
    numub1.Draw("p")
    numub2.Draw("p same")
    canvas.Print(output + ")")

def plot_combined_flux_ratios(dir1, dir2, output):
    canvas = ROOT.TCanvas("canvas", "canvas", 1024, 768)
    filename1 = dir1 + "/" + fhc_all[0] + "/sk_tuned13av2_13anom_" + fhc_all[0] + "_numode_fine.root"
    nue1, nueb1, numu1, numub1 = get_flux_weight_histograms(filename1)
    filename2 = dir2 + "/" + fhc_all[1] + "/sk_tuned13av3_13anom_" + fhc_all[1][:-3] + "_numode_fine.root"
    nue2, nueb2, numu2, numub2 = get_flux_weight_histograms(filename2)
    nue_ratio = get_flux_ratios(nue2, nue1, fhc_all[1] + "/" + fhc_all[0] + " (FHC) #nu_{e}")
    nue_ratio.Draw("p")
    canvas.Print(output + "(")
    nueb_ratio = get_flux_ratios(nueb2, nueb1, fhc_all[1] + "/" + fhc_all[0] + " (FHC) #bar{#nu}_{e}")
    nueb_ratio.Draw("p")
    canvas.Print(output)
    numu_ratio = get_flux_ratios(numu2, numu1, fhc_all[1] + "/" + fhc_all[0] + " (FHC) #nu_{#mu}")
    numu_ratio.Draw("p")
    canvas.Print(output)
    numub_ratio = get_flux_ratios(numub2, numub1, fhc_all[1] + "/" + fhc_all[0] + " (FHC) #bar{#nu_{#mu}}")
    numub_ratio.Draw("p")
    canvas.Print(output)

    filename1 = dir1 + "/" + rhc_all[0] + "/sk_tuned13av2_13anom_" + rhc_all[0] + "_antinumode_fine.root"
    nue1, nueb1, numu1, numub1 = get_flux_weight_histograms(filename1)
    filename2 = dir2 + "/" + rhc_all[1] + "/sk_tuned13av3_13anom_" + rhc_all[1][:-3] + "_antinumode_fine.root"
    nue2, nueb2, numu2, numub2 = get_flux_weight_histograms(filename2)
    nue_ratio = get_flux_ratios(nue2, nue1, rhc_all[1] + "/" + rhc_all[0] + " (FHC) #nu_{e}")
    nue_ratio.Draw("p")
    canvas.Print(output)
    nueb_ratio = get_flux_ratios(nueb2, nueb1, rhc_all[1] + "/" + rhc_all[0] + " (FHC) #bar{#nu}_{e}")
    nueb_ratio.Draw("p")
    canvas.Print(output)
    numu_ratio = get_flux_ratios(numu2, numu1, rhc_all[1] + "/" + rhc_all[0] + " (FHC) #nu_{#mu}")
    numu_ratio.Draw("p")
    canvas.Print(output)
    numub_ratio = get_flux_ratios(numub2, numub1, rhc_all[1] + "/" + rhc_all[0] + " (FHC) #bar{#nu}_{#mu}")
    numub_ratio.Draw("p")
    canvas.Print(output + ")")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Check the ratio of two sets of flux tuning files")
    parser.add_argument("--dir-1", dest="dir1", type=str, help="A directory containing flux tuning files")
    parser.add_argument("--dir-2", dest="dir2", type=str, help="A directory containing flux tuning files")
    parser.add_argument("-o", "--output", type=str, help="The output file")
    args = parser.parse_args()
    plot_flux_ratios(args.dir1, args.dir2, args.output)
    plot_special_flux_ratios(args.dir1, args.dir2, "flux_special.pdf")
    plot_combined_flux_ratios(args.dir1, args.dir2, "flux_combined.pdf")
  
