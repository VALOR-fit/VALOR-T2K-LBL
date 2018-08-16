#!/usr/bin/env python

import ROOT
import argparse
import re
import sys
import numpy as np

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

def get_2d_contours(f, source, target):
    keys = [key.GetName() for key in f.GetListOfKeys()]
    contours = []
    i = 0
    while (source + str(i) in keys):
        contours.append(f.Get(source + str(i)))
        contours[-1].SetNameTitle(target + str(i), target + str(i))
        i += 1
    return contours

def get_2d_bestfit(f, source, target):
    bf = f.Get(source).Clone(target)
    return bf

def get_1d_contours(f_nh, f_ih, source, target):
    chisq_nh = f_nh.Get(source)
    dchisq_nh = chisq_nh.Clone(target + "nh")
    dchisq_nh.SetNameTitle(target + "nh", target + "nh")
    dchisq_nh.SetDirectory(0)

    chisq_ih = f_ih.Get(source)
    dchisq_ih = chisq_ih.Clone(target + "ih")
    dchisq_ih.SetNameTitle(target + "ih", target + "ih")
    dchisq_ih.SetDirectory(0)

    minimum = min(dchisq_nh.GetMinimum(), dchisq_ih.GetMinimum())

    for i in range(1, chisq_nh.GetNbinsX() + 1):
        dchisq_nh.SetBinContent(i, chisq_nh.GetBinContent(i) - minimum)
        dchisq_ih.SetBinContent(i, chisq_ih.GetBinContent(i) - minimum)

    return (dchisq_nh, dchisq_ih)

def get_appearance(dir):
    subdir =  "mixing13_deltacp/data/"
    # Appearance wrc normal
    f_app_wrc_nh = ROOT.TFile(dir + "/" + subdir + "/" + \
        "fit_oa2018_normal_reactor_mixing13_deltacp_nexp_1_0.022_-1.601-job0_-_0_0_contour.root", "OPEN")
    app_68_nh = get_2d_contours(f_app_wrc_nh, "g68_valor_", "app_wrc_nh_68_")
    app_90_nh = get_2d_contours(f_app_wrc_nh, "g90_valor_", "app_wrc_nh_90_")
    app_bf_nh = get_2d_bestfit(f_app_wrc_nh, "best_fit_point_marker", "app_wrc_nh_bf")
    # Appearance wrc inverted
    f_app_wrc_ih = ROOT.TFile(dir + "/" + subdir + "/" + \
        "fit_oa2018_inverted_reactor_mixing13_deltacp_nexp_1_0.022_-1.601-job0_-_0_0_contour.root", "OPEN")
    app_68_ih = get_2d_contours(f_app_wrc_ih, "g68_valor_", "app_wrc_ih_68_")
    app_90_ih = get_2d_contours(f_app_wrc_ih, "g90_valor_", "app_wrc_ih_90_")
    app_bf_ih = get_2d_bestfit(f_app_wrc_ih, "best_fit_point_marker", "app_wrc_ih_bf")
    # Appearance worc normal
    f_app_worc_nh = ROOT.TFile(dir + "/" + subdir + "/" + \
        "fit_oa2018_normal_t2k_mixing13_deltacp_nexp_1_0.022_-1.601-job0_-_0_0_contour.root", "OPEN")
    app_worc_68_nh = get_2d_contours(f_app_worc_nh, "g68_valor_", "app_worc_nh_68_")
    app_worc_90_nh = get_2d_contours(f_app_worc_nh, "g90_valor_", "app_worc_nh_90_")
    app_worc_bf_nh = get_2d_bestfit(f_app_worc_nh, "best_fit_point_marker", "app_worc_nh_bf")
    # Appearance worc inverted
    f_app_worc_ih = ROOT.TFile(dir + "/" + subdir + "/" + \
        "fit_oa2018_inverted_t2k_mixing13_deltacp_nexp_1_0.022_-1.601-job0_-_0_0_contour.root", "OPEN")
    app_worc_68_ih = get_2d_contours(f_app_worc_ih, "g68_valor_", "app_worc_ih_68_")
    app_worc_90_ih = get_2d_contours(f_app_worc_ih, "g90_valor_", "app_worc_ih_90_")
    app_worc_bf_ih = get_2d_bestfit(f_app_worc_ih, "best_fit_point_marker", "app_worc_ih_bf")

    out_file = ROOT.TFile("t2k_appearance_fits.root", "RECREATE")
    for contour in app_68_nh:
        contour.Write()
    for contour in app_90_nh:
        contour.Write()
    app_bf_nh.Write("app_wrc_nh_bf")
    for contour in app_68_ih:
        contour.Write()
    for contour in app_90_ih:
        contour.Write()
    app_bf_ih.Write("app_wrc_ih_bf")
    for contour in app_worc_68_nh:
        contour.Write()
    for contour in app_worc_90_nh:
        contour.Write()
    app_worc_bf_nh.Write("app_worc_nh_bf")
    for contour in app_worc_68_ih:
        contour.Write()
    for contour in app_worc_90_ih:
        contour.Write()
    app_worc_bf_ih.Write("app_worc_ih_bf")

    out_file.Close()
    f_app_wrc_nh.Close()
    f_app_wrc_ih.Close()

def get_atmospheric(dir):
    subdir =  "mixing23_dmsq23/data/"
    # Atmospheric wrc normal
    f_atm_wrc_nh = ROOT.TFile(dir + "/" + subdir + "/" + \
        "fit_oa2018_normal_reactor_mixing23_dmsq23_data_smeared_contour.root", "OPEN")
    atmos_68_nh = get_2d_contours(f_atm_wrc_nh, "g68_valor_", "atm_wrc_nh_68_")
    atmos_90_nh = get_2d_contours(f_atm_wrc_nh, "g90_valor_", "atm_wrc_nh_90_")
    atmos_bf_nh = get_2d_bestfit(f_atm_wrc_nh, "best_fit_point_marker", "atm_wrc_nh_bf")
    # Atmospheric wrc inverted
    f_atm_wrc_ih = ROOT.TFile(dir + "/" + subdir + "/" + \
        "fit_oa2018_inverted_reactor_mixing23_dmsq23_data_smeared_contour.root", "OPEN")
    atmos_68_ih = get_2d_contours(f_atm_wrc_ih, "g68_valor_", "atm_wrc_ih_68_")
    atmos_90_ih = get_2d_contours(f_atm_wrc_ih, "g90_valor_", "atm_wrc_ih_90_")
    atmos_bf_ih = get_2d_bestfit(f_atm_wrc_ih, "best_fit_point_marker", "atm_wrc_ih_bf")
    # theta23 wrc
    subdir =  "mixing23/data"
    f_theta23_wrc_nh = ROOT.TFile(dir + "/" + subdir + "/" + \
        "fit_oa2018_normal_reactor_mixing23_nexp_1_0.528_-999.000-job0_-_0_0_contour.root", "OPEN")
    f_theta23_wrc_ih = ROOT.TFile(dir + "/" + subdir + "/" + \
        "fit_oa2018_inverted_reactor_mixing23_nexp_1_0.528_-999.000-job0_-_0_0_contour.root", "OPEN")
    theta23_wrc_nh, theta23_wrc_ih = get_1d_contours(f_theta23_wrc_nh, f_theta23_wrc_ih, \
        "h_chisq", "theta23_wrc_")
    out_file = ROOT.TFile("t2k_atmospheric_fits.root", "RECREATE")
    for contour in atmos_68_nh:
        contour.Write()
    for contour in atmos_90_nh:
        contour.Write()
    atmos_bf_nh.Write("atm_wrc_nh_bf")
    for contour in atmos_68_ih:
        contour.Write()
    for contour in atmos_90_ih:
        contour.Write()
    atmos_bf_ih.Write("atm_wrc_ih_bf")
    theta23_wrc_nh.Write()
    theta23_wrc_ih.Write()
    out_file.Close()
    f_theta23_wrc_nh.Close()
    f_theta23_wrc_ih.Close()
    f_atm_wrc_nh.Close()
    f_atm_wrc_ih.Close()

def get_deltacp(dir):
    # deltacp wrc
    subdir =  "deltacp/data"
    f_dcp_wrc_nh = ROOT.TFile(dir + "/" + subdir + "/" + \
        "fit_oa2018_normal_reactor_deltacp_nexp_1_-1.601_-999.000-job0_-_0_0_contour.root", "OPEN")
    f_dcp_wrc_ih = ROOT.TFile(dir + "/" + subdir + "/" + \
        "fit_oa2018_inverted_reactor_deltacp_nexp_1_-1.601_-999.000-job0_-_0_0_contour.root", "OPEN")
    dcp_wrc_nh, dcp_wrc_ih = get_1d_contours(f_dcp_wrc_nh, f_dcp_wrc_ih, \
        "h_chisq", "deltacp_wrc_")
    out_file = ROOT.TFile("t2k_deltacp_fits.root", "RECREATE")
    dcp_wrc_nh.Write()
    dcp_wrc_ih.Write()
    out_file.Close()
    f_dcp_wrc_nh.Close()
    f_dcp_wrc_ih.Close()

def get_contours(dir):
    get_atmospheric(dir)
    get_appearance(dir)
    get_deltacp(dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create the official data files")
    parser.add_argument("-d", "--dir", type=str, dest="dir", help="The top-level directory containing all fits")
    args = parser.parse_args()

    if not args.dir is None:
        get_contours(args.dir)
