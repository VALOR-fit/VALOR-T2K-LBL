#!/usr/bin/env python

import ROOT
import math
import sys
import os

from array import array
import numpy as np

ROOT.PyConfig.IgnoreCommandLineOptions = True

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)

def make_full_hybrid(outfile, datafile):
    df = ROOT.TFile(datafile)

    mcf1 = ROOT.TFile("spectra_osc_asimova_1Rmu.root")
    h1Rmu = mcf1.Get("curr_ereco_Vs_q2_1Rmu")
    mcf2 = ROOT.TFile("spectra_osc_asimova_RHC1Rmu.root")
    hRHC1Rmu = mcf2.Get("curr_ereco_Vs_q2_RHC1Rmu")
    mcf3 = ROOT.TFile("spectra_osc_asimova_1Re.root")
    h1Re = mcf3.Get("curr_ereco_Vs_q2_1Re")
    mcf4 = ROOT.TFile("spectra_osc_asimova_RHC1Re.root")
    hRHC1Re = mcf4.Get("curr_ereco_Vs_q2_RHC1Re")
    mcf5 = ROOT.TFile("spectra_osc_asimova_MultiRe.root")
    hMultiRe = mcf5.Get("curr_ereco_Vs_q2_MultiRe")

    hybrid = ROOT.TFile(outfile, "RECREATE")
    h1Rmu.Write("mock_data_1Rmu_0")
    hRHC1Rmu.Write("mock_data_RHC1Rmu_0")
    h1Re.Write("mock_data_1Re_0")
    hRHC1Re.Write("mock_data_RHC1Re_0")
    hMultiRe.Write("mock_data_MultiRe_0")


def make_hybrid(outfile, datafile, mcFHC1Rmu=None, mcRHC1Rmu=None):
    df = ROOT.TFile(datafile)

    if not mcFHC1Rmu is None:
        mcf1 = ROOT.TFile(mcFHC1Rmu)
        h1Rmu = mcf1.Get("curr_ereco_Vs_q2_1Rmu")
        h1Rmu.SetName("mock_data_1Rmu_0")
        h1Rmu.SetTitle("mock_data_1Rmu_0")
    if not mcRHC1Rmu is None:
        mcf2 = ROOT.TFile(mcRHC1Rmu)
        hRHC1Rmu = mcf2.Get("curr_ereco_Vs_q2_RHC1Rmu")
        hRHC1Rmu.SetName("mock_data_RHC1Rmu_0")
        hRHC1Rmu.SetTitle("mock_data_RHC1Rmu_0")
    h1Re = df.Get("mock_data_1Re_0")
    hRHC1Re = df.Get("mock_data_RHC1Re_0")
    hMultiRe = df.Get("mock_data_MultiRe_0")

    hybrid = ROOT.TFile(outfile, "RECREATE")
    h1Rmu.Write()
    hRHC1Rmu.Write()
    h1Re.Write()
    hRHC1Re.Write()
    hMultiRe.Write()
    hybrid.Close()
    df.Close()
    mcf1.Close()
    mcf2.Close()

if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Plot spectra")
    parser.add_option("-a", "--fhc-mu-mc-file", dest="mcfile1", help='Input file containing Asimov 1Rmu spectra')
    parser.add_option("-b", "--rhc-mu-mc-file", dest="mcfile2", help='Input file containing Asimov RHC1Rmu spectra')
    parser.add_option("-o", "--output-file", dest="outfile", help='Output file for hyrbid spectra')
    parser.add_option("-d", "--data-file", dest="datafile", help='Input file containing data spectra')
    
    v = optparse.Values()
    parser.parse_args(values = v)
    o, other = parser.parse_args()
  
    make_full_hybrid(o.outfile, o.datafile)
  
    if other:
        print "Warning - ignored arguments: ", other
