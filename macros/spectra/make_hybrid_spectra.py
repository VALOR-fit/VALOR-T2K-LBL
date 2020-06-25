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

samples=["1Rmu", "RHC1Rmu", "1Re", "RHC1Re", "MultiRe"]

def make_full_hybrid(outfile, datafile, mcfileprefix, mcsamples=[]):
    df = ROOT.TFile(datafile)
    hybrid = ROOT.TFile(outfile, "RECREATE")
    print "Output", outfile

    for sample in samples:
        if sample in mcsamples:
            mcf = ROOT.TFile(mcfileprefix + sample + ".root")
            hsample = mcf.Get("curr_ereco_Vs_q2_" + sample)
            hybrid.cd()
            hsample.Write("mock_data_" + sample + "_0")
            mcf.Close()
            print sample, "from MC"
        else:
            hsample = df.Get("mock_data_" + sample + "_0")
            hybrid.cd()
            hsample.Write("mock_data_" + sample + "_0")
            print sample, "from data"

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
    parser.add_option("-o", "--output-file-prefix", dest="outfile", help='Output file prefix for hybrid spectra')
    parser.add_option("-d", "--data-file", dest="datafile", help='Input file containing data spectra')
    parser.add_option("-f", "--fake-data-prefix", dest="prefix", help='Full path and prefit to fake data spectra')
    
    v = optparse.Values()
    parser.parse_args(values = v)
    o, other = parser.parse_args()
  
    # Mostly fake data
    make_full_hybrid(o.outfile + "fhc1rmu_data.root", o.datafile, o.prefix, ["1Re","RHC1Re","RHC1Rmu","MultiRe"])
    make_full_hybrid(o.outfile + "fhc1re_data.root", o.datafile, o.prefix, ["1Rmu","RHC1Re","RHC1Rmu","MultiRe"])
    make_full_hybrid(o.outfile + "rhc1rmu_data.root", o.datafile, o.prefix, ["1Re","1Rmu","RHC1Re","MultiRe"])
    make_full_hybrid(o.outfile + "rhc1re_data.root", o.datafile, o.prefix, ["1Re","1Rmu","RHC1Rmu","MultiRe"])
    make_full_hybrid(o.outfile + "fhc1pi_data.root", o.datafile, o.prefix, ["1Re","1Rmu","RHC1Re","RHC1Rmu"])
    make_full_hybrid(o.outfile + "bestfit_asimov.root", o.datafile, o.prefix, ["1Re","1Rmu","RHC1Re","RHC1Rmu","MultiRe"])
    make_full_hybrid(o.outfile + "disapp_1pi_data.root", o.datafile, o.prefix, ["1Re","RHC1Re"])

    # Mostly real data
#    make_full_hybrid(o.outfile + "fhc1rmu_data.root", o.datafile, o.prefix, ["1Rmu"])
#    make_full_hybrid(o.outfile + "fhc1re_data.root", o.datafile, o.prefix, ["1Re"])
#    make_full_hybrid(o.outfile + "rhc1rmu_data.root", o.datafile, o.prefix, ["RHC1Rmu"])
#    make_full_hybrid(o.outfile + "rhc1re_data.root", o.datafile, o.prefix, ["RHC1Re"])
#    make_full_hybrid(o.outfile + "fhc1pi_data.root", o.datafile, o.prefix, ["MultiRe"])

    if other:
        print "Warning - ignored arguments: ", other
