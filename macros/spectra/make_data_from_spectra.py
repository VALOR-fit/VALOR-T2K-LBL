#!/usr/bin/env python

import ROOT
import math
import sys
import os

from array import array

ROOT.PyConfig.IgnoreCommandLineOptions = True

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)

samples=["1Rmu", "RHC1Rmu", "1Re", "RHC1Re", "MultiRe"]

def make_data(outfile, mcfileprefix, mcsamples=[]):
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
            print "Skipping", sample

if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Plot spectra")
    parser.add_option("-o", "--output-file", dest="outfile", help='Output file prefix for data spectra')
    parser.add_option("-f", "--fake-data-prefix", dest="prefix", help='Full path and prefit to fake data spectra')
    
    v = optparse.Values()
    parser.parse_args(values = v)
    o, other = parser.parse_args()
  
    make_data(o.outfile, o.prefix, ["1Re","1Rmu","RHC1Re","RHC1Rmu","MultiRe"])

    if other:
        print "Warning - ignored arguments: ", other
