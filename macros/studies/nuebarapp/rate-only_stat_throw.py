import ROOT
import copy
import math
import array
from array import array
import sys
import os
import subprocess

ROOT.PyConfig.IgnoreCommandLineOptions = True

env_var = "T2KOSCANA"
try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");

ROOT.gROOT.SetBatch(1)

valor   = ROOT.valor
util    = valor.util
conv    = valor.conv
opt     = valor.Opt.Instance()
t2kopt  = valor.T2KOpt.Instance()
env     = valor.Env.Instance()
rnd     = valor.Rnd.Instance()
arridx  = util.ArrIdx.Instance()

def main(infilename, infilename_data, ntoys, nthrows, tag, itoy_start):
    # Set up outputs
    Nevents_max = 100
    outname = tag + ".root"
    fout = ROOT.TFile(outname, "recreate")
    hist_nevents_weighted   = ROOT.TH1D("Nevents_weighted"  , "Nevents_weighted"  , Nevents_max, 0, Nevents_max)
    hist_nevents_unweighted = ROOT.TH1D("Nevents_unweighted", "Nevents_unweighted", Nevents_max, 0, Nevents_max)
    tree_out = ROOT.TTree("output_tree","output_tree")
    pval_left_weighted    = array('f', [-1.])   # For why an array is used see https://root.cern.ch/how/how-write-ttree-python
    pval_left_unweighted  = array('f', [-1.])
    pval_right_weighted   = array('f', [-1.])
    pval_right_unweighted = array('f', [-1.])
    Nevents_fakedata      = array('i', [-1 ])
    tree_out.Branch("pvalue_left_weighted"   , pval_left_weighted   , "pval_left_weighted/F"   )
    tree_out.Branch("pvalue_left_unweighted" , pval_left_unweighted , "pval_left_unweighted/F" )
    tree_out.Branch("pvalue_right_weighted"  , pval_right_weighted  , "pval_right_weighted/F"  )
    tree_out.Branch("pvalue_right_unweighted", pval_right_unweighted, "pval_right_unweighted/F")
    tree_out.Branch("Eventrate_fakedata"     , Nevents_fakedata     , "Eventrate_fakedata/I"  )

    # Read inputs
    fin = ROOT.TFile(infilename, "read")
    if not fin:
        print "Cannot open fake-data file", infilename, ". Exiting."
        util.Die(-1)
    fin_data = ROOT.TFile(infilename_data, "read")
    if not fin_data:
        print "Cannot open data file", infilename_data, ". Exiting."
        util.Die(-1)
    Nevents_data_hist = fin_data.Get("exp_data_RHC1Re_0")
    Nevents_data = Nevents_data_hist.Integral()
    for itoy in range(itoy_start, itoy_start + ntoys):
        if itoy % 1000 == 0:
            print "itoy =", itoy - itoy_start, "of", ntoys
        chisq = 0
        Nevents_nuebar = 0
        hist_nevents_weighted_itoy   = ROOT.TH1D("Nevents_weighted_" + str(itoy)  , "Nevents_weighted_" + str(itoy)  , Nevents_max, 0, Nevents_max)
        hist_nevents_unweighted_itoy = ROOT.TH1D("Nevents_unweighted_" + str(itoy), "Nevents_unweighted_" + str(itoy), Nevents_max, 0, Nevents_max)
        for isample in range(0, conv.kns):
            s = util.GetSKSample(isample)
            sample_name = util.SampleAsStringCompact(s)
            if s == conv.kSKsUndefined or sample_name == "":
                continue
            hname_exp_read =  "exp_data_" + sample_name + "_" + str(itoy)
            #print "*** hname exp = ", hname_exp_read
            hist_exp_sample = fin.Get(hname_exp_read)
            nexp_sample = hist_exp_sample.Integral()
            hname_obs_read =  "mock_data_" + util.SampleAsStringCompact(s) + "_" + str(itoy)
            hist_obs_sample = fin.Get(hname_obs_read)
            nobs_sample = hist_obs_sample.Integral()
            if(s == conv.kSKsRHCFCFV1RingELike):
                Nevents_nuebar = nexp_sample
            else:
                chisq += util.ChisqPoisson(nobs_sample, nexp_sample)
        # Fluctuate the nuebar event rate
        for ithrow in range(0, nthrows):
            Nevents_fluct = rnd.RndToyMCstat().Poisson(Nevents_nuebar)
            hist_nevents_weighted.Fill(Nevents_fluct, ROOT.TMath.Exp(-chisq*0.5))
            hist_nevents_weighted_itoy.Fill(Nevents_fluct, ROOT.TMath.Exp(-chisq*0.5))
            hist_nevents_unweighted.Fill(Nevents_fluct)
            hist_nevents_unweighted_itoy.Fill(Nevents_fluct)
        # Calculate p-values for this fake-data toy
        Nevents_fakedata     [0] = int(Nevents_nuebar)
        pval_left_unweighted [0] = p_value(hist_nevents_unweighted_itoy, math.floor(Nevents_data), "left" )
        pval_right_unweighted[0] = p_value(hist_nevents_unweighted_itoy, math.ceil(Nevents_data) , "right")
        pval_left_weighted   [0] = p_value(hist_nevents_weighted_itoy  , math.floor(Nevents_data), "left" )
        pval_right_weighted  [0] = p_value(hist_nevents_weighted_itoy  , math.ceil(Nevents_data) , "right")
        tree_out.Fill()
        del hist_nevents_weighted_itoy
        del hist_nevents_unweighted_itoy

    # Write out
    fout.cd()
    hist_nevents_weighted.Write("Nevents_nuebar_weighted")
    hist_nevents_unweighted.Write("Nevents_nuebar_unweighted")
    tree_out.Write("fake-data_tree")
    fout.Close()
    fin_data.Close()
    print "Done!"


def p_value(distribution, data, tail):
    ntot = distribution.Integral()
    ibin_data = distribution.FindBin(data)
    pval = -999
    if tail == "right":
        pval = distribution.Integral(ibin_data, distribution.GetNbinsX())
    elif tail == "left":
        pval = distribution.Integral(1, ibin_data)
    else:
        print "p_value: Unrecognised argument", tail, "exiting."
        util.Die(-1)
    return pval/ntot


if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Perform Poisson throws (default 10k) on the number of nuebar events for each of the fake-data toys in input file and store the weighted and unweighted flucts")
    parser.add_option('--infilename', dest='infilename', default="", help='Root file containing the fake-data toys to be fluctuated')
    parser.add_option('--infilename-data', dest='infilename_data', default="", help='Root file containing the data or the Asimov, in order to compute p-values')
    parser.add_option('--ntoys', dest='ntoys', default=100000, type=int, help='Number of fake-data toys in the input file')
    parser.add_option('--nthrows', dest='nthrows', default=10000, type=int, help='Number of Poisson throws')
    parser.add_option('--itoy-start', dest='itoy_start', default=0, type=int, help='Start reading the input file at this toy number')
    parser.add_option('--outfilename', dest='tag', default="nuebarapp_deltachisq", help='Output file name')


    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o, other = parser.parse_args()

    main(o.infilename, o.infilename_data, o.ntoys, o.nthrows, o.tag, o.itoy_start)

    if other:
        print "Warning - ignored arguments: ",other
                                                   
