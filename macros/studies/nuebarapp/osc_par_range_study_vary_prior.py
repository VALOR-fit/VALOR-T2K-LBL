#!/usr/bin/env python

# Input a bunch of files (each output by generate_toymc_ensemble_nuebarapp.cxx) corresponding to different
# oscillation parameter prior ranges. This macro then compares the likelihood weights distributions and chisq
# TODO: look at nh and ih priors as well

import ROOT
import copy
from copy import deepcopy
import collections
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



# Define the prior ranges to test here

# dCP, special case, use two priors to cut out the Dchisq peak. So 1st prior is [-pi, priors[0]], 2nd peak is [priors[1], +pi]  (peak is at 1.5radians)
#priors = [
#    [1.2, 1.6],
#    [1.1, 1.8],
#    [0.9, 2.0],
#    [0.7, 2.2],
#    [0.5, 2.4],
#    [0.3, 2.6],
#    [0.1, 2.8]
#    [-0.1, 3.0],
#]

## Dmsq23:
#priors = [
#    [2.3e-3, 2.65e-3],
#    [2.3e-3, 2.67e-3],
#    [2.3e-3, 2.69e-3],
#    [2.3e-3, 2.71e-3]
#]


# mixing23:
priors = [
    [0.300, 0.700],
    [0.39, 0.65],
    [0.40, 0.63],
    [0.41, 0.625],
    [0.42, 0.62]
]

hists_itoy_names = [
    "exp_data_1Rmu_",
    "exp_data_1Re_",
    "exp_data_RHC1Rmu_",
    "exp_data_MultiRe_"
]


def main(nominal_filename, tag, chisq_min, param, target, ntoys, data_filename):
    # Read the data file and get the control sample spectra
    # TODO: want to do this properly by looping over conv.KSKSamples, but this is 0 with python (instead of the expected array). Why?!
    data = []
    data_file = ROOT.TFile(data_filename, "read")
    data.append(data_file.Get("mock_data_1Rmu_0"))
    data.append(data_file.Get("mock_data_1Re_0"))
    data.append(data_file.Get("mock_data_RHC1Rmu_0"))
    data.append(data_file.Get("mock_data_MultiRe_0"))
    
    # Read the nominal file and extract the probability of accepting the toy and the mixing parameter value
    nominal_file = ROOT.TFile(nominal_filename, "read")
    nominal_values = []
    nominal_exp_accepted = 0
    min_throw =  999999999
    max_throw = -999999999
    for itoy in range(0, ntoys):
        if (itoy + 1) % 1000 == 0:
            print "Reading toy", itoy
        chisq = 0
        for isample in range(0, len(hists_itoy_names)):
            hname = hists_itoy_names[isample] + str(itoy)
            hist = nominal_file.Get(hname)
            if not hist:
                print "Could not find a hist named", hname, "in the input file", nominal_filename
                util.Die(1)
            chisq += util.Chisq_Penalty_TH2D(data[isample],hist)    
            del hist
        prob = ROOT.TMath.Exp(-0.5*(chisq - chisq_min))
        fpl = nominal_file.Get("iparm_"+str(itoy))
        param_value = fpl.Value(util.PhysParamIdFromRawString(param))
        if param_value < min_throw:
            min_throw = param_value
        elif param_value > max_throw:
            max_throw = param_value
        nominal_values.append([prob, param_value])
        nominal_exp_accepted += prob    # incrememt the exprected number of accepted toys by toy prob * number of toys with this prob (which is one, as we're reading event-by-event)
    #print nominal_values    
    throws = ROOT.TH1D("throws","throws", 100, min_throw, max_throw)
    canv = ROOT.TCanvas("throws_canv", "throws_canv", 800, 600)
    for pair in nominal_values:
        throws.Fill(pair[1])
    throws.Draw()
    canv.Update()
    print "throws integral =", throws.Integral()
    canv.SaveAs(param+"_throws_"+str(ntoys)+"_toys.pdf")    
    print "For the nominal priors, we expect to accept", nominal_exp_accepted, "toys out of", ntoys

    for prior in priors:
        exp_accepted = 0
        exp_excluded = 0
        for pair in nominal_values:
            if param == "delta_cp":
                if (pair[1] >= -ROOT.TMath.Pi() and pair[1] <= prior[0]) or (pair[1] >= prior[1] and pair[1] <= ROOT.TMath.Pi()):
                    # The mixing parameter from this nominal toy is within the prior range
                    exp_accepted += pair[0]
                else:
                    # The nominal toy is outside of this prior range, so we would have excluded this toy
                    exp_excluded += pair[0]
            else:
                if pair[1] >= prior[0] and pair[1] <= prior[1]:
                    # The mixing parameter from this nominal toy is within the prior range
                    exp_accepted += pair[0]
                else:
                    # The nominal toy is outside of this prior range, so we would have excluded this toy
                    exp_excluded += pair[0]
        if exp_accepted == 0:
            print "For prior [", prior[0], ",", prior[1], "], there were no accepted toys"
            continue
        exp_excluded_scaled = exp_excluded * float(target)/float(exp_accepted)   # Scale the number of excluded toys to the number expected if we'd generated target number of accepted toys
        # Calculate the maximum change in p-value. At worst, all excluded toys could have been in the tail of interest and the toys that replace them are not, or vice-vera
        dpval = exp_excluded/float(exp_accepted)   
        print "For prior [", prior[0], ",", prior[1], "] with", ntoys, "total toys (accepted + rejected), we expect to accept", exp_accepted, "toys and exclude", exp_excluded, "toys from the nominal set of toys. Scaling to", target, "accepted toys, expect to reject", exp_excluded_scaled, "nominal toys. This gives a maximum change in p-value of", dpval


     
                
                

    del throws
    del canv
    print "Done!"


if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="See the effects of different priors on nuebar appearance toy generation.")
    parser.add_option('--tag', dest='tag', default="", help='Tag used when naming output files')
    parser.add_option('--param', dest='param', default="", help='The parameter being varied. e.g. sin^2(theta_23)')
    parser.add_option('--chisq-min', dest='chisq_min', default=0, type=float, help='Minimum chisq summed over the control samples')
    parser.add_option('--target-accepted', dest='target', default=20000, type=int, help='The target number of accepted toys you intend to generate after prior range optimisation')
    parser.add_option('--nominal-file', dest='nominal_filename', default="", help='The file containing toys generated with the nominal set of priors, includes the rejected toys too')
    parser.add_option('--data-file', dest='data_filename', default="", help='The file containing data or fake data')
    parser.add_option('--ntoys', dest='ntoys', default=0, type=int, help='Total number of toys (accepted + rejected) in the input nominal file')

    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o, other = parser.parse_args()

    main(o.nominal_filename, o.tag, o.chisq_min, o.param, o.target, o.ntoys, o.data_filename)

    if other:
        print "Warning - ignored arguments: ",other
