#!/usr/bin/env python

import ROOT
from copy import deepcopy
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

def main(infile_pattern, nfiles, ntoys, outfile_name, count_failed_toys, ignore_spectra):
    # Each file contains the following:
    #   10 TH2Ds per toy:
    #       mock_data_1Rmu_[itoy], exp_data_1Rmu_[itoy], mock_data_1Re_[itoy],exp_data_1Re_[itoy],
    #       mock_data_RHC1Rmu_[itoy], exp_data_RHC1Rmu_[itoy], mock_data_RHC1Re_[itoy], exp_data_RHC1Re_[itoy],
    #       mock_data_MultiRe_[itoy], exp_data_MultiRe_[itoy]
    #   6 Osc par TH1Ds * 2 for hierarchies * 2 for barred and non-barred params:
    #       hToyHisto_dcp_nh, hToyHisto_Mix13_nh, hToyHisto_Mix23_nh, hToyHisto_Dm23_nh, hToyHisto_Mix12_nh,
    #       hToyHisto_Dm12_nh,
    #       hToyHisto_dcp_ih etc.
    #       hToyHisto_dcpBar_nh etc.
    #   1 FitParamList per toy:
    #       iparm_[itoy]
    #   1 chisq distribution
    # Here [itoy] is the toy MC number and each file follows on from the previous, e.g. for 10 toys across 2 files
    # you'd have [itoy] = 0 to 4 in the 1st file and [itoy] = 5 to 9 in the 2nd file.

    hists_itoy_names = [
        "mock_data_1Rmu_[itoy]",
        "exp_data_1Rmu_[itoy]",
        "mock_data_1Re_[itoy]",
        "exp_data_1Re_[itoy]",
        "mock_data_RHC1Rmu_[itoy]",
        "exp_data_RHC1Rmu_[itoy]",
        "mock_data_RHC1Re_[itoy]",
        "exp_data_RHC1Re_[itoy]",
        "mock_data_MultiRe_[itoy]",
        "exp_data_MultiRe_[itoy]"
    ]

    pars = [
        "dcp",
        "Mix13",
        "Mix23",
        "Dm23",
        "Mix12",
        "Dm12"
    ]

    hists_ipar_names = []
    hists_ipar_names_accepted = []
    for bar in ["", "Bar"]:
        for mh in ["nh", "ih"]:
            for par in pars:
                hists_ipar_names.append("hToyHisto_" + par + bar + "_" + mh)
                if count_failed_toys == 1:
                    hists_ipar_names_accepted.append("hToyHisto_" + par + bar + "_" + mh + "_accepted")

    print "Initialise the array of osc par hists"
    hists_ipar = {}
    for name in hists_ipar_names:
        # Open the first file to get the properties of the hist
        file0_name = infile_pattern.replace("[x]", str(0))
        file0 = ROOT.TFile(file0_name, "read")
        file0_hist = file0.Get(name)
        if not file0_hist:
            print "Could not find a hist named", name, "in the input file", file0_name
            util.Die(1)
        hists_ipar[name] = deepcopy(file0_hist)

    hists_ipar_accepted = {}    
    if count_failed_toys == 1:
        for name in hists_ipar_names_accepted:
            # Open the first file to get the properties of the hist
            file0_name = infile_pattern.replace("[x]", str(0))
            file0 = ROOT.TFile(file0_name, "read")
            file0_hist = file0.Get(name)
            if not file0_hist:
                print "Could not find a hist named", name, "in the input file", file0_name
                util.Die(1)
            hists_ipar_accepted[name] = deepcopy(file0_hist)

    hist_chisq = 0
    hist_chisq_accepted = 0
    file0_name = infile_pattern.replace("[x]", str(0))
    file0 = ROOT.TFile(file0_name, "read")
    file0_hist = file0.Get("chisq_distribution")
    hist_chisq = deepcopy(file0_hist)
    if count_failed_toys == 1:
        file0_hist_accepted = file0.Get("chisq_distribution_accepted")
        hist_chisq_accepted = deepcopy(file0_hist_accepted)

    outfile = ROOT.TFile(outfile_name, "recreate")

    print "Loop over each input file and combine all objects"
    for ifile in range(0,nfiles):
        print "Opening file", ifile
        filename = infile_pattern.replace("[x]", str(ifile))
        file = ROOT.TFile(filename, "read")
        if not file:
            print "Could not find the input file", filename
            util.Die(1)
        outfile.cd()
        # Add the contents of the osc par hists to the corresponding hists in the output file
        # Skip the first file, as we've already read the osc par throw and chisq hists from it
        if ifile != 0:
            for hname in hists_ipar_names:
                hist = file.Get(hname)
                if not hist:
                    print "Could not find a hist named", hname, "in the input file", filename
                    util.Die(1)
                out_hist = hists_ipar[hname]
                out_hist.Add(hist)
            hist_chisq_in = file.Get("chisq_distribution")
            hist_chisq.Add(hist_chisq_in)
            if count_failed_toys == 1:
                hist_chisq_accepted_in = file.Get("chisq_distribution_accepted")
                hist_chisq_accepted.Add(hist_chisq_accepted_in)
                for hname in hists_ipar_names_accepted:
                    hist = file.Get(hname)
                    if not hist:
                        print "Could not find a hist named", hname, "in the input file", filename
                        util.Die(1)
                    out_hist = hists_ipar_accepted[hname]
                    out_hist.Add(hist)
        # Add each spectrum hist and FitParamList to the output file
        if ignore_spectra == 0:
            itoy_start = ifile*ntoys
            itoy_end = itoy_start + ntoys   # Omit the -1 here because the range() function will give numbers up to itoy_end -1
            for itoy in range(itoy_start, itoy_end):
                if itoy % 1000 == 0 and not itoy == 0:
                    print "     toy #", itoy
                for name in hists_itoy_names:
                    # Spectrum
                    hname = name.replace("[itoy]", str(itoy))
                    hist = file.Get(hname)
                    if not hist:
                        print "Could not find a hist named", hname, "in the input file", filename
                        util.Die(1)
                    hist.Write(hname)
                # FitParamList
                pname = "iparm_" + str(itoy)
                fpl = file.Get(pname)
                fpl.Write(pname)


    hist_chisq.Write("chisq_distribution")
    for name in hists_ipar_names:
        hists_ipar[name].Write(name)
    if count_failed_toys == 1:
        hist_chisq_accepted.Write("chisq_distribution_accepted")
        for name in hists_ipar_names_accepted:
            hists_ipar_accepted[name].Write(name)

    
    print "Done!"



if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Combine many fake-data toy root files into one")
    parser.add_option('--infile-pattern', dest='infile_pattern', default="", help='The names of the input files with the run number replaced with \'[x]\'. e.g. if input files are called file0.root, file1.root etc. this input would be file[x].root')
    parser.add_option('--nfiles', dest='nfiles', default=1, type=int, help='The number of input files')
    parser.add_option('--ntoys', dest='ntoys', default=1, type=int, help='The number of toys per input file')
    parser.add_option('--outfile-name', dest='outfile_name', default="", help='Name of the combined output root file')
    parser.add_option('--count-failed-toys', dest='count_failed_toys', default=0, type=int, help='If MC sampling was used to generate toys and the rejected toys were stored, use this option')
    parser.add_option('--ignore-spectra', dest='ignore_spectra', default=0, type=int, help='When e.g. optimising the range of a prior, the spectra are not needed, so to save a few hours, ignore them when combining files')

    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o, other = parser.parse_args()

    main(o.infile_pattern, o.nfiles, o.ntoys, o.outfile_name, o.count_failed_toys, o.ignore_spectra)

    if other:
        print "Warning - ignored arguments: ",other


