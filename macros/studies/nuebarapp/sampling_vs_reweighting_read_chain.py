#!/usr/bin/env python

import ROOT
import collections
from collections import namedtuple
import re

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

#
# Histogram parameters to tweak
#
nbins_dchisq    = 160 #700
low_dchisq      = -25. #-20
high_dchisq     = 55.  #50
low_nevents     = -1   #-1
high_nevents    = 60   #50

tree_name_shape = "fitnt"
tree_name_rate = "fake-data_tree"
rand = ROOT.TRandom3()


def main(infilename, tag, infilename_data, test_statistic, read_nevents, start_nevents, chisq_min, fakefilename_pattern, ntoys_fake_file, beta):
    print "reading from file(s):", infilename

    read_nevents += start_nevents
    if not (test_statistic == "rate+shape" or test_statistic == "rate-only"):
        print "Test statistic", test_statistic, "not understood. Exiting."
        util.Die(1)

    data = read_data(infilename_data, test_statistic)

    # Set up the TChains and add all input toy files
    chain = setup_chain(infilename, test_statistic)

    if test_statistic == "rate+shape":
        # Loop over the trees to read the Delta chisq values and fill the histograms
        print "Reading delta chisq values from the input trees"
        to_save = read_from_tree(chain, beta, start_nevents, read_nevents, data, "rate+shape", data.test_stat, chisq_min,
                                 "hDchisq_beta0", nbins_dchisq, low_dchisq, high_dchisq, fakefilename_pattern, ntoys_fake_file)
    else:
        to_save = read_from_tree(chain, beta, start_nevents, read_nevents, data, "rate-only", data.test_stat, chisq_min,
                               "hNevents_beta0", high_nevents - low_nevents, low_nevents, high_nevents, fakefilename_pattern, ntoys_fake_file)

    warn_over_under_flow(to_save.hist_accepted)
    warn_over_under_flow(to_save.hist_rejected)
    warn_over_under_flow(to_save.hist_all_toys)

    # Convert list of weights & tails into TVectorD
    vec_wghts_accepted = ROOT.TVectorD(len(to_save.wghts_accepted))
    vec_tails_accepted = ROOT.TVectorD(len(to_save.wghts_accepted))
    for i in range(0, len(to_save.wghts_accepted)):
        pair = to_save.wghts_accepted[i]
        vec_wghts_accepted[i] = pair[0]
        tail = pair[1]
        if tail == "left":
            tail = 0
        if tail == "right":  # Not elif because we only store info on the tail in the rate+shape analysis. For rate-only, pair[1] is the p-value
            tail = 1
        vec_tails_accepted[i] = tail

    vec_wghts_rejected = ROOT.TVectorD(len(to_save.wghts_rejected))
    vec_tails_rejected = ROOT.TVectorD(len(to_save.wghts_rejected))
    for i in range(0, len(to_save.wghts_rejected)):
        pair = to_save.wghts_rejected[i]
        vec_wghts_rejected[i] = pair[0]
        tail = pair[1]
        if tail == "left":
            tail = 0
        if tail == "right":  # Not elif because we only store info on the tail in the rate+shape analysis. For rate-only, pair[1] is the p-value
            tail = 1
        vec_tails_rejected[i] = tail

    vec_wghts_all_toys = ROOT.TVectorD(len(to_save.wghts_all_toys))
    vec_tails_all_toys = ROOT.TVectorD(len(to_save.wghts_all_toys))
    for i in range(0, len(to_save.wghts_all_toys)):
        pair = to_save.wghts_all_toys[i]
        vec_wghts_all_toys[i] = pair[0]
        tail = pair[1]
        if tail == "left":
            tail = 0
        if tail == "right":  # Not elif because we only store info on the tail in the rate+shape analysis. For rate-only, pair[1] is the p-value
            tail = 1
        vec_tails_all_toys[i] = tail

    print "Number of accepted toys =", to_save.hist_accepted.Integral()
    print "Sum weights of rejected toys =", to_save.hist_rejected.Integral()

    # Save outputs
    save_name = "sampling_vs_reweighting_input_hists_" + test_statistic + "_beta=" + str(beta) + "_" + tag + ".root"
    vec_save_name = "pvals"
    if test_statistic == "rate+shape":
        vec_save_name = "tails"
    print "Saving to file", save_name
    fout = ROOT.TFile(save_name, "recreate")
    to_save.hist_accepted.Write("accepted")
    to_save.hist_rejected.Write("rejected")
    to_save.hist_all_toys.Write("all_toys")
    vec_wghts_accepted.Write("vec_wghts_accepted")
    vec_wghts_rejected.Write("vec_wghts_rejected")
    vec_wghts_all_toys.Write("vec_wghts_all_toys")
    vec_tails_accepted.Write("vec_"+ vec_save_name +"_accepted")
    vec_tails_rejected.Write("vec_"+ vec_save_name +"_rejected")
    vec_tails_all_toys.Write("vec_"+ vec_save_name +"_all_toys")
    fout.Close()

    print "Done!"




def warn_over_under_flow(hist):
    overflows = hist.GetBinContent(hist.GetNbinsX()+1)
    underflows = hist.GetBinContent(0)
    text1 = "WARNING: histogram \""+ hist.GetTitle()+ "\" has"
    text2 = "points. Consider extending the histogram range (p-val calculation won't be accurate)."
    if overflows != 0:
        print text1, "overflow", text2
    if underflows != 0:
        print text1, "underflow", text2


def read_data(infilename_data, test_stat_name):
    if not (test_stat_name == "rate-only" or test_stat_name == "rate+shape"):
        print "Test statistic", test_stat_name, "not understood. Exiting."
        util.Die(-1)
    data_tuple = namedtuple('data', ['FHC1Rmu', 'FHC1Re', 'RHC1Rmu', 'FHCMultiRe', 'test_stat'])     # only interested in chisq w.r.t. the control samples, so exclude RHC1Re
    data = data_tuple(FHC1Rmu=0, FHC1Re=0, RHC1Rmu=0, FHCMultiRe=0, test_stat=0)
    tree_name = "fitnt"
    tree_data = ROOT.TTree()
    if infilename_data != "":
        file_data  = ROOT.TFile(infilename_data ,"READ")
        tree_data  = file_data.Get(tree_name)
        if tree_data.GetEntries() != 1:
            print "Input data tree has", tree_data.GetEntries(), "Expected one entry. Exiting."
            util.Die(-1)
        for event in tree_data:     # Reading trees in pyroot is weird, but there'll be a better way to do this than a for loop with only one iteration...
            data_stat = 0
            if test_stat_name == "rate+shape":
                data_stat = event.deltachisq01
                print "# nuebar data Dchisq =", data_stat
            else:
                data_stat = event.hErecoToy_RHC1Re.Integral()
                print "# nuebar events =", data_stat
            data = data._replace(
                FHC1Rmu    = event.hErecoToy_1Rmu,
                FHC1Re     = event.hErecoToy_1Re,
                RHC1Rmu    = event.hErecoToy_RHC1Rmu,
                FHCMultiRe = event.hErecoToy_MultiRe,
                test_stat  = data_stat
            )
    return data


def setup_chain(infilenames, test_statistic):
    tname = tree_name_shape
    if test_statistic == "rate-only":
        tname = tree_name_rate

    chain = ROOT.TChain(tname)
    for name in infilenames:
        chain.Add(name)

    return chain


def read_from_tree(chain, beta, start_nevents, read_nevents, data_hists, test_stat_name, test_stat_data, chisq_min, hname, nbins, xlow, xhigh, fakefilename_pattern, ntoys_fake_file):
    if not (test_stat_name == "rate-only" or test_stat_name == "rate+shape"):
        print "Test statistic", test_stat_name, "not understood. Exiting."
        util.Die(-1)
    hist_accepted  = ROOT.TH1D(hname+"_accepted", hname+"_accepted", nbins, xlow, xhigh)
    hist_rejected  = ROOT.TH1D(hname+"_rejected", hname+"_rejected", nbins, xlow, xhigh)
    hist_all_toys  = ROOT.TH1D(hname+"_all_toys", hname+"_all_toys", nbins, xlow, xhigh)
    wghts_accepted = []
    wghts_rejected = []
    wghts_all_toys = []

    #print "*** test: start_nevents =", start_nevents, "; read_nevents =", read_nevents      # TODO: remove when done testing
    itoy_fake = 0
    fake_file_name = ""
    fake_file = 0
    batch_num = -999
    ievent = 0
    for event in chain:
        if ievent < start_nevents:
            ievent += 1
            itoy_fake += 1
            continue
        if read_nevents > 0 and ievent >= read_nevents:
            break

        if ievent % 5000 == 0:
            print "reading event", ievent

        # When using re-weighting, then we'll need the control sample chisq to calculate the weights, which requires the
        # spectra from the fake data toy file.
        # Work out which fake data toy file we should read based on the name of the file containing marginalised toys
        # Note: Currently this is dependant on the way I've generated the fake data toys & marginalised them - I have n
        # fake data toy files (named batch 0 ... batch n) and each of these was converted into m files of marginalised
        # toys. So here I need to search for the number after the string "batch" in the marginalised toy file
        # TODO: If the way the toys are generated or files have been named changes from the above description, this code will need to be changed
        current_marg_file_name = chain.GetFile().GetPath()
        #print "*** test: reading marginalisation output file", current_marg_file_name       # TODO: remove when done testing
        partial = re.sub("^.*batch", "", current_marg_file_name)
        batch_num_str = re.sub("-.*", "", partial)

        fake_file_name_try_partial = re.sub("\[b\]", str(beta), fakefilename_pattern)
        fake_file_name_try = re.sub("\[i\]", batch_num_str, fake_file_name_try_partial)
        if fake_file_name != fake_file_name_try or not fake_file:
            if fake_file:
                fake_file.Close()
            fake_file_name = fake_file_name_try
            fake_file = ROOT.TFile(fake_file_name, "read")
            print "*** test: reading fake data toy file", fake_file_name       # TODO: remove when done testing

        # If we're about to read from a new fake data file, work out the starting index
        if int(batch_num_str) != batch_num:
            batch_num = int(batch_num_str)
            #print "**** test: batch_num =", batch_num, "; ntoys_fake_file =", ntoys_fake_file       # TODO: remove when done testing
            itoy_fake = batch_num * ntoys_fake_file
            if itoy_fake < start_nevents:
                itoy_fake += start_nevents - itoy_fake

        # Get spectra & calculate the chisq
        #print "*** test: reading fake specta with index", itoy_fake       # TODO: remove when done testing
        hexp_1Rmu    = fake_file.Get("exp_data_1Rmu_"    + str(itoy_fake))
        hexp_1Re     = fake_file.Get("exp_data_1Re_"     + str(itoy_fake))
        hexp_RHC1Rmu = fake_file.Get("exp_data_RHC1Rmu_" + str(itoy_fake))
        hexp_RHC1Re  = fake_file.Get("exp_data_RHC1Re_"  + str(itoy_fake))
        hexp_MultiRe = fake_file.Get("exp_data_MultiRe_" + str(itoy_fake))

        chisq = 0
        chisq += util.Chisq_Penalty_TH2D(data_hists.FHC1Rmu   , hexp_1Rmu)
        chisq += util.Chisq_Penalty_TH2D(data_hists.FHC1Re    , hexp_1Re)
        chisq += util.Chisq_Penalty_TH2D(data_hists.RHC1Rmu   , hexp_RHC1Rmu)
        chisq += util.Chisq_Penalty_TH2D(data_hists.FHCMultiRe, hexp_MultiRe)
        prob = ROOT.TMath.Exp(-0.5*(chisq-chisq_min))

        del hexp_1Rmu
        del hexp_1Re
        del hexp_RHC1Rmu
        del hexp_MultiRe

        accepted = False
        rand_throw = rand.Uniform(0,1)
        if rand_throw <= prob:
            accepted = True

        if test_stat_name == "rate+shape":
            test_stat_toy = event.deltachisq01
            hist_all_toys.Fill(test_stat_toy, prob)     # Fill with prob as the bin content so that we can compare all toys to the rejected toys (the latter requires the weights to be applied)
            if test_stat_toy <= test_stat_data:
                wghts_all_toys.append([prob,"left"])
            if test_stat_toy >= test_stat_data: # Not elif, because if deltachisq01 = Delta_chisq_data, then this value is used in both tails
                wghts_all_toys.append([prob,"right"])

            if accepted == True or accepted == "True":
                hist_accepted.Fill(test_stat_toy)   # Don't fill with prob = bin content, as we count the accepted toys as having been sampled so the weight is not explicitly needed
                if  test_stat_toy <= test_stat_data:
                    wghts_accepted.append([prob,"left"])
                if  test_stat_toy >= test_stat_data:
                    wghts_accepted.append([prob,"right"])
            else:
                hist_rejected.Fill(test_stat_toy, prob)
                if  test_stat_toy <= test_stat_data:
                    wghts_rejected.append([prob,"left"])
                if  test_stat_toy >= test_stat_data:
                    wghts_rejected.append([prob,"right"])

            #print "*** test: ievent =", ievent, "; Accepted =", accepted, "; Delta chisq =", test_stat_toy, "\n"      # TODO: remove when done testing
        else:
            test_stat_toy = hexp_RHC1Re.Integral()
            #hist_all_toys.Fill(test_stat_toy, prob)     # Fill with prob as the bin content so that we can compare all toys to the rejected toys (the latter requires the weights to be applied)   #TODO: delete
            rate_only_fill_from_hist(event, hist_all_toys, prob)    # Read rate-only Nevents hist from the input file into whatever histogram binning is required here
            rate_only_pval_toy = 0
            if beta == 0:
                rate_only_pval_toy = event.pvalue_right
            else:
                rate_only_pval_toy = event.pvalue_left
            wghts_all_toys.append([prob, rate_only_pval_toy])    # In the rate-only analysis we don't associate the weight with a tail, as we apply the weights as a weighted mean p-value. And sd we already have p-values for each toy, store the pval instead of the tail
            if accepted == True or accepted == "True":
                #hist_accepted.Fill(test_stat_toy)  #TODO: delete
                rate_only_fill_from_hist(event, hist_accepted, 1.)
                wghts_accepted.append([prob, rate_only_pval_toy])
            else:
                #hist_rejected.Fill(test_stat_toy, prob)  #TODO: delete
                rate_only_fill_from_hist(event, hist_rejected, prob)
                wghts_rejected.append([prob, rate_only_pval_toy])

        if ievent % 10000 == 0:
            print "reading the beta =", beta, "chain, event #", ievent
        ievent += 1
        itoy_fake += 1

    returns = namedtuple('read_tree', ['hist_accepted' , 'hist_rejected' , 'hist_all_toys',
                                       'wghts_accepted', 'wghts_rejected', 'wghts_all_toys'])
    to_return = returns(
        hist_accepted  = hist_accepted,
        hist_rejected  = hist_rejected,
        hist_all_toys  = hist_all_toys,
        wghts_accepted = wghts_accepted,
        wghts_rejected = wghts_rejected,
        wghts_all_toys = wghts_all_toys
    )
    return to_return


def rate_only_fill_from_hist(tree_event, hist_to_fill, prob):
    hNevents_fluct = tree_event.hist_nevents_fluct
    for ibin in range(1, hNevents_fluct.GetNbinsX()+1):
        read_center  = hNevents_fluct.GetBinCenter(ibin)
        read_content = hNevents_fluct.GetBinContent(ibin)
        hist_to_fill.Fill(read_center, prob * read_content)     # The weight is only dependant on the control samples, so is the same for each fluctuation of the event RHC1Re rate
    del hNevents_fluct


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(usage="%prog [options]..", description="Select which marginalised toys are labelled as accepted or rejected and store them in hists. Use the output of this macro as the input to sampling_vs_reweighting_pvals.py")
    parser.add_argument('--infilename', dest='infilename', nargs="+", help='Root file(s) output by marginalisation')
    parser.add_argument('--infilename_data', dest='infilename_data', default="", help='Root file output by marginalisation of the data')
    parser.add_argument('--fakefilename-pattern', dest='fakefilename_pattern', default="", help='Naming pattering of the fake data toy file(s) used to produce the --infilename_beta files. Replace beta with [b] and the batch number with [i]. e.g. tag_0_[b]_[i]-[i].root')
    parser.add_argument('--tag', dest='tag', default="nuebarapp_deltachisq", help='Tag used as the output file name')
    parser.add_argument('--test-statistic', dest='test_statistic', default="rate+shape", help='The test statistic used in the marginalisation and fit')
    parser.add_argument('--read-nevents', dest='read_nevents', default=-1, type=int, help='Read only the first n events (staring from --start-nevents) from the input trees. This works for the p-values of both analyses, but the rate-only graph will be plot will all events. Set to <= 0 to read all events (default)')
    parser.add_argument('--start-nevents', dest='start_nevents', default=0, type=int, help='Start reading the input trees from this event number')
    parser.add_argument('--beta', dest='beta', default=0, type=int, help='Value of beta corresponding to --infilename')
    parser.add_argument('--ntoys-fakedata-file', dest='ntoys_fake_file', default=50000, type=int, help='The number of toys in each fake data file in --fakefilename_pattern')
    parser.add_argument('--chisq-min', dest='chisq_min', default=0, type=float, help='Minimum control sample chisq')

    o = parser.parse_args()
    main(o.infilename, o.tag, o.infilename_data, o.test_statistic, o.read_nevents, o.start_nevents, o.chisq_min,
         o.fakefilename_pattern, o.ntoys_fake_file, o.beta)

    other = argparse.REMAINDER  # TODO: Test that this correctly identifies the remaining args
    if other:
        print "Warning - ignored arguments: ",other

 
 
