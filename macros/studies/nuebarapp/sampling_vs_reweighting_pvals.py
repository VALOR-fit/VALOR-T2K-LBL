#!/usr/bin/env python

import ROOT
import collections
from collections import namedtuple

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
nbins_dchisq    = 250 #700
low_dchisq      = -20. #-20
high_dchisq     = 45.  #50
low_nevents     = -1   #-1
high_nevents    = 60   #50
scale_max_rate  = 5
scale_max_shape = 2
leg_text_size   = 0.03
line_thickness  = 2
draw_prelim     = 1
prelim_text     = "T2K Run 1-9 Preliminary"

tree_name = "fitnt"
rand = ROOT.TRandom3()


def main(infilename_beta0, infilename_beta1, tag, infilename_data, data_label, test_statistic, read_nevents, start_nevents, chisq_min, fakefilename_beta0, fakefilename_beta1):
    read_nevents += start_nevents
    if not (test_statistic == "rate+shape" or test_statistic == "rate-only"):
        print "Test statistic", test_statistic, "not understood. Exiting."
        util.Die(1)
    file_beta0 = ROOT.TFile(infilename_beta0,"READ")
    file_beta1 = ROOT.TFile(infilename_beta1,"READ")

    if test_statistic == "rate+shape":
        data = read_data(infilename_data, "rate+shape")

        tree_beta0 = file_beta0.Get(tree_name)
        tree_beta1 = file_beta1.Get(tree_name)
        n_entries = tree_beta0.GetEntries()
        if n_entries != tree_beta1.GetEntries():
            print "Input fake-data trees have", n_entries, "and", tree_beta1.GetEntries(), "for beta = 0,1, respectively. Expected to be equal. Exiting."
            util.Die(-1)

        # Loop over the trees to read the Delta chisq values and fill the histrograms
        print "Reading delta chisq values from the input trees"
        beta0 = read_from_tree(tree_beta0, 0, start_nevents, read_nevents, data, "rate+shape", data.test_stat,
                               chisq_min, "hDchisq_beta0", nbins_dchisq, low_dchisq, high_dchisq, fakefilename_beta0)
        beta1 = read_from_tree(tree_beta1, 1, start_nevents, read_nevents, data, "rate+shape", data.test_stat,
                               chisq_min, "hDchisq_beta1", nbins_dchisq, low_dchisq, high_dchisq, fakefilename_beta1)


        # Plot the Delta chisq distribuitons
        canv = draw_test_statistic_distributions(beta0, beta1, data, "#Delta #chi^{2} = #chi^{2}(#beta=0) - #chi^{2}(#beta=1)",
                                          tag, "Delta chisq", "#Delta #chi^{2}", data_label)

        # Calculate p-values
        print "\n\n***** p-values & Bayes factor *****"
        print "nuebar data Dchisq =", data.test_stat
        p_val = -999
        err_pval = -999

        p_val    = p_value_from_hist_unweighted(beta0.hist_accepted, data.test_stat, "right")
        err_pval = err_pvalue_unweighted(p_val, beta0.hist_accepted.Integral())
        print "Rate+shape (true beta = 0) p-value (right tail) for accepted MC sampling toys (from hist) =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val    = p_value_from_list_unweighted("right", beta0.wghts_accepted) # no need to calculate the error from the list, it would be exactly the same as from the hist
        print "Rate+shape (true beta = 0) p-value (right tail) for accepted MC sampling toys (from list) =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val    = p_value_reweighting("right", beta0.wghts_accepted)
        err_pval = err_pvalue_reweighting(p_val, beta0.wghts_accepted)
        print "Rate+shape (true beta = 0) p-value (right tail) for accepted reweighted toys =",  p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val    = p_value_reweighting("right", beta0.wghts_all_toys)
        err_pval = err_pvalue_reweighting(p_val, beta0.wghts_all_toys)
        print "Rate+shape (true beta = 0) p-value (right tail) for all reweighted toys =",       p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val    = p_value_reweighting("right", beta0.wghts_rejected)
        err_pval = err_pvalue_reweighting(p_val, beta0.wghts_rejected)
        print "Rate+shape (true beta = 0) p-value (right tail) for rejected reweighted toys =",  p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val = p_value_from_hist_unweighted(beta1.hist_accepted, data.test_stat, "left")
        err_pval = err_pvalue_unweighted(p_val, beta1.hist_accepted.Integral())
        print "Rate+shape (true beta = 1) p-value (left tail) for accepted MC sampling toys (from hist) =",  p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val    = p_value_from_list_unweighted("right", beta1.wghts_accepted) # no need to calculate the error from the list, it would be exactly the same as from the hist
        print "Rate+shape (true beta = 1) p-value (left tail) for accepted MC sampling toys (from list) =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val    = p_value_reweighting("left", beta1.wghts_accepted)
        err_pval = err_pvalue_reweighting(p_val, beta1.wghts_accepted)
        print "Rate+shape (true beta = 1) p-value (left tail) for accepted reweighted toys =",   p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val    = p_value_reweighting("left", beta1.wghts_all_toys)
        err_pval = err_pvalue_reweighting(p_val, beta1.wghts_all_toys)
        print "Rate+shape (true beta = 1) p-value (left tail) for all reweighted toys =",        p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val    = p_value_reweighting("left", beta1.wghts_rejected)
        err_pval = err_pvalue_reweighting(p_val, beta1.wghts_rejected)
        print "Rate+shape (true beta = 1) p-value (left tail) for rejected reweighted toys =",   p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        draw_control_chisq_hist(beta0, 0, tag)
        draw_control_chisq_hist(beta1, 1, tag)

        print "Bayes factor = exp(- Delta chisq / 2) =", ROOT.TMath.Exp(-0.5*data.test_stat)
        print "\n\n"
    else:
        #
        # rate-only
        #
        data = read_data(infilename_data, "rate-only")

        tree_beta0  = file_beta0.Get("fake-data_tree")
        tree_beta1  = file_beta1.Get("fake-data_tree")

        beta0 = read_from_tree(tree_beta0, 0, start_nevents, read_nevents, data, "rate-only", data.test_stat,chisq_min,
                               "hNevents_beta0", high_nevents - low_nevents, low_nevents, high_nevents, fakefilename_beta0)
        beta1 = read_from_tree(tree_beta1, 1, start_nevents, read_nevents, data, "rate-only", data.test_stat, chisq_min,
                               "hNevents_beta1", high_nevents - low_nevents, low_nevents, high_nevents, fakefilename_beta1)

        # Plot the Delta chisq distribuitons
        canv = draw_test_statistic_distributions(beta0, beta1, data, "Number of observed RHC 1Re-like events", tag,
                                                 "Nevents", "Events", data_label)

        # Calculate p-values
        print "\n\n***** p-values *****"
        print "# nuebar data events =", data.test_stat

        p_val = mean_p_value_rate_only(beta0.wghts_accepted)
        err_pval = err_pvalue_unweighted(p_val, len(beta0.wghts_accepted))
        print "Rate-only (true beta = 0) p-value (right tail) for accepted MC sampling toys =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val = mean_p_value_rate_only_reweighting(beta0.wghts_accepted)
        err_pval = err_pvalue_unweighted(p_val, len(beta0.wghts_accepted))  # For the rate-only analysis, the error on the weighted mean p-value uses the same equation as the unweighted case (the weights cancel out)
        print "Rate-only (true beta = 0) p-value (right tail) for accepted reweighted toys =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val = mean_p_value_rate_only_reweighting(beta0.wghts_all_toys)
        err_pval = err_pvalue_unweighted(p_val, len(beta0.wghts_all_toys))
        print "Rate-only (true beta = 0) p-value (right tail) for all reweighted toys =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val = mean_p_value_rate_only_reweighting(beta0.wghts_rejected)
        err_pval = err_pvalue_unweighted(p_val, len(beta0.wghts_rejected))
        print "Rate-only (true beta = 0) p-value (right tail) for rejected reweighted toys =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val = mean_p_value_rate_only(beta1.wghts_accepted)
        err_pval = err_pvalue_unweighted(p_val, len(beta1.wghts_accepted))
        print "Rate-only (true beta = 1) p-value (right tail) for accepted MC sampling toys =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val = mean_p_value_rate_only_reweighting(beta1.wghts_accepted)
        err_pval = err_pvalue_unweighted(p_val, len(beta1.wghts_accepted))  # For the rate-only analysis, the error on the weighted mean p-value uses the same equation as the unweighted case (the weights cancel out)
        print "Rate-only (true beta = 1) p-value (right tail) for accepted reweighted toys =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val = mean_p_value_rate_only_reweighting(beta1.wghts_all_toys)
        err_pval = err_pvalue_unweighted(p_val, len(beta1.wghts_all_toys))
        print "Rate-only (true beta = 1) p-value (right tail) for all reweighted toys =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"

        p_val = mean_p_value_rate_only_reweighting(beta1.wghts_rejected)
        err_pval = err_pvalue_unweighted(p_val, len(beta1.wghts_rejected))
        print "Rate-only (true beta = 1) p-value (right tail) for rejected reweighted toys =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"
        print "\n\n"

    # Save outputs
    fout = ROOT.TFile(tag+"_" + test_statistic + "_MCsampling_vs_reweighting_hists.root","recreate")
    beta0.hist_accepted.Write("accepted_beta0")
    beta1.hist_accepted.Write("accepted_beta1")
    beta0.hist_rejected.Write("rejected_beta0")
    beta1.hist_rejected.Write("rejected_beta1")
    canv.Write(test_statistic + "_canvas")
    fout.Close()

    print "Done!"


# p-value calculation when re-weighting is not used
def p_value_from_hist_unweighted(distribution, data, tail):
    #print "Calculating p-value under the assumption that MC sampling was used to generate toys"
    ntot = distribution.Integral()
    ibin_data = distribution.FindBin(data)
    pval = -999
    if tail == "right":
        pval = distribution.Integral(ibin_data,distribution.GetNbinsX())
    elif tail == "left":
        pval = distribution.Integral(1,ibin_data)
    else:
        print "p_value: Unrecognised argument", tail, "exiting."
        util.Die(-1)
    return pval/ntot


# p-value calculation similar to p_value_from_hist_unweighted(), but the p-value is calculated from a list instead
def p_value_from_list_unweighted(tail, wght_tail_list):
    # Format of wght_tail_list is [[weight1, tail1], [weight2, tail2], ...]
    sum_toys_tail = sum(1.0 if pair[1] == tail else 0 for pair in wght_tail_list)
    total_toys = len(wght_tail_list)
    return sum_toys_tail / total_toys


# p-value calculation when re-weighting is used (so we calculate with the weights)
def p_value_reweighting(tail, wght_tail_list):
    # Format of wght_tail_list is [[weight1, tail1], [weight2, tail2], ...]
    sum_wghts_tail = sum(pair[0] if pair[1] == tail else 0 for pair in wght_tail_list)
    total_wghts = sum(pair[0] for pair in wght_tail_list)
    return sum_wghts_tail / total_wghts


def mean_p_value_rate_only(wght_pval_list):
    # Format of wght_pval_list is [[weight1, toy_pval1], [weight2, toy_pval2], ...]
    sum_pval = sum(pair[1] for pair in wght_pval_list)
    sum_wghts = len(wght_pval_list) # Using MC sampling, all weights are treated as being 1
    return sum_pval / sum_wghts


def mean_p_value_rate_only_reweighting(wght_pval_list):
    # Format of wght_pval_list is [[weight1, toy_pval1], [weight2, toy_pval2], ...]
    numerator = sum(pair[0] * pair[1] for pair in wght_pval_list)
    sum_wghts = sum(pair[0] for pair in wght_pval_list)
    return numerator / sum_wghts


# p-value error calculation when re-weighting is not used
def err_pvalue_unweighted(pval, Ntoys):
    #print "Calculating p-value errors under the assumption that MC sampling was used to generate toys"
    return ROOT.TMath.Sqrt(pval*(1-pval)/Ntoys)


# p-value error calculation when re-weighting is used (so we calculate with the weights)
def err_pvalue_reweighting(pval, wght_tail_list):
    #print "Calculating p-value errors under the assumption that re-weighting was used to generate toys"
    sum_wghts = sum(pair[0] for pair in wght_tail_list)
    sumsq_wghts_left  = sum(pair[0]*pair[0] if pair[1] == "left"  else 0 for pair in wght_tail_list)
    sumsq_wghts_right = sum(pair[0]*pair[0] if pair[1] == "right" else 0 for pair in wght_tail_list)
    return ROOT.TMath.Sqrt(ROOT.TMath.Power((1 - pval)/sum_wghts,2) * sumsq_wghts_right + ROOT.TMath.Power(pval/sum_wghts,2) * sumsq_wghts_left)


def sigma(p_val):
    return ROOT.TMath.Sqrt(2.) * ROOT.TMath.ErfInverse(1 - p_val)


def warn_over_under_flow(hist):
    overflows = hist.GetBinContent(hist.GetNbinsX()+1)
    underflows = hist.GetBinContent(0)
    text1 = "WARNING: histogram \""+ hist.GetTitle()+ "\" has"
    text2 = "points. Consider extending the histogram range (p-val calculation won't be accurate)."
    if overflows != 0:
        print text1, "overflow", text2
    if underflows != 0:
        print text1, "underflow", text2

def area_norm(hist1, hist2, range_low = 0, range_high = 0):
    # Normalise hist1 according to the area of hist2
    integral1 = 0
    integral2 = 0
    if range_low == 0 and range_high == 0:
        integral1 = hist1.Integral("width")
        integral2 = hist2.Integral("width")
    else:
        hist1_bin_low  = hist1.FindBin(range_low)
        hist2_bin_low  = hist2.FindBin(range_low)
        hist1_bin_high = hist1.FindBin(range_high)
        hist2_bin_high = hist2.FindBin(range_high)
        integral1 = hist1.Integral(hist1_bin_low, hist1_bin_high)
        integral2 = hist2.Integral(hist2_bin_low, hist2_bin_high)
    hist1.Scale(integral2/integral1)


def plot_pval_dist(hist, beta, tag):    # TODO: use this
    scale_max = 5
    canv = ROOT.TCanvas("p-values", "p-values", 800, 600)
    canv.SetLogy(1)
    #canv.SetLogx(1)
    leg = ROOT.TLegend(0.58,0.6,0.9,0.95)
    leg.SetTextSize(0.03)

    # TODO: what range do we want? (set to 60 and read every 5th as a test)
    nproj = 60
    maxi = scale_max * hist.ProjectionY().GetMaximum()
    hist_proj_tot = hist.ProjectionX("_px_tot", 1, hist.GetNbinsX())
    hist_proj_tot.Draw("hist")
    hist_proj_tot.GetXaxis().SetTitle("p-value")
    hist_proj_tot.GetYaxis().SetTitle("Entries")
    hist_proj_tot.SetMaximum(maxi)
    Nplot = 0
    for ibin in range(1, nproj+1):
        if (ibin-1) % 5 != 0:
            continue
        hist_proj = hist.ProjectionX("_px_"+str(ibin-1), ibin, ibin+1)  # hist containing the p-value distribution corresponding to fake-data toys with Nevents = ibin-1
        hist_proj.SetFillColor(Nplot+2)
        hist_proj.SetFillStyle(3001+Nplot)
        draw_opt = "hist,same"
        hist_proj.DrawCopy(draw_opt)
        leg_text = "# Events = " + str(ibin-1)
        leg.AddEntry(hist_proj, leg_text, "F")
        canv.Update()
        Nplot += 1
    leg.AddEntry(hist_proj_tot, "Total", "l")
    leg.Draw("same") 
    canv.Update()
    tag_beta = "_pval_beta=0"
    if beta == 1:
        tag_beta = "_pval_beta=1"
    outfile_name = tag + tag_beta + "_nevents.pdf"
    outfile_name2 = tag + tag_beta + "_nevents.png"
    outfile_name3 = tag + tag_beta + "_nevents.eps"
    print "saving p-val image to file", outfile_name
    canv.SaveAs(outfile_name, "recreate")
    canv.SaveAs(outfile_name2, "recreate")
    canv.SaveAs(outfile_name3, "recreate")


def count_pvals_below_threshold(hist, threshold):    # TODO: use this
    count = 0
    threshold_bin = hist.FindBin(threshold) - 1 # -1 because we want pvals < threshold
    frac = hist.Integral(1, threshold_bin) / hist.Integral(1, hist.GetNbinsX())
    return frac


def check_lengths(wghts_list, tails_list, name):
    if len(wghts_list) != len(tails_list):
        print "length of weights list named", name, "=", len(wghts_list)
        print "length of tails list named", name,  "=", len(tails_list)
        print "These should be equal. Exiting."
        util.Die(1)


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


def read_from_tree(tree, beta, start_nevents, read_nevents, data_hists, test_stat_name, test_stat_data, chisq_min, hname, nbins, xlow, xhigh, fakefilename):
    if not (test_stat_name == "rate-only" or test_stat_name == "rate+shape"):
        print "Test statistic", test_stat_name, "not understood. Exiting."
        util.Die(-1)
    hist_accepted  = ROOT.TH1D(hname+"_accepted", hname+"_accepted", nbins, xlow, xhigh)
    hist_rejected  = ROOT.TH1D(hname+"_rejected", hname+"_rejected", nbins, xlow, xhigh)
    hist_all_toys  = ROOT.TH1D(hname+"_all_toys", hname+"_all_toys", nbins, xlow, xhigh)
    wghts_accepted = []
    wghts_rejected = []
    wghts_all_toys = []

    # In the rate+shape weighted analysis, we need the fake data toys to get the control sample chisq
    fake_file = ROOT.TFile(fakefilename, "read")

    ievent = 0
    for event in tree:
        if ievent < start_nevents:
            ievent += 1
            continue
        if read_nevents > 0 and ievent >= read_nevents:
            break

        # When using re-weighting, then we'll need the control sample chisq to calculate the weights
        hexp_1Rmu    = fake_file.Get("exp_data_1Rmu_"    + str(ievent))
        hexp_1Re     = fake_file.Get("exp_data_1Re_"     + str(ievent))
        hexp_RHC1Rmu = fake_file.Get("exp_data_RHC1Rmu_" + str(ievent))
        hexp_RHC1Re  = fake_file.Get("exp_data_RHC1Re_"  + str(ievent))
        hexp_MultiRe = fake_file.Get("exp_data_MultiRe_" + str(ievent))

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
        else:
            test_stat_toy = hexp_RHC1Re.Integral()
            hist_all_toys.Fill(test_stat_toy, prob)     # Fill with prob as the bin content so that we can compare all toys to the rejected toys (the latter requires the weights to be applied)
            rate_only_pval_toy = 0
            if beta == 0:
                rate_only_pval_toy = event.pvalue_right
            else:
                rate_only_pval_toy = event.pvalue_left
            wghts_all_toys.append([prob, rate_only_pval_toy])    # In the rate-only analysis we don't associate the weight with a tail, as we apply the weights as a weighted mean p-value. And sd we already have p-values for each toy, store the pval instead of the tail
            if accepted == True or accepted == "True":
                hist_accepted.Fill(test_stat_toy)
                wghts_accepted.append([prob, rate_only_pval_toy])
            else:
                hist_rejected.Fill(test_stat_toy, prob)
                wghts_rejected.append([prob, rate_only_pval_toy])

        if ievent % 10000 == 0:
            print "reading the beta =", beta, "tree, event #", ievent
        ievent += 1

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


def set_draw_opt(hist_accepted, hist_rejected, beta, max, xtitle):
    line_col_acc = 1
    line_col_rej = 1
    if beta == 0:
        line_col_acc = 2
        line_col_rej = 1
    else:
        line_col_acc = 3
        line_col_rej = 4
    hist_accepted.SetLineColor(line_col_acc)
    hist_rejected.SetLineColor(line_col_rej)
    hist_accepted.SetLineWidth(line_thickness)
    hist_rejected.SetLineWidth(line_thickness)
    hist_accepted.GetXaxis().SetTitle(xtitle)
    hist_accepted.GetYaxis().SetTitle("Entries")
    hist_accepted.SetMaximum(max)
    hist_accepted.SetMinimum(0.5)   # TODO: change this min if we plot weighted hists only?


def draw_control_chisq_hist(from_tree, beta, tag):
    nbins = 100
    max_chisq = 1.1*max(-2.0*ROOT.TMath.Log(pair[0]) for pair in from_tree.wghts_all_toys)
    hist_accepted = ROOT.TH1D("chisq_control_beta"+str(beta)+"_accepted", "chisq_control_beta"+str(beta)+"_accepted", nbins, 0, max_chisq)
    hist_rejected = ROOT.TH1D("chisq_control_beta"+str(beta)+"_rejected", "chisq_control_beta"+str(beta)+"_rejected", nbins, 0, max_chisq)
    hist_all_toys = ROOT.TH1D("chisq_control_beta"+str(beta)+"_all_toys", "chisq_control_beta"+str(beta)+"_all_toys", nbins, 0, max_chisq)
    for pair in from_tree.wghts_accepted:
        hist_accepted.Fill(-2.0*ROOT.TMath.Log(pair[0]))
    for pair in from_tree.wghts_rejected:
        hist_rejected.Fill(-2.0*ROOT.TMath.Log(pair[0]))
    for pair in from_tree.wghts_all_toys:
        hist_all_toys.Fill(-2.0*ROOT.TMath.Log(pair[0]))    # Just to get the maximum, don't actually need to plot it

    hist_accepted.SetLineColor(2)
    hist_rejected.SetLineColor(3)
    hist_accepted.SetFillColor(2)
    hist_rejected.SetFillColor(3)

    canv = ROOT.TCanvas("chisq_control_canv_beta"+str(beta), "chisq_control_canv_beta"+str(beta), 800, 600)
    stack = ROOT.THStack("chisq_control_stack_beta"+str(beta), "chisq_control_stack_beta"+str(beta))
    stack.Add(hist_rejected)
    stack.Add(hist_accepted)
    stack.SetMaximum(1.1*hist_all_toys.GetMaximum())
    stack.SetMinimum(0.5)
    stack.Draw()
    stack.GetXaxis().SetTitle("Control Sample #chi^{2} - #chi^{2}_{min.} for #beta = "+str(beta))
    stack.GetYaxis().SetTitle("Entries")

    leg = ROOT.TLegend(0.75, 0.75, 0.9, 0.95)
    leg.AddEntry(hist_rejected, "Rejected toys", "f")
    leg.AddEntry(hist_accepted, "Accepted toys", "f")
    leg.Draw()

    canv.Update()
    canv.SaveAs(tag+"_control_sample_chisq_beta"+str(beta)+".pdf")


def draw_test_statistic_distributions(beta0, beta1, data, xtitle, tag, canv_name, data_line_name, data_label):
    # Set plotting options
    maximum = 1.1 * max(beta0.hist_accepted.GetMaximum(), beta0.hist_rejected.GetMaximum(),
                               beta1.hist_accepted.GetMaximum(), beta1.hist_rejected.GetMaximum())
    set_draw_opt(beta0.hist_accepted, beta0.hist_rejected, 0, maximum, xtitle)
    set_draw_opt(beta1.hist_accepted, beta1.hist_rejected, 1, maximum, xtitle)

    # Setup vertical line for data or Asimov
    data_line = ROOT.TLine(data.test_stat, 0.5, data.test_stat, maximum)
    data_line.SetLineStyle(2)
    data_line.SetLineWidth(line_thickness)

    # Setup legend
    leg_ylow = 0.73
    if data.test_stat != 0:     # If we read no data, the data test statistic will be 0
        leg_ylow = 0.78
    leg = ROOT.TLegend(0.6,leg_ylow,0.9,0.95)
    leg.SetTextSize(leg_text_size)
    leg.AddEntry(beta0.hist_accepted, "Accepted (#beta = 0)","l")
    leg.AddEntry(beta0.hist_rejected, "Rejected (#beta = 0)","l")
    leg.AddEntry(beta1.hist_accepted, "Accepted (#beta = 1)","l")
    leg.AddEntry(beta1.hist_rejected, "Rejected (#beta = 1)","l")
    if data.test_stat != 0:
        leg.AddEntry(data_line, data_line_name + " ("+data_label+") = "+str("{0:.2f}".format(data.test_stat)),"l")

    canv = ROOT.TCanvas(canv_name, data_line_name, 800, 600)
    beta0.hist_accepted.Draw()
    beta0.hist_rejected.Draw("same")
    beta1.hist_accepted.Draw("same")
    beta1.hist_rejected.Draw("same")
    if data.test_stat != 0:
        data_line.Draw("same")
    leg.Draw("same")
    if draw_prelim == 1:
        prelim = ROOT.TLatex()
        prelim.SetTextSize(0.05)
        prelim.DrawLatexNDC(0.5, 0.96, prelim_text)
    canv.Update()

    save_name = tag + "_" + canv_name.replace(" ", "")
    print "saving", canv_name, "image to file", save_name + ".pdf"
    canv.SaveAs(save_name+".pdf","recreate")
    #canv.SaveAs(save_name+".eps","recreate")
    #canv.SaveAs(save_name+".png","recreate")

    # Plot the histograms on a log scale
    beta0.hist_accepted.SetMaximum(scale_max_shape*maximum)
    data_line.SetY2(scale_max_shape*maximum)
    canv.SetLogy(1)
    canv.Update()
    save_name += "_logscale"
    print "saving", canv_name, "(log scale) image to file", save_name + ".pdf"
    canv.SaveAs(save_name+".pdf","recreate")
    #canv.SaveAs(save_name+".eps","recreate")
    #canv.SaveAs(save_name+".png","recreate")

    # Check for any underflow or overflow points that won't have been plot
    warn_over_under_flow(beta0.hist_accepted)
    warn_over_under_flow(beta0.hist_rejected)
    warn_over_under_flow(beta1.hist_accepted)
    warn_over_under_flow(beta1.hist_rejected)

    return canv


if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Compare the MC sampling to reweighting Delta chisq and event rate graphs and calculate p-values")
    parser.add_option('--infilename_beta0', dest='infilename_beta0', default="", help='Root file output by marginalisation, with fake data toys (distribution toys) for beta = 0')
    parser.add_option('--infilename_beta1', dest='infilename_beta1', default="", help='Root file output by marginalisation, with fake data toys (distribution toys) for beta = 1')
    parser.add_option('--infilename_data', dest='infilename_data', default="", help='Root file output by marginalisation of the data')
    parser.add_option('--fakefilename_beta0', dest='fakefilename_beta0', default="", help='Fake data toy file used to produce --infilename_beta0')
    parser.add_option('--fakefilename_beta1', dest='fakefilename_beta1', default="", help='Fake data toy file used to produce --infilename_beta1')
    parser.add_option('--data_label', dest='data_label', default="data", help='Label to be put in the legend entry corresponding to the data')
    parser.add_option('--tag', dest='tag', default="nuebarapp_deltachisq", help='Tag used as the output file name')
    parser.add_option('--test-statistic', dest='test_statistic', default="rate+shape", help='The test statistic used in the marginalisation and fit')
    parser.add_option('--read-nevents', dest='read_nevents', default=-1, type=int, help='Read only the first n events (staring from --start-nevents) from the input trees. This works for the p-values of both analyses, but the rate-only graph will be plot will all events. Set to <= 0 to read all events (default)')
    parser.add_option('--start-nevents', dest='start_nevents', default=0, type=int, help='Start reading the input trees from this event number')
    parser.add_option('--chisq-min', dest='chisq_min', default=0, type=float, help='Minimum control sample chisq')


    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o, other = parser.parse_args()

    main(o.infilename_beta0, o.infilename_beta1, o.tag, o.infilename_data, o.data_label, o.test_statistic,
         o.read_nevents, o.start_nevents, o.chisq_min, o.fakefilename_beta0, o.fakefilename_beta1)

    if other:
        print "Warning - ignored arguments: ",other

 
 
