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


# Define the filenames of the input files here, with the prior range
# (format example: names["nominal_0_[x]_0-0.root"] = [0.3, 0.7])
# The first file should use the nominal priors
# Here the [x] in the filename replaces the value of beta,
# this is so that it knows to look at e.g. "nominal_0_0_0-0.root" and "nominal_0_1_0-0.root"
names = collections.OrderedDict()
names["../nominal/nominal_specify_priors_90000_toys_0_[x]_0-combined_no_spectra.root"] = [0.300, 0.700]
names[""]
#names["mixing23_0.355_0.68_90000_toys_0_[x]_0-combined.root"] = [0.355, 0.68]
#names["mixing23_0.375_0.66_90000_toys_0_[x]_0-combined.root"] = [0.375, 0.66]
#names["mixing23_0.39_0.65_90000_toys_0_[x]_0-combined.root"] = [0.39, 0.65]
#names["mixing23_0.41_0.625_90000_toys_0_[x]_0-combined.root"] = [0.41, 0.625]





def main(tag, chisq_min, param, target):
    canv_chisq = ROOT.TCanvas("chisq", "#chi^{2}", 800, 600)
    leg_chisq = ROOT.TLegend(0.50,0.8,0.9,0.95)
    max_chisq = -999
    hists_chisq = []

    canv_log_like_weights = ROOT.TCanvas("log_like_weights", "0.5 #times (#chi^{2} - #chi^{2}_{min.})", 800, 600)
    leg_log_like_weights = ROOT.TLegend(0.50,0.8,0.9,0.95)
    max_log_like_weights = -999
    hists_log_like_weights = []

    canv_prior_nh = ROOT.TCanvas("prior_nh", "prior_nh", 800, 600)
    leg_prior_nh = ROOT.TLegend(0.50,0.8,0.9,0.95)
    max_prior_nh = -999
    hists_prior_nh = []

    canv_prior_ih = ROOT.TCanvas("prior_ih", "prior_ih", 800, 600)
    leg_prior_ih = ROOT.TLegend(0.50,0.8,0.9,0.95)
    max_prior_ih = -999
    hists_prior_ih = []

    weights_nominal_fine_beta0 = ROOT.TH1D("nominal_weights_fine_binning_beta0", "nominal_weights_fine_binning_beta0", 1000000, 0, 200)
    weights_nominal_fine_beta1 = ROOT.TH1D("nominal_weights_fine_binning_beta1", "nominal_weights_fine_binning_beta1", 1000000, 0, 200)

    beta_vals = ["0", "1"]
    for beta in beta_vals:
        iter = 0
        for filename in names:
            filename_read = filename.replace("[x]",beta)
            print "Reading file", filename_read
            file = ROOT.TFile(filename_read, "read")

            # Draw the chisq distribution
            chisq_distribution_in = file.Get("chisq_distribution")
            chisq_distribution_accepted_in = file.Get("chisq_distribution_accepted")
            # Rebin, as they are stored using 100000 bins, which is not good for plotting
            chisq_distribution = chisq_distribution_in.Rebin(1000, "chisq_distribution_rebin")    # new hist should have 100 bins
            chisq_distribution_accepted = chisq_distribution_accepted_in.Rebin(1000, "chisq_distribution_accepted_rebin")    # new hist should have 100 bins
            append_hists(chisq_distribution, chisq_distribution_accepted, beta, hists_chisq, iter)
            max_chisq = set_draw(hists_chisq, leg_chisq, iter, names[filename], "#chi^{2} (#beta = "+beta+")", max_chisq, beta, param)

            # Draw the log likelihood weights distribution
            weights = ROOT.TH1D("weights", "weights", 100, 0, 200)
            weights_accepted = ROOT.TH1D("weights_accepted", "weights_accepted", 100, 0, 200)
            for ibin in range(1, chisq_distribution_in.GetNbinsX() + 1):
                wght = 0.5 * (chisq_distribution_in.GetBinCenter(ibin) - chisq_min)
                weights.Fill(wght, chisq_distribution_in.GetBinContent(ibin))
                weights_accepted.Fill(wght, chisq_distribution_accepted_in.GetBinContent(ibin))
                if iter == 0:
                    if beta == "0":
                        weights_nominal_fine_beta0.Fill(wght, chisq_distribution_in.GetBinContent(ibin))
                    else:
                        weights_nominal_fine_beta1.Fill(wght, chisq_distribution_in.GetBinContent(ibin))
            append_hists(weights, weights_accepted, beta, hists_log_like_weights, iter)
            max_log_like_weights = set_draw(hists_log_like_weights, leg_log_like_weights, iter, names[filename], "-ln(toy acceptance probability) (#beta = "+beta+")", max_log_like_weights, beta, param)

            # Draw the NH and IH priors
            prior_nh_in = file.Get("hToyHisto_" + param + "_nh")
            prior_nh_accepted_in = file.Get("hToyHisto_" + param + "_nh_accepted")
            prior_nh = prior_nh_in.Rebin(2, "hToyHisto_" + param + "_nh_rebin")
            prior_nh_accepted = prior_nh_accepted_in.Rebin(2, "hToyHisto_" + param + "_nh_accepted_rebin")
            append_hists(prior_nh, prior_nh_accepted, beta, hists_prior_nh, iter)
            max_prior_nh = set_draw(hists_prior_nh, leg_prior_nh, iter, names[filename], rootex_name(param) + " NH (#beta = "+beta+")", max_prior_nh, beta, param)

            prior_ih_in = file.Get("hToyHisto_" + param + "_ih")
            prior_ih_accepted_in = file.Get("hToyHisto_" + param + "_ih_accepted")
            prior_ih = prior_ih_in.Rebin(2, "hToyHisto_" + param + "_ih_rebin")
            prior_ih_accepted = prior_ih_accepted_in.Rebin(2, "hToyHisto_" + param + "_ih_accepted_rebin")
            append_hists(prior_ih, prior_ih_accepted, beta, hists_prior_ih, iter)
            max_prior_ih = set_draw(hists_prior_ih, leg_prior_ih, iter, names[filename], rootex_name(param) + " IH (#beta = "+beta+")", max_prior_ih, beta, param)               

            iter += 1

    # To group all graphs in the same document whilst ensuring that they're in the correct order, draw outside the loop
    draw(hists_chisq, canv_chisq, leg_chisq, max_chisq, tag, True, False, True)
    draw(hists_log_like_weights, canv_log_like_weights, leg_log_like_weights, max_log_like_weights, tag, False, False, True)
    draw(hists_prior_nh, canv_prior_nh, leg_prior_nh, max_prior_nh, tag, False, False, True)
    draw(hists_prior_ih, canv_prior_ih, leg_prior_ih, max_prior_ih, tag, False, True, True)

    summarise(hists_log_like_weights, names, weights_nominal_fine_beta0, weights_nominal_fine_beta1, target, param)


    print "Done!"


def append_hists(hist1, hist2, beta, list, iter):
    if beta == "0":
        list.append([copy.deepcopy(hist1), 0])
        list.append([copy.deepcopy(hist2), 0])
    elif beta == "1":
        list[2*iter][1] = copy.deepcopy(hist1)
        list[2*iter+1][1] = copy.deepcopy(hist2)
    else:
        print "Unrecognised value of beta", beta, ". Exiting."
        util.Die(1)


def set_draw(hists, leg, iter, prior_range, xtitle, maximum, beta, param):
    hist_all = hists[2*iter][int(beta)]
    hist_accepted = hists[2*iter+1][int(beta)]
    hist_all.SetLineColor(iter+1)
    hist_accepted.SetLineColor(iter+1)
    hist_accepted.SetLineStyle(2)
    if beta == "0":
        # Use the same legend for both beta = 0 and beta = 1
        if iter == 0:
            leg.AddEntry(hist_accepted, "Accepted toys", "l")
        leg.AddEntry(hist_all, rootex_name(param) + " [" + str(prior_range[0]) + ", "  + str(prior_range[1]) + "] all toys", "l")    
    if iter == 0:
        hist_all.GetXaxis().SetTitle(xtitle)
        hist_all.GetYaxis().SetTitle("Number of toys")
    return max(maximum, hist_all.GetMaximum())


def draw(hists, canv_beta0, leg, maximum, tag, first_to_save, last_to_save, logy):
    canv_beta1 = ROOT.TCanvas(canv_beta0.GetName()+"_beta1", canv_beta0.GetTitle()+"_beta1", canv_beta0.GetWindowWidth(), canv_beta0.GetWindowHeight())
    # Draw the beta = 0 hists on one page and the beta = 1 hists on the next page
    canv_active = 0
    for beta in range(0,2):
        first = True
        if beta == 0:
            canv_active = canv_beta0
        else:
            canv_active = canv_beta1
        canv_active.cd()
        if logy:
            canv_active.SetLogy(1)
        for pair in hists:
            if first:
                set_max = 1.1
                if logy:
                    pair[beta].SetMinimum(0.5)
                    set_max = 10
                pair[beta].SetMaximum(set_max*maximum)    
                pair[beta].DrawCopy()
                first = False
            else:
                pair[beta].DrawCopy("same")
        leg.Draw()  # The legends should be the same for beta = 0 and beta = 1
        if first_to_save and beta == 0:
            canv_active.SaveAs(tag+".pdf(")
        elif last_to_save and beta == 1:
            canv_active.SaveAs(tag+".pdf)")
        else:
            canv_active.SaveAs(tag+".pdf")


def summarise(hists, names, weights_nominal_fine_beta0, weights_nominal_fine_beta1, ntoys_accepted_target, param):
    for beta in range(0,2):
        print "For beta =", beta
        hist_wghts_fine_nominal = 0
        if beta == 0:
            hist_wghts_fine_nominal = weights_nominal_fine_beta0
        if beta == 1:
            hist_wghts_fine_nominal = weights_nominal_fine_beta1
        for i in range(0,len(hists)):
            if i % 2 == 0:
                print "  For prior =",rootex_name(param) + " [" + str(names.values()[i/2][0]) + ", "  + str(names.values()[i/2][1]) + "]"
                print "    total toys =", hists[i][beta].Integral()
            else:
                print "    accepted toys =", hists[i][beta].Integral()
            

def rootex_name(param):
    if param == "dcp":
        return "#delta_{CP}"
    elif param == "Mix13":
        return "sin^{2}2#theta_{13}"
    elif param == "Mix23":
        return "sin^{2}#theta_{23}"
    elif param == "Mix12":
        return "sin^{2}2#theta_{12}"
    elif param == "Dm12":
        return "#Delta m^{2}_{12}"
    elif param == "Dm23":
        return "#Delta m^{2}_{23}"
    else:
        print "Unexpected parameter name. Exiting."
        util.Die(1)
    return ""



if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="See the effects of different priors on nuebar appearance toy generation.")
    parser.add_option('--tag', dest='tag', default="", help='Tag used when naming output files')
    parser.add_option('--param', dest='param', default="", help='The parameter being varied. Choices: dcp, Mix13, Mix23, Mix12, Dm12, Dm23')
    parser.add_option('--chisq-min', dest='chisq_min', default=0, type=float, help='Minimum chisq summed over the control samples')
    parser.add_option('--target-accepted', dest='target', default=20000, type=float, help='The target number of accepted toys you intend to generate after prior range optimisation')

    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o, other = parser.parse_args()

    main(o.tag, o.chisq_min, o.param, o.target)

    if other:
        print "Warning - ignored arguments: ",other
