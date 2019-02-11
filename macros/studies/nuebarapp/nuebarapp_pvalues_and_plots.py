#!/usr/bin/env python

import ROOT
import collections
from collections import namedtuple
import copy
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

#
# Histogram parameters to tweak
#
nbins_dchisq    = 160
low_dchisq      = -25. 
high_dchisq     = 55.  
low_nevents     = -1   
high_nevents    = 60   
scale_max_rate  = 5
scale_max_shape = 2
leg_text_size   = 0.02 #0.025 is good for posters/presentations
line_thickness  = 2
draw_prelim     = 1
prelim_text     = "T2K Run 1-9 Preliminary" #"T2K Run 1-9 Work in Progress"

# global variables
tree_name = "fitnt"
rand = ROOT.TRandom3()
maximum_Dchisq  = 0
maximum_nevents = 0
label = ""
data_label = ""
label_comp = ""
tag = ""
read_nevents = 0
start_nevents = 0
test_statistic = ""


def main(infilename_beta0, infilename_beta1, tag_save, infilename_data, data_label_save,
            infilename_comp_nobs_beta0, infilename_comp_nobs_beta1,
            infilename_comp_Dchisq_beta0, infilename_comp_Dchisq_beta1, test_statistic_save,
            reweighting, valor_comp, label_save, label_comp_save, read_nevents_save, start_nevents_save,
            nev_comp, Dchisq_comp, chisq_min, fakefilename_beta0, fakefilename_beta1):
    global label; label = label_save
    global data_label; data_label = data_label_save
    global label_comp; label_comp = label_comp_save
    global tag; tag = tag_save
    global read_nevents; read_nevents = read_nevents_save
    global start_nevents; start_nevents = start_nevents_save
    global test_statistic; test_statistic = test_statistic_save
    global tree_name

    if not (test_statistic == "rate+shape" or test_statistic == "rate-only"):
        print "Test statistic [", test_statistic, "] not understood. Exiting."
        util.Die(1)

    # Read fake data trees from the input files
    file_beta0 = ROOT.TFile(infilename_beta0,"READ")
    file_beta1 = ROOT.TFile(infilename_beta1,"READ")
    if test_statistic == "rate-only":
        tree_name = "fake-data_tree"
    tree_beta0 = file_beta0.Get(tree_name)
    tree_beta1 = file_beta1.Get(tree_name)
    n_entries = tree_beta0.GetEntries()
    if n_entries != tree_beta1.GetEntries():
        print "WARNING: Input fake-data trees have", n_entries, "and", tree_beta1.GetEntries(), "for beta = 0,1, respectively. Expected to be equal"
        #util.Die(-1)

    data = read_data(infilename_data)

    if test_statistic == "rate+shape":
        # Loop over the trees to read the Delta chisq values and fill the histrograms
        print "Reading delta chisq values from the input trees"
        beta0 = read_from_tree(tree_beta0, 0, data, data.test_stat, chisq_min, "hDchisq_beta0", nbins_dchisq,
                               low_dchisq, high_dchisq, fakefilename_beta0, reweighting)
        beta1 = read_from_tree(tree_beta1, 1, data, data.test_stat, chisq_min, "hDchisq_beta1", nbins_dchisq,
                               low_dchisq, high_dchisq, fakefilename_beta1, reweighting)


        # Plot the Delta chisq distribuitons on a linear scale
        canvs = draw_test_statistic_distributions(beta0, beta1, data, "#Delta #chi^{2} = #chi^{2}(#beta=0) - #chi^{2}(#beta=1)",
                                                 "Delta chisq", "#Delta #chi^{2}")

        # Plot comparison hists on a linear scale
        plot_comparison_both(beta0.hist, infilename_comp_Dchisq_beta0, 0, beta1.hist, infilename_comp_Dchisq_beta1, 1,
                             canvs.linear, canvs.leg, valor_comp, Dchisq_comp, "linear")

        # Plot comparison hists on a log scale
        #  TODO: having to call plot_comparison_both() again isn't ideal, especially if we're getting comparison data from a tree. Try re-using the same canvas and getting the data line from somewhere (to adjust it's height)
        plot_comparison_both(beta0.hist, infilename_comp_Dchisq_beta0, 0, beta1.hist, infilename_comp_Dchisq_beta1, 1,
                             canvs.log, canvs.leg, valor_comp, Dchisq_comp, "log")


        # Calculate p-values
        print "\n\n***** p-values & Bayes factor *****"
        print "nuebar data Dchisq =", data.test_stat
        p_val = -999
        err_pval = -999

        p_val    = p_value_reweighting("right", beta0.wghts)    # Use the reweighting p-val method even for unweighted toys, because the unweighted toys have all weights = 1
        err_pval = err_pvalue_reweighting(p_val, beta0.wghts)   # The reweighting p-val error function should be fine for unweighted too, as all the weights are stored as 1 in the unweighted case
        print "Rate+shape (true beta = 0) p-value (right tail) =", p_val, "+/-", err_pval, "= (", sigma(p_val), err_sigma(p_val, err_pval), ") sigma"
        print "*** test: from hist, p-val = ", p_value_from_hist_unweighted(beta0.hist, data.test_stat, "right")

        p_val    = p_value_reweighting("left", beta1.wghts)
        err_pval = err_pvalue_reweighting(p_val, beta1.wghts)
        print "Rate+shape (true beta = 1) p-value (left tail) =", p_val, "+/-", err_pval, "= (", sigma(p_val), err_sigma(p_val, err_pval), ")sigma"
        print "*** test: from hist, p-val = ", p_value_from_hist_unweighted(beta1.hist, data.test_stat, "left")

        print "Bayes factor = exp(- Delta chisq / 2) =", ROOT.TMath.Exp(-0.5*data.test_stat)
        print "\n\n"
    else:
        # Loop over the trees to read the RHC1Re event-rate values and fill the histrograms
        print "Reading Nevents values from the input trees"
        beta0 = read_from_tree(tree_beta0, 0, data, data.test_stat, chisq_min, "hNevents_beta0",
                               high_nevents - low_nevents, low_nevents, high_nevents, fakefilename_beta0, reweighting)
        beta1 = read_from_tree(tree_beta1, 1, data, data.test_stat, chisq_min, "hNevents_beta1",
                               high_nevents - low_nevents, low_nevents, high_nevents, fakefilename_beta1, reweighting)

        # Plot the Delta chisq distribuitons on a linear scale
        canvs = draw_test_statistic_distributions(beta0, beta1, data, "Number of observed RHC 1Re-like events",
                                                  "Nevents", "Events")

        # Plot comparison hists on a linear scale
        plot_comparison_both(beta0.hist, infilename_comp_nobs_beta0, 0, beta1.hist, infilename_comp_nobs_beta1, 1,
                             canvs.linear, canvs.leg, valor_comp, nev_comp, "linear")

        # Plot comparison hists on a log scale
        #  TODO: having to call plot_comparison_both() again isn't ideal, especially if we're getting comparison data from a tree. Try re-using the same canvas and getting the data line from somewhere (to adjust it's height)
        plot_comparison_both(beta0.hist, infilename_comp_nobs_beta0, 0, beta1.hist, infilename_comp_nobs_beta1, 1,
                             canvs.log, canvs.leg, valor_comp, nev_comp, "log")


        # Calculate p-values
        print "\n\n***** p-values *****"
        print "# nuebar data events =", data.test_stat
        p_val = -999
        err_pval = -999

        p_val = mean_p_value_rate_only_reweighting(beta0.wghts)    # Use the reweighting p-val method even for unweighted toys, because the unweighted toys have all weights = 1
        err_pval = err_pvalue_unweighted(p_val, len(beta0.wghts))  # TODO: error for the reweighted rate-only p-values
        print "Rate-only (true beta = 0) p-value (right tail) =", p_val, "+/-", err_pval, "= (", sigma(p_val), err_sigma(p_val, err_pval), ") sigma"

        p_val = mean_p_value_rate_only_reweighting(beta1.wghts)
        err_pval = err_pvalue_unweighted(p_val, len(beta1.wghts))
        print "Rate-only (true beta = 1) p-value (left tail) =", p_val, "+/-", err_pval, "= (", sigma(p_val), err_sigma(p_val, err_pval), ") sigma"
        print "\n\n"

    # Save outputs
    fout = ROOT.TFile(tag+"_" + test_statistic + "_hists.root","recreate")
    beta0.hist.Write("beta0")
    beta1.hist.Write("beta1")
    canvs.log.Write(test_statistic + "_canvas")
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
    return float(sum_toys_tail) / float(total_toys)


# p-value calculation when re-weighting is used (so we calculate with the weights)
def p_value_reweighting(tail, wght_tail_list):
    # Format of wght_tail_list is [[weight1, tail1], [weight2, tail2], ...]
    sum_wghts_tail = sum(pair[0] if pair[1] == tail else 0 for pair in wght_tail_list)
    total_wghts = sum(pair[0] for pair in wght_tail_list)
    return float(sum_wghts_tail) / float(total_wghts)


def mean_p_value_rate_only(wght_pval_list):
    # Format of wght_pval_list is [[weight1, toy_pval1], [weight2, toy_pval2], ...]
    sum_pval = sum(pair[1] for pair in wght_pval_list)
    sum_wghts = len(wght_pval_list) # Using MC sampling, all weights are treated as being 1
    return float(sum_pval) / float(sum_wghts)


def mean_p_value_rate_only_reweighting(wght_pval_list):
    # Format of wght_pval_list is [[weight1, toy_pval1], [weight2, toy_pval2], ...]
    numerator = sum(pair[0] * pair[1] for pair in wght_pval_list)
    sum_wghts = sum(pair[0] for pair in wght_pval_list)
    return float(numerator) / float(sum_wghts)


# p-value error calculation when re-weighting is not used
def err_pvalue_unweighted(pval, Ntoys):
    #print "Calculating p-value errors under the assumption that MC sampling was used to generate toys"
    return ROOT.TMath.Sqrt(pval*(1-pval)/Ntoys)


# p-value error calculation when re-weighting is used (so we calculate with the weights)
def err_pvalue_reweighting(pval, wght_tail_list):
    sum_wghts = sum(pair[0] for pair in wght_tail_list)
    sumsq_wghts_left  = sum(pair[0]*pair[0] if pair[1] == "left"  else 0 for pair in wght_tail_list)
    sumsq_wghts_right = sum(pair[0]*pair[0] if pair[1] == "right" else 0 for pair in wght_tail_list)
    print "sum_wghts =", sum_wghts
    print "sumsq_wghts_left =", sumsq_wghts_left
    print "sumsq_wghts_right =", sumsq_wghts_right
    return ROOT.TMath.Sqrt(ROOT.TMath.Power((1 - pval)/sum_wghts,2) * sumsq_wghts_right + ROOT.TMath.Power(pval/sum_wghts,2) * sumsq_wghts_left)


# Convert an error on a p-value to an asymmetric error on the number of sigma
def err_sigma(p_val, err_pval):
    sigma_up   = sigma(p_val - err_pval)    # small p-val = large sigma
    sigma_down = sigma(p_val + err_pval)
    err_up     = sigma_up - sigma(p_val)
    err_down   = sigma(p_val) - sigma_down
    return "+ " + str(err_up )+ " - " + str(err_down)


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


def plot_pval_dist(hist, beta):    # TODO: use this?
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


def count_pvals_below_threshold(hist, threshold):    # TODO: use this?
    count = 0
    threshold_bin = hist.FindBin(threshold) - 1 # -1 because we want pvals < threshold
    frac = hist.Integral(1, threshold_bin) / hist.Integral(1, hist.GetNbinsX())
    return frac


def make_data_comparison_line(tstat_comp, leg):
    if tstat_comp != -999:
        max = 0
        if test_statistic == "rate+shape":
            max = maximum_Dchisq
        else:
            max = maximum_nevents
        data_line_comp = ROOT.TLine(tstat_comp, 0, tstat_comp, max)
        data_line_comp.SetLineStyle(5)
        data_line_comp.SetLineWidth(line_thickness)
        #data_line_comp.Draw("same")
        #leg.AddEntry(data_line_comp,"#Delta #chi^{2} ("+data_label+", "+label_comp+") = "+str("{0:.2f}".format(tstat_comp)),"l")
        return copy.deepcopy(data_line_comp)


def plot_comparison_both(valor_hist1, infile1name_comp, beta_file1, valor_hist2, infile2name_comp, beta_file2, canv, leg, valor_comp, tstat_comp, save_tag):
    outtag = ""
    label_save = label_comp + " (#beta = "
    draw_opt = "same"
    statistic_name = "#Delta #chi^{2}"
    if test_statistic == "rate-only":
        outtag = "_comparison_nevents"
        draw_opt = draw_opt + ",hist"
        statistic_name = "Events"
    elif test_statistic == "rate+shape":
        outtag = "_comparison_Deltachisq"
    else:
        print "Test statistic", test_statistic, "not understood in plot_comparison_both(). Exiting."
        util.Die(-1)
    
    new_canv = canv.DrawClone("new_comp_canv")
    new_canv.cd()
    new_leg = leg.DrawClone()
    
    # Draw data comparison line
    data_line_comp = make_data_comparison_line(tstat_comp, new_leg)
    if tstat_comp != -999:
        new_leg.AddEntry(data_line_comp, statistic_name + " ("+data_label+", "+label_comp+") = "+str("{0:.2f}".format(tstat_comp)),"l")
        data_line_comp.Draw("same")

    comp_h1 = ROOT.TH1D()
    comp_h2 = ROOT.TH1D()
    if infile1name_comp != "":
        comp_h1 = plot_comparison_single(valor_hist1, infile1name_comp, beta_file1, new_leg, valor_comp, label_save)
        comp_h1.Draw(draw_opt)
    if infile2name_comp != "":
        comp_h2 = plot_comparison_single(valor_hist2, infile2name_comp, beta_file2, new_leg, valor_comp, label_save)
        comp_h2.Draw(draw_opt)

    if infile1name_comp != "" or infile2name_comp != "":
        if test_statistic == "rate+shape":
            new_leg.SetX1NDC(0.9*leg.GetX1NDC())
        new_leg.SetY1NDC(0.92*leg.GetY1NDC())
        new_canv.Update()
        ext1 = ".pdf"
        ext2 = ".png"
        ext3 = ".eps"
        outfile_name  = tag + outtag + "_" + save_tag + ext1
        outfile_name2 = tag + outtag + "_" + save_tag + ext2
        outfile_name3 = tag + outtag + "_" + save_tag + ext3
        print "saving comparison image to file", outfile_name
        new_canv.SaveAs(outfile_name,"recreate")
        #new_canv.SaveAs(outfile_name2,"recreate")
        #new_canv.SaveAs(outfile_name3,"recreate")


def plot_comparison_single(valor_hist, infilename_comp, beta, leg, valor_comp, label_in):
    label_save = label_in + str(beta) + ")"
    file_comp = ROOT.TFile(infilename_comp,"READ")
    hist_comp = 0
    if test_statistic == "rate-only":
        hist_comp = ROOT.TH1D("comparison_"+str(beta), "comparison_"+str(beta), high_nevents - low_nevents, low_nevents, high_nevents)
    else:
        hist_comp = ROOT.TH1D("comparison_"+str(beta), "comparison_"+str(beta), nbins_dchisq, low_dchisq, high_dchisq)
    if valor_comp == False or valor_comp == "False" or valor_comp == "":
        # If reading a digitized P-Theta histogram, directly get the hist then re-bin it to the same as binning as the VALOR hist
        hname = "h0"
        hist_comp_raw = file_comp.Get(hname)
        for ibin in range(1, hist_comp_raw.GetNbinsX()+1):
            centre = hist_comp_raw.GetBinCenter(ibin)
            nentries = hist_comp_raw.GetBinContent(ibin)
            hist_comp.Fill(centre, nentries)
    else:
        # If comparing to another VALOR result, read the Dchisq or Nevents values into a hist
        tree = 0
        input_comp = 0
        if test_statistic == "rate+shape":
            tree = file_comp.Get("fitnt")
            # TODO: Currently there's no support for reweighted valor comparison hists. Put this in
            input_comp = read_from_tree(tree, beta, 0, 0, 0, "comparison_"+str(beta), nbins_dchisq, low_dchisq, high_dchisq, "", False)
        else:
            tree = file_comp.Get("fake-data_tree")
            input_comp = read_from_tree(tree, beta, 0, 0, 0, "comparison_"+str(beta), high_nevents - low_nevents, low_nevents, high_nevents, "", False)
        hist_comp = input_comp.hist


    area_norm(hist_comp, valor_hist)
    max = 0
    if test_statistic == "rate+shape":
        max = maximum_Dchisq
    else:
        max = maximum_nevents
    if beta == 0:
        set_draw_opt(hist_comp, 0, max, "", True)
    else:
        set_draw_opt(0, hist_comp, max, "", True)

    print "\n\n ***** Compare RMS values and integral between the VALOR and P-Theta histograms *****"
    if test_statistic == "rate-only":
        print "RMS value (axis 1) for rate-only distribution (beta =", str(beta), ", valor) =", valor_hist.GetRMS(1)
        print "RMS value (axis 1) for rate-only distribution (beta =", str(beta), ", P-theta) =", hist_comp.GetRMS(1)
        print "mean value (axis 1) for rate-only distribution (beta =", str(beta), ", valor) =", valor_hist.GetMean(1)
        print "mean value (axis 1) for rate-only distribution (beta =", str(beta), ", P-theta) =", hist_comp.GetMean(1)
        print "peak value for rate-only distribution (beta =", str(beta), ", valor) =", valor_hist.GetBinCenter(valor_hist.GetMaximumBin())
        print "peak value for rate-only distribution (beta =", str(beta), ", P-theta) =", hist_comp.GetBinCenter(hist_comp.GetMaximumBin())
    elif test_statistic == "rate+shape":
        print "RMS value (axis 1) for rate+shape distribution (beta =", str(beta), ", valor) =", valor_hist.GetRMS(1)
        print "RMS value (axis 1) for rate+shape distribution (beta =", str(beta), ", P-theta) =", hist_comp.GetRMS(1)
        print "mean value (axis 1) for rate+shape distribution (beta =", str(beta), ", valor) =", valor_hist.GetMean(1)
        print "mean value (axis 1) for rate+shape distribution (beta =", str(beta), ", P-theta) =", hist_comp.GetMean(1)
        print "peak value for rate+shape distribution (beta =", str(beta), ", valor) =", valor_hist.GetBinCenter(valor_hist.GetMaximumBin())
        print "peak value for rate+shape distribution (beta =", str(beta), ", P-theta) =", hist_comp.GetBinCenter(hist_comp.GetMaximumBin())

    print "valor_Hist integral =", valor_hist.Integral(), "; comp_hist integral =",  hist_comp.Integral()
    print "\n\n"
    new_hist = copy.deepcopy(hist_comp)
    # Only add the label if it does not already exist (this can happen as plot_comparison_both() is called twice on the same canvas [with log scale on and off])
    add_label = True
    leg_entries = leg.GetListOfPrimitives()
    entry = leg_entries.First()
    while(entry):
        if label_save == entry.GetLabel():
            add_label = False
            break
        entry = leg_entries.After(entry)
    if add_label:
        leg.AddEntry(new_hist,label_save,"l")
    return new_hist


def check_lengths(wghts_list, tails_list, name):
    if len(wghts_list) != len(tails_list):
        print "length of weights list named", name, "=", len(wghts_list)
        print "length of tails list named", name,  "=", len(tails_list)
        print "These should be equal. Exiting."
        util.Die(1)


def read_data(infilename_data):
    if not (test_statistic == "rate-only" or test_statistic == "rate+shape"):
        print "Test statistic", test_statistic, "not understood in read_data(). Exiting."
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
            if test_statistic == "rate+shape":
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


def read_from_tree(tree, beta, data_hists, test_stat_data, chisq_min, hname, nbins, xlow, xhigh, fakefilename, reweighting):
    if not (test_statistic == "rate-only" or test_statistic == "rate+shape"):
        print "Test statistic", test_statistic, "not understood in read_from_tree(). Exiting."
        util.Die(-1)
    hist = ROOT.TH1D(hname, hname, nbins, xlow, xhigh)
    wghts = []

    # Note on the input format of rate+shape and rate-only files:
    # rate+shape is tree of marginalised events (with no explicit weights)
    # rate-only is tree of unweighted p-values (output from stat fluct macro) and hists containing all stat throws for that toy
    # In both cases the inputs do not have explicit weight applied, instead we apply the explicit weights here if the reweighting option is used
    ievent = 0
    for event in tree:
        if ievent < start_nevents:
            ievent += 1
            continue
        if read_nevents > 0 and ievent >= read_nevents:
            break

        prob = 1
        if data_hists != 0 and reweighting == 1:
            # When using re-weighting, then we'll need the control sample chisq to calculate the weights
            fake_file = ROOT.TFile(fakefilename, "read")
            if not fake_file:
                print "ERROR: fake data file named [", fakefilename, "] not found. Exiting."
                util.Die(1)
            hexp_1Rmu    = fake_file.Get("exp_data_1Rmu_"    + str(ievent))
            hexp_1Re     = fake_file.Get("exp_data_1Re_"     + str(ievent))
            hexp_RHC1Rmu = fake_file.Get("exp_data_RHC1Rmu_" + str(ievent))
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

        if test_statistic == "rate+shape":
            test_stat_toy = event.deltachisq01
            hist.Fill(test_stat_toy, prob)
            if test_stat_toy <= test_stat_data:
                wghts.append([prob,"left"])
            if test_stat_toy >= test_stat_data: # Not elif, because if deltachisq01 = Delta_chisq_data, then this value is used in both tails
                wghts.append([prob,"right"])
        else:
            # rate-only
            if beta == 0:
                rate_only_pval_toy = event.pvalue_right
            else:
                rate_only_pval_toy = event.pvalue_left
            # In the rate-only analysis we don't associate the weight with a tail, as we apply the weights as a weighted mean p-value. And as we already have p-values for each toy, store the pval instead of the tail
            wghts.append([prob, rate_only_pval_toy])

            # Read rate-only Nevents hist from the input file into whatever histogram binning is required here
            hNevents_fluct = event.hist_nevents_fluct
            for ibin in range(1, hNevents_fluct.GetNbinsX()+1):
                read_center  = hNevents_fluct.GetBinCenter(ibin)
                read_content = hNevents_fluct.GetBinContent(ibin)
                hist.Fill(read_center, prob * read_content)     # The weight is only dependant on the control samples, so is the same for each fluctuation of the event RHC1Re rate
            del hNevents_fluct

        if ievent % 1000 == 0:
            print "reading the beta =", beta, "tree, event #", ievent
        ievent += 1

    returns = namedtuple('read_tree', ['hist', 'wghts'])
    to_return = returns(
        hist  = hist,
        wghts = wghts
    )
    return to_return


def set_draw_opt(hist_beta0, hist_beta1, max, xtitle, is_comparison):
    line_col_beta0 = 8
    line_col_beta1 = 2
    line_style_beta0 = 1
    line_style_beta1 = 1
    if is_comparison:
        line_col_beta0 = ROOT.TColor.kMagenta+1     #40   # colour 40 and 41 look good on monitors, but don't show up at all on projectors
        line_col_beta1 = 4  #41
        line_style_beta0 = 2    #3
        line_style_beta1 = 2    #3
    if hist_beta0:
        hist_beta0.SetLineColor(line_col_beta0)
        hist_beta0.SetLineStyle(line_style_beta0)
        hist_beta0.SetLineWidth(line_thickness)
        hist_beta0.GetXaxis().SetTitle(xtitle)
        hist_beta0.GetYaxis().SetTitle("Entries")
        hist_beta0.SetMaximum(max)
        hist_beta0.SetMinimum(0.5)   # TODO: change this min if we plot weighted hists only?
    if hist_beta1:
        hist_beta1.SetLineColor(line_col_beta1)
        hist_beta1.SetLineStyle(line_style_beta1)
        hist_beta1.SetLineWidth(line_thickness)
        hist_beta1.GetXaxis().SetTitle(xtitle)
        hist_beta1.GetYaxis().SetTitle("Entries")
        hist_beta1.SetMaximum(max)
        hist_beta1.SetMinimum(0.5)   # TODO: change this min if we plot weighted hists only?


def draw_test_statistic_distributions(beta0, beta1, data, xtitle, canv_name, data_line_name):
    # Set plotting options
    global maximum_Dchisq
    global maximum_nevents
    maximum = 1.1 * max(beta0.hist.GetMaximum(), beta1.hist.GetMaximum())
    if test_statistic == "rate+shape":   # Store the max for use in other methods
        maximum_Dchisq = maximum
    else:
        maximum_nevents = maximum

    set_draw_opt(beta0.hist, beta1.hist, maximum, xtitle, False)

    # Setup vertical line for data or Asimov
    data_line = ROOT.TLine(data.test_stat, 0.5, data.test_stat, maximum)
    data_line.SetLineStyle(2)
    data_line.SetLineWidth(line_thickness)

    # Setup legend
    leg_ylow = 0.8
    if data.test_stat != 0:     # If we read no data, the data test statistic will be 0
        leg_ylow = 0.78
    leg = ROOT.TLegend(0.65,leg_ylow,0.9,0.95)
    leg.SetTextSize(leg_text_size)
    leg.AddEntry(beta0.hist, label + " (#beta = 0)", "l")
    leg.AddEntry(beta1.hist, label + " (#beta = 1)", "l")
    if data.test_stat != 0:
        leg.AddEntry(data_line, data_line_name + " ("+data_label+") = "+str("{0:.2f}".format(data.test_stat)),"l")

    canv = ROOT.TCanvas(canv_name, data_line_name, 800, 600)
    beta0.hist.Draw()
    beta1.hist.Draw("same")
    if data.test_stat != 0:
        data_line.Draw("same")
    leg.Draw("same")
    if draw_prelim == 1:
        prelim = ROOT.TLatex()
        prelim.SetTextSize(0.05)
        prelim.DrawLatexNDC(0.45, 0.96, prelim_text)
    beta0.hist.SetMaximum(1.1*maximum)
    data_line.SetY2(1.1*maximum)
    canv.Update()

    save_name = tag + "_" + canv_name.replace(" ", "")
    print "saving", canv_name, "image to file", save_name + ".pdf"
    canv.SaveAs(save_name+".pdf","recreate")
    #canv.SaveAs(save_name+".eps","recreate")
    #canv.SaveAs(save_name+".png","recreate")

    # Plot the histograms on a log scale
    canv_log = canv.Clone(canv_name + "_logscale")
    beta0.hist.SetMaximum(scale_max_shape*maximum)
    data_line.SetY2(scale_max_shape*maximum)
    canv_log.SetLogy(1)
    canv_log.Update()
    save_name += "_logscale"
    print "saving", canv_name, "(log scale) image to file", save_name + ".pdf"
    canv_log.SaveAs(save_name+".pdf","recreate")
    #canv.SaveAs(save_name+".eps","recreate")
    #canv.SaveAs(save_name+".png","recreate")

    # Check for any underflow or overflow points that won't have been plot
    warn_over_under_flow(beta0.hist)
    warn_over_under_flow(beta1.hist)

    return_canvs = namedtuple('canvases', ['linear', 'log', 'leg'])     # TODO: probably better to only return the linear canvas and then in the draw comparison method I can re-apply the log scale there
    to_return = return_canvs(
        linear = canv,
        log = canv_log,
        leg = leg
    )
    return to_return


if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Plot the Delta chisq and event rate graphs and calculate p-values for the nuebar appearance analysis")
    parser.add_option('--infilename_beta0', dest='infilename_beta0', default="", help='Root file output by marginalisation, with fake data toys (distribution toys) for beta = 0')
    parser.add_option('--infilename_beta1', dest='infilename_beta1', default="", help='Root file output by marginalisation, with fake data toys (distribution toys) for beta = 1')
    parser.add_option('--infilename_data', dest='infilename_data', default="", help='Root file output by marginalisation of the data')
    parser.add_option('--infilename_comp_nobs_beta0', dest='infilename_comp_nobs_beta0', default="", help='Root file containing histogram of a fake data distribution to be compared with ours for beta = 0')
    parser.add_option('--infilename_comp_nobs_beta1', dest='infilename_comp_nobs_beta1', default="", help='Root file containing histogram of a fake data distribution to be compared with ours for beta = 1')
    parser.add_option('--infilename_comp_Dchisq_beta0', dest='infilename_comp_Dchisq_beta0', default="", help='Root file containing histogram of a fake data distribution to be compared with ours for beta = 0')
    parser.add_option('--infilename_comp_Dchisq_beta1', dest='infilename_comp_Dchisq_beta1', default="", help='Root file containing histogram of a fake data distribution to be compared with ours for beta = 1')
    parser.add_option('--fakefilename_beta0', dest='fakefilename_beta0', default="", help='Fake data toy file used to produce --infilename_beta0')
    parser.add_option('--fakefilename_beta1', dest='fakefilename_beta1', default="", help='Fake data toy file used to produce --infilename_beta1')
    parser.add_option('--data_label', dest='data_label', default="data", help='Label to be put in the legend entry corresponding to the data')
    parser.add_option('--tag', dest='tag', default="nuebarapp_deltachisq", help='Tag used as the output file name')
    parser.add_option('--label', dest='label', default="Fake-data", help='Label used for the legend')
    parser.add_option('--label-comp', dest='label_comp', default="P-theta Events", help='Label used for the legend for comparison hists')
    parser.add_option('--nev-comp', dest='nev_comp', default=-999, type=float, help='Data or Asimov number of RHC 1Re-like events from the results we\'re comparing to')
    parser.add_option('--Dchisq-comp', dest='Dchisq_comp', default=-999, type=float, help='Data or Asimov Delta chisq from the results we\'re comparing to')
    parser.add_option('--test-statistic', dest='test_statistic', default="rate+shape", help='The test statistic used in the marginalisation and fit')
    parser.add_option('--reweighting', dest='reweighting', action='store_true', default="False", help='Include this flag to explicitly use control sample weights in p-value calculations. Only use if fake data toys were not generated by MC sampling')
    parser.add_option('--valor-comp', dest='valor_comp', action='store_true', default="False", help='Include this flag to compare to another VALOR file isntead of the default P-Theta file')
    parser.add_option('--read-nevents', dest='read_nevents', default=-1, type=int, help='Read only the first n events (staring from --start-nevents) from the input trees. This works for the p-values of both analyses, but the rate-only graph will be plot will all events. Set to <= 0 to read all events (default)')
    parser.add_option('--start-nevents', dest='start_nevents', default=0, type=int, help='Start reading the input trees from this event number')
    parser.add_option('--chisq-min', dest='chisq_min', default=0, type=float, help='Minimum control sample chisq')


    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o, other = parser.parse_args()

    main(o.infilename_beta0, o.infilename_beta1, o.tag, o.infilename_data, o.data_label,
         o.infilename_comp_nobs_beta0, o.infilename_comp_nobs_beta1,
         o.infilename_comp_Dchisq_beta0, o.infilename_comp_Dchisq_beta1, o.test_statistic,
         o.reweighting, o.valor_comp, o.label, o.label_comp, o.read_nevents, o.start_nevents,
         o.nev_comp, o.Dchisq_comp, o.chisq_min, o.fakefilename_beta0, o.fakefilename_beta1)

    if other:
        print "Warning - ignored arguments: ",other



