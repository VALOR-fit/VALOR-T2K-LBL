#!/usr/bin/env python

import ROOT
import copy
import math
import array
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

#
# Histogram parameters to tweak
#
nbins_dchisq  = 1100#50
low_dchisq    = -30.
high_dchisq   = 80.
low_nevents   = -1
high_nevents  = 65
scale_max     = 5
leg_text_size = 0.03

def main(infilename_beta0, infilename_beta1, tag, infilename_data, data_label, infilename_comp_nobs_beta0, infilename_comp_nobs_beta1, infilename_comp_Dchisq_beta0, infilename_comp_Dchisq_beta1, test_statistic, prior_predictive, valor_comp, label, label_comp):    
    tstat = True
    if not test_statistic == "rate+shape":
        if test_statistic == "rate-only":
            tstat = False
        else:
            print "Test statistic", test_statistic, "not understood. Exiting."
            util.Die(1)
         
    #
    # Setup histograms
    #
    hist_deltachisq01_beta0 = ROOT.TH1D("deltachisq01_beta0", "deltachisq01_beta0", nbins_dchisq              , low_dchisq , high_dchisq )
    hist_deltachisq01_beta1 = ROOT.TH1D("deltachisq01_beta1", "deltachisq01_beta1", nbins_dchisq              , low_dchisq , high_dchisq )
    hist_nevents_beta0      = ROOT.TH1D("nevents_beta0"     , "nevents_beta0"     , high_nevents - low_nevents, low_nevents, high_nevents)
    hist_nevents_beta1      = ROOT.TH1D("nevents_beta1"     , "nevents_beta1"     , high_nevents - low_nevents, low_nevents, high_nevents)
    nevents_data     = -999
    Delta_chisq_data = -999

    file_beta0 = ROOT.TFile(infilename_beta0,"READ")
    file_beta1 = ROOT.TFile(infilename_beta1,"READ")

    tree_name = "fitnt"
    tree_data = ROOT.TTree()
    if infilename_data != "":
        file_data  = ROOT.TFile(infilename_data ,"READ")
        tree_data  = file_data.Get(tree_name)
        if tree_data.GetEntries() != 1:
            print "Input data tree has", tree_data.GetEntries(), "Expected one entry. Exiting."
            util.Die(-1)


    if tstat == True:
        tree_beta0 = file_beta0.Get(tree_name)
        tree_beta1 = file_beta1.Get(tree_name)
        n_entries = tree_beta0.GetEntries()
        if n_entries != tree_beta1.GetEntries():
            print "Input fake-data trees have", n_entries, "and", tree_beta1.GetEntries(), "for beta = 0,1, respectively. Expected to be equal. Exiting."
            util.Die(-1)
                    
        #
        # Loop over the trees to read the Delta chisq values and fill the histrograms
        #
        print "Reading delta chisq values from the input trees"
        # Fill with Delta chisq = chisq(beta = 0 for marginalisation toys) - chisq(beta = 1 for marginalisation toys) ...
        # For distribution toys generated with beta = 0
        nevents = tree_beta0.GetEntries()
        ievent = 0
        for event in tree_beta0:
            hist_deltachisq01_beta0.Fill(event.deltachisq01)
            if ievent % 10000 == 0 and not ievent == 0:
                print "reading the beta = 0 tree, event #", ievent
            ievent += 1

        # And for distribution toys generated with beta = 1
        ievent = 0
        for event in tree_beta1:
            hist_deltachisq01_beta1.Fill(event.deltachisq01)
            if ievent % 10000 == 0 and not ievent == 0:
                print "reading the beta = 1 tree, event #", ievent
            ievent += 1

        # And for the data
        if infilename_data != "":
            for event in tree_data:
                Delta_chisq_data = event.deltachisq01   # There'll be a better way to do this than a for loop with only one iteration...
                print "# nuebar data Dchisq =", Delta_chisq_data

        #
        # Plot the Delta chisq histograms
        #
        hist_deltachisq01_beta0.SetLineColor(3)
        hist_deltachisq01_beta1.SetLineColor(2)
        hist_deltachisq01_beta0.GetXaxis().SetTitle("#Delta #chi^{2} = #chi^{2}(#beta=0) - #chi^{2}(#beta=1)")
        hist_deltachisq01_beta0.GetYaxis().SetTitle("Entries")
        maximum_dchisq = 1.1 * max(hist_deltachisq01_beta0.GetMaximum(), hist_deltachisq01_beta1.GetMaximum())
        hist_deltachisq01_beta0.SetMaximum(maximum_dchisq)
        data_line_dchisq = ROOT.TLine(Delta_chisq_data, 0, Delta_chisq_data, maximum_dchisq)
        data_line_dchisq.SetLineStyle(2)

        leg_ylow = 0.73
        if infilename_data != "":
            leg_ylow = 0.78
        leg = ROOT.TLegend(0.6,leg_ylow,0.9,0.95)
        leg.SetTextSize(leg_text_size)
        leg.AddEntry(hist_deltachisq01_beta0, "Fake-data #beta = 0","l")
        leg.AddEntry(hist_deltachisq01_beta1, "Fake-data #beta = 1","l")
        if infilename_data != "":
            leg.AddEntry(data_line_dchisq,"#Delta #chi^{2} ("+data_label+") = "+str("{0:.2f}".format(Delta_chisq_data)),"l")

        canv = ROOT.TCanvas("Delta chisq", "#Delta #chi^{2}", 800, 600)
        hist_deltachisq01_beta0.SetMinimum(0.5)
        hist_deltachisq01_beta0.Draw()
        hist_deltachisq01_beta1.Draw("same")
        if infilename_data != "":
            data_line_dchisq.Draw("same")
        leg.Draw("same")
        canv.Update()

        outfile_name = tag + "_Deltachisq.pdf"
        outfile_name2 = tag + "_Deltachisq.png"
        print "saving Delta chisq image to file", outfile_name
        canv.SaveAs(outfile_name,"recreate")
        canv.SaveAs(outfile_name2,"recreate")

        #
        # Plot the Delta chisq histograms on a log scale
        #
        hist_deltachisq01_beta0.SetMaximum(scale_max*maximum_dchisq)
        data_line_dchisq.SetY2(scale_max*maximum_dchisq)
        canv.SetLogy(1)
        canv.Update()
        outfile_name = tag + "_Deltachisq_logscale.pdf"
        outfile_name2 = tag + "_Deltachisq_logscale.png"
        print "saving Delta chisq log scale image to file", outfile_name
        canv.SaveAs(outfile_name,"recreate")
        canv.SaveAs(outfile_name2,"recreate")

        # If we have comparison data, plot it too
        plot_comparison_both(hist_deltachisq01_beta0, infilename_comp_Dchisq_beta0, 0, hist_deltachisq01_beta1,
                             infilename_comp_Dchisq_beta1, 1, "rate+shape", canv, leg, tag, valor_comp,
                             prior_predictive, label_comp)
        
        #
        # Check for any overflow points that won't have been plot
        #
        warn_over_under_flow(hist_deltachisq01_beta0)
        warn_over_under_flow(hist_deltachisq01_beta1)
        
        #
        # Calculate p-values
        #
        print "\n\n***** p-values & Bayes factor *****"
        print "nuebar data Dchisq =", Delta_chisq_data
        p_val = p_value(hist_deltachisq01_beta0, Delta_chisq_data, "right")
        err_pval = err_pvalue(p_val, hist_deltachisq01_beta0, prior_predictive)
        print "Rate+shape (true beta = 0) p-value (right tail) =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"
        p_val = p_value(hist_deltachisq01_beta1, Delta_chisq_data, "left" )
        err_pval = err_pvalue(p_val, hist_deltachisq01_beta1, prior_predictive)
        print "Rate+shape (true beta = 1) p-value (left  tail) =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"
        print "Bayes factor = exp(- Delta chisq / 2) =", ROOT.TMath.Exp(-0.5*Delta_chisq_data)
        print "\n\n"

        # Save outputs
        fout = ROOT.TFile(tag+"_rate+shape_hists.root","recreate")
        hist_deltachisq01_beta0.Write("rate+shape_stat_beta0")
        hist_deltachisq01_beta1.Write("rate+shape_stat_beta1")
        canv.Write("rate+shape_canvas")
        fout.Close()

    else:
        #
        # rate-only
        #
        hname = "Nevents_nuebar_weighted"
        if prior_predictive:
            hname = "Nevents_nuebar_unweighted"
        histfill_nevents_beta0 = file_beta0.Get(hname)
        histfill_nevents_beta1 = file_beta1.Get(hname)
        bin_width = histfill_nevents_beta0.GetBinCenter(2) - histfill_nevents_beta0.GetBinCenter(1)
        for ibin in range(0, histfill_nevents_beta0.GetNbinsX()+2):           
            nev = histfill_nevents_beta0.GetBinCenter(ibin) - bin_width/2.
            nentries = histfill_nevents_beta0.GetBinContent(ibin)
            hist_nevents_beta0.Fill(nev, nentries)
        for ibin in range(0, histfill_nevents_beta1.GetNbinsX()+2):           
            nev = histfill_nevents_beta1.GetBinCenter(ibin) - bin_width/2.
            nentries = histfill_nevents_beta1.GetBinContent(ibin)
            hist_nevents_beta1.Fill(nev, nentries)
        
        if infilename_data != "":
            for event in tree_data:
                nevents_data = event.hErecoToy_RHC1Re.Integral()   # There'll be a better way to do this than a for loop with only one iteration...
                print "# nuebar events =", nevents_data

        tree_beta0  = file_beta0.Get("fake-data_tree")
        tree_beta1  = file_beta1.Get("fake-data_tree")
        hist_pval_beta0 = ROOT.TH2D("pval_beta0", "pval_beta0", 200, -3, 3, high_nevents - low_nevents, low_nevents, high_nevents)    # 2D because want to separate out p-vals from fake-data toys with different event rates when plotting
        hist_pval_beta1 = ROOT.TH2D("pval_beta1", "pval_beta1", 200, -3, 3, high_nevents - low_nevents, low_nevents, high_nevents)
        hist_pval_nolog_beta0 = ROOT.TH2D("pval_nolog_beta0", "pval_nolog_beta0", 200, 0, 2, high_nevents - low_nevents, low_nevents, high_nevents)
        hist_pval_nolog_beta1 = ROOT.TH2D("pval_nolog_beta1", "pval_nolog_beta1", 200, 0, 2, high_nevents - low_nevents, low_nevents, high_nevents)
        ievent = 0
        for event in tree_beta0:
            if prior_predictive:
                hist_pval_beta0.Fill(ROOT.TMath.Log10(event.pvalue_right_unweighted), event.Eventrate_fakedata)
                hist_pval_nolog_beta0.Fill(event.pvalue_right_unweighted, event.Eventrate_fakedata)
            else:
                hist_pval_beta0.Fill(ROOT.TMath.Log10(event.pvalue_right_weighted), event.Eventrate_fakedata)
                hist_pval_nolog_beta0.Fill(event.pvalue_right_weighted, event.Eventrate_fakedata)
            if ievent % 10000 == 0 and not ievent == 0:
                print "reading the beta = 0 tree, event #", ievent
            ievent += 1
        ievent = 0
        for event in tree_beta1:
            if prior_predictive:
                hist_pval_beta1.Fill(ROOT.TMath.Log10(event.pvalue_left_unweighted), event.Eventrate_fakedata)
                hist_pval_nolog_beta1.Fill(event.pvalue_left_unweighted, event.Eventrate_fakedata)
            else:
                hist_pval_beta1.Fill(ROOT.TMath.Log10(event.pvalue_left_unweighted), event.Eventrate_fakedata)
                hist_pval_nolog_beta1.Fill(event.pvalue_left_unweighted, event.Eventrate_fakedata)
            if ievent % 10000 == 0 and not ievent == 0:
                print "reading the beta = 1 tree, event #", ievent
            ievent += 1

        #
        # Plot the nevents histograms
        #
        hist_nevents_beta0.SetLineColor(3)
        hist_nevents_beta1.SetLineColor(2)
        #hist_nevents_beta0.SetFillColor(3)
        #hist_nevents_beta1.SetFillColor(2)
        #hist_nevents_beta0.SetBarWidth(0.1)
        #hist_nevents_beta1.SetBarWidth(0.1)
        #hist_nevents_beta0.SetBarOffset(0.0)
        #hist_nevents_beta1.SetBarOffset(0.1)
        hist_nevents_beta0.SetBarWidth(0.2)
        hist_nevents_beta1.SetBarWidth(0.2)
        #hist_nevents_beta0.SetFillStyle(1001)
        #hist_nevents_beta1.SetFillStyle(1001)
        hist_nevents_beta0.SetLineStyle(1)
        hist_nevents_beta1.SetLineStyle(1)
        hist_nevents_beta1.SetBarOffset(0.2)
        hist_nevents_beta0.GetXaxis().SetTitle("Number of observed events")
        hist_nevents_beta0.GetYaxis().SetTitle("Entries")
        maximum_nevents = scale_max * max(hist_nevents_beta0.GetMaximum(), hist_nevents_beta1.GetMaximum())
        hist_nevents_beta0.SetMaximum(maximum_nevents)
        data_line_nevents = ROOT.TLine(nevents_data, 0, nevents_data, maximum_nevents)
        data_line_nevents.SetLineStyle(2)

        leg_nevents = ROOT.TLegend(0.50,0.8,0.9,0.95)
        leg_nevents.SetTextSize(leg_text_size)
        leg_nevents.AddEntry(hist_nevents_beta0, "Fake-data #beta = 0","l")
        leg_nevents.AddEntry(hist_nevents_beta1, "Fake-data #beta = 1","l")
        if infilename_data != "":
            leg_nevents.AddEntry(data_line_nevents,"Events ("+data_label+") = "+str("{0:.2f}".format(nevents_data)),"l")

        canv_nevents = ROOT.TCanvas("nevents", "## events", 800, 600)
        canv_nevents.SetLogy(1)
        #hist_nevents_beta0.Draw("B")
        #hist_nevents_beta1.Draw("same,B")
        hist_nevents_beta0.SetMinimum(0.5)
        hist_nevents_beta0.Draw("hist")
        hist_nevents_beta1.Draw("same,hist")
        leg_nevents.Draw("same")
        if infilename_data != "":
            data_line_nevents.Draw("same")
        canv_nevents.Update()

        outfile_name_nevents = tag + "_nevents.pdf"
        outfile_name_nevents2 = tag + "_nevents.png"
        print "saving nevents image to file", outfile_name_nevents
        canv_nevents.SaveAs(outfile_name_nevents,"recreate")
        canv_nevents.SaveAs(outfile_name_nevents2,"recreate")

        # If we have comparison data, plot it too
        plot_comparison_both(hist_nevents_beta0, infilename_comp_nobs_beta0, 0, hist_nevents_beta1,
                             infilename_comp_nobs_beta1, 1, "rate-only", canv_nevents, leg_nevents,
                             tag, valor_comp, prior_predictive, label_comp)
        

        #
        # Check for any overflow points that won't have been plot
        #
        warn_over_under_flow(hist_nevents_beta0)
        warn_over_under_flow(hist_nevents_beta1)
        
        #
        # Calculate p-values
        #
        print "\n\n***** p-values *****"
        print "# nuebar data events =", nevents_data
        p_val = hist_pval_nolog_beta0.ProjectionX().GetMean()
        err_pval = err_pvalue(p_val, hist_pval_nolog_beta0, prior_predictive)
        print "Rate only (true beta = 0) p-value (right tail) =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"    # TODO: error on sigma?
        p_val = hist_pval_nolog_beta1.ProjectionX().GetMean()
        err_pval = err_pvalue(p_val, hist_pval_nolog_beta1, prior_predictive)
        print "Rate only (true beta = 1) p-value (left tail) =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"
        print "\n\n"


        #
        # Plot the p-value distributions
        # TODO: Test this method
        plot_pval_dist(hist_pval_beta0, 0, "log_"+tag)
        plot_pval_dist(hist_pval_beta1, 1, "log_"+tag)
        plot_pval_dist(hist_pval_nolog_beta0, 0, tag)
        plot_pval_dist(hist_pval_nolog_beta1, 1, tag)

        # Save outputs
        fout = ROOT.TFile(tag+"_rate-only_hists.root","recreate")
        hist_nevents_beta0.Write("rate-only_stat_beta0")
        hist_nevents_beta1.Write("rate-only_stat_beta1")
        canv_nevents.Write("rate-only_canvas")
        fout.Close()

    print "Done!"



def p_value(distribution, data, tail):
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


def err_pvalue(pval, distribution, prior_predictive):
    if prior_predictive:
        return ROOT.TMath.Sqrt(pval*(1-pval)/(distribution.Integral()))
    else:
        print "err_pvalue(): Posterior predictive errors not yet supported! Exiting."
        util.Die(1)


def sigma(p_val):
    return ROOT.TMath.Sqrt(2.) * ROOT.TMath.ErfInverse(1 - p_val)


def warn_over_under_flow(hist):
    overflows = hist.GetBinContent(hist.GetNbinsX()+1)
    underflows = hist.GetBinContent(0)
    text2 = "points. Consider extending the histogram range (p-val calculation won't be accurate)."
    if overflows != 0:
        print text1, "overflow", text2
    if underflows != 0:
        print text1, "underflow", text2


def plot_comparison_both(valor_hist1, infile1name_comp, beta_file1, valor_hist2, infile2name_comp, beta_file2, statistic, canv, leg, tag, valor_comp, prior_predictive, label_comp):
    outtag = ""
    label = label_comp + " (Fake-data #beta = "
    hname  = ""
    draw_opt = "same"
    if statistic == "rate-only":
        outtag = "_comparison_nevents"
        hname = "Nobs"
        #draw_opt = draw_opt + ",B"
        draw_opt = draw_opt + ",hist"
    elif statistic == "rate+shape":
        outtag = "_comparison_Deltachisq_logscale"
        hname = "Dchisq"
    else:
        print "Test statistic", statistic, "not understood. Exiting."
        util.Die(-1)

    comp_h1 = ROOT.TH1D()
    comp_h2 = ROOT.TH1D()
    if infile1name_comp != "":
        comp_h1 = plot_comparison_single(valor_hist1, infile1name_comp, beta_file1, statistic, leg, hname, label, valor_comp, prior_predictive)
        comp_h1.Draw(draw_opt)
    if infile2name_comp != "":
        comp_h2 = plot_comparison_single(valor_hist2, infile2name_comp, beta_file2, statistic, leg, hname, label, valor_comp, prior_predictive)
        comp_h2.Draw(draw_opt)
        
    if infile1name_comp != "" or infile2name_comp != "":
        if statistic == "rate+shape":
            leg.SetX1NDC(0.9*leg.GetX1NDC())
        leg.SetY1NDC(0.92*leg.GetY1NDC())
        canv.Update()
        ext1 = ".pdf"
        ext2 = ".png"
        outfile_name = tag + outtag + ext1
        outfile_name2 = tag + outtag + ext2
        print "saving comparison image to file", outfile_name
        canv.SaveAs(outfile_name,"recreate")
        canv.SaveAs(outfile_name2,"recreate")


def plot_comparison_single(valor_hist, infilename_comp, beta, statistic, leg, hname, label, valor_comp, prior_predictive):
    label = label + str(beta) + ")"
    col = 0
    style = 3
    if beta == 0:
        col = 40
    elif beta == 1:
        col = 41
    else:
        print "beta = ", beta, "not understood. Exiting."
        util.Die(-1)
    
    file_comp = ROOT.TFile(infilename_comp,"READ")
    hist_comp = 0
    if valor_comp == False or valor_comp == "False" or valor_comp == "":
        # If reading a digitized P-Theta histogram, directly get the hist
        hname = "h0"
        hist_comp = file_comp.Get(hname)
    else:
        # If comparing to another VALOR result, read the Dchisq or Nevents values into a hist
        if statistic == "rate-only":
            hist_comp = ROOT.TH1D("comprison_"+str(beta), "comparison_"+str(beta), high_nevents - low_nevents, low_nevents, high_nevents)
            hname = "Nevents_nuebar_weighted"
            if prior_predictive:
                hname = "Nevents_nuebar_unweighted"
            histfill_nevents = file_comp.Get(hname)
            bin_width = histfill_nevents.GetBinCenter(2) - histfill_nevents.GetBinCenter(1)
            for ibin in range(0, histfill_nevents.GetNbinsX()+2):           
                nev = histfill_nevents.GetBinCenter(ibin) - bin_width/2.
                nentries = histfill_nevents.GetBinContent(ibin)
                hist_comp.Fill(nev, nentries)
        else:
            hist_comp = ROOT.TH1D("comparison_"+str(beta), "comparison_"+str(beta), nbins_dchisq, low_dchisq, high_dchisq)
            tree_name = "fitnt"
            tree = file_comp.Get(tree_name)
            ievent = 0
            for event in tree:
                hist_comp.Fill(event.deltachisq01)
                if ievent % 10000 == 0 and not ievent == 0:
                    print "reading the beta =", str(beta), "comparison tree, event #", ievent
                ievent += 1
    area_norm(hist_comp, valor_hist)
    hist_comp.SetLineColor(col)
    hist_comp.SetLineStyle(2)

    print "\n\n ***** Compare RMS values and integral between the VALOR and P-Theta histograms *****"
    if statistic == "rate-only":
        # Draw as a bar chart
        #hist_comp.SetFillColor(col)
        ##hist_comp.SetFillStyle(3003)
        #hist_comp.SetBarWidth(0.1)
        #hist_comp.SetFillStyle(3003)
        #hist_comp.SetFillColor(col)
        #hist_comp.SetBarWidth(0.2)
        #if beta == 0:
            #hist_comp.SetBarOffset(0.2)
            #hist_comp.SetBarOffset(0.)
        #elif beta ==1:
            #hist_comp.SetBarOffset(0.3)
            #hist_comp.SetBarOffset(0.2)
        print "RMS value (axis 0) for rate-only distribution (beta =", str(beta), ", valor) =", valor_hist.GetRMS(0)
        print "RMS value (axis 0) for rate-only distribution (beta =", str(beta), ", P-theta) =", hist_comp.GetRMS(0)
        print "RMS value (axis 1) for rate-only distribution (beta =", str(beta), ", valor) =", valor_hist.GetRMS(1)
        print "RMS value (axis 1) for rate-only distribution (beta =", str(beta), ", valor) =", valor_hist.GetRMS(1)
        print "RMS value (axis 1) for rate-only distribution (beta =", str(beta), ", P-theta) =", hist_comp.GetRMS(1)
    elif statistic == "rate+shape":
        # draw as a hist
        hist_comp.SetLineStyle(style)
        print "RMS value (axis 0) for rate+shape distribution (beta =", str(beta), ", valor) =", valor_hist.GetRMS(0)
        print "RMS value (axis 0) for rate+shape distribution (beta =", str(beta), ", P-theta) =", hist_comp.GetRMS(0)
        print "RMS value (axis 1) for rate+shape distribution (beta =", str(beta), ", valor) =", valor_hist.GetRMS(1)
        print "RMS value (axis 1) for rate+shape distribution (beta =", str(beta), ", P-theta) =", hist_comp.GetRMS(1)

    print "valor_Hist integral =", valor_hist.Integral(), "; comp_hist integral =",  hist_comp.Integral()
    print "\n\n"
    new_hist = copy.deepcopy(hist_comp)
    leg.AddEntry(new_hist,label,"l")
    return new_hist

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


def plot_pval_dist(hist, beta, tag):
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
    print "saving p-val image to file", outfile_name
    canv.SaveAs(outfile_name, "recreate")
    canv.SaveAs(outfile_name2, "recreate")


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
    parser.add_option('--data_label', dest='data_label', default="data", help='Label to be put in the legend entry corresponding to the data')
    parser.add_option('--tag', dest='tag', default="nuebarapp_deltachisq", help='Tag used as the output file name')
    parser.add_option('--label', dest='label', default="Events", help='Label used for the legend')
    parser.add_option('--label-comp', dest='label_comp', default="P-theta Events", help='Label used for the legend for comparison hists')
    parser.add_option('--test-statistic', dest='test_statistic', default="rate+shape", help='The test statistic used in the marginalisation and fit')
    parser.add_option('--prior-predictive', dest='prior_predictive', action='store_true', default="False", help='Include this flag to make the rate-only analysis prior-predictive (for rate+shape, the input file determines this)')
    parser.add_option('--valor-comp', dest='valor_comp', action='store_true', default="False", help='Include this flag to compare to another VALOR file isntead of the default P-Theta file')


    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o, other = parser.parse_args()

    main(o.infilename_beta0, o.infilename_beta1, o.tag, o.infilename_data, o.data_label,
         o.infilename_comp_nobs_beta0, o.infilename_comp_nobs_beta1,
         o.infilename_comp_Dchisq_beta0, o.infilename_comp_Dchisq_beta1, o.test_statistic, o.prior_predictive, o.valor_comp, o.label, o.label_comp)

    if other:
        print "Warning - ignored arguments: ",other
