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
nbins_dchisq    = 250 #700
low_dchisq      = -20. #-20
high_dchisq     = 45.  #50
low_nevents     = -1   #-1
high_nevents    = 60   #50
scale_max_rate  = 5
scale_max_shape = 2
leg_text_size   = 0.02 #0.03
line_thickness  = 2 # 2 is good for projectors?
draw_prelim     = 1
prelim_text     = "T2K Run 1-9 Preliminary"

def main(infilename_beta0, infilename_beta1, tag, infilename_data, data_label, infilename_comp_nobs_beta0, infilename_comp_nobs_beta1, infilename_comp_Dchisq_beta0, infilename_comp_Dchisq_beta1, test_statistic, prior_predictive, valor_comp, label, label_comp, read_nevents, start_nevents, infilename_comp_data, nev_comp, Dchisq_comp):
    read_nevents += start_nevents
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

        # For the data
        hDataErecoToy_1Rmu = 0
        hDataErecoToy_1Re = 0
        hDataErecoToy_RHC1Rmu = 0
        hDataErecoToy_MultiRe = 0
        if infilename_data != "":
            for event in tree_data:
                Delta_chisq_data = event.deltachisq01   # Reading trees in pyroot is weird, but there'll be a better way to do this than a for loop with only one iteration...
                print "# nuebar data Dchisq =", Delta_chisq_data
                hDataErecoToy_1Rmu = event.hErecoToy_1Rmu
                hDataErecoToy_1Re = event.hErecoToy_1Re
                hDataErecoToy_RHC1Rmu = event.hErecoToy_RHC1Rmu
                hDataErecoToy_MultiRe = event.hErecoToy_MultiRe

        # Fill with Delta chisq = chisq(beta = 0 for marginalisation toys) - chisq(beta = 1 for marginalisation toys) ...
        # For distribution toys generated with beta = 0
        nevents = tree_beta0.GetEntries()
        wghts_beta0 = []
        tail_beta0 = []
        ievent = 0
        for event in tree_beta0:
            if ievent < start_nevents:
                ievent += 1
                continue
            if read_nevents > 0 and ievent >= read_nevents:
                break
            # If we're using re-weighting, then we'll need the control sample chisq to calculate the weights
            chisq = 0
            chisq += util.Chisq_Penalty_TH2D(hDataErecoToy_1Rmu, event.hErecoExp_1Rmu)
            chisq += util.Chisq_Penalty_TH2D(hDataErecoToy_1Re, event.hErecoExp_1Re)
            chisq += util.Chisq_Penalty_TH2D(hDataErecoToy_RHC1Rmu, event.hErecoExp_RHC1Rmu)
            chisq += util.Chisq_Penalty_TH2D(hDataErecoToy_MultiRe, event.hErecoExp_MultiRe)
            wghts_beta0.append(ROOT.TMath.Exp(-0.5*chisq))  # No need to divide by exp(-0.5*chisq_min) as we do when generating toys, as this factor is canceled when calculating p-values their errors
            if prior_predictive == True or prior_predictive == "True":
                hist_deltachisq01_beta0.Fill(event.deltachisq01)
            else:
                hist_deltachisq01_beta0.Fill(event.deltachisq01, ROOT.TMath.Exp(-0.5*chisq))
                print "***test: Filling beta = 0 hist with", event.deltachisq01, ",", ROOT.TMath.Exp(-0.5*chisq)
                print "event.hErecoExp_1Rmu.Integral() =", event.hErecoExp_1Rmu.Integral(), "; data =", hDataErecoToy_1Rmu.Integral()
                print "event.hErecoExp_1Re.Integral() =", event.hErecoExp_1Re.Integral(), "; data =", hDataErecoToy_1Re.Integral()
                print "event.hErecoExp_RHC1Rmu.Integral() =", event.hErecoExp_RHC1Rmu.Integral(), "; data =", hDataErecoToy_RHC1Rmu.Integral()
                print "event.hErecoExp_MultiRe.Integral() =", event.hErecoExp_MultiRe.Integral(), "; data =", hDataErecoToy_MultiRe.Integral(), "\n\n"
            if event.deltachisq01 <= Delta_chisq_data:
                tail_beta0.append("left")
            # Not elif, because if deltachisq01 = Delta_chisq_data, then this value is used in both tails    
            if event.deltachisq01 >= Delta_chisq_data:    
                tail_beta0.append("right")
            if ievent % 1000 == 0:
                print "reading the beta = 0 tree, event #", ievent
            ievent += 1
        if len(wghts_beta0) != len(tail_beta0):
            print "length of wghts_beta0 =", len(wghts_beta0)
            print "length of tail_beta0 =", len(tail_beta0)
            print "These should be equal. Exiting."
            util.Die(1)

        # And for distribution toys generated with beta = 1
        ievent = 0
        wghts_beta1 = []
        tail_beta1 = []
        for event in tree_beta1:
            if ievent < start_nevents:
                ievent += 1
                continue
            if read_nevents > 0 and ievent >= read_nevents:
                break
            chisq = 0
            chisq += util.Chisq_Penalty_TH2D(hDataErecoToy_1Rmu, event.hErecoExp_1Rmu)
            chisq += util.Chisq_Penalty_TH2D(hDataErecoToy_1Re, event.hErecoExp_1Re)
            chisq += util.Chisq_Penalty_TH2D(hDataErecoToy_RHC1Rmu, event.hErecoExp_RHC1Rmu)
            chisq += util.Chisq_Penalty_TH2D(hDataErecoToy_MultiRe, event.hErecoExp_MultiRe)
            wghts_beta1.append(ROOT.TMath.Exp(-0.5*chisq))  # No need to divide by exp(-0.5*chisq_min) as we do when generating toys, as this factor is canceled when calculating p-values their errors
            if prior_predictive == True or prior_predictive == "True":
                hist_deltachisq01_beta1.Fill(event.deltachisq01)
            else:
                hist_deltachisq01_beta1.Fill(event.deltachisq01, ROOT.TMath.Exp(-0.5*chisq))
            if event.deltachisq01 <= Delta_chisq_data:
                tail_beta1.append("left")
            # Not elif, because if deltachisq01 = Delta_chisq_data, then this value is used in both tails
            if event.deltachisq01 >= Delta_chisq_data:
                tail_beta1.append("right")
            if ievent % 1000 == 0:
                print "reading the beta = 1 tree, event #", ievent
            ievent += 1
        if len(wghts_beta1) != len(tail_beta1):
            print "length of wghts_beta1 =", len(wghts_beta1)
            print "length of tail_beta1 =", len(tail_beta1)
            print "These should be equal. Exiting."
            util.Die(1)

        #
        # Plot the Delta chisq histograms
        #
        hist_deltachisq01_beta0.SetLineColor(8)
        hist_deltachisq01_beta0.SetLineWidth(line_thickness)
        hist_deltachisq01_beta1.SetLineColor(2)
        hist_deltachisq01_beta1.SetLineWidth(line_thickness)
        hist_deltachisq01_beta0.GetXaxis().SetTitle("#Delta #chi^{2} = #chi^{2}(#beta=0) - #chi^{2}(#beta=1)")
        hist_deltachisq01_beta0.GetYaxis().SetTitle("Entries")
        maximum_dchisq = 1.1 * max(hist_deltachisq01_beta0.GetMaximum(), hist_deltachisq01_beta1.GetMaximum())
        hist_deltachisq01_beta0.SetMaximum(maximum_dchisq)
        data_line_dchisq = ROOT.TLine(Delta_chisq_data, 0, Delta_chisq_data, maximum_dchisq)
        data_line_dchisq.SetLineStyle(2)
        data_line_dchisq.SetLineWidth(line_thickness)

        leg_ylow = 0.8
        if infilename_data != "":
            leg_ylow = 0.83
        leg = ROOT.TLegend(0.72,leg_ylow,0.9,0.95)
        leg.SetTextSize(leg_text_size)
        leg.AddEntry(hist_deltachisq01_beta0, label + " (#beta = 0)","l")
        leg.AddEntry(hist_deltachisq01_beta1, label + " (#beta = 1)","l")
        if infilename_data != "":
            data_line_label = data_label
            if Dchisq_comp != -999:
                data_line_label = data_label + ", " + label
            leg.AddEntry(data_line_dchisq,"#Delta #chi^{2} ("+data_line_label+") = "+str("{0:.2f}".format(Delta_chisq_data)),"l")

        canv = ROOT.TCanvas("Delta chisq", "#Delta #chi^{2}", 800, 600)
        if prior_predictive == True or prior_predictive == "True":
            hist_deltachisq01_beta0.SetMinimum(0.5)
        else:
            hist_deltachisq01_beta0.SetMinimum(hist_deltachisq01_beta0.GetMinimum()*0.9)
        hist_deltachisq01_beta0.Draw()
        hist_deltachisq01_beta1.Draw("same")
        if infilename_data != "":
            data_line_dchisq.Draw("same")
        leg.Draw("same")
        if draw_prelim == 1:
            prelim = ROOT.TLatex()
            prelim.SetTextSize(0.05)
            prelim.DrawLatexNDC(0.5, 0.96, prelim_text)
        canv.Update()

        outfile_name = tag + "_Deltachisq.pdf"
        outfile_name2 = tag + "_Deltachisq.png"
        outfile_name3 = tag + "_Deltachisq.eps"
        print "saving Delta chisq image to file", outfile_name
        canv.SaveAs(outfile_name,"recreate")
        canv.SaveAs(outfile_name2,"recreate")
        canv.SaveAs(outfile_name3,"recreate")

        # If we have comparison data, plot it too    TODO: test: does this disrupt the logscale call to plot_comparison_both?
        if Dchisq_comp != -999:
            data_line_comp = ROOT.TLine(Dchisq_comp, 0, Dchisq_comp, maximum_dchisq)
            data_line_comp.SetLineStyle(5)
            data_line_comp.SetLineWidth(line_thickness)
            data_line_comp.Draw("same")
            leg.AddEntry(data_line_comp,"#Delta #chi^{2} ("+data_label+", "+label_comp+") = "+str("{0:.2f}".format(Dchisq_comp)),"l")
        plot_comparison_both(hist_deltachisq01_beta0, infilename_comp_Dchisq_beta0, 0, hist_deltachisq01_beta1,
                             infilename_comp_Dchisq_beta1, 1, "rate+shape", canv, leg, tag, valor_comp,
                             prior_predictive, label_comp, read_nevents, infilename_comp_data)

        #
        # Plot the Delta chisq histograms on a log scale
        #
        hist_deltachisq01_beta0.SetMaximum(scale_max_shape*maximum_dchisq)
        data_line_dchisq.SetY2(scale_max_shape*maximum_dchisq)
        if Dchisq_comp != -999:
            data_line_comp.SetY2(scale_max_shape*maximum_dchisq)
        canv.SetLogy(1)
        canv.Update()
        outfile_name = tag + "_Deltachisq_logscale.pdf"
        outfile_name2 = tag + "_Deltachisq_logscale.png"
        outfile_name3 = tag + "_Deltachisq_logscale.eps"
        print "saving Delta chisq log scale image to file", outfile_name
        canv.SaveAs(outfile_name,"recreate")
        canv.SaveAs(outfile_name2,"recreate")
        canv.SaveAs(outfile_name3,"recreate")

        # If we have comparison data, plot it too
        plot_comparison_both(hist_deltachisq01_beta0, infilename_comp_Dchisq_beta0, 0, hist_deltachisq01_beta1,
                             infilename_comp_Dchisq_beta1, 1, "rate+shape", canv, leg, tag+"_logscale", valor_comp,
                             prior_predictive, label_comp, read_nevents, infilename_comp_data)
        
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
        p_val = -999
        err_pval = -999
        if prior_predictive == True or prior_predictive == "True":
            p_val = p_value(hist_deltachisq01_beta0, Delta_chisq_data, "right")
            err_pval = err_pvalue(p_val, hist_deltachisq01_beta0)
        else:    
            p_val = p_value_reweighting("right", wghts_beta0, tail_beta0)
            err_pval = err_pvalue_reweighting(p_val, wghts_beta0, tail_beta0)
        print "Rate+shape (true beta = 0) p-value (right tail) =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"
        if prior_predictive == True or prior_predictive == "True":
            p_val = p_value(hist_deltachisq01_beta1, Delta_chisq_data, "left")
            err_pval = err_pvalue(p_val, hist_deltachisq01_beta1)
        else:
            p_val = p_value_reweighting("left", wghts_beta1, tails_beta1)
            err_pval = err_pvalue_reweighting(p_val, wghts_beta1, tails_beta1)
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
        if prior_predictive == True or prior_predictive == "True":
            hname = "Nevents_nuebar_unweighted"
        print "*** test: hname main dist =", hname
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
            if ievent < start_nevents:
                ievent += 1
                continue
            if read_nevents > 0 and ievent >= read_nevents:
                break
            if prior_predictive == True or prior_predictive == "True":
                hist_pval_beta0.Fill(ROOT.TMath.Log10(event.pvalue_right_unweighted), event.Eventrate_fakedata)
                hist_pval_nolog_beta0.Fill(event.pvalue_right_unweighted, event.Eventrate_fakedata)
            else:
                hist_pval_beta0.Fill(ROOT.TMath.Log10(event.pvalue_right_weighted), event.Eventrate_fakedata)
                hist_pval_nolog_beta0.Fill(event.pvalue_right_weighted, event.Eventrate_fakedata)
            if ievent % 1000 == 0 and not ievent == 0:
                print "reading the beta = 0 tree, event #", ievent
            ievent += 1
        ievent = 0
        for event in tree_beta1:
            if ievent < start_nevents:
                ievent += 1
                continue
            if read_nevents > 0 and ievent >= read_nevents:
                break
            if prior_predictive == True or prior_predictive == "True":
                hist_pval_beta1.Fill(ROOT.TMath.Log10(event.pvalue_left_unweighted), event.Eventrate_fakedata)
                hist_pval_nolog_beta1.Fill(event.pvalue_left_unweighted, event.Eventrate_fakedata)
            else:
                hist_pval_beta1.Fill(ROOT.TMath.Log10(event.pvalue_left_weighted), event.Eventrate_fakedata)
                hist_pval_nolog_beta1.Fill(event.pvalue_left_weighted, event.Eventrate_fakedata)
            if ievent % 1000 == 0 and not ievent == 0:
                print "reading the beta = 1 tree, event #", ievent
            ievent += 1

        #
        # Plot the nevents histograms
        #
        hist_nevents_beta0.SetLineColor(8)
        hist_nevents_beta0.SetLineWidth(line_thickness)
        hist_nevents_beta1.SetLineColor(2)
        hist_nevents_beta1.SetLineWidth(line_thickness)
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
        hist_nevents_beta0.GetXaxis().SetTitle("Number of observed RHC 1Re-like events")
        hist_nevents_beta0.GetYaxis().SetTitle("Entries")
        maximum_nevents = 1.1 * max(hist_nevents_beta0.GetMaximum(), hist_nevents_beta1.GetMaximum())
        hist_nevents_beta0.SetMaximum(maximum_nevents)
        data_line_nevents = ROOT.TLine(nevents_data, 0, nevents_data, maximum_nevents)
        data_line_nevents.SetLineStyle(2)
        data_line_nevents.SetLineWidth(line_thickness)

        leg_nevents = ROOT.TLegend(0.57,0.85,0.9,0.95)
        leg_nevents.SetTextSize(leg_text_size)
        leg_nevents.AddEntry(hist_nevents_beta0, label + " (#beta = 0)","l")
        leg_nevents.AddEntry(hist_nevents_beta1, label + " (#beta = 1)","l")
        if infilename_data != "":
            data_line_label = data_label
            if nev_comp != -999:
                data_line_label = data_label + ", " + label
            leg_nevents.AddEntry(data_line_nevents,"Events ("+data_line_label+") = "+str("{0:.2f}".format(nevents_data)),"l")

        canv_nevents = ROOT.TCanvas("nevents", "## events", 800, 600)
        #hist_nevents_beta0.Draw("B")
        #hist_nevents_beta1.Draw("same,B")
        hist_nevents_beta0.SetMinimum(0.5)
        hist_nevents_beta0.Draw("hist")
        hist_nevents_beta1.Draw("same,hist")
        leg_nevents.Draw("same")
        if infilename_data != "":
            data_line_nevents.Draw("same")
        if draw_prelim == 1:
            prelim = ROOT.TLatex()
            prelim.SetTextSize(0.05)
            prelim.DrawLatexNDC(0.5, 0.96, prelim_text)
        canv_nevents.Update()

        outfile_name_nevents = tag + "_nevents.pdf"
        outfile_name_nevents2 = tag + "_nevents.png"
        outfile_name_nevents3 = tag + "_nevents.eps"
        print "saving nevents image to file", outfile_name_nevents
        canv_nevents.SaveAs(outfile_name_nevents,"recreate")
        canv_nevents.SaveAs(outfile_name_nevents2,"recreate")
        canv_nevents.SaveAs(outfile_name_nevents3,"recreate")

        # If we have comparison data, plot it too     TODO: test: does this disrupt the logscale call to plot_comparison_both?
        if nev_comp != -999:
            data_line_comp = ROOT.TLine(nev_comp, 0, nev_comp, maximum_nevents)
            data_line_comp.SetLineStyle(5)
            data_line_comp.SetLineWidth(line_thickness)
            data_line_comp.Draw("same")
            leg_nevents.AddEntry(data_line_comp,"Events ("+data_label+", "+label_comp+") = "+str("{0:.2f}".format(nev_comp)),"l")
        plot_comparison_both(hist_nevents_beta0, infilename_comp_nobs_beta0, 0, hist_nevents_beta1,
                             infilename_comp_nobs_beta1, 1, "rate-only", canv_nevents, leg_nevents,
                             tag, valor_comp, prior_predictive, label_comp, read_nevents, infilename_comp_data)

        # Also plot on a log scale
        maximum_nevents = scale_max_rate * max(hist_nevents_beta0.GetMaximum(), hist_nevents_beta1.GetMaximum())
        hist_nevents_beta0.SetMaximum(maximum_nevents)
        data_line_nevents.SetY2(maximum_nevents)
        data_line_comp.SetY2(maximum_nevents)
        canv_nevents.SetLogy(1)
        canv_nevents.Update()
        outfile_name_nevents = tag + "_nevents_logscale.pdf"
        outfile_name_nevents2 = tag + "_nevents_logscale.png"
        outfile_name_nevents3 = tag + "_nevents_logscale.eps"
        print "saving logscale nevents image to file", outfile_name_nevents
        canv_nevents.SaveAs(outfile_name_nevents,"recreate")
        canv_nevents.SaveAs(outfile_name_nevents2,"recreate")
        canv_nevents.SaveAs(outfile_name_nevents3,"recreate")


        # If we have comparison data, plot it too
        plot_comparison_both(hist_nevents_beta0, infilename_comp_nobs_beta0, 0, hist_nevents_beta1,
                             infilename_comp_nobs_beta1, 1, "rate-only", canv_nevents, leg_nevents,
                             tag+"_logscale", valor_comp, prior_predictive, label_comp, read_nevents, infilename_comp_data)
        

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
        err_pval = err_pvalue(p_val, hist_pval_nolog_beta0)
        print "Rate only (true beta = 0) p-value (right tail) =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"    # TODO: error on sigma?
        p_val = hist_pval_nolog_beta1.ProjectionX().GetMean()
        err_pval = err_pvalue(p_val, hist_pval_nolog_beta1)
        print "Rate only (true beta = 1) p-value (left tail) =", p_val, "+/-", err_pval, "=", sigma(p_val), "sigma"
        print "\n\n"


        #
        # Plot the p-value distributions
        # TODO: Test this method
        plot_pval_dist(hist_pval_beta0, 0, "log_"+tag)
        plot_pval_dist(hist_pval_beta1, 1, "log_"+tag)
        plot_pval_dist(hist_pval_nolog_beta0, 0, tag)
        plot_pval_dist(hist_pval_nolog_beta1, 1, tag)
        print "For beta = 0, the fraction of toys with pval < 0.1 =",  count_pvals_below_threshold(hist_pval_nolog_beta0.ProjectionX(), 0.1)
        print "For beta = 1, the fraction of toys with pval < 0.1 =",  count_pvals_below_threshold(hist_pval_nolog_beta1.ProjectionX(), 0.1)
        print "For beta = 0, the fraction of toys with pval < 0.05 =", count_pvals_below_threshold(hist_pval_nolog_beta0.ProjectionX(), 0.05)
        print "For beta = 1, the fraction of toys with pval < 0.05 =", count_pvals_below_threshold(hist_pval_nolog_beta1.ProjectionX(), 0.05)

        # Save outputs
        fout = ROOT.TFile(tag+"_rate-only_hists.root","recreate")
        hist_nevents_beta0.Write("rate-only_stat_beta0")
        hist_nevents_beta1.Write("rate-only_stat_beta1")
        canv_nevents.Write("rate-only_canvas")
        fout.Close()

    print "Done!"


# p-value calculation when re-weighting is not used
def p_value(distribution, data, tail):
    print "Calculating p-value under the assumption that MC sampling was used to generate toys"
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


# p-value calculation when re-weighting is used (so we calculate with the weights)
def p_value_reweighting(tail, wghts, tails):
    print "Calculating p-value under the assumption that re-weighting was used to generate toys"
    sum_wghts_tail = 0
    for i in range(0, len(wghts)):
        if tails[i] == tail:
            sum_wghts_tail += wghts[i]
    return sum_wghts_tail / sum(wghts)        


# p-value error calculation when re-weighting is not used
def err_pvalue(pval, distribution):
    print "Calculating p-value errors under the assumption that MC sampling was used to generate toys"
    return ROOT.TMath.Sqrt(pval*(1-pval)/(distribution.Integral()))


# p-value error calculation when re-weighting is used (so we calculate with the weights)
def err_pvalue_reweighting(pval, wghts, tails):
    print "Calculating p-value errors under the assumption that re-weighting was used to generate toys"
    sum_wghts = sum(wghts)
    sumsq_wghts_left = 0
    sumsq_wghts_right = 0
    for i in range(0, len(wghts)):
        if tails[i] == "left":
            sumsq_wghts_left += wghts[i] * wghts[i]
        else:
            sumsq_wghts_right += wghts[i] * wghts[i]
    return TMath.Sqrt(TMath.Power((1 - pval)/sum_wghts,2) * sumsq_wghts_right + TMath.Power(pval/sum_wghts,2) * sumsq_wghts_left)        


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


def plot_comparison_both(valor_hist1, infile1name_comp, beta_file1, valor_hist2, infile2name_comp, beta_file2, statistic, canv, leg, tag, valor_comp, prior_predictive, label_comp, read_nevents, infilename_comp_data):
    outtag = ""
    label = label_comp + " (#beta = "
    hname  = ""
    draw_opt = "same"
    if statistic == "rate-only":
        outtag = "_comparison_nevents"
        hname = "Nobs"
        #draw_opt = draw_opt + ",B"
        draw_opt = draw_opt + ",hist"
    elif statistic == "rate+shape":
        outtag = "_comparison_Deltachisq"
        hname = "Dchisq"
    else:
        print "Test statistic", statistic, "not understood. Exiting."
        util.Die(-1)

    if infilename_comp_data != "" and not (valor_comp == False or valor_comp == "False" or valor_comp == ""):
        data_comp = -999
        file_comp_data = ROOT.TFile(infilename_comp_data,"READ")
        tree_name = "fitnt"
        tree_data = ROOT.TTree()
        tree_data  = file_comp_data.Get(tree_name)
        if tree_data.GetEntries() != 1:
            print "Input comparison data tree has", tree_data.GetEntries(), "Expected one entry. Exiting."
            util.Die(-1)
        for event in tree_data:
            if statistic == "rate-only":
                data_comp = event.hErecoToy_RHC1Re.Integral()   # There'll be a better way to do this than a for loop with only one iteration...
            else:
                data_comp = event.deltachisq01
        data_line = ROOT.TLine(data_comp, 0, data_comp, valor_hist1.GetMaximum())
        data_line.SetLineStyle(6)
        data_line.SetLineWidth(line_thickness)
        data_type = "#Delta #chi^{2}"
        if statistic == "rate-only":
            data_type = "Events"
        leg.AddEntry(data_line,data_type+" (Asimov A) = "+str("{0:.2f}".format(data_comp)),"l")
        data_line.Draw("same")

    comp_h1 = ROOT.TH1D()
    comp_h2 = ROOT.TH1D()
    if infile1name_comp != "":
        comp_h1 = plot_comparison_single(valor_hist1, infile1name_comp, beta_file1, statistic, leg, hname, label, valor_comp, prior_predictive, read_nevents)
        comp_h1.Draw(draw_opt)
    if infile2name_comp != "":
        comp_h2 = plot_comparison_single(valor_hist2, infile2name_comp, beta_file2, statistic, leg, hname, label, valor_comp, prior_predictive, read_nevents)
        comp_h2.Draw(draw_opt)
        
    if infile1name_comp != "" or infile2name_comp != "":
        if statistic == "rate+shape":
            leg.SetX1NDC(0.9*leg.GetX1NDC())
        leg.SetY1NDC(0.92*leg.GetY1NDC())
        canv.Update()
        ext1 = ".pdf"
        ext2 = ".png"
        ext3 = ".eps"
        outfile_name = tag + outtag + ext1
        outfile_name2 = tag + outtag + ext2
        outfile_name3 = tag + outtag + ext3
        print "saving comparison image to file", outfile_name
        canv.SaveAs(outfile_name,"recreate")
        canv.SaveAs(outfile_name2,"recreate")
        canv.SaveAs(outfile_name3,"recreate")


def plot_comparison_single(valor_hist, infilename_comp, beta, statistic, leg, hname, label, valor_comp, prior_predictive, read_nevents):
    label = label + str(beta) + ")"
    col = 0
    style = 3
    if beta == 0:
        col = ROOT.TColor.kMagenta+1 #40   # colour 40 and 41 look good on monitors, but don't show up at all on projectors
    elif beta == 1:
        col = 4 #41
    else:
        print "beta = ", beta, "not understood. Exiting."
        util.Die(-1)
    
    file_comp = ROOT.TFile(infilename_comp,"READ")
    hist_comp = 0
    if statistic == "rate-only":
        hist_comp = ROOT.TH1D("comprison_"+str(beta), "comparison_"+str(beta), high_nevents - low_nevents, low_nevents, high_nevents)
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
        if statistic == "rate-only":
            hname = "Nevents_nuebar_unweighted"
            #hname = "Nevents_nuebar_weighted"
            #if prior_predictive == True or prior_predictive == "True":
            #    hname = "Nevents_nuebar_unweighted"
            #print "*** test: hname comparison =", hname    
            histfill_nevents = file_comp.Get(hname)
            bin_width = histfill_nevents.GetBinCenter(2) - histfill_nevents.GetBinCenter(1)
            for ibin in range(0, histfill_nevents.GetNbinsX()+2):           
                nev = histfill_nevents.GetBinCenter(ibin) - bin_width/2.
                nentries = histfill_nevents.GetBinContent(ibin)
                hist_comp.Fill(nev, nentries)
        else:
            tree_name = "fitnt"
            tree = file_comp.Get(tree_name)
            ievent = 0
            for event in tree:
                if read_nevents > 0 and ievent >= read_nevents and statistic == "rate+shape":
                    break
                hist_comp.Fill(event.deltachisq01)
                if ievent % 10000 == 0:
                    print "reading the beta =", str(beta), "comparison tree, event #", ievent
                ievent += 1
    area_norm(hist_comp, valor_hist)
    hist_comp.SetLineColor(col)
    hist_comp.SetLineWidth(line_thickness)
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
        print "RMS value (axis 1) for rate-only distribution (beta =", str(beta), ", valor) =", valor_hist.GetRMS(1)
        print "RMS value (axis 1) for rate-only distribution (beta =", str(beta), ", P-theta) =", hist_comp.GetRMS(1)
        print "mean value (axis 1) for rate-only distribution (beta =", str(beta), ", valor) =", valor_hist.GetMean(1)
        print "mean value (axis 1) for rate-only distribution (beta =", str(beta), ", P-theta) =", hist_comp.GetMean(1)
        print "peak value for rate-only distribution (beta =", str(beta), ", valor) =", valor_hist.GetBinCenter(valor_hist.GetMaximumBin())
        print "peak value for rate-only distribution (beta =", str(beta), ", P-theta) =", hist_comp.GetBinCenter(hist_comp.GetMaximumBin())
    elif statistic == "rate+shape":
        # draw as a hist
        hist_comp.SetLineStyle(style)
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
        if label == entry.GetLabel():
            add_label = False
            break
        entry = leg_entries.After(entry)
    if add_label:    
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
    outfile_name3 = tag + tag_beta + "_nevents.eps"
    print "saving p-val image to file", outfile_name
    canv.SaveAs(outfile_name, "recreate")
    canv.SaveAs(outfile_name2, "recreate")
    canv.SaveAs(outfile_name3, "recreate")


def count_pvals_below_threshold(hist, threshold):
    count = 0
    threshold_bin = hist.FindBin(threshold) - 1 # -1 because we want pvals < threshold
    frac = hist.Integral(1, threshold_bin) / hist.Integral(1, hist.GetNbinsX())
    return frac


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
    parser.add_option('--infilename_comp_data', dest='infilename_comp_data', default="", help='Root file output by marginalisation of the comparison data')
    parser.add_option('--data_label', dest='data_label', default="data", help='Label to be put in the legend entry corresponding to the data')
    parser.add_option('--tag', dest='tag', default="nuebarapp_deltachisq", help='Tag used as the output file name')
    parser.add_option('--label', dest='label', default="Fake-data", help='Label used for the legend')
    parser.add_option('--label-comp', dest='label_comp', default="P-theta Events", help='Label used for the legend for comparison hists')
    parser.add_option('--nev-comp', dest='nev_comp', default=-999, type=float, help='Data or Asimov number of RHC 1Re-like events from the results we\'re comparing to')
    parser.add_option('--Dchisq-comp', dest='Dchisq_comp', default=-999, type=float, help='Data or Asimov Delta chisq from the results we\'re comparing to')
    parser.add_option('--test-statistic', dest='test_statistic', default="rate+shape", help='The test statistic used in the marginalisation and fit')
    parser.add_option('--prior-predictive', dest='prior_predictive', action='store_true', default="False", help='Include this flag to make the rate-only analysis prior-predictive (for rate+shape, the input file determines this)')
    parser.add_option('--valor-comp', dest='valor_comp', action='store_true', default="False", help='Include this flag to compare to another VALOR file isntead of the default P-Theta file')
    parser.add_option('--read-nevents', dest='read_nevents', default=-1, type=int, help='Read only the first n events (staring from --start-nevents) from the input trees. This works for the p-values of both analyses, but the rate-only graph will be plot will all events. Set to <= 0 to read all events (default)')
    parser.add_option('--start-nevents', dest='start_nevents', default=0, type=int, help='Start reading the input trees from this event number')


    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o, other = parser.parse_args()

    main(o.infilename_beta0, o.infilename_beta1, o.tag, o.infilename_data, o.data_label,
         o.infilename_comp_nobs_beta0, o.infilename_comp_nobs_beta1,
         o.infilename_comp_Dchisq_beta0, o.infilename_comp_Dchisq_beta1, o.test_statistic, 
         o.prior_predictive, o.valor_comp, o.label, o.label_comp, o.read_nevents, o.start_nevents, o.infilename_comp_data, o.nev_comp, o.Dchisq_comp)

    if other:
        print "Warning - ignored arguments: ",other
