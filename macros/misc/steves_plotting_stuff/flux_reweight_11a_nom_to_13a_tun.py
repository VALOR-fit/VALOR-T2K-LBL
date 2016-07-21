#!/usr/bin/env python

import ROOT
import os
import math

R = ROOT

for suff in ["","_fine"]:

  # Laura's files - get 13a tuning from 13a nominal.
  input_files_13atun = dict()
  input_files_13atun["FHC"] = os.environ["T2KOSCANA"]+"inputs/official/2014/flux/sk_tuned13av1.0_13anom_run1-5b"+suff+".root"
  input_files_13atun["RHC"] = os.environ["T2KOSCANA"]+"inputs/official/2014/flux/sk_tuned13av1.0_13anom_run5c_antinumode"+suff+".root"

  # Sukuki-san's files - get 11bv32 from 11a.
  input_files_11anom = dict()
  input_files_11anom["FHC"] = os.environ["T2KOSCANA"]+"inputs/official/2013/flux/11b_v3.2_run1-4a/sk_tuned11bv3.2_11anom_run1-run4"+suff+".root"
  input_files_11anom["RHC"] = os.environ["T2KOSCANA"]+"inputs/official/2013/flux/11b_v3.2_antinu/sk_tuned11bv3.2_11anom_m250ka"+suff+".root"

  species = ["numu","numub","nue","nueb"]

  for beam_mode in ["FHC","RHC"]:
    infile_13atun = R.TFile(input_files_13atun[beam_mode])
    infile_11anom = R.TFile(input_files_11anom[beam_mode])
    
    assert(infile_13atun)
    assert(infile_11anom)
    
    outfile = R.TFile(os.environ["T2KOSCANA"]+"inputs/flux_reweight/reweight_13aMC_to_13aFlux_"+beam_mode+suff+".root","RECREATE")
    for s in species:
      
      tun_hname = "enu_sk_tuned13a_"+s
      tun_hist = infile_13atun.Get(tun_hname)
      assert(tun_hist)
          
      nom_hname = "enu_sk_11a_norm_"+s
      nom_hist = infile_11anom.Get(nom_hname)
      assert(nom_hist)
      
      out_hist = R.TH1D(tun_hist)
      out_hist.Reset()
      out_hist.GetSumw2().Set(0)
      out_hist.SetTitle("13b tuned flux / 11a nominal flux")
      out_hist.SetName(tun_hname+"_11a_ratio")
      
      for ibin in xrange(1,out_hist.GetXaxis().GetNbins()+1):
        center = tun_hist.GetBinCenter(ibin)
        if center > 10:
          out_hist.Fill(center,1)
          out_hist.SetBinError(ibin,0)
          continue
          
        assert(abs(center-nom_hist.GetBinCenter(ibin))<1E-10)
        
        tun_val = tun_hist.GetBinContent(ibin)
        tun_err = tun_hist.GetBinError  (ibin)
        nom_val = nom_hist.GetBinContent(ibin)
        nom_err = nom_hist.GetBinError  (ibin)
        
        wght_11a_nom_to_13a_tun =  tun_val / nom_val
        frac_error = math.sqrt( (tun_err/tun_val)**2 + (nom_err/nom_val)**2 )
        
        out_hist.Fill(center, wght_11a_nom_to_13a_tun)
        out_hist.SetBinError(ibin,frac_error * wght_11a_nom_to_13a_tun)
      
      outfile.cd()
      out_hist.Write()

print "Done"
