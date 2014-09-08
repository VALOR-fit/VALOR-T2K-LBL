
import ROOT
import os

R = ROOT

# Megan's files - get 13a nominal from 11bv32.
input_files_13anom = dict()
input_files_13anom["FHC"] = os.environ["T2KOSCANA"]+"inputs/official/2014/flux/sk_11btune_13anom_250ka_ratios.root"
input_files_13anom["RHC"] = os.environ["T2KOSCANA"]+"inputs/official/2014/flux/sk_11btune_13anom_m250ka_ratios.root"

# Laura's files - get 13a tuning from 13a nominal.
input_files_13atun = dict()
input_files_13atun["FHC"] = os.environ["T2KOSCANA"]+"inputs/official/2014/flux/sk_tuned13av1.0_13anom_run1-5b.root"
input_files_13atun["RHC"] = os.environ["T2KOSCANA"]+"inputs/official/2014/flux/sk_tuned13av1.0_13anom_run1-5b.root"

# Sukuki-san's files - get 11bv32 from 11a.
input_files_11btun = dict()
input_files_11btun["FHC"] = os.environ["T2KOSCANA"]+"inputs/official/2013/flux/11b_v3.2_run1-4a/sk_tuned11bv3.2_11anom_run1-run4.root"
input_files_11btun["RHC"] = os.environ["T2KOSCANA"]+"inputs/official/2013/flux/11b_v3.2_antinu/sk_tuned11bv3.2_11anom_m250ka.root"

input_files_11btun

species = ["numu","numub","nue","nueb"]

for beam_mode in ["FHC","RHC"]:
#~ for beam_mode in ["nu"]:
  infile_n    = R.TFile(input_files_13anom[beam_mode])
  infile_t    = R.TFile(input_files_13atun[beam_mode])
  infile_11bt = R.TFile(input_files_11btun[beam_mode])
  
  assert(infile_n)
  assert(infile_t)
  assert(infile_11bt)
  
  outfile = R.TFile(os.environ["T2KOSCANA"]+"inputs/flux_reweight/reweight_13aMC_to_13aFlux_"+beam_mode+".root","RECREATE")
  for s in species:
    
    in_hname = "enu_sk_tuned11b_"+s+"_ratio"
    in_hist = infile_n.Get(in_hname)
    assert(in_hist)
        
    tun_hname = "enu_sk_tuned11b_"+s+"_ratio"
    tun_hist = infile_11bt.Get(tun_hname)
    
    #~ weight_hname = "enu_sk_tuned13a_"+s+"_ratio"
    weight_hname = "enu_sk_tuned13a_"+s+"_ratio"
    weight_hist = infile_t.Get(weight_hname)
    assert(weight_hist)
    
    out_hist = R.TH1D(in_hist)
    out_hist.Reset()
    out_hist.GetSumw2().Set(0)
    out_hist.SetTitle("13b tuned flux / 11a nominal flux")
    
    #~ weight_hist.Draw()
    #~ in_hist.SetLineColor(R.kRed)
    #~ in_hist.Draw("same")
    #~ raw_input()
    #~ print in_hist.GetXaxis().GetNbins(),weight_hist.GetXaxis().GetNbins()
    #~ assert(in_hist.GetXaxis().GetNbins()==weight_hist.GetXaxis().GetNbins())
    
    for ibin in xrange(1,out_hist.GetXaxis().GetNbins()+1):
      center = in_hist.GetBinCenter(ibin)
      if center > 10:
        out_hist.Fill(center,1)
        break
      #~ print center,weight_hist.GetBinCenter(ibin)
      assert(abs(center-weight_hist.GetBinCenter(ibin))<1E-10)
      
      wght_11a_nom_to_11b_tun = tun_hist   .GetBinContent(ibin)
      wght_13a_nom_to_11b_tun = in_hist    .GetBinContent(ibin)
      wght_13a_nom_to_13a_tun = weight_hist.GetBinContent(ibin)
      
      wght_11b_tun_to_13a_nom = 1./wght_13a_nom_to_11b_tun
      
      wght_11a_nom_to_13a_tun =  wght_11a_nom_to_11b_tun *  wght_11b_tun_to_13a_nom * wght_13a_nom_to_13a_tun 
      
      out_hist.Fill(center, wght_11a_nom_to_13a_tun)
    
    outfile.cd()
    out_hist.Write(in_hname)

print "Done"
