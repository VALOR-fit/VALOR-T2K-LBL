
import ROOT
import os

R = ROOT

input_files = dict()
input_files["nu" ] = os.environ["T2KOSCANA"]+"inputs/official/2013/flux/11b_v3.2_run1-4a/sk_tuned11bv3.2_11anom_run1-run4_fine.root"
input_files["anu"] = os.environ["T2KOSCANA"]+"inputs/official/2013/flux/11b_v3.2_antinu/sk_tuned11bv3.2_11anom_m250ka_fine.root"

angles = [ "2p0","2p25" ]
species = ["numu","numub","nue","nueb"]

weight_file = R.TFile(os.environ["T2KOSCANA"]+"inputs/test/flux_offaxis_weights.root")
assert(weight_file)

for beam_mode in input_files.iterkeys():
  infile = R.TFile(input_files[beam_mode])
  for angle in angles:
    outfile = R.TFile("reweighted_flux_"+angle+"_"+beam_mode+".root","RECREATE")
    for s in species:
      
      in_hname = "enu_sk_tuned11b_"+s+"_ratio"
      in_hist = infile.Get(in_hname)
      assert(in_hist)
      
      weight_hname = "weights_oa_"+angle+"_"+beam_mode+"mode_"+s
      weight_hist = weight_file.Get(weight_hname)
      assert(weight_hist)
      
      out_hist = R.TH1D(in_hist)
      out_hist.Reset()
      
      for ibin in xrange(1,out_hist.GetXaxis().GetNbins()+1):
        center = in_hist.GetBinCenter(ibin)
        assert(abs(center-weight_hist.GetBinCenter(ibin))<1E-10)
        out_hist.Fill(center,in_hist.GetBinContent(ibin)*weight_hist.GetBinContent(ibin))
      
      outfile.cd()
      out_hist.Write(in_hname)

print "Done"
