//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for various tweaks of
// nd280 efficiency curves
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_sk_ereco_spectrum_twk_nd280_eff
(
  const char * hypothesis         = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12     = 0.86,                           /* */
  double       input_mixing13     = 0.10,                           /* */
  double       input_mixing23     = 1.00,                           /* */
  double       input_deltacp      = 0.00,                           /* */
  double       input_dm12sq       = 7.58E-5,                        /* */
  double       input_dm23sq       = 2.32E-3,                        /* */
  const char * config_file        = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  TRandom3 rand(0);

  // n tweaks
  const int n = 3;

  // book histograms
  TH1D * hEreco[n];

  // set physics hypothesis
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  }
  cout << "Physics hypothesis set: " << endl;
  cout << ph;
  
  // create & config superk spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;
  SuperkExpectation skexp;
  skexp.Read(config);

  // generate predictions
  for(int i = 0; i < n; i++) {
     cout << "Generating SuperK spectrum number: " << i << endl; 
     // tweak spectrum generation model
     if(i>0) {
          double twkdial = (i==1) ? +1 : -1;
          skexp.GetND280Measurement()->TweakEfficiency(twkdial);
     }
     // calculate spectra
     skexp.Generate(ph);  
     // store spectra
     hEreco[i] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
     hEreco[i]->SetDirectory(0);
  }//i

  // calculate err envelopes
  TH1D * alo = new TH1D("alo",  "spectrum err envelope: low  boundary",            conv::kNErecoBins1Rmu,conv::kErecoBinEdges1Rmu);
  TH1D * ahi = new TH1D("ahi",  "spectrum err envelope: high boundary",            conv::kNErecoBins1Rmu,conv::kErecoBinEdges1Rmu);
  TH1D * flo = new TH1D("flo",  "spectrum fractional err envelope: low  boundary", conv::kNErecoBins1Rmu,conv::kErecoBinEdges1Rmu);
  TH1D * fhi = new TH1D("fhi",  "spectrum fractional err envelope: high boundary", conv::kNErecoBins1Rmu,conv::kErecoBinEdges1Rmu);
  for(int ibin=1; ibin <= conv::kNErecoBins1Rmu; ibin++) {
    double nlo =  999999;
    double nhi = -999999;
    double def = hEreco[0]->GetBinContent(ibin);
    for(int i = 0; i < n; i++) {    
       double ncurr = hEreco[i]->GetBinContent(ibin);
       nlo = TMath::Min(nlo, ncurr);
       nhi = TMath::Max(nhi, ncurr);
    }
    alo->SetBinContent(ibin,nlo);
    ahi->SetBinContent(ibin,nhi);
    flo->SetBinContent(ibin, ((def>0) ? (nlo-def)/def : 0) );
    fhi->SetBinContent(ibin, ((def>0) ? (nhi-def)/def : 0) );
  }

  // draw
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  gStyle->SetOptStat(0);
  c->Divide(1,2);
  c->cd(1);
  hEreco[0]->Draw();
  hEreco[0]->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hEreco[0]->GetYaxis()->SetTitle("Number of events");
  alo->SetLineStyle(kDashed);
  ahi->SetLineStyle(kDashed);
  alo->Draw("same");
  ahi->Draw("same");
  c->cd(2);
  TH1F * hframe = (TH1F*) c->GetPad(2)->DrawFrame(0,-0.15,conv::kErecoBinEdges1Rmu[conv::kNErecoBins1Rmu],0.15);
  hframe->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hframe->Draw();
  flo->Draw("same");
  fhi->Draw("same");

  // write-out in file
  string oscparams = util::AsString(hypothesis,
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-twk_nd280_eff-%s.root",oscparams.c_str()),"recreate");
  hEreco[0]->Write("sk_ereco_spectrum_twk_nd280_eff_nominal");
  hEreco[1]->Write("sk_ereco_spectrum_twk_nd280_eff_plus_one_sigma");
  hEreco[2]->Write("sk_ereco_spectrum_twk_nd280_eff_minus_one_sigma");
  f.Close();

}
