//
// ------------------------------------------------------------------------------------------------
// Calculate SuperK Ereco spectrum error envelope due to flux shape systematic.
// Flux systematic variations generated using the Cholseky method.
// The error in each Ereco bin is 1 standard deviation from N random fluxes
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_sk_ereco_spectrum_twk_beam_shape
(
  const char * hypothesis              = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_sinsq2theta12     = 0.86,                           /* */
  double       input_sinsq2theta13     = 0.10,                           /* */
  double       input_sinsq2theta23     = 1.00,                           /* */
  double       input_deltacp           = 0.00,                           /* */
  double       input_dm12sq            = 7.58E-5,                        /* */
  double       input_dm23sq            = 2.32E-3,                        /* */
  const char * config_file             = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  const int nexp       = 1000; // number of 'flux tweaks'
  const int netruebins = 19;   // number of Etrue (variable size) bins in flux covariance matrix
  // beam WG binning:
  double bin_edges[netruebins+1] = {
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 7.0, 10.0 };

  // output histograms and graphs
  TH1D * hEreco [nexp];  //
  TH1D * hErecoMax = 0;  //
  TH1D * hErecoMin = 0;  //
  TH1D * hErecoRatio1EtrueBinVarMax[netruebins+1]; //
  TH1D * hErecoRatio1EtrueBinVarMin[netruebins+1]; //
  TGraph * grNumEvtVersusToyExpId = 0;
  double nev  [nexp];
  double expid[nexp];

  // set physics hypothesis
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypSinsq2Theta23,  input_sinsq2theta23);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypSinsq2Theta12,  input_sinsq2theta12);
    ph.SetCurrValue(conv::kPhypSinsq2Theta13,  input_sinsq2theta13);
    ph.SetCurrValue(conv::kPhypSinsq2Theta23,  input_sinsq2theta23);
    ph.SetCurrValue(conv::kPhypDeltaCP,        input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,         input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq);
  }
  cout << "Physics hypothesis set: " << endl;
  cout << ph;
  
  // create & config superk spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;
  SuperkExpectation skexp;
  skexp.Read(config);
  
  // generate nexp random spectra
  for(int iexp = 0; iexp < nexp; iexp++) {
     // tweak spectrum generation model
     cout << "** iteration: " << iexp << endl;
     if(iexp>0) {
       // use 0^th entry to store the untweaked SK reco energy spectrum
       skexp.GetNuFlux()->TweakShape();       
     }
     // generate predictions
     cout << "Generating SK predictions " << endl;
     skexp.Generate(ph);
     // store spectrum
     hEreco[iexp] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
     hEreco[iexp]->SetDirectory(0);
  }

  // find min/max and calculate fractional deviation from nominal
  hErecoMax = new TH1D(*hEreco[0]);
  hErecoMin = new TH1D(*hEreco[0]);
  hErecoMax->SetDirectory(0);
  hErecoMin->SetDirectory(0);
  for(int ibin=1; ibin <= conv::kNErecoBins1Rmu; ibin++) {
    double nom   = hEreco[0]->GetBinContent(ibin);
    double sumsq = 0;
    for(int iexp = 1; iexp < nexp; iexp++) {
       double N = hEreco[iexp]->GetBinContent(ibin);
       sumsq += TMath::Power(N-nom,2);	
    }//iexp
    double sig = TMath::Sqrt(sumsq/(nexp-1));
    double min = nom-sig;
    double max = nom+sig; 
    if(nom>0){
      hErecoMax -> SetBinContent (ibin,  (max-nom)/nom);
      hErecoMin -> SetBinContent (ibin,  (min-nom)/nom);
    }
  }//ibin
  hErecoMax->GetYaxis()->SetRangeUser(-1,1);
  hErecoMin->GetYaxis()->SetRangeUser(-1,1);

  // get number of events vs toy experiment id to illustrate that there is no
  // normalization change
  for(int iexp = 1; iexp < nexp; iexp++) {
    nev  [iexp] = hEreco[iexp]->Integral();
    expid[iexp] = iexp;
  }
  grNumEvtVersusToyExpId = new TGraph(nexp,expid,nev);

  //
  // check the effect in Ereco distributions by modifying each of the variable Etrue bins
  // (in the Beam WG variable binning scheme used in the flux covariance matrix) by +/-100%
  //

  bool force_input_twk_dials = true;
  bool correlated            = true;
  bool fractional            = true;

  for(int ietrue_bin = 1; ietrue_bin <= netruebins; ietrue_bin++)
  {
       cout << "Number of events (nominal) = " << hEreco[0]->Integral() << endl;

       // init
       TVectorD twkv(netruebins);
       for(int k = 0; k < netruebins; k++) {
         twkv[k] = 0.;
       }

       // tweak given Etrue bin only
       twkv[ietrue_bin-1] = +1; // +100%
       skexp.GetNuFlux()->TweakShape(force_input_twk_dials,&twkv,correlated,fractional);       

       // generate predictions
       skexp.Generate(ph);

       // store spectrum
       hErecoRatio1EtrueBinVarMax[ietrue_bin] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
       hErecoRatio1EtrueBinVarMax[ietrue_bin]->SetDirectory(0);
       hErecoRatio1EtrueBinVarMax[ietrue_bin]->SetTitle(Form("Flux shape, %.3f - %.3f GeV, +100 percent change",bin_edges[ietrue_bin-1],bin_edges[ietrue_bin]));
       cout << "Number of events (Etrue variable bin " << ietrue_bin << ", +100%) = " << hErecoRatio1EtrueBinVarMax[ietrue_bin]->Integral() << endl;
       hErecoRatio1EtrueBinVarMax[ietrue_bin]->Divide(hEreco[0]);

       // tweak given Etrue bin only
       twkv[ietrue_bin-1] = -1; // -100%
       skexp.GetNuFlux()->TweakShape(force_input_twk_dials,&twkv,correlated,fractional);       

       // generate predictions
       skexp.Generate(ph);

       // store spectrum
       hErecoRatio1EtrueBinVarMin[ietrue_bin] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
       hErecoRatio1EtrueBinVarMin[ietrue_bin]->SetDirectory(0);
       hErecoRatio1EtrueBinVarMin[ietrue_bin]->SetTitle(Form("Flux shape, %.3f - %.3f GeV, -100 percent change",bin_edges[ietrue_bin-1],bin_edges[ietrue_bin]));
       cout << "Number of events (Etrue variable bin " << ietrue_bin << ", -100%) = " << hErecoRatio1EtrueBinVarMin[ietrue_bin]->Integral() << endl;
       hErecoRatio1EtrueBinVarMin[ietrue_bin]->Divide(hEreco[0]);
  }

  //
  // draw
  //
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  gStyle->SetOptStat(0);
  c->Divide(1,2);
  c->cd(1);
  for(int iexp = 0; iexp < nexp; iexp++) {
    if(iexp==0) {
      hEreco[iexp]->GetXaxis()->SetTitle("E_{reco} (GeV)");
      hEreco[iexp]->GetYaxis()->SetTitle("Number of events");
      hEreco[iexp]->SetLineWidth(2);
      hEreco[iexp]->Draw();
    } else {
      hEreco[iexp]->Draw("same");
   }
  }//iexp
  c->cd(2);
  hErecoMax->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoMax->GetYaxis()->SetTitle("Relative error");
  hErecoMax->Draw();
  hErecoMin->Draw("same");
  c->Update();
  TCanvas * c2 = new TCanvas("c2","",20,20,500,500);
  c2->Divide(5,4);
  for(int ietrue_bin = 1; ietrue_bin <= netruebins; ietrue_bin++)
  {
     c2->cd(ietrue_bin);
     hErecoRatio1EtrueBinVarMax[ietrue_bin]->SetLineColor(kBlue);
     hErecoRatio1EtrueBinVarMin[ietrue_bin]->SetLineColor(kRed);
     hErecoRatio1EtrueBinVarMax[ietrue_bin]->Draw();
     hErecoRatio1EtrueBinVarMin[ietrue_bin]->Draw("same");
  }
  
  // write-out in file
  string oscparams = util::AsString(hypothesis,
     input_sinsq2theta12,input_sinsq2theta13,input_sinsq2theta23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-twk_beam_shape-%s.root",oscparams.c_str()),"recreate");
  hErecoMax->Write("error_envelope_upper_boundary");
  hErecoMin->Write("error_envelope_lower_boundary");
  for(int iexp = 0; iexp < nexp; iexp++) {
     const char * title = Form("sk_ereco_spectrum_twk_beam_shape_try%d",iexp);
     hEreco[iexp]->Write(title);
  }
  grNumEvtVersusToyExpId->Write("num_of_events_vs_toy_experiment");
  for(int ietrue_bin = 1; ietrue_bin <= netruebins; ietrue_bin++)
  {
     hErecoRatio1EtrueBinVarMax[ietrue_bin]->Write(Form("ereco_spectrum_ratio_twk_etrue_varbin_%d_by_plus_100pct", ietrue_bin));
     hErecoRatio1EtrueBinVarMin[ietrue_bin]->Write(Form("ereco_spectrum_ratio_twk_etrue_varbin_%d_by_minus_100pct",ietrue_bin));
  }

  f.Close();
}
