//
// ------------------------------------------------------------------------------------------------
// Calculate SuperK Ereco spectrum error envelope due to BANFF parameter variations.
// BANFF parameter variations generated using the Cholseky method.
// The error in each Ereco bin is 1 standard deviation from N random fluxes
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_sk_ereco_spectrum_twk_banff
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

  const int nexp = 200; // number of BANFF param tweaks

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
    
  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Set physics hypothesis
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

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Physics hypothesis: \n" << ph << endl;
  cout << *Opt::Instance() << endl;

  // Output histograms and graphs
  TH1D * hEreco [nexp];  //
  TH1D * hErecoMax = 0;  //
  TH1D * hErecoMin = 0;  //
  TH1D * hErecoRatio1EtrueBinVarMax[conv::kNumOfParamsInBANFFfit]; //
  TH1D * hErecoRatio1EtrueBinVarMin[conv::kNumOfParamsInBANFFfit]; //

  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);
  
  // Generate nexp random spectra
  for(int iexp = 0; iexp < nexp; iexp++) {
     // tweak spectrum generation model
     cout << " ** iteration: " << iexp << endl;
     if(iexp>0) {
       // use 0^th entry to store the untweaked SK reco energy spectrum
       skexp.GetBANFF()->Tweak();       
     }
     skexp.PrintSystematicVariations();
     // generate predictions
     skexp.Generate(ph);
     // store spectrum
     hEreco[iexp] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
     hEreco[iexp]->SetDirectory(0);
  }

  // Find min/max and calculate fractional deviation from nominal
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
  hErecoMax->GetYaxis()->SetRangeUser(-2,2);
  hErecoMin->GetYaxis()->SetRangeUser(-2,2);

  //
  // Check the effect in Ereco distributions by modifying each of the variable Etrue bins
  // (in the Beam WG variable binning scheme used in the flux covariance matrix) by +/-100%
  //

  bool force_input_twk_dials = true;
  bool correlated            = true;

  for(int iparam = 0; iparam < conv::kNumOfParamsInBANFFfit; iparam++)
  {
       //
       // init
       //
       cout << "- Reseting all BANFF tweaking fials " << endl;
       TVectorD twkv(conv::kNumOfParamsInBANFFfit);
       for(int k = 0; k < conv::kNumOfParamsInBANFFfit; k++) {
         twkv[k] = 0.;
       }
       skexp.GetBANFF()->Tweak(force_input_twk_dials,0); //reset
       //skexp.GetBANFF()->Tweak(force_input_twk_dials,&twkv,correlated);       

       // tweak iparam only neglecting correlations
       conv::PhysParamId_t p = (conv::PhysParamId_t) ( (int)conv::kNuisBANFF00 + iparam );
       cout << "- Tweaking BANFF parameter: " << iparam << " [" << util::AsString(p) << "] only (neglecting correlations)" << endl;

       //
       // tweak current param only by +100% - neglect correlations
       //
       cout << "- Tweaking param by +100%" << endl;
       twkv[iparam] = +1; // +100%
       skexp.GetBANFF()->Tweak(force_input_twk_dials,&twkv,correlated);  
       skexp.PrintSystematicVariations();
       // generate predictions
       skexp.Generate(ph);
       // store spectrum
       hErecoRatio1EtrueBinVarMax[iparam] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
       hErecoRatio1EtrueBinVarMax[iparam]->SetDirectory(0);
       hErecoRatio1EtrueBinVarMax[iparam]->SetTitle(Form("BANFF_%d, +100 percent change",iparam));
       hErecoRatio1EtrueBinVarMax[iparam]->Divide(hEreco[0]);

       //
       // tweak current param only by -100% - neglect correlations
       //
       cout << "- Tweaking by -100%" << endl;
       twkv[iparam] = -0.5; // -100%
       skexp.GetBANFF()->Tweak(force_input_twk_dials,&twkv,correlated);       
       skexp.PrintSystematicVariations();
       // generate predictions
       skexp.Generate(ph);
       // store spectrum
       hErecoRatio1EtrueBinVarMin[iparam] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
       hErecoRatio1EtrueBinVarMin[iparam]->SetDirectory(0);
       hErecoRatio1EtrueBinVarMin[iparam]->SetTitle(Form("BANFF_%d, -100 percent change",iparam));
       hErecoRatio1EtrueBinVarMin[iparam]->Divide(hEreco[0]);
  }

  //
  // Draw
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
  c2->Divide(4,7);
  c2->Update();
  for(int iparam = 0; iparam < conv::kNumOfParamsInBANFFfit; iparam++)
  {
     c2->cd(iparam+1);
     hErecoRatio1EtrueBinVarMax[iparam]->SetLineColor(kBlue);
     hErecoRatio1EtrueBinVarMin[iparam]->SetLineColor(kRed);
     hErecoRatio1EtrueBinVarMax[iparam]->Draw();
     hErecoRatio1EtrueBinVarMin[iparam]->Draw("same");
  }
  
  // Write-out in file
  string oscparams = util::AsString(ph.Hypothesis(),
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-twk_banff-%s.root",oscparams.c_str()),"recreate");
  hErecoMax->Write("error_envelope_upper_boundary");
  hErecoMin->Write("error_envelope_lower_boundary");
  for(int iexp = 0; iexp < nexp; iexp++) {
     const char * title = Form("sk_ereco_spectrum_twk_beam_shape_try%d",iexp);
     hEreco[iexp]->Write(title);
  }
  for(int iparam = 0; iparam < conv::kNumOfParamsInBANFFfit; iparam++)
  {
     hErecoRatio1EtrueBinVarMax[iparam]->Write(Form("ereco_spectrum_ratio_twk_banff_%d_by_plus_100pct", iparam));
     hErecoRatio1EtrueBinVarMin[iparam]->Write(Form("ereco_spectrum_ratio_twk_banff_%d_by_minus_100pct",iparam));
  }

  f.Close();
}
