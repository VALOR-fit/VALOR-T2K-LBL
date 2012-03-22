//
// ------------------------------------------------------------------------------------------------
// Calculate value of chisq penalty term due to BANFF parameter variations.
// BANFF parameter variations generated using the Cholseky method.
// A second method is intended to be implemented as cross check.
// based on the macro test_sk_ereco_spectrum_twk_banff
// ------------------------------------------------------------------------------------------------
//

void test_banff_chisq_penalty_term
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

  const int nexp = 5000; // number of BANFF param tweaks

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
    
  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;
  Opt::Instance()->gOscAnaTypeChoice     = conv::kAnaCombined; //options: kAnaCombined, kAnaNuEAppOnly, kAnaNuMuDisappOnly

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

  // Output histograms 
  TH2D * BANFFcov  = new TH2D("BANFFcov", "BANFF Covariance Matrix", conv::kNumOfParamsInBANFFfit, 0, conv::kNumOfParamsInBANFFfit-1, conv::kNumOfParamsInBANFFfit, 0, conv::kNumOfParamsInBANFFfit-1  );  //
  TH1D * chisqpenalty  = new TH1D("chisqpenalty", "Chisq Penalty Term", 100, 0, 50);  //
  
  TH1D * hdeviations    [conv::kNumOfParamsInBANFFfit];  //
  const int np = conv::kNumOfParamsInBANFFfit;
  for(int i=0; i<np; i++) {
    hdeviations[i] = new TH1D(Form("hdeviation%d", i), Form("hdeviation%d", i), 100, -10, 10);
  }

  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  //Checking loaded inputs in BANFF class
  TMatrixD cov = skexp.GetBANFF()->CovMatrix(); 
  cov.Print();
  for(int i=0; i<np; i++) {
    for(int j=0; j<np; j++) {
      BANFFcov->SetBinContent(i+1,j+1,cov[i][j]);
    }
  }

  cout << "Best fit values: " << endl;
  TVectorD bestfit(conv::kNumOfParamsInBANFFfit);
  bestfit = skexp.GetBANFF()->BestFitValues(); 
  bestfit.Print();

  //
  // Generate nexp random spectra
  //
  for(int iexp = 0; iexp < nexp; iexp++) {
    // tweak spectrum generation model
    cout << " ** iteration: " << iexp << endl;
    // use 0^th entry to store the untweaked SK reco energy spectrum
    skexp.GetBANFF()->Tweak();       
    double chisq = skexp.GetBANFF()->ChisqPenalty();       
    chisqpenalty -> Fill(chisq);
  }
  
  //
  // Check the effect in Ereco distributions by modifying each of the variable Etrue bins
  // (in the Beam WG variable binning scheme used in the flux covariance matrix) by +/-100%
  //

  bool force = true;
  TVectorD chisq_sigma       (conv::kNumOfParamsInBANFFfit);//vector for chisq after +100% forced deviations
  TVectorD chisq_minus_sigma (conv::kNumOfParamsInBANFFfit);//vector for chisq after -100% forced deviations
  TVectorD chisq_sigma2      (conv::kNumOfParamsInBANFFfit);//vector for chisq after +100% forced deviations - second method
  TVectorD chisq_minus_sigma2(conv::kNumOfParamsInBANFFfit);//vector for chisq after -100% forced deviations - second method

  for(int iparam = 0; iparam < conv::kNumOfParamsInBANFFfit; iparam++)
  {
       //
       // init
       //
       cout << "- Reseting all BANFF parameters to best-fit ones " << endl;
       skexp.GetBANFF()->Tweak(force,0); //reset

       TVectorD v(conv::kNumOfParamsInBANFFfit);//vector for forced deviations
       v = skexp.GetBANFF()->BestFitValues();

       conv::PhysParamId_t p = (conv::PhysParamId_t) ( (int)conv::kNuisBANFF00 + iparam );
       cout << "- Tweaking BANFF parameter: " << iparam << " [" << util::AsString(p) << "] only (neglecting correlations)" << endl;

       //
       // tweak current param only by +100% - neglect correlations
       //
       cout << "- Tweaking param by +100%" << endl;
       v[iparam] = 2 * skexp.GetBANFF()->BestFitValue(p); // +100%
       skexp.GetBANFF()->Tweak(force,&v);  
       double chisq = skexp.GetBANFF()->ChisqPenalty();
       chisq_sigma [iparam] = chisq;

       //
       // tweak current param only by -100% - neglect correlations
       //
       cout << "- Tweaking by -100%" << endl;
       v[iparam] = 0.5 * skexp.GetBANFF()->BestFitValue(p); // -100%
       skexp.GetBANFF()->Tweak(force,&v);  
       double chisq_minus = skexp.GetBANFF()->ChisqPenalty();
       chisq_minus_sigma [iparam] = chisq_minus;
  }
  
  //
  // Print out
  //

  cout << " iparam | \t +100% | \t -100% " << endl;
  for(int iparam = 0; iparam < conv::kNumOfParamsInBANFFfit; iparam++)
    {
      cout << iparam << "   |   " << chisq_sigma [iparam] << "   |   " << chisq_minus_sigma [iparam]  << endl;
    }

  //
  // Draw
  //
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  gStyle->SetOptStat(0);
  //  gStyle->SetPalette(1);
  BANFFcov->Draw("colz");

  TCanvas * c1 = new TCanvas("c1","",20,20,500,500);
  gStyle->SetOptStat(0);
  gStyle->SetFillColor(kWhite);
  chisqpenalty->Draw();

  // Write-out in file
  TFile f("banff_chisq_penalty_term.root","recreate");
  BANFFcov->Write();
  chisqpenalty->Write();
  f.Close();

}
