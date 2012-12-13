void chisq_for_n_systematics
(
  double       input_mixing12     = -99,                            /* */
  double       input_mixing13     = -99,                            /* */
  double       input_mixing23     = -99,                            /* */
  double       input_deltacp      = -99,                            /* */
  double       input_dm12sq       = -99,                            /* */
  double       input_dm23sq       = -99,                            /* */
  bool         save_toy_expts     = false,                          /* save the ereco spectra? */
  const char * config_file        = "analysis_configuration.xml",   /* */
  const int    nexp               = 100                           /* number of experiments to make */

) 
{
  using namespace std;
  using namespace t2kana;

  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Set physics hypothesis
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;
  FitParamList * fp = FitParamList::NamedParameterList("nuosc3f_numudisapp_fix4");
  if(fp->PhysicsHypothesis() == conv::kHypNuMuDisapp2F) {
    fp->Update(conv::kPhypMixing23,  input_mixing23);
    fp->Update(conv::kPhypDm23sq,    input_dm23sq);
  }
  else if(fp->PhysicsHypothesis() == conv::kHypOsc3F) {
    fp->Update(conv::kPhypMixing12,  input_mixing12);
    fp->Update(conv::kPhypMixing13,  input_mixing13);
    fp->Update(conv::kPhypMixing23,  input_mixing23);
    fp->Update(conv::kPhypDeltaCP,   input_deltacp);
    fp->Update(conv::kPhypDm12sq,    input_dm12sq);
    fp->Update(conv::kPhypDm23sq,    input_dm23sq);
  }

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Physics parameter list: \n" << *fp << endl;
  cout << *Opt::Instance() << endl;

  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  SuperkObservation * skobs = 0;
  FitInputs * fit_inputs = 0;

  // Output histograms
  TH1D * hNsk    [conv::kNumOfNuisanceParams1Rmu] = {0};
  TH1D * hNskExp [conv::kNumOfNuisanceParams1Rmu] = {0};
  TH1D * hChisq  [conv::kNumOfNuisanceParams1Rmu] = {0};
  TH1D * hChisq73[conv::kNumOfNuisanceParams1Rmu] = {0};
  TH1D * hNskSummary = new TH1D("hNskSummary","hNskSummary",conv::kNumOfNuisanceParams1Rmu+1,0,conv::kNumOfNuisanceParams1Rmu);
  TH1D * hNskExpSummary = new TH1D("hNskExpSummary","hNskExpSummary",conv::kNumOfNuisanceParams1Rmu+1,0,conv::kNumOfNuisanceParams1Rmu);
  TH1D * hChisqSummary = new TH1D("hChisqSummary","hChisqSummary",conv::kNumOfNuisanceParams1Rmu+1,0,conv::kNumOfNuisanceParams1Rmu);
  TH1D * hChisq73Summary = new TH1D("hChisq73Summary","hChisq73Summary",conv::kNumOfNuisanceParams1Rmu+1,0,conv::kNumOfNuisanceParams1Rmu);

  // Open output files
  string oscparams = util::AsString(fp->PhysicsHypothesis(),
				    input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("chisq_for_n_systematics_%s.root",oscparams.c_str()),"RECREATE");

  //remove the [] so can use it in SaveAs
  while(oscparams.find_first_of("[]") != string::npos) {
    oscparams.erase(oscparams.find_first_of("[]"), 1);
  }

  TRandom3 rand(0);

  //loop over modes
  for(int ip = 0; ip <= conv::kNumOfNuisanceParams1Rmu; ip++) {

    string modestr = "";
    switch(ip) {
    case 0: modestr = "stat"; break;
    case 1: modestr = ""; break;

    default: modestr = ""; break;
    }

    hNsk[ip]     = new TH1D(Form("Nsk_%s",modestr.c_str()),Form("Nsk_%s",modestr.c_str()),3000,0,300);
    hNskExp[ip]  = new TH1D(Form("NskExp_%s",modestr.c_str()),Form("NskExp_%s",modestr.c_str()),3000,0,300);
    hChisq[ip]   = new TH1D(Form("Chisq_%s",modestr.c_str()),Form("Chisq_%s",modestr.c_str()),1000,0,50);
    hChisq73[ip] = new TH1D(Form("Chisq73_%s",modestr.c_str()),Form("Chisq73_%s",modestr.c_str()),2500,0,150);

    //Loop over experiments
    for(int ie = 0; ie < nexp; ie++) {
      
      skobs = new SuperkObservation;
      skexp->ResetSystematicVariations();

      //setup correlated tweaks
      TVectorD vb(conv::kNumOfParamsInBANFFfit);
      vb = skexp->GetBANFF()->BestFitValues();
      int ib = -1;
      TVectorD vs(3); vs.Zero();
      TVectorD vf(conv::kNumOfFSIParams); vf.Zero();

      bool applied_negevents_condition = true;
      while(applied_negevents_condition) {

	//do the tweaks
	if(ip >= 1) {
	  skexp->GetSuperkDetector()->TweakEfficiency1RmuCCnonQE    (rand.Gaus(0.,1.));
	}
	if(ip >= 2) {
	  skexp->GetSuperkDetector()->TweakEfficiency1RmuNC         (rand.Gaus(0.,1.));
	}
	if(ip >= 3) {
	  skexp->GetSuperkDetector()->TweakEnergyScale              (rand.Gaus(0.,1.));
	}
	if(ip >= 4) {
	  skexp->GetNuXSec()->TweakSpecFunc      (rand.Uniform(0.,1.));
	}
	if(ip >= 5) {
	  skexp->GetNuXSec()->TweakPiLessDelta   (rand.Gaus(0.,1.));
	}
	if(ip >= 6) {
	  skexp->GetNuXSec()->TweakBindE         (rand.Gaus(0.,1.));
	}
	if(ip >= 7) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFFMAQE);
	  v[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 8) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFFCCQEE0);
	  v[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 9) {
	  skexp->GetSuperkDetector()->TweakEfficiency1RmuCCQEOther  (rand.Gaus(0.,1.));
	}
	if(ip >= 10) {
	  vs[0] = rand.Gaus(0.,1.);
	}
	if(ip >= 11) {
	  vs[1] = rand.Gaus(0.,1.);
	}
	if(ip >= 12) {
	  vs[2] = rand.Gaus(0.,1.);
	}
	if(ip >= 13) {
	  skexp->GetSuperkDetector()->TweakEfficiency1RmuNuECC      (rand.Gaus(0.,1.));
	}
	if(ip >= 14) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF00);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 15) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF01);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 16) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF02);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 17) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF03);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 18) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF04);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 19) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF05);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 20) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF06);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 21) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF07);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 22) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF08);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 23) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF09);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 24) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF10);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 25) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF11);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 26) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF12);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 27) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF13);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 28) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF14);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 29) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFF15);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 30) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFFMARES);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 31) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFFCCQEE1);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 32) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFFCCQEE2);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 33) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFFCC1piE0);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 34) {
	  ib = util::ArrIdx::Instance()->BANFFSystIdx(conv::kNuisBANFFCC1piE1);
	  vb[ib] = skexp->GetBANFF()->BestFitValues()[ib] + (tweak * sqrt(skexp->GetBANFF()->CovMatrix()[ib][ib]));
	}
	if(ip >= 35) {
	  skexp->GetNuXSec()->TweakFermiMom      (rand.Gaus(0.,1.));
	}
	if(ip >= 36) {
	  skexp->GetNuXSec()->TweakWshape        (rand.Gaus(0.,1.));
	}
	if(ip >= 37) {
	  skexp->GetNuXSec()->TweakShapeCCoth    (rand.Gaus(0.,1.));
	}
	if(ip >= 38) {
	  skexp->GetNuXSec()->TweakNormCCcoh     (rand.Gaus(0.,1.));
	}
	if(ip >= 39) {
	  skexp->GetNuXSec()->TweakNormNC1piC    (rand.Gaus(0.,1.));
	}
	if(ip >= 40) {
	  skexp->GetNuXSec()->TweakNormNCoth     (rand.Gaus(0.,1.));
	}
	if(ip >= 41) {
	  skexp->GetNuXSec()->TweakNormNuEToNuMu (rand.Gaus(0.,1.));
	}
	if(ip >= 42) {
	  skexp->GetNuXSec()->TweakNormAntiNuToNu(rand.Gaus(0.,1.));
	}
	if(ip >= 43) {
	  vf[0] = rand.Gaus(0.,1.);
	}
	if(ip >= 44) {
	  vf[1] = rand.Gaus(0.,1.);
	}
	if(ip >= 45) {
	  vf[2] = rand.Gaus(0.,1.);
	}
	if(ip >= 46) {
	  vf[3] = rand.Gaus(0.,1.);
	}
	if(ip >= 47) {
	  vf[4] = rand.Gaus(0.,1.);
	}
	if(ip >= 48) {
	  vf[5] = rand.Gaus(0.,1.);
	}

	//do correllated tweaks
	skexp->GetBANFF()->Tweak(true,&vb,1);
	skexp->GetSuperkDetector()->TweakEfficiency1RmuCCQERingCounting(true, vs[0], vs[1], vs[2]);
	skexp->GetFSI()->Tweak(true,&vf);

	//get the observation
	//	skexp->PrintSystematicVariations();
	applied_negevents_condition = skobs->CreateMockData(*skexp, *fp, false);
      }

      //get the nsk
      double nobs = 0, nex = 0;
      for(int ibin = 1; ibin <= conv::kNErecoBins1Rmu; ibin++) {
	nobs += skobs->ErecoSpectrumArray             (s)[ibin-1];
	nex  += skexp->CurrPredictedErecoSpectrumArray(s)[ibin-1];
      }
      hNsk[ip]->Fill(nobs);
      hNskExp[ip]->Fill(nex);

      //skexp->PrintSystematicVariations();

      //reset expectation to find the chisq
      skexp->ResetSystematicVariations();
      
      //get the chisqs
      fit_inputs = new FitInputs(skobs,skexp,fp);
      double chisq73 = fit_inputs->Chisquare("poisson", s);
      hChisq73[ip]->Fill(chisq73);

      fit_inputs->SetChisqDefinition(true);
      double chisq = fit_inputs->Chisquare("poisson", s);
      hChisq[ip]->Fill(chisq);


      cout << " observed:  " << nobs
	   << " expected:  " << nexp
	   << " chisq(5):  " << chisq
	   << " chisq(73): " << chisq73
	   << endl;

      delete skobs;
      delete fit_inputs;

    }//iexp

    hNsk    [ip]->Write();
    hNskExp [ip]->Write();
    hChisq  [ip]->Write();
    hChisq73[ip]->Write();

    hNskSummary    ->SetBinContent(ip,hNsk[ip]->GetRMS()     / hNsk[ip]->GetMean());
    hNskExpSummary ->SetBinContent(ip,hNskExp[ip]->GetRMS()  / hNskExp[ip]->GetMean());
    hChisqSummary  ->SetBinContent(ip,hChisq[ip]->GetRMS()   / hChisq[ip]->GetMean());
    hChisq73Summary->SetBinContent(ip,hChisq73[ip]->GetRMS() / hChisq73[ip]->GetMean());

    if(ip == 2) break;
  }//ip


  hNskSummary    ->Write();
  hNskExpSummary ->Write();
  hChisqSummary  ->Write();
  hChisq73Summary->Write();

  f.Close();

  delete skexp;

}
