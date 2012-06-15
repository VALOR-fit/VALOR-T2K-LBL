//
// ------------------------------------------------------------------------------------------------
// Calculate SuperK Ereco spectrum error envelope. Tweaking all systematics.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

//#include "Conventions.h"

void generate_sk_ereco_spectrum_err_envelope
(
  const char * hypothesis         = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  const char * mixing             = "sinsq2theta",                  /* options: "sinsq2theta", "sinsqtheta" */
  const char * hierarchy          = "normal",                       /* options: "normal", "inverted" */
  double       input_mixing12     = -99,                            /* */
  double       input_mixing13     = -99,                            /* */
  double       input_mixing23     = -99,                            /* */
  double       input_deltacp      = -99,                            /* */
  double       input_dm12sq       = -99,                            /* */
  double       input_dm23sq       = -99,                            /* */
  const char * config_file        = "analysis_configuration.xml"    /* */
)
{
  gROOT->Reset();

  using namespace std;
  using namespace t2kana;

  const int nexp = 5; // number of toy MC experiments

  // Analysis confifuration
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  // Global options
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixing);
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  cout << *Opt::Instance() << endl;

  // Output histograms and graphs
  TH1D * hEreco [nexp];        //
  TH1D * hErecoFracErrP  = 0;  //
  TH1D * hErecoFracErrM = 0;   //
  TH1D * hErecoDef = 0;        //
  TH1D * hErecoMax = 0;        //
  TH1D * hErecoMin = 0;        //

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
  
  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // Generate nexp random spectra
  for(int iexp = 0; iexp < nexp; iexp++) {
     // tweak spectrum generation model
     if(iexp%50==0) {
        cout << "* Current toy MC experiment: " << iexp << endl;
     }
     if(iexp>0) {
       // use 0^th entry to store the untweaked SK reco energy spectrum
       skexp.TweakAllSystematics();
     }
     // generate predictions
     skexp.Generate(ph);
     // store spectrum
     hEreco[iexp] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike));
     hEreco[iexp]->SetDirectory(0);
  }

  // Store nominal
  hErecoDef = new TH1D(*hEreco[0]);
  hErecoDef->SetDirectory(0);

  // Find min/max and calculate fractional deviation from nominal
  hErecoFracErrP = new  TH1D(*hEreco[0]);
  hErecoFracErrM = new  TH1D(*hEreco[0]);
  hErecoMax      = new TH1D(*hEreco[0]);
  hErecoMin      = new TH1D(*hEreco[0]);
  hErecoFracErrP -> SetDirectory(0);
  hErecoFracErrM -> SetDirectory(0);
  hErecoMax      -> SetDirectory(0);
  hErecoMin      -> SetDirectory(0);
  
  for(int ibin=1; ibin <= conv::kNErecoBins1Rmu; ibin++) {
    double nom   = hErecoDef->GetBinContent(ibin);
    double sumsq = 0;
    for(int iexp = 1; iexp < nexp; iexp++) {
       double N = hEreco[iexp]->GetBinContent(ibin);
       sumsq += TMath::Power(N-nom,2);	
    }//iexp
    double sig  = TMath::Sqrt(sumsq/(nexp-1));
    double min  = nom-sig;
    double max  = nom+sig; 
    double errp = 0;
    double errm = 0;
    if(nom>0){
      errp = (max-nom)/nom;
      errm = (min-nom)/nom;
    }
    hErecoFracErrP -> SetBinContent (ibin,  errp);
    hErecoFracErrM -> SetBinContent (ibin,  errm);
    hErecoMax      -> SetBinContent (ibin,  max);
    hErecoMin      -> SetBinContent (ibin,  min);
  }//ibin
  hErecoFracErrP->GetYaxis()->SetRangeUser(-1,1);
  hErecoFracErrM->GetYaxis()->SetRangeUser(-1,1);

  // Draw
  hErecoDef->SetLineStyle(kDashed);
  hErecoMin->SetLineStyle(kSolid);
  hErecoMax->SetLineStyle(kSolid);
  hErecoMin->SetLineWidth(2);
  hErecoMax->SetLineWidth(2);
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  gStyle->SetOptStat(0);
  c->Divide(1,2);
  c->cd(1);
  hErecoMax  -> Draw();
  hErecoMin  -> Draw("same");
  hErecoDef  -> Draw("same");
  hErecoMax->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoMax->GetYaxis()->SetTitle("Number of events");
  c->cd(2);
  hErecoFracErrP->SetLineWidth(2);
  hErecoFracErrM->SetLineWidth(2);
  hErecoFracErrP -> Draw();
  hErecoFracErrM -> Draw("same");
  hErecoFracErrP->GetXaxis()->SetTitle("E_{reco} (GeV)");
  hErecoFracErrP->GetYaxis()->SetTitle("Fractional error");
  c->Update();

  // Write-out in file
  string oscparams = util::AsString(ph.Hypothesis(),
     input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);
  TFile f(Form("sk_ereco_spectrum-err_envelope-%s.root",oscparams.c_str()),"recreate");
  hErecoFracErrP -> Write("spectrum_frac_err_plus");
  hErecoFracErrM -> Write("spectrum_frac_err_minus");
  hErecoMax      -> Write("spectrum_error_envelope_upper_boundary");
  hErecoMin      -> Write("spectrum_error_envelope_lower_boundary");
  hErecoDef      -> Write("spectrum_nominal");
  f.Close();
}
