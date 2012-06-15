//
// ------------------------------------------------------------------------------------------------
// Plot number of events predicted at SuperK (NSK) vs oscillation params mixing13, deltacp
// fixing the other 4 oscillation parameters.
// Also, plot fractional change wrt the 2-flavour hypothesis.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

//#include "Conventions.h"

void test_nsk_twk_osc_mixing13_deltacp_2
(
  double       input_mixing12     = -99,                            /* */
  double       input_mixing23     = -99,                            /* */
  double       input_dm12sq       = -99,                            /* */
  double       input_dm23sq       = -99,                            /* */
  const char * config_file        = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;

  const double pi = TMath::Pi();

  double mixing13_min  =  0.0;
  double mixing13_max  =  0.25;
  double mixing13_step =  0.01;
  double deltacp_min   = -pi;
  double deltacp_max   =  pi;
  double deltacp_step  =  2*pi/90.;

  // make sure that min and max are at the bin centre
  int mixing13_nbins = 1 + (int) ((mixing13_max - mixing13_min) / mixing13_step);
  int deltacp_nbins  = 1 + (int) ((deltacp_max  - deltacp_min ) / deltacp_step );

  TH2D * nsk_vs_osc_params = 
   new TH2D("nsk_vs_osc_params", "number of 1-ring mu-like events",
       mixing13_nbins, 
       mixing13_min - 0.5*mixing13_step, 
       mixing13_max + 0.5*mixing13_step,
       deltacp_nbins, 
       deltacp_min - 0.5*deltacp_step, 
       deltacp_max + 0.5*deltacp_step);

  TH2D * frac_dnsk_vs_osc_params = 
   new TH2D("frac_dnsk_vs_osc_params", "% change in number of 1-ring mu-like events wrt 2-flavour case",
       mixing13_nbins, 
       mixing13_min - 0.5*mixing13_step, 
       mixing13_max + 0.5*mixing13_step,
       deltacp_nbins, 
       deltacp_min - 0.5*deltacp_step, 
       deltacp_max + 0.5*deltacp_step);

  nsk_vs_osc_params       -> GetXaxis() -> SetTitle("sin^{2}2#theta_{13}");
  nsk_vs_osc_params       -> GetYaxis() -> SetTitle("#delta_{CP}");
  frac_dnsk_vs_osc_params -> GetXaxis() -> SetTitle("sin^{2}2#theta_{13}");
  frac_dnsk_vs_osc_params -> GetYaxis() -> SetTitle("#delta_{CP}");

  cout << "Generating predictions for mixing23 = " << input_mixing23 << ", dm23sq = " << input_dm23sq
       << " eV^2, mixing12 = " << input_mixing12 << ", dm12sq = " << input_dm12sq << " eV^2" << endl;

  // create & config superk spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;
  SuperkExpectation skexp;
  skexp.Read(config);

  // Set physics hypothesis
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
  }
  if(input_dm23sq < -1)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1)  input_dm12sq = conv::kGlobFitDm12sq;

  // reference (2-flavour) physics hypothesis
  PhysHypothesis phref(t2kana::conv::kHypNuMuDisapp2F);
  phref.SetCurrValue(conv::kPhypMixing23, input_mixing23);
  phref.SetCurrValue(conv::kPhypDm23sq,   input_dm23sq  );

  // prediction for reference hypothesis
  skexp.Generate(phref);
  double nskref = skexp.CurrPredictedErecoSpectrum(s)->Integral();
  cout << "Nsk[1Rmu,2-flavour] = " << nskref << endl;

  // 3-flavour physics hypothesis, setting fixed params
  PhysHypothesis ph(conv::kHypOsc3F);
  ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
  ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq  );
  ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
  ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq  );
   
  // loop over varied params and generate predictions
  for(int i = 1; i <= nsk_vs_osc_params->GetNbinsX(); i++) {
    for(int j = 1; j <= nsk_vs_osc_params->GetNbinsY(); j++) {
       double mixing13 = nsk_vs_osc_params->GetXaxis()->GetBinCenter(i);
       double deltacp       = nsk_vs_osc_params->GetYaxis()->GetBinCenter(j);
       ph.SetCurrValue(conv::kPhypMixing13,  mixing13);
       ph.SetCurrValue(conv::kPhypDeltaCP,        deltacp      );
       // generate predictions
       skexp.Generate(ph);
       // integrate Ereco distribution 
       double nsk = skexp.CurrPredictedErecoSpectrum(s)->Integral();
       double frac_dnsk = 100*(nsk-nskref)/nskref;
       // store values
       nsk_vs_osc_params->SetBinContent(i,j,nsk);
       frac_dnsk_vs_osc_params->SetBinContent(i,j,frac_dnsk);
       cout << "Nsk[1Rmu,3-flavour](mixing13 = " << mixing13 << ", deltacp = " << deltacp << ") = " << nsk 
            << ", change from the 2-flavour hypothesis = " << frac_dnsk << "%" << endl;
    }//j
  }//i

  // draw
  gStyle->SetOptStat(0);
  TCanvas * c1 = new TCanvas;
  nsk_vs_osc_params       -> Draw("colz");
  TCanvas * c2 = new TCanvas;
  frac_dnsk_vs_osc_params -> Draw("colz");

  // save
  string oscparams = util::AsString("3f_osc",   
     input_mixing12,-1,input_mixing23,-1,input_dm12sq,input_dm23sq,"no-theta13-deltacp");
  TFile f(Form("nsk_vs_mixing13_deltacp-%s.root",oscparams.c_str()),"recreate");
  nsk_vs_osc_params       -> Write("nsk");
  frac_dnsk_vs_osc_params -> Write("frac_dnsk");
  f.Close();
}
