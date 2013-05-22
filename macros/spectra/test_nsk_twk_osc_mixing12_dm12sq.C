//
// ------------------------------------------------------------------------------------------------
// Plot number of events predicted at SuperK (NSK) vs oscillation params sinq2theta12, dm12sq
// fixing the other 4 oscillation parameters.
// Also, plot fractional change wrt the 2-flavour hypothesis.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

//#include "Conventions.h"

void test_nsk_twk_osc_mixing12_dm12sq
(
  double       input_mixing13     = -99,                            /* */
  double       input_mixing23     = -99,                            /* */
  double       input_deltacp      = -99,                            /* */
  double       input_dm23sq       = -99,                            /* */
  const char * config_file        = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;

  double mixing12_min  = 0, mixing12_max = 0;
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    mixing12_min  = conv::kGlobFitSinsq2Theta12 * (1 - 3*conv::kGlobFitErrSinsq2Theta12);
    mixing12_max  = conv::kGlobFitSinsq2Theta12 * (1 + 3*conv::kGlobFitErrSinsq2Theta12);
  }
  else{
    mixing12_min  = conv::kGlobFitSinsqTheta12 * (1 - 3*conv::kGlobFitErrSinsqTheta12);
    mixing12_max  = conv::kGlobFitSinsqTheta12 * (1 + 3*conv::kGlobFitErrSinsqTheta12);
  }
  double mixing12_step = (mixing12_max-mixing12_min)/50.;
  double dm12sq_min    = conv::kGlobFitDm12sq * (1 - 3*conv::kGlobFitErrDm12sq);
  double dm12sq_max    = conv::kGlobFitDm12sq * (1 + 3*conv::kGlobFitErrDm12sq);
  double dm12sq_step   = (dm12sq_max-dm12sq_min)/50.;

  // make sure that min and max are at the bin centre
  int mixing12_nbins = 1 + TMath::CeilNint ((mixing12_max - mixing12_min) / mixing12_step);
  int dm12sq_nbins   = 1 + TMath::CeilNint ((dm12sq_max   - dm12sq_min  ) / dm12sq_step  );

  TH2D * nsk_vs_osc_params = 
   new TH2D("nsk_vs_osc_params", "number of 1-ring mu-like events",
       mixing12_nbins, 
       mixing12_min - 0.5*mixing12_step, 
       mixing12_max + 0.5*mixing12_step,
       dm12sq_nbins, 
       dm12sq_min - 0.5*dm12sq_step, 
       dm12sq_max + 0.5*dm12sq_step);

  TH2D * frac_dnsk_vs_osc_params = 
   new TH2D("frac_dnsk_vs_osc_params", "% change in number of 1-ring mu-like events wrt 2-flavour case",
       mixing12_nbins, 
       mixing12_min - 0.5*mixing12_step, 
       mixing12_max + 0.5*mixing12_step,
       dm12sq_nbins, 
       dm12sq_min - 0.5*dm12sq_step, 
       dm12sq_max + 0.5*dm12sq_step);

  nsk_vs_osc_params       -> GetXaxis() -> SetTitle("sin^{2}2#theta_{12}");
  nsk_vs_osc_params       -> GetYaxis() -> SetTitle("#Deltam_{12}^{2}");
  frac_dnsk_vs_osc_params -> GetXaxis() -> SetTitle("sin^{2}2#theta_{12}");
  frac_dnsk_vs_osc_params -> GetYaxis() -> SetTitle("#Deltam_{12}^{2}");

  cout << "Generating predictions for mixing23 = " << input_mixing23 << ", dm23sq = " << input_dm23sq
       << " eV^2, mixing13 = " << input_mixing13 << ", deltacp = " << input_deltacp << endl;

  // create & config superk spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;
  SuperkExpectation skexp;
  skexp.Read(config);

  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsq2Theta23;
  }
  else {
    if(input_mixing13 < -1) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1) input_mixing23 = conv::kGlobFitSinsqTheta23;
  }
  if(input_dm23sq < -1)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_deltacp < -90) input_deltacp = conv::kGlobFitDeltaCP;

  // reference (2-flavour) physics hypothesis
  PhysHypothesis phref(t2kana::conv::kHypNuMuDisapp2F);
  phref.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
  phref.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq       );

  // prediction for reference hypothesis
  skexp.Generate(phref);
  double nskref = skexp.CurrPredictedErecoSpectrum(s)->Integral();
  cout << "Nsk[1Rmu,2-flavour] = " << nskref << endl;

  // 3-flavour physics hypothesis, setting fixed params
  PhysHypothesis ph(conv::kHypOsc3F);
  ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
  ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq       );
  ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
  ph.SetCurrValue(conv::kPhypDeltaCP,        input_deltacp      );
   
  // loop over varied params and generate predictions
  for(int i = 1; i <= nsk_vs_osc_params->GetNbinsX(); i++) {
    for(int j = 1; j <= nsk_vs_osc_params->GetNbinsY(); j++) {
       double mixing12 = nsk_vs_osc_params->GetXaxis()->GetBinCenter(i);
       double dm12sq        = nsk_vs_osc_params->GetYaxis()->GetBinCenter(j);
       ph.SetCurrValue(conv::kPhypMixing12,  mixing12);
       ph.SetCurrValue(conv::kPhypDm12sq,         dm12sq       );
       // generate predictions
       skexp.Generate(ph);
       // integrate Ereco distribution 
       double nsk = skexp.CurrPredictedErecoSpectrum(s)->Integral();
       double frac_dnsk = 100*(nsk-nskref)/nskref;
       // store values
       nsk_vs_osc_params->SetBinContent(i,j,nsk);
       frac_dnsk_vs_osc_params->SetBinContent(i,j,frac_dnsk);
       cout << "Nsk[1Rmu,3-flavour](mixing12 = " << mixing12 << ", dm12sq = " << dm12sq << " eV^2) = " << nsk 
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
  string oscparams = util::AsString(conv::kHypOsc3F,   
     -1,input_mixing13,input_mixing23,input_deltacp,-1,input_dm23sq);
  TFile f(Form("nsk_vs_mixing12_dm12sq-%s.root",oscparams.c_str()),"recreate");
  nsk_vs_osc_params       -> Write("nsk");
  frac_dnsk_vs_osc_params -> Write("frac_dnsk");
  f.Close();
}
