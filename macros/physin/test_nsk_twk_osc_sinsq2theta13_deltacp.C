//
// ------------------------------------------------------------------------------------------------
// Plot number of events predicted at SuperK (NSK) vs oscillation params sinsq2theta13, deltacp
// fixing the other 4 oscillation parameters.
// Also, plot fractional change wrt the 2-flavour hypothesis.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_nsk_twk_osc_sinsq2theta13_deltacp
(
  double       input_sinsq2theta23     = 1.00,                           /* */
  double       input_dm23sq            = 2.32E-3,                        /* */
  double       input_sinsq2theta12     = 0.86,                           /* */
  double       input_dm12sq            = 7.58E-5,                        /* */
  const char * config_file             = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  const double pi = TMath::Pi();

  double sinsq2theta13_min  =  0.0;
  double sinsq2theta13_max  =  0.25;
  double sinsq2theta13_step =  0.01;
  double deltacp_min        = -pi;
  double deltacp_max        =  pi;
  double deltacp_step       =  2*pi/90.;

  // make sure that min and max are at the bin centre
  int sinsq2theta13_nbins = 1 + (int) ((sinsq2theta13_max - sinsq2theta13_min) / sinsq2theta13_step);
  int deltacp_nbins       = 1 + (int) ((deltacp_max       - deltacp_min      ) / deltacp_step      );

  TH2D * nsk_vs_osc_params = 
   new TH2D("nsk_vs_osc_params", "number of 1-ring mu-like events",
       sinsq2theta13_nbins, 
       sinsq2theta13_min - 0.5*sinsq2theta13_step, 
       sinsq2theta13_max + 0.5*sinsq2theta13_step,
       deltacp_nbins, 
       deltacp_min - 0.5*deltacp_step, 
       deltacp_max + 0.5*deltacp_step);

  TH2D * frac_dnsk_vs_osc_params = 
   new TH2D("frac_dnsk_vs_osc_params", "% change in number of 1-ring mu-like events wrt 2-flavour case",
       sinsq2theta13_nbins, 
       sinsq2theta13_min - 0.5*sinsq2theta13_step, 
       sinsq2theta13_max + 0.5*sinsq2theta13_step,
       deltacp_nbins, 
       deltacp_min - 0.5*deltacp_step, 
       deltacp_max + 0.5*deltacp_step);

  nsk_vs_osc_params       -> GetXaxis() -> SetTitle("sin^{2}2#theta_{13}");
  nsk_vs_osc_params       -> GetYaxis() -> SetTitle("#delta_{CP}");
  frac_dnsk_vs_osc_params -> GetXaxis() -> SetTitle("sin^{2}2#theta_{13}");
  frac_dnsk_vs_osc_params -> GetYaxis() -> SetTitle("#delta_{CP}");

  cout << "Generating predictions for sinsq2theta23 = " << input_sinsq2theta23 << ", dm23sq = " << input_dm23sq
       << " eV^2, sinsq2theta12 = " << input_sinsq2theta12 << ", dm12sq = " << input_dm12sq << " eV^2" << endl;

  // create & config superk spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;
  SuperkExpectation skexp;
  skexp.Read(config);

  // reference (2-flavour) physics hypothesis
  PhysHypothesis phref(t2kana::conv::kHypNuMuDisapp2F);
  phref.SetCurrValue(conv::kPhypSinsq2Theta23,  input_sinsq2theta23);
  phref.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq       );

  // prediction for reference hypothesis
  skexp.Generate(phref);
  double nskref = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
  cout << "Nsk[1Rmu,2-flavour] = " << nskref << endl;

  // 3-flavour physics hypothesis, setting fixed params
  PhysHypothesis ph(conv::kHypOsc3F);
  ph.SetCurrValue(conv::kPhypSinsq2Theta23,  input_sinsq2theta23);
  ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq       );
  ph.SetCurrValue(conv::kPhypSinsq2Theta12,  input_sinsq2theta12);
  ph.SetCurrValue(conv::kPhypDm12sq,         input_dm12sq       );
   
  // loop over varied params and generate predictions
  for(int i = 1; i <= nsk_vs_osc_params->GetNbinsX(); i++) {
    for(int j = 1; j <= nsk_vs_osc_params->GetNbinsY(); j++) {
       double sinsq2theta13 = nsk_vs_osc_params->GetXaxis()->GetBinCenter(i);
       double deltacp       = nsk_vs_osc_params->GetYaxis()->GetBinCenter(j);
       ph.SetCurrValue(conv::kPhypSinsq2Theta13,  sinsq2theta13);
       ph.SetCurrValue(conv::kPhypDeltaCP,        deltacp      );
       // generate predictions
       skexp.Generate(ph);
       // integrate Ereco distribution 
       double nsk = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
       double frac_dnsk = 100*(nsk-nskref)/nskref;
       // store values
       nsk_vs_osc_params->SetBinContent(i,j,nsk);
       frac_dnsk_vs_osc_params->SetBinContent(i,j,frac_dnsk);
       cout << "Nsk[1Rmu,3-flavour](sinsq2theta13 = " << sinsq2theta13 << ", deltacp = " << deltacp << ") = " << nsk 
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
     input_sinsq2theta12,-1,input_sinsq2theta23,-1,input_dm12sq,input_dm23sq,"no-theta13-deltacp");
  TFile f(Form("nsk_vs_sinsq2theta13_deltacp-%s.root",oscparams.c_str()),"recreate");
  nsk_vs_osc_params       -> Write("nsk");
  frac_dnsk_vs_osc_params -> Write("frac_dnsk");
  f.Close();
}
