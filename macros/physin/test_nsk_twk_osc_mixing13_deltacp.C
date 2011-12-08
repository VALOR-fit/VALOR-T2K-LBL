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

void test_nsk_twk_osc_mixing13_deltacp
(
  double       input_mixing23          = 1.00,                           /* */
  double       input_dm23sq            = 2.4E-3,                        /* */
  double       input_mixing12          = 0.8704,                           /* */
  double       input_dm12sq            = 7.6E-5,                        /* */
  const char * config_file             = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  const double pi = TMath::Pi();

  double mixing13_min       =  0.0;
  double mixing13_max       =  0.25;
  double mixing13_step      =  0.01;
  double deltacp_min        = -pi;
  double deltacp_max        =  pi;
  double deltacp_step       =  2*pi/90.;

  // make sure that min and max are at the bin centre
  int mixing13_nbins      = 1 + (int) ((mixing13_max - mixing13_min) / mixing13_step);
  int deltacp_nbins       = 1 + (int) ((deltacp_max       - deltacp_min      ) / deltacp_step      );

  TH2D * nsk_vs_osc_params = 
   new TH2D("nsk_vs_osc_params", "number of 1-ring e-like events",
       mixing13_nbins, 
       mixing13_min - 0.5*mixing13_step, 
       mixing13_max + 0.5*mixing13_step,
       deltacp_nbins, 
       deltacp_min - 0.5*deltacp_step, 
       deltacp_max + 0.5*deltacp_step);

  TH2D * frac_dnsk_vs_osc_params = 
   new TH2D("frac_dnsk_vs_osc_params", "% change in number of 1-ring e-like events wrt 2-flavour case",
       mixing13_nbins, 
       mixing13_min - 0.5*mixing13_step, 
       mixing13_max + 0.5*mixing13_step,
       deltacp_nbins, 
       deltacp_min - 0.5*deltacp_step, 
       deltacp_max + 0.5*deltacp_step);

  //lorena - added nsk histos for signal and background independently
  TH2D * nsk_vs_osc_params_sig = 
   new TH2D("nsk_vs_osc_params_sig", "number of 1-ring e-like events (SIGNAL)",
       mixing13_nbins, 
       mixing13_min - 0.5*mixing13_step, 
       mixing13_max + 0.5*mixing13_step,
       deltacp_nbins, 
       deltacp_min - 0.5*deltacp_step, 
       deltacp_max + 0.5*deltacp_step);

  TH2D * nsk_vs_osc_params_back = 
   new TH2D("nsk_vs_osc_params_back", "number of 1-ring e-like events (BACKGROUND)",
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
  nsk_vs_osc_params_sig   -> GetXaxis() -> SetTitle("sin^{2}2#theta_{13}");//lorena
  nsk_vs_osc_params_sig   -> GetYaxis() -> SetTitle("#delta_{CP}");
  nsk_vs_osc_params_back  -> GetXaxis() -> SetTitle("sin^{2}2#theta_{13}");
  nsk_vs_osc_params_back  -> GetYaxis() -> SetTitle("#delta_{CP}");


  cout << "Generating predictions for mixing23 = " << input_mixing23 << ", dm23sq = " << input_dm23sq
       << " eV^2, mixing12 = " << input_mixing12 << ", dm12sq = " << input_dm12sq << " eV^2" << endl;

  // create & config superk spectrum predictor
  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

 // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  SuperkExpectation skexp;
  skexp.Read(config);

  // reference (2-flavour) physics hypothesis
  PhysHypothesis phref(t2kana::conv::kHypNuMuDisapp2F);
  phref.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
  phref.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq       );

  // prediction for reference hypothesis
  skexp.Generate(phref);
  double nskref = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
  cout << "Nsk[1Rmu,2-flavour] = " << nskref << endl;

  // 3-flavour physics hypothesis, setting fixed params
  PhysHypothesis ph(conv::kHypOsc3F);
  ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
  ph.SetCurrValue(conv::kPhypDm23sq,         input_dm23sq       );
  ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
  ph.SetCurrValue(conv::kPhypDm12sq,         input_dm12sq       );

  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << "- Physics hypothesis: " << ph << endl;
  cout << *Opt::Instance() << endl;


  
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
       double nsk = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingELike)->Integral();
       double frac_dnsk = 100*(nsk-nskref)/nskref;
       //lorena - now separate signal and background
       double nsksig = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingELike, conv::kSKmSigNuECCQE)->Integral();
       nsksig += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingELike, conv::kSKmSigNuECCoth)->Integral();
       nsksig += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingELike, conv::kSKmSigNuECC1pi)->Integral();
       double nskback = nsk - nsksig;
       // store values
       nsk_vs_osc_params->SetBinContent(i,j,nsk);
       nsk_vs_osc_params_sig->SetBinContent(i,j,nsksig);//lorena
       nsk_vs_osc_params_back->SetBinContent(i,j,nskback);
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
  TCanvas * c3 = new TCanvas;
  nsk_vs_osc_params_sig       -> Draw("colz");
  TCanvas * c4 = new TCanvas;
  nsk_vs_osc_params_back       -> Draw("colz");


  // save
  string oscparams = util::AsString(ph.Hypothesis(),   
		     input_mixing12,-1, input_mixing23,-1, input_dm12sq,input_dm23sq);
  TFile f(Form("nsk_vs_mixing13_deltacp-%s.root",oscparams.c_str()),"recreate");
  nsk_vs_osc_params       -> Write("nsk");
  frac_dnsk_vs_osc_params -> Write("frac_dnsk");
  nsk_vs_osc_params_sig       -> Write("nsk_sig");
  nsk_vs_osc_params_back       -> Write("nsk_back");
  //f.Close();//lorena


  //____________________________________

  //lorena - now repeat the analysis changing deltasqm23 and fixing deltaCP (may go in different place...)
  mixing13_min  =  0.0;
  mixing13_max  =  0.25;
  mixing13_step =  0.01;
  double dm23sq_min        = 1.0E-3;
  double dm23sq_max        = 4.0E-3;
  double dm23sq_step       =  0.1E-3;

  // make sure that min and max are at the bin centre
  int mixing13_nbins = 1 + (int) ((mixing13_max - mixing13_min) / mixing13_step);
  int dm23sq_nbins       = 1 + (int) ((dm23sq_max       - dm23sq_min      ) / dm23sq_step      );

  TH2D * nsk_vs_osc_params_bis = 
   new TH2D("nsk_vs_osc_params_bis", "number of 1-ring e-like events",
       mixing13_nbins, 
       mixing13_min - 0.5*mixing13_step, 
       mixing13_max + 0.5*mixing13_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);

  TH2D * frac_dnsk_vs_osc_params_bis = 
   new TH2D("frac_dnsk_vs_osc_params_bis", "% change in number of 1-ring e-like events wrt 2-flavour case",
       mixing13_nbins, 
       mixing13_min - 0.5*mixing13_step, 
       mixing13_max + 0.5*mixing13_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);

  //lorena - added nsk histos for signal and background independently
  TH2D * nsk_vs_osc_params_bis_sig = 
   new TH2D("nsk_vs_osc_params_bis_sig", "number of 1-ring e-like events (SIGNAL)",
       mixing13_nbins, 
       mixing13_min - 0.5*mixing13_step, 
       mixing13_max + 0.5*mixing13_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);
  TH2D * nsk_vs_osc_params_bis_back = 
   new TH2D("nsk_vs_osc_params_bis_back", "number of 1-ring e-like events (BACKGROUND)",
       mixing13_nbins, 
       mixing13_min - 0.5*mixing13_step, 
       mixing13_max + 0.5*mixing13_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);


  nsk_vs_osc_params_bis       -> GetXaxis() -> SetTitle("sin^{2}2#theta_{13}");
  nsk_vs_osc_params_bis       -> GetYaxis() -> SetTitle("#Delta m^{2}_{23}");
  frac_dnsk_vs_osc_params_bis -> GetXaxis() -> SetTitle("sin^{2}2#theta_{13}");
  frac_dnsk_vs_osc_params_bis -> GetYaxis() -> SetTitle("#Delta m^{2}_{23}");
  nsk_vs_osc_params_bis_sig   -> GetXaxis() -> SetTitle("sin^{2}2#theta_{13}");//lorena
  nsk_vs_osc_params_bis_sig   -> GetYaxis() -> SetTitle("#Delta m^{2}_{23}");
  nsk_vs_osc_params_bis_back  -> GetXaxis() -> SetTitle("sin^{2}2#theta_{13}");
  nsk_vs_osc_params_bis_back  -> GetYaxis() -> SetTitle("#Delta m^{2}_{23}");


  cout << "Generating predictions for mixing23 = " << input_mixing23 << ", dm23sq = " << input_dm23sq
       << " eV^2, mixing12 = " << input_mixing12 << ", dm12sq = " << input_dm12sq << " eV^2" << endl;

  
  
  // prediction for reference hypothesis
  skexp.Generate(phref);
  double nskref = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike)->Integral();
  cout << "Nsk[1Rmu,2-flavour] = " << nskref << endl;

  double input_deltaCP   = 0.0 ;

  // 3-flavour physics hypothesis, setting fixed params
  // PhysHypothesis ph(conv::kHypOsc3F);
  ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
  ph.SetCurrValue(conv::kPhypDeltaCP,        input_deltaCP       );
  ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
  ph.SetCurrValue(conv::kPhypDm12sq,         input_dm12sq       );
   

  // loop over varied params and generate predictions
  for(int i = 1; i <= nsk_vs_osc_params_bis->GetNbinsX(); i++) {
    for(int j = 1; j <= nsk_vs_osc_params_bis->GetNbinsY(); j++) {
       double mixing13 = nsk_vs_osc_params_bis->GetXaxis()->GetBinCenter(i);
       double deltam23sq       = nsk_vs_osc_params_bis->GetYaxis()->GetBinCenter(j);
       ph.SetCurrValue(conv::kPhypMixing13,  mixing13);
       ph.SetCurrValue(conv::kPhypDm23sq ,      deltam23sq    );
       // generate predictions
       skexp.Generate(ph);
       // integrate Ereco distribution 
       double nsk = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingELike)->Integral();
       double frac_dnsk = 100*(nsk-nskref)/nskref;
       //lorena - now separate signal and background
       double nsksig = skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingELike, conv::kSKmSigNuECCQE)->Integral();
       nsksig += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingELike, conv::kSKmSigNuECCoth)->Integral();
       nsksig += skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingELike, conv::kSKmSigNuECC1pi)->Integral();
       double nskback = nsk - nsksig;
       // store values
       nsk_vs_osc_params_bis->SetBinContent(i,j,nsk);
       nsk_vs_osc_params_bis_sig->SetBinContent(i,j,nsksig);//lorena
       nsk_vs_osc_params_bis_back->SetBinContent(i,j,nskback);
       frac_dnsk_vs_osc_params_bis->SetBinContent(i,j,frac_dnsk);
       cout << "Nsk[1Rmu,3-flavour](mixing13 = " << mixing13 << ", deltam23sq  = " << deltam23sq << ") = " << nsk 
            << ", change from the 2-flavour hypothesis = " << frac_dnsk << "%" << endl;
    }//j
  }//i

  // draw
  gStyle->SetOptStat(0);
  TCanvas * c5 = new TCanvas;
  nsk_vs_osc_params_bis       -> Draw("colz");
  TCanvas * c6 = new TCanvas;
  frac_dnsk_vs_osc_params_bis -> Draw("colz");
  TCanvas * c7 = new TCanvas;
  nsk_vs_osc_params_bis_sig       -> Draw("colz");
  TCanvas * c8 = new TCanvas;
  nsk_vs_osc_params_bis_back       -> Draw("colz");

  // save
  //string oscparams = util::AsString("3f_osc",   
  //input_mixing12,-1,input_mixing23,-1,input_dm12sq,input_dm23sq,"no-theta13-deltacp");
  //  TFile f(Form("nsk_vs_mixing13_deltacp-%s.root",oscparams.c_str()),"recreate");
  nsk_vs_osc_params_bis       -> Write("nsk_bis");
  frac_dnsk_vs_osc_params_bis -> Write("frac_dnsk_bis");
  nsk_vs_osc_params_sig_bis       -> Write("nsk_bis_sig");
  nsk_vs_osc_params_back_bis       -> Write("nsk_bis_back");
  f.Close();


}
