//
// ------------------------------------------------------------------------------------------------
// Plot number of events predicted at SuperK (NSK) vs oscillation params
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

//#include "Conventions.h"

void test_nsk_twk_osc
(
  bool         thrownon23  = false,
  const int    nthrows     = 100,
  const char * hypothesis  = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  const char * mixparamchoice       = "sinsqtheta",
  const char * hierarchy            = "normal",
  const char * config_file = "analysis_configuration.xml"
)
{
  using namespace std;
  using namespace valor;

  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  double mixing23_min  =  0.00;
  double mixing23_max  =  1.00;
  double mixing23_step =  0.005;
  double dm23sq_min    =  1.0E-3;
  double dm23sq_max    =  4.0E-3;
  double dm23sq_step   =  5E-5;

 // make sure that min and max are at the bin centre
  int mixing23_nbins = 1 + (int)( ((mixing23_max - mixing23_min) / mixing23_step) + 0.01);
  int dm23sq_nbins   = 1 + (int)( ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  ) + 0.01);

 TH2D * nsk_vs_osc_params = 
   new TH2D("nsk_vs_osc_params", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // fit params
  PhysHypothesis ph(hypothesis);

  // generate predictions
  for(int i = 1; i <= nsk_vs_osc_params->GetNbinsX(); i++) {
    for(int j = 1; j <= nsk_vs_osc_params->GetNbinsY(); j++) {
       double mixing23 = nsk_vs_osc_params->GetXaxis()->GetBinCenter(i);
       double dm23sq   = nsk_vs_osc_params->GetYaxis()->GetBinCenter(j);
       ph.SetCurrValue(valor::conv::kPhypMixing23, mixing23);
       ph.SetCurrValue(valor::conv::kPhypDm23sq,   dm23sq);
       if(ph.Hypothesis() == conv::kHypOsc3F) {
	 if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
	   ph.SetCurrValue(conv::kPhypMixing12,  conv::kGlobFitSinsq2Theta12);
	   ph.SetCurrValue(conv::kPhypMixing13,  conv::kGlobFitSinsq2Theta13);
	 }
	 else {
	   ph.SetCurrValue(conv::kPhypMixing12,  conv::kGlobFitSinsqTheta12);
	   ph.SetCurrValue(conv::kPhypMixing13,  conv::kGlobFitSinsqTheta13);
	 }
	 ph.SetCurrValue(conv::kPhypCPV    ,   conv::kGlobFitDeltaCP);
	 ph.SetCurrValue(conv::kPhypDm12sq,    conv::kGlobFitDm12sq);
       }

       // generate predictions
       skexp.Generate(ph);

       // integrate Ereco distribution
       double nsk = skexp.CurrPredictedErecoSpectrum(s)->Integral();

       if(thrownon23) {
	 double nskmin = 99999;
	 double nskmax = 0;
	 for(int t = 0; t < nthrows; t++) {
	   ph.SetCurrValue(conv::kPhypMixing12,  util::RandomizeMixing12());
	   ph.SetCurrValue(conv::kPhypMixing13,  util::RandomizeMixing13());
	   ph.SetCurrValue(conv::kPhypCPV    ,   util::RandomizeCPV());
	   ph.SetCurrValue(conv::kPhypDm12sq,    util::RandomizeDm12sq());
	   skexp.Generate(ph);
	   double thisnsk = skexp.CurrPredictedErecoSpectrum(s)->Integral();
	   //nsksq += pow(thisnsk - nsk, 2);
	   nskmin = TMath::Min(thisnsk, nskmin);
	   nskmax = TMath::Max(thisnsk, nskmax);
	 }
	 nsk = (TMath::Abs(nsk-nskmin) > TMath::Abs(nsk-nskmax)) ? nskmin : nskmax;
       }

       nsk_vs_osc_params->SetBinContent(i,j,nsk);
    }//j
    cout << "x bin " << i << " of " << nsk_vs_osc_params->GetNbinsX() << endl;
  }//i

  gStyle->SetOptStat(0);

  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)
    nsk_vs_osc_params->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
  else
    nsk_vs_osc_params->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  nsk_vs_osc_params->GetYaxis()->SetTitle("|#Deltam^{2}_{32}| (eV^{2}/c^{4})");

  nsk_vs_osc_params->Draw("colz");

  TFile f("./nsk_vs_osc_params.root","recreate");
  nsk_vs_osc_params->Write("nsk_vs_osc_params");
  f.Close();
}
