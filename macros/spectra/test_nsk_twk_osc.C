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
  bool         thrownon23  = true,
  const int    nthrows     = 10,
  const char * sk_sample   = "RHC1Re",          /* options: "1Rmu", "1Re", "MultiRmu" */
  const char * hypothesis  = "3f_osc",                 /* options: "2f_numu_disap", "3f_osc", "3f_oscbar" */
  const char * mixparamchoice       = "sinsqtheta",
  const char * hierarchy            = "normal",
  const char * config_file = "analysis_configuration_run1to6_test.xml",

  double       input_mixing12       = -99,                            /* */
  double       input_mixing13       = -99,                            /* */
  double       input_mixing23       = -99,                            /* */
  double       input_deltacp        = -99,                            /* */
  double       input_dm12sq         = -99,                            /* */
  double       input_dm23sq         = -99,                            /* */
)
{
  using namespace std;
  using namespace valor;

  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample

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
  double dm23sq_min    =  2.00001E-4;
  double dm23sq_max    =  3.99999E-2;
  double dm23sq_step   =  5E-4;
  // make sure that min and max are at the bin centre
  int mixing23_nbins = 1 + (int)( ((mixing23_max - mixing23_min) / mixing23_step) + 0.01);
  int dm23sq_nbins   = 1 + (int)( ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  ) + 0.01);


  double mixing13_min  =  0.00;
  double mixing13_max  =  0.20;
  double mixing13_step =  0.001;
  //double mixing13_step =  0.05;
  double dcp_min       =  -conv::kpi;
  double dcp_max       =  conv::kpi;
  double dcp_step      =  conv::kpi / 100;
  
  // make sure that min and max are at the bin centre
  int mixing13_nbins = 1 + (int)( ((mixing13_max - mixing13_min) / mixing13_step) + 0.01);
  int dcp_nbins   = 1 + (int)( ((dcp_max   - dcp_min  ) / dcp_step  ) + 0.01);

  TH2D * nsk_vs_osc_params;
  TString sk_sample_str = sk_sample;

  if(sk_sample_str.Contains("1Rmu")){

    cout << endl;
    cout << "Make plot mixing23 Vs dm23sq" << endl;
    cout << endl;

    nsk_vs_osc_params = 
      new TH2D("nsk_vs_osc_params", "",
	       mixing23_nbins, 
	       mixing23_min - 0.5*mixing23_step, 
	       mixing23_max + 0.5*mixing23_step,
	       dm23sq_nbins, 
	       dm23sq_min - 0.5*dm23sq_step, 
	       dm23sq_max + 0.5*dm23sq_step);
  }
  if(sk_sample_str.Contains("1Re")){

    cout << endl;
    cout << "Make plot mixing13 Vs dcp" << endl;
    cout << endl;

    nsk_vs_osc_params = 
      new TH2D("nsk_vs_osc_params", "",
	       mixing13_nbins, 
	       mixing13_min - 0.5*mixing13_step, 
	       mixing13_max + 0.5*mixing13_step,
	       dcp_nbins, 
	       dcp_min - 0.5*dcp_step, 
	       dcp_max + 0.5*dcp_step);
  }

  // superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // fit params
  PhysHypothesis ph(hypothesis);

  // generate predictions
  for(int i = 1; i <= nsk_vs_osc_params->GetNbinsX(); i++) {
    for(int j = 1; j <= nsk_vs_osc_params->GetNbinsY(); j++) {
      
      // double mixing23 = nsk_vs_osc_params->GetXaxis()->GetBinCenter(i);
      // double dm23sq   = nsk_vs_osc_params->GetYaxis()->GetBinCenter(j);
      double par1 = nsk_vs_osc_params->GetXaxis()->GetBinCenter(i);
      double par2 = nsk_vs_osc_params->GetYaxis()->GetBinCenter(j);
      	
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
	 ph.SetCurrValue(valor::conv::kPhypMixing23,  input_mixing23);
	 ph.SetCurrValue(valor::conv::kPhypDm23sq,    input_dm23sq);
       } else
	 if(ph.Hypothesis() == conv::kHypOsc3F) {

	   ph.SetCurrValue(valor::conv::kPhypMixing12,  input_mixing12);
	   ph.SetCurrValue(valor::conv::kPhypDm12sq,    input_dm12sq);

	   if(sk_sample_str.Contains("1Rmu")){	   
	     ph.SetCurrValue(valor::conv::kPhypMixing13,  input_mixing13);
	     ph.SetCurrValue(valor::conv::kPhypCPV    ,   input_deltacp);
	     ph.SetCurrValue(valor::conv::kPhypMixing23, par1);
	     ph.SetCurrValue(valor::conv::kPhypDm23sq,   par2);
	   }
	   else if(sk_sample_str.Contains("1Re")){
	     ph.SetCurrValue(valor::conv::kPhypMixing13,  par1);
	     ph.SetCurrValue(valor::conv::kPhypCPV    ,   par2);
	     ph.SetCurrValue(valor::conv::kPhypMixing23, input_mixing23);
	     ph.SetCurrValue(valor::conv::kPhypDm23sq,   input_dm23sq);
	   }

	   
	 } else
	   if(ph.Hypothesis() == conv::kHypOsc3FBar){
	     ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
	     ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
	     ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
	     ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
	     ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
	     ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);

	     ph.SetCurrValue(conv::kPhypMixing12Bar,  input_mixing12);
	     ph.SetCurrValue(conv::kPhypDm12sqBar,    input_dm12sq);

	     if(sk_sample_str.Contains("1Rmu")){	   
	       ph.SetCurrValue(conv::kPhypMixing13Bar,  input_mixing13);
	       ph.SetCurrValue(conv::kPhypCPVBar    ,   input_deltacp);
	       ph.SetCurrValue(valor::conv::kPhypMixing23Bar, par1);
	       ph.SetCurrValue(valor::conv::kPhypDm23sqBar,   par2);
	     }
	     else if(sk_sample_str.Contains("1Re")){	   
	       ph.SetCurrValue(conv::kPhypMixing13Bar,  par1);
	       ph.SetCurrValue(conv::kPhypCPVBar    ,   par2);
	       ph.SetCurrValue(valor::conv::kPhypMixing23Bar, input_mixing23);
	       ph.SetCurrValue(valor::conv::kPhypDm23sqBar,   input_dm23sq);
	     }	  
   
	   }


       // generate predictions
       skexp.Generate(ph);

       // integrate Ereco distribution
       double nsk = skexp.CurrPredictedErecoSpectrum(s)->Integral();

       if(thrownon23) {
	 double nskmin = 99999;
	 double nskmax = 0;
	 for(int t = 0; t < nthrows; t++) {

	   if(ph.Hypothesis() == conv::kHypOsc3F) {
	     ph.SetCurrValue(conv::kPhypMixing12,  util::RandomizeMixing12());
	     ph.SetCurrValue(conv::kPhypDm12sq,    util::RandomizeDm12sq());
	     
	     if(sk_sample_str.Contains("1Rmu")){
	       ph.SetCurrValue(conv::kPhypMixing13,  util::RandomizeMixing13());
	       ph.SetCurrValue(conv::kPhypCPV     ,  util::RandomizeCPV());
	     }
	     else if(sk_sample_str.Contains("1Re")){
	       ph.SetCurrValue(conv::kPhypMixing23,  util::RandomizeMixing23());
	       ph.SetCurrValue(conv::kPhypDm23sq  ,  util::RandomizeDm23sq());
	     }
	   }
	   else
	    if(ph.Hypothesis() == conv::kHypOsc3FBar){
	      ph.SetCurrValue(conv::kPhypMixing12Bar,  util::RandomizeMixing12Bar());
	      ph.SetCurrValue(conv::kPhypDm12sq,    util::RandomizeDm12sqBar());
	      
	      if(sk_sample_str.Contains("1Rmu")){
		ph.SetCurrValue(conv::kPhypMixing13Bar,  util::RandomizeMixing13Bar());
		ph.SetCurrValue(conv::kPhypCPVBar     ,  util::RandomizeCPVBar());
	      }
	      else if(sk_sample_str.Contains("1Re")){
		ph.SetCurrValue(conv::kPhypMixing23Bar,  util::RandomizeMixing23Bar());
		ph.SetCurrValue(conv::kPhypDm23sqBar  ,  util::RandomizeDm23sqBar());
	      }
	    }
	   
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

  if(sk_sample_str.Contains("1Rmu")){
    if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta)
      nsk_vs_osc_params->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
    else
      nsk_vs_osc_params->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
    nsk_vs_osc_params->GetYaxis()->SetTitle("|#Deltam^{2}_{32}| (eV^{2}/c^{4})");
  }
  else if(sk_sample_str.Contains("1Re")){
    nsk_vs_osc_params->GetXaxis()->SetTitle("sin^{2}#theta_{13}");
    nsk_vs_osc_params->GetYaxis()->SetTitle("#Delta_{cp}");
  }


  nsk_vs_osc_params->Draw("colz");
    
  TFile f("./nsk_vs_osc_params.root","recreate");
  nsk_vs_osc_params->Write("nsk_vs_osc_params");
  f.Close();
}
