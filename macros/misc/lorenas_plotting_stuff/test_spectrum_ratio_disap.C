//
// ------------------------------------------------------------------------------------------------
// Plot ratio at minimum for disappearance spectrum at different values of theta23
//
// Lorena Escudero
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

void test_spectrum_ratio_disap
(
  const char * hypothesis           = "3f_osc",                       /* options: "2f_numu_disap", "3f_osc" */
  double       input_mixing12       = 0.306,                            /* */
  double       input_mixing13       = 0.0251,                            /* */
  double       input_mixing23       = 0.5,                            /* */
  double       input_deltacp        = 0,
  double       input_dm23sq         = 0.0024375,                            /* */
  double       input_dm12sq         = 7.5E-5,                            /* */
  const int    nbins_mixing23       = 41,
  double       min_mixing23         = 0.3,
  double       max_mixing23         = 0.7,
  const char * config_file          = "analysis_configuration_Run4a.xml"  /* */
)
{
  using namespace std;
  using namespace t2kana;
  
  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsqTheta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;
  //  Opt::Instance()->gOscAnaTypeChoice     = conv::kAnaNuMuDisappOnly;
  Opt::Instance()->gOscDmsqScheme        = conv::kNuDmsqFogliLisi;
  
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
  if(input_dm23sq < -1)  {
    conv::NuOscDmsqScheme_t odmsq = Opt::Instance()->gOscDmsqScheme;                                                                                                                                    
    switch (odmsq) {                                                                                                                                                                                        
    case (conv::kNuDmsqFogliLisi): 
      input_dm23sq  =  conv::kGlobFitDmsqFL;
      break;                                                                                                                                                                                                
    case (conv::kNuDmsq12andAbs23):            
      input_dm23sq = conv::kGlobFitDm23sq;
      break;                                                                                                                                                                                                
    }                       
  }
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
  
  //first find the spectrum with no oscillations
  ph.SetCurrValue(conv::kPhypMixing12,  0.0);
  ph.SetCurrValue(conv::kPhypMixing13,  0.0);
  ph.SetCurrValue(conv::kPhypMixing23,  0.0);
  ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp);
  ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
  ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  
  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  cout << "Creating skexp" << endl;
  skexp->Read(config);
  skexp->Generate(ph);
  
  //output file
  const char * outfile0 = "disap_spectum_ratio.root";
  TFile f0(outfile0,"recreate");
  
  TH1D * hErecoTotNoOsc[conv::kns];
  
  for(int is=0; is<conv::kns; is++) {
    conv::SuperkSample_t s = conv::kSKSamples[is];
    if(s==conv::kSKmUndefined) continue;
    const int nbins = util::NumErecoBins(s);
    const double * ereco_bin_edges = util::ErecoBinEdges(s);
    hErecoTotNoOsc[is] = new TH1D("","",nbins,ereco_bin_edges);
    for(int ir = 0; ir < nbins; ir++){
      double nev = skexp->CurrPredictedErecoSpectrumArray(s)[ir];
      double energy = 0.5*(util::ErecoBinEdges(s,ir) + util::ErecoBinEdges(s,ir+1));
      hErecoTotNoOsc[is]->Fill(energy,nev);
    }
    
    //    cout << "Total =" << hErecoTotNoOsc[is]->Integral() << " events. " << endl;
    
    f0.cd();
    hErecoTotNoOsc[is] -> Write(Form("hErecoTotNoOsc_%s",util::AsStringCompact(s).c_str()));
  } //end of first loop over samples, to save the spectrum with no oscillations
  
  
  TH1D * hErecoTot[conv::kns][nbins_mixing23];
  TH1D * hErecoTotRatio[conv::kns][nbins_mixing23];
  
  TH1D * min_dist = new TH1D("min_dist","Ratio at minimum",nbins_mixing23,0.3,0.7);
  TH1D * events_dist = new TH1D("events_dist","Number of events",nbins_mixing23,0.3,0.7);
  TH1D * min_energy_dist = new TH1D("min_energy_dist","Energy at minimum",nbins_mixing23,0.3,0.7);
  
  double step_mixing23 = (max_mixing23-min_mixing23)/nbins_mixing23;
  
  for(int index = 0; index<nbins_mixing23; index++){
    
    double mix23 = min_mixing23 + (index)*step_mixing23; 
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  mix23);
    ph.SetCurrValue(conv::kPhypDeltaCP,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    
    skexp->Generate(ph);
    
    
    //loop over samples
    for(int is=0; is<conv::kns; is++) {
      conv::SuperkSample_t s = conv::kSKSamples[is];
      if(s==conv::kSKmUndefined) continue;
      const int nbins = util::NumErecoBins(s);
      const double * ereco_bin_edges = util::ErecoBinEdges(s);
      hErecoTot[is][index] = new TH1D("","",nbins,ereco_bin_edges);
      for(int ir = 0; ir < nbins; ir++){
	double nev = skexp->CurrPredictedErecoSpectrumArray(s)[ir];
	double energy = 0.5*(util::ErecoBinEdges(s,ir) + util::ErecoBinEdges(s,ir+1));
	hErecoTot[is][index]->Fill(energy,nev);
      }
      
      //cout << "Total =" << hErecoTot[is][index]->Integral() << " events. " << endl;
      
      hErecoTotRatio[is][index] = new TH1D(*hErecoTot[is][index]);
      hErecoTotRatio[is][index]->Divide(hErecoTotNoOsc[is]);
      
      double min = 1;
      double min_energy = 0;
      
      for(int ix = 1; ix <= nbins; ix++){
	double content = hErecoTotRatio[is][index]->GetBinContent(ix);
	if((content < min)&&(content != 0.0)){
	  min = content;
	  min_energy =  hErecoTotRatio[is][index]->GetBinCenter(ix);
	}
      }
      
      min_dist -> SetBinContent(index+1,min);
      min_energy_dist -> SetBinContent(index+1,min_energy);
      events_dist->SetBinContent(index+1,hErecoTot[is][index]->Integral());
      
      f0.cd();
      hErecoTot[is][index] -> Write(Form("hErecoTot_%s_%d",util::AsStringCompact(s).c_str(),index));
      hErecoTotRatio[is][index] -> Write(Form("hErecoTotRatio_%s_%d",util::AsStringCompact(s).c_str(),index));
    } //end of first loop over samples, to save the spectrum histograms
    
    //  delete skexp;
    
  }
  
  f0.cd();
  min_dist->GetYaxis()->SetTitle("Ratio at minimum");
  min_dist->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  min_dist->Write("min_dist");
  min_energy_dist->Write("min_energy_dist");
  events_dist -> Write("events_dist");
  
  
}
