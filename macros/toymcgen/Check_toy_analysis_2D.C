
using namespace std;
using namespace valor;
using namespace valor::conv;

void Check_toy_analysis_2D
(
 const char * sk_sample = "RHC1Re",   /* "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */
  
 int first_experiment  = 0, //0, //4404, //4739,                       /* */
 int last_experiment   = 80, //9999, //4404, //4739,                       /* */
  
 const char * xmlfile = "",
 const char * rootfile = "Toys_fit_submit_fitmarg_joint_dcp_nexp_100_2.000_-999.000-job96-_0_99.root",

 const char * histoname = "mock_data_RHC1Re_%d",
 //const char * histoname = "mock_data_RHC1Rmu_%d",
 //const char * histoname = "mock_data_1Re_%d",
 //const char * histoname = "mock_data_1Rmu_%d",
 
 bool using_xml_file = false,

 // Default
 double       input_mixing12           = -99,           
 double       input_mixing13           = -99, 
 double       input_mixing23           = -99, 
 double       input_deltacp            = -99, 
 double       input_dm12sq             = -99,
 double       input_dm23sq             = -99,
 double       input_beta_bar           = 1,

 //const char * config_file      = "analysis_configuration_run1to6_joint_Nu6.914E20POT_NuBar4.0108E20POT.xml",
 //const char * config_file      = "analysis_configuration_run1to4_joint_Nu6.57E20POT.xml",
 //const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA_Nu6.57E20POT.xml",
 const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA.xml",
 //const char * config_file      = "analysis_configuration_test_theta_elike_RFGRPA_FullStat.xml",

 const char * fit_param_list        = "nuosc3f+nuisance", // dcp

 const char * dmsqscheme               = "auto",   /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * mixparamchoice           = "sinsqtheta", /* sinsq2theta or sinsqtheta */
 const char * hierarchy                = "normal",     /* normal or inverted */
 const char * mixparamoctant           = "two",       /* one==kNuOctantOne two==kNuOctantTwo */

 bool use_new_definition_of_dn_over_n = true
 ){

  bool RHCFlag = false;
  TString str_sk_sample(sk_sample);
  if(str_sk_sample.Contains("RHC"))
    RHCFlag = true;

  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }

  // Analysis configuration
  TString config  (Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file));
  TString xmlfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile));

  FileStat_t dummy;

  // Check for config file, and XML data file.
  // If not found in inputs or config, see if it works as a normal path
  if (gSystem->GetPathInfo(config.Data(), dummy)) 
    {
      config = config_file;
      if (gSystem->GetPathInfo(config.Data(), dummy)) {
	cerr<<"Config file not found: "<<config_file<<endl;
	exit(-1);
      }
    }
  
  // Read toy files (xml or root file)
  TFile *frootfile;
  if ( using_xml_file && gSystem->GetPathInfo(xmlfile2.Data(), dummy)) 
    {
      xmlfile2 = xmlfile;
      if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) {
	cerr<<"Data input file not found: "<<xmlfile<<endl;
	exit(-1);
      }
    }
  else{
    frootfile = new TFile(rootfile,"READ");
  }
  
  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy);
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(dmsqscheme);
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(mixparamchoice);
  Opt::Instance()->gOscOctantChoice      = util::OctantFromString(mixparamoctant);
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }
  Opt::Instance()->SetDefaultParameterValues();

  //
  
  //Set to defaults
  if(Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsq2Theta) {
    if(input_mixing13 < -1 && input_mixing13 > -900) input_mixing13 = conv::kGlobFitSinsq2Theta13;
    if(input_mixing23 < -1 && input_mixing23 > -900) input_mixing23 = conv::kGlobFitSinsq2Theta23;
    if(input_mixing12 < -1 && input_mixing12 > -900) input_mixing12 = conv::kGlobFitSinsq2Theta12;
  }
  else {
    if(input_mixing13 < -1 && input_mixing13 > -900) input_mixing13 = conv::kGlobFitSinsqTheta13;
    if(input_mixing23 < -1 && input_mixing23 > -900) input_mixing23 = conv::kGlobFitSinsqTheta23;
    if(input_mixing12 < -1 && input_mixing12 > -900) input_mixing12 = conv::kGlobFitSinsqTheta12;
  }
  if(input_dm23sq < -1 && input_dm23sq > -900)  input_dm23sq = conv::kGlobFitDm23sq;
  if(input_dm12sq < -1 && input_dm12sq > -900)  input_dm12sq = conv::kGlobFitDm12sq;
  if(input_deltacp < -90 && input_deltacp > -900) input_deltacp = conv::kGlobFitDeltaCP;

  if(input_beta_bar < -90 && input_beta_bar > -900) input_beta_bar = conv::kGlobFitBetaBar;


  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());
 
  // Fit param list for fakedata sets
  FitParamList * fparm  = util::NamedParameterList(fit_param_list); 
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used for fakedata sets: " << util::AsString(hyp) << endl;

  // set oscillation parameter inputs as necessary (used to generate toy-MC expt)
  if(hyp == kHypOsc3F) {
    fparm->SetInitValue(kPhypMixing12, input_mixing12);
    fparm->SetInitValue(kPhypMixing13, input_mixing13);
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
    fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
  }
  else if(hyp == kHypOsc3FBar) {
    fparm->SetInitValue(kPhypMixing12, input_mixing12);
    fparm->SetInitValue(kPhypMixing13, input_mixing13);
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
    fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
    
    fparm->SetInitValue(kPhypMixing12Bar, input_mixing12);
    fparm->SetInitValue(kPhypMixing13Bar, input_mixing13);
    fparm->SetInitValue(kPhypMixing23Bar, input_mixing23);
    fparm->SetInitValue(kPhypCPVBar    ,  input_deltacp);
    fparm->SetInitValue(kPhypDm12sqBar,   input_dm12sq);
    fparm->SetInitValue(kPhypDm23sqBar,   input_dm23sq);
  }
  else if(hyp == kHypOsc3FBarNormBar) {
    fparm->SetInitValue(kPhypMixing12, input_mixing12);
    fparm->SetInitValue(kPhypMixing13, input_mixing13);
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
    fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
    
    fparm->SetInitValue(kPhypMixing12Bar, input_mixing12);
    fparm->SetInitValue(kPhypMixing13Bar, input_mixing13);
    fparm->SetInitValue(kPhypMixing23Bar, input_mixing23);
    fparm->SetInitValue(kPhypCPVBar    ,  input_deltacp);
    fparm->SetInitValue(kPhypDm12sqBar,   input_dm12sq);
    fparm->SetInitValue(kPhypDm23sqBar,   input_dm23sq);

    fparm->SetInitValue(kPhypBetaBar,     input_beta_bar);
  }
  else if(hyp == kHypOsc3FNormBar) {
    fparm->SetInitValue(kPhypMixing12, input_mixing12);
    fparm->SetInitValue(kPhypMixing13, input_mixing13);
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypCPV    ,  input_deltacp);
    fparm->SetInitValue(kPhypDm12sq,   input_dm12sq);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
    fparm->SetInitValue(kPhypBetaBar,  input_beta_bar);
  }
  else if(hyp == kHypNuMuDisapp2F) {
    fparm->SetInitValue(kPhypMixing23, input_mixing23);
    fparm->SetInitValue(kPhypDm23sq,   input_dm23sq);
  }

  // Get physics hypothesis

  PhysHypothesis ph(*fparm);
  
  // Get the nominal

  skexp->Generate(ph);
  // TH1D *hErecoNominal = skexp->CurrPredictedErecoSpectrum(s);
  // hErecoNominal->SetDirectory(0);
  // hErecoNominal->GetXaxis()->SetTitle("E_{reco} (GeV)");
  // hErecoNominal->GetYaxis()->SetTitle("Events per bin");
  // hErecoNominal->SetLineWidth(2);
  TH2D *hErecoVsQ2Nominal = skexp->CurrPredictedErecoVsQ2Spectrum(s);
  hErecoVsQ2Nominal->SetDirectory(0);
  hErecoVsQ2Nominal->GetYaxis()->SetTitle("E_{reco} (GeV)");
  hErecoVsQ2Nominal->GetXaxis()->SetTitle("Q^{2}");
  hErecoVsQ2Nominal->GetZaxis()->SetTitle("Events per bin");
  hErecoVsQ2Nominal->SetLineWidth(2);
   
  double nominal_total_nsk = hErecoVsQ2Nominal->Integral();

  //

  int nexperiments = last_experiment - first_experiment + 1; 
  int n_ereco_bins = util::NumErecoBins(s);  
  double * sumsq_per_bin = (double*)calloc(n_ereco_bins,sizeof(double));
  VALORArray * nsk_per_toy = new VALORArray(nexperiments);
  double * nsk_per_toy_arr = nsk_per_toy->GetArray();

  TH1D * hNsk = new TH1D("hNsk","hNsk",10000,0,1000);

  //

  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;

  cout << "- Toy-MC" << endl;
  cout << "   - Reading experiments from: " << xmlfile2 << endl; 
  cout << "   - Using experiments in range: [" << first_experiment << ", " << last_experiment << "]" << endl; 
  cout << "   - Number of toy experiments: " << nexperiments << endl; 
  cout << "   - True oscillation parameters: " << endl;
  cout << "       > " << util::AsString(kPhypMixing12) << " = " << input_mixing12 << endl;
  cout << "       > " << util::AsString(kPhypMixing13) << " = " << input_mixing13 << endl;
  cout << "       > " << util::AsString(kPhypMixing23) << " = " << input_mixing23 << endl;
  cout << "       > " << util::AsString(kPhypCPV    )  << " = " << input_deltacp << endl;
  cout << "       > " << util::AsString(kPhypDm12sq)   << " = " << input_dm12sq << " eV^2/c^4" << endl;
  cout << "       > " << util::AsString(kPhypDm23sq)   << " = " << input_dm23sq << " eV^2/c^4" << endl;
  if(hyp == kHypOsc3FBarNormBar || hyp == kHypOsc3FNormBar){
    cout << "       > " << util::AsString(kPhypBetaBar)     << " = " << input_beta_bar << endl;
  }
  else{
  }
  
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << fit_param_list << endl;
  cout << "- Outputs" << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;
  cout << *Rnd::Instance() << endl;

  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;

  // Loop over toy-MC experiments

  for(int iexperiment = first_experiment; iexperiment <= last_experiment; iexperiment++) {
    cout << "\n\n *** Pseudo-experiment = " << iexperiment << endl;

    TH2D *hErecoVsQ2Toy;
    TH1D *hErecoToy;
    SuperkObservation * skobs_fakedata = new SuperkObservation;
    
    if(!using_xml_file){
      char name[200];
      sprintf(name,histoname,iexperiment);
      
      cout << name << endl;

      hErecoVsQ2Toy = (TH2D*) frootfile->Get(name);
      //hErecoToy = (TH1D*) frootfile->Get(name)->ProjectionY();
    }
    else{
      // get toy MC with/without systematics
      
      skobs_fakedata->Read(xmlfile2.Data(),iexperiment,using_xml_file);
      hErecoVsQ2Toy = (TH2D*) skobs_fakedata->TwoDSpectrum(s);
      //hErecoToy = (TH1D*) skobs_fakedata->ErecoSpectrum(s);
    }
    
    cout <<  "# of events in " << util::AsStringCompact(s) << " sample: " 
	 << hErecoVsQ2Toy->Integral() << endl;
    
    //hErecoToy->DrawCopy();
    
    //store Nsk
    double nskint = hErecoVsQ2Toy->Integral();
    hNsk->Fill(nskint);
    nsk_per_toy_arr[iexperiment] = nskint;
    
    // Clean-up
    delete skobs_fakedata;  
    delete hErecoVsQ2Toy;

  } // loop over toys

  double dn_over_n = hNsk->GetRMS() / hNsk->GetMean();
  cout << endl;
  cout << "RMS / Mean  = " << dn_over_n << endl;
  cout << endl;

  if (use_new_definition_of_dn_over_n) {
    double sum_nsk = 0;
    for (int iexp = 0 ; iexp < nexperiments; iexp++) {
      sum_nsk += nsk_per_toy_arr[iexp];
    }
    double mean_nsk = sum_nsk / nexperiments;
    
    double sumsq = 0;
    for (int iexp = 0 ; iexp < nexperiments; iexp++) {
      //cerr<<nsk_per_toy_arr[iexp]<<"\t"<<mean_nsk<<endl;
      sumsq += pow(nsk_per_toy_arr[iexp]-mean_nsk,2);
    }
    
    double expected_deviation = sqrt(sumsq/nexperiments);
    cerr<<"Expected deviation: "<<expected_deviation<<endl;
    cerr<<"Nominal total SK: "<<nominal_total_nsk<<endl;
    
    dn_over_n = expected_deviation/nominal_total_nsk;
    
  }

  cout << "dn / n = " << dn_over_n << endl;



  // // Find min/max and calculate fractional deviation from nominal
  // TH1D * hErecoMax = new TH1D(*hErecoNominal);
  // TH1D * hErecoMin = new TH1D(*hErecoNominal);
  // hErecoMax->SetDirectory(0);
  // hErecoMin->SetDirectory(0);
  // hErecoMax->GetYaxis()->SetTitle("Ratio of events at 1 #sigma to nominal");
  // hErecoMin->GetYaxis()->SetTitle("Ratio of events at 1 #sigma to nominal");
  // TH1D * hErecoMaxEnv = new TH1D(*hErecoNominal);
  // TH1D * hErecoMinEnv = new TH1D(*hErecoNominal);
  // hErecoMaxEnv->SetDirectory(0);
  // hErecoMinEnv->SetDirectory(0);
  // for(int ibin=1; ibin <= util::NumErecoBins(s); ibin++) {
  //   double nom   = hErecoNominal->GetBinContent(ibin);
  //   double sumsq = 0;
  //   double sig = TMath::Sqrt(sumsq_per_bin[ibin-1]/nexp);
  //   double min = nom-sig;
  //   double max = nom+sig;
  //   if(nom>0){
  //     hErecoMax -> SetBinContent (ibin,  (max-nom)/nom);
  //     hErecoMin -> SetBinContent (ibin,  (min-nom)/nom);
  //     hErecoMaxEnv -> SetBinContent (ibin,  max);
  //     hErecoMinEnv -> SetBinContent (ibin,  min);
  //   }
  // }//ibin



  return;
}
