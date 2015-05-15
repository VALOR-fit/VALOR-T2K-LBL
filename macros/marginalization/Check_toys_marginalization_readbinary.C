
//
// ------------------------------------------------------------------------------------------------
//
// Check the extraction of the systematic and oscillation parameters done for marginalization
// as well as the Ereco distributions.
//
// Davide Sgalaberna
//
// ------------------------------------------------------------------------------------------------



using namespace std;
using namespace valor;
using namespace valor::conv;

void Check_toys_marginalization_readbinary
(
 const char * config_file              = "analysis_configuration_run1to6_test_app.xml",

 //const char * fit_param_list       = "nuosc3fbar_normbar_fixosc+nuisance",
 //const char * fit_param_list       = "nuosc3fbar_normbar_fitdcpbar+nuisance",
 //const char * fit_param_list       = "nuosc3f_normbar_fitdcp+nuisance",
 //const char * fit_param_list       = "nuosc3f_normbar_fitdcp+nuisance",
 const char * fit_param_list       = "nuosc3f_normbar_fit23_fitdcp+nuisance",

 const char * sk_sample = "RHC1Re",   /* "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */

 const char * hypothesis           = "3f_oscbar_normbar",                       /* options: "2f_numu_disap", "3f_osc, 3f_oscbar, 3f_oscbar_normbar" */
 double       input_mixing12       = -99,                         
 double       input_mixing13       = -99,                         
 double       input_mixing23       = -99,                         
 double       input_deltacp        = -99,                         
 double       input_dm12sq         = -99,                         
 double       input_dm23sq         = -99,                       
 double       input_betabar        = -99,

 const char * dmsqscheme          = "abs23",   /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy           = "normal"   ,
 const char * mixparamchoice      = "sinsqtheta",
 const char * mixparamoctant      = "two",       /* one==kNuOctantOne two==kNuOctantTwo */

 const int InputBin0  = 1,
 const int InputBin1  = 0,

 ////////////////////// MARGINALIZATION PARAMETERS ///////////////////

 const char * infilename           = "../../inputs/marginalization/toys/nuebarapp_1D_T2Kreact_betabar_ToydcpJoint_10kToys/BetaBar_4p0E20POT_toydcp", // (beta=0,1), dcp T2K+react joint

 //const char * infilename           = "../../inputs/marginalization/toys/nuebarapp_1D_T2Kreact_betabar_statonly_fixosc_10kToys/buildmarg_betabar",  // T2K+reactor (betabar=0,1), stat only
 //const char * infilename           = "../../inputs/marginalization/toys/nuebarapp_1D_T2Kreact_betabar_fixosc_10kToys/BetaBar_4p0E20POT",  // T2K+reactor (betabar=0,1)
 //const char * infilename           = "../../inputs/marginalization/toys/nuebarapp_1D_T2Kreact_betabar_fixosc_10kToys_grid/build_marg",  // T2K+reactor (betabar grid)
 //const char * infilename           = "../../inputs/marginalization/toys/nuebarapp_1D_T2KFullStat_T2Kreact_betabar_fixosc_10kToys/BetaBar_7p8E21POT",  // T2K+reactor T2K Full Stat (betabar=0,1) 

 //const int    nfiles               = 201,         // # of files to read
 const int    nfiles               = 1,         // # of files to read

 const int    NtoysMarg            = 10000,   // # toys marginalization 
 
 const char * freeparname0        = "betaBar",  // free parameter
 const double freeparmin0         = 0.,         // minimum value of free parameter
 const double freeparmax0         = 1.,         // maxinum value of free parameter
 const int    NptScan0            = 2,          // # pts of free parameter to scan
 
 const char * freeparname1        = "",  // free parameter
 const double freeparmin1         = -3.14159265358979312e+00, // minimum value of free parameter
 const double freeparmax1         = +3.14159265358979312e+00, // maxinum value of free parameter
 const int    NptScan1            = 1, //201, //21,      // # pts of free parameter to scan

 const bool  doMargMH             = true,    // marginalize Mass Hierarchy
 
 )
{
  
  FileStat_t dummy;
  
  TString opt(fit_param_list); 
  
  // Analysis configuration
  TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  
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

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());
 
  // Fit param list for fakedata sets
  FitParamList * fparm  = util::NamedParameterList(fit_param_list); 
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used for fakedata sets: " << util::AsString(hyp) << endl;

  // Set Marginalization

  PhysParamId_t freepar0;
  string sfreeparname0(freeparname0);  
  // NuBar mode
  if(sfreeparname0=="mixing23Bar")      freepar0 = conv::kPhypMixing23Bar;
  else if(sfreeparname0=="mixing13Bar") freepar0 = conv::kPhypMixing13Bar;
  else if(sfreeparname0=="dmsq23Bar")   freepar0 = conv::kPhypDm23sqBar;
  else if(sfreeparname0=="dmsq12Bar")   freepar0 = conv::kPhypDm12sqBar;
  else if(sfreeparname0=="mixing12Bar") freepar0 = conv::kPhypMixing12Bar;
  else if(sfreeparname0=="deltacpBar")  freepar0 = conv::kPhypCPVBar;
  else if(sfreeparname0=="betaBar")     freepar0 = conv::kPhypBetaBar;
  // Nu mode
  else if(sfreeparname0=="mixing23")    freepar0 = conv::kPhypMixing23;
  else if(sfreeparname0=="mixing13")    freepar0 = conv::kPhypMixing13;
  else if(sfreeparname0=="dmsq23")      freepar0 = conv::kPhypDm23sq;
  else if(sfreeparname0=="dmsq12")      freepar0 = conv::kPhypDm12sq;
  else if(sfreeparname0=="mixing12")    freepar0 = conv::kPhypMixing12;
  else if(sfreeparname0=="deltacp")     freepar0 = conv::kPhypCPV;
  //
  else{
    freepar0 = conv::kPhypUndefined;
    cout << "Param 0 is not defined" << endl;
    exit(1);
  }
  
  PhysParamId_t freepar1 = conv::kPhypUndefined;
  string sfreeparname1(freeparname1);  
  // NuBar mode
  if(sfreeparname1=="mixing23Bar")      freepar1 = conv::kPhypMixing23Bar;
  else if(sfreeparname1=="mixing13Bar") freepar1 = conv::kPhypMixing13Bar;
  else if(sfreeparname1=="dmsq23Bar")   freepar1 = conv::kPhypDm23sqBar;
  else if(sfreeparname1=="dmsq12Bar")   freepar1 = conv::kPhypDm12sqBar;
  else if(sfreeparname1=="mixing12Bar") freepar1 = conv::kPhypMixing12Bar;
  else if(sfreeparname1=="deltacpBar")  freepar1 = conv::kPhypCPVBar;
  else if(sfreeparname1=="betaBar")     freepar1 = conv::kPhypBetaBar;
  // Nu mode
  else if(sfreeparname1=="mixing23")    freepar1 = conv::kPhypMixing23;
  else if(sfreeparname1=="mixing13")    freepar1 = conv::kPhypMixing13;
  else if(sfreeparname1=="dmsq23")      freepar1 = conv::kPhypDm23sq;
  else if(sfreeparname1=="dmsq12")      freepar1 = conv::kPhypDm12sq;
  else if(sfreeparname1=="mixing12")    freepar1 = conv::kPhypMixing12;
  else if(sfreeparname1=="deltacp")     freepar1 = conv::kPhypCPV;
  //
  else{
    freepar1 = conv::kPhypUndefined;
    cout << "1D fit" << endl;
  }
  
  
  // Set free param, ntoys, # pt to scan, min par, max par (if min and max not passed, use default ones for that par)
  
  // free parameter 
  PhysParamId_t parfree[2];
  parfree[0] = freepar0;
  parfree[1] = freepar1;
  
  // # of toys for marginalization
  int ntoys = NtoysMarg;
  
  // # of free parameter pts to scan 
  int nbinspar[2];
  nbinspar[0] = NptScan0;
  nbinspar[1] = NptScan1;
  
  // free parameter range
  double minpar[2];
  double maxpar[2];
  minpar[0] = freeparmin0;
  maxpar[0] = freeparmax0;
  minpar[1] = freeparmin1;
  maxpar[1] = freeparmax1;
  
  // free parameter sub-range
  
  const int    first_bin_freepar0   = -999; // first freepar of the sub-range
  const int    last_bin_freepar0    = -999; // last freepar of the sub-range
  const int    first_bin_freepar1   = -999; // first freepar of the sub-range
  const int    last_bin_freepar1    = -999; // last freepar of the sub-range
  
  int subminbin[2];
  int submaxbin[2];
  subminbin[0] = first_bin_freepar0; submaxbin[0] = last_bin_freepar0;
  subminbin[1] = first_bin_freepar1; submaxbin[1] = last_bin_freepar1;
  
  // Implement the constructor giving the following input parameters:
  FitMarginalization *fitmarg = new FitMarginalization();   
  fitmarg->SetPhysParFree(parfree,ntoys,nbinspar,minpar,maxpar);
  fitmarg->SetSubRange(subminbin,submaxbin); // By default it scans all the range
  fitmarg->Init(skexp); // Before need SetPhysParFree() to initialize the arrays
  
  // Set Mass Hierarchy marginalization
  fitmarg->SetMargMassHierarchy(doMargMH); // false by default 
    
  // Set fit parameter list
  fitmarg->SetFitParamList(fit_param_list);
  fitmarg->PrintPhysHypothesis();
  
  
  // Build the map fileID - sub-range
  
  fitmarg->SetMaxIterationSingleFit(2);  // Max # of iterations per fit
  fitmarg->SetMaxFitNumber(3); // Max # of fits max # of iterations
  
  fitmarg->SetTagFileName(infilename);
  fitmarg->SetFileID(nfiles);
  
  
  // Get file ID corresponding to the input free parameters
  
  int fileid = fitmarg->GetFileIDRead(0,0);
  fitmarg->ReadToysFromBinary(fileid,infilename);



  // Get Ereco toys
  
  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }

  //fitmarg->GetToyMargHisto(s,0,0,0,0)->Draw();
  //cout << fitmarg->GetToyMargHisto(s,0,0,0,1)->Integral() << endl;


  
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
  
  int ipar0 = InputBin0;
  int ipar1 = InputBin1;
  
  cout << "par0 = " << ipar0 << endl;
  cout << "par1 = " << ipar1 << endl;

  TH1F* h_marg_toy;
  TH1F* h_up = new TH1F("h_up","h_up",util::NumErecoBins(s),util::ErecoBinEdges(s));
  TH1F* h_down = new TH1F("h_down","h_down",util::NumErecoBins(s),util::ErecoBinEdges(s));
  TH1D* h_nominal;
  TH2F* h_2d = new TH2F("h_2d","Envelope",util::NumErecoBins(s),util::ErecoBinEdges(s),1000,0,6);
  
  for(int itoy=0;itoy<NtoysMarg;itoy++){
    
    //Normal Hierarchy
    h_marg_toy = fitmarg->GetToyMargHisto(s,itoy,ipar0,ipar1,0);
    for(int i=1;i<=util::NumErecoBins(s);i++){
      h_2d->Fill(h_marg_toy->GetBinCenter(i),h_marg_toy->GetBinContent(i));
    }
    
    //cout << "NH: " << h_marg_toy->Integral() << endl;
      
    //Inverted Hierarchy
    h_marg_toy = fitmarg->GetToyMargHisto(s,itoy,ipar0,ipar1,1);
    for(int i=1;i<=util::NumErecoBins(s);i++){
      h_2d->Fill(h_marg_toy->GetBinCenter(i),h_marg_toy->GetBinContent(i));
    }
    
    //cout << "IH: " << h_marg_toy->Integral() << endl;
    //cout << endl;
  }
  
  // Get tuned MC - with BANFF tuning

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
  if(input_betabar < -1) input_betabar = conv::kGlobFitBetaBar;
 
  PhysHypothesis ph(hypothesis);
  if(ph.Hypothesis() == conv::kHypNuMuDisapp2F) {
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3F) {
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3FBar){
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    ph.SetCurrValue(conv::kPhypMixing12Bar,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13Bar,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23Bar,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPVBar    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sqBar,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sqBar,    input_dm23sq);
  } else
  if(ph.Hypothesis() == conv::kHypOsc3FBarNormBar){
    ph.SetCurrValue(conv::kPhypMixing12,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPV    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sq,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sq,    input_dm23sq);
    ph.SetCurrValue(conv::kPhypMixing12Bar,  input_mixing12);
    ph.SetCurrValue(conv::kPhypMixing13Bar,  input_mixing13);
    ph.SetCurrValue(conv::kPhypMixing23Bar,  input_mixing23);
    ph.SetCurrValue(conv::kPhypCPVBar    ,   input_deltacp);
    ph.SetCurrValue(conv::kPhypDm12sqBar,    input_dm12sq);
    ph.SetCurrValue(conv::kPhypDm23sqBar,    input_dm23sq);
    ph.SetCurrValue(conv::kPhypBetaBar,    input_betabar);
  }
 
  skexp->Generate(ph);
  h_nominal = new TH1D(*skexp->CurrPredictedErecoSpectrum(s));
  //cout << h_nominal->Integral() << endl;
 



  for(int i=1;i<=util::NumErecoBins(s);i++){
    TH1D* h_proj = h_2d->ProjectionY("",i,i);
    double average = average(h_proj);
    double rms = rms(h_proj,average);
    h_up->SetBinContent(i,average+ rms);
    h_down->SetBinContent(i,average-rms);  
  }
  
  TCanvas * c = new TCanvas("c","",20,20,1000,1000);
  h_2d->Draw("colz");
  h_up->Draw("same");
  h_down->Draw("same");
  h_nominal->SetLineColor(kYellow);
  h_nominal->Draw("same");
  
  //c->SaveAs("draw_marg_toys.root"); 
  // // Store toy distributions in output file  
  // char outname[200];
  // sprintf(outname,"check_toy_throws.root",outfilename);
  // //~ cout << outname << endl;
  // TFile *fout = new TFile(outname,"RECREATE");
  // fout->cd();
  // TObjArray * arr = fitmarg->GetToyParamThrowHistograms();
  // arr->Write();
  // fout->Close();
    
  return;
}

double average(TH1D* h1){
  double avg(0);
  
  for(int i=1; i<=h1->GetNbinsX();i++){
    avg += h1->GetBinCenter(i)*h1->GetBinContent(i);
  }
  avg = avg/h1->GetEntries();
  
  return avg;
}

double rms(TH1D* h1,double avg){
  double rms(0);
  
  for(int i=1; i<=h1->GetNbinsX();i++){
    rms += ((h1->GetBinCenter(i) - avg) **2) * (h1->GetBinContent(i));
  }
  rms = rms/h1->GetEntries();
  rms = rms**0.5;
  return rms;
}
