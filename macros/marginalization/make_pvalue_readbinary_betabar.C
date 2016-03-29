
// ------------------------------------------------------------------------------------------------
//
This macro takes the marg toys produced by buid_marg_toys and
Raj Shah
//
// ------------------------------------------------------------------------------------------------




//
// TODO:
// - set correct oscillation parameters for reactor constraint
//

using namespace std;
using namespace valor;
using namespace valor::conv;

void make_pvalue_readbinary_betabar
(
 
 ////////////////////////// GENERAL PARAMETERS /////////////////////
 
 const char * tag                      = "fitmarg_betabar_DEV",    /* tag for output file */
 const char * xmlfile                  = "../outputs/test/toymc_prova_joint-0.xml",
 const char * rootfile                 = "../outputs/test/toymc_prova_joint-0.root",
 bool using_xml_file                   = false,
  bool         dotoys                   = 1,
 // Default
 double       input_mixing12           = -99,
 double       input_mixing13           = -99,
 double       input_mixing23           = -99,
 double       input_deltacp            = -99,
 double       input_dm12sq             = -99,
 double       input_dm23sq             = -99,
 double       input_beta_bar           = 1., //0.2,
 
 const char * config_file              = "analysis_configuration_run1to6_test_app_4.0108E20POT_fluxJune.xml",
 
 bool         stat_fluctuations        = false,    /* include statistical fluctuation in toy MC expts? (not for marginalization) */
 int          tweak_systematics        = 0,       /* generate toy (both fakedata sets and marginalization) MC with randomized systematics? */

 
 //const char * fit_param_list       = "nuosc3fbar_normbar_fixall+nuisance",
 //const char * fit_param_list       = "nuosc3fbar_normbar_fixosc+nuisance",
 //const char * fit_param_list       = "nuosc3fbar_normbar_fitdcpbar+nuisance",
 //const char * fit_param_list       = "nuosc3f_normbar_fitdcp+nuisance",
 //const char * fit_param_list       = "nuosc3f_normbar_fit23_fitdcp+nuisance",
 const char * fit_param_list           = "nuosc3f_normbar_fix12+nuisance",
 
 const char * dmsqscheme               = "auto",   /* auto==kNuDmsqAuto abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * mixparamchoice           = "sinsqtheta", /* sinsq2theta or sinsqtheta */
 const char * hierarchy                = "normal",     /* normal or inverted */
 const char * mixparamoctant           = "two",       /* one==kNuOctantOne two==kNuOctantTwo */
 
 ////////////////////// MARGINALIZATION PARAMETERS ///////////////////
 
 //const char * infilename = "../../inputs/marginalization/toys/EPS_RELEASE/BetaBar_Toy_23_th13dcp_Corr_Mach3_fluxJune_50k_DEV",
 const char * infilename = "../../inputs/marginalization/toys/EPS_RELEASE/BetaBar_Toy_23_th13dcp_Corr_Mach3_fluxJune_10k_DEV",
 
 const int    nfiles               = 2, //401, 1        // # of files to read
 
 const int    NtoysMarg            = 10000,  // # toys marginalization
 
 const bool  doMargMH             = true,   // marginalize Mass Hierarchy
)
{

  //Set the first parameter BetaBar
  const char * freeparname0 = "betaBar";  // free parameter
  PhysParamId_t freepar0    = conv::kPhypBetaBar;
  const double freeparmin0  = 0;
  const double freeparmax0  = 1;
  const int NptScan0        = 0;
  
  // Set the second parameter as undefined
  const char * freeparname1        = "";  // free parameter
  PhysParamId_t freepar1      = conv::kPhypUndefined;
  const double freeparmin1    = 0.; // minimum value of free parameter
  const double freeparmax1    = 0.; // maxinum value of free parameter
  const int    NptScan1       = 1; //201, //21,      // # pts of free parameter to scan
                                        //
  
  TString opt(fit_param_list);
  
  // Analysis configuration
  TString config  (Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file));
  TString xmlfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile));
  TString rootfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), rootfile));

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
  
  if (!dotoys && gSystem->GetPathInfo(xmlfile2.Data(), dummy))
  {
    xmlfile2 = xmlfile;
    if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) {
      cerr<<"Data input file not found: "<<xmlfile<<endl;
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
  
  // Check whether some osc. params need to be randomized
  bool randomize_mixing12 = (TMath::Abs(input_mixing12 + 999) < kASmallNum);
  bool randomize_mixing13 = (TMath::Abs(input_mixing13 + 999) < kASmallNum);
  bool randomize_mixing23 = (TMath::Abs(input_mixing23 + 999) < kASmallNum);
  bool randomize_deltacp  = (TMath::Abs(input_deltacp  + 999) < kASmallNum);
  bool randomize_dm12sq   = (TMath::Abs(input_dm12sq   + 999) < kASmallNum);
  bool randomize_dm23sq   = (TMath::Abs(input_dm23sq   + 999) < kASmallNum);
  
  
  bool randomize_beta_bar   = (TMath::Abs(input_beta_bar   + 999) < kASmallNum);
  
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
  
  
  
  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << tag << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  
  
  
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << fit_param_list << endl;
  cout << "- Outputs" << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;
  cout << *Rnd::Instance() << endl;
  
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;
  
  // Object for storing input oscillation and nuisance params in each toy-MC
  PhysParamList * iparm = new PhysParamList;
  
  SuperkObservation * skobs_fakedata = new SuperkObservation;

  if(!dotoys){
    if(using_xml_file){
      cout << "Use xml file" << endl;
      skobs_fakedata->Read(xmlfile2.Data(),iexperiment,using_xml_file);
    }
    else{
      cout << "Use root file" << endl;
      skobs_fakedata->Read(rootfile2.Data(),iexperiment,using_xml_file);
    }
  }
  else{
    bool applied_negevents_condition = true;
    bool unphysical_tweak = true;
    bool rethrow_if_parameters_unphysical = true;// both marginalization and fakedata sets
    while(applied_negevents_condition || unphysical_tweak) {
      if(tweak_systematics) {
        // tweak all sources of systematic error
        if(tweak_systematics == 1) {
          cout << "- Randomizing all systematic params in pdf..." << endl;
          skexp->TweakAllSystematics();
        }
        else if(tweak_systematics == 2) {
          cout << "- Randomizing BANFF systematic params in pdf..." << endl;
          skexp->TweakAllBANFFSystematics();
        }
        else if(tweak_systematics == 3) {
          cout << "- Randomizing uncorrelated-xsec systematic params in pdf..." << endl;
          skexp->TweakAllUncorrelatedXSecSystematics();
        }
        else if(tweak_systematics == 4) {
          cout << "- Randomizing SKDetFSI systematic params in pdf..." << endl;
          skexp->TweakAllSKDetectorFSISystematics();
        }
        skexp->PrintSystematicVariations();
        
        PhysParamList true_nuisance_params = skexp->Tweaks();
        if (rethrow_if_parameters_unphysical) {
          unphysical_tweak = !(valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(true_nuisance_params));
          if (unphysical_tweak) {cout<<"Parameters are unphysical - rethrowing"<<std::endl;}
        }
        else {
          unphysical_tweak = false;
        }
        
        // generate mock data using the current PDF (incl. systematic variations)
        // + statistical fluctuations if requested (Poisson fluctuation in each bin)
        cout << "- Generating toy-MC..." << endl;
        //note that fparm does not contain the tweaks, so use only the PhysHyp here
        applied_negevents_condition = skobs_fakedata->CreateMockData(*skexp, PhysHypothesis(*fparm), stat_fluctuations);
        cout << "applied_negevents_condition = " << applied_negevents_condition << endl;
        
        
        if(!(applied_negevents_condition || unphysical_tweak)) {
          // keep track of all nuisance params set so as to add in output tree
          iparm->Append(true_nuisance_params);
        }
        // reset all systematics so that the PDF used in the fitter is the nominal one
        cout << "- Reseting all systematic params in pdf..." << endl;
        skexp->ResetSystematicVariations();
        //skexp->PrintSystematicVariations();
      }
      else {
        //make sure systematics are reset
        skexp->ResetSystematicVariations();
        // generate mock data using the current PDF (without systematic variations)
        // with statistical fluctuations if requested (Poisson fluctuation in each bin)
        cout << "- Generating toy-MC..." << endl;
        applied_negevents_condition = skobs_fakedata->CreateMockData(*skexp, PhysHypothesis(*fparm), stat_fluctuations);
        cout << "applied_negevents_condition = " << applied_negevents_condition << endl;
        
        unphysical_tweak = false; // Otherwise it doesn't exit the loop
      }
    }//applied_negevents_condition
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
  fitmarg->SetTestStatistic("rate+shape");     // Set Test Statistic ("rate-only" or "rate+shape")
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
  
  
  //betabar bin 0
  fitmarg->FillExpectedNuEBarAppHisto(skobs_fakedata, 0, 1); // Fill expected NueBar app histo, betabin=0
  
  //betabar bin 1
  fitmarg->FillExpectedNuEBarAppHisto(skobs_fakedata, 1, 1); // Fill expected NueBar app histo, betabin=1
  
  
  // Create file to store results
  TFile f(Form("%s.root",tag), "recreate");

  
  fitmarg->GetNuEBarHisto(0)->Write("Beta0");
  fitmarg->GetNuEBarHisto(1)->Write("Beta1");
  
  f.Close();

  // Clean-up
  delete skobs_fakedata;
  
  fitmarg->CleanUp();
  delete fitmarg;
  
  
  
  
  
  
  
  // Store ntuple with best-fit results
  
  
  return;
}
