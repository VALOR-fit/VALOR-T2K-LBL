
//
// ------------------------------------------------------------------------------------------------           
//
// Run toys for marginalization and store into binary files. 
// Toys for different sub-ranges of the free paramters can be run in parallel jobs.
// A file .bin contains the array, while a twin file .bin contains the sub-range of the free parameter
// Each sub-range contains the same number of toys, coming from the same random seed.
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------

using namespace std;
using namespace valor;
using namespace valor::conv;

void build_marg_toys
(
 const int    fileid               = 2,   // file id

 const int    first_bin_freepar0   = 18,   // first freepar of the sub-range
 const int    last_bin_freepar0    = 20, // last freepar of the sub-range
 const int    first_bin_freepar1   =-999,   // first freepar of the sub-range
 const int    last_bin_freepar1    =-999, // last freepar of the sub-range

 const char * tag                  = "MH_FitMarg_mixing23_1p7E20_prefit",
 //const char * tag                  = "submit_mixing23_1p7E20_postfit_NoRethrow",

 //const char * config_file          = "analysis_configuration_1.7E20_prefit.xml",   
 const char * config_file          = "analysis_configuration_1.7E20.xml",   

 //const char * fit_param_list     = "nuosc3fbar_numudisapp_constrainall+nuisance",  
 const char * fit_param_list       = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",

 const char * freeparname0          = "mixing23Bar",  // free parameter
 const double freeparmin0           = 0,     // minimum value of free parameter
 const double freeparmax0           = 1,     // maxinum value of free parameter
  const int   NptScan0              = 21,  // # pts of free parameter to scan

 const char * freeparname1          = "dmsq23Bar",  // free parameter
 //const char * freeparname1          = "",  // free parameter
 const double freeparmin1           = 0.0002,     // minimum value of free parameter
 const double freeparmax1           = 0.02,     // maxinum value of free parameter
 const int    NptScan1              = 21,  // # pts of free parameter to scan

 const int    NtoysMarg            = 10,   // # toys marginalization 

 const int tweak_systematics       = 1,      /* generate toy MC with randomized systematics? */

 const bool   doMargMH             = true,    // marginalize Mass Hierarchy
 const bool   doRethrowUnphysical  = false,    // rethrow unphysical systematic toys

 const char * dmsqscheme              = "abs23",       /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy               = "normal",    /* normal or inverted */
 const char * mixparamchoice          = "sinsqtheta",  /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant          = "one"         /* one==kNuOctantOne two==kNuOctantTwo */
 )
{
  FileStat_t dummy;
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
 
  PhysParamId_t freepar0;
  string sfreeparname0(freeparname0);  
  // NuBar mode
  if(sfreeparname0=="mixing23Bar")      freepar0 = conv::kPhypMixing23Bar;
  else if(sfreeparname0=="mixing13Bar") freepar0 = conv::kPhypMixing13Bar;
  else if(sfreeparname0=="dmsq23Bar")   freepar0 = conv::kPhypDm23sqBar;
  else if(sfreeparname0=="dmsq12Bar")   freepar0 = conv::kPhypDm12sqBar;
  else if(sfreeparname0=="mixing12Bar") freepar0 = conv::kPhypMixing12Bar;
  else if(sfreeparname0=="deltacpBar")  freepar0 = conv::kPhypCPVBar;
  // Nu mode
  else if(sfreeparname0=="mixing23")    freepar0 = conv::kPhypMixing23;
  else if(sfreeparname0=="mixing13")    freepar0 = conv::kPhypMixing13;
  else if(sfreeparname0=="dmsq23")      freepar0 = conv::kPhypDm23sq;
  else if(sfreeparname0=="dmsq12")      freepar0 = conv::kPhypDm12sq;
  else if(sfreeparname0=="mixing12")    freepar0 = conv::kPhypMixing12;
  else if(sfreeparname0=="deltacp")     freepar0 = conv::kPhypCPV;
  //
  else{
    // Parameter 0 must be always defined
    cout << "Parameter 0 is not defined" << endl;
    exit(1);
  }

  PhysParamId_t freepar1;
  string sfreeparname1(freeparname1);  
  // NuBar mode
  if(sfreeparname1=="mixing23Bar")      freepar1 = conv::kPhypMixing23Bar;
  else if(sfreeparname1=="mixing13Bar") freepar1 = conv::kPhypMixing13Bar;
  else if(sfreeparname1=="dmsq23Bar")   freepar1 = conv::kPhypDm23sqBar;
  else if(sfreeparname1=="dmsq12Bar")   freepar1 = conv::kPhypDm12sqBar;
  else if(sfreeparname1=="mixing12Bar") freepar1 = conv::kPhypMixing12Bar;
  else if(sfreeparname1=="deltacpBar")  freepar1 = conv::kPhypCPVBar;
  // Nu mode
  else if(sfreeparname1=="mixing23")    freepar1 = conv::kPhypMixing23;
  else if(sfreeparname1=="mixing13")    freepar1 = conv::kPhypMixing13;
  else if(sfreeparname1=="dmsq23")      freepar1 = conv::kPhypDm23sq;
  else if(sfreeparname1=="dmsq12")      freepar1 = conv::kPhypDm12sq;
  else if(sfreeparname1=="mixing12")    freepar1 = conv::kPhypMixing12;
  else if(sfreeparname1=="deltacp")     freepar1 = conv::kPhypCPV;
  //
  else{
    cout << "Parameter 1 is not defined" << endl;
    cout << "1D fit" << endl;
    freepar1 = conv::kPhypUndefined;
  }



  // Implement the constructor giving the following input parameters:
  // free param, ntoys, # pt to scan, min par, max par.
  // If min and max not passed, use default ones for that parameter (see FitMarginalization.h)
  
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
  int subminbin[2];
  int submaxbin[2];
  subminbin[0] = first_bin_freepar0; submaxbin[0] = last_bin_freepar0;
  subminbin[1] = first_bin_freepar1; submaxbin[1] = last_bin_freepar1;

  // Implement the constructor giving the following input parameters:
  FitMarginalization *fitmarg = new FitMarginalization();   
  fitmarg->SetPhysParFree(parfree,ntoys,nbinspar,minpar,maxpar);
  fitmarg->SetSubRange(subminbin,submaxbin); // By default it scans all the range
  fitmarg->Init(); // Before need SetPhysParFree() to initialize the arrays


  // Set fit parameter list
  fitmarg->SetFitParamList(fit_param_list);
  fitmarg->PrintPhysHypothesis();      
  
  // Set Mass Hierarchy marginalization
  fitmarg->SetMargMassHierarchy(doMargMH); // false by default 
  
  // Decide whether to rethrow or not the unphysical systematic toys 
  fitmarg->SetRethrowSystUnphysical(doRethrowUnphysical);

  // Set priors for not fixed oscillation parameters
  fitmarg->SetPrior(kPhypMixing23Bar,"flat",0,1);
  fitmarg->SetPrior(kPhypMixing13Bar,"flat",0,1);
  fitmarg->SetPrior(kPhypMixing12Bar,"flat",0,1);
  fitmarg->SetPrior(kPhypDm23sqBar,"flat",conv::kMinFitDm23sq,conv::kMaxFitDm23sq);
  fitmarg->SetPrior(kPhypDm12sqBar,"flat",conv::kMinFitDm12sq,conv::kMaxFitDm12sq);
  fitmarg->SetPrior(kPhypCPVBar,"flat",-TMath::Pi(),+TMath::Pi());
  

  // Make toys and store in binary files
    
  // Produce toys for different sub-ranges of the free parameter and store in different files.
  // Useful for running on cluster when many toys are needed (e.g. 10k x 200)

  // int min_bin[2];
  // int max_bin[2];    
  // min_bin[0] = first_bin_freepar0;
  // max_bin[0] = last_bin_freepar0;
  // min_bin[1] = first_bin_freepar1; 
  // max_bin[1] = last_bin_freepar1;
  //fitmarg->MakeToysMarginalization(tweak_systematics,skexp,min_bin,max_bin); 
  fitmarg->MakeToysMarginalization(tweak_systematics,skexp); 
  fitmarg->WriteToys2Binary(fileid,tag);

  return;
}
