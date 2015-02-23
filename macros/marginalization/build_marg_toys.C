
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
 const int    fileid               = 0,   // file id
 const int    first_bin_freepar    =-999,   // first freepar of the sub-range
 const int    last_bin_freepar     =-999, // last freepar of the sub-range

 //const char * tag                  = "DEV_FitMarg_mixing23_1p7E20_prefit",
 const char * tag                  = "MH_FitMarg_mixing23_1p7E20_prefit",

 //const char * config_file          = "analysis_configuration_1.7E20_prefit.xml",   
 const char * config_file          = "analysis_configuration_1.7E20.xml",   

 //const char * fit_param_list     = "nuosc3fbar_numudisapp_constrainall+nuisance",  
 const char * fit_param_list       = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",

 const char * freeparname          = "mixing23Bar",  // free parameter
 const int    freeparmin           = 0,     // minimum value of free parameter
 const int    freeparmax           = 1,     // maxinum value of free parameter

 const int    NtoysMarg            = 10,   // # toys marginalization 
 const int    NptScan              = 201,  // # pts of free parameter to scan

 const int tweak_systematics       = 1,      /* generate toy MC with randomized systematics? */

 const bool   doMargMH             = true,    // marginalize Mass Hierarchy
 const bool   doRethrowUnphysical  = true,    // rethrow unphysical systematic toys

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
 
  PhysParamId_t freepar;
  string sfreeparname(freeparname);
  
  // NuBar mode
  if(sfreeparname=="mixing23Bar")      freepar = conv::kPhypMixing23Bar;
  else if(sfreeparname=="mixing13Bar") freepar = conv::kPhypMixing13Bar;
  else if(sfreeparname=="dmsq23Bar")   freepar = conv::kPhypDm23sqBar;
  else if(sfreeparname=="dmsq12Bar")   freepar = conv::kPhypDm12sqBar;
  else if(sfreeparname=="mixing12Bar") freepar = conv::kPhypMixing12Bar;
  else if(sfreeparname=="deltacpBar")  freepar = conv::kPhypCPVBar;
  // Nu mode
  else if(sfreeparname=="mixing23")    freepar = conv::kPhypMixing23;
  else if(sfreeparname=="mixing13")    freepar = conv::kPhypMixing13;
  else if(sfreeparname=="dmsq23")      freepar = conv::kPhypDm23sq;
  else if(sfreeparname=="dmsq12")      freepar = conv::kPhypDm12sq;
  else if(sfreeparname=="mixing12")    freepar = conv::kPhypMixing12;
  else if(sfreeparname=="deltacp")     freepar = conv::kPhypCPV;
  //
  else{
    cout << "Wrong PhysParamId" << endl;
    exit(1);
  }

  // Implement the constructor giving the following input parameters:
  // free param, ntoys, # pt to scan, min par, max par.
  // If min and max not passed, use default ones for that parameter (see FitMarginalization.h)
  FitMarginalization *fitmarg = new FitMarginalization(freepar,NtoysMarg,NptScan,freeparmin,freeparmax); 

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
  
  fitmarg->MakeToysMarginalization(tweak_systematics,skexp,first_bin_freepar,last_bin_freepar); 
  fitmarg->WriteToys2Binary(fileid,tag);

  return;
}
