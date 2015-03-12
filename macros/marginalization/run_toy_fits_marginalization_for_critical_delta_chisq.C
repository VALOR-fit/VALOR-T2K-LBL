
// ------------------------------------------------------------------------------------------------           
//
// Intended to calculate a critical delta Chi-square surface for the marginalized toys.
// Currently wildly unfinished.
// 
// Steve Dennis, March 2015
//                                                                                           
// ------------------------------------------------------------------------------------------------


using namespace std;
using namespace valor;
using namespace valor::conv;

void run_toy_fits_marginalization_for_critical_delta_chisq
(
  const bool   readtoys             = true,
  const char * outfile_tag          = "marginalization_toys",
  
  const int    Nfakedata           = 10000,
  
  const int    IDoutfile            = 0,
  const int    first_bin_freepar0   = -999,   // first freepar of the sub-range
  const int    last_bin_freepar0    = -999, // last freepar of the sub-range
  const int    first_bin_freepar1   = -999,   // first freepar of the sub-range
  const int    last_bin_freepar1    = -999, // last freepar of the sub-range
  
  const char * infilename           = "toys",
  
  const int    nfiles               = 1,         // # of files to read
  
  const char * config_file          = "analysis_configuration_2.5E20.xml",  
  
  const char * fit_param_list       = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",
  
  const int    NtoysMarg            = 10000,   // # toys marginalization 
  
  const char * freeparname0        = "mixing23Bar",  // free parameter
  const double freeparmin0         = 0.,             // minimum value of free parameter
  const double freeparmax0         = 1.,             // maxinum value of free parameter
  const int    NptScan0            = 21,             // # pts of free parameter to scan
  
  const char * freeparname1        = "dmsq23Bar",  // free parameter
  const double freeparmin1         = 0.0002, // minimum value of free parameter
  const double freeparmax1         = 0.02, // maxinum value of free parameter
  const int    NptScan1            = 21,      // # pts of free parameter to scan
  
  int         tweak_systematics    = 1,       /* generate toy MC with randomized systematics? */
  
  const bool  doMargMH             = true,    // marginalize Mass Hierarchy
  const bool  doRethrowUnphysical  = true,    // rethrow unphysical systematic toys
  
  const char * dmsqscheme          = "abs23",      /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
  const char * hierarchy           = "normal",     /* normal or inverted */
  const char * mixparamchoice      = "sinsqtheta", /* sinsq2theta or sinsqtheta */
  const char * mixparamoctant      = "one"         /* one==kNuOctantOne two==kNuOctantTwo */
)
{
  
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
  
  FileStat_t dummy;

  TString config    = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

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
  
  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());
  
  
  PhysParamId_t freepar0 = util::PhysParamIdFromString(freeparname0);
  if (freepar0 < 0 || util::IsSystematic(freepar0))
  {
    cout << "Cannot handle parameter 0 "<< freeparname0 <<", p= "<<freepar0 << "\t"<<util::AsString(freepar0) << endl;
    exit(1);
  }

  PhysParamId_t freepar1;
  string sfreeparname1(freeparname1);  
  if ( sfreeparname1=="" || sfreeparname1=="none" )
  {
    cout << "1D fit" << endl;
    freepar1 = conv::kPhypUndefined;
  }
  else {
    freepar1 = util::PhysParamIdFromString(freeparname1);
    
    if (freepar1 < 1 || util::IsSystematic(freepar1))
    {
      cout << "Cannot handle parameter 1 "<< freeparname1 <<", p= "<<freepar1 << "\t"<<util::AsString(freepar1) << endl;
      exit(1);
    }
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
  
  const FitParamList fparm = fitmarg->FitParamList();
  
  // Start generating fake datasets
  SuperkObservation * skobs = new SuperkObservation;
  
  VALORArrayF   * fFakeDataEreco[conv::kns][NFakedata];
  PhysParamList * fTrueValues   [NFakedata];
  double        * fTrueChisq    [NFakedata];
  
  PhysHypothesis ph(*fparm);
  
  for(int ifakedata = 0 ; ifakedata < Nfakedata; ifakedata++) {
    bool got_good_tweaks = false;
    while (true) {
      skexp->RandomizeSystematics();
      PhysParamList true_nuisance_params = skexp->Tweaks();
      got_good_tweaks = (not doRethrowUnphysical) || T2KPhysicalParameterLimits::Instance()->AllWithinLimits(true_nuisance_params);
      if (got_good_tweaks) {
        fTrueValues[ifakedata] = new PhysParamList(true_nuisance_params);
        break;
      }
    }
    
    for (ip = 0 ; ip < conv::kNuis000; i++) {// loop over physics parameters
      PhysParamId_t p = (PhysParamId_t)ip;
      if (not fparm->Exists(p)) {continue}
      if (fparm->IsFixed(p) and (p!=freepar0) (p!=freepar1) ) {continue}
      double value = fparm->Randomize(p);
      fTrueValues[ifakedata]->SetValue(p,value);
      ph.SetValue(p,value);
    }
    
    skobs->Generate(skexp,ph);
    
    for (int is = 0; is < conv::kns ; is++) {
        SuperkSample_t s = conv::kSKsamples[is];
        if (s==conv::kSKsUndefined) {continue;}
        fFakeDataEreco[is][ifakedata] = new VALORArrayF(skobs->ErecoSpectrumVArray(s));
    }
    
    fTrueChisq[NFakeData] = skobs->Statistic(skexp,fparm);
  }
  
  // 
  
  // Read toys from binary files and fill marginalized chisq
  for(int ifile=0;ifile<nfiles;ifile++){
  
    cout << endl;
    cout << "////////////////////////////////" << endl;
    cout << endl;
    cout << "Read file ID " << ifile << endl; 
    cout << endl;

    fitmarg->ReadToysFromBinary(ifile,infilename);

    fitmarg->FillChisqMarg(skobs);
    fitmarg->Minimize();

  }  
  
  
  
  return;
}
