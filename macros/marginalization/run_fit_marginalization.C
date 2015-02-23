
using namespace std;
using namespace valor;
using namespace valor::conv;

void run_fit_marginalization
(
 const bool   readtoys             = true,

 //const char * infilename           = "../../inputs/marginalization/toys/margtoys_cluster",
 //const char * infilename           = "DEV_FitMarg_mixing23_1p7E20_prefit",
 //const char * infilename           = "MH_FitMarg_mixing23_1p7E20_prefit",
 const char * infilename           = "ciao",

 const int    nfiles               = 21,         // # of files to read

 const char * config_file          = "analysis_configuration_1.7E20_prefit.xml",   /* analysis configuration file */

 //const char * fit_param_list      = "nuosc3fbar_numudisapp_constrainall+nuisance",  
 const char * fit_param_list       = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",

 const char * freeparname          = "mixing23Bar",  // free parameter
 const int    freeparmin           = 0,     // minimum value of free parameter
 const int    freeparmax           = 1,     // maxinum value of free parameter

 const int    NtoysMarg            = 10,   // # toys marginalization 
 const int    NptScan              = 21,   // # pts of free parameter to scan

 int          tweak_systematics    = 1,       /* generate toy MC with randomized systematics? */

 const bool   doMargMH             = true,    // marginalize Mass Hierarchy
 const bool   doRethrowUnphysical  = true,    // rethrow unphysical systematic toys

 const char * dmsqscheme           = "abs23",      /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy            = "inverted",     /* normal or inverted */
 const char * mixparamchoice       = "sinsqtheta", /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant        = "one",       /* one==kNuOctantOne two==kNuOctantTwo */

 const char * xmlfile = "fakedata/201411.fakedata/fakedata_201411_0p170e21_A.xml",  // input XML file 
 bool using_xml_file  = true
)
{

  FileStat_t dummy;

  TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  TString xmlfile2 = Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile);

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

  if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) 
  {
    xmlfile2 = xmlfile;
    if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) {
      cerr<<"Data input file not found: "<<xmlfile<<endl;
      exit(-1);
    }
  }
  
  int ifakedata = 0;
  SuperkObservation * skobs_fakedata = new SuperkObservation;
  skobs_fakedata->Read(xmlfile2.Data(),ifakedata,using_xml_file);

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
  // free param, ntoys, # pt to scan, min par, max par (if min and max not passed, use default ones for that par)
  FitMarginalization *fitmarg = new FitMarginalization(freepar,NtoysMarg,NptScan,freeparmin,freeparmax); 

  // Set fit parameter list
  fitmarg->SetFitParamList(fit_param_list);
  fitmarg->PrintPhysHypothesis();
      
  // Set Mass Hierarchy marginalization
  fitmarg->SetMargMassHierarchy(doMargMH); // false by default 
  
  if(!readtoys){

    // Decide whether to rethrow or not the unphysical systematic toys 
    fitmarg->SetRethrowSystUnphysical(doRethrowUnphysical);
    
    // Produce toys and fill marginalized chisq

    // Set priors for not fixed oscillation parameters
    fitmarg->SetPrior(kPhypMixing23Bar,"flat",0,1);
    fitmarg->SetPrior(kPhypMixing13Bar,"flat",0,1);
    fitmarg->SetPrior(kPhypMixing12Bar,"flat",0,1);
    fitmarg->SetPrior(kPhypDm23sqBar,"flat",conv::kMinFitDm23sq,conv::kMaxFitDm23sq);
    fitmarg->SetPrior(kPhypDm12sqBar,"flat",conv::kMinFitDm12sq,conv::kMaxFitDm12sq);
    fitmarg->SetPrior(kPhypCPVBar,"flat",-TMath::Pi(),+TMath::Pi());

    // Make toys and store in binary files
    fitmarg->MakeToysMarginalization(tweak_systematics,skexp);   // scan all free parameter range
    fitmarg->FillChisqMarg(skobs_fakedata);
  }

  else{
       
    // Read toys from binary files and fill marginalized chisq

    for(int ifile=0;ifile<nfiles;ifile++){
      
      cout << endl;
      cout << "////////////////////////////////" << endl;
      cout << endl;
      cout << "Read file ID " << ifile << endl; 
      cout << endl;

      fitmarg->ReadToysFromBinary(ifile,infilename);
      fitmarg->FillChisqMarg(skobs_fakedata); 
    }
  
  }
    
  // Find minimum of marginalized chisq and find best-fit values
  fitmarg->Minimize(); // Find the minimum of the TH1F chisqmarg
  fitmarg->PrintBestFit();
  fitmarg->GetChisqMarg()->Draw();
  
  return;
}
