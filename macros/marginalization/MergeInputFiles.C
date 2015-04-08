

//
// ------------------------------------------------------------------------------------------------           
//
// Read all the input binary files with the toys and save in a single binary file.
// Useful since running toys in parallel jobs takes less time but reading more files
// when performing the fit takes more time.
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------


using namespace std;
using namespace valor;
using namespace valor::conv;

void MergeInputFiles(

 const int    first_bin_freepar0   = -999,   // first freepar of the sub-range
 const int    last_bin_freepar0    = -999, // last freepar of the sub-range
 const int    first_bin_freepar1   = -999,   // first freepar of the sub-range
 const int    last_bin_freepar1    = -999, // last freepar of the sub-range

 //const char * tag                  = "mixing13_Vs_dcp_4p0E20_MHmarg_postfit",
 const char * tag                  = "submit_mixing13_4p0E20_postfit",

 const int    nfiles               = 201,         // # of files to read

 //const char * config_file          = "analysis_configuration_1.7E20_prefit.xml",   
 //const char * config_file          = "analysis_configuration_1.7E20.xml",   
 const char * config_file          = "analysis_configuration_run1to6_test_app.xml",

 //const char * fit_param_list       = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",
 const char * fit_param_list       = "nuosc3fbar_nueapp_fixnonbar_fix4+nuisance",
 //const char * fit_param_list       = "nuosc3fbar_fixnonbar_12fix+nuisance",

 const int    NtoysMarg            = 10000,   // # toys marginalization 

 const char * freeparname0        = "mixing13Bar",  // free parameter
 const double freeparmin0         = 0,              // minimum value of free parameter
 const double freeparmax0         = 0.3,              // maxinum value of free parameter
 //const char * freeparname0        = "mixing23Bar",  // free parameter
 //const double freeparmin0         = 0,              // minimum value of free parameter
 //const double freeparmax0         = 1,              // maxinum value of free parameter
 const int    NptScan0            = 201,            // # pts of free parameter to scan

 //const char * freeparname1          = "",  // free parameter
 const char * freeparname1          = "deltacpBar",  // free parameter
 const double freeparmin1           = -3.14159265358979312e+00,  // minimum value of free parameter
 const double freeparmax1           = +3.14159265358979312e+00,  // maxinum value of free parameter
 //const char * freeparname1          = "dmsq23Bar",  // free parameter
 //const double freeparmin1           = 0.0002,     // minimum value of free parameter
 //const double freeparmax1           = 0.02,     // maxinum value of free parameter
 //const int    NptScan1              = 201,  // # pts of free parameter to scan
 const int    NptScan1              = 1,  // # pts of free parameter to scan

 const bool   doMargMH             = true,    // marginalize Mass Hierarchy
 const bool   doRethrowUnphysical  = true,    // rethrow unphysical systematic toys
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
  fitmarg->Init(skexp); // Before need SetPhysParFree() to initialize the arrays

  
  // Set fit parameter list
  fitmarg->SetFitParamList(fit_param_list);
  fitmarg->PrintPhysHypothesis();
      
  // Set Mass Hierarchy marginalization
  fitmarg->SetMargMassHierarchy(doMargMH); // false by default 

  for(int ifile=0;ifile<nfiles;ifile++){
  //for(int ifile=10;ifile<22;ifile++){
  //for(int ifile=0;ifile<10;ifile++){
    
    cout << endl;
    cout << "////////////////////////////////" << endl;
    cout << endl;
    cout << "Read file ID " << ifile << endl; 
    cout << endl;
    
    fitmarg->ReadToysFromBinary(ifile,tag); // Read each binary file
    fitmarg->MergeInputFiles(); // Store toys of current binary file
  }
  
  fitmarg->WriteToys2BinaryMerge(tag); // Write all the stored toys in a final binary file

  return;
}
