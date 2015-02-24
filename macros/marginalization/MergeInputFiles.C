

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
 const bool   readtoys             = false,

 const char * tag                  = "ciao",
 //const char * tag                  = "MH_FitMarg_mixing23_1p7E20_prefit",
 
 const int    nfiles               = 11,         // # of files to read

 const char * config_file          = "analysis_configuration_1.7E20_prefit.xml",   /* analysis configuration file */

 const char * fit_param_list       = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",

 const int    NtoysMarg            = 10,   // # toys marginalization 

 const char * freeparname_1        = "mixing23Bar",  // free parameter
 const int    NptScan_1            = 201,            // # pts of free parameter to scan
 const int    freeparmin_1         = 0,              // minimum value of free parameter
 const int    freeparmax_1         = 1,              // maxinum value of free parameter

 const char * freeparname_2        = "dmsq23Bar",  // free parameter
 const int    NptScan_2            = 201,          // # pts of free parameter to scan
 const int    freeparmin_2         = conv::kMinFitDm23sq, // minimum value of free parameter
 const int    freeparmax_2         = conv::kMaxFitDm23sq, // maxinum value of free parameter

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


  
  PhysParamId_t freepar_1;
  string sfreeparname_1(freeparname_1);
  
  // NuBar mode
  if(sfreeparname_1=="mixing23Bar")      freepar_1 = conv::kPhypMixing23Bar;
  else if(sfreeparname_1=="mixing13Bar") freepar_1 = conv::kPhypMixing13Bar;
  else if(sfreeparname_1=="dmsq23Bar")   freepar_1 = conv::kPhypDm23sqBar;
  else if(sfreeparname_1=="dmsq12Bar")   freepar_1 = conv::kPhypDm12sqBar;
  else if(sfreeparname_1=="mixing12Bar") freepar_1 = conv::kPhypMixing12Bar;
  else if(sfreeparname_1=="deltacpBar")  freepar_1 = conv::kPhypCPVBar;
  // Nu mode
  else if(sfreeparname_1=="mixing23")    freepar_1 = conv::kPhypMixing23;
  else if(sfreeparname_1=="mixing13")    freepar_1 = conv::kPhypMixing13;
  else if(sfreeparname_1=="dmsq23")      freepar_1 = conv::kPhypDm23sq;
  else if(sfreeparname_1=="dmsq12")      freepar_1 = conv::kPhypDm12sq;
  else if(sfreeparname_1=="mixing12")    freepar_1 = conv::kPhypMixing12;
  else if(sfreeparname_1=="deltacp")     freepar_1 = conv::kPhypCPV;
  //
  else{
    cout << "Wrong PhysParamId" << endl;
    exit(1);
  }
  
  // Implement the constructor giving the following input parameters:
  // free param, ntoys, # pt to scan, min par, max par (if min and max not passed, use default ones for that par)
  FitMarginalization *fitmarg = new FitMarginalization(freepar_1,NtoysMarg,NptScan_1,freeparmin_1,freeparmax_1); 

  // Set fit parameter list
  fitmarg->SetFitParamList(fit_param_list);
  fitmarg->PrintPhysHypothesis();
      
  // Set Mass Hierarchy marginalization
  fitmarg->SetMargMassHierarchy(doMargMH); // false by default 


  for(int ifile=0;ifile<nfiles;ifile++){
    
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
