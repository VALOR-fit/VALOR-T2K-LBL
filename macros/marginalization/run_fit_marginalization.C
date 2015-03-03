
using namespace std;
using namespace valor;
using namespace valor::conv;

void run_fit_marginalization
(
 const bool   readtoys             = false,

 const char * outfilename          = "ChisqMarg",
 
 const int    IDoutfile            = 0,
 const int    first_bin_freepar0   = -999,   // first freepar of the sub-range
 const int    last_bin_freepar0    = -999, // last freepar of the sub-range
 const int    first_bin_freepar1   = -999,   // first freepar of the sub-range
 const int    last_bin_freepar1    = -999, // last freepar of the sub-range

 //const char * infilename           = "../../inputs/marginalization/toys/mixing23_1p7E20_postfit_NoRethrow_Merged",
 //const char * infilename           = "MH_FitMarg_mixing23_1p7E20_prefit",
 //const char * infilename           = "MH_FitMarg_mixing23_1p7E20_prefit_Merged",
 //const char * infilename           = "submit_mixing23_1p7E20_postfit_NoRethrow",
 const char * infilename           = "submit_mixing23_1p7E20_postfit_NoRethrow_Merged",
 //const char * infilename           = "../../inputs/marginalization/toys/mixing23_dmsq32_1p7E20_postfit_NoRethrow_Merged",

 const int    nfiles               = 1,         // # of files to read

 //const char * config_file          = "analysis_configuration_1.7E20_prefit.xml",  
 //const char * config_file          = "analysis_configuration_1.7E20.xml",  
 const char * config_file          = "analysis_configuration_run1to6_test_app.xml",
 
 //const char * fit_param_list       = "nuosc3fbar_fixnonbar_fix5+nuisance",
 const char * fit_param_list       = "nuosc3fbar_fixnonbar_12fix+nuisance",
 //const char * fit_param_list      = "nuosc3fbar_numudisapp_constrainall+nuisance",  
 //const char * fit_param_list       = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",

 const int    NtoysMarg            = 1,   // # toys marginalization 
 
 const char * freeparname0        = "mixing13Bar",  // free parameter
 //const char * freeparname0        = "mixing23Bar",  // free parameter
 const double freeparmin0         = 0.,              // minimum value of free parameter
 const double freeparmax0         = 0.3,            // maxinum value of free parameter
 const int    NptScan0            = 21,            // # pts of free parameter to scan

 const char * freeparname1        = "deltacpBar",  // free parameter
 const double freeparmin1         = -TMath::Pi(), // minimum value of free parameter
 const double freeparmax1         = +TMath::Pi(), // maxinum value of free parameter
 //const char * freeparname1        = "dmsq23Bar",  // free parameter
 //const double freeparmin1         = 0.0002, // minimum value of free parameter
 //const double freeparmax1         = 0.02, // maxinum value of free parameter
 const int    NptScan1            = 21,          // # pts of free parameter to scan

 int         tweak_systematics    = 1,       /* generate toy MC with randomized systematics? */

 const bool  doMargMH             = false,    // marginalize Mass Hierarchy
 const bool  doRethrowUnphysical  = false,    // rethrow unphysical systematic toys

 const char * dmsqscheme          = "abs23",      /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy           = "normal",     /* normal or inverted */
 const char * mixparamchoice      = "sinsqtheta", /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant      = "one",       /* one==kNuOctantOne two==kNuOctantTwo */

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
  
  if(!readtoys){

    // Decide whether to rethrow or not the unphysical systematic toys 
    fitmarg->SetRethrowSystUnphysical(doRethrowUnphysical);
    
    // Produce toys and fill marginalized chisq

    // Set priors for not fixed oscillation parameters
    fitmarg->SetPrior(kPhypMixing23Bar,"flat",0,1);
    fitmarg->SetPrior(kPhypMixing13Bar,"flat",0.07,0.15);
    fitmarg->SetPrior(kPhypMixing12Bar,"flat",0,1);
    fitmarg->SetPrior(kPhypDm23sqBar,"flat",conv::kMinFitDm23sq,conv::kMaxFitDm23sq);
    fitmarg->SetPrior(kPhypDm12sqBar,"flat",conv::kMinFitDm12sq,conv::kMaxFitDm12sq);
    fitmarg->SetPrior(kPhypCPVBar,"flat",-TMath::Pi(),+TMath::Pi());

    fitmarg->MakeToysMarginalization(tweak_systematics,skexp); 

    //fitmarg->GetToyMargHisto(conv::kSKsRHCFCFV1RingMuLike,0,0,0,0)->Draw();
    //fitmarg->GetToyMargHisto(conv::kSKsRHCFCFV1RingELike,0,0,0,0)->Draw();
    
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
  fitmarg->Minimize(); 
  fitmarg->PrintBestFit();
  
  if(fitmarg->IsFit2D()){
    fitmarg->GetChisqMarg2D()->Draw("colz");
  }
  else{
    fitmarg->GetChisqMarg()->Draw();
  }
  
  
  // Store ChisqMarg (full or partial) in output file
  
  char outname[200];
  sprintf(outname,"%s_%d.root",outfilename,IDoutfile);
  
  cout << outname << endl;

  TFile *fout = new TFile(outname,"RECREATE");
  fitmarg->GetChisqMarg()->Write();
  fitmarg->GetChisqMarg2D()->Write();
  fout->Close();
  
  return;
}
