
// ------------------------------------------------------------------------------------------------           
//
// Macro to perform fit marginalizing the systematic and oscillation parameters
// Run toys for marginalization for all or sub-range of the free parameter and fill the Chisq.
// In the end the minimum of the marginalized Chisq distribution is found.
// Toys for different sub-ranges of the free paramters can be run in parallel jobs.
// The other possibility is to produce toys with build_marg_toys.C and store in binary files.
// Then this macro reads the toys from binary files, fill the marginalized Chisq and minimize. 
//
// Davide Sgalaberna                                                                                              
//                                                                                           
// ------------------------------------------------------------------------------------------------


using namespace std;
using namespace valor;
using namespace valor::conv;

void run_fit_marginalization
(
 const bool   readtoys             = true,

 // Parameter for chisq minimization
 const bool   minimizeHisto        = true, // if not all the chisq histo is filled
 const int    MaxIterFit           = 2, 
 const int    MaxNfit              = 3,
 
 const int    strategy             = 0,    // 0 --> use steps, 1 --> use intervals

 // Step strategy
 const double start_X              = 0.102,
 const double start_Y              = -1.66504,
 const string StrategyStartPt      = "GridFix", //"Single", // "GridRnd", "GridFix"
 const int    NStartPt_X           = 12,  // # of throws for starting pt extraction
 const int    NStartPt_Y           = 4,  // # of throws for starting pt extraction

 // Interval strategy
 const bool   doUncorrXY           = false, // DON'T USE FOR sin2th13 Vs dcp
 const bool   RedoBorderFit        = true,
 const int    MaxNfitBorder        = 1,

 //const double MinStep0             = 0.,     // min x value for fit steps
 //const double MaxStep0             = 1.,     // max x value for fit steps
 const double MinStep0             = 0.,     // min x value for fit steps
 const double MaxStep0             = 0.3,     // max x value for fit steps
 
 //const double MinStep1             = 0.0002, // min y value for fit steps
 //const double MaxStep1             = 0.02,   // max y value for fit steps
 //const double MinStep1             = -TMath::Pi(), // min y value for fit steps
 const double MinStep1             = -TMath::Pi(), // min y value for fit steps
 const double MaxStep1             = +TMath::Pi(),   // max y value for fit steps
 //

 //const char * outfilename          = "ChisqMarg_mixing23_Vs_dmsq23_Nbins201x201",
 const char * outfilename          = "ChisqMarg_mixing23_Vs_dmsq23_Nbins21x21",
 //const char * outfilename          = "ChisqMarg_prova",
 
 const int    IDoutfile            = 0,
 const int    first_bin_freepar0   = -999,   // first freepar of the sub-range
 const int    last_bin_freepar0    = -999, // last freepar of the sub-range
 const int    first_bin_freepar1   = -999,   // first freepar of the sub-range
 const int    last_bin_freepar1    = -999, // last freepar of the sub-range

 //const char * infilename           = "../../inputs/marginalization/toys/mixing23_1p7E20_postfit_NoRethrow_Merged",
 //const char * infilename           = "submit_mixing23_1p7E20_postfit_NoRethrow",
 //const char * infilename           = "submit_mixing23_1p7E20_postfit_NoRethrow_Merged", // test
 //const char * infilename           = "../../inputs/marginalization/toys/mixing23_dmsq32_1p7E20_postfit_NoRethrow_Merged",
 //const char * infilename           = "MH_FitMarg_mixing13Vsdeltacp_1p7E20_prefit",
 const char * infilename           = "MH_FitMarg_mixing13Vsdeltacp_1p7E20_201x201_prefit",

 const int    nfiles               = 1,         // # of files to read

 //const char * config_file          = "analysis_configuration_1.7E20_prefit.xml",  
 //const char * config_file          = "analysis_configuration_1.7E20.xml",  
 const char * config_file          = "analysis_configuration_run1to6_test_app.xml",
 
 //const char * fit_param_list       = "nuosc3fbar_fixnonbar_fix5+nuisance",
 const char * fit_param_list       = "nuosc3fbar_fixnonbar_12fix+nuisance",
 //const char * fit_param_list      = "nuosc3fbar_numudisapp_constrainall+nuisance",  
 //const char * fit_param_list       = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",

 const int    NtoysMarg            = 1,   // # toys marginalization 
 //const int    NtoysMarg            = 2,   // # toys marginalization 
 
 //const char * freeparname0        = "mixing23Bar",  // free parameter
 //const double freeparmin0         = 0.,              // minimum value of free parameter
 //const double freeparmax0         = 1.,            // maxinum value of free parameter
 //const int    NptScan0            = 21, //2, //21,        // # pts of free parameter to scan
 const char * freeparname0        = "mixing13Bar",  // free parameter
 const double freeparmin0         = 0.,              // minimum value of free parameter
 const double freeparmax0         = 0.3, //0.05, //0.3,   // maxinum value of free parameter
 const int    NptScan0            = 201, //2, //21,        // # pts of free parameter to scan
 
 //const char * freeparname1        = "",  // free parameter
 const char * freeparname1        = "deltacpBar",  // free parameter
 const double freeparmin1         = -TMath::Pi(), // minimum value of free parameter
 const double freeparmax1         = +TMath::Pi(), // maxinum value of free parameter
 const int    NptScan1            = 201, //21,      // # pts of free parameter to scan
 //const char * freeparname1        = "dmsq23Bar",  // free parameter
 //const double freeparmin1         = 0.0002, // minimum value of free parameter
 //const double freeparmax1         = 0.02, // maxinum value of free parameter
 //const int    NptScan1            = 21, //21,      // # pts of free parameter to scan

 int         tweak_systematics    = 1,       /* generate toy MC with randomized systematics? */

 const bool  doMargMH             = true,    // marginalize Mass Hierarchy
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

    if(!minimizeHisto){
      fitmarg->FillChisqMarg(skobs_fakedata);
      fitmarg->MinimizeHisto();
    }
    else{
      // Max # of iterations per fit
      fitmarg->SetMaxIterationSingleFit(MaxIterFit); 
      // Max # of fits max # of iterations
      fitmarg->SetMaxFitNumber(MaxNfit);   

      // Use chisq minimization with steps
      if(strategy==0){
	fitmarg->SetMinimInterval(MinStep0,MaxStep0,MinStep1,MaxStep1);
	fitmarg->SetStartBin(start_X,start_Y); 
	fitmarg->SetNumStartBin(NStartPt_X,NStartPt_Y);
	fitmarg->SetExtractStartPt(StrategyStartPt);
	fitmarg->MinimizeChisqStep(skobs_fakedata);      
      }

      // Use chisq minimization with intervals
      else if(strategy==1){
	// Redo the fit if minimum on the border  
	fitmarg->SetBorderDoFit(RedoBorderFit); // true by default
	// Max # of fits max # of iterations if the first minimum is on the border
	fitmarg->SetBorderMaxFitNumber(MaxNfitBorder); // 2*Ndim by default
	fitmarg->SetMinimInterval(MinStep0,MaxStep0,MinStep1,MaxStep1);
	fitmarg->SetMinimUncorrXY(doUncorrXY);
	fitmarg->MinimizeChisqInterval(skobs_fakedata);      
      }

      cout << "# of calls of CalcChisqMargPt(): " 
	   << fitmarg->GetNCallsCalcChisqMargPt() << endl;
      
      TCanvas *cChisqStep = new TCanvas("cChisqStep","cChisqStep");
      if(fitmarg->IsFit2D())
	fitmarg->GetChisqStep_2D()->DrawCopy("colz");
      else
	fitmarg->GetChisqStep_1D()->DrawCopy("");
    }
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

      if(!minimizeHisto){
	fitmarg->FillChisqMarg(skobs_fakedata);
	fitmarg->MinimizeHisto();
      }
      else{
	// Max # of iterations per fit
	fitmarg->SetMaxIterationSingleFit(MaxIterFit); 
	// Max # of fits max # of iterations
	fitmarg->SetMaxFitNumber(MaxNfit);   
	
	// Use chisq minimization with steps
	if(strategy==0){
	  fitmarg->SetMinimInterval(MinStep0,MaxStep0,MinStep1,MaxStep1);
	  fitmarg->SetStartBin(start_X,start_Y); 
	  fitmarg->SetNumStartBin(NStartPt_X,NStartPt_Y);
	  fitmarg->SetExtractStartPt(StrategyStartPt);
	  fitmarg->MinimizeChisqStep(skobs_fakedata);      
	}

	// Use chisq minimization with intervals
	else if(strategy==1){
	  // Redo the fit if minimum on the border  
	  fitmarg->SetBorderDoFit(RedoBorderFit); // true by default
	  // Max # of fits max # of iterations if the first minimum is on the border
	  fitmarg->SetBorderMaxFitNumber(MaxNfitBorder); // 2*Ndim by default
	  fitmarg->SetMinimInterval(MinStep0,MaxStep0,MinStep1,MaxStep1);
	  fitmarg->SetMinimUncorrXY(doUncorrXY);
	  fitmarg->MinimizeChisqInterval(skobs_fakedata);      
	}
	
	cout << "# of calls of CalcChisqMargPt(): " 
	     << fitmarg->GetNCallsCalcChisqMargPt() << endl;

	TCanvas *cChisqStep = new TCanvas("cChisqStep","cChisqStep");
	if(fitmarg->IsFit2D())
	  fitmarg->GetChisqStep_2D()->DrawCopy("colz");
	else
	  fitmarg->GetChisqStep_1D()->DrawCopy("");
      }
      
    }  
  }
 
  // Find minimum of marginalized chisq and find best-fit values
  fitmarg->PrintBestFit();
  
  TCanvas *cChisqMarg = new TCanvas("cChisqMarg","cChisqMarg");
  if(fitmarg->IsFit2D()){
    fitmarg->GetChisqMarg2D()->DrawCopy("colz");
    //fitmarg->GetChisqMarg2D()->DrawCopy("colz text");
    //fitmarg->GetChisqMarg2D()->DrawCopy("LEGO");
  }
  else{
    double min = fitmarg->GetChisqBestFit();
    fitmarg->GetChisqMarg()->SetMinimum(min);
    fitmarg->GetChisqMarg()->DrawCopy();
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
