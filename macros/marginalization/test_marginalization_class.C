
using namespace std;
using namespace valor;
using namespace valor::conv;

void test_marginalization_class
(
 const bool   readtoys = true;
 const char * outfilename            = "framework",
 const char * infilename             = "framework",

 const char * config_file      = "analysis_configuration_1.7E20_prefit.xml",   /* analysis configuration file */

 const char * fit_param_list    = "nuosc3fbar_numudisapp_fixnonbar_fix4+nuisance",

 int          tweak_systematics       = 1,                                /* generate toy MC with randomized systematics? */

 const char * dmsqscheme              = "abs23",                             /* abs23==kNuDmsq12andAbs23 abs31==kNuDmsq12andAbs31 fogli==kNuDmsqFogliLisi */
 const char * hierarchy               = "inverted",                            /* normal or inverted */
 const char * mixparamchoice          = "sinsqtheta",                        /* sinsq2theta or sinsqtheta */
 const char * mixparamoctant          = "one",                               /* one==kNuOctantOne two==kNuOctantTwo */

 const char * xmlfile = "fakedata/201411.fakedata/fakedata_201411_0p170e21_A.xml",  // input XML file 
 bool using_xml_file  = true,

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
  //skobs_fakedata->ErecoSpectrum(kSKsRHCFCFV1RingMuLike)->Draw();
  //return;







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







  FitMarginalization *fitmarg = new FitMarginalization(conv::kPhypMixing23Bar,10,21); // free param, ntoys, npt param scan
  //fitmarg->Init();
  //fitmarg->CleanUp();
  
  // Set fit parameter list
  fitmarg->SetFitParamList(fit_param_list);
  //cout << *fitmarg->GetFitParamList() << endl;

  fitmarg->PrintPhysHypothesis();




  if(!readtoys){
    
    // Make toys and store in binary files
    
    // Produce toys for different sub-ranges of the free parameter and store in different files.
    // Useful for running on cluster when many toys are needed (e.g. 10k x 200)
    
    fitmarg->MakeToysMarginalization(tweak_systematics,skexp,0,10); 
    
    // Toys for marginalization are stored in binary files (reading and writing very fast)
    // Several files can be written for different sub-ranges of the free parameter each
    // Sub-range bins are given as input in MakeToysMarginalization() and stored in a second binary file for now
    fitmarg->WriteToys2Binary(0,infilename);
    
    fitmarg->MakeToysMarginalization(tweak_systematics,skexp,11,14);
    fitmarg->WriteToys2Binary(1,infilename);
    
    fitmarg->MakeToysMarginalization(tweak_systematics,skexp,15,20);
    fitmarg->WriteToys2Binary(2,infilename);
    
    // Once all the toys are extracted (attention not to make any overlap in the indices) 
    // the marginalized chisq is computed by averaging over all the stored toys. 
    // The fake data set is taken as input.
    // This method doesn't minimize yet
    fitmarg->FillChisqMarg(skobs_fakedata);

  }

  else{
    
    // Read toys from output binary files 
    // When reading different files for free param sub-ranges, ChisqMarg must be filled after reading each file 
    // Not needed when making toys
    
    cout << "/////////////////////////" << endl;
    fitmarg->ReadToysFromBinary(0,infilename);

    // Fill the ChisqMarg after reading each toy file. Otherwise we need to copy each element to another array
    // increasing the number of loops. A TH1F chisqmarg is filled as well.
    fitmarg->FillChisqMarg(skobs_fakedata); 
  
    cout << "/////////////////////////" << endl;
    fitmarg->ReadToysFromBinary(1,infilename);
    fitmarg->FillChisqMarg(skobs_fakedata);

    cout << "/////////////////////////" << endl;
    fitmarg->ReadToysFromBinary(2,infilename);
    fitmarg->FillChisqMarg(skobs_fakedata); 
    
  }
  
  
  // Find minimum of marginalized chisq and find best-fit values
  fitmarg->Minimize(); // Find the minimum of the TH1F chisqmarg
  fitmarg->PrintBestFit();
  fitmarg->GetChisqMarg()->Draw();
  






  return;
}
