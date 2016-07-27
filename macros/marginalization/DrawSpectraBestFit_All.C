
//
// ------------------------------------------------------------------------------------------------
//
// Make the best-fit spectra with marginalization
// Takes as input the following files:
// - binary file with likelihood weights from the fit w/ marginalization
// - binary file with marginalization toys
//
// Davide Sgalaberna
//
// ------------------------------------------------------------------------------------------------
//

#include "Build.h"
#include "AnalysisOptions.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

//#define USE_UNCOMPRESSED_VARRAY_FORMAT // Disable compressed files in $VALOR_FMWK/Utils/src/Utilities/VALORArray.cxx


void DrawSpectraBestFit_All
(
 const char *tag = "fitmarg_prova",

 const bool draw_histobin = false,

 const char * xmlfile      = "",
 const char * rootfile     = "data_run1to6/SKdata_Run16.root",
 const bool using_xml_file = false,
 int   iexperiment         = 0, // data (or fake data) event to show                    
 
 const char * config_file  = "analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_fitbias.xml",

 const char * sk_sample    = "RHC1Re",   /* "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */

 bool SystMode = 1,  // 0 --> GetMode, 1 --> GetMean
 bool WeightMode = 1, // 0 --> All weights are 1, 1 --> Likelihoods are weights
 
 ////////////////////// MARGINALIZATION PARAMETERS ///////////////////

 const char * likeinfilename = "liketoy_fit_fit_marg_PROVA_nexp_1_0.450_-1.601-job0_toy1_x_oscpar81_oscpar51_MHmarg.bin",
 //const char * likeinfilename = "liketoy_fit_fitmarg_joint_betaBar_prova_nexp_1_1.000_-999.000-job0_toy10000_x_oscpar2_oscpar1_MHmarg.bin",
 
 const double InputPar0 = 0.405,  // Input value of parameter 0 
 const double InputPar1 = -1.2566371, // Input value of parameter 1

 //const char * fit_param_list = "nuosc3f_normbar+nuisance", // betabar
 const char * fit_param_list = "nuosc3f+nuisance", // dcp

 //const char * infilename = "../TECH_NOTE_14_3_16/BETABAR/TOYS/joint_betabar_react_MargMH_10kToys/joint_betabar_niwg2015c_react_MargMH_10kToys",
 const char * infilename = "prova",
 
 const int    nfiles               = 1, //201, // 1         // # of files to read

 const int    NtoysMarg            = 1, //10000,   // # toys marginalization 
 
 // const char * freeparname0        = "deltacp",  // free parameter
 // const double freeparmin0         = -3.14159265358979312e+00,         // minimum value of free parameter
 // const double freeparmax0         = +3.14159265358979312e+00,         // maxinum value of free parameter
 // const int    NptScan0            = 201,                              // # pts of free parameter to scan
 // const char * freeparname1        = "",  // free parameter
 // const double freeparmin1         = 0.,  // minimum value of free parameter
 // const double freeparmax1         = 0.,  // maxinum value of free parameter
 // const int    NptScan1            = 1,   // # pts of free parameter to scan

 // const char * freeparname0        = "betaBar",  // free parameter
 // const double freeparmin0         = 0.,         // minimum value of free parameter
 // const double freeparmax0         = 1.,         // maxinum value of free parameter
 // const int    NptScan0            = 2,          // # pts of free parameter to scan
 // const char * freeparname1        = "",         // free parameter
 // const double freeparmin1         = 0.,         // minimum value of free parameter
 // const double freeparmax1         = 0.,         // maxinum value of free parameter
 // const int    NptScan1            = 1,          // # pts of free parameter to scan

 const char * freeparname0        = "mixing23",  // free parameter
 const double freeparmin0         = 0.3,         // minimum value of free parameter
 const double freeparmax0         = 0.7,         // maxinum value of free parameter
 const int    NptScan0            = 81,          // # pts of free parameter to scan
 const char * freeparname1        = "deltacp",         // free parameter
 const double freeparmin1         = -3.14159265358979312e+00,         // minimum value of free parameter
 const double freeparmax1         = +3.14159265358979312e+00 ,         // maxinum value of free parameter
 const int    NptScan1            = 51,          // # pts of free parameter to scan

 const bool  doMargMH             = true,    // marginalize Mass Hierarchy
 const char * hierarchy_marg      = "normal",     /* normal or inverted */

 //const char * fileosctoys         = "fit_fitmarg_joint_dcp_react_nh_nexp_1_-1.601_-999.000-job0-_0_0.root"
 )
{
  char name[256];

  FileStat_t dummy;
  
  TString opt(fit_param_list); 
  
  // Analysis configuration
  TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  
  TString xmlfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile));
  TString rootfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), rootfile));

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

  bool dotoys = false;
  if (!dotoys && gSystem->GetPathInfo(xmlfile2.Data(), dummy)){
    xmlfile2 = xmlfile;
    if (gSystem->GetPathInfo(xmlfile2.Data(), dummy)) {
      cerr<<"Data input file not found: "<<xmlfile<<endl;
      exit(-1);
    }
  }

  cout << "- Toy-MC" << endl;
  if(using_xml_file) cout << "   - Reading experiments from: " << xmlfile2 << endl; 
  else cout << "   - Reading experiments from: " << rootfile2 << endl; 
  cout << "   - Using experiment: " << iexperiment << endl;
  

  SuperkObservation * skobs_fakedata = new SuperkObservation;
  if(using_xml_file){
    cout << "Use xml file" << endl;
    skobs_fakedata->Read(xmlfile2.Data(),iexperiment,using_xml_file);
  }
  else{
    cout << "Use root file" << endl;
    skobs_fakedata->Read(rootfile2.Data(),iexperiment,using_xml_file);
  }


  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config.Data());
 
  // Fit param list for fakedata sets
  FitParamList * fparm  = util::NamedParameterList(fit_param_list); 
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used for fakedata sets: " << util::AsString(hyp) << endl;

  // Set Marginalization

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
  else if(sfreeparname0=="betaBar")     freepar0 = conv::kPhypBetaBar;
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
  else if(sfreeparname1=="betaBar")     freepar1 = conv::kPhypBetaBar;
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
  fitmarg->SetSubRange(subminbin,submaxbin); // By default it scans all the range
  fitmarg->Init(skexp); // Before need SetPhysParFree() to initialize the arrays

  // Set if to store full likelihood for each toy (after fitmarg->Init())
  fitmarg->SetStoreLikeWeights(true); 
   
  // Set Mass Hierarchy marginalization
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(hierarchy_marg); // needed to build the right binary file name 
  fitmarg->SetMargMassHierarchy(doMargMH); // false by default 
    
  // Set fit parameter list
  fitmarg->SetFitParamList(fit_param_list);
  fitmarg->PrintPhysHypothesis();
  
  
  // Build the map fileID - sub-range
  
  fitmarg->SetMaxIterationSingleFit(2);  // Max # of iterations per fit
  fitmarg->SetMaxFitNumber(3); // Max # of fits max # of iterations
  
  fitmarg->SetTagFileName(infilename);
  fitmarg->SetFileID(nfiles);

  // Read fit toy likelihoods from binary file
  fitmarg->ReadLikeToysFromBinary(likeinfilename);


  // Get input oscillation parameters
  
  int ipar0 = conv::kBadResult;
  int ipar1 = conv::kBadResult;
  double iparbincenter0 = conv::kBadResult;
  double iparbincenter1 = conv::kBadResult;
  
  ipar0 = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->FindBin(InputPar0)-1; // bin of the parameter (histo bin = param bin + 1)
  iparbincenter0 = fitmarg->GetMapChisqRangeFileID_2D()->GetXaxis()->GetBinCenter(ipar0+1);
  
  if(fitmarg->IsFit2D()){
    ipar1 = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->FindBin(InputPar1)-1; // bin of the parameter (histo bin = param bin + 1)
    iparbincenter1 = fitmarg->GetMapChisqRangeFileID_2D()->GetYaxis()->GetBinCenter(ipar1+1);
  }
  else{
    ipar1 = 0;
    iparbincenter1 = conv::kBadResult;
  }
  
  cout << endl;
  cout << "Input oscillation parameters: " << endl;
  cout << " - parameter 0: " << freeparname0 << " --> " << InputPar0 << endl;
  cout << "   grid bin is " << ipar0 << endl;
  cout << "   grid bin center = " << iparbincenter0 << endl;
  cout << " - parameter 1: " << freeparname1 << " --> " << InputPar1 << endl;
  cout << "   grid bin is " << ipar1 << endl;
  cout << "   grid bin center = " << iparbincenter1 << endl;
  cout << endl;
  
  // Set variables to read the fit toy likelihoods
  
  double truepar0 = -999;
  double truepar1 = -999;
  truepar0 = fparm->InitValue(freepar0);
  if(fitmarg->IsFit2D()) truepar1 = fparm->InitValue(freepar1);	

  // Declare histograms
  
  TH1D *h_bestfit_marg[conv::kns];
  TH1D *hErecoToy[conv::kns];
  TH1D **h_bin[conv::kns];

  // Loop over the SK samples
  for(int ikns=0;ikns<conv::kns;ikns++){
    
    //conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
    SuperkSample_t s = kSKSamples[ikns];
    if(s == kSKsUndefined) continue;
  
    cout << endl;
    cout << "SK sample: " << util::AsStringCompact(s) << endl;
    cout << endl;

    // Get Ereco toys
    
    // Declare the histograms
    const int n_ereco_bins = util::NumErecoBins(s);  
    const int n_q2_bins = util::NumQ2Bins(s);  
    TH1F* h_marg_toy;
    
    //TH1D *h_bestfit_marg = new TH1D("h_bestfit_marg","h_bestfit_marg",util::NumErecoBins(s),util::ErecoBinEdges(s));
    sprintf(name,"h_bestfit_marg_%s",util::AsStringCompact(s).c_str());
    h_bestfit_marg[ikns] = new TH1D(name,name,util::NumErecoBins(s),util::ErecoBinEdges(s));
    
    //
    
    






    //int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
        
    //
    
    const int NErecBins = util::NumErecoBins(s);
    double **NentriesErecoToy = new double*[NErecBins];
    for(int i=0;i<NErecBins;i++) NentriesErecoToy[i] = new double[NtoysMarg];
    
    double *LikeToy = new double[NtoysMarg];
    double SumLike = 0.;
    
    // Loop over marginalization toys
    
    for(int itoy=0;itoy<NtoysMarg;itoy++){
      
#ifndef __USE_2D_TEMPLATES__
      h_marg_toy = (TH1F*)fitmarg->GetToyMargHisto(s,itoy,ipar0,ipar1);
#else
      h_marg_toy = (TH1F*)fitmarg->GetToyMargHisto2D(s,itoy,ipar0,ipar1)->ProjectionY()->Clone("h_marg_toy");
#endif
      
      for(int i=1;i<=util::NumErecoBins(s);i++){ 
	NentriesErecoToy[i-1][itoy] = (double) h_marg_toy->GetBinContent(i);
	//cout << h_marg_toy->GetBinCenter(i) << endl;
      }
      
      //if(fitmarg->IsFit2D())
      //LikeToy[itoy] = (double) fitmarg->GetLikeToyRead(InputPar0,truepar1,itoy,true); 
      //else
      //LikeToy[itoy] = (double) fitmarg->GetLikeToyRead(InputPar0,InputPar1,itoy,true); 

      LikeToy[itoy] = (double) fitmarg->GetLikeToyRead(InputPar0,InputPar1,itoy,true); 
      
      //cout << "toy " << itoy << " - like = " <<  LikeToy[itoy] << endl;
      
      if(LikeToy[itoy] < 0){
	cerr << endl;
	cerr << "Marginalization toy likelihood is negative!" << endl;
	cerr << "Toy " << itoy << " --> Likelihood = " << LikeToy[itoy] << endl;
	cerr << endl;
	exit(1);
      }
      
      SumLike += LikeToy[itoy];
      
    } // loop itoy
    
    
    // Set and declare bin histograms 
    
    double *lowedge = new double[n_ereco_bins];
    for(int i=0;i<n_ereco_bins;i++) lowedge[i] = 9999999.;
    double *upedge = new double[n_ereco_bins];
    for(int i=0;i<n_ereco_bins;i++) upedge[i] = 0.;
    
    for (int ibin = 0 ; ibin < n_ereco_bins ; ibin++) {
      for(int itoy=0;itoy<NtoysMarg;itoy++){
	double cont = NentriesErecoToy[ibin][itoy]; 
	if(cont < lowedge[ibin]){
	  lowedge[ibin] = cont;
	}
      }        
      for(int itoy=0;itoy<NtoysMarg;itoy++){
	double cont = NentriesErecoToy[ibin][itoy]; 
	if(cont > upedge[ibin]){
	  upedge[ibin] = cont;
	}      
      }
      // histograms must contain the first and last entries 
      if(lowedge[ibin]==0.) lowedge[ibin] -= 1.E-3;
      upedge[ibin] += 1.E-3;
      //cout << "bin # " << ibin << " --> lowedge = " << lowedge[ibin] << " - upedge = " << upedge[ibin] << endl;
    
    } // loop over Ereco bins
    
    
    
    // Fill bin histograms
  
    h_bin[ikns] = new TH1D*[util::NumErecoBins(s)];
    for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){
      sprintf(name,"h_bin_%s_%d",util::AsStringCompact(s).c_str(),ibin);
      h_bin[ikns][ibin] = new TH1D(name,name,50,lowedge[ibin],upedge[ibin]);
    }
    
    for(int i=0;i<util::NumErecoBins(s);i++){
      for(int itoy=0;itoy<NtoysMarg;itoy++){
	double weight = LikeToy[itoy]/SumLike;
	if(!WeightMode) weight = 1.;
	h_bin[ikns][i]->Fill(NentriesErecoToy[i][itoy],weight);
      }
    }
    
    
    // Find the mode/mean of each Ereco bin
    
    double *bin_mode = new double[n_ereco_bins];
    
    for(int iereco=0;iereco<util::NumErecoBins(s);iereco++){
      
      double y = h_bin[ikns][iereco]->GetMean();
      double ey = h_bin[ikns][iereco]->GetRMS();
      
      if(SystMode==0){ // Get the mode
	cerr << endl;
	cerr << "SystMode==1 is not implemented" << endl;
	cerr << endl;
	exit(1);
	
	// 	cout << "bin # " << iereco << " --> Integral = " << h_bin[ikns][iereco]->Integral()
	// 	     << " - Mean = " << y << " - RMS = " << ey << endl;
	// 	TF1 *f1 = new TF1("f1","gaus", y-ey,y+ey);  
	// 	TFitResultPtr r = h_bin[ikns][iereco]->Fit("f1","S","goff",y-ey,y+ey);       
	// 	bin_mode[iereco] = r.Parameter(1);
	// 	cout << bin_mode[iereco] << " - ";
	// 	bin_mode[iereco] = f1->GetMaximumX();
	// 	cout << bin_mode[iereco] << endl;
      } 
      else if(SystMode==1){ // Get the mean
	bin_mode[iereco] = y;
      }
            
      if(bin_mode[iereco]<0.) bin_mode[iereco] = 0.;      
      h_bestfit_marg[ikns]->SetBinContent(iereco+1,bin_mode[iereco]);
      
      //cout << "# of events in Ereco " << iereco << " = " << bin_mode[iereco] << endl;
    }
    
    sprintf(name,"hErecoToy_%s",util::AsStringCompact(s).c_str());
#ifndef __USE_2D_TEMPLATES__
    hErecoToy[ikns] = (TH1D*)skobs_fakedata->ErecoSpectrum(s)->Clone(name);
#else
    hErecoToy[ikns] = (TH1D*)skobs_fakedata->TwoDSpectrum(s)->ProjectionY()->Clone(name);
#endif  
   
  } // end loop over sk samples
  
   
  // Draw bin histograms

  TCanvas **c_bin[conv::kns];
  if(draw_histobin){
    for(int ikns=0;ikns<conv::kns;ikns++){
      SuperkSample_t s = kSKSamples[ikns];
      if(s == kSKsUndefined) continue;      
      c_bin[ikns] = new TCanvas*[util::NumErecoBins(s)];
      for(int iereco=0;iereco<util::NumErecoBins(s);iereco++){
	sprintf(name,"c_bin_%s_%d",util::AsStringCompact(s).c_str(),iereco);    
	c_bin[ikns][iereco] = new TCanvas(name,name);
	h_bin[ikns][iereco]->Draw();
      }
    }
  }
  
  
  
  // Draw best-fit histogram marginalizing the systematics

  TCanvas *c_bestfit_marg[conv::kns];   
  for(int ikns=0;ikns<conv::kns;ikns++){
    SuperkSample_t s = kSKSamples[ikns];
    if(s == kSKsUndefined) continue;
    sprintf(name,"c_bestfit_marg_%s",util::AsStringCompact(s).c_str());
    c_bestfit_marg[ikns] = new TCanvas(name,name);
    hErecoToy[ikns]->GetXaxis()->SetTitle("Reconstructed neutrino energy (GeV)");
    hErecoToy[ikns]->GetYaxis()->SetTitle("Entries");
    hErecoToy[ikns]->Draw("pe");
    h_bestfit_marg[ikns]->Draw("same");
  }
  
  // Store best-fit spectrum in a file
  sprintf(name,"spectrum_bestfit_%s.root",tag);
  TFile *fout = new TFile(name,"RECREATE");

  for(int ikns=0;ikns<conv::kns;ikns++){
    SuperkSample_t s = kSKSamples[ikns];
    if(s == kSKsUndefined) continue;
    hErecoToy[ikns]->Write();
    h_bestfit_marg[ikns]->Write();
  }
  fout->Close();
    
    
    
    



  //FitParamList * fparm_osc  = util::NamedParameterList("nuosc3f"); 
  // TFile *filetreein = new TFile(fileosctoys,"READ");
  // if(!filetreein->IsOpen()){
  //   cerr << endl;
  //   cerr << "File " << infilename << " is not open" << endl;
  //   cerr << endl;
  // }
  // TTree *treein = (TTree*)filetreein->Get("treemargtoys");
  // if(!treein){
  //   cerr << endl;
  //   cerr << "Tree object is empty" << endl;
  //   cerr << endl;
  // }
  //treein->SetBranchAddress("best_fit_params",&fparm_osc);
  //int Nentries = treein->GetEntries();
  //cout << "Nentries = " << Nentries << endl;
  // if(Nentries!=NtoysMarg){
  //   cerr << endl;
  //   cerr << "Nentries osc parameter tree = " << Nentries << endl;
  //   cerr << "NtoysMarg of marg binary files = " << NtoysMarg << endl;
  //   cerr << "They must be identical!" << endl;
  //   cerr << endl;
  //   exit(1);
  // }
  // double BFProf_Mix13=0.;
  // double BFProf_Mix23=0.;
  // double BFProf_Mix12=0.;
  // double BFProf_Dm23=0.;
  // double BFProf_Dm12=0.;
  // double BFProf_DeltaCP=0.;
  // double LikeMax=0.;
  // double BFMarg_Mix13=0.;
  // double BFMarg_Mix23=0.;
  // double BFMarg_Mix12=0.;
  // double BFMarg_Dm23=0.;
  // double BFMarg_Dm12=0.;
  // double BFMarg_DeltaCP=0.;
  // double *ToyMix13   = new double[NtoysMarg];
  // double *ToyMix23   = new double[NtoysMarg];
  // double *ToyMix12   = new double[NtoysMarg];
  // double *ToyDm23    = new double[NtoysMarg];
  // double *ToyDm12    = new double[NtoysMarg];
  // double *ToyDeltaCP = new double[NtoysMarg];
  // for(int ientry=0;ientry<Nentries;ientry++){
  //   treein->GetEntry(ientry);
  //   cout << endl;
  //   cout << "Toy " << ientry << endl;
  //   //cout << *fparm_osc << endl;
  //   int np = fparm_osc->GetNumberParameters();
  //   for(int ip = 0; ip < np; ip++) {
  //     conv::PhysParamId_t p = fparm_osc->WhichParam(ip);
  //     //cout << "Param: " << util::AsStringCompact(p) << " ";
  //     //if(util::IsSystematic(p)) cout << " --> Systematic "; // not oscillation parameters 
  //     //if(fparm_osc->IsFixed(p)) cout << " --> Fixed "; 
  //     //if(fparm_osc->IsNuisance(p)) cout << " --> Nuisance ";
  //     //cout << " --> " << fparm_osc->CurrValue(p);
  //     //cout << endl;
  //     if(util::IsSystematic(p)) continue; // not oscillation parameters 
  //     if(fparm_osc->IsFixed(p)) continue; 
  //     if(fparm_osc->IsNuisance(p)) continue;
  //     cout << "Param: " << util::AsStringCompact(p) << " = "
  // 	   << fparm_osc->CurrValue(p)
  // 	   << endl;
  //     if     (p==kPhypMixing13)  ToyMix13[ientry] = fparm_osc->CurrValue(p); 
  //     else if(p==kPhypMixing23)  ToyMix23[ientry] = fparm_osc->CurrValue(p); 
  //     else if(p==kPhypMixing12)  ToyMix12[ientry] = fparm_osc->CurrValue(p); 
  //     else if(p==kPhypDm23sq)    ToyDm23[ientry]  = fparm_osc->CurrValue(p); 
  //     else if(p==kPhypDm12sq)    ToyDm12[ientry]  = fparm_osc->CurrValue(p); 
  //     else if(p==kPhypCPV)       ToyDm12[ientry]  = fparm_osc->CurrValue(p); 
  //     else if(p==kPhypMixing13Bar)  ToyMix13[ientry] = fparm_osc->CurrValue(p); 
  //     else if(p==kPhypMixing23Bar)  ToyMix23[ientry] = fparm_osc->CurrValue(p); 
  //     else if(p==kPhypMixing12Bar)  ToyMix12[ientry] = fparm_osc->CurrValue(p); 
  //     else if(p==kPhypDm23sqBar)    ToyDm23[ientry]  = fparm_osc->CurrValue(p); 
  //     else if(p==kPhypDm12sqBar)    ToyDm12[ientry]  = fparm_osc->CurrValue(p); 
  //     else if(p==kPhypCPVBar)       ToyDm12[ientry]  = fparm_osc->CurrValue(p); 
  //     //if(util::IsSystematic(p)) continue; // not oscillation parameters 
  //   }
  // }
  // // Get best-fit oscillation parameters
  // for(int itoy=0;itoy<NtoysMarg;itoy++){
  //   double weight = LikeToy[itoy]/SumLike;
  //   if(!WeightMode) weight = 1.;
  //   // Take best-fit of marginalized osc params with marginalization
  //   BFMarg_Mix13   += ToyMix13[itoy]   * weight;
  //   BFMarg_Mix23   += ToyMix23[itoy]   * weight;
  //   BFMarg_Mix12   += ToyMix12[itoy]   * weight;
  //   BFMarg_Dm23    += ToyDm23[itoy]    * weight;
  //   BFMarg_Dm12    += ToyDm12[itoy]    * weight;
  //   BFMarg_DeltaCP += ToyDeltaCP[itoy] * weight;
  //   // Take best-fit of marginalized osc params with profiling (approximation)
  //   if(LikeMax<LikeToy[itoy]){ 
  //     BFProf_Mix13   = ToyMix13[itoy]; 
  //     BFProf_Mix23   = ToyMix23[itoy]; 
  //     BFProf_Mix12   = ToyMix12[itoy]; 
  //     BFProf_Dm23    = ToyDm23[itoy]; 
  //     BFProf_Dm12    = ToyDm12[itoy]; 
  //     BFProf_DeltaCP = ToyDeltaCP[itoy]; 
  //     LikeMax        = LikeToy[itoy];
  //   }
  //   // for(int i=0;i<util::NumErecoBins(s);i++){
  //   //   double NTotEvt = 0;
  //   //   #ifndef __USE_2D_TEMPLATES__
  //   //   NTotEvt = (double)fitmarg->GetToyMargHisto(s,itoy,ipar0,ipar1)->Integral();
  //   //   #else
  //   //   NTotEvt = (double)fitmarg->GetToyMargHisto2D(s,itoy,ipar0,ipar1)->Integral();
  //   //   #endif
  //   //   cout << util::AsStringCompact(s) 
  //   // 	   << " --> # of events = " << NTotEvt  
  //   // 	   << endl;
  //   // }
  //   // cout << " - Like = " << LikeToy[itoy] << endl;
  //   // cout << " - Mix13   = " << ToyMix13[itoy]   << endl;
  //   // cout << " - Mix23   = " << ToyMix23[itoy]   << endl;
  //   // cout << " - Mix12   = " << ToyMix12[itoy]   << endl;
  //   // cout << " - Dm23    = " << ToyDm23[itoy]    << endl;
  //   // cout << " - Dm12    = " << ToyDm12[itoy]    << endl;
  //   // cout << " - DeltaCP = " << ToyDeltaCP[itoy] << endl;
  //   // cout << endl;
  // }
  // // Draw best-fit oscillation parameters
  // TCanvas *c_bestfit_mix13 = new TCanvas("c_bestfit_mix13","c_bestfit_mix13");
  // fparm_osc->GetToyHisto(conv::kPhypMixing13)->DrawClone();
  // TLine *lMargMix13 = new TLine(BFMarg_Mix13,0,BFMarg_Mix13,fparm_osc->GetToyHisto(conv::kPhypMixing13)->GetMaximum());
  // lMargMix13->SetLineColor(kRed);
  // lMargMix13->SetLineWidth(2);
  // lMargMix13->Draw();
  // TLine *lProfMix13 = new TLine(BFProf_Mix13,0,BFProf_Mix13,fparm_osc->GetToyHisto(conv::kPhypMixing13)->GetMaximum());
  // lProfMix13->SetLineColor(kBlue);
  // lProfMix13->SetLineWidth(2);
  // lProfMix13->Draw();
  // TLegend *legoscpar = new TLegend(0.6, 0.4, 0.9, 0.9);  
  // legoscpar->SetFillColor(0);
  // legoscpar->SetBorderSize(0);
  // legoscpar->AddEntry(lMargMix13,"Marginalized","l");
  // legoscpar->AddEntry(lProfMix13,"Profiled","l");
  // legoscpar->Draw();
  // TCanvas *c_bestfit_mix23 = new TCanvas("c_bestfit_mix23","c_bestfit_mix23");
  // fparm_osc->GetToyHisto(conv::kPhypMixing23)->DrawClone();
  // TLine *lMargMix23 = new TLine(BFMarg_Mix23,0,BFMarg_Mix23,fparm_osc->GetToyHisto(conv::kPhypMixing23)->GetMaximum());
  // lMargMix23->SetLineColor(kRed);
  // lMargMix23->SetLineWidth(2);
  // lMargMix23->Draw();
  // TLine *lProfMix23 = new TLine(BFProf_Mix23,0,BFProf_Mix23,fparm_osc->GetToyHisto(conv::kPhypMixing23)->GetMaximum());
  // lProfMix23->SetLineColor(kBlue);
  // lProfMix23->SetLineWidth(2);
  // lProfMix23->Draw();
  // legoscpar->Draw();
  // TCanvas *c_bestfit_mix12 = new TCanvas("c_bestfit_mix12","c_bestfit_mix12");
  // fparm_osc->GetToyHisto(conv::kPhypMixing12)->DrawClone();
  // TLine *lMargMix12 = new TLine(BFMarg_Mix12,0,BFMarg_Mix12,fparm_osc->GetToyHisto(conv::kPhypMixing12)->GetMaximum());
  // lMargMix12->SetLineColor(kRed);
  // lMargMix12->SetLineWidth(2);
  // lMargMix12->Draw();
  // TLine *lProfMix12 = new TLine(BFProf_Mix12,0,BFProf_Mix12,fparm_osc->GetToyHisto(conv::kPhypMixing12)->GetMaximum());
  // lProfMix12->SetLineColor(kBlue);
  // lProfMix12->SetLineWidth(2);
  // lProfMix12->Draw();
  // legoscpar->Draw();
  // TCanvas *c_bestfit_dm23 = new TCanvas("c_bestfit_dm23","c_bestfit_dm23");
  // fparm_osc->GetToyHisto(conv::kPhypDm23sq)->DrawClone();
  // TLine *lMargDm23 = new TLine(BFMarg_Dm23,0,BFMarg_Dm23,fparm_osc->GetToyHisto(conv::kPhypDm23sq)->GetMaximum());
  // lMargDm23->SetLineColor(kRed);
  // lMargDm23->SetLineWidth(2);
  // lMargDm23->Draw();
  // TLine *lProfDm23 = new TLine(BFProf_Dm23,0,BFProf_Dm23,fparm_osc->GetToyHisto(conv::kPhypDm23sq)->GetMaximum());
  // lProfDm23->SetLineColor(kBlue);
  // lProfDm23->SetLineWidth(2);
  // lProfDm23->Draw();
  // legoscpar->Draw();
  // TCanvas *c_bestfit_dm12 = new TCanvas("c_bestfit_dm12","c_bestfit_dm12");
  // fparm_osc->GetToyHisto(conv::kPhypDm12sq)->DrawClone();
  // TLine *lMargDm12 = new TLine(BFMarg_Dm12,0,BFMarg_Dm12,fparm_osc->GetToyHisto(conv::kPhypDm12sq)->GetMaximum());
  // lMargDm12->SetLineColor(kRed);
  // lMargDm12->SetLineWidth(2);
  // lMargDm12->Draw();
  // TLine *lProfDm12 = new TLine(BFProf_Dm12,0,BFProf_Dm12,fparm_osc->GetToyHisto(conv::kPhypDm12sq)->GetMaximum());
  // lProfDm12->SetLineColor(kBlue);
  // lProfDm12->SetLineWidth(2);
  // lProfDm12->Draw();
  // legoscpar->Draw();
  // TCanvas *c_bestfit_dcp = new TCanvas("c_bestfit_dcp","c_bestfit_dcp");
  // fparm_osc->GetToyHisto(conv::kPhypCPV)->DrawClone();
  // TLine *lMargDcp = new TLine(BFMarg_DeltaCP,0,BFMarg_DeltaCP,fparm_osc->GetToyHisto(conv::kPhypCPV)->GetMaximum());
  // lMargDcp->SetLineColor(kRed);
  // lMargDcp->SetLineWidth(2);
  // lMargDcp->Draw();
  // TLine *lProfDcp = new TLine(BFProf_DeltaCP,0,BFProf_DeltaCP,fparm_osc->GetToyHisto(conv::kPhypCPV)->GetMaximum());
  // lProfDcp->SetLineColor(kBlue);
  // lProfDcp->SetLineWidth(2);
  // lProfDcp->Draw();
  // legoscpar->Draw();
  // // Print resutls
  // cout << endl;
  // cout << "Best-fit of marginalized oscillation parameters w/ marginalization" << endl;
  // cout << "(weighted average of osc params with likelihoods from all marginalization toys):" << endl;
  // cout << " - Mix13   = " << BFMarg_Mix13 << endl;
  // cout << " - Mix23   = " << BFMarg_Mix23 << endl;
  // cout << " - Mix12   = " << BFMarg_Mix12 << endl;
  // cout << " - Dm23    = " << BFMarg_Dm23 << endl;
  // cout << " - Dm12    = " << BFMarg_Dm12 << endl;
  // cout << " - DeltaCP = " << BFMarg_DeltaCP << endl;
  // cout << endl;
  // cout << "Best-fit of marginalized oscillation parameters w/ profiling" << endl;
  // cout << "(w/ highest likelihood from all marginalization toys):" << endl;
  // cout << " - Mix13   = " << BFProf_Mix13 << endl;
  // cout << " - Mix23   = " << BFProf_Mix23 << endl;
  // cout << " - Mix12   = " << BFProf_Mix12 << endl;
  // cout << " - Dm23    = " << BFProf_Dm23 << endl;
  // cout << " - Dm12    = " << BFProf_Dm12 << endl;
  // cout << " - DeltaCP = " << BFProf_DeltaCP << endl;
  // cout << endl;




  return;
}

TH1D * spectrum_to_per_50mev(TH1D * h) {
  TH1D * hout = new TH1D(*h);
  for(int ix = 1; ix <= hout->GetNbinsX(); ix++) {
    hout->SetBinContent(ix, 0.05 * h->GetBinContent(ix) / h->GetBinWidth(ix));
  }
  hout->GetYaxis()->SetTitle("Events per 0.050 GeV");
  return hout;
}

