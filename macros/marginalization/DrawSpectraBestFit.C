
//
// ------------------------------------------------------------------------------------------------
//
// Check the extraction of the systematic and oscillation parameters done for marginalization
// as well as the Ereco distributions.
//
// Davide Sgalaberna
//
// ------------------------------------------------------------------------------------------------

#include "Build.h"
#include "AnalysisOptions.h"

using namespace std;
using namespace valor;
using namespace valor::conv;

void DrawSpectraBestFit
(
 const char *tag = "sinsqth13Vsdcp_T2Konly_Run14_DataNH",

 const char * xmlfile      = "data_run1to4/data_newstruct.xml",
 //const char * xmlfile      = "data_run1to6/data.xml",
 const char * rootfile     = "../outputs/test/toymc_prova_joint-0.root",
 const bool using_xml_file = true,
 int   iexperiment         = 0, // data (or fake data) event to show                    
 
 //const char * config_file  = "analysis_configuration_test_theta.xml",
 //const char * config_file  = "analysis_configuration_run1to6_joint_Nu6.914E20POT_NuBar4.0108E20POT.xml",
 const char * config_file  = "analysis_configuration_run1to4_joint_Nu6.57E20POT.xml",

 const char * sk_sample    = "1Re",   /* "1Rmu", "1Re", "RHC1Rmu", "RHC1Re", "MultiRmu" */

 bool SystMode = 1,  // 0 --> GetMode, 1 --> GetMean
 bool WeightsAreOne = false, // (for nominal or bestfit spectra) true --> All weights are 1, false --> Likelihoods are weights
 
 ////////////////////// MARGINALIZATION PARAMETERS ///////////////////

 //const char * likeinfilename = "/mnt/lnec/davides/VALOR/results/fit/liketoy_fitmarg_sinsqth13Vsdcp_T2Konly_Run14_Data_IH_FixMEC_toy10000_x_oscpar81_oscpar51_MHIH.bin", 
 //const char * likeinfilename = "/mnt/lnec/davides/VALOR/results/fit/liketoy_fitmarg_sinsqth13Vsdcp_T2Konly_Run14_Data_NH_StatOnly_toy1000_x_oscpar81_oscpar51_MHNH.bin",
 const char * likeinfilename = "/mnt/lnec/davides/VALOR/results/fit/liketoy_fitmarg_sinsqth13Vsdcp_T2Konly_Run14_Data_NH_toy10000_x_oscpar81_oscpar51_MHNH.bin", 

 const double InputPar0 = 0.048025,  // Input value of parameter 0 
 const double InputPar1 = 1.25664, // Input value of parameter 1

 const char * fit_param_list        = "nuosc3f+nuisance", // dcp

 //const char * infilename = "/mnt/lnec/davides/VALOR/results/joint_sinsqth13_Vs_dcp_t2konly_IH_10kToys_FixMEC/joint_sinsqth13_Vs_dcp_t2konly_constrall_IH_10kToys_FixMEC",
 //const char * infilename = "/mnt/lnec/davides/VALOR/results/joint_sinsqth13_Vs_dcp_t2konly_IH_10kToys_StatOnly/joint_sinsqth13_Vs_dcp_t2konly_constrall_IH_1kToys_StatOnly",
 const char * infilename = "/mnt/lnec/davides/VALOR/results/joint_sinsqth13_Vs_dcp_t2konlyrange_NH_10kToys_NEWTOYS/joint_dcp_t2konly_constrall_NH_10kToys",

 //const char * infilename = "../../inputs/marginalization/toys/JOINT/joint_dcp_react_NH_1kToys_StatOnly/joint_dcp_t2kreact_constrall_NH_1kToys_StatOnly",
 //const char * infilename = "../../inputs/marginalization/toys/JOINT/joint_dcp_react_IH_1kToys_StatOnly/joint_dcp_t2kreact_constrall_IH_1kToys_StatOnly",
 //const char * infilename = "../../inputs/marginalization/toys/JOINT/joint_dcp_react_NH_1kToys_OldJoint_StatOnly/joint_dcp_t2kreact_constrall_NH_1kToys_StatOnly",
 //const char * infilename = "../../inputs/marginalization/toys/JOINT/joint_dcp_react_IH_1kToys_OldJoint_StatOnly/joint_dcp_t2kreact_constrall_IH_1kToys_StatOnly",
 //const char * infilename = "../../inputs/marginalization/toys/JOINT/joint_dcp_react_NH_10kToys/joint_dcp_t2kreact_constrall_NH_10kToys",
 //const char * infilename = "../../inputs/marginalization/toys/JOINT/joint_dcp_react_IH_10kToys/joint_dcp_t2kreact_constrall_IH_10kToys",
 //const char * infilename = "../../inputs/marginalization/toys/JOINT/joint_dcp_react_NH_10kToys_OldJoint/joint_dcp_t2kreact_constrall_NH_10kToys_OldJoint",
 //const char * infilename = "../../inputs/marginalization/toys/JOINT/joint_dcp_react_IH_10kToys_OldJoint/joint_dcp_t2kreact_constrall_IH_10kToys_OldJoint",
 //const char * infilename = "submit_prova",
 //const char * infilename = "PROVA_SystAll",
 
 //const int    nfiles               = 201, //201, // 1         // # of files to read
 const int    nfiles               = 546, //201, // 1         // # of files to read

 const int    NtoysMarg            = 10000,   // # toys marginalization 
 
 const char * freeparname0        = "mixing13",  // free parameter
 const double freeparmin0         = 0., // minimum value of free parameter
 const double freeparmax0         = 0.113, // maxinum value of free parameter
 const int    NptScan0            = 81, //81, //201, //21,      // # pts of free parameter to scan
 const char * freeparname1        = "deltacp",  // free parameter
 const double freeparmin1         = -3.14159265358979312e+00,         // minimum value of free parameter
 const double freeparmax1         = +3.14159265358979312e+00, // 10., 1         // maxinum value of free parameter
 const int    NptScan1            = 51, //51, //201, //3, //21, // 801, 2          // # pts of free parameter to scan 
 
 // const char * freeparname0        = "deltacp",  // free parameter
 // const double freeparmin0         = -3.14159265358979312e+00,         // minimum value of free parameter
 // const double freeparmax0         = +3.14159265358979312e+00, // 10., 1         // maxinum value of free parameter
 // const int    NptScan0            = 201, //201, // 801, 2          // # pts of free parameter to scan
 // const char * freeparname1        = "",  // free parameter
 // const double freeparmin1         = 0., // minimum value of free parameter
 // const double freeparmax1         = 0., // maxinum value of free parameter
 // const int    NptScan1            = 1, //201, //21,      // # pts of free parameter to scan

 const bool  doMargMH             = false,    // marginalize Mass Hierarchy
 const char * hierarchy_marg      = "normal"     /* normal or inverted */
 )
{
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

  // Get Ereco toys
  
  conv::SuperkSample_t s = util::SKSampleFromString(sk_sample); // requested sample
  if(s == conv::kSKsUndefined) {
    cerr << "** Undefined sample!" << endl;
    return;
  }
  
  // Declare the histograms
  const int n_ereco_bins = util::NumErecoBins(s);  
  TH1F* h_marg_toy;
  TH1D *h_bestfit_marg = new TH1D("h_bestfit_marg","h_bestfit_marg",util::NumErecoBins(s),util::ErecoBinEdges(s));
  
  //
  
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s); // index of specified sample
  
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
    h_marg_toy = (TH1F*)fitmarg->GetToyMargHisto2D(s,itoy,ipar0,ipar1)->ProjectionY();
#endif
    
    for(int i=1;i<=util::NumErecoBins(s);i++) NentriesErecoToy[i-1][itoy] = h_marg_toy->GetBinContent(i);
    
    if(fitmarg->IsFit2D())
      LikeToy[itoy] = fitmarg->GetLikeToyRead(InputPar0,truepar1,itoy,true); 
    else
      LikeToy[itoy] = fitmarg->GetLikeToyRead(InputPar0,InputPar1,itoy,true); 
    //cout << "toy " << itoy << " - like = " << like << endl;

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
    cout << "bin # " << ibin << " --> lowedge = " << lowedge[ibin] << " - upedge = " << upedge[ibin] << endl;
    
  } // loop over Ereco bins
  

  
  // Fill bin histograms
  
  TH1D **h_bin = new TH1D*[util::NumErecoBins(s)];
  char name[200];
  for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){
    sprintf(name,"h_bin_%d",ibin);
    h_bin[ibin] = new TH1D(name,name,50,lowedge[ibin],upedge[ibin]);
  }
  
  for(int i=0;i<util::NumErecoBins(s);i++){
    for(int itoy=0;itoy<NtoysMarg;itoy++){
      double weight = LikeToy[itoy]/SumLike;
      if(WeightsAreOne) weight = 1.;
      h_bin[i]->Fill(NentriesErecoToy[i][itoy],weight);
    }
  }
    
  // Find the mode/mean of each Ereco bin
  
  double *bin_mode = new double[n_ereco_bins];
  
  for(int iereco=0;iereco<util::NumErecoBins(s);iereco++){
    
    double y = h_bin[iereco]->GetMean();
    double ey = h_bin[iereco]->GetRMS();
    
    if(SystMode==0){ // Get the mode
      cout << "bin # " << iereco << " --> Integral = " << h_bin[iereco]->Integral()
	   << " - Mean = " << y << " - RMS = " << ey << endl;
      TF1 *f1 = new TF1("f1","gaus", y-ey,y+ey);  
      //TF1 *f1 = new TF1("f1","pol2", y-ey,y+ey);  
      //TF1 *f1 = new TF1("f1","landau",y-ey,y+ey);  
      //TF1 *f1 = new TF1("f1","gaus+landau",h_bin[iereco]->GetXaxis()->GetBinLowEdge(1),h_bin[iereco]->GetXaxis()->GetBinUpEdge(h_bin[iereco]->GetNbinsX()));  
      //TF1 *f1 = new TF1("f1","gaus",h_bin[iereco]->GetXaxis()->GetBinLowEdge(1),h_bin[iereco]->GetXaxis()->GetBinUpEdge(h_bin[iereco]->GetNbinsX()));        
      TFitResultPtr r = h_bin[iereco]->Fit("f1","S","goff",y-ey,y+ey);       
      bin_mode[iereco] = r->Parameter(1);
      cout << bin_mode[iereco] << " - ";
      bin_mode[iereco] = f1->GetMaximumX();
      cout << bin_mode[iereco] << endl;
    } 
    else if(SystMode==1){ // Get the mean
      bin_mode[iereco] = y;
    }
    
    if(bin_mode[iereco]<0.) bin_mode[iereco] = 0.;
    
    h_bestfit_marg->SetBinContent(iereco+1,bin_mode[iereco]);
  
    cout << "# of events in Ereco " << iereco << " = " << bin_mode[iereco] << endl;
  }
  
  // Draw bin histograms
  char name[200];
  TCanvas **c_bin = new TCanvas*[util::NumErecoBins(s)];
  for(int iereco=0;iereco<util::NumErecoBins(s);iereco++){
    sprintf(name,"c_bin_%d",iereco);    
    c_bin[iereco] = new TCanvas(name,name);
    h_bin[iereco]->Draw();
  }
  
  TH1D *hErecoToy;
#ifndef __USE_2D_TEMPLATES__
  hErecoToy = (TH1D*)skobs_fakedata->ErecoSpectrum(s)->Clone("hErecoToy");
#else
  hErecoToy = (TH1D*)skobs_fakedata->TwoDSpectrum(s)->ProjectionX()->Clone("hErecoToy");
#endif  
  
  // Draw best-fit histogram marginalizing the systematics
  TCanvas *c_bestfit_marg = new TCanvas("c_bestfit_marg","c_bestfit_marg");
  hErecoToy->Draw("pe");
  h_bestfit_marg->Draw("same");

  // Store best-fit spectrum in a file
  char name[200];
  sprintf(name,"spectrum_bestfit_%s.root",tag);
  TFile *fout = new TFile(name,"RECREATE");
  hErecoToy->Write();
  h_bestfit_marg->Write();
  fout->Close();

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

