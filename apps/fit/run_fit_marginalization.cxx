
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


// How to run the macro:


//
// TODO:
// - set correct oscillation parameters for reactor constraint
// 


//
// TODO:
// - Add string of freeparname0 and freeparname1
//




#include "SuperkObservation.h"
#include "SuperkExpectation.h"
#include "T2KPhysicalParameterLimits.h"

#include "Utilities.h"
#include "Binning.h"

#include "Build.h"
#include "UtlForMacrosT2XK.h"
#include "UtlFitParamListFactory.h"
#include "AnalysisOptions.h"
#include "FitMarginalization.h"

#include "Utilities/UtlForMacros.h"
#include "Utilities/OscParams.h"
#include "Utilities/RunOpt.h"
#include "Utilities/Env.h"
#include "Conventions/PhysHypothesisId.h"

#include "TROOT.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TStopwatch.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TMath.h"


#include <sstream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace valor;
using namespace valor::conv;

class MyRunOpt : public RunOpt
{
public:
  MyRunOpt(int argc, char ** argv, int parseopt) {
    Init();
    ReadFromCommandLine(argc, argv, parseopt);
    Print(cout);
  }
  ~MyRunOpt() {}
  void   Print (ostream & stream, const int w = 40) const {
    RunOpt::Print(stream, w);
  }
private:

  void Init() {
    RunOpt::Init();
  }
  void ReadFromCommandLine(int argc, char ** argv, int parseopt) {
    RunOpt::ReadFromCommandLine(argc, argv, parseopt);
  }
  void Usage(const int w = 40) {
    RunOpt::Usage(w);
  }
};



int main(int argc, char ** argv)
{
  TStopwatch sw;

  Env::Instance()->AutoSet("t2xk");
  cout << *Env::Instance() << endl;
  gROOT->SetBatch(true);
  MyRunOpt runopt(argc, argv, kParseCommon | kParseToys | kParseGlobal | kParseFitMarg | kParse3FOscParams);
  
  TString opt(runopt.Fit_param_list()); 

  // Analysis configuration
  TString config  (Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), runopt.ConfigFile().c_str() ));
  TString xmlfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), runopt.Xmlfile().c_str() ));
  TString rootfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), runopt.Rootfile().c_str() ));

  FileStat_t dummy;
  
  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(runopt.GlobalMassHierarchy());
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(runopt.GlobalDmsqScheme());
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(runopt.GlobalMixParamScheme());
  Opt::Instance()->gOscOctantChoice      = util::OctantFromString(runopt.GlobalMixParamOctant());
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }
  Opt::Instance()->SetDefaultParameterValues();

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(runopt.ConfigFile().c_str());

  // Fit param list for fakedata sets
  FitParamList * fparm  = util::NamedParameterList(runopt.Fit_param_list().c_str()); 
  
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used for fakedata sets: " << util::AsString(hyp) << endl;

  int first_experiment = runopt.FirstExperiment();
  int last_experiment = runopt.LastExperiment();
  int nexperiments = last_experiment - first_experiment + 1; 


  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << runopt.FileTag() << endl;
  cout << "- Reading analysis configuration from: " << config << endl;

  if(runopt.DoToys()){
    cout << "- Toy-MC generation" << endl;
    cout << "   - Number of toy experiments: " << nexperiments << endl; 
    cout << "   - Randomize systematics : "      << ((runopt.SystTweakMode()) ?      "YES " : "NO") << endl;
    cout << "   - Statistical fluctuations: "    << ((runopt.StatFluctuations()) ?      "YES " : "NO") << endl;
  }
  else{
    cout << "- Reading analysis configuration from: " << config << endl;
    cout << "- Toy-MC" << endl;
    if(runopt.Using_xml_file()) cout << "   - Reading experiments from: " << xmlfile2 << endl; 
    else cout << "   - Reading experiments from: " << rootfile2 << endl; 
    cout << "   - Using experiments in range: [" << first_experiment << ", " << last_experiment << "]" << endl; 
  }
  
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << runopt.Fit_param_list() << endl;
  cout << "- Outputs" << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;
  //cout << *Rnd::Instance() << endl;

  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;


  // Create file to store results

  string oscstr;// = util::AsString(hyp,input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);  

  //TFile *ftoys = new TFile(Form("Toys_%s-%s_%d_%d.root", runopt.FileTag().c_str(), oscstr.c_str(), first_experiment, last_experiment), "recreate");
  TFile *f = new TFile(Form("%s-%s_%d_%d.root", runopt.FileTag().c_str(), oscstr.c_str(), first_experiment, last_experiment), "recreate");  

  // Set Tree toy/fit
  int run=-999;
  string tag="";
  int nexpts=-999;
  int ToyID=-999;
  double truepar0=-999;
  double truepar1=-999;    
  int nfits=-999;
  int nchisqcalls=-999;
  int niterations=-999;
  int fitstatus=-999;
  double chisqmin=conv::kALargeNum;
  double chisqtrue0=conv::kALargeNum;
  double chisqtrue1=conv::kALargeNum;
  double deltachisqtrue=conv::kALargeNum;
  double deltachisq0BF=conv::kALargeNum;
  TH2F *hchisq2D_out = new TH2F("hchisq2D_out","hchisq2D_out",1,0,1,1,0,1);
  TH1F *hchisq1D_out = new TH1F("hchisq1D_out","hchisq1D_out",1,0,1);
  TH2F *hdeltachisq2D_out = new TH2F("hdeltachisq2D_out","hdeltachisq2D_out",1,0,1,1,0,1);
  TH1F *hdeltachisq1D_out = new TH1F("hdeltachisq1D_out","hdeltachisq1D_out",1,0,1);
  TH2F *hstep2D_out = new TH2F("hstep2D_out","hstep2D_out",1,0,1,1,0,1);
  TH1F *hstep1D_out = new TH1F("hstep1D_out","hstep1D_out",1,0,1);
  double bestfit0=-999;
  double bestfit1=-999;
  double bestfiterrNeg0=-999;
  double bestfiterrNeg1=-999;
  double bestfiterrPos0=-999;
  double bestfiterrPos1=-999;
  double bestfiterrSym0=-999;
  double bestfiterrSym1=-999;
  double bestfitFullChisq0=-999;
  double bestfitFullChisq1=-999;
  double chisqminFullChisq=conv::kALargeNum;
#ifndef __USE_2D_TEMPLATES__
  TH1F *hErecoToy_1Rmu = new TH1F("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1);
  TH1F *hErecoToy_1Re = new TH1F("hErecoToy_1Re","hErecoToy_1Re",1,0,1);
  TH1F *hErecoToy_RHC1Rmu = new TH1F("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1);
  TH1F *hErecoToy_RHC1Re = new TH1F("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1);
  TH1F *hErecoToy_sample[conv::kns] = {};
#else
  TH2F *hErecoToy_1Rmu = new TH2F("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1,1,0,1);
  TH2F *hErecoToy_1Re = new TH2F("hErecoToy_1Re","hErecoToy_1Re",1,0,1,1,0,1);
  TH2F *hErecoToy_RHC1Rmu = new TH2F("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1,1,0,1);
  TH2F *hErecoToy_RHC1Re = new TH2F("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1,1,0,1);
  TH2F *hErecoToy_sample[conv::kns] = {};
#endif
  hErecoToy_sample[conv::kSKsFCFV1RingMuLike   ] = hErecoToy_1Rmu;
  hErecoToy_sample[conv::kSKsFCFV1RingELike    ] = hErecoToy_1Re ;
  hErecoToy_sample[conv::kSKsRHCFCFV1RingMuLike] = hErecoToy_RHC1Rmu;
  hErecoToy_sample[conv::kSKsRHCFCFV1RingELike ] = hErecoToy_RHC1Re ;
 
  // Object for storing input oscillation and nuisance params in each toy-MC
  PhysParamList * iparm = new PhysParamList;

  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  fitnt->Branch("best_fit_params",      "valor::FitParamList",  &fparm,      1280000,0);
  fitnt->Branch("input_params",         "valor::PhysParamList", &iparm,      1280000,0);
  fitnt->Branch("truepar0", &truepar0 );
  fitnt->Branch("truepar1", &truepar1 );
  fitnt->Branch("bestfit0", &bestfit0 );
  fitnt->Branch("bestfit1", &bestfit1 );

  fitnt->Branch("hchisq1D",     "TH1F", &hchisq1D_out,     32000,0 );
  fitnt->Branch("hstep1D",      "TH1F", &hstep1D_out,      32000,0 );
  fitnt->Branch("hdeltachisq1D","TH1F", &hdeltachisq1D_out,32000,0 );
  fitnt->Branch("hchisq2D",     "TH2F", &hchisq2D_out,     32000,0 );
  fitnt->Branch("hdeltachisq2D","TH2F", &hdeltachisq2D_out,32000,0 );
  fitnt->Branch("hstep2D",      "TH2F", &hstep2D_out,      32000,0 );
  
  fitnt->Branch("fitstatus",&fitstatus,"fitstatus/I");
  fitnt->Branch("chisqmin",&chisqmin,"chisqmin/D");
  fitnt->Branch("chisqtrue0",&chisqtrue0,"chisqtrue0/D");
  fitnt->Branch("chisqtrue1",&chisqtrue1,"chisqtrue1/D");
  fitnt->Branch("deltachisqtrue",&deltachisqtrue,"deltachisqtrue/D");
  fitnt->Branch("nfits",&nfits,"nfits/I");
  fitnt->Branch("nchisqcalls",&nchisqcalls,"nchisqcalls/I");
  fitnt->Branch("niterations",&niterations,"niterations/I");
  fitnt->Branch("bestfiterrNeg0",&bestfiterrNeg0);
  fitnt->Branch("bestfiterrNeg1",&bestfiterrNeg1);
  fitnt->Branch("bestfiterrPos0",&bestfiterrPos0);
  fitnt->Branch("bestfiterrPos1",&bestfiterrPos1);
  fitnt->Branch("bestfiterrSym0",&bestfiterrSym0);
  fitnt->Branch("bestfiterrSym1",&bestfiterrSym1);
  fitnt->Branch("bestfitFullChisq0", &bestfitFullChisq0 );
  fitnt->Branch("bestfitFullChisq1", &bestfitFullChisq1 );
  fitnt->Branch("chisqminFullChisq",&chisqminFullChisq);
  fitnt->Branch("ToyID",&ToyID);
  fitnt->Branch("Job",&run);
  fitnt->Branch("Tag",&tag);
  fitnt->Branch("nexpts",&nexpts);
  #ifndef __USE_2D_TEMPLATES__
  fitnt->Branch("hErecoToy_1Rmu",     "TH1F",&hErecoToy_1Rmu,     32000,0);
  fitnt->Branch("hErecoToy_1Re",      "TH1F",&hErecoToy_1Re,      32000,0);
  fitnt->Branch("hErecoToy_RHC1Rmu",  "TH1F",&hErecoToy_RHC1Rmu,  32000,0);
  fitnt->Branch("hErecoToy_RHC1Re",   "TH1F",&hErecoToy_RHC1Re,   32000,0);
  #else  
  fitnt->Branch("hErecoToy_1Rmu",     "TH2F",&hErecoToy_1Rmu,     32000,0);
  fitnt->Branch("hErecoToy_1Re",      "TH2F",&hErecoToy_1Re,      32000,0);
  fitnt->Branch("hErecoToy_RHC1Rmu",  "TH2F",&hErecoToy_RHC1Rmu,  32000,0);
  fitnt->Branch("hErecoToy_RHC1Re",   "TH2F",&hErecoToy_RHC1Re,   32000,0);
  #endif
 
  // Update physics fitparamlist with chosen value of oscillation parameters
  bool set_init = false;
  set_init = true;
  util::OscParams::Instance()->UpdateFitParamList(fparm,set_init);
  set_init = false;
  util::OscParams::Instance()->UpdateFitParamList(fparm,set_init);

  int weird_applied_neg      = 0;
  int weird_unphysical_tweak = 0;
  last_experiment = runopt.LastExperiment();
  for(int iexperiment = runopt.FirstExperiment(); iexperiment <= last_experiment; iexperiment++) {
     cout << "\n\n*** Pseudo-experiment = " << iexperiment << endl;

     // Generate mock data
     cout << "- Generating toy-MC..." << endl;
     SuperkObservation * skobs_fakedata = new SuperkObservation;
     PhysHypothesis ph(runopt.Hypothesis());

     //util::OscParams::Instance()->Print(cout);
     
     // Update physics hypothesis and fitparamlist with chosen value of oscillation parameters
     util::OscParams::Instance()->UpdatePhysHypothesis(ph);
     
     if(runopt.DoToys()){
       
       // Randomize specified osc params for current toy MC expt
       util::OscParams::Instance()->Randomize();
       
       // Set the mass hierarchy
       util::SetMassHierarchy(skexp);
       
       // Update physics hypothesis and fitparamlist with chosen value of oscillation parameters
       util::OscParams::Instance()->UpdatePhysHypothesis(ph);
       set_init = false;
       util::OscParams::Instance()->UpdateFitParamList(fparm,set_init);
     
       // Fill toy histograms of osc params
       util::OscParams::Instance()->FillOscToyHisto(ph);
       
       // Reset input true params and add true oscillation params
       iparm->Clear();
       PhysParamList true_oscillation_params = fparm->AsPhysParamList("current");
       iparm->Append(true_oscillation_params);
       
       // Loop over toys
       bool applied_negevents_condition = true;
       bool unphysical_tweak = true;
       while(applied_negevents_condition || unphysical_tweak) {
         
         util::PerformSysThrow(runopt.SystTweakMode(), skexp, runopt.VerboseSystThrows());
         
         // Rethrow unphysical parameters
         PhysParamList true_nuisance_params = skexp->Tweaks();
         if (runopt.RethrowIfUnphysical()) {
           unphysical_tweak = !(valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(true_nuisance_params));
           if (unphysical_tweak) {
             cout << "Parameters are unphysical - rethrowing" << std::endl;
             weird_unphysical_tweak++;
             continue;
           }
         }
         else {
           unphysical_tweak = false;
         }
         
         applied_negevents_condition = skobs_fakedata->CreateMockData(*skexp, ph, runopt.StatFluctuations());
         cout << "applied_negevents_condition = " << applied_negevents_condition << endl;
         if(applied_negevents_condition)
           weird_applied_neg++;
       }//applied_negevents_condition
       
       cout << *fparm << endl;
       
     }
     
     else { // Read xml file
       
       if(runopt.Using_xml_file()){
         cout << "Use xml file" << endl;
         skobs_fakedata->Read(xmlfile2.Data(),iexperiment,runopt.Using_xml_file());
       }
       else{
         cout << "Use root file" << endl;
         skobs_fakedata->Read(rootfile2.Data(),iexperiment,runopt.Using_xml_file());
       }
     }
     
     /////////////////// SET MARGINALIZATION FIT ////////////////
  
     FitMarginalization *fitmarg = new FitMarginalization();   
     
     // Reset all the variables
     truepar0=-999;
     truepar1=-999;    
     nfits=-999;
     nchisqcalls=-999;
     niterations=-999;
     fitstatus=-999;
     chisqmin=conv::kALargeNum;
     chisqtrue0=conv::kALargeNum;
     chisqtrue1=conv::kALargeNum;
     deltachisqtrue=conv::kALargeNum;
     bestfit0=-999;
     bestfit1=-999;
     bestfiterrNeg0=-999;
     bestfiterrNeg1=-999;
     bestfiterrPos0=-999;
     bestfiterrPos1=-999;
     bestfiterrSym0=-999;
     bestfiterrSym1=-999;  
     bestfitFullChisq0=-999;
     bestfitFullChisq1=-999;
     chisqminFullChisq=conv::kALargeNum;
   
     // Perform the fit w/ marginalization
     util::PerformMargFit(argc,argv,fitmarg,skexp,skobs_fakedata); 
     
     nchisqcalls = fitmarg->GetNCallsCalcChisqMargPt();

     bestfiterrNeg0 = fitmarg->GetParamBestFitErrNeg(0);
     bestfiterrNeg1 = fitmarg->GetParamBestFitErrNeg(1);
     bestfiterrPos0 = fitmarg->GetParamBestFitErrPos(0);
     bestfiterrPos1 = fitmarg->GetParamBestFitErrPos(1);
     
     bestfiterrSym0 = fitmarg->GetParamBestFitErrSym(0);
     bestfiterrSym1 = fitmarg->GetParamBestFitErrSym(1);
     
     // Get true chisq 
     
     PhysParamId_t freepar0 = conv::kPhypUndefined;
     PhysParamId_t freepar1 = conv::kPhypUndefined;
    
     freepar0 = util::GetFreePar(runopt.Freeparname0());
     if(fitmarg->IsFit2D()) freepar1 = util::GetFreePar(runopt.Freeparname1());
     
     truepar0 = fparm->InitValue(freepar0);
     if(fitmarg->IsFit2D()) truepar1 = fparm->InitValue(freepar1);
     
     chisqtrue0 = (double)fitmarg->GetChisqMargTrueRead(truepar0,truepar1,skobs_fakedata,true);  // TODO
     
     cout << endl;
     cout << "Calculate Chisq for true point:" << endl;
     cout << "truepar0 = " << truepar0 << endl;
     cout << "truepar1 = " << truepar1 << endl;
     cout << endl;

     // Get Best Fit point
     bestfit0 = fitmarg->GetParamBestFit(0); 
     bestfit1 = fitmarg->GetParamBestFit(1);
     chisqmin = fitmarg->GetChisqBestFit();

     // Draw Full Chisq
     if(runopt.DrawFullChisq()){
       
       cout << endl;
       cout << "Draw the full Chisq" << endl;
       cout << endl;

       fitmarg->FillChisqMarg(skobs_fakedata,true);
       fitmarg->MinimizeHisto();
       bestfitFullChisq0 = fitmarg->GetParamBestFit(0);
       if(fitmarg->IsFit2D()) bestfitFullChisq1 = fitmarg->GetParamBestFit(1); 
       chisqminFullChisq = (double)fitmarg->GetChisqBestFit();
       fitmarg->BuildDeltaChisqMarg();
     }

          
     // Get delta chisq
     deltachisq0BF = chisqtrue0 - chisqmin;
     
     // Get the fakedata histogram
     for(int is = 0; is < kns; is++) {
       SuperkSample_t s = kSKSamples[is];
       if(s == kSKsUndefined) continue;
    
#ifndef __USE_2D_TEMPLATES__
       TH1D * hspec = skobs_fakedata->ErecoSpectrum(s);
#else
       TH2D * hspec = skobs_fakedata->TwoDSpectrum(s);
#endif
       hspec->Copy(*(hErecoToy_sample[s]));
       delete hspec;
     }
     
     // Fill the output tree
     
     if(fitmarg->IsFit2D()){
       fitmarg->GetChisqMarg2D     ()->Copy(*hchisq2D_out      );
       fitmarg->GetDeltaChisqMarg2D()->Copy(*hdeltachisq2D_out );
       fitmarg->GetChisqStep_2D    ()->Copy(*hstep2D_out       );
     
       int npt_chisq_x = hchisq2D_out->GetNbinsX();
       int npt_chisq_y = hchisq2D_out->GetNbinsY();
       for(int iy=0;iy<npt_chisq_y;iy++){
         for(int ix=0;ix<npt_chisq_x;ix++){
           cout << hchisq2D_out->GetBinContent(ix+1,iy+1) << " ";
         }
         cout << endl;
       }
       
     }
     else{
       fitmarg->GetChisqMarg1D     ()->Copy(*hchisq1D_out      );
       fitmarg->GetDeltaChisqMarg1D()->Copy(*hdeltachisq1D_out);
       fitmarg->GetChisqStep_1D    ()->Copy(*hstep1D_out);

       int npt_chisq_x = hchisq1D_out->GetNbinsX();
       for(int ix=0;ix<npt_chisq_x;ix++){
         cout << hchisq1D_out->GetBinContent(ix+1) << " ";
       }
       cout << endl;
       for(int ix=0;ix<npt_chisq_x;ix++){
         cout << hdeltachisq1D_out->GetBinContent(ix+1) << " ";
       }
       
     }
     
     // Fill toy informations (toy id, run, # of experiments, tag)
     ToyID = iexperiment;
     tag = runopt.FileTag();
     nexpts = runopt.NExperiments();
     run = runopt.Run();

     // Fill the tree 
     int iflag_tree = fitnt->Fill();
     cout << endl;
     cout << "iflag_tree = " << iflag_tree << endl;
     cout << endl;
     
     // Print output
     
     cout << endl;
     cout << "Minimum Marginalized Chisq = " << chisqmin << endl;        
     cout << "bestfit 0 = " << bestfit0;
     cout << " +" << bestfiterrPos0 
          << " -" << bestfiterrNeg0 
          << " (MINOS)"
          << endl;
     cout << "bestfit 0 = " << bestfit0
          << " +" << bestfiterrSym0 
          << " -" << bestfiterrSym0 
          << " (HESSE)"
          << endl;
     if(fitmarg->IsFit2D()){
       cout << "bestfit 1 = " << bestfit1; 
       cout << " +" << bestfiterrPos1
            << " -" << bestfiterrNeg1 
            << " (MINOS)"
            << endl;
       cout << "bestfit 1 = " << bestfit1
            << " +" << bestfiterrSym1
            << " -" << bestfiterrSym1 
            << " (HESSE)"
            << endl;
     }
     cout << endl;
     cout << "# of calls of CalcChisqMargPt(): " 
          << nchisqcalls << endl;
     cout << endl;
     cout << "Chisq True (" << truepar0 << "," << truepar1 << ") = " << chisqtrue0 << endl;
     cout << "Delta Chisq = Chisq_true - Chisq_bf = " << deltachisq0BF << endl;
     cout << endl;
     cout << "Minimum of full Chisq histogram:" << endl;
     cout << "  - chisqmin = " << chisqminFullChisq << endl; 
     cout << "  - bestfit0 = " << bestfitFullChisq0 << endl; 
     if(fitmarg->IsFit2D())
       cout << "  - bestfit1 = " << bestfitFullChisq1 << endl; 
     cout << endl;
     
     for(int is = 0; is < kns; is++) {
       SuperkSample_t s = kSKSamples[is];
       if(s == kSKsUndefined) continue;
       cout <<  "# of events in " << util::AsStringCompact(s) << " sample: "; 
       
#ifndef __USE_2D_TEMPLATES__
       cout << skobs_fakedata->ErecoSpectrum(s)->Integral() << endl;
#else
       cout << skobs_fakedata->TwoDSpectrum(s)->Integral() << endl;
#endif
       
     }
     cout << endl;
     
     
     //      // Store mock data in separated root file 
     //      // as when generate_toymc_ensemble.cxx is run
     //      ftoys->cd();
     // #ifndef __USE_2D_TEMPLATES__ //1D templates
     //      TH1D * mock_data[kns];
     // #else
     //      TH2D * mock_data[kns];
     // #endif
     //      for (int is = 0 ; is < conv::kns ; is++) {
     //        SuperkSample_t s = conv::kSKSamples[is];
     //        if (s == conv::kSKsUndefined) continue;
     // #ifndef __USE_2D_TEMPLATES__ 
     //        mock_data[is] = skobs_fakedata->ErecoSpectrum(s);
     // #else
     //        mock_data[is] = skobs_fakedata->TwoDSpectrum(s);
     // #endif
     //        mock_data[is]->SetDirectory(0);     
     //        mock_data[is]->Write(Form("mock_data_%s_%d", util::AsStringCompact(s).c_str(), iexperiment));
     //        delete mock_data[is];
     //}//is

     

     int npt_chisq_x = hchisq1D_out->GetNbinsX();
     for(int ix=0;ix<npt_chisq_x;ix++){
       cout << hchisq1D_out->GetBinContent(ix+1) << " ";
     }
     cout << endl;



     // // Write chisq
     // f->cd();
     // hchisq1D_out->Write();



     // Clean-up
     delete skobs_fakedata;
     fitmarg->CleanUp();
     delete fitmarg;


     hchisq1D_out->GetNbinsX();
     for(int ix=0;ix<npt_chisq_x;ix++){
       cout << hchisq1D_out->GetBinContent(ix+1) << " ";
     }
     cout << endl;

     
  } // loop over toys
  
  delete hErecoToy_1Rmu;
  delete hErecoToy_1Re;
  delete hErecoToy_RHC1Rmu;
  delete hErecoToy_RHC1Re;

  delete hchisq2D_out;
  delete hdeltachisq2D_out;
  delete hstep2D_out;
  delete hchisq1D_out;
  delete hdeltachisq1D_out;
  delete hstep1D_out;

  // Store ntuple with best-fit results
  f->cd();
  fitnt->Write(); // store the tree
  util::OscParams::Instance()->Write(f); // store the toy osc param histograms
  f->Close();
  delete f;
  
  sw.Stop();
  sw.Print();
}
