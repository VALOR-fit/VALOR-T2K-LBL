
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
// IMPORTANT:
//
// An issue is the mass hierarchy used in the marginalization FitMarginalization uses gMassHierarchyAssumed (the global one). 
// This variable is the same that is set for the toys to generate the fake datasets and there could be a conflict.           
// For this reasone util::SetMassHierarchy(skexp) is used to set MH for fake datasets at each toy iteration.
// Furthermore when making fits of NH and IH, gMassHierarchyAssumed changes its value.
// It's important to draw the full chisq distribution using FitMarginalization::FillChisqMarg() just after the 
// minimization has been performed, otherwise gMassHierarchyAssumed changes (NH->IH) and the marginalization toy files get wrong path and it gives error.
// The fit must be done with the following order:
// - Build fake dataset with given gMassHierarchyAssumed
// - Minimization NH (gMassHierarchyAssumed=normal)
// - Draw Full Chisq NH (gMassHierarchyAssumed=normal)
// - Minimization IH (gMassHierarchyAssumed=inverted)
// - Draw Full Chisq IH (gMassHierarchyAssumed=inverted)
// At the next event iteration gMassHierarchyAssumed is set again to the defined one for the fake datasets.
// 
// NEED TO BE UPDATED THE WAY WE PASS MASS HIERARCHY TO FitMarginalization in order to avoid this conflict 
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
  TStopwatch sw_all;

  Env::Instance()->AutoSet("t2xk");
  cout << *Env::Instance() << endl;
  gROOT->SetBatch(true);
  //MyRunOpt runopt(argc, argv, kParseCommon | kParseToys | kParseGlobal | kParseFitMarg | kParse3FOscParams);
  MyRunOpt *runopt = new MyRunOpt(argc, argv, kParseCommon | kParseToys | kParseGlobal | kParseFitMarg | kParse3FOscParams);
  
  TString opt(runopt->Fit_param_list()); 

  // Analysis configuration
  TString config  (Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), runopt->ConfigFile().c_str() ));
  TString xmlfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), runopt->Xmlfile().c_str() ));
  TString rootfile2(Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), runopt->Rootfile().c_str() ));

  FileStat_t dummy;
  
  // Global options
  Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString(runopt->GlobalMassHierarchy());
  Opt::Instance()->gOscDmsqScheme        = util::DmsqSchemeFromString(runopt->GlobalDmsqScheme());
  Opt::Instance()->gOscMixingParamChoice = util::MixingTypeFromString(runopt->GlobalMixParamScheme());
  Opt::Instance()->gOscOctantChoice      = util::OctantFromString(runopt->GlobalMixParamOctant());
  if(Opt::Instance()->Undefined()) {
    cerr << "Global options not defined correctly. Exiting..." << endl;
    cerr << *Opt::Instance() << endl;
    exit(1);
  }
  Opt::Instance()->SetDefaultParameterValues();

  // Get superk spectrum predictor & configure
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(runopt->ConfigFile().c_str());

  // Fit param list for fakedata sets
  FitParamList * fparm  = util::NamedParameterList(runopt->Fit_param_list().c_str()); 
  
  PhysHypothesis_t hyp = fparm->PhysicsHypothesis();
  cout << "- The following physics hypothesis is used for fakedata sets: " << util::AsString(hyp) << endl;

  int first_experiment = runopt->FirstExperiment();
  int last_experiment = runopt->LastExperiment();
  int nexperiments = last_experiment - first_experiment + 1; 


  // Print settings
  cout << "\n[Job inputs]" << endl;
  cout << "- Job: " << runopt->FileTag() << endl;
  cout << "- Reading analysis configuration from: " << config << endl;

  if(runopt->DoToys()){
    cout << "- Toy-MC generation" << endl;
    cout << "   - Number of toy experiments: " << nexperiments << endl; 
    cout << "   - Randomize systematics : "      << ((runopt->SystTweakMode()) ?      "YES " : "NO") << endl;
    cout << "   - Statistical fluctuations: "    << ((runopt->StatFluctuations()) ?      "YES " : "NO") << endl;
  }
  else{
    cout << "- Reading analysis configuration from: " << config << endl;
    cout << "- Toy-MC" << endl;
    if(runopt->Using_xml_file()) cout << "   - Reading experiments from: " << xmlfile2 << endl; 
    else cout << "   - Reading experiments from: " << rootfile2 << endl; 
    cout << "   - Using experiments in range: [" << first_experiment << ", " << last_experiment << "]" << endl; 
  }
  
  cout << "- Fit" << endl;
  cout << "   - Fit Parameter List: " << runopt->Fit_param_list() << endl;
  cout << "- Outputs" << endl;
  cout << endl;
  cout << *Opt::Instance() << endl;
  //cout << *Rnd::Instance() << endl;

  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(hyp) << endl;


  // Create file to store results

  string oscstr;// = util::AsString(hyp,input_mixing12,input_mixing13,input_mixing23,input_deltacp,input_dm12sq,input_dm23sq);  

  TFile *ftoys = new TFile(Form("Toys_%s-%s_%d_%d.root", runopt->FileTag().c_str(), oscstr.c_str(), first_experiment, last_experiment), "recreate");
  TFile *f = new TFile(Form("%s-%s_%d_%d.root", runopt->FileTag().c_str(), oscstr.c_str(), first_experiment, last_experiment), "recreate");  

  // Set Tree toy/fit
  
  double CPUtime_all=conv::kBadResult;
  double REALtime_all=conv::kBadResult;
  double CPUtime_Minim_nh=conv::kBadResult;
  double CPUtime_Minim_ih=conv::kBadResult;
  double CPUtime_FullChisq_nh=conv::kBadResult;
  double CPUtime_FullChisq_ih=conv::kBadResult;
  double REALtime_Minim_nh=conv::kBadResult;
  double REALtime_Minim_ih=conv::kBadResult;
  double REALtime_FullChisq_nh=conv::kBadResult;
  double REALtime_FullChisq_ih=conv::kBadResult;
  double CPUtime_end=conv::kBadResult;
  double REALtime_end=conv::kBadResult;

  bool IsFit2D_nh=false;
  bool IsFit2D_ih=false;
  int run=-999;
  string tag="";
  int nexpts=-999;
  int ToyID=-999;
  
  // true variables
  double truepar0=-999;
  double truepar1=-999;    
  string trueMH = "";
  string trueparname0 = "";
  string trueparname1 = "";
    
  // global fit output
  string MH_glob = "";
  double deltachisq0BF_glob=conv::kALargeNum;  
  double deltachisq01_glob=conv::kALargeNum;  
  double chisqmin_glob=conv::kALargeNum;
  double chisqmin_nh=conv::kALargeNum;
  double chisqmin_ih=conv::kALargeNum;
  
#ifndef __USE_2D_TEMPLATES__
  TH1D *hErecoToy_1Rmu = new TH1D("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1);
  TH1D *hErecoToy_1Re = new TH1D("hErecoToy_1Re","hErecoToy_1Re",1,0,1);
  TH1D *hErecoToy_RHC1Rmu = new TH1D("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1);
  TH1D *hErecoToy_RHC1Re = new TH1D("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1);
  TH1D *hErecoToy_sample[conv::kns] = {};
  TH1D *hErecoExp_1Rmu = new TH1D("hErecoExp_1Rmu","hErecoExp_1Rmu",1,0,1);
  TH1D *hErecoExp_1Re = new TH1D("hErecoExp_1Re","hErecoExp_1Re",1,0,1);
  TH1D *hErecoExp_RHC1Rmu = new TH1D("hErecoExp_RHC1Rmu","hErecoExp_RHC1Rmu",1,0,1);
  TH1D *hErecoExp_RHC1Re = new TH1D("hErecoExp_RHC1Re","hErecoExp_RHC1Re",1,0,1);
  TH1D *hErecoExp_sample[conv::kns] = {};
#else
  TH2D *hErecoToy_1Rmu = new TH2D("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1,1,0,1);
  TH2D *hErecoToy_1Re = new TH2D("hErecoToy_1Re","hErecoToy_1Re",1,0,1,1,0,1);
  TH2D *hErecoToy_RHC1Rmu = new TH2D("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1,1,0,1);
  TH2D *hErecoToy_RHC1Re = new TH2D("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1,1,0,1);
  TH2D *hErecoToy_sample[conv::kns] = {};
  TH2D *hErecoExp_1Rmu = new TH2D("hErecoExp2D_1Rmu","hErecoExp2D_1Rmu",1,0,1,1,0,1);
  TH2D *hErecoExp_1Re = new TH2D("hErecoExp2D_1Re","hErecoExp2D_1Re",1,0,1,1,0,1);
  TH2D *hErecoExp_RHC1Rmu = new TH2D("hErecoExp2D_RHC1Rmu","hErecoExp2D_RHC1Rmu",1,0,1,1,0,1);
  TH2D *hErecoExp_RHC1Re = new TH2D("hErecoExp2D_RHC1Re","hErecoExp2D_RHC1Re",1,0,1,1,0,1); 
  TH2D *hErecoExp_sample[conv::kns] = {};
#endif

  hErecoToy_sample[conv::kSKsFCFV1RingMuLike   ] = hErecoToy_1Rmu;
  hErecoToy_sample[conv::kSKsFCFV1RingELike    ] = hErecoToy_1Re ;
  hErecoToy_sample[conv::kSKsRHCFCFV1RingMuLike] = hErecoToy_RHC1Rmu;
  hErecoToy_sample[conv::kSKsRHCFCFV1RingELike ] = hErecoToy_RHC1Re ;
  hErecoExp_sample[conv::kSKsFCFV1RingMuLike   ] = hErecoExp_1Rmu;
  hErecoExp_sample[conv::kSKsFCFV1RingELike    ] = hErecoExp_1Re ;
  hErecoExp_sample[conv::kSKsRHCFCFV1RingMuLike] = hErecoExp_RHC1Rmu;
  hErecoExp_sample[conv::kSKsRHCFCFV1RingELike ] = hErecoExp_RHC1Re ;
 
  // Object for storing input oscillation and nuisance params in each toy-MC
  PhysParamList * iparm = new PhysParamList;

  // Create tree to store outputs
  TTree * fitnt = new TTree("fitnt","Results of VaLOR data and/or toy-MC experiment fits");
  fitnt->Branch("best_fit_params",      "valor::FitParamList",  &fparm,      1280000,0);
  fitnt->Branch("input_params",         "valor::PhysParamList", &iparm,      1280000,0);

  fitnt->Branch("CPUtime_Minim_nh",&CPUtime_Minim_nh);
  fitnt->Branch("CPUtime_Minim_ih",&CPUtime_Minim_ih);
  fitnt->Branch("CPUtime_FullChisq_nh",&CPUtime_FullChisq_nh);
  fitnt->Branch("CPUtime_FullChisq_ih",&CPUtime_FullChisq_ih);
  fitnt->Branch("REALtime_Minim_nh",&REALtime_Minim_nh);
  fitnt->Branch("REALtime_Minim_ih",&REALtime_Minim_ih);
  fitnt->Branch("REALtime_FullChisq_nh",&REALtime_FullChisq_nh);
  fitnt->Branch("REALtime_FullChisq_ih",&REALtime_FullChisq_ih);
  fitnt->Branch("CPUtime_end",&CPUtime_end);
  fitnt->Branch("REALtime_end",&REALtime_end);

  fitnt->Branch("ToyID",&ToyID);
  fitnt->Branch("Job",&run);
  fitnt->Branch("Tag",&tag);
  fitnt->Branch("nexpts",&nexpts);

  // true variables
  fitnt->Branch("truepar0", &truepar0 );
  fitnt->Branch("truepar1", &truepar1 );
  fitnt->Branch("trueMH",   &trueMH   );
  fitnt->Branch("trueparname0",&trueparname0);
  fitnt->Branch("trueparname1",&trueparname1);
  
  // global fit output
  fitnt->Branch("deltachisq0BF_glob",&deltachisq0BF_glob,"deltachisq0BF_glob/D");
  fitnt->Branch("deltachisq01_glob",&deltachisq01_glob,"deltachisq01_glob/D");
  fitnt->Branch("MH_glob",   &MH_glob   );
  fitnt->Branch("chisqmin_glob",&chisqmin_glob,"chisqmin_glob/D");
  fitnt->Branch("chisqmin_nh",&chisqmin_nh,"chisqmin_nh/D");
  fitnt->Branch("chisqmin_ih",&chisqmin_ih,"chisqmin_ih/D");
    
  // histograms 
#ifndef __USE_2D_TEMPLATES__
  fitnt->Branch("hErecoToy_1Rmu",     "TH1D",&hErecoToy_1Rmu,     32000,0);
  fitnt->Branch("hErecoToy_1Re",      "TH1D",&hErecoToy_1Re,      32000,0);
  fitnt->Branch("hErecoToy_RHC1Rmu",  "TH1D",&hErecoToy_RHC1Rmu,  32000,0);
  fitnt->Branch("hErecoToy_RHC1Re",   "TH1D",&hErecoToy_RHC1Re,   32000,0);
  fitnt->Branch("hErecoExp_1Rmu",     "TH1D",&hErecoExp_1Rmu,     32000,0);
  fitnt->Branch("hErecoExp_1Re",      "TH1D",&hErecoExp_1Re,      32000,0);
  fitnt->Branch("hErecoExp_RHC1Rmu",  "TH1D",&hErecoExp_RHC1Rmu,  32000,0);
  fitnt->Branch("hErecoExp_RHC1Re",   "TH1D",&hErecoExp_RHC1Re,   32000,0);
#else  
  fitnt->Branch("hErecoToy_1Rmu",     "TH2D",&hErecoToy_1Rmu,     32000,0);
  fitnt->Branch("hErecoToy_1Re",      "TH2D",&hErecoToy_1Re,      32000,0);
  fitnt->Branch("hErecoToy_RHC1Rmu",  "TH2D",&hErecoToy_RHC1Rmu,  32000,0);
  fitnt->Branch("hErecoToy_RHC1Re",   "TH2D",&hErecoToy_RHC1Re,   32000,0);
  fitnt->Branch("hErecoExp_1Rmu",     "TH2D",&hErecoExp_1Rmu,     32000,0);
  fitnt->Branch("hErecoExp_1Re",      "TH2D",&hErecoExp_1Re,      32000,0);
  fitnt->Branch("hErecoExp_RHC1Rmu",  "TH2D",&hErecoExp_RHC1Rmu,  32000,0);
  fitnt->Branch("hErecoExp_RHC1Re",   "TH2D",&hErecoExp_RHC1Re,   32000,0);
#endif
 
  // Update physics fitparamlist with chosen value of oscillation parameters
  bool set_init = false;
  set_init = true;
  util::OscParams::Instance()->UpdateFitParamList(fparm,set_init);
  set_init = false;
  util::OscParams::Instance()->UpdateFitParamList(fparm,set_init);

  int weird_applied_neg      = 0;
  int weird_unphysical_tweak = 0;
  last_experiment = runopt->LastExperiment();
  for(int iexperiment = runopt->FirstExperiment(); iexperiment <= last_experiment; iexperiment++) {
    
    TStopwatch sw_fit;

    cout << endl;
    cout << "CPU time" << endl;
    cout << endl;
    sw_fit.Stop();
    sw_fit.Print();
    sw_fit.Start(0);
    
    cout << "\n\n*** Pseudo-experiment = " << iexperiment << endl;
    
    // Generate mock data
     cout << "- Generating toy-MC..." << endl;
     SuperkObservation * skobs_fakedata = new SuperkObservation;
     PhysHypothesis ph(runopt->Hypothesis());

     util::OscParams::Instance()->Print(cout);
     
     // Update physics hypothesis and fitparamlist with chosen value of oscillation parameters
     util::OscParams::Instance()->UpdatePhysHypothesis(ph);
     
     if(runopt->DoToys()){
       
       // Randomize specified osc params for current toy MC expt
       util::OscParams::Instance()->Randomize();
       
       // Set the mass hierarchy (to avoid conflicts with fitmarginalization mass hierarchy)
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
         
         util::PerformSysThrow(runopt->SystTweakMode(), skexp, runopt->VerboseSystThrows());
         
         // Rethrow unphysical parameters
         PhysParamList true_nuisance_params = skexp->Tweaks();
         if (runopt->RethrowIfUnphysical()) {
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
         
         applied_negevents_condition = skobs_fakedata->CreateMockData(*skexp, ph, runopt->StatFluctuations());
         cout << "applied_negevents_condition = " << applied_negevents_condition << endl;
         if(applied_negevents_condition)
           weird_applied_neg++;
       }//applied_negevents_condition
       
       //cout << *fparm << endl;
       
     }
     
     else { // Read xml file
       
       if(runopt->Using_xml_file()){
         cout << "Use xml file" << endl;
         skobs_fakedata->Read(xmlfile2.Data(),iexperiment,runopt->Using_xml_file());
       }
       else{
         cout << "Use root file" << endl;
         skobs_fakedata->Read(rootfile2.Data(),iexperiment,runopt->Using_xml_file());
       }
     }

     
     /////////////////// SET MARGINALIZATION FIT ////////////////
  
     // Reset all the variables

     truepar0=-999;
     truepar1=-999;  
     trueMH="";
     trueparname0="";
     trueparname1="";
     deltachisq0BF_glob=conv::kALargeNum;  
     deltachisq01_glob=conv::kALargeNum;  
     MH_glob = "";
     chisqmin_glob=conv::kALargeNum;
     chisqmin_nh=conv::kALargeNum;
     chisqmin_ih=conv::kALargeNum;

     //
     
     trueparname0 = runopt->Freeparname0();
     trueparname1 = runopt->Freeparname1();     
     if( trueparname0 != "betaBar" ){
       cerr << endl;
       cerr << "The true parameter name is " << trueparname0 << endl;
       cerr << "You can fit the Mass Hierarchy only using betaBar as free parameter" << endl;
       cerr << endl;
       exit(1);
     }

     PhysParamId_t freepar0 = conv::kPhypUndefined;
     PhysParamId_t freepar1 = conv::kPhypUndefined;
     
     freepar0 = util::GetFreePar(runopt->Freeparname0());     
     truepar0 = fparm->InitValue(freepar0);

     //freepar1 = util::GetFreePar(runopt->Freeparname1());     
     //truepar1 = fparm->InitValue(freepar1);
        
     FitMarginalization *fitmarg_nh = new FitMarginalization();   
     FitMarginalization *fitmarg_ih = new FitMarginalization();   
      
     //
     // Perform the fit w/ marginalization
     //

     //
     // IMPORTANT:
     //
     // An issue is the mass hierarchy used in the marginalization FitMarginalization uses gMassHierarchyAssumed (the global one). 
     // This variable is the same that is set for the toys to generate the fake datasets and there could be a conflict.           
     // For this reasone util::SetMassHierarchy(skexp) is used to set MH for fake datasets at each toy iteration.
     // Furthermore when making fits of NH and IH, gMassHierarchyAssumed changes its value.
     // It's important to draw the full chisq distribution using FitMarginalization::FillChisqMarg() just after the 
     // minimization has been performed, otherwise gMassHierarchyAssumed changes (NH->IH) and the marginalization toy files get wrong path and it gives error.
     // The fit must be done with the following order:
     // - Build fake dataset with given gMassHierarchyAssumed
     // - Minimization NH (gMassHierarchyAssumed=normal)
     // - Draw Full Chisq NH (gMassHierarchyAssumed=normal)
     // - Minimization IH (gMassHierarchyAssumed=inverted)
     // - Draw Full Chisq IH (gMassHierarchyAssumed=inverted)
     // At the next event iteration gMassHierarchyAssumed is set again to the defined one for the fake datasets.
     // 

     if(truepar0 != 1.){
       
       cerr << endl;
       cerr << "truepar0 = " << truepar0 << endl;
       cerr << "truepar1 = " << truepar1 << endl;
       cerr << endl;
       cerr << "////////////////////////////////////////////////////////////" << endl;
       cerr << "//                                                        //" << endl;
       cerr << "// WARNING: the true parameter (betaBar) must be set to 1 //" << endl;
       cerr << "//          for Mass Hierarchy fit                        //" << endl;
       cerr << "//                                                        //" << endl;
       cerr << "//          Force it to be 1!!!                           //" << endl;
       cerr << "//                                                        //" << endl;
       cerr << "////////////////////////////////////////////////////////////" << endl;
       cerr << endl;
       
       truepar0 = 1.; // only for betabar = 1
     }

     cout << endl;
     cout << "//  Start fit of Normal Hierarchy" << endl;
     cout << endl;

     util::PerformMargFit(argc,argv,runopt,fitmarg_nh,skexp,skobs_fakedata,false,"normal"); // Normal Hierarchy marginalization toys
     Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal"); // needed to build the right binary file name 
     chisqmin_nh = (double)fitmarg_nh->GetChisqMargTrueRead(truepar0,truepar1,skobs_fakedata,true);  
     
     cout << endl;
     cout << "Fit of Normal Hierarchy is done!" << endl;
     cout << endl;
     
     /////

     cout << endl;
     cout << "//  Start fit of Inverted Hierarchy" << endl;
     cout << endl;
     
     util::PerformMargFit(argc,argv,runopt,fitmarg_ih,skexp,skobs_fakedata,false,"inverted"); // Inverted Hierarchy marginalization toys
     Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted"); // needed to build the right binary file name 		       
     chisqmin_ih = (double)fitmarg_ih->GetChisqMargTrueRead(truepar0,truepar1,skobs_fakedata,true);  
     
     cout << endl;
     cout << "Fit of Inverted Hierarchy is done!" << endl;
     cout << endl;

  
     // Checks

     IsFit2D_nh = fitmarg_nh->IsFit2D();
     IsFit2D_ih = fitmarg_ih->IsFit2D();
         
     if(IsFit2D_nh || IsFit2D_ih){
       cerr << endl;
       cerr << "Cannot use 2D fits for MH fit" << endl;
       cerr << "Do it only with 1D fits" << endl;
       cerr << endl;
       exit(1);
     }


     trueMH = runopt->GlobalMassHierarchy(); 

         
     // Get values of global minimum
     
     cout << "Get values of global minimum" << endl;

     if(chisqmin_nh<chisqmin_ih){ // best-fit is normal hierarchy
       MH_glob = "normal";
       chisqmin_glob = chisqmin_nh;
     }
     else{ // best-fit is inverted hierarchy
       MH_glob = "inverted";
       chisqmin_glob = chisqmin_ih;
     }
     
     deltachisq01_glob = chisqmin_nh - chisqmin_ih;
     if(trueMH == "normal")        deltachisq0BF_glob = chisqmin_nh - chisqmin_glob;
     else if(trueMH == "inverted") deltachisq0BF_glob = chisqmin_ih - chisqmin_glob;

     
     // Get the fakedata histogram

     cout << "Get the fakedata histogram" << endl;

     for(int is = 0; is < kns; is++) {
       SuperkSample_t s = kSKSamples[is];
       if(s == kSKsUndefined) continue;
    
#ifndef __USE_2D_TEMPLATES__
       TH1D * hspec_obs = (TH1D*)skobs_fakedata->ErecoSpectrum(s);
       TH1D * hspec_exp = (TH1D*)skexp->CurrPredictedErecoSpectrum(s);
#else
       TH2D * hspec_obs = (TH2D*)skobs_fakedata->TwoDSpectrum(s);
       TH2D * hspec_exp = (TH2D*)skexp->CurrPredictedErecoVsQ2Spectrum(s);  
#endif
       hspec_obs->Copy(*(hErecoToy_sample[s]));
       hspec_exp->Copy(*(hErecoExp_sample[s]));
       delete hspec_obs;
       delete hspec_exp;
     }
          
     // Fill toy informations (toy id, run, # of experiments, tag)
     ToyID = iexperiment;
     tag = runopt->FileTag();
     nexpts = runopt->NExperiments();
     run = runopt->Run();

     // Fill the tree 
     int iflag_tree = fitnt->Fill();
     cout << endl;
     cout << "Output fill tree flag = " << iflag_tree << endl;
     cout << endl;
     
     
     // Print output
 
     cout << endl;
     cout << "/// True parameters ///" << endl;
     cout << " - truepar0: " << runopt->Freeparname0() << " = " << truepar0 << endl;
     cout << " - truepar1: " << runopt->Freeparname1() << " = " << truepar1 << endl;
     cout << " - Mass Hierarchy: " << trueMH << endl;
     cout << endl;

     cout << "/// Normal Hierarchy ///" << endl;
     cout << " - Minimum Marginalized Chisq = " << chisqmin_nh << endl;        
     cout << endl;
     
     cout << "/// Inverted Hierarchy ///" << endl;
     cout << " - Minimum Marginalized Chisq = " << chisqmin_ih << endl;        
     cout << endl;
     
     cout << "/// Global minimum ///" << endl;
     cout << " - DeltaChisq = chisq(true) - chisq(bf) = " << deltachisq0BF_glob << endl;
     cout << " - DeltaChisq = chisq(NH) - chisq(IH) = " << deltachisq01_glob << endl;
     cout << " - Minimum Marginalized Chisq = " << chisqmin_glob << endl;  
     cout << " - Mass Hierarchy: " << MH_glob << endl;
     cout << endl;
     
    //Print out # events for sample
    
     for(int is = 0; is < kns; is++) {
       SuperkSample_t s = kSKSamples[is];
       if(s == kSKsUndefined) continue;
       
       cout << "Sample " << util::AsStringCompact(s) << endl;
       cout << " - # of observed events = ";       
#ifndef __USE_2D_TEMPLATES__
       cout << skobs_fakedata->ErecoSpectrum(s)->Integral() << endl;
#else
       cout << skobs_fakedata->TwoDSpectrum(s)->Integral() << endl;
#endif
       cout << " - # of expected events = ";       
#ifndef __USE_2D_TEMPLATES__
       cout << skexp->CurrPredictedErecoSpectrum(s)->Integral() << endl;
#else
       cout << skexp->CurrPredictedErecoVsQ2Spectrum(s)->Integral() << endl;
#endif
       cout << endl;
     }
     cout << endl;
     
     // Store mock data in separated root file 
     // as when generate_toymc_ensemble.cxx is run
     ftoys->cd();
#ifndef __USE_2D_TEMPLATES__ //1D templates
     TH1D * mock_data[kns];
#else
     TH2D * mock_data[kns];
#endif
     for (int is = 0 ; is < conv::kns ; is++) {
       SuperkSample_t s = conv::kSKSamples[is];
       if (s == conv::kSKsUndefined) continue;
#ifndef __USE_2D_TEMPLATES__ 
       mock_data[is] = skobs_fakedata->ErecoSpectrum(s);
#else
       mock_data[is] = skobs_fakedata->TwoDSpectrum(s);
#endif
       mock_data[is]->SetDirectory(0);     
       mock_data[is]->Write(Form("mock_data_%s_%d", util::AsStringCompact(s).c_str(), iexperiment));
       delete mock_data[is];
     }//is
     
     // Clean-up
     delete skobs_fakedata;
     fitmarg_nh->CleanUp();
     fitmarg_ih->CleanUp();
     delete fitmarg_nh;
     delete fitmarg_ih;
 
  } // loop over toys
  
  delete hErecoToy_1Rmu;
  delete hErecoToy_1Re;
  delete hErecoToy_RHC1Rmu;
  delete hErecoToy_RHC1Re;
  delete hErecoExp_1Rmu;
  delete hErecoExp_1Re;
  delete hErecoExp_RHC1Rmu;
  delete hErecoExp_RHC1Re;

  cout << endl;
  cout << "Total CPU time" << endl;
  cout << endl;
  sw_all.Stop();
  sw_all.Print();
  sw_all.Start(0);

  // Store ntuple with best-fit results
  f->cd();
  fitnt->Write(); // store the tree
  util::OscParams::Instance()->Write(f); // store the toy osc param histograms
  f->Close();
  delete f;  
  delete runopt;
}
