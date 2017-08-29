
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

#include "valor/Utilities/UtlForMacros.h"
#include "valor/Oscillation/OscParams.h"
#include "valor/CmdLine/RunOpt_TO_DELETE.h"
#include "valor/Utilities/Env.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Conventions/PhysHypothesisId.h"

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

/*  // For testing purposes only
  const char* data_file = Form("%s/inputs/data_run1to8/SKdata_Run18.root",
     Env::Instance()->TopLevelDir().c_str());
  SuperkObservation new_data = SuperkObservation();
  new_data.Read(data_file, 0, false);
  // End testing*/
  
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
  
  // fit status variables
  int nfits_nh=-999;
  int nfits_ih=-999;
  int nchisqcalls_nh=-999;
  int nchisqcalls_ih=-999;
  int niterations_nh=-999;
  int niterations_ih=-999;
  int fitstatus_nh=-999;
  int fitstatus_ih=-999;
  
  // true variables
  double truepar0=-999;
  double truepar1=-999;
  string trueMH = "";
  string trueparname0 = "";
  string trueparname1 = "";
  double chisqtrue0=conv::kALargeNum;
  
  // global fit output
  string MH_glob = "";
  double deltachisq0BF_glob=conv::kALargeNum;
  double chisqmin_glob=conv::kALargeNum;
  double bestfit0_glob=-999;
  double bestfit1_glob=-999;
  double bestfiterrNeg0_glob=-999;
  double bestfiterrNeg1_glob=-999;
  double bestfiterrPos0_glob=-999;
  double bestfiterrPos1_glob=-999;
  double bestfiterrSym0_glob=-999;
  double bestfiterrSym1_glob=-999;
  double bestfitFullChisq0_glob=-999;
  double bestfitFullChisq1_glob=-999;
  double chisqminFullChisq_glob=conv::kALargeNum;
  
  // normal hierarchy fit output
  double chisqmin_nh=conv::kALargeNum;
  double bestfit0_nh=-999;
  double bestfit1_nh=-999;
  double bestfiterrNeg0_nh=-999;
  double bestfiterrNeg1_nh=-999;
  double bestfiterrPos0_nh=-999;
  double bestfiterrPos1_nh=-999;
  double bestfiterrSym0_nh=-999;
  double bestfiterrSym1_nh=-999;
  double bestfitFullChisq0_nh=-999;
  double bestfitFullChisq1_nh=-999;
  double chisqminFullChisq_nh=conv::kALargeNum;
  
  // inverted hierarchy fit output
  double chisqmin_ih=conv::kALargeNum;
  double bestfit0_ih=-999;
  double bestfit1_ih=-999;
  double bestfiterrNeg0_ih=-999;
  double bestfiterrNeg1_ih=-999;
  double bestfiterrPos0_ih=-999;
  double bestfiterrPos1_ih=-999;
  double bestfiterrSym0_ih=-999;
  double bestfiterrSym1_ih=-999;
  double bestfitFullChisq0_ih=-999;
  double bestfitFullChisq1_ih=-999;
  double chisqminFullChisq_ih=conv::kALargeNum;
  
  // histograms
  TH2D *hchisq2D_out_nh = new TH2D("hchisq2D_out_nh","hchisq2D_out_nh",1,0,1,1,0,1);
  TH2D *hchisq2D_out_ih = new TH2D("hchisq2D_out_ih","hchisq2D_out_ih",1,0,1,1,0,1);
  TH1D *hchisq1D_out_nh = new TH1D("hchisq1D_out_nh","hchisq1D_out_nh",1,0,1);
  TH1D *hchisq1D_out_ih = new TH1D("hchisq1D_out_ih","hchisq1D_out_ih",1,0,1);
  TH2D *hdeltachisq2D_out_nh = new TH2D("hdeltachisq2D_out_nh","hdeltachisq2D_out_nh",1,0,1,1,0,1);
  TH2D *hdeltachisq2D_out_ih = new TH2D("hdeltachisq2D_out_ih","hdeltachisq2D_out_ih",1,0,1,1,0,1);
  TH1D *hdeltachisq1D_out_nh = new TH1D("hdeltachisq1D_out_nh","hdeltachisq1D_out_nh",1,0,1);
  TH1D *hdeltachisq1D_out_ih = new TH1D("hdeltachisq1D_out_ih","hdeltachisq1D_out_ih",1,0,1);
  TH2D *hstep2D_out_nh = new TH2D("hstep2D_out_nh","hstep2D_out_nh",1,0,1,1,0,1);
  TH2D *hstep2D_out_ih = new TH2D("hstep2D_out_ih","hstep2D_out_ih",1,0,1,1,0,1);
  TH1D *hstep1D_out_nh = new TH1D("hstep1D_out_nh","hstep1D_out_nh",1,0,1);
  TH1D *hstep1D_out_ih = new TH1D("hstep1D_out_ih","hstep1D_out_ih",1,0,1);
  
#ifndef __USE_2D_TEMPLATES__
  TH1D *hErecoToy_1Rmu = new TH1D("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1);
  TH1D *hErecoToy_1Re = new TH1D("hErecoToy_1Re","hErecoToy_1Re",1,0,1);
  TH1D *hErecoToy_RHC1Rmu = new TH1D("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1);
  TH1D *hErecoToy_RHC1Re = new TH1D("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1);
  TH1D *hErecoToy_MultiRe = new TH1D("hErecoToy_MultiRe","hErecoToy_MultiRe",1,0,1);
  TH1D *hErecoToy_sample[conv::kns] = {};
  TH1D *hErecoExp_1Rmu = new TH1D("hErecoExp_1Rmu","hErecoExp_1Rmu",1,0,1);
  TH1D *hErecoExp_1Re = new TH1D("hErecoExp_1Re","hErecoExp_1Re",1,0,1);
  TH1D *hErecoExp_RHC1Rmu = new TH1D("hErecoExp_RHC1Rmu","hErecoExp_RHC1Rmu",1,0,1);
  TH1D *hErecoExp_RHC1Re = new TH1D("hErecoExp_RHC1Re","hErecoExp_RHC1Re",1,0,1);
  TH1D *hErecoExp_MultiRe = new TH1D("hErecoExp_MultiRe","hErecoExp_MultiRe",1,0,1);
  TH1D *hErecoExp_sample[conv::kns] = {};
#else
  TH2D *hErecoToy_1Rmu = new TH2D("hErecoToy_1Rmu","hErecoToy_1Rmu",1,0,1,1,0,1);
  TH2D *hErecoToy_1Re = new TH2D("hErecoToy_1Re","hErecoToy_1Re",1,0,1,1,0,1);
  TH2D *hErecoToy_RHC1Rmu = new TH2D("hErecoToy_RHC1Rmu","hErecoToy_RHC1Rmu",1,0,1,1,0,1);
  TH2D *hErecoToy_RHC1Re = new TH2D("hErecoToy_RHC1Re","hErecoToy_RHC1Re",1,0,1,1,0,1);
  TH2D *hErecoToy_MultiRe = new TH2D("hErecoToy_MultiRe","hErecoToy_MultiRe",1,0,1,1,0,1);
  TH2D *hErecoToy_sample[conv::kns] = {};
  TH2D *hErecoExp_1Rmu = new TH2D("hErecoExp2D_1Rmu","hErecoExp2D_1Rmu",1,0,1,1,0,1);
  TH2D *hErecoExp_1Re = new TH2D("hErecoExp2D_1Re","hErecoExp2D_1Re",1,0,1,1,0,1);
  TH2D *hErecoExp_RHC1Rmu = new TH2D("hErecoExp2D_RHC1Rmu","hErecoExp2D_RHC1Rmu",1,0,1,1,0,1);
  TH2D *hErecoExp_RHC1Re = new TH2D("hErecoExp2D_RHC1Re","hErecoExp2D_RHC1Re",1,0,1,1,0,1);
  TH2D *hErecoExp_MultiRe = new TH2D("hErecoExp2D_MultiRe","hErecoExp2D_MultiRe",1,0,1,1,0,1);
  TH2D *hErecoExp_sample[conv::kns] = {};
#endif
  
  hErecoToy_sample[conv::kSKsFCFV1RingMuLike   ] = hErecoToy_1Rmu;
  hErecoToy_sample[conv::kSKsFCFV1RingELike    ] = hErecoToy_1Re ;
  hErecoToy_sample[conv::kSKsRHCFCFV1RingMuLike] = hErecoToy_RHC1Rmu;
  hErecoToy_sample[conv::kSKsRHCFCFV1RingELike ] = hErecoToy_RHC1Re ;
  hErecoToy_sample[conv::kSKsFCFVMultiRingELike ] = hErecoToy_MultiRe ;
  
  hErecoExp_sample[conv::kSKsFCFV1RingMuLike   ] = hErecoExp_1Rmu;
  hErecoExp_sample[conv::kSKsFCFV1RingELike    ] = hErecoExp_1Re ;
  hErecoExp_sample[conv::kSKsRHCFCFV1RingMuLike] = hErecoExp_RHC1Rmu;
  hErecoExp_sample[conv::kSKsRHCFCFV1RingELike ] = hErecoExp_RHC1Re ;
  hErecoExp_sample[conv::kSKsFCFVMultiRingELike ] = hErecoExp_MultiRe ;
  
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
  
  // fit status
  fitnt->Branch("IsFit2D_nh",&IsFit2D_nh);
  fitnt->Branch("IsFit2D_ih",&IsFit2D_ih);
  fitnt->Branch("fitstatus_nh",&fitstatus_nh,"fitstatus_nh/I");
  fitnt->Branch("fitstatus_ih",&fitstatus_ih,"fitstatus_ih/I");
  fitnt->Branch("nfits_nh",&nfits_nh,"nfits_nh/I");
  fitnt->Branch("nfits_ih",&nfits_ih,"nfits_ih/I");
  fitnt->Branch("nchisqcalls_nh",&nchisqcalls_nh,"nchisqcalls_nh/I");
  fitnt->Branch("nchisqcalls_ih",&nchisqcalls_ih,"nchisqcalls_ih/I");
  fitnt->Branch("niterations_nh",&niterations_nh,"niterations_nh/I");
  fitnt->Branch("niterations_ih",&niterations_ih,"niterations_ih/I");
  
  // true variables
  fitnt->Branch("truepar0", &truepar0 );
  fitnt->Branch("truepar1", &truepar1 );
  fitnt->Branch("trueMH",   &trueMH   );
  fitnt->Branch("chisqtrue0",&chisqtrue0,"chisqtrue0/D");
  fitnt->Branch("trueparname0",&trueparname0);
  fitnt->Branch("trueparname1",&trueparname1);
  
  
  // global fit output
  fitnt->Branch("deltachisq0BF_glob",&deltachisq0BF_glob,"deltachisq0BF_glob/D");
  fitnt->Branch("MH_glob",   &MH_glob   );
  fitnt->Branch("chisqmin_glob",&chisqmin_glob,"chisqmin_glob/D");
  fitnt->Branch("bestfit0_glob", &bestfit0_glob );
  fitnt->Branch("bestfit1_glob", &bestfit1_glob );
  fitnt->Branch("bestfiterrNeg0_glob",&bestfiterrNeg0_glob);
  fitnt->Branch("bestfiterrNeg1_glob",&bestfiterrNeg1_glob);
  fitnt->Branch("bestfiterrPos0_glob",&bestfiterrPos0_glob);
  fitnt->Branch("bestfiterrPos1_glob",&bestfiterrPos1_glob);
  fitnt->Branch("bestfiterrSym0_glob",&bestfiterrSym0_glob);
  fitnt->Branch("bestfiterrSym1_glob",&bestfiterrSym1_glob);
  fitnt->Branch("bestfitFullChisq0_glob", &bestfitFullChisq0_glob);
  fitnt->Branch("chisqminFullChisq_glob",&chisqminFullChisq_glob);
  
  // normal hierarchy fit output
  fitnt->Branch("chisqmin_nh",&chisqmin_nh,"chisqmin_nh/D");
  fitnt->Branch("bestfit0_nh", &bestfit0_nh );
  fitnt->Branch("bestfit1_nh", &bestfit1_nh );
  fitnt->Branch("bestfiterrNeg0_nh",&bestfiterrNeg0_nh);
  fitnt->Branch("bestfiterrNeg1_nh",&bestfiterrNeg1_nh);
  fitnt->Branch("bestfiterrPos0_nh",&bestfiterrPos0_nh);
  fitnt->Branch("bestfiterrPos1_nh",&bestfiterrPos1_nh);
  fitnt->Branch("bestfiterrSym0_nh",&bestfiterrSym0_nh);
  fitnt->Branch("bestfiterrSym1_nh",&bestfiterrSym1_nh);
  fitnt->Branch("bestfitFullChisq0_nh", &bestfitFullChisq0_nh );
  fitnt->Branch("chisqminFullChisq_nh",&chisqminFullChisq_nh);
  
  // inverted hierarchy fit output
  fitnt->Branch("chisqmin_ih",&chisqmin_ih,"chisqmin_ih/D");
  fitnt->Branch("bestfit0_ih", &bestfit0_ih );
  fitnt->Branch("bestfit1_ih", &bestfit1_ih );
  fitnt->Branch("bestfiterrNeg0_ih",&bestfiterrNeg0_ih);
  fitnt->Branch("bestfiterrNeg1_ih",&bestfiterrNeg1_ih);
  fitnt->Branch("bestfiterrPos0_ih",&bestfiterrPos0_ih);
  fitnt->Branch("bestfiterrPos1_ih",&bestfiterrPos1_ih);
  fitnt->Branch("bestfiterrSym0_ih",&bestfiterrSym0_ih);
  fitnt->Branch("bestfiterrSym1_ih",&bestfiterrSym1_ih);
  fitnt->Branch("bestfitFullChisq0_ih", &bestfitFullChisq0_ih );
  fitnt->Branch("chisqminFullChisq_ih",&chisqminFullChisq_ih);
  
  // histograms
  fitnt->Branch("hchisq1D_nh",     "TH1D", &hchisq1D_out_nh,     32000,0 );
  fitnt->Branch("hchisq1D_ih",     "TH1D", &hchisq1D_out_ih,     32000,0 );
  fitnt->Branch("hstep1D_nh",      "TH1D", &hstep1D_out_nh,      32000,0 );
  fitnt->Branch("hstep1D_ih",      "TH1D", &hstep1D_out_ih,      32000,0 );
  fitnt->Branch("hdeltachisq1D_nh","TH1D", &hdeltachisq1D_out_nh,32000,0 );
  fitnt->Branch("hdeltachisq1D_ih","TH1D", &hdeltachisq1D_out_ih,32000,0 );
  fitnt->Branch("hchisq2D_nh",     "TH2D", &hchisq2D_out_nh,     32000,0 );
  fitnt->Branch("hchisq2D_ih",     "TH2D", &hchisq2D_out_ih,     32000,0 );
  fitnt->Branch("hdeltachisq2D_nh","TH2D", &hdeltachisq2D_out_nh,32000,0 );
  fitnt->Branch("hdeltachisq2D_ih","TH2D", &hdeltachisq2D_out_ih,32000,0 );
  fitnt->Branch("hstep2D_nh",      "TH2D", &hstep2D_out_nh,      32000,0 );
  fitnt->Branch("hstep2D_ih",      "TH2D", &hstep2D_out_ih,      32000,0 );
  
#ifndef __USE_2D_TEMPLATES__
  fitnt->Branch("hErecoToy_1Rmu",     "TH1D",&hErecoToy_1Rmu,     32000,0);
  fitnt->Branch("hErecoToy_1Re",      "TH1D",&hErecoToy_1Re,      32000,0);
  fitnt->Branch("hErecoToy_RHC1Rmu",  "TH1D",&hErecoToy_RHC1Rmu,  32000,0);
  fitnt->Branch("hErecoToy_RHC1Re",   "TH1D",&hErecoToy_RHC1Re,   32000,0);
  fitnt->Branch("hErecoToy_MultiRe",   "TH1D",&hErecoToy_MultiRe,   32000,0);
  fitnt->Branch("hErecoExp_1Rmu",     "TH1D",&hErecoExp_1Rmu,     32000,0);
  fitnt->Branch("hErecoExp_1Re",      "TH1D",&hErecoExp_1Re,      32000,0);
  fitnt->Branch("hErecoExp_RHC1Rmu",  "TH1D",&hErecoExp_RHC1Rmu,  32000,0);
  fitnt->Branch("hErecoExp_RHC1Re",   "TH1D",&hErecoExp_RHC1Re,   32000,0);
  fitnt->Branch("hErecoExp_MultiRe",   "TH1D",&hErecoExp_MultiRe,   32000,0);
#else
  fitnt->Branch("hErecoToy_1Rmu",     "TH2D",&hErecoToy_1Rmu,     32000,0);
  fitnt->Branch("hErecoToy_1Re",      "TH2D",&hErecoToy_1Re,      32000,0);
  fitnt->Branch("hErecoToy_RHC1Rmu",  "TH2D",&hErecoToy_RHC1Rmu,  32000,0);
  fitnt->Branch("hErecoToy_RHC1Re",   "TH2D",&hErecoToy_RHC1Re,   32000,0);
  fitnt->Branch("hErecoToy_MultiRe",   "TH2D",&hErecoToy_MultiRe,   32000,0);
  fitnt->Branch("hErecoExp_1Rmu",     "TH2D",&hErecoExp_1Rmu,     32000,0);
  fitnt->Branch("hErecoExp_1Re",      "TH2D",&hErecoExp_1Re,      32000,0);
  fitnt->Branch("hErecoExp_RHC1Rmu",  "TH2D",&hErecoExp_RHC1Rmu,  32000,0);
  fitnt->Branch("hErecoExp_RHC1Re",   "TH2D",&hErecoExp_RHC1Re,   32000,0);
  fitnt->Branch("hErecoExp_MultiRe",   "TH2D",&hErecoExp_MultiRe,   32000,0);
#endif
  
  valor::conv::PhysParamId_t x_penalty = kPhypUndefined;
  valor::conv::PhysParamId_t y_penalty = kPhypUndefined;
  if(runopt->ToyPenaltyTerm2D()){
    x_penalty = valor::util::PhysParamIdFromString(runopt->ToyPenaltyTerm2DHistX());
    y_penalty = valor::util::PhysParamIdFromString(runopt->ToyPenaltyTerm2DHistY());
  }
  
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
    SuperkObservation * skobs_fakedata = new SuperkObservation();
    PhysHypothesis ph(runopt->Hypothesis());
    
    util::OscParams::Instance()->Print(cout);
    
    // Update physics hypothesis and fitparamlist with chosen value of oscillation parameters
    util::OscParams::Instance()->UpdatePhysHypothesis(ph);
    
    if(runopt->DoToys()){
      if(!runopt->ToyPenaltyTerm2D()){
        // Randomize specified osc params for current toy MC expt
        util::OscParams::Instance()->Randomize();
        
        // Set the mass hierarchy (to avoid conflicts with fitmarginalization mass hierarchy)
        util::SetMassHierarchy(skexp);
        
        // Update physics hypothesis and fitparamlist with chosen value of oscillation parameters
        util::OscParams::Instance()->UpdatePhysHypothesis(ph);
        set_init = false;
        util::OscParams::Instance()->UpdateFitParamList(fparm,set_init);
      }
      else{
        bool passed = false;
        while (passed == false)
        {
          // Randomize specified osc params for current toy MC expt
          util::OscParams::Instance()->Randomize();
          
          // Set the mass hierarchy (to avoid conflicts with fitmarginalization mass hierarchy)
          util::SetMassHierarchy(skexp);
          
          // Update physics hypothesis and fitparamlist with chosen value of oscillation parameters
          util::OscParams::Instance()->UpdatePhysHypothesis(ph);
          set_init = false;
          util::OscParams::Instance()->UpdateFitParamList(fparm,set_init);
          
          double val_x = fparm->CurrValue(x_penalty);
          double val_y = fparm->CurrValue(y_penalty);
          double random_number = Rnd::Instance()->RndOther().Uniform(0,1);
          int xbin = runopt->ToyPenaltyTerm2DHist()->GetXaxis()->FindBin(val_x);
          int ybin = runopt->ToyPenaltyTerm2DHist()->GetYaxis()->FindBin(val_y);
          
          if(runopt->ToyPenaltyTerm2DHist()->GetBinContent(xbin,ybin) > random_number) passed = true;
        }
      }
      
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
     
/*      // Replace some fake data sample spectra with real spectra for testing
      bool replace_samples[conv::kns];
      for(int is = 0; is < conv::kns; ++is)
      {
         SuperkSample_t s = kSKSamples[is];
//         if(util::IsMultiReSample(s) // Keep CC1pi MC)
//         if(util::Is1ReSample(s) && util::IsRHCSample(s))   // Keep RHC1Re MC
         if(util::Is1RMuSample(s) && util::IsRHCSample(s))   // Keep FHC1Rmu MC
         {
            replace_samples[is] = false;
         }
         else
         {  // Replace MC with data
            replace_samples[is] = true;
         }
      }

      // For testing purposes only, replace selected samples with real data
      skobs_fakedata->ReplaceSamplesWithData(new_data, replace_samples);
      // End testing*/
      //  cout << *fparm << endl;
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
    
    FitMarginalization *fitmarg_nh = new FitMarginalization();
    FitMarginalization *fitmarg_ih = new FitMarginalization();
    
    // Reset all the variables
    
    nfits_nh=-999;
    nfits_ih=-999;
    nchisqcalls_nh=-999;
    nchisqcalls_ih=-999;
    niterations_nh=-999;
    niterations_ih=-999;
    fitstatus_nh=-999;
    fitstatus_ih=-999;
    
    truepar0=-999;
    truepar1=-999;
    trueMH="";
    chisqtrue0=conv::kALargeNum;
    trueparname0="";
    trueparname1="";
    
    deltachisq0BF_glob=conv::kALargeNum;
    MH_glob = "";
    chisqmin_glob=conv::kALargeNum;
    bestfit0_glob=-999;
    bestfit1_glob=-999;
    bestfiterrNeg0_glob=-999;
    bestfiterrNeg1_glob=-999;
    bestfiterrPos0_glob=-999;
    bestfiterrPos1_glob=-999;
    bestfiterrSym0_glob=-999;
    bestfiterrSym1_glob=-999;
    bestfitFullChisq0_glob=-999;
    chisqminFullChisq_glob=conv::kALargeNum;
    
    chisqmin_nh=conv::kALargeNum;
    bestfit0_nh=-999;
    bestfit1_nh=-999;
    bestfiterrNeg0_nh=-999;
    bestfiterrNeg1_nh=-999;
    bestfiterrPos0_nh=-999;
    bestfiterrPos1_nh=-999;
    bestfiterrSym0_nh=-999;
    bestfiterrSym1_nh=-999;
    bestfitFullChisq0_nh=-999;
    chisqminFullChisq_nh=conv::kALargeNum;
    
    chisqmin_ih=conv::kALargeNum;
    bestfit0_ih=-999;
    bestfit1_ih=-999;
    bestfiterrNeg0_ih=-999;
    bestfiterrNeg1_ih=-999;
    bestfiterrPos0_ih=-999;
    bestfiterrPos1_ih=-999;
    bestfiterrSym0_ih=-999;
    bestfiterrSym1_ih=-999;
    bestfitFullChisq0_ih=-999;
    bestfitFullChisq1_ih=-999;
    chisqminFullChisq_ih=conv::kALargeNum;
    
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
    
    cout << endl;
    cout << "/////////////////////////////////////" << endl;
    cout << "//                                 //" << endl;
    cout << "//  Start fit of Normal Hierarchy  //" << endl;
    cout << "//                                 //" << endl;
    cout << "/////////////////////////////////////" << endl;
    cout << endl;

    util::PerformMargFit(argc,argv,runopt,fitmarg_nh,skexp,skobs_fakedata,false,"normal"); // Normal Hierarchy marginalization toys

    nchisqcalls_nh = fitmarg_nh->GetNCallsCalcChisqMargPt();
    
    cout << endl;
    cout << "Fit of Normal Hierarchy is done!" << endl;
    cout << endl;
    
    cout << endl;
    cout << "CPU time" << endl;
    cout << endl;
    sw_fit.Stop();
    sw_fit.Print();
    CPUtime_Minim_nh = sw_fit.CpuTime();
    sw_fit.Start(0);
    
    if(runopt->DrawFullChisq()){
      
      cout << "Fill Full Chisq Marg for Normal Hierarchy" << endl;
      fitmarg_nh->FillChisqMarg(skobs_fakedata,true);
      
      cout << endl;
      cout << "CPU time" << endl;
      cout << endl;
      sw_fit.Stop();
      sw_fit.Print();
      CPUtime_FullChisq_nh = sw_fit.CpuTime();
      REALtime_FullChisq_nh = sw_fit.RealTime();
      sw_fit.Start(0);
    }
    
    cout << endl;
    cout << "///////////////////////////////////////" << endl;
    cout << "//                                   //" << endl;
    cout << "//  Start fit of Inverted Hierarchy  //" << endl;
    cout << "//                                   //" << endl;
    cout << "///////////////////////////////////////" << endl;
    cout << endl;
    
    util::PerformMargFit(argc,argv,runopt,fitmarg_ih,skexp,skobs_fakedata,false,"inverted"); // Inverted Hierarchy marginalization toys
    
    nchisqcalls_ih = fitmarg_ih->GetNCallsCalcChisqMargPt();
    
    cout << endl;
    cout << "Fit of Inverted Hierarchy is done!" << endl;
    cout << endl;
    
    cout << endl;
    cout << "CPU time" << endl;
    cout << endl;
    sw_fit.Stop();
    sw_fit.Print();
    CPUtime_Minim_ih = sw_fit.CpuTime();
    sw_fit.Start(0);
    
    if(runopt->DrawFullChisq()){
      
      cout << "Fill Full Chisq Marg for Inverted Hierarchy" << endl;
      fitmarg_ih->FillChisqMarg(skobs_fakedata,true);
      
      cout << endl;
      cout << "CPU time" << endl;
      cout << endl;
      sw_fit.Stop();
      sw_fit.Print();
      CPUtime_FullChisq_ih = sw_fit.CpuTime();
      REALtime_FullChisq_ih = sw_fit.RealTime();
      sw_fit.Start(0);
    }
    
    // Checks
    
    IsFit2D_nh = fitmarg_nh->IsFit2D();
    IsFit2D_ih = fitmarg_ih->IsFit2D();
    
    if(IsFit2D_nh != IsFit2D_ih){
      cerr << endl;
      cerr << "Fits for NH and IH are different!" << endl;
      
      if(IsFit2D_nh) cerr << "NH fit is 2D" << endl;
      else cerr << "NH fit is not 2D" << endl;
      
      if(IsFit2D_ih) cerr << "IH fit is 2D" << endl;
      else cerr << "IH fit is not 2D" << endl;
      
      cerr << endl;
      exit(1);
    }
    
    // Get Best Fit point for each mass hierachy
    
    cout << "Get Best Fit points" << endl;
    bestfit0_nh = fitmarg_nh->GetParamBestFit(0);
    bestfit0_ih = fitmarg_ih->GetParamBestFit(0);
    bestfit1_nh = fitmarg_nh->GetParamBestFit(1);
    bestfit1_ih = fitmarg_ih->GetParamBestFit(1);
    chisqmin_nh = fitmarg_nh->GetChisqBestFit();
    chisqmin_ih = fitmarg_ih->GetChisqBestFit();
    
    // Get true chisq (only for true mass hierarchy)
    
    PhysParamId_t freepar0 = conv::kPhypUndefined;
    PhysParamId_t freepar1 = conv::kPhypUndefined;
    
    freepar0 = util::GetFreePar(runopt->Freeparname0());
    if(fitmarg_nh->IsFit2D()) freepar1 = util::GetFreePar(runopt->Freeparname1());
    
    truepar0 = fparm->InitValue(freepar0);
    if(fitmarg_nh->IsFit2D()) truepar1 = fparm->InitValue(freepar1);
    trueMH = runopt->GlobalMassHierarchy();
    
    trueparname0 = runopt->Freeparname0();
    trueparname1 = runopt->Freeparname1();
    
    
    // Not for BetaBar rate-only analysis with likelihood reweight method
    
    cout << "Calculate Chisq for true point (mass hierarchy as well):" << endl;
    
    if(trueMH=="normal"){
      Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal"); // needed to build the right binary file name
      chisqtrue0 = (double)fitmarg_nh->GetChisqMargTrueRead(truepar0,truepar1,skobs_fakedata,true);
    }
    else if(trueMH=="inverted"){
      Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted"); // needed to build the right binary file name
      chisqtrue0 = (double)fitmarg_ih->GetChisqMargTrueRead(truepar0,truepar1,skobs_fakedata,true);
    }
    else{
      cerr << endl;
      cerr << "The true Mass Hierarchy is not defined" << endl;
      cerr << endl;
      exit(1);
    }
    
    
    
    
    // Draw Full Chisq
    
    if(runopt->DrawFullChisq()){
      
      // cout << "Fill Full Chisq Marg for Normal Hierarchy" << endl;
      // fitmarg_nh->FillChisqMarg(skobs_fakedata,true);
      // cout << endl;
      // cout << "CPU time" << endl;
      // cout << endl;
      // sw_fit.Stop();
      // sw_fit.Print();
      // CPUtime_FullChisq_nh = sw_fit.CpuTime();
      // REALtime_FullChisq_nh = sw_fit.RealTime();
      // sw_fit.Start(0);
      // cout << "Fill Full Chisq Marg for Inverted Hierarchy" << endl;
      // fitmarg_ih->FillChisqMarg(skobs_fakedata,true);
      // cout << endl;
      // cout << "CPU time" << endl;
      // cout << endl;
      // sw_fit.Stop();
      // sw_fit.Print();
      // CPUtime_FullChisq_ih = sw_fit.CpuTime();
      // REALtime_FullChisq_ih = sw_fit.RealTime();
      // sw_fit.Start(0);
      
      cout << "Minimize Histo" << endl;
      fitmarg_nh->MinimizeHisto();
      fitmarg_ih->MinimizeHisto();
      
      bestfitFullChisq0_nh = fitmarg_nh->GetParamBestFit(0);
      bestfitFullChisq0_ih = fitmarg_ih->GetParamBestFit(0);
      
      if(fitmarg_nh->IsFit2D()) bestfitFullChisq1_nh = fitmarg_nh->GetParamBestFit(1);
      if(fitmarg_ih->IsFit2D()) bestfitFullChisq1_ih = fitmarg_ih->GetParamBestFit(1);
      
      chisqminFullChisq_nh = (double)fitmarg_nh->GetChisqBestFit();
      chisqminFullChisq_ih = (double)fitmarg_ih->GetChisqBestFit();
      
      cout << endl;
      cout << "Full Chisq Minima:" << endl;
      cout << endl;
      cout << " - Normal Hierarchy:" << endl;
      cout << " - bestfitFullChisq0 = " << bestfitFullChisq0_nh << endl;
      cout << " - bestfitFullChisq1 = " << bestfitFullChisq1_nh << endl;
      cout << " - chisqminFullChisq = " << chisqminFullChisq_nh << endl;
      cout << endl;
      cout << " - Inverted Hierarchy:" << endl;
      cout << " - bestfitFullChisq0 = " << bestfitFullChisq0_ih << endl;
      cout << " - bestfitFullChisq1 = " << bestfitFullChisq1_ih << endl;
      cout << " - chisqminFullChisq = " << chisqminFullChisq_ih << endl;
      cout << endl;
      
      cout << "Build deltaChisqMarg distribution" << endl;
      fitmarg_nh->BuildDeltaChisqMarg();
      fitmarg_ih->BuildDeltaChisqMarg();
      
    }
    
    
    
    
    
    
    // // Store toy likelihoods (needed to draw best-fit spectra w/ marginalization)
    // if(runopt->StoreLikeWeights()){
    //   fitmarg_nh->WriteLikeToys2Binary(TString::Format("%s_nh",runopt->FileTag().c_str()));
    //   fitmarg_ih->WriteLikeToys2Binary(TString::Format("%s_ih",runopt->FileTag().c_str()));
    // }
    
    
    
    
    
    //
    // Get fit errors
    //
    // Since the chisq was minimized reading binary files at each step you
    // must adopt the same approach for error computations as well.
    //
    
    if(runopt->TypeError().Contains("MINOS")){
      fitmarg_nh->CalcBestFitErrors(skobs_fakedata,"MINOS",true);
      fitmarg_ih->CalcBestFitErrors(skobs_fakedata,"MINOS",true);
    }
    if(runopt->TypeError().Contains("HESSE")){
      fitmarg_nh->CalcBestFitErrors(skobs_fakedata,"HESSE",true);
      fitmarg_ih->CalcBestFitErrors(skobs_fakedata,"HESSE",true);
    }
    
    bestfiterrNeg0_nh = fitmarg_nh->GetParamBestFitErrNeg(0);
    bestfiterrNeg0_ih = fitmarg_ih->GetParamBestFitErrNeg(0);
    
    bestfiterrNeg1_nh = fitmarg_nh->GetParamBestFitErrNeg(1);
    bestfiterrNeg1_ih = fitmarg_ih->GetParamBestFitErrNeg(1);
    
    bestfiterrPos0_nh = fitmarg_nh->GetParamBestFitErrPos(0);
    bestfiterrPos0_ih = fitmarg_ih->GetParamBestFitErrPos(0);
    
    bestfiterrPos1_nh = fitmarg_nh->GetParamBestFitErrPos(1);
    bestfiterrPos1_ih = fitmarg_ih->GetParamBestFitErrPos(1);
    
    bestfiterrSym0_nh = fitmarg_nh->GetParamBestFitErrSym(0);
    bestfiterrSym0_ih = fitmarg_ih->GetParamBestFitErrSym(0);
    
    bestfiterrSym1_nh = fitmarg_nh->GetParamBestFitErrSym(1);
    bestfiterrSym1_ih = fitmarg_ih->GetParamBestFitErrSym(1);
    
    
    // Get values of global minimum
    
    if(chisqmin_nh<chisqmin_ih){ // best-fit is normal hierarchy
      MH_glob = "normal";
      deltachisq0BF_glob = chisqtrue0 - chisqmin_nh; // chisqtrue0 is for true MH
      chisqmin_glob = chisqmin_nh;
      bestfit0_glob = bestfit0_nh;
      bestfit1_glob = bestfit1_nh;
      bestfiterrNeg0_glob = bestfiterrNeg0_nh;
      bestfiterrNeg1_glob = bestfiterrNeg1_nh;
      bestfiterrPos0_glob = bestfiterrPos0_nh;
      bestfiterrPos1_glob = bestfiterrPos1_nh;
      bestfiterrSym0_glob = bestfiterrSym0_nh;
      bestfiterrSym1_glob = bestfiterrSym1_nh;
      bestfitFullChisq0_glob = bestfitFullChisq0_nh;
      chisqminFullChisq_glob = chisqminFullChisq_nh;
    }
    else{ // best-fit is inverted hierarchy
      MH_glob = "inverted";
      deltachisq0BF_glob = chisqtrue0 - chisqmin_ih;
      chisqmin_glob = chisqmin_ih;
      bestfit0_glob = bestfit0_ih;
      bestfit1_glob = bestfit1_ih;
      bestfiterrNeg0_glob = bestfiterrNeg0_ih;
      bestfiterrNeg1_glob = bestfiterrNeg1_ih;
      bestfiterrPos0_glob = bestfiterrPos0_ih;
      bestfiterrPos1_glob = bestfiterrPos1_ih;
      bestfiterrSym0_glob = bestfiterrSym0_ih;
      bestfiterrSym1_glob = bestfiterrSym1_ih;
      bestfitFullChisq0_glob = bestfitFullChisq0_ih;
      chisqminFullChisq_glob = chisqminFullChisq_ih;
    }
    
    
    
    // Get the fakedata histogram
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
    
    // Fill the output tree
    
    if(fitmarg_nh->IsFit2D()){
      fitmarg_nh->GetChisqMarg2D     ()->Copy(*hchisq2D_out_nh      );
      fitmarg_ih->GetChisqMarg2D     ()->Copy(*hchisq2D_out_ih      );
      fitmarg_nh->GetDeltaChisqMarg2D()->Copy(*hdeltachisq2D_out_nh );
      fitmarg_ih->GetDeltaChisqMarg2D()->Copy(*hdeltachisq2D_out_ih );
      fitmarg_nh->GetChisqStep_2D    ()->Copy(*hstep2D_out_nh       );
      fitmarg_ih->GetChisqStep_2D    ()->Copy(*hstep2D_out_ih       );
      
      if(runopt->VerboseFit()){
        // Print out the full chisq
        cout << endl;
        cout << "Full chisq values: Normal Hierarchy:" << endl;
        int npt_chisq_x = hchisq2D_out_nh->GetNbinsX();
        int npt_chisq_y = hchisq2D_out_nh->GetNbinsY();
        for(int iy=0;iy<npt_chisq_y;iy++){
          for(int ix=0;ix<npt_chisq_x;ix++){
            cout << hchisq2D_out_nh->GetBinContent(ix+1,iy+1) << " ";
          }
          cout << endl;
        }
        cout << endl;
        cout << "Full chisq values: Inverted Hierarchy:" << endl;
        npt_chisq_x = hchisq2D_out_ih->GetNbinsX();
        npt_chisq_y = hchisq2D_out_ih->GetNbinsY();
        for(int iy=0;iy<npt_chisq_y;iy++){
          for(int ix=0;ix<npt_chisq_x;ix++){
            cout << hchisq2D_out_ih->GetBinContent(ix+1,iy+1) << " ";
          }
          cout << endl;
        }
      }
    }
    else{
      fitmarg_nh->GetChisqMarg1D     ()->Copy(*hchisq1D_out_nh      );
      fitmarg_ih->GetChisqMarg1D     ()->Copy(*hchisq1D_out_ih      );
      fitmarg_nh->GetDeltaChisqMarg1D()->Copy(*hdeltachisq1D_out_nh);
      fitmarg_ih->GetDeltaChisqMarg1D()->Copy(*hdeltachisq1D_out_ih);
      fitmarg_nh->GetChisqStep_1D    ()->Copy(*hstep1D_out_nh);
      fitmarg_ih->GetChisqStep_1D    ()->Copy(*hstep1D_out_ih);
      
      cout << endl;
      if(runopt->VerboseFit()){
        cout << "Full chisq values: Normal Hierarchy:" << endl;
        cout << " - Chisq:" << endl;
        int npt_chisq_x = hchisq1D_out_nh->GetNbinsX();
        for(int ix=0;ix<npt_chisq_x;ix++){
          cout << hchisq1D_out_nh->GetBinContent(ix+1) << " ";
        }
        cout << endl;
        cout << " - DeltaChisq:" << endl;
        for(int ix=0;ix<npt_chisq_x;ix++){
          cout << hdeltachisq1D_out_nh->GetBinContent(ix+1) << " ";
        }
        cout << endl;
        
        cout << endl;
        cout << "Full chisq values: Inverted Hierarchy:" << endl;
        cout << " - Chisq:" << endl;
        npt_chisq_x = hchisq1D_out_ih->GetNbinsX();
        for(int ix=0;ix<npt_chisq_x;ix++){
          cout << hchisq1D_out_ih->GetBinContent(ix+1) << " ";
        }
        cout << endl;
        cout << " - DeltaChisq:" << endl;
        for(int ix=0;ix<npt_chisq_x;ix++){
          cout << hdeltachisq1D_out_ih->GetBinContent(ix+1) << " ";
        }
        cout << endl;
      }
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
    cout << " - true chisq: " << chisqtrue0 << endl;
    cout << endl;
    
    cout << "/// Normal Hierarchy ///" << endl;
    cout << " - Minimum Marginalized Chisq = " << chisqmin_nh << endl;
    cout << " - bestfit 0 = " << bestfit0_nh;
    cout << " +" << bestfiterrPos0_nh
    << " -" << bestfiterrNeg0_nh
    << " (MINOS)"
    << endl;
    cout << " - bestfit 0 = " << bestfit0_nh
    << " +" << bestfiterrSym0_nh
    << " -" << bestfiterrSym0_nh
    << " (HESSE)"
    << endl;
    if(fitmarg_nh->IsFit2D()){
      cout << " - bestfit 1 = " << bestfit1_nh;
      cout << " +" << bestfiterrPos1_nh
      << " -" << bestfiterrNeg1_nh
      << " (MINOS)"
      << endl;
      cout << " - bestfit 1 = " << bestfit1_nh
      << " +" << bestfiterrSym1_nh
      << " -" << bestfiterrSym1_nh
      << " (HESSE)"
      << endl;
    }
    cout << " - # of calls of CalcChisqMargPt(): " << nchisqcalls_nh << endl;
    cout << "Minimum of full Chisq histogram:" << endl;
    cout << " - chisqmin = " << chisqminFullChisq_nh << endl;
    cout << " - bestfit0 = " << bestfitFullChisq0_nh << endl;
    if(fitmarg_nh->IsFit2D())
      cout << "  - bestfit1 = " << bestfitFullChisq1_nh << endl;
    cout << endl;
    
    cout << "/// Inverted Hierarchy ///" << endl;
    cout << " - Minimum Marginalized Chisq = " << chisqmin_ih << endl;
    cout << " - bestfit 0 = " << bestfit0_ih;
    cout << " +" << bestfiterrPos0_ih
    << " -" << bestfiterrNeg0_ih
    << " (MINOS)"
    << endl;
    cout << " - bestfit 0 = " << bestfit0_ih
    << " +" << bestfiterrSym0_ih
    << " -" << bestfiterrSym0_ih
    << " (HESSE)"
    << endl;
    if(fitmarg_ih->IsFit2D()){
      cout << " - bestfit 1 = " << bestfit1_ih;
      cout << " +" << bestfiterrPos1_ih
      << " -" << bestfiterrNeg1_ih
      << " (MINOS)"
      << endl;
      cout << " - bestfit 1 = " << bestfit1_ih
      << " +" << bestfiterrSym1_ih
      << " -" << bestfiterrSym1_ih
      << " (HESSE)"
      << endl;
    }
    cout << " - # of calls of CalcChisqMargPt(): " << nchisqcalls_ih << endl;
    cout << "Minimum of full Chisq histogram:" << endl;
    cout << " - chisqmin = " << chisqminFullChisq_ih << endl;
    cout << " - bestfit0 = " << bestfitFullChisq0_ih << endl;
    if(fitmarg_ih->IsFit2D())
      cout << " - bestfit1 = " << bestfitFullChisq1_ih << endl;
    cout << endl;
    
    cout << "/// Global minimum ///" << endl;
    cout << " - DeltaChisq = " << deltachisq0BF_glob << endl;
    cout << " - Minimum Marginalized Chisq = " << chisqmin_glob << endl;
    cout << " - Mass Hierarchy: " << MH_glob << endl;
    cout << " - bestfit 0 = " << bestfit0_glob;
    cout << " +" << bestfiterrPos0_glob
    << " -" << bestfiterrNeg0_glob
    << " (MINOS)"
    << endl;
    cout << " - bestfit 0 = " << bestfit0_glob
    << " +" << bestfiterrSym0_glob
    << " -" << bestfiterrSym0_glob
    << " (HESSE)"
    << endl;
    if(fitmarg_nh->IsFit2D()){
      cout << " - bestfit 1 = " << bestfit1_glob;
      cout << " +" << bestfiterrPos1_glob
      << " -" << bestfiterrNeg1_glob
      << " (MINOS)"
      << endl;
      cout << " - bestfit 1 = " << bestfit1_glob
      << " +" << bestfiterrSym1_glob
      << " -" << bestfiterrSym1_glob
      << " (HESSE)"
      << endl;
    }
    cout << "Minimum of full Chisq histogram:" << endl;
    cout << " - chisqmin = " << chisqminFullChisq_glob << endl;
    cout << " - bestfit0 = " << bestfitFullChisq0_glob << endl;
    if(fitmarg_nh->IsFit2D())
      cout << " - bestfit1 = " << bestfitFullChisq1_glob << endl;
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
    
    //
    
    // int npt_chisq_x = hchisq1D_out->GetNbinsX();
    // for(int ix=0;ix<npt_chisq_x;ix++){
    //   cout << hchisq1D_out->GetBinContent(ix+1) << " ";
    // }
    // cout << endl;
    
    // Clean-up
    delete skobs_fakedata;
    fitmarg_nh->CleanUp();
    fitmarg_ih->CleanUp();
    delete fitmarg_nh;
    delete fitmarg_ih;
    
    // hchisq1D_out->GetNbinsX();
    // for(int ix=0;ix<npt_chisq_x;ix++){
    //   cout << hchisq1D_out->GetBinContent(ix+1) << " ";
    // }
    // cout << endl;
    
  } // loop over toys
  
  delete hErecoToy_1Rmu;
  delete hErecoToy_1Re;
  delete hErecoToy_RHC1Rmu;
  delete hErecoToy_RHC1Re;
  delete hErecoToy_MultiRe;
  
  delete hErecoExp_1Rmu;
  delete hErecoExp_1Re;
  delete hErecoExp_RHC1Rmu;
  delete hErecoExp_RHC1Re;
  delete hErecoExp_MultiRe;
  
  delete hchisq2D_out_nh;
  delete hchisq2D_out_ih;
  delete hdeltachisq2D_out_nh;
  delete hdeltachisq2D_out_ih;
  delete hstep2D_out_nh;
  delete hstep2D_out_ih;
  delete hchisq1D_out_nh;
  delete hchisq1D_out_ih;
  delete hdeltachisq1D_out_nh;
  delete hdeltachisq1D_out_ih;
  delete hstep1D_out_nh;
  delete hstep1D_out_ih;
  
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
