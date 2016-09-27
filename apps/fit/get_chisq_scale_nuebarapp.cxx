//
// ------------------------------------------------------------------------------------------------
////Find minimum chisq of all toy throws
//This macro is used to find the minimum chisq for lots of toy throws, this is needed to make the
//nuebarapp and postrior predictive toy generation more efficient.
//The ouput of this app is a number, the minimum chisq.
//Give the app the setup eg. what to randomise and ranges of parametres and the "data"
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//


// Run in the following way:

//$T2KOSCANA/bin/generate_toymc_ensemble --enable-stat-fluct --syst-tweak-mode 1 --enable-rethrow-unphysical-toys --randomize-all-oscparams --prior-reactors --cpv-nh -3.14 --cpv-ih -3.14 --enable-marg-mass-hierarchy --normal-mass-hierarchy-prior 0.5 --enable-synch-osc-ih-to-nh --file-tag "ciao" --config-file analysis_configuration_test_theta_elike_RFGRPA_Nu6.57E20POT.xml --run 0 --enable-verbose-syst-throws 0 --dmsq-scheme "auto" --mix-param-scheme "sinsqtheta" --cpv-phase-scheme "delta" --mass-hierarchy "normal" --mix-param-octant "two" --hypothesis "3f_osc" --prior 1 --first-experiment 0 --n-experiments 10

#include "SuperkObservation.h"
#include "SuperkExpectation.h"
#include "T2KPhysicalParameterLimits.h"

#include "Utilities.h"
#include "Binning.h"
#include "AnalysisOptions.h"

#include "Build.h"
#include "UtlForMacrosT2XK.h"

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
    stream << "Macro-specific running options:" << endl
	   << " - " << setw(w) << "Random seed for position inside bin:" <<  fToySeed << endl
	   << "Output: Minimum chisq from all toy throws"  << endl;
    stream << endl;
  }
  int      ToySeed        () { return fToySeed; }
private:
  int    fToySeed;
  void Init() {
    RunOpt::Init();
    fToySeed        = 0;
  }
  void ReadFromCommandLine(int argc, char ** argv, int parseopt) {
    RunOpt::ReadFromCommandLine(argc, argv, parseopt);
    if(fParser->OptionExists("toy-random-seed"))
    fToySeed = fParser->ArgAsInt("toy-random-seed");
  }
  
  void Usage(const int w = 40) {
    RunOpt::Usage(w);
    cout << "Macro-specific running options:" << endl
    << " " << setw(w) << "--toy-random-seed"           << "Random seed used to specify where inside the bin the toys are generated" << endl
    << endl;
  }
};

int main(int argc, char ** argv)
{
  TStopwatch sw;
  TH1D* chisq_distribution = new TH1D("chisq_distribution","chisq_distribution",100000,0,1000);
  double chisq_min = 9999999; //Starting value of minimum chisq. This macro will fill this double with the min chisq from all of the toys
  Env::Instance()->AutoSet("t2xk");
  cout << *Env::Instance() << endl;
  gROOT->SetBatch(true);
  MyRunOpt *runopt = new MyRunOpt(argc, argv, kParseCommon | kParseToys | kParseGlobal | kParseFitMarg | kParse3FOscParams);
  
  TRandom3 rand(runopt->ToySeed()); //random number to see where inside a bin a toy is generated
  
  TString config  (Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), runopt->ConfigFile().c_str() ));
  TString xmlfile2(runopt->Xmlfile().c_str() );
  TString rootfile2(runopt->Rootfile().c_str());
  
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
  
  if(!runopt->DoToys()){
    cout << "- Reading analysis configuration from: " << config << endl;
    cout << "- Toy-MC" << endl;
    if(runopt->Using_xml_file()) cout << "   - Reading experiments from: " << xmlfile2 << endl;
    else cout << "   - Reading experiments from: " << rootfile2 << endl;
    cout << "   - Using experiments in range: [" << 0 << ", " << 0 << "]" << endl;
  }
  
  // Create & config superk spectrum predictor
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(runopt->ConfigFile().c_str());
  
  SuperkObservation * skobs_data = new SuperkObservation();
  
  if(!runopt->DoToys()){
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined || s == conv::kSKsRHCFCFV1RingELike) continue;
      
      if(runopt->Using_xml_file()){
        cout << "Use xml file" << endl; // Only read 1st toy
        skobs_data->ReadSample(xmlfile2.Data(),0,runopt->Using_xml_file(),s);
      }
      else{
        cout << "Use root file" << endl; // Only read 1st toy
        skobs_data->ReadSample(rootfile2.Data(),0,runopt->Using_xml_file(),s);
      }
    }
  }
  
  //Compare input data file with each expectation toy to weight the fake data
  TH2D * real_data[kns];
  
  for (int is = 0 ; is < conv::kns ; is++) {
    SuperkSample_t s = conv::kSKSamples[is];
    if (s == conv::kSKsUndefined || s == conv::kSKsRHCFCFV1RingELike) continue;
    real_data[is] = skobs_data->TwoDSpectrum(s);
  }
  
  int weird_applied_neg      = 0;
  int weird_unphysical_tweak = 0;
  for(int iexperiment = runopt->FirstExperiment(), last_experiment = runopt->LastExperiment(); iexperiment <= last_experiment; iexperiment++) {
    cout << "\n\n*** Pseudo-experiment = " << iexperiment << endl;
    
    // Generate mock data
    cout << "- Generating toy-MC..." << endl;
    SuperkObservation * skobs = new SuperkObservation();
    PhysHypothesis ph(runopt->Hypothesis());
    
    
    //util::OscParams::Instance()->Print(cout);
    
    // Randomize specified osc params for current toy MC expt
    util::OscParams::Instance()->Randomize();
    
    // Set the mass hierarchy
    util::SetMassHierarchy(skexp);
    
    // Update physics hypothesis with chosen value of oscillation parameters
    util::OscParams::Instance()->UpdatePhysHypothesis(ph);
    
    // Fill toy histograms of osc params
    util::OscParams::Instance()->FillOscToyHisto(ph);
    
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
      
      applied_negevents_condition = skobs->CreateMockData(*skexp, ph, runopt->StatFluctuations());
      cout << "applied_negevents_condition = " << applied_negevents_condition << endl;
      if(applied_negevents_condition)
      weird_applied_neg++;
      
      
    }//applied_negevents_condition
    
    
    //Get exp data
#ifndef __USE_2D_TEMPLATES__ //1D templates
    TH1D * exp_data[kns];
#else
    TH2D * exp_data[kns];
#endif
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
      
#ifndef __USE_2D_TEMPLATES__
      exp_data[is] = skexp->CurrPredictedErecoSpectrum(s);
#else
      exp_data[is] = skexp->CurrPredictedErecoVsQ2Spectrum(s);
#endif
      exp_data[is]->SetDirectory(0);
    }//is
    
    //Calculate the chisq for 3 samples
    double chisq = 0;
    
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined || s == conv::kSKsRHCFCFV1RingELike) continue;
      chisq += valor::util::Chisq_Penalty_TH2D(real_data[is],exp_data[is]);
    }
    chisq_distribution->Fill(chisq);
    if(chisq < chisq_min) chisq_min = chisq;
    
    for (int is = 0 ; is < conv::kns ; is++) {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
      delete exp_data[is];
    }
    
    delete skobs;
  }//iexperiment
  
  
  cout << endl
  << " - Summary - " << endl
  << "    * Created " << runopt->NExperiments() << " toys" << endl
  << "    * " << weird_unphysical_tweak << " unphysical_tweaks rethrown" << endl
  << "    * " << weird_applied_neg << " applied_negevents_condition rethrown" << endl
  << endl
  << " ************************************" << endl
  << "      chisq min = " << chisq_min << endl
  << " ************************************" << endl
  << endl;
  
  TFile output("output.root","recreate");
  chisq_distribution->Write();
  output.Close();
  
  sw.Stop();
  sw.Print();
  delete skexp;
  delete skobs_data;
}


