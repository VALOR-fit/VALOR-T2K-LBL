//____________________________________________________________________________ 
/*

 >> ChangeLog:

*/ 
//____________________________________________________________________________

#define __show_fit_init_info__ 
#define __show_fit_result__ 
#define __show_fit_steps__ 
#define __show_fit_cpu_time__

#include <cassert> 
#include <cstdlib> 
#include <iostream> 
#include <sstream>

#include <TStopwatch.h>
#include <Math/GSLMinimizer.h>
#include <Math/GSLNLSMinimizer.h>
#include <Math/GSLSimAnMinimizer.h>

#include "FitterGSL.h" 
#include "SuperkObservation.h" 
#include "SuperkExpectation.h" 
#include "FitInputs.h" 
#include "FitParamList.h" 
#include "Utilities.h"
#include "Conventions.h"
#include "BANFF.h"
#include "FSI.h"
#include "SuperkDetector.h"

using std::cerr; 
using std::cout; 
using std::endl; 
using std::ostringstream;

using namespace t2kana; 
using namespace t2kana::conv;

   
//____________________________________________________________________________ 
FitterGSL::FitterGSL() :
FitterABC()
{
  this->Init();
}
//____________________________________________________________________________ 
FitterGSL::FitterGSL(FitInputs * fit_inputs, Option_t * option) :
FitterABC()
{
  this->Init(fit_inputs, option);
}
//____________________________________________________________________________ 
FitterGSL::~FitterGSL() 
{

}
//____________________________________________________________________________ 
void FitterGSL::Init(void) 
{
  fFitInputs      = 0;
  fOption         = "";
  fNFreeFitParams = 0;
  fNFitBins       = 0;

  fFunc      = 0;
  fMinimizer = 0;
}
//____________________________________________________________________________ 
void FitterGSL::Init(FitInputs * fit_inputs, Option_t * option) 
{
#ifdef __show_fit_init_info__
  cout << "FitterGSL/Init >> Fitter initialization..." << endl;
#endif

  this->Init();

  fFitInputs = fit_inputs;
  if(!fFitInputs) {
     cerr << "FitterGSL/Init >> ** Null fit inputs!" << endl;
     exit(1);
  }

  fOption = option;

  FitParamList & fit_params = fFitInputs->FitParams();
  int param_list_size = fit_params.Size();

#ifdef __show_fit_init_info__
  cout << "FitterGSL/Init >> Fit parameter list size = " 
       << param_list_size << endl;
#endif

  if(fOption.Contains("SimAnneal")) {
    fMinimizer = new ROOT::Math::GSLSimAnMinimizer;
  }
  else
  if(fOption.Contains("LevMarq")) {
    fMinimizer = new ROOT::Math::GSLNLSMinimizer;
  }
  else
  if(fOption.Contains("ConjugateFR")) {
    fMinimizer = new ROOT::Math::GSLMinimizer(ROOT::Math::kConjugateFR);
  }
  else
  if(fOption.Contains("ConjugatePR")) {
    fMinimizer = new ROOT::Math::GSLMinimizer(ROOT::Math::kConjugatePR);
  }
  else
  if(fOption.Contains("VectorBFGS")) {
    fMinimizer = new ROOT::Math::GSLMinimizer(ROOT::Math::kVectorBFGS2);
  }
  else
  if(fOption.Contains("SteepestDescent")) {
    fMinimizer = new ROOT::Math::GSLMinimizer(ROOT::Math::kSteepestDescent);
  }
  else {
     cout << "FitterGSL/Init >> Can not handle option: " << option << endl;
     exit(1);
  }

  fMinimizer->SetMaxFunctionCalls(1000000);
  fMinimizer->SetMaxIterations(100000);
  fMinimizer->SetTolerance(0.001);

  // calculate degrees of freedom
  fNFitBins       = kNErecoBins1Rmu; // HERE, FIX
  fNFreeFitParams = fit_params.GetNumberFreeParameters();

  if (fOption.Contains("LevMarq") ){ 
    fFitMethodFunc = fFitInputs->AsFitMethodFunction(); 
  }
  else{
    fFunc = fFitInputs->AsFunctor(); 
    fMinimizer->SetFunction(*fFunc);
  }
  fit_params.SetTo(fMinimizer);

  
#ifdef __show_fit_init_info__
  cout << "FitterGSL/Init >> NDoF = " << this->NDOF()
       << " (Number of bins = " << this->NBinsInFit()
       << ", Number of free parameters = " << this->NFitParams() << ")"
       << endl;
#endif
}
//____________________________________________________________________________ 
bool FitterGSL::DoTheFit(void) 
{
#ifdef __show_fit_cpu_time__
  TStopwatch w;
  w.Start();
#endif
  fFitTime = 0;

  if (fOption.Contains("LevMarq") ){
    ROOT::Fit::Fitter fitter;
    fFitInputs->FitParams().SetTo(&fitter);
    //fitter.Config().MinimizerOptions().SetPrintLevel(1);
    fitter.Config().MinimizerOptions().SetMinimizerType("GSLMultiFit");
    
    fitter.FitFCN(*fFitMethodFunc,0,fFitInputs->GetNPoints(),true);
    ROOT::Fit::FitResult result = fitter.Result();
    result.Print(std::cout);
  }
  else{
    if(!fMinimizer) return false;
    // Run the fitter once 
    fMinimizer->Minimize();
  }
  
#ifdef __show_fit_cpu_time__
  w.Stop();
  cout << "FitterGSL/DoTheFit >> Fit performed in " 
       << w.CpuTime() << " sec." << endl;
  fFitTime = w.CpuTime();
#endif
  
  fFitInputs->FitParams().CopyCurrToBestFitValues();
  return true;
}
//____________________________________________________________________________ 
