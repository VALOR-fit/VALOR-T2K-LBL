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

#include "FitterLevembergMarquardt.h" 
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
FitterLevembergMarquardt::FitterLevembergMarquardt() :
FitterABC()
{
  this->Init();
}
//____________________________________________________________________________ 
FitterLevembergMarquardt::FitterLevembergMarquardt(
  FitInputs * fit_inputs, Option_t * option) :
FitterABC()
{
  this->Init(fit_inputs, option);
}
//____________________________________________________________________________ 
FitterLevembergMarquardt::~FitterLevembergMarquardt() 
{

}
//____________________________________________________________________________ 
void FitterLevembergMarquardt::Init(void) 
{
  fFitInputs      = 0;
  fOption         = "";
  fNFreeFitParams = 0;
  fNFitBins       = 0;

  fFunc      = 0;
  fMinimizer = 0;
}
//____________________________________________________________________________ 
void FitterLevembergMarquardt::Init(FitInputs * fit_inputs, Option_t * /*option*/) 
{
#ifdef __show_fit_init_info__
  cout << "FitterLevembergMarquardt/Init >> Fitter initialization..." << endl;
#endif

  this->Init();

  fFitInputs = fit_inputs;
  if(!fFitInputs) {
     cerr << "FitterLevembergMarquardt/Init >> ** Null fit inputs!" << endl;
     exit(1);
  }

  FitParamList & fit_params = fFitInputs->FitParams();
  int param_list_size = fit_params.Size();

#ifdef __show_fit_init_info__
  cout << "FitterLevembergMarquardt/Init >> Fit parameter list size = " 
       << param_list_size << endl;
#endif

  fMinimizer = new ROOT::Math::GSLNLSMinimizer;

  fMinimizer->SetMaxFunctionCalls(1000000);
  fMinimizer->SetMaxIterations(100000);
  fMinimizer->SetTolerance(0.001);

  fFunc = fFitInputs->AsFunctor(); 
  fMinimizer->SetFunction(*fFunc);

  // set params
  fit_params.SetTo(fMinimizer);

  // calculate degrees of freedom
  fNFitBins       = kNErecoBins1Rmu; // HERE, FIX
  fNFreeFitParams = fit_params.GetNumberFreeParameters();

#ifdef __show_fit_init_info__
  cout << "FitterLevembergMarquardt/Init >> NDoF = " << this->NDOF()
       << " (Number of bins = " << this->NBinsInFit()
       << ", Number of free parameters = " << this->NFitParams() << ")"
       << endl;
#endif
}
//____________________________________________________________________________ 
bool FitterLevembergMarquardt::DoTheFit(void) 
{
#ifdef __show_fit_cpu_time__
  TStopwatch w;
  w.Start();
#endif
  fFitTime = 0;

  if(!fMinimizer) return false;

  // Run the fitter once 
  fMinimizer->Minimize();

#ifdef __show_fit_cpu_time__
  w.Stop();
  cout << "FitterLevembergMarquardt/DoTheFit >> Fit performed in " 
       << w.CpuTime() << " sec." << endl;
  fFitTime = w.CpuTime();
#endif

  return true;
}
//____________________________________________________________________________ 
