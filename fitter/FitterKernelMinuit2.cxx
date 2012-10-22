//____________________________________________________________________________ 
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/ 
//____________________________________________________________________________

#define __show_fit_init_info__ 
#define __show_fit_init_debug__ 
#define __show_fit_result__ 
#define __show_fit_steps__ 
#define __show_fit_steps_debug__
#define __show_fit_cpu_time__
#define __get_hesse_errors__
//#define __get_minos_errors__

#include <cassert> 
#include <cstdlib> 
#include <iostream> 
#include <sstream>

#include <TMath.h> 
#include <TGraph.h> 
#include <TH2D.h> 
#include <THStack.h> 
#include <TMarker.h> 
#include <TLegend.h> 
#include <TString.h> 
#include <TGaxis.h> 
#include <TCanvas.h> 
#include <TFile.h>
#include <TSystem.h>
#include <TStopwatch.h>

#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnUserCovariance.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnApplication.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnHesse.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnPrint.h>

#include "FitterKernelMinuit2.h" 
#include "SuperkObservation.h" 
#include "SuperkExpectation.h" 
#include "FitInputs.h" 
#include "FitParamList.h" 
#include "Utilities.h"
#include "Conventions.h"

using std::cerr; 
using std::cout; 
using std::endl; 
using std::ostringstream;

using namespace t2kana; 
using namespace t2kana::conv;
   
//____________________________________________________________________________ 
FitterKernelMinuit2::FitterKernelMinuit2() :
FitterABC()
{
  this->Init();
}
//____________________________________________________________________________ 
FitterKernelMinuit2::FitterKernelMinuit2(FitInputs * fit_inputs, Option_t * option) :
FitterABC()
{
  this->Init(fit_inputs, option);
}
//____________________________________________________________________________ 
FitterKernelMinuit2::~FitterKernelMinuit2() 
{
  if(fResults)    { delete fResults;    }
}
//____________________________________________________________________________ 
void FitterKernelMinuit2::Init(void) 
{
  fFitInputs      = 0;
  fOption         = "";
  fNFreeFitParams = 0;
  fNFitBins       = 0;
  fResults        = 0;
}
//____________________________________________________________________________ 
void FitterKernelMinuit2::Init(FitInputs * fit_inputs, Option_t * option) 
{
#ifdef __show_fit_init_info__
  cout << "FitterKernelMinuit2/Init >> Fitter initialization..." << endl;
#endif

  this->Init();

  fResults = new TObjArray();

  fFitInputs = fit_inputs;
  if(!fFitInputs) {
     cerr << "FitterKernelMinuit2/Init >> ** Null fit inputs!" << endl;
     exit(1);
  }

  FitParamList & fit_params = fFitInputs->FitParams();
  int param_list_size = fit_params.Size();

#ifdef __show_fit_init_info__
  cout << "FitterKernelMinuit2/Init >> Fit parameter list size = " 
       << param_list_size << endl;
#endif
  
  // set params
  fit_params.SetTo(fUserParams);

  cout << fUserParams << endl;

#ifdef __show_fit_init_info__
  cout << "FitterKernelMinuit2/Init >> Constructing covariance matrix " << endl;
#endif


#ifdef __show_fit_init_debug__
  cout << "FitterKernelMinuit2/Init >> Covariance matrix:" << endl;
  //covariance_matrix.Print();
#endif

  // calculate degrees of freedom

  fNFitBins       = kNErecoBins1Rmu;
  fNFreeFitParams = fit_params.GetNumberFreeParameters();

#ifdef __show_fit_init_info__
  cout << "FitterKernelMinuit2/Init >> NDoF = " << this->NDOF()
       << " (Number of bins = " << this->NBinsInFit()
       << ", Number of free parameters = " << this->NFitParams() << ")"
       << endl;
#endif
}
//____________________________________________________________________________ 
bool FitterKernelMinuit2::DoTheFit(void) 
{
#ifdef __show_fit_cpu_time__
  TStopwatch w;
  w.Start();
#endif

  // Run the fitter once - we have the option here to check the status code
  // and if the fitter failed, re-seed the fitter and try again (not needed
  // so far).
  bool ok = this->RunOnce();

#ifdef __show_fit_cpu_time__
  w.Stop();
  cout << "FitterKernelMinuit2/DoTheFit >> Fit performed in " 
       << w.CpuTime() << " sec." << endl;
#endif

  return ok;
}
//____________________________________________________________________________ 
bool FitterKernelMinuit2::RunOnce(void) 
{
  // a few arrays to store parameters
  unsigned int nfreeparams = (unsigned int)fFitInputs->FitParams().GetNumberFreeParameters();
  unsigned int nparams     = (unsigned int)fFitInputs->FitParams().Size();
  vector<double> bestfit_free;
  vector<double> err_migrad_cov;
  vector<double> err_migrad_err;
#ifdef __get_hesse_errors__
  vector<double> err_hesse_cov;
  vector<double> err_hesse_err;
#endif
#ifdef __get_minos_errors__
  vector<double> err_minosplus;
  vector<double> err_minosminus;
#endif

  // Set fitter strategy
  // In general, low values of level mean fewer function calls and high values mean more reliable minimization.
  // Currently allowed values are 0 (low), 1 (default), and 2 (high).
  ROOT::Minuit2::MnStrategy strategy(1);

  // Get the fitter & pass it the inputs
#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Getting the fitter & giving it inputs (data & pdf)..." << endl;
#endif

  ROOT::Minuit2::MnMigrad migrad(*fFitInputs, fUserParams);
  
  // Run MINUIT's minimization step
#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Running MIGRAD..." << endl;
#endif

  // max number of function calls
  // if zero, max_calls = 200 + 100 * npar + 5 * npar**2
  int max_calls = 0;

  // value used for terminating iteration procedure.
  // For example, MIGRAD will stop iterating when edm (expected distance from minimum) will be:
  // edm < tolerance * 10**-3
  double tolerance = 0.1;

  //Minimise!
  ROOT::Minuit2::FunctionMinimum bestfit = migrad(max_calls, tolerance);

#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> MIGRAD minimisation complete" << endl;
#endif

  // Done! Copy current --> best-fit values
  fFitInputs->FitParams().CopyCurrToBestFitValues();

  //get the best fit vector, for use in this function
  std::vector<ROOT::Minuit2::MinuitParameter> bestfit_params = bestfit.UserParameters().Parameters();
  for(unsigned int ip = 0; ip < bestfit_params.size(); ip++) {
    if(!bestfit_params[ip].IsFixed()) {
      bestfit_free.push_back(bestfit_params[ip].Value());
    }
#ifdef __show_fit_steps_debug__
    cout << "FitterKernelMinuit2/RunOnce >> Bestfit for par#" << ip << " = " << bestfit_params[ip].Value() << endl;
#endif
  }

#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Getting fit success status..." << endl;
#endif
  //Check if fit succeeded
  if(!bestfit.IsValid()){
    std::cout << "FitterKernelMinuit2/RunOnce >> Fit Failed" << std::endl;
    return false;
  }
  else {
    std::cout << "FitterKernelMinuit2/RunOnce >> Fit succeeded" << std::endl;
    fStatusCode = 3;
  }

  // Get errors
#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Getting Migrad errors..." << endl;
#endif

  // migrad covariance
  ROOT::Minuit2::MnUserCovariance cov_migrad(bestfit.UserCovariance());
#ifdef __show_fit_steps_debug__
  cout << cov_migrad << endl;
#endif
  for(unsigned int ip = 0; ip < nfreeparams; ip++) {
    double tmp = cov_migrad(ip,ip);
    err_migrad_cov.push_back(tmp);
#ifdef __show_fit_steps_debug__
    cout << "FitterKernelMinuit2/RunOnce >> Migrad covariance for par#" << ip << " = " << bestfit_free[ip] << " +/- " << err_migrad_cov[ip] << endl;
#endif
  }

  // migrad error
  err_migrad_err =  bestfit.UserParameters().Errors(); 
  for(unsigned int ip = 0; ip < nparams; ip++) {
#ifdef __show_fit_steps_debug__
    cout << "FitterKernelMinuit2/RunOnce >> Migrad error for par#" << ip << " = " << bestfit_free[ip] << " +/- " << err_migrad_err[ip] << endl;
#endif
  }
  fFitInputs->FitParams().SetMigradErrors(err_migrad_err);

#ifdef __get_hesse_errors__
  //Get Hesse errors
#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Getting Hesse errors..." << endl;
#endif
  // argument is strategy
  ROOT::Minuit2::MnHesse hesse(2);
  // perform hesse error calculation
  hesse(*fFitInputs, bestfit);

  // hesse covariance
  ROOT::Minuit2::MnUserCovariance cov_hesse(bestfit.UserCovariance());
#ifdef __show_fit_steps_debug__
  cout << cov_hesse << endl;
#endif
  for(unsigned int ip = 0; ip < nfreeparams; ip++) {
    err_hesse_cov.push_back(cov_hesse(ip,ip));
#ifdef __show_fit_steps_debug__
    cout << "FitterKernelMinuit2/RunOnce >> Hesse covariance for par#" << ip << " = " << bestfit_free[ip] << " +/- " << err_hesse_cov[ip] << endl;
#endif
  }

  // hesse error
  err_hesse_err =  bestfit.UserParameters().Errors(); 
  fFitInputs->FitParams().SetHesseErrors(err_hesse_err);
  for(unsigned int ip = 0; ip < nfreeparams; ip++) {
#ifdef __show_fit_steps_debug__
    cout << "FitterKernelMinuit2/RunOnce >> Hesse error for par#" << ip << " = " << bestfit_free[ip] << " +/- " << err_hesse_err[ip] << endl;
#endif
  }
#endif

#ifdef __get_minos_errors__
#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Getting Minos errors..." << endl;
#endif

  //Get Minos Errors
  ROOT::Minuit2::MnMinos minos(*fFitInputs,bestfit);
  for(unsigned int ip = 0; ip < nparams; ip++) {
    if(!bestfit_params[ip].IsFixed())
      continue;
    double up = minos.Upper(ip);
    double lo = minos.Lower(ip);
    err_minosminus.push_back(pow(lo,2));
    err_minosplus. push_back (pow(up,2));
#ifdef __show_fit_steps_debug__
    cout << "FitterKernelMinuit2/RunOnce >> par#" << ip << " = " << bestfit_free[ip] << " +" << up << " -" << lo << " (minos)" << endl;
#endif
  }
  fFitInputs->FitParams().SetMinosErrors(err_minosplus, err_minosminus);
#endif

  // Print fit result
#ifdef __show_fit_result__
#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Printing fit results..." << endl;
#endif
  cout << "FitterKernelMinuit2/RunOnce >> fcn =  " << bestfit.Fval() << endl;
  cout << "FitterKernelMinuit2/RunOnce >> chi2 = " << fFitInputs->Chisquare() << endl;
  cout << bestfit << endl;
#endif

#ifdef __show_fit_steps__
  cout << "FitterKernelMinuit2/RunOnce >> Done" << endl;
#endif

  return (fStatusCode == 3);
}
//____________________________________________________________________________ 
