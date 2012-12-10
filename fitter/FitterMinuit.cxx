//____________________________________________________________________________ 
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/ 
//____________________________________________________________________________

//#define __show_fit_init_info__ 
#define __show_fit_result__ 
//#define __show_fit_steps__ 
#define __show_fit_cpu_time__

#include <cassert> 
#include <cstdlib> 
#include <iostream> 
#include <sstream>

#include <TMath.h> 
#include <TRandom3.h> 
#include <TVirtualFitter.h> 
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

#include "FitterMinuit.h" 
#include "SuperkObservation.h" 
#include "SuperkExpectation.h" 
#include "FitInputs.h" 
#include "FitParamList.h" 
#include "Utilities.h"
#include "Conventions.h"
#include "Opt.h"

using std::cerr; 
using std::cout; 
using std::endl; 
using std::ostringstream;

using namespace t2kana; 
using namespace t2kana::conv;
   
//____________________________________________________________________________ 
FitterMinuit::FitterMinuit() :
FitterABC()
{
  this->Init();
}
//____________________________________________________________________________ 
FitterMinuit::FitterMinuit(FitInputs * fit_inputs, Option_t * option) :
FitterABC()
{
  this->Init(fit_inputs, option);
}
//____________________________________________________________________________ 
FitterMinuit::~FitterMinuit() 
{

}
//____________________________________________________________________________ 
void FitterMinuit::Init(void) 
{
  fFitInputs      = 0;
  fOption         = "";
  fNFreeFitParams = 0;
  fNFitBins       = 0;
  fGetHesseErrors = true;
}
//____________________________________________________________________________ 
void FitterMinuit::Init(FitInputs * fit_inputs, Option_t * option) 
{
#ifdef __show_fit_init_info__
  cout << "FitterMinuit/Init >> Fitter initialization..." << endl;
#endif

  this->Init();

  fFitInputs = fit_inputs;
  if(!fFitInputs) {
     cerr << "FitterMinuit/Init >> ** Null fit inputs!" << endl;
     exit(1);
  }

  FitParamList & fit_params = fFitInputs->FitParams();
  int param_list_size = fit_params.Size();

#ifdef __show_fit_init_info__
  cout << "FitterMinuit/Init >> Fit parameter list size = " 
       << param_list_size << endl;
#endif

  fOption = option;

  if(fOption.Contains("minuit")) {
    // get fitter
    TVirtualFitter::SetDefaultFitter("Minuit");
    TVirtualFitter * fitter = TVirtualFitter::Fitter(0,param_list_size);

    // set print level
    Double_t arglist[5];
    arglist[0] = -1;
    fitter->ExecuteCommand("SET PRINT",arglist,1);

    // set params
    fit_params.SetTo(fitter);

    // set minimization function
    fitter->SetFCN(FitFunc);

    /*
    // free parameters
    arglist[0] = 4; //deltacp
    fitter->ExecuteCommand("SET LIMits",arglist, 1);
    */
  }

  // calculate degrees of freedom

  fNFitBins       = kNErecoBins1Rmu; // HERE, FIX
  fNFreeFitParams = fit_params.GetNumberFreeParameters();

#ifdef __show_fit_init_info__
  cout << "FitterMinuit/Init >> NDoF = " << this->NDOF()
       << " (Number of bins = " << this->NBinsInFit()
       << ", Number of free parameters = " << this->NFitParams() << ")"
       << endl;
#endif
}
//____________________________________________________________________________ 
bool FitterMinuit::DoTheFit(void) 
{
#ifdef __show_fit_cpu_time__
  TStopwatch w;
  w.Start();
#endif
  fFitTime = 0;

  // Run the fitter once - we have the option here to check the status code
  // and if the fitter failed, re-seed the fitter and try again
  bool ok = this->RunOnce(false);

  // 4: migrad not converged
  int ntries = 0;
  while(fMigradStatusCode == 4 && ntries < 10) {
    cout << "Original fit failed; reseeding & trying again. Attempt: " << ntries << endl;
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    TRandom3 rand(0);
    switch(omp) {
    case(conv::kNuMixSinsq2Theta):
      fFitInputs->FitParams().SetInitValue(conv::kPhypMixing23, rand.Uniform(0.80,0.98));
      break;
    case(conv::kNuMixSinsqTheta):
      fFitInputs->FitParams().SetInitValue(conv::kPhypMixing23, rand.Uniform(0.35,0.55));
      break;
    case(conv::kNuMixUndefined):
      cerr << "Mixing parameter type undefinded. Exiting..." << endl;
      exit(1);
      break;
    }
    fFitInputs->FitParams().SetInitValue(conv::kPhypDm23sq, rand.Uniform(0.0020,0.0030));
    this->RunOnce(false);
    ntries++;
  }

  if(fMigradStatusCode == 4) {
    cout << "Fits still failed; trying SIMPLEX" << endl;
    this->RunOnce(true);
  }

#ifdef __show_fit_cpu_time__
  w.Stop();
  cout << "FitterMinuit/DoTheFit >> Fit performed in " 
       << w.CpuTime() << " sec." << endl;
  fFitTime = w.CpuTime();
#endif

  return ok;
}
//____________________________________________________________________________ 
bool FitterMinuit::RunOnce(bool use_simplex) 
{
  fHesseStatusCode = -1;

  // Get the fitter
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Getting the fitter..." << endl;
#endif
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();

  // Pass fit inputs (data, pdf) to MINUIT
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Passing fit inputs (data & pdf)..." << endl;
#endif
  gMinuit->SetObjectFit(fFitInputs);

  // Run MINUIT's minimization step
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Running MIGRAD..." << endl;
#endif
  Double_t arglist[100];
  arglist[0] = 50000; // num of func calls
  arglist[1] = 0.01;  // tolerance
  if(!use_simplex)
    fMigradStatusCode = fitter->ExecuteCommand("MIGRAD",arglist, 2);
  //SIMPLEX
  else {
    arglist[1] = 0.1; // higher tolerance
    fMigradStatusCode = fitter->ExecuteCommand("SIMPLEX",arglist, 2);
  }

  //return on not converged
  if(fMigradStatusCode == 4 && !use_simplex)
    return false;
  else if(fMigradStatusCode == 4 && use_simplex) {
    //copy current to best fit, to stop horrible crashes
    fFitInputs->FitParams().CopyCurrToBestFitValues();
    return false;
  }

  // Print fit result
#ifdef __show_fit_result__
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Printing MIGRAD fit results..." << endl;
#endif
  Double_t aminn = 0;
  fitter->PrintResults(3,aminn);
#endif

  // Get the error matrix status
  double ha(0.0);     // dummy vars
  double edm, errdef; // dummy vars
  int nvpar, nparx;   // dummy vars
  int errMatrixStatus = fitter->GetStats(ha, edm, errdef, nvpar, nparx);
  fStatusCode = errMatrixStatus;

  //Use Hesse to make the errors better
  if(fGetHesseErrors || fStatusCode != 3)
    fHesseStatusCode = fitter->ExecuteCommand("HESSE",arglist, 1);

  // Done! Copy current --> best-fit values
  fFitInputs->FitParams().CopyCurrToBestFitValues();

  // Get error matrix status
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Getting error matrix status..." << endl;
#endif
  errMatrixStatus = fitter->GetStats(ha, edm, errdef, nvpar, nparx);
  fStatusCode = errMatrixStatus;

  // Get errors
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Getting errors..." << endl;
#endif
  const int nparamsmax = knppmax;
  double errplus [nparamsmax] = {0.};
  double errminus[nparamsmax] = {0.};
  double errpara [nparamsmax] = {0.};
  for(int ip = 0; ip < nparamsmax; ip++) {
    double glob_corr; // dummy vars
    fitter->GetErrors(ip,errplus[ip],errminus[ip],errpara[ip],glob_corr);
  }
  fFitInputs->FitParams().SetErrors(errplus, errminus, errpara);

  // Print fit result
#ifdef __show_fit_result__
#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Printing HESSE fit results..." << endl;
#endif
  if(fGetHesseErrors) {
    Double_t amin = 0;
    fitter->PrintResults(3,amin);
    //fitter->PrintResults(1,amin); // prints limits
  }
  cout << "FitterMinuit/RunOnce >> Final error matrix status = " 
       << fStatusCode << endl; 
#endif

#ifdef __show_fit_steps__
  cout << "FitterMinuit/RunOnce >> Done" << endl;
#endif

  return (fStatusCode == 3);
}
//____________________________________________________________________________ 
