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
#define __show_fit_result__ 
//#define __show_fit_steps__ 
#define __show_fit_cpu_time__

#include <cassert> 
#include <cstdlib> 
#include <iostream> 
#include <sstream>

#include <TMath.h> 
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

#include "FitterKernel.h" 
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
FitterKernel::FitterKernel() :
FitterABC()
{
  this->Init();
}
//____________________________________________________________________________ 
FitterKernel::FitterKernel(FitInputs * fit_inputs, Option_t * option) :
FitterABC()
{
  this->Init(fit_inputs, option);
}
//____________________________________________________________________________ 
FitterKernel::~FitterKernel() 
{
//  if(fResults) { delete fResults; }
}
//____________________________________________________________________________ 
void FitterKernel::Init(void) 
{
  fFitInputs      = 0;
  fOption         = "";
  fNFreeFitParams = 0;
  fNFitBins       = 0;
  fResults        = 0;
}
//____________________________________________________________________________ 
void FitterKernel::Init(FitInputs * fit_inputs, Option_t * option) 
{
#ifdef __show_fit_init_info__
  cout << "FitterKernel/Init >> Fitter initialization..." << endl;
#endif

  this->Init();

  fResults = new TObjArray();

  fFitInputs = fit_inputs;
  if(!fFitInputs) {
     cerr << "FitterKernel/Init >> ** Null fit inputs!" << endl;
     exit(1);
  }

  FitParamList & fit_params = fFitInputs->FitParams();
  int param_list_size = fit_params.Size();

#ifdef __show_fit_init_info__
  cout << "FitterKernel/Init >> Fit parameter list size = " 
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

  fNFitBins       = kNErecoBins1Rmu;
  fNFreeFitParams = fit_params.GetNumberFreeParameters();

#ifdef __show_fit_init_info__
  cout << "FitterKernel/Init >> NDoF = " << this->NDOF()
       << " (Number of bins = " << this->NBinsInFit()
       << ", Number of free parameters = " << this->NFitParams() << ")"
       << endl;
#endif
}
//____________________________________________________________________________ 
bool FitterKernel::DoTheFit(void) 
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
  cout << "FitterKernel/DoTheFit >> Fit performed in " 
       << w.CpuTime() << " sec." << endl;
#endif

  return ok;
}
//____________________________________________________________________________ 
bool FitterKernel::RunOnce(void) 
{
  // Get the fitter
#ifdef __show_fit_steps__
  cout << "FitterKernel/RunOnce >> Getting the fitter..." << endl;
#endif
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();

  // Pass fit inputs (data, pdf) to MINUIT
#ifdef __show_fit_steps__
  cout << "FitterKernel/RunOnce >> Passing fit inputs (data & pdf)..." << endl;
#endif
  gMinuit->SetObjectFit(fFitInputs);

  // Run MINUIT's minimization step
#ifdef __show_fit_steps__
  cout << "FitterKernel/RunOnce >> Running MIGRAD..." << endl;
#endif
  Double_t arglist[100];
  arglist[0] = 5000; // num of func calls
  arglist[1] = 0.01; // tolerance
  fitter->ExecuteCommand("MIGRAD",arglist, 2);

  // Print fit result
#ifdef __show_fit_result__
#ifdef __show_fit_steps__
  cout << "FitterKernel/RunOnce >> Printing MIGRAD fit results..." << endl;
#endif
  Double_t aminn = 0;
  fitter->PrintResults(3,aminn);
#endif

  //Use Hesse to make the errors better
  fitter->ExecuteCommand("HESSE",arglist, 1);

  // Done! Copy current --> best-fit values
  fFitInputs->FitParams().CopyCurrToBestFitValues();

  // Get error matrix status
#ifdef __show_fit_steps__
  cout << "FitterKernel/RunOnce >> Getting error matrix status..." << endl;
#endif
  double ha(0.0);     // dummy vars
  double edm, errdef; // dummy vars
  int nvpar, nparx;   // dummy vars
  int errMatrixStatus = fitter->GetStats(ha, edm, errdef, nvpar, nparx);
  fStatusCode = errMatrixStatus;

  // Get errors
#ifdef __show_fit_steps__
  cout << "FitterKernel/RunOnce >> Getting errors..." << endl;
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
  cout << "FitterKernel/RunOnce >> Printing HESSE fit results..." << endl;
#endif
  Double_t amin = 0;
  fitter->PrintResults(3,amin);
  //fitter->PrintResults(1,amin); // prints limits
  cout << "FitterKernel/RunOnce >> Final error matrix status = " 
       << fStatusCode << endl;
 
//  TString smpl = fFitInputs->FitParams().GetTheSample(); 
//  cout << "FitterKernel/RunOnce >> Final chisq = ";
//  if (smpl =="numu"){
//        cout<<fFitInputs->Chisquare("poisson",conv::kSKsFCFV1RingMuLike)<< endl;
//  }else if(smpl =="nue"){
//        cout<<fFitInputs->Chisquare("poisson",conv::kSKsFCFV1RingELike)<< endl;
//  }else{
//        cout<<fFitInputs->Chisquare() << endl;
//  }

#endif

#ifdef __show_fit_steps__
  cout << "FitterKernel/RunOnce >> Done" << endl;
#endif

  return (fStatusCode == 3);
}
//____________________________________________________________________________ 
