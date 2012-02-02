//____________________________________________________________________________ 
/*
 Combined fitter, combining Superk1RMuErecoFitter and Superk1ReNFitter
  These are chagnes in SuperkCombinedFitter.cxx:

 Changes since v1.0.0:

 @ 01/02/2012 - LE
   First added in v2.2.1
*/ 
//____________________________________________________________________________

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

#include "Superk1ReNFitter.h" 
#include "Superk1RmuErecoFitter.h" 
#include "SuperkCombinedFitter.h" 
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
SuperkCombinedFitter::SuperkCombinedFitter() :
FitterKernel() 
{
  fitter1Rmu = new Superk1RmuErecoFitter();
  fitter1Re  = new Superk1ReNFitter();
}
//____________________________________________________________________________ 
SuperkCombinedFitter::SuperkCombinedFitter(FitInputs * fi, Option_t * opt) :
FitterKernel(fi,opt)
{
  fitter1Rmu = new Superk1RmuErecoFitter(fi, opt);
  fitter1Re  = new Superk1ReNFitter(fi, opt);
}
//____________________________________________________________________________ 
SuperkCombinedFitter::~SuperkCombinedFitter() 
{

}
//____________________________________________________________________________ 
void SuperkCombinedFitter::Analyze(Option_t * option) 
{
  util::SetDefaultStyle();

  TString opt = option;
  cout << "SuperkCombinedFitter/Analyze >> Option: " << opt.Data() << endl;

  //
  // Print best-fit results used for generating the plots
  //
  fFitInputs->FitParams().ShowOnlyBestFitValuesInNextPrint();
  cout << "SuperkCombinedFitter/Analyze >> "
       << "Using the following best-fit results: " << fFitInputs->FitParams();

  //call each independent fitter to do the analysis
  fitter1Rmu->Analyze(option);
  fitter1Re->Analyze(option);
  
}
//____________________________________________________________________________ 
void SuperkCombinedFitter::Save(const char * filename, const char * opt_file) 
{
  Save(filename);
}
//____________________________________________________________________________ 
void SuperkCombinedFitter::Save(const char * filename) 
{
  util::SetDefaultStyle();

  cout << "SuperkCombinedFitter/Save >> Writing fit results in: " 
       << filename << endl;

  fitter1Rmu->Save(filename, "recreate");
  fitter1Re ->Save(filename, "update");
}
//____________________________________________________________________________ 
