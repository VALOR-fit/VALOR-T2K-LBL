//________________________________________________________________________________________
/*
 Changes since v1.0.0:

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__  
#define __show_init_mesg_debug__   
#define __show_cholesky_mesg_info__
//#define __show_cholesky_mesg_debug__
//#define __show_tweak_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TVectorD.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TSystem.h>

#include "NuFlux.h" 
#include "Conventions.h" 
#include "Utilities.h" 

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;
using namespace t2kana::util;

ClassImp(BANFF)

//________________________________________________________________________________________
BANFF::BANFF()
{
  this->Init();
}
//________________________________________________________________________________________
BANFF::~BANFF()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void NuFlux::Init(void)
{
  fBestFitValues = 0; 
  fCovMtrx       = 0; 
  fCovMtrxChoL   = 0; 
  fIsTweaked     = 0; 
}
//________________________________________________________________________________________
void BANFF::CleanUp(void)
{
  if ( fBestFitValues ) delete fBestFitValues;
  if ( fCovMtrx       ) delete fCovMtrx;
  if ( fCovMtrxChoL   ) delete fCovMtrxChoL;
}
//________________________________________________________________________________________
bool NuFlux::BANFF(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read BANFF inputs to oscillation analysis
//
#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Reading BANFF inputs from XML configuration file." << endl;
#endif

  // Get BANFF fit central values
  fBestFitValues = util::XML2TVectorD(xml,root_node,"banff/best_fit_values");
  if(!fBestFitValues) {
    cerr << "Exp/BANFF >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }

  // Get BANFF fit covariance matrix
  fCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/covariance_matrix");
  if(!fCovMtrx) {
    cerr << "Exp/BANFF >> No BANFF covariance matrix was found." << endl;
    exit(1);
  }

#ifdef __show_init_mesg_info__
  cout << "Exp/BANFF >> Done reading beam inputs!" << endl;
#endif

  return true;
}
//________________________________________________________________________________________
void BANFF::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  if(opt.Contains("??")) { 

    //
    // ...
    //

    c->Update();
  }
}
//________________________________________________________________________________________
void NuFlux::Tweak(
  bool       force_input_twkdials,  /* force the input twk dials or generate twk dials internally */
  TVectorD * input_twkv,            /* input twk dials (if forced) - set 0 to reset */
  bool       correlated,            /* twk dials are the correlated or uncorrelated ones? */
  bool       fractional             /* twk dials expressing a fractional or absolute flux change? */
)    
{ 
#ifdef __show_cholesky_mesg_debug__
  cout << "Exp/BANFF/Tweak >> Tweaking SuperK nu_mu flux shape" << endl;
#endif

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //
  if(fCovMtrxChoL == 0) 
  {
#ifdef __show_cholesky_mesg_info__  
    cout << "Exp/BANFF/Tweak >> "
         << "Running Cholesky's decomposition method on the BANFF covariance matrix." << endl;
#endif

    TMatrixD cov(this->CovMatrix());
#ifdef __show_shape_cholesky_mesg_debug__
    cout << "Exp/BANFF/Tweak >> BANFF covariance matrix: " << endl;
    cov.Print();
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fCovMtrxChoL = new TMatrixD(L);
#ifdef __show_cholesky_mesg_debug__
    cout << "Exp/BANFF/Tweak >> Cholesky decomposition of BANFF covariance matrix: " << endl;
    fCovMtrxChoL->Print();
#endif
  }

/*
  //
  // Forcing null tweaking dials.
  // Reset tweaked flux to nominal flux.
  //
  if(force_input_twkdials && input_twkv==0) {
#ifdef __show_tweak_mesg_debug__
    cout << "Exp/BANFF/Tweak >> Reset tweaked numu flux histogram to nominal values " << endl;
#endif
     this->FillInTweakedNuMuFluxHisto(0);
     return;
  }

  //
  // Check whether the values of the tweaking dials are forced.
  // If the forced tweaking dials are supposed to be the correlated ones then we don't
  // need to use the Cholesky decomposition of the flux covariance matrix.
  //
  if(force_input_twkdials && correlated) {
#ifdef __show_tweak_mesg_debug__
     cout << "Exp/BANFF/Tweak >> Forcing input twk dials (correlated)! " << endl;
#endif
     this->FillInTweakedNuMuFluxHisto(input_twkv,fractional);
     return;
  }

  //
  // Check whether the values of the tweaking dials are forced.
  // If the forced tweaking dials are supposed to be the un-correlated ones then use
  // the Cholesky decomposition of the flux covariance matrix to generate
  // correlated tweaking dials.
  //
  if(force_input_twkdials && !correlated) {
#ifdef __show_tweak_mesg_debug__
     cout << "Exp/BANFF/Tweak >> Forcing input twk dials (uncorrelated)! " << endl;
#endif
     TMatrixD & L = *fSuperkNuMuFluxCovCholeskyL;
     TVectorD twkv = util::CholeskyCalculateCorrelatedParamVariations(L,*input_twkv);
     this->FillInTweakedNuMuFluxHisto(&twkv,false);
     return;
  }

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  // Use Cholseky method to generate random fluxes
  //
  if(!force_input_twkdials) {
#ifdef __show_tweak_mesg_debug__
    cout << "Exp/BANFF/Tweak >> Generating random twk dials (correlated)! " << endl;
#endif
    TMatrixD & L = *fSuperkNuMuFluxCovCholeskyL;
    TVectorD twkv = util::CholeskyGenerateCorrelatedParamVariations(L);
    this->FillInTweakedNuMuFluxHisto(&twkv,false);
    return;
  }
*/

}
//________________________________________________________________________________________
void NuFlux::PrintSystematicVariations(void) const
{
  cout << "BANFF/Tweak >> ... " << endl;
}
//________________________________________________________________________________________
