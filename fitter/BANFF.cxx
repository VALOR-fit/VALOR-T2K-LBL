//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis) 
 see code in the corresponding SVN branch.
 
 Changes since 2011_v0.9.0:

 @

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

#include "BANFF.h" 
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
void BANFF::Init(void)
{
  fBestFitValues = 0; 
  fCurrentValues = 0; 
  fCovMtrx       = 0; 
  fCovMtrxChoL   = 0; 
}
//________________________________________________________________________________________
void BANFF::CleanUp(void)
{
  if ( fBestFitValues ) delete fBestFitValues;
  if ( fCurrentValues ) delete fCurrentValues;
  if ( fCovMtrx       ) delete fCovMtrx;
  if ( fCovMtrxChoL   ) delete fCovMtrxChoL;
}
//________________________________________________________________________________________
bool BANFF::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
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

  // init current values to best-fit values
  fCurrentValues = new TVectorD(*fBestFitValues);

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
double BANFF::BestFitValue(PhysParamId_t p) const
{
  const TVectorD & bestfit = this->BestFitValues();
  int idx = this->Param2Idx(p);
  return bestfit[idx];
}
//________________________________________________________________________________________
double BANFF::CurrentValue(PhysParamId_t p) const
{
  const TVectorD & current_values = *fCurrentValues;
  int idx = this->Param2Idx(p);
  return current_values[idx];
}
//________________________________________________________________________________________
void BANFF::Tweak(
  bool       force_input_twkv,  /* force the input twk dials or generate twk dials internally */
  TVectorD * input_twkv,        /* input twk dials (if forced) - set 0 to reset */
  bool       correlated         /* twk dials are the correlated or uncorrelated ones? */
)    
{ 
#ifdef __show_cholesky_mesg_debug__
  cout << "Exp/BANFF/Tweak >> Tweaking BANFF fit parameters" << endl;
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

  //
  // Forcing null tweaking dials. This resets the BANFF params to their best-fit values.
  //

  if(force_input_twkv && input_twkv==0) {
#ifdef __show_tweak_mesg_debug__
    cout << "Exp/BANFF/Tweak >> Reset BANFF params to best-fit values." << endl;
#endif
     this->DoTweak(0);
     return;
  }

  //
  // Check whether the values of the tweaking dials are forced.
  // If the forced tweaking dials are supposed to be the correlated ones then we don't
  // need to use the Cholesky decomposition of the flux covariance matrix.
  //

  if(force_input_twkv && correlated) {
#ifdef __show_tweak_mesg_debug__
     cout << "Exp/BANFF/Tweak >> Forcing input twk dials (correlated)! " << endl;
#endif
     this->DoTweak(input_twkv);
     return;
  }

  //
  // Check whether the values of the tweaking dials are forced.
  // If the forced tweaking dials are supposed to be the un-correlated ones then use
  // the Cholesky decomposition of the flux covariance matrix to generate
  // correlated tweaking dials.
  //

  if(force_input_twkv && !correlated) {
#ifdef __show_tweak_mesg_debug__
     cout << "Exp/BANFF/Tweak >> Forcing input twk dials (uncorrelated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL;
     TVectorD twkv = util::CholeskyCalculateCorrelatedParamVariations(L,*input_twkv);
     this->DoTweak(&twkv);
     return;
  }

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  // Use Cholseky method to generate random fluxes
  //

  if(!force_input_twkv) {
#ifdef __show_tweak_mesg_debug__
    cout << "Exp/BANFF/Tweak >> Generating random twk dials (correlated)! " << endl;
#endif
    TMatrixD & L = *fCovMtrxChoL;
    TVectorD twkv = util::CholeskyGenerateCorrelatedParamVariations(L);
    this->DoTweak(&twkv);
    return;
  }

}
//________________________________________________________________________________________
void BANFF::PrintSystematicVariations(void) const
{
  cout << "BANFF/Tweak >> ... " << endl;
}
//________________________________________________________________________________________
int BANFF::Param2Idx(PhysParamId_t p) const
{
  int iparam = (int) p;
  int first  = (int) kNuisBANFF00;
  int last   = (int) kNuisBANFF36;
  assert(iparam >= first && iparam <= last);
  int idx = iparam-first;
  return idx;  
}
//________________________________________________________________________________________
void BANFF::DoTweak(TVectorD * input_twkv)
{
// Apply input_twkv to fBestFitValues and get fCurrentValues
//

  TVectorD & current_values = *fCurrentValues;
  TVectorD & bestfit_values = *fBestFitValues;
  TVectorD & twk_dials      = *input_twkv;

//
// TODO
// TODO
// TODO
//

}
//________________________________________________________________________________________
