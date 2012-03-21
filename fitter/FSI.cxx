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
//#define __show_init_mesg_debug__
#define __show_fsi_cholesky_mesg_info__
#define __show_fsi_cholesky_mesg_debug__
#define __show_fsi_tweak_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <TMath.h>

#include "FSI.h"
#include "Utilities.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;

using namespace t2kana;
using namespace t2kana::conv;
using namespace t2kana::util;

ClassImp(FSI)

//________________________________________________________________________________________
FSI::FSI()
{
  this->Init();
}
//________________________________________________________________________________________
FSI::~FSI()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void FSI::Init(void)
{
  fCurrValues    = 0;
  fCovMtrx       = 0;
  fCovMtrxChoL   = 0;

/*
  for(int is=0; is<kns; is++) {
    for(int im=0; im<knm; im++) {
      fFracTun [is][im] = 0;
    }
  }
*/
}
//________________________________________________________________________________________
void FSI::CleanUp(void)
{
  if ( fCurrValues    ) delete fCurrValues;
  if ( fCovMtrx       ) delete fCovMtrx;
  if ( fCovMtrxChoL   ) delete fCovMtrxChoL;
/*
  for(int is=0; is<kns; is++) {
    for(int im=0; im<knm; im++) {
      if (fFracTun [is][im]) delete fFracTun [is][im];
    }
  }
*/
}
//________________________________________________________________________________________
bool FSI::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/FSI >> Reading FSI inputs from input XML file." << endl;
#endif

/*
  const int np = conv::kNumOfFSIParams;

  // Get FSI covariance matrix
  fCovMtrx = util::XML2TMatrixD(xml,root_node,"xsec/fsi/covariance_matrix");
  if(!fCovMtrx) {
    cerr << "Exp/FSI >> No FSI covariance matrix was found." << endl;
    exit(1);
  }
#ifdef __show_init_mesg_info__
  cout << "Exp/FSI >> Input FSI covariance matrix: " << endl;
  util::PrintCompact(*fCovMtrx,7,4);
#endif

  // check cov matrix dimension against no
  // TODO

  // Create/init current FSI params
  fCurrValues = new TVectorD(np);
  for(int ip = 0; ip < np; ip++) 
  {
     (*fCurrValues)[ip] = 1.;
  }

*/

/*
  //
  //
  //

  cout << "Exp/FSI >> Booking fractional MC templates..." << endl;

  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     const int      n_ereco_bins    = util::NumErecoBins (s);
     const double * ereco_bin_edges = util::ErecoBinEdges(s);
     for(int im = 0; im < knm; im++) {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       fFracTun[is][im] = new TH2D(
          Form("nsk_tun_%d_%d",is,im), util::AsString(m).c_str(),
               n_ereco_bins, ereco_bin_edges, kNEtrueBins, kEtrueMin, kEtrueMax);
       fFracTun[is][im]->SetDirectory(0);
       fFracTun[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fFracTun[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
     }
  }
*/

#ifdef __show_init_mesg_info__
  cout << "Exp/FSI >> Done reading neutrino cross-section inputs!" << endl;
#endif

  return true;  
}
//________________________________________________________________________________________
double FSI::SystWght(
  SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
  double tot_wght = 1.;




  return tot_wght;
}
//________________________________________________________________________________________
double FSI::CurrentValue(PhysParamId_t p) const
{
  return 0;
}
//________________________________________________________________________________________
void FSI::PrintSystematicVariations(void) const
{
  //TVectorD & fsi_current_values = *fCurrValues;
  //cout << "Exp/FSI >> FSI parameters: " << endl;
  //util::PrintCompact(fsi_current_values,12,6,"FSI params");
}
//________________________________________________________________________________________
void FSI::Tweak(
  bool       force_inputs,    /* force the input params or generate random ones internally  */
  TVectorD * input_fsi_params /* input vector of FSI params - Set to NULL to reset to nominal ones */
)
{
#ifdef __show_fsi_cholesky_mesg_debug__
  cout << "Exp/FSI/Tweak >> Tweaking FSI parameters" << endl;
#endif

  //
  // If this is the first time this function is called, then run the Cholesky
  // decomposition of the flux covariance matrix
  //

  if(fCovMtrxChoL == 0)
  {
#ifdef __show_fsi_cholesky_mesg_info__
    cout << "Exp/FSI/Tweak >> "
         << "Running Cholesky's decomposition method on the FSI covariance matrix." << endl;
#endif
    
    TMatrixD cov(this->CovMatrix());
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/FSI/Tweak >> FSI covariance matrix: " << endl;
    util::PrintCompact(cov,7,4);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fCovMtrxChoL = new TMatrixD(L);
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/FSI/Tweak >> Cholesky decomposition of FSI covariance matrix: " << endl;
    util::PrintCompact(*fCovMtrxChoL,7,4);
#endif
  }

  TVectorD & current_values = *fCurrValues;

  //
  // Force input FSI parameters
  //
  if(force_inputs)
  {
     // Input vector is NULL.
     // This resets the FSI params to their best-fit values.
     if(input_fsi_params==0) {
#ifdef __show_fsi_tweak_mesg_debug__
        cout << "Exp/FSI/Tweak >> Reset FSI params to nominal values." << endl;
#endif
        const int np = conv::kNumOfFSIParams;
        for(int ip = 0; ip < np; ip++) {
           current_values[ip] = 1.;
        }
        return;
     }
     else
     {
#ifdef __show_fsi_tweak_mesg_debug__
        cout << "Exp/FSI/Tweak >> Forcing input FSI parameters. " << endl;
#endif
        TVectorD & vec = *input_fsi_params;
        current_values = vec;
        return;
     }
  }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
  {   
#ifdef __show_fsi_tweak_mesg_debug__
     cout << "Exp/FSI/Tweak >> Generating random FSI params (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL; 
     TVectorD vec = util::CholeskyGenerateCorrelatedParamVariations(L);
     current_values = vec;
     return;
  }

}
//________________________________________________________________________________________
