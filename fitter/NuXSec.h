//____________________________________________________________________________
/*!
   
\class    t2kana::NuXSec
      
\brief    Neutrino generator / cross section info used in the oscillation fit
   
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________
   
#ifndef _t2kana_nu_xsec_h_
#define _t2kana_nu_xsec_h_

#include <TObject.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TSpline.h>
#include <TXMLEngine.h>

#include "Conventions.h"

using t2kana::conv::kns;
using t2kana::conv::knm;
using t2kana::conv::kNumOfXSecSystWithSpl;
using t2kana::conv::kXsNErecoBins;
using t2kana::conv::kXsNEtrueBins;

namespace t2kana {

class BANFF;

class NuXSec: public TObject
{
public:

  NuXSec();
 ~NuXSec();

  // Reweight to BANFF best-fit params
  double SuperkXSecFinTuningWght (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, double Ereco, double Etrue) const;

  // Calculate weight in response to tweaked systematics relative to the BANFF best-fit (not the nominal MC)
  double SuperkXSecSystWght (t2kana::conv::SuperkSample_t s, t2kana::conv::SuperkMode_t m, double Ereco, double Etre) const;

  //

  // Get current value of FSI params
  double FSICurrentValue (conv::PhysParamId_t p) const;
  // Get input FSI covariance matrix and its Cholesky decomposition
  const TMatrixD & FSICovMatrix         (void) const { return *fFSICovMtrx;       }
  const TMatrixD & FSICovMatrixCholesky (void) const { return *fFSICovMtrxChoL;   }

  // Tweak 
  void TweakFSI (bool force_inputs = false, TVectorD * input_vec = 0);

  // Print
  void PrintSystematicVariations (void) const;

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

  // BANFF fit info
  void TellMeAboutTheBANFFfit(BANFF * banff);
   
private:

  void Init    (void);
  void CleanUp (void);

  BANFF *     fBANFF;               ///< BANFF fit info
  TVectorD *  fFSICurrentValues;    ///< Current values of FSI params
  TMatrixD *  fFSICovMtrx;          ///< FSI covariance matrix
  TMatrixD *  fFSICovMtrxChoL;      ///< Triangular matrix from Cholesky decomp. of FSI covar. matrix

  TSpline3 *  fXSecWghtSpl[kns][knm][kNumOfXSecSystWithSpl][kXsNErecoBins][kXsNEtrueBins];

ClassDef(NuXSec, 4);
};

} // t2kana

#endif
