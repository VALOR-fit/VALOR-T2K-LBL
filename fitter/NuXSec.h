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

using t2kana::conv::knwspl;

namespace t2kana {

class BANFF;

class NuXSec: public TObject
{
public:

  NuXSec();
 ~NuXSec();

  // Reweight to BANFF best-fit params
  double SuperkXSecFinTuningWght (conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;

  // Calculate weight in response to tweaked systematics relative to the BANFF best-fit (not the nominal MC)
  double SuperkXSecSystWght (conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;

  // Tweak 
  void TweakFSI (bool force_inputs = false, TVectorD * input_vec = 0);

  // Get current value of systematic params handled by this class
  double CurrentValue (conv::PhysParamId_t p) const;

  // Print
  void PrintSystematicVariations (void) const;

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

  //
  // Access data held by this object
  //

  // Get input FSI covariance matrix and its Cholesky decomposition
  const TMatrixD & FSICovMatrix         (void) const { return *fFSICovMtrx;       }
  const TMatrixD & FSICovMatrixCholesky (void) const { return *fFSICovMtrxChoL;   }

  // Get weight spline for each sample / mode / parameters /Ereco bin / Etrue bin
  TSpline3 * FindWghtSpline(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;

  //
  // Set interdepencies with other classes
  //

  // BANFF fit info
  void TellMeAboutTheBANFFfit(BANFF * banff);
   
private:

  void Init    (void);
  void CleanUp (void);

  int  WghtSplineIdx        (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  bool ApplySyst            (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  bool ApplySystToSKMode    (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m) const;
  bool ApplySystToEnergyBin (conv::PhysParamId_t p, double Ereco, double Etrue) const;

  //
  // Private data members
  //

  BANFF *     fBANFF;               ///< BANFF fit info
  TVectorD *  fFSICurrentValues;    ///< Current values of FSI params
  TMatrixD *  fFSICovMtrx;          ///< FSI covariance matrix
  TMatrixD *  fFSICovMtrxChoL;      ///< Triangular matrix from Cholesky decomp. of FSI covar. matrix
  TSpline3 *  fXSecWghtSpl[knwspl]; ///< Weight splines. Index is function of sample, mode, systematic param, Ereco bin, Etrue bin

ClassDef(NuXSec, 4);
};

} // t2kana

#endif
