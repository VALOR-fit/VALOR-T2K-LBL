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
#include <TXMLEngine.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TH2.h>
#include <TSpline.h>

#include "Conventions.h"

using t2kana::conv::knwspl;

namespace t2kana {

class BANFF;

class NuXSec: public TObject
{
public:

  NuXSec();
 ~NuXSec();

  //
  // Calculate weights
  //

  // Reweight to best-fit value of input BANFF param
  double FinTuningWght(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;

  // As above, but compute the combined weight for all BANFF cross-section parameters to the input sample/mode and 
  // energy bin. The input MC template is used for calculating how to share between the MC template bins weights computed
  // over different var.-sized Ereco,Etrue bins
  double FinTuningWght(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue, const TH2D * mctempl) const;

  // Calculate weight in response to tweaking the input parameter
  // For the BANFF parameters, the weight is relative to the best-fit (tuned MC) parameter, not the nominal MC one
  double SystWght(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;

  // As above, but compute the combined weight for all cross-section parameters.
  double SystWght(conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue, const TH2D * mctempl) const;

  //
  // Tweak systematic parameters held in this object.
  // Note, BANFF cross-section systematic parameters are tweaked via the BANFF object and they are 
  // correlated with the flux systematic params.
  //
  void TweakFSI (bool force_inputs = false, TVectorD * input_vec = 0);

  //
  // Access data held by this object
  //

  // Get current value of systematic params handled by this class
  double CurrentValue (conv::PhysParamId_t p) const;

  // Get input FSI covariance matrix and its Cholesky decomposition
  const TMatrixD & FSICovMatrix         (void) const { return *fFSICovMtrx;       }
  const TMatrixD & FSICovMatrixCholesky (void) const { return *fFSICovMtrxChoL;   }

  // Get weight spline for each sample / mode / parameters / Ereco bin / Etrue bin
  TSpline3 * FindWghtSpline(conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;

  // Print
  void PrintSystematicVariations (void) const;

  //
  // Configure
  //

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

  //
  // Set interdepencies with other classes
  //

  // BANFF fit info
  void TellMeAboutTheBANFFfit(BANFF * banff);
   
private:

  void Init    (void);
  void CleanUp (void);

  // Calculate weight spline index
  int  WghtSplineIdx        (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const;
  // Apply systematic to inmput sample/mode and energy bins?
  bool ApplySyst            (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue) const; 
  bool ApplySystToSKMode    (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m) const;
  bool ApplySystToEnergyBin (conv::PhysParamId_t p, double Ereco, double Etrue) const;
  //
  double FractionOfEvents   (conv::PhysParamId_t p, conv::SuperkSample_t s, conv::SuperkMode_t m, double Ereco, double Etrue, const TH2D * mctempl) const;

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
