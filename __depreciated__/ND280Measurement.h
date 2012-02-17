//____________________________________________________________________________
/*!

\class    t2kana::ND280Measurement

\brief    Encapsulate the ND280 measurements

\author   STFC, Rutherford Appleton Laboratory / Oxford T2K groups

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          Jelena Ilic         <Jelena.Ilic \at stfc.ac.uk>
          Nick Grant          <nick.grant \at stfc.ac.uk>
          Geoff Pearce        <geoff.pearce \at stfc.ac.uk>
          Alfons Weber        <alfons.weber \at stfc.ac.uk>
   
\created  October 20, 2010
   
*/
//____________________________________________________________________________
    
#ifndef _t2kana_nd280_measurement_h_
#define _t2kana_nd280_measurement_h_

#include <TObject.h>
#include <TXMLEngine.h>

#include "Conventions.h"

class TH1D;
class TH1F;
class TVirtualPad;

using t2kana::conv::knmnd;

namespace t2kana {

class NuXSec;

class ND280Measurement: public TObject
{
public:
  ND280Measurement();
 ~ND280Measurement();

  // Get nd280 oscillation analysis inputs
  double CorrectionFactor (void) const;
  double NObsCCLike       (bool nominal=false) const; ///< observed number of CC-like events (by default report tweaked value, if tweaked)
  double NExpNCBkg        (bool nominal=false) const; ///< expected NC contamination of CC-like sample (by default report tweaked value, if tweaked)
  double NExpCC           (bool nominal=false) const; ///< expected number of CC events  (by default report tweaked value, if tweaked)
  double NumOfEvents      (t2kana::conv::ND280Mode_t mode, bool nominal=false) const;
  double NumOfEvents      (t2kana::conv::ND280Mode_t mode, int bin, bool nominal=false) const;
  double NumOfEvents      (t2kana::conv::ND280Mode_t mode, double Etrue, bool nominal=false) const;
  double EfficiencyRelErr (double Et, int sign) const; ///< relative uncertainty on efficiency

  void   UpdateEtrueSpectra (void);
  TH1D*  CurrNuMuEtrue      (void) { return fNuMuEt;    } ///< current (allowing for applied weights) Etrue distribution for numu
  TH1D*  CurrNuMuBarEtrue   (void) { return fNuMuBarEt; } ///<      -//-  for numubar
  TH1D*  CurrNuEEtrue       (void) { return fNuEEt;     } ///<      -//-  for nue

  // Draw inputs: 
  // "Saclay/EtrueSpectrum/NuMu/All"
  // "Saclay/EtrueSpectrum/NuMu/CC"
  // "Saclay/EtrueSpectrum/NuMu/NC"
  // "Saclay/EtrueSpectrum/NuMuBar/CC"
  // "Saclay/EtrueSpectrum/NuMuBar/NC"
  // "Saclay/EtrueSpectrum/NuMuBar/CC"
  // "Saclay/EtrueSpectrum/NuMuBar/NC"
  void Draw(TVirtualPad * c, Option_t * opt = "") const;
      
  // Systematic variations
  void TweakEfficiency (double twkdial);

  // Statistical fluctuations
  void FluctuateNObsCCLike(double twkdial);

  // Print systematic variations
  void PrintSystematicVariations(void) const;

  //
  void TellMeAboutCrossSections(NuXSec * geninp) { fGenInputs = geninp; }
   
  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);
   
private:

  void   Init    (void);
  void   CleanUp (void);

  double  fPOTData;                ///< number of POT (data)
  double  fNumObsCCLike;           ///< Observed number of CC-like events
  double  fNumExpNCBkgDeflt;       ///<
  double  fNumNExpCCDeflt;         ///<
  double  fPOTMC;                  ///< number of POT (MC)
  TH1F *  fNumSelEt[knmnd];        ///< number of events passing cuts as function of Etrue for various modes
  double  fSelEffTwkDial;          ///< selection efficiency tweaking dial
  double  fSelEffRelErrPlus;       ///< relative error on selection efficiency
  double  fSelEffRelErrMinus;      ///< relative error on selection efficiency
  double  fNObsCCLikeTwkDial;      ///< tweak observed number of CC-like events within statistical error
  double  fRelStatErr;             ///< relative statistical error on observed number of CC-like events

  TH1D *  fNuMuEt;                 ///< fEtrue distribution for numu    - needed for beam norm err calculation
  TH1D *  fNuMuBarEt;              ///< fEtrue distribution for numubar - needed for beam norm err calculation
  TH1D *  fNuEEt;                  ///< fEtrue distribution for nue     - needed for beam norm err calculation

  double  fNomCorr;                ///< nominal correction factor (without efficiency, cross section tweaks)
   
  NuXSec * fGenInputs;

ClassDef(ND280Measurement, 2)
};

} // t2kana

#endif
