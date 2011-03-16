//____________________________________________________________________________
/*!

\class    t2kana::NuFlux

\brief    Neutrino flux inputs to the oscillation analysis

\author   STFC, Rutherford Appleton Laboratory / Oxford T2K groups

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          Jelena Ilic         <Jelena.Ilic \at stfc.ac.uk>
          Nick Grant          <nick.grant \at stfc.ac.uk>
          Geoff Pearce        <geoff.pearce \at stfc.ac.uk>
          Alfons Weber        <alfons.weber \at stfc.ac.uk>

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_nu_flux_h_
#define _t2kana_nu_flux_h_

#include <TObject.h>
#include <TMatrixD.h>
#include <TMatrixTSym.h>
#include <TXMLEngine.h>

class TH1D;
class TVirtualPad;

class HadErrCalc;

namespace t2kana {

class NuFlux: public TObject
{
public:
  NuFlux();   
 ~NuFlux();

  // Get beam oscillation analysis inputs
  double              SuperkNuMuFlux10aTo10dWght    (double Etrue) const;
  double              SuperkNuMuBarFlux10aTo10dWght (double Etrue) const;
  double              SuperkNuEFlux10aTo10dWght     (double Etrue) const;
  double              SuperkNuMuFlux                (double Etrue, bool wght_to_10d=true) const;
  double              SuperkNuMuFluxPOT             (void) const;
  const TMatrixD &    SuperkNuMuFluxCov             (void) const; ///< Superk flux covariance matrix
  double              SuperkNuMuFluxShapeWeight     (double Etrue);
  double              SuperkFluxNormErr             (TH1D * hnumu1Rmu, TH1D * hnumubar1Rmu, TH1D * hnue1Rm,
                                                     TH1D * hnumuND, TH1D * hnumubarND, TH1D * hnueND) const; ///< normalization err
  double              SuperkFluxNorm                (TH1D * hnumu1Rmu, TH1D * hnumubar1Rmu, TH1D * hnue1Rm,
                                                     TH1D * hnumuND, TH1D * hnumubarND, TH1D * hnueND) const; ///< normalization
  bool                UsingKuboSansCode             (void) { return (fHEC!=0); }
  bool                NormTweaked                   (void) { return (fSuperkNuMuFluxNormTwkDial!=0); }
  bool                ShapeTweaked                  (void) { return fEnableFluxShapeTweaking; };
   
  // Drawing option:
  // - "SKFlux/10a/NuMu"              : nominal SuperK numu flux (10a)
  // - "SKFlux/10d/NuMu"              : weighted SuperK numu flux (10d)
  // - "SKFlux/Both10a10d/NuMu"       : nominal (10a) and weighted SuperK numu flux (10d)
  // - "SKFlux/Weights/10aTo10d/NuMu" : 10a->10d numu flux weights
  // - "SKFlux/CovarianceMatrix/NuMu" : SK numu flux covariance matrix
  void Draw(TVirtualPad * c, Option_t * opt = "SKFlux/10d") const;
   
  // Systematic variations   
  void TweakNorm  (double twkdial) { fSuperkNuMuFluxNormTwkDial = twkdial; }
  void TweakShape (bool enable_tweaking=true);
   
  // Validation-mode only systematic variations (tweak flux in a single etrue bin; Bin numbers inn Beam WG's variable binning scheme)
  void TweakShapeSetValidationMode         (bool mode)               { fTweakShapeVldMode = mode; }
  void TweakShapeValidationSetBinVariation (int bin, double twkdial) { fTweakShapeVldModeBin = bin; fTweakShapeVldModeDial = twkdial; }
   
  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

private:

  void Init    (void);
  void CleanUp (void);

  TH1D *        fSuperkNuMuFlux10a;            ///< 10a SuperK numu flux estimation (50 MeV binning)
  TH1D *        fSuperkNuMuFlux10aVarBin;      ///< 10a SuperK numu flux estimation (beam group var binning scheme)
  double        fSuperkNuMuFluxPOT;            ///< integrated POT for numu flux simulation
  TMatrixD *    fSuperkNuMuFluxCov;            ///< SuperK numu flux estimation
  TH1D *        fSuperkNuMuWght10aTo10d;       ///< SuperK numu    flux weights (10a -> 10d)
  TH1D *        fSuperkNuMuBarWght10aTo10d;    ///< SuperK numubar flux weights (10a -> 10d)
  TH1D *        fSuperkNuEWght10aTo10d;        ///< SuperK nue     flux weights (10a -> 10d)
  double        fSuperkNuMuFluxNormTwkDial;    ///<
  bool          fSuperkNuMuFluxHaveCholeskyL;  ///<
  TMatrixD *    fSuperkNuMuFluxCovCholeskyL;   ///<
  TH1D *        fSuperkNuMuFluxTweaked;        ///< SuperK numu flux estimation (tweaked / 50 MeV binning)
  bool          fEnableFluxShapeTweaking;      ///<
  HadErrCalc *  fHEC;                          ///< Kubo-san / BEAM WG norm err calc code
  double        fSuperkNormRelErr_simple;      ///< a fixed norm rel error (to use when not using Kubo-san's code)
      
  bool   fTweakShapeVldMode;     ///< validation-mode systematic variation enabled
  int    fTweakShapeVldModeBin;  ///< Etrue variable size bin tweaked
  double fTweakShapeVldModeDial; ///< flux bin tweaking dial	

ClassDef(NuFlux, 2)
};

} //t2kana

#endif
