//____________________________________________________________________________
/*!

\class    t2kana::NuFlux

\brief    Neutrino flux inputs to the oscillation analysis

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

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
class TH2D;
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
  double*             SuperkFluxNormErr             (TH1D * hnumu1Rmu, TH1D * hnumubar1Rmu, TH1D * hnue1Rmu,
                                                     TH1D * hnumu1Re, TH1D * hnumubar1Re, TH1D * hnue1Re, TH1D * hnuesig1Re,
                                                     TH1D * hnumuND, TH1D * hnumubarND, TH1D * hnueND) const; ///< normalization err
  double*             SuperkFluxNorm                (TH1D * hnumu1Rmu, TH1D * hnumubar1Rmu, TH1D * hnue1Rmu,
                                                     TH1D * hnumu1Re, TH1D * hnumubar1Re, TH1D * hnue1Re, TH1D * hnuesig1Re,
                                                     TH1D * hnumuND, TH1D * hnumubarND, TH1D * hnueND) const; ///< normalization
  bool                UsingKuboSansCode             (void) { return (fHEC!=0); }
  bool                NormTweaked                   (void) { return (fSuperkFluxNormTwkDial!=0); }
  bool                ShapeTweaked                  (void) { return (fFluxShapeWasTweaked && !fIgnoreFluxShapeTweaks); }
   
  // Systematic variations (norm + shape)
  void TweakNorm  (double twkdial);
  void TweakShape (bool force_input_twkdials = false, TVectorD * input_twkv = 0, bool correlated = true, bool fractional = false);

  void IgnoreShapeTweaking(bool ignore) { fIgnoreFluxShapeTweaks = ignore; }
   
  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

  // Drawing option:
  // - "SKFlux/10a/NuMu"              : nominal SuperK numu flux (10a)
  // - "SKFlux/10d/NuMu"              : weighted SuperK numu flux (10d)
  // - "SKFlux/Both10a10d/NuMu"       : nominal (10a) and weighted SuperK numu flux (10d)
  // - "SKFlux/Weights/10aTo10d/NuMu" : 10a->10d numu flux weights
  // - "SKFlux/CovarianceMatrix/NuMu" : SK numu flux covariance matrix
  void Draw(TVirtualPad * c, Option_t * opt = "SKFlux/10d") const;

  // Print
  void PrintSystematicVariations (void) const;

private:

  void Init    (void);
  void CleanUp (void);

  void BuildNuMuFluxCovMatrix     (TH2D * fractional_covariance);
  void FillInTweakedNuMuFluxHisto (TVectorD * input_twkv, bool fractional = false);

  TH1D *        fSuperkNuMuFlux10a;            ///< 10a SuperK numu flux estimation (50 MeV binning)
  TH1D *        fSuperkNuMuFlux10aVarBin;      ///< 10a SuperK numu flux estimation (Beam WG var binning scheme)
  double        fSuperkNuMuFluxPOT;            ///< integrated POT for files used for numu flux simulation
  TMatrixD *    fSuperkNuMuFluxCov;            ///< SuperK numu flux covariance matrix (Beam WG var binning scheme)
  TH1D *        fSuperkNuMuWght10aTo10d;       ///< SuperK numu    flux weights (10a -> 10d)
  TH1D *        fSuperkNuMuBarWght10aTo10d;    ///< SuperK numubar flux weights (10a -> 10d)
  TH1D *        fSuperkNuEWght10aTo10d;        ///< SuperK nue     flux weights (10a -> 10d)
  double        fSuperkFluxNormTwkDial;        ///< flux normalization tweaking dial (1Re and 1Rmu tweaked coherently)
  bool          fSuperkNuMuFluxHaveCholeskyL;  ///<
  TMatrixD *    fSuperkNuMuFluxCovCholeskyL;   ///<
  TH1D *        fSuperkNuMuFluxTweaked;        ///< Tweaked SuperK numu flux estimation (50 MeV binning)
  bool          fFluxShapeWasTweaked;          ///< Was there a flux shape tweaking?
  bool          fIgnoreFluxShapeTweaks;        ///< Ignore flux shape tweaking?
  HadErrCalc *  fHEC;                          ///< Kubo-san's / Beam WG norm err calc code
  double        fSuperkNormRelErr_simple;      ///< a fixed norm rel error (to use when not using Kubo-san's code)
      
ClassDef(NuFlux, 3)
};

} //t2kana

#endif
