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

namespace t2kana {

class NuFlux: public TObject
{
public:
  NuFlux();   
 ~NuFlux();

  // Get beam oscillation analysis inputs
  double  SuperkNuMuFlux10aTo10dWght    (double Etrue) const;
  double  SuperkNuMuBarFlux10aTo10dWght (double Etrue) const;
  double  SuperkNuEFlux10aTo10dWght     (double Etrue) const;
  double  SuperkNuMuFlux                (double Etrue, bool wght_to_10d=true) const;
      
  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

  // Drawing option:
  // - "SKFlux/10a/NuMu"              : nominal SuperK numu flux (10a)
  // - "SKFlux/10d/NuMu"              : weighted SuperK numu flux (10d)
  // - "SKFlux/Both10a10d/NuMu"       : nominal (10a) and weighted SuperK numu flux (10d)
  // - "SKFlux/Weights/10aTo10d/NuMu" : 10a->10d numu flux weights
  // - "SKFlux/CovarianceMatrix/NuMu" : SK numu flux covariance matrix
  void Draw(TVirtualPad * c, Option_t * opt = "SKFlux/10d") const;

private:

  void Init    (void);
  void CleanUp (void);

  TH1D *        fSuperkNuMuFlux10a;            ///< 10a SuperK numu flux estimation (50 MeV binning)
  TH1D *        fSuperkNuMuFlux10aVarBin;      ///< 10a SuperK numu flux estimation (Beam WG var binning scheme)
  double        fSuperkNuMuFluxPOT;            ///< integrated POT for files used for numu flux simulation
  TH1D *        fSuperkNuMuWght10aTo10d;       ///< SuperK numu    flux weights (10a -> 10d)
  TH1D *        fSuperkNuMuBarWght10aTo10d;    ///< SuperK numubar flux weights (10a -> 10d)
  TH1D *        fSuperkNuEWght10aTo10d;        ///< SuperK nue     flux weights (10a -> 10d)
  double        fSuperkFluxNormTwkDial;        ///< flux normalization tweaking dial (1Re and 1Rmu tweaked coherently)
      
ClassDef(NuFlux, 4)
};

} //t2kana

#endif
