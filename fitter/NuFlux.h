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

class BANFF;

class NuFlux: public TObject
{
public:
  NuFlux();   
 ~NuFlux();

  //      
  double SuperkFluxPOT (void) const;

  // Weights for nominal SK flux -> tuned SK flux
  double  SuperkNuMuFluxTuneWght    (double Etrue) const; 
  double  SuperkNuMuBarFluxTuneWght (double Etrue) const; 
  double  SuperkNuEFluxTuneWght     (double Etrue) const; 

  // Tuned (default) or nominal SK flux
  double  SuperkNuMuFlux    (double Etrue, bool tuned=true) const;
  double  SuperkNuMuBarFlux (double Etrue, bool tuned=true) const;
  double  SuperkNuEFlux     (double Etrue, bool tuned=true) const;

  // SK flux weights (ratios of tuned SK fluxes) as a result of BANFF tweaks
  double  SuperkNuMuFluxWght    (double Etrue) const; 
  double  SuperkNuMuBarFluxWght (double Etrue) const; 
  double  SuperkNuEFluxWght     (double Etrue) const; 

  // Read inputs from XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

  // Drawing option:
  // - "SKFlux/10a/NuMu"              : nominal SuperK numu flux (10a)
  // - "SKFlux/10d/NuMu"              : weighted SuperK numu flux (10d)
  // - "SKFlux/Both10a10d/NuMu"       : nominal (10a) and weighted SuperK numu flux (10d)
  // - "SKFlux/Weights/10aTo10d/NuMu" : 10a->10d numu flux weights
  // - "SKFlux/CovarianceMatrix/NuMu" : SK numu flux covariance matrix
  void Draw(TVirtualPad * c, Option_t * opt = "SKFlux/10d") const;

  // 
  void TellMeAboutTheBANFFfit(BANFF * banff);

private:

  void Init    (void);
  void CleanUp (void);

  double Flux         (double Etrue, bool tuned, TH1D * nominal_flux, TH1D * tuning_weights) const;
  double FluxTuneWght (double Etrue, TH1D * tuning_weights) const;

  double   fSuperkFluxPOT;                ///< integrated POT for files used for flux simulation
  TH1D *   fSuperkNuMuFluxNom;            ///< nominal SuperK numu    flux (50 MeV binning)
  TH1D *   fSuperkNuMuBarFluxNom;         ///< nominal SuperK numubar flux (50 MeV binning)
  TH1D *   fSuperkNuEFluxNom;             ///< nominal SuperK nue     flux (50 MeV binning)
  TH1D *   fSuperkNuMuFluxTuneWght;       ///< SuperK numu    flux tuning weights (10a -> 10d)
  TH1D *   fSuperkNuMuBarFluxTuneWght;    ///< SuperK numubar flux tuning weights (10a -> 10d)
  TH1D *   fSuperkNuEFluxTuneWght;        ///< SuperK nue     flux tuning weights (10a -> 10d)
  BANFF *  fBANFF;                        ///< BANFF fit info
      
ClassDef(NuFlux, 5)
};

} //t2kana

#endif
