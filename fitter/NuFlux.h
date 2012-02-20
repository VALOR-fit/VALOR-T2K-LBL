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
#include <TH1D.h>

class TVirtualPad;

namespace t2kana {

class BANFF;

class NuFlux: public TObject
{
public:
  NuFlux();   
 ~NuFlux();

  // POT used for flux simulation
  double SuperkFluxPOT (void) const;

  // Tuned (default) or nominal SK flux
  double SuperkFlux (int neutrino_pdg, double Etrue, bool tuned=true, bool tweaked=true) const;

  // Weights for nominal SK flux -> tuned SK flux
  double SuperkFluxTuneWght (int neutrino_pdg, double Etrue) const; 

  // SK flux weights (ratios of tuned SK fluxes) as a result of BANFF tweaks
  double SuperkFluxWght (int neutrino_pdg, double Etrue) const; 

  // Access nominal flux and weight histograms held in this object
  TH1D * SuperkFluxSpectrum    (int neutrino_pdg, bool tuned=true, bool tweaked=true) const;
  TH1D * SuperkFluxTuningWghts (int neutrino_pdg) const;

  // Read configuration and data from input XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);

  // BANFF fit info
  void TellMeAboutTheBANFFfit(BANFF * banff);

  // Drawing options:
  // - "Flux/SK/NuMu/TuneWghts"    :
  // - "Flux/SK/NuMuBar/TuneWghts" :
  // - "Flux/SK/NuE/TuneWghts"     :
  // - "Flux/SK/NuEBar/TuneWghts"  :
  // - "Flux/SK/NuMu"              :
  // - "Flux/SK/NuMuBar"           :
  // - "Flux/SK/NuE"               :
  // - "Flux/SK/NuEBar"            :
  void Draw(TVirtualPad * c, Option_t * opt) const;

private:

  void Init               (void);
  void CleanUp            (void);
  void BookAndInitNomFlux (void);
  void BookAndInitTwkFlux (void);
  void CalcTwkFlux        (void) const; ///< Calculate tweaked flux as a result of BANFF param variations

  //
  // Note :
  //  - nominal flux : refers to the Superk flux used for MC production
  //  - tuned flux   : nominal flux + flux tuning weights
  //  - tweaked flux : tuned flux + BANFF fit parameter variations
  //
  TH1D *           fSuperkNuMuFluxNom;            ///< nominal SuperK numu    flux (50 MeV binning)
  TH1D *           fSuperkNuMuBarFluxNom;         ///< nominal SuperK numubar flux (50 MeV binning)
  TH1D *           fSuperkNuEFluxNom;             ///< nominal SuperK nue     flux (50 MeV binning)
  TH1D *           fSuperkNuEBarFluxNom;          ///< nominal SuperK nuebar  flux (50 MeV binning)
  TH1D *           fSuperkNuMuFluxTun;            ///< tuned   SuperK numu    flux (50 MeV binning)
  TH1D *           fSuperkNuMuBarFluxTun;         ///< tuned   SuperK numubar flux (50 MeV binning)
  TH1D *           fSuperkNuEFluxTun;             ///< tuned   SuperK nue     flux (50 MeV binning)
  TH1D *           fSuperkNuEBarFluxTun;          ///< tuned   SuperK nuebar  flux (50 MeV binning)
  mutable TH1D *   fSuperkNuMuFluxTwk;            ///< tweaked SuperK numu    flux (50 MeV binning)
  mutable TH1D *   fSuperkNuMuBarFluxTwk;         ///< tweaked SuperK numubar flux (50 MeV binning)
  mutable TH1D *   fSuperkNuEFluxTwk;             ///< tweaked SuperK nue     flux (50 MeV binning)
  mutable TH1D *   fSuperkNuEBarFluxTwk;          ///< tweaked SuperK nuebar  flux (50 MeV binning)
  TH1D *           fSuperkNuMuFluxTunWght;        ///< SuperK numu    flux tuning weights
  TH1D *           fSuperkNuMuBarFluxTunWght;     ///< SuperK numubar flux tuning weights
  TH1D *           fSuperkNuEFluxTunWght;         ///< SuperK nue     flux tuning weights
  TH1D *           fSuperkNuEBarFluxTunWght;      ///< SuperK nue     flux tuning weights
  BANFF *          fBANFF;                        ///< BANFF fit info
      
ClassDef(NuFlux, 5)
};

} //t2kana

#endif
