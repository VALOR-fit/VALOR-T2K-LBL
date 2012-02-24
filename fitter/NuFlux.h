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

#include "Conventions.h"

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

  // SK flux (current, tuned/fin, tuned/ext, nominal - depending on input flux type)
  double SuperkFlux (int neutrino_pdg, double Etrue, conv::FluxType_t type = conv::kFlxCurrent) const;

  // SK flux tuning and systematic variation weights
  double SuperkFluxExtTuningWght  (int neutrino_pdg, double Etrue) const;  // nominal --> tuned (ext. data only)
  double SuperkFluxFinTuningWght  (int neutrino_pdg, double Etrue) const;  // tuned (ext. data only) --> tuned (ext. + T2K data)
  double SuperkFluxSystWght       (int neutrino_pdg, double Etrue) const;  // weight due to systematic variation of tuned (ext. + T2K data) flux

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

  // Initialize and clean-up
  void Init    (void);
  void CleanUp (void);

  // Book and initialize/compute spectra not loaded from config files
  void BookAndInitAuxSpectra(void);

  // Apply BANFF weights
  void ApplyBANFFWeights (bool rew_to_banff_bestfit) const; 

  // Access flux and weight histograms held in this object
  TH1D * SuperkFluxSpectrum       (int neutrino_pdg, conv::FluxType_t type = conv::kFlxCurrent) const;
  TH1D * SuperkFluxExtTuningWghts (int neutrino_pdg) const;

  //
  // Note on fluxes (see also Conventions.h):
  //  - nominal          : flux used for MC production
  //  - tuned (external) : nominal flux + flux tuning weights based on external data (NA61,...) only
  //  - tuned (final)    : nominal flux + flux tuning weights based on external data (NA61,...) and T2K detector data (ND280 data fit by BANFF,...)
  //  - current          : current flux (current systematic variation of combined-tuning flux)
  //
  TH1D *           fSKNuMuFluxNominal;           ///< nominal SuperK numu    flux (50 MeV binning)
  TH1D *           fSKNuMuBarFluxNominal;        ///< nominal SuperK numubar flux (50 MeV binning)
  TH1D *           fSKNuEFluxNominal;            ///< nominal SuperK nue     flux (50 MeV binning)
  TH1D *           fSKNuEBarFluxNominal;         ///< nominal SuperK nuebar  flux (50 MeV binning)
  TH1D *           fSKNuMuFluxTunedExt;          ///< tuned   SuperK numu    flux (50 MeV binning) based on external (NA61,...) data only 
  TH1D *           fSKNuMuBarFluxTunedExt;       ///< tuned   SuperK numubar flux (50 MeV binning) based on external (NA61,...) data only
  TH1D *           fSKNuEFluxTunedExt;           ///< tuned   SuperK nue     flux (50 MeV binning) based on external (NA61,...) data only
  TH1D *           fSKNuEBarFluxTunedExt;        ///< tuned   SuperK nuebar  flux (50 MeV binning) based on external (NA61,...) data only
  TH1D *           fSKNuMuFluxTunedFinal;        ///< tuned   SuperK numu    flux (50 MeV binning) based on external (NA61,...) + T2K detector (ND280,...) data
  TH1D *           fSKNuMuBarFluxTunedFinal;     ///< tuned   SuperK numubar flux (50 MeV binning) based on external (NA61,...) + T2K detector (ND280,...) data
  TH1D *           fSKNuEFluxTunedFinal;         ///< tuned   SuperK nue     flux (50 MeV binning) based on external (NA61,...) + T2K detector (ND280,...) data
  TH1D *           fSKNuEBarFluxTunedFinal;      ///< tuned   SuperK nuebar  flux (50 MeV binning) based on external (NA61,...) + T2K detector (ND280,...) data
  mutable TH1D *   fSKNuMuFluxCurrent;           ///< current SuperK numu    flux (50 MeV binning)
  mutable TH1D *   fSKNuMuBarFluxCurrent;        ///< current SuperK numubar flux (50 MeV binning)
  mutable TH1D *   fSKNuEFluxCurrent;            ///< current SuperK nue     flux (50 MeV binning)
  mutable TH1D *   fSKNuEBarFluxCurrent;         ///< current SuperK nuebar  flux (50 MeV binning)
  TH1D *           fSKNuMuFluxExtTuningWght;     ///< SuperK numu    flux tuning weights
  TH1D *           fSKNuMuBarFluxExtTuningWght;  ///< SuperK numubar flux tuning weights
  TH1D *           fSKNuEFluxExtTuningWght;      ///< SuperK nue     flux tuning weights
  TH1D *           fSKNuEBarFluxExtTuningWght;   ///< SuperK nue     flux tuning weights
  BANFF *          fBANFF;                       ///< BANFF fit info
      
ClassDef(NuFlux, 5)
};

} //t2kana

#endif
