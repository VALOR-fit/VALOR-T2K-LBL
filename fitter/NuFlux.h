//____________________________________________________________________________
/*!

\class    valor::NuFlux

\brief    Neutrino flux inputs to the oscillation analysis

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _valor_nu_flux_h_
#define _valor_nu_flux_h_

#include <TObject.h>
#include <TMatrixD.h>
#include <TMatrixTSym.h>
#include <TXMLEngine.h>
#include <TH1D.h>

#include "Conventions.h"
#include "Enumerations.h"

class TVirtualPad;

namespace valor {

class BANFF;

class NuFlux: public TObject
{
public:
  NuFlux();   
 ~NuFlux();
 
  using TObject::Read; // Suppress clang warnings
  using TObject::Draw; // Suppress clang warnings

  // SK flux tuning and systematic variation weights
  double SuperkFluxExtTuningWght     (int neutrino_pdg, double Etrue) const;  // wght: nominal --> ext/tuned flux
  double SuperkFluxFinTuningWght     (int neutrino_pdg, double Etrue) const;  // wght: ext/tuned --> fin/tuned flux
  double SuperkFluxSystWght          (int neutrino_pdg, double Etrue, std::string option = "") const;  // wght: fin/tuned flux -> current (incl. syst. variation)
  double SuperkFlux2p2hSystWght       (int neutrino_pdg, double Etrue) const;  // wght: ext/tuned flux -> current (incl. syst. variation)
  double SuperkRHCFluxExtTuningWght  (int neutrino_pdg, double Etrue) const;  // wght: nominal --> ext/tuned flux
  double SuperkRHCFluxFinTuningWght  (int neutrino_pdg, double Etrue) const;  // wght: ext/tuned --> fin/tuned flux
  double SuperkRHCFluxSystWght       (int neutrino_pdg, double Etrue, std::string option = "") const;  // wght: fin/tuned flux -> current (incl. syst. variation)
  double SuperkRHCFlux2p2hSystWght    (int neutrino_pdg, double Etrue) const;  // wght: ext/tuned flux -> current (incl. syst. variation)  

  // ... as above, slighty different interface
  // note that these check for FHC/RHC so there is only a single copy of each
  double SuperkFluxExtTuningWght  (conv::SuperkSample_t s, conv::SuperkMode_t m, double Etrue) const;  
  double SuperkFluxFinTuningWght  (conv::SuperkSample_t s, conv::SuperkMode_t m, double Etrue) const;  
  double SuperkFluxSystWght       (conv::SuperkSample_t s, conv::SuperkMode_t m, double Etrue, std::string option = "") const;

  // SK flux (current, tuned/fin, tuned/ext, nominal - depending on input flux type)
  double SuperkFlux (conv::SuperkSample_t s,int neutrino_pdg, double Etrue, conv::FluxType_t type = conv::kFlxCurrent) const;
  double SuperkFlux (int neutrino_pdg, double Etrue, conv::FluxType_t type = conv::kFlxCurrent) const;
  double RHCSuperkFlux (int neutrino_pdg, double Etrue, conv::FluxType_t type = conv::kFlxCurrent) const;

  // Read configuration and data from input XML file
  bool Read (TXMLEngine* xml, XMLNodePointer_t root_node);
  bool ReadRHC (TXMLEngine* xml, XMLNodePointer_t root_node);

  // NEW FLUX FHC+RHC
  //bool ReadFHCpRHC (TXMLEngine* xml, XMLNodePointer_t root_node);
  //

  // BANFF fit info
  void TellMeAboutTheBANFFfit(BANFF * banff);
  void TellMeAboutTheRHCBANFFfit(BANFF * banff);

  // NEW FLUX FHC+RHC
  //void TellMeAboutTheFHCpRHCBANFFfit(BANFF * banff);
  //

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
  // Drawing options:
  // - "Flux/SK/NuMuRHC/TuneWghts"    :
  // - "Flux/SK/NuMuBarRHC/TuneWghts" :
  // - "Flux/SK/NuERHC/TuneWghts"     :
  // - "Flux/SK/NuEBarRHC/TuneWghts"  :
  // - "Flux/SK/NuMuRHC"              :
  // - "Flux/SK/NuMuBarRHC"           :
  // - "Flux/SK/NuERHC"               :
  // - "Flux/SK/NuEBarRHC"            :
  void DrawRHC(TVirtualPad * c, Option_t * opt) const;

private:

  // Check that BANFFSystConv.h is setup correctly
  void CheckListOfBANFFParams(int nparams, int nfluxparams, const conv::PhysParamId_t * params, const char * vname, bool allow_gap);
  void CheckIsFluxParam(int nparams, int nfluxparams, const conv::PhysParamId_t * params, const bool   * flux,   const char * vname);
  void CheckFluxPdg    (int nparams, int nfluxparams, const conv::PhysParamId_t * params, const int    * pdg,    const char * vname);
  void CheckFluxERange (int nparams, int nfluxparams, const conv::PhysParamId_t * params, const double * erange, const char * vname);

  // Initialize and clean-up
  void Init    (void);
  void CleanUp (void);

  // Book and initialize/compute spectra not loaded from config files
  void BookAndInitAuxSpectra(void);
  void BookAndInitRHCAuxSpectra(void);

  // NEW FLUX FHC+RHC
  //void BookAndInitFHCpRHCAuxSpectra(void);
  //

  // Apply BANFF weights
  void ApplyBANFFWeights (bool get_current) const; 
  void ApplyRHCBANFFWeights (bool get_current) const;


  // NEW FLUX FHC+RHC
  //void ApplyFHCpRHCBANFFWeights (bool get_current) const;
  //


  // Access flux and weight histograms held in this object
  TH1D * SuperkFluxSpectrum       (int neutrino_pdg, conv::FluxType_t type = conv::kFlxCurrent) const;
  TH1D * SuperkFluxExtTuningWghts (int neutrino_pdg) const;
  //TH1D * SuperkRHCFluxSpectrum       (int neutrino_pdg, conv::FluxType_t type = conv::kRHCFlxCurrent) const;
  TH1D * SuperkRHCFluxExtTuningWghts (int neutrino_pdg) const;

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
  TH1D *           fSKRHCNuMuFluxNominal;           ///< nominal SuperK numu    flux (50 MeV binning)
  TH1D *           fSKRHCNuMuBarFluxNominal;        ///< nominal SuperK numubar flux (50 MeV binning)
  TH1D *           fSKRHCNuEFluxNominal;            ///< nominal SuperK nue     flux (50 MeV binning)
  TH1D *           fSKRHCNuEBarFluxNominal;         ///< nominal SuperK nuebar  flux (50 MeV binning)
  TH1D *           fSKRHCNuMuFluxTunedExt;          ///< tuned   SuperK numu    flux (50 MeV binning) based on external (NA61,...) data only 
  TH1D *           fSKRHCNuMuBarFluxTunedExt;       ///< tuned   SuperK numubar flux (50 MeV binning) based on external (NA61,...) data only
  TH1D *           fSKRHCNuEFluxTunedExt;           ///< tuned   SuperK nue     flux (50 MeV binning) based on external (NA61,...) data only
  TH1D *           fSKRHCNuEBarFluxTunedExt;        ///< tuned   SuperK nuebar  flux (50 MeV binning) based on external (NA61,...) data only
  TH1D *           fSKRHCNuMuFluxTunedFinal;        ///< tuned   SuperK numu    flux (50 MeV binning) based on external (NA61,...) + T2K detector (ND280,...) data
  TH1D *           fSKRHCNuMuBarFluxTunedFinal;     ///< tuned   SuperK numubar flux (50 MeV binning) based on external (NA61,...) + T2K detector (ND280,...) data
  TH1D *           fSKRHCNuEFluxTunedFinal;         ///< tuned   SuperK nue     flux (50 MeV binning) based on external (NA61,...) + T2K detector (ND280,...) data
  TH1D *           fSKRHCNuEBarFluxTunedFinal;      ///< tuned   SuperK nuebar  flux (50 MeV binning) based on external (NA61,...) + T2K detector (ND280,...) data
  mutable TH1D *   fSKRHCNuMuFluxCurrent;           ///< current SuperK numu    flux (50 MeV binning)
  mutable TH1D *   fSKRHCNuMuBarFluxCurrent;        ///< current SuperK numubar flux (50 MeV binning)
  mutable TH1D *   fSKRHCNuEFluxCurrent;            ///< current SuperK nue     flux (50 MeV binning)
  mutable TH1D *   fSKRHCNuEBarFluxCurrent;         ///< current SuperK nuebar  flux (50 MeV binning)
  TH1D *           fSKRHCNuMuFluxExtTuningWght;     ///< SuperK numu    flux tuning weights
  TH1D *           fSKRHCNuMuBarFluxExtTuningWght;  ///< SuperK numubar flux tuning weights
  TH1D *           fSKRHCNuEFluxExtTuningWght;      ///< SuperK nue     flux tuning weights
  TH1D *           fSKRHCNuEBarFluxExtTuningWght;   ///< SuperK nue     flux tuning weights
      
ClassDef(NuFlux, 6)
};

} //valor

#endif
