//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @ 18/2/13 - LE
 Added back the selection of dmsq schema,from branch 2010_v3: in the 3-flavour analyses, take the global NuOscDmsqScheme_t option into account. 
 Doubled arrays (independent ones for 1Rmu and 1Re)
 Number of Etrue, Ereco bins can be different for 1Rmu and 1Re.

 31-12-2013  Davide
 Added samples and modes for reverse horn current running

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__
//#define __show_init_mesg_debug__
//#define __ereco_scale_debug__
//#define __show_init_mesg_debug2__
//#define __show_update_mesg_info__
//#define __show_hypothesis__ 
//#define __show_hypothesis_probabilities__
//#define __show_banfftun_syst_var_wght_mesg_debug__
//#define __show_syst_var_wght_mesg_debug__ 
//#define __show_spectra_update_mesg_debug__
//#define __show_skipped_templates__
//#define __show_isanevent_mesg_debug__
//#define __show_negative_events_mesg_debug__
//#define __save_templates_debug__   // saves to files, and cout's number of events, at various stages of applying tunes to the templates:
// 1 number of events
// 2 number of events, weighted by pot & mc normalisation
// 3 number of events, weighted by pot & mc normalisation, flux tuning
// 4 number of events, weighted by pot & mc normalisation, flux tuning, banff tuning
// 5 number of events, weighted by pot & mc normalisation, flux tuning, banff tuning, multiplicative systematic variations
// 6 number of events, weighted by pot & mc normalisation, flux tuning, banff tuning, multiplicative systematic variations, nonmultiplicative systematic variations
// 7 number of events, weighted by pot & mc normalisation, flux tuning, banff tuning, multiplicative systematic variations, nonmultiplicative systematic variations, physics hypothesis
#ifdef __save_templates_debug__
//#define __save_templates_debug_mcstat__ //do the extra MCStatErrDebug
#endif
//#define __apply_minevents_to_individual_modes__ //special debug mode that applies kSetNegEventsTo to all modes. Used for comparisons with the KDI 1 sigma spectra
//#define NO_EVENTS_CAN_BE_NEGATIVE // if defined, will stop events being negative at all stages (not just building the Ereco spectra)

//#define POTWEIGHTSAREONE
//#define FLUXWEIGHTSAREONE
//#define BANFFWEIGHTSAREONE
//#define SYSTWEIGHTSAREONE
//#define PHYSHYPWEIGHTSAREONE
//#define PHYSHYPWEIGHTAREONEANDNUEISZERO //needs PHYSHYPWEIGHTSAREONE to work
#define INCLUDEEVENTSOVERTHIRTYGEV // include Ereco overflow events in the last bin?

#define METHOD_FOR_ESCALE 1 // 1 = new linear method, 0 = old (should be deprecated and removed once new method fully validated)

#define USE_FLUX_EXT_WEIGHT // Instead of using 'wgt' branch in files, use '_ratio' flux plots.

//#define __INCLUDE_MC_STAT_FLUCTUATIONS__ // fluctuate bin contents in final MC templates based on MC stat error in each bin

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include <TTree.h>
#include <TH2D.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TMath.h>
#include <TCanvas.h>

#include "SuperkMC.h"
#include "SuperkDetectorFSI.h"
#include "NuXSec.h"
#include "NuFlux.h"
#include "Utilities.h"
#include "Conventions/Opt.h"
#include "Utilities/Rnd.h"
#include "Utilities/NuOscSterile.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace valor;
using namespace valor::conv;
using namespace valor::util;

ClassImp(SuperkMC)

//________________________________________________________________________________________
SuperkMC::SuperkMC()
{
  this->Init();
}
//________________________________________________________________________________________
SuperkMC::~SuperkMC()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkMC::Init(void)
{
  fPOT = 0;
  fPOTRHC = 0;
  fXMEC = -1;
  fUnphysicalOscProb = false;

  for(int is=0; is<kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    int n_etrue_bins = kNEtrueBins[is];
    int n_ereco_bins = kNErecoBins[is];
    int n_mc_etrue_bins = kMCNEtrueBins;//not sample dependent
    int n_xs_ereco_bins = kNXsErecoBins[is];
    
    fMCVariableTemplArray   .push_back( VALORArray(knm, n_mc_etrue_bins, n_xs_ereco_bins) );
    fMCTemplArrayOri        .push_back( VALORArray(knm, n_etrue_bins   , n_ereco_bins   ) );
    fMCTemplArrayNom        .push_back( VALORArray(knm, n_etrue_bins   , n_ereco_bins   ) );
    fMCTemplArrayTun        .push_back( VALORArray(knm, n_etrue_bins   , n_ereco_bins   ) );
    fMCTemplArrayFin        .push_back( VALORArray(knm, n_etrue_bins   , n_ereco_bins   ) );
    fMCStatFracErr          .push_back( VALORArray(knm, n_etrue_bins   , n_ereco_bins   ) );
    fMCEtrueSpectrumTotArray.push_back( VALORArray(n_etrue_bins) );
    fMCErecoSpectrumTotArray.push_back( VALORArray(n_ereco_bins) );
    fMCErecoSpectrumArray   .push_back( VALORArray(knm, n_ereco_bins) );
    
    fMCVariableTemplArray   [is].SetName( Form("mc_templ_variable_%s"    , util::AsStringCompact(s).c_str()));
    fMCTemplArrayOri        [is].SetName( Form("mc_templ_ori_%s"         , util::AsStringCompact(s).c_str()));
    fMCTemplArrayNom        [is].SetName( Form("mc_templ_nom_%s"         , util::AsStringCompact(s).c_str()));
    fMCTemplArrayTun        [is].SetName( Form("mc_templ_tun_%s"         , util::AsStringCompact(s).c_str()));
    fMCTemplArrayFin        [is].SetName( Form("mc_templ_fin_%s"         , util::AsStringCompact(s).c_str()));
    fMCStatFracErr          [is].SetName( Form("mc_staterr_%s"           , util::AsStringCompact(s).c_str()));
    fMCEtrueSpectrumTotArray[is].SetName( Form("mc_spectrum_etrue_tot_%s", util::AsStringCompact(s).c_str()));
    fMCErecoSpectrumTotArray[is].SetName( Form("mc_spectrum_obs_tot_%s"  , util::AsStringCompact(s).c_str()));
    fMCErecoSpectrumArray   [is].SetName( Form("mc_spectrum_obs_%s"      , util::AsStringCompact(s).c_str()));
    
    fMCVariableTemplArray   [is].SetTitle( Form("MC variable template  %s"     , util::AsString(s).c_str()));
    fMCTemplArrayOri        [is].SetTitle( Form("MC template ori %s"           , util::AsString(s).c_str()));
    fMCTemplArrayNom        [is].SetTitle( Form("MC template nom %s"           , util::AsString(s).c_str()));
    fMCTemplArrayTun        [is].SetTitle( Form("MC template tun %s"           , util::AsString(s).c_str()));
    fMCTemplArrayFin        [is].SetTitle( Form("MC template fin %s"           , util::AsString(s).c_str()));
    fMCStatFracErr          [is].SetTitle( Form("MC Fractional StatErr %s"     , util::AsString(s).c_str()));
    fMCEtrueSpectrumTotArray[is].SetTitle( Form("MC Spectrum Etrue_Total %s"   , util::AsString(s).c_str()));
    fMCErecoSpectrumTotArray[is].SetTitle( Form("MC Spectrum Ereco Total %s"   , util::AsString(s).c_str()));
    fMCErecoSpectrumArray   [is].SetTitle( Form("MC Spectrum Ereco per mode %s", util::AsString(s).c_str()));
    
    fMCErecoSpectrumTot[is] = 0;
    fMCEtrueSpectrumTot[is] = 0;
    for(int im=0; im<knm; im++) {
      fMCTemplOri      [is][im] = 0;
      fMCTemplNom      [is][im] = 0;
      
      fMCErecoSpectrum [is][im] = 0;
      fMCVariableTempl [is][im] = 0;
      
      fNNonEmptyVariableBins [is][im] = 0;
      for (int ivb = 0; ivb < kMaxVarBins; ivb++) {
        fNNonEmptyFixedBins    [is][im][ivb] = 0;
        fNonEmptyFixedBins     [is][im][ivb] = NULL;
        fNonEmptyVariableBins  [is][im][ivb] = 0;
      }
    }//im
    fMCTreeNuMu     .push_back(0);
    fMCTreeNuMuBar  .push_back(0);
    fMCTreeNuE      .push_back(0);
    fMCTreeNuEBar   .push_back(0);
    fMCTreeOscNuE   .push_back(0);
    fMCTreeOscNuEBar.push_back(0);
  }//is
  
  ////  

  // FHC
  fMCPOTNuMu           = 0;
  fMCPOTNuMuBar        = 0;
  fMCPOTNuE            = 0;
  fMCPOTNuEBar         = 0;
  fMCPOTOscNuE         = 0;
  fMCPOTOscNuEBar      = 0;

  // RHC
  fMCPOTRHCNuMu           = 0;
  fMCPOTRHCNuMuBar        = 0;
  fMCPOTRHCNuE            = 0;
  fMCPOTRHCNuEBar         = 0;
  fMCPOTRHCOscNuE         = 0;
  fMCPOTRHCOscNuEBar      = 0;

  // NEW MULTIRING
  fMCPOTMultiRingNuMu           = 0;
  fMCPOTMultiRingNuMuBar        = 0;
  fMCPOTMultiRingNuE            = 0;
  fMCPOTMultiRingNuEBar         = 0;
  fMCPOTMultiRingOscNuE         = 0;
  fMCPOTMultiRingOscNuEBar      = 0;

  ////

  fBrEreco             = 0;
  fBrGenCode           = 0;
  
  fNuFlux              = 0;
  fNuXSec              = 0;  
  fSKDetFSI            = 0;
  
  arr_index = util::ArrIdx::Instance();

  cout << "Exp/SKMC >> Init complete" << endl;

}
//________________________________________________________________________________________
void SuperkMC::CleanUp(void)
{
  for(int is = 0; is < kns; is++) {
    for(int im = 0; im < knm; im++) {
      //already deleted Ori & Nom (the if doesn't work!)
      //if( fMCTemplOri     [is][im] ) delete fMCTemplOri     [is][im];
      //if( fMCTemplNom     [is][im] ) delete fMCTemplNom     [is][im];
      
      if( fMCErecoSpectrum[is][im] ) delete fMCErecoSpectrum[is][im];
      if( fMCVariableTempl[is][im] ) delete fMCVariableTempl[is][im];
      
      for (int ivb = 0; ivb < fNNonEmptyVariableBins[is][im] ; ivb++) {
        if (fNonEmptyFixedBins[is][im][ivb]) delete fNonEmptyFixedBins[is][im][ivb];
      }
    } //im
    if( fMCErecoSpectrumTot [is] ) delete fMCErecoSpectrumTot [is];
    if( fMCEtrueSpectrumTot [is] ) delete fMCEtrueSpectrumTot [is];
  }
}
//________________________________________________________________________________________
bool SuperkMC::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
  this->CleanUp();

  // instantiate all necessary histograms
  this->Book();

#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> Reading SuperK MC inputs from input XML file." << endl;
#endif

  //
  // Get data POT
  // 
  fPOT = util::XML2Dbl(xml,root_node,"data_pot/superk/fhc");
#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> " << "POT_{data} (FHC) = " << fPOT << endl;
#endif
  //
  // Get data POT RHC
  // 
  fPOTRHC = util::XML2Dbl(xml,root_node,"data_pot/superk/rhc");
#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> " << "POT_{data} (RHC) = " << fPOTRHC << endl;
#endif

  // NEW MULTIRING 

  //
  // Get data POT FHC MultiRing
  // 

#ifdef __MULTIRING_FHC_RUNNING__
  fPOTMultiRing = util::XML2Dbl(xml,root_node,"data_pot/superk/multiring_fhc");
#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> " << "POT_{data} (MultiRing) = " << fPOTMultiRing << endl;
#endif
#endif
  ///
  
  //
  // Get MC trees and MC POT
  //

  // FHC
  fMCPOTHistNuMu      = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numu/pot_fhc"         );
  fMCPOTHistNuMuBar   = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numubar/pot_fhc"      );
  fMCPOTHistNuE       = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nue/pot_fhc"          );
  fMCPOTHistNuEBar    = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nuebar/pot_fhc"       );
  fMCPOTHistOscNuE    = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/signal_nue/pot_fhc"   );
  fMCPOTHistOscNuEBar = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/signal_nuebar/pot_fhc");

  fMCPOTNuMu       = fMCPOTHistNuMu     ->GetBinContent(2) / fMCPOTHistNuMu     ->GetBinContent(1);
  fMCPOTNuMuBar    = fMCPOTHistNuMuBar  ->GetBinContent(2) / fMCPOTHistNuMuBar  ->GetBinContent(1);
  fMCPOTNuE        = fMCPOTHistNuE      ->GetBinContent(2) / fMCPOTHistNuE      ->GetBinContent(1);
  fMCPOTNuEBar     = fMCPOTHistNuEBar   ->GetBinContent(2) / fMCPOTHistNuEBar   ->GetBinContent(1);
  fMCPOTOscNuE     = fMCPOTHistOscNuE   ->GetBinContent(2) / fMCPOTHistOscNuE   ->GetBinContent(1);
  fMCPOTOscNuEBar  = fMCPOTHistOscNuEBar->GetBinContent(2) / fMCPOTHistOscNuEBar->GetBinContent(1);

  // RHC
  fMCPOTHistRHCNuMu      = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numu/pot_rhc"         );
  fMCPOTHistRHCNuMuBar   = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numubar/pot_rhc"      );
  fMCPOTHistRHCNuE       = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nue/pot_rhc"          );
  fMCPOTHistRHCNuEBar    = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nuebar/pot_rhc"       );
  fMCPOTHistRHCOscNuE    = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/signal_nue/pot_rhc"   );
  fMCPOTHistRHCOscNuEBar = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/signal_nuebar/pot_rhc");

  fMCPOTRHCNuMu    = fMCPOTHistRHCNuMu     ->GetBinContent(2) / fMCPOTHistRHCNuMu     ->GetBinContent(1);
  fMCPOTRHCNuMuBar = fMCPOTHistRHCNuMuBar  ->GetBinContent(2) / fMCPOTHistRHCNuMuBar  ->GetBinContent(1);
  fMCPOTRHCNuE     = fMCPOTHistRHCNuE      ->GetBinContent(2) / fMCPOTHistRHCNuE      ->GetBinContent(1);
  fMCPOTRHCNuEBar  = fMCPOTHistRHCNuEBar   ->GetBinContent(2) / fMCPOTHistRHCNuEBar   ->GetBinContent(1);
  fMCPOTRHCOscNuE  = fMCPOTHistRHCOscNuE   ->GetBinContent(2) / fMCPOTHistRHCOscNuE   ->GetBinContent(1);
  fMCPOTRHCOscNuE  = fMCPOTHistRHCOscNuEBar->GetBinContent(2) / fMCPOTHistRHCOscNuEBar->GetBinContent(1);

  // NEW MULTIRING
#ifdef __MULTIRING_FHC_RUNNING__
  fMCPOTHistMultiRingNuMu       = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numu/pot");
  fMCPOTHistMultiRingNuMuBar    = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numubar/pot");
  fMCPOTHistMultiRingNuE        = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nue/pot");
  fMCPOTHistMultiRingNuEBar     = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nuebar/pot");
  fMCPOTHistMultiRingOscNuE     = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/signal_nue/pot");
  fMCPOTHistMultiRingOscNuEBar  = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/signal_nuebar/pot");
  
  fMCPOTMultiRingNuMu       = fMCPOTHistMultiRingNuMu     ->GetBinContent(2) / fMCPOTHistMultiRingNuMu     ->GetBinContent(1);
  fMCPOTMultiRingNuMuBar    = fMCPOTHistMultiRingNuMuBar  ->GetBinContent(2) / fMCPOTHistMultiRingNuMuBar  ->GetBinContent(1);
  fMCPOTMultiRingNuE        = fMCPOTHistMultiRingNuE      ->GetBinContent(2) / fMCPOTHistMultiRingNuE      ->GetBinContent(1);
  fMCPOTMultiRingNuEBar     = fMCPOTHistMultiRingNuEBar   ->GetBinContent(2) / fMCPOTHistMultiRingNuEBar   ->GetBinContent(1);
  fMCPOTMultiRingOscNuE     = fMCPOTHistMultiRingOscNuE   ->GetBinContent(2) / fMCPOTHistMultiRingOscNuE   ->GetBinContent(1);
  fMCPOTMultiRingOscNuEBar  = fMCPOTHistMultiRingOscNuEBar->GetBinContent(2) / fMCPOTHistMultiRingOscNuEBar->GetBinContent(1);
#endif
  ////
  
  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if (s==kSKsUndefined) continue;
     fMCTreeNuMu     [is] = util::XML2TTree (xml,root_node,Form("mc/superk/nominal/numu/%s"         ,util::GetXMLMCTreeName(s).c_str()));
     fMCTreeNuMuBar  [is] = util::XML2TTree (xml,root_node,Form("mc/superk/nominal/numubar/%s"      ,util::GetXMLMCTreeName(s).c_str()));
     fMCTreeNuE      [is] = util::XML2TTree (xml,root_node,Form("mc/superk/nominal/nue/%s"          ,util::GetXMLMCTreeName(s).c_str()));
     fMCTreeNuEBar   [is] = util::XML2TTree (xml,root_node,Form("mc/superk/nominal/nuebar/%s"       ,util::GetXMLMCTreeName(s).c_str()));
     fMCTreeOscNuE   [is] = util::XML2TTree (xml,root_node,Form("mc/superk/nominal/signal_nue/%s"   ,util::GetXMLMCTreeName(s).c_str()));
     fMCTreeOscNuEBar[is] = util::XML2TTree (xml,root_node,Form("mc/superk/nominal/signal_nuebar/%s",util::GetXMLMCTreeName(s).c_str()));
     assert (fMCTreeNuMu      [is]);
     assert (fMCTreeNuMuBar   [is]);
     assert (fMCTreeNuE       [is]);
     assert (fMCTreeNuEBar    [is]);
     assert (fMCTreeOscNuE    [is]);
     assert (fMCTreeOscNuEBar [is]);
  }




#ifdef __show_init_mesg_info__
  // FHC
  cout << "Exp/SKMC >> " << "POT {MC; numu}    = "    << fMCPOTNuMu       << endl;
  cout << "Exp/SKMC >> " << "POT {MC; numubar} = "    << fMCPOTNuMuBar    << endl;
  cout << "Exp/SKMC >> " << "POT {MC; nue}     = "    << fMCPOTNuE        << endl;
  cout << "Exp/SKMC >> " << "POT {MC; nuebar}  = "    << fMCPOTNuEBar     << endl;
  cout << "Exp/SKMC >> " << "POT {MC; osc_nue} = "    << fMCPOTOscNuE     << endl;
  cout << "Exp/SKMC >> " << "POT {MC; osc_nuebar} = " << fMCPOTOscNuEBar  << endl;
  
  // RHC
  cout << "RHC: Exp/SKMC >> " << "POT {MC; numu}    = "    << fMCPOTRHCNuMu       << endl;
  cout << "RHC: Exp/SKMC >> " << "POT {MC; numubar} = "    << fMCPOTRHCNuMuBar    << endl;
  cout << "RHC: Exp/SKMC >> " << "POT {MC; nue}     = "    << fMCPOTRHCNuE        << endl;
  cout << "RHC: Exp/SKMC >> " << "POT {MC; nuebar}  = "    << fMCPOTRHCNuEBar     << endl;
  cout << "RHC: Exp/SKMC >> " << "POT {MC; osc_nue} = "    << fMCPOTRHCOscNuE     << endl;
  cout << "RHC: Exp/SKMC >> " << "POT {MC; osc_nuebar} = " << fMCPOTRHCOscNuEBar  << endl;

  // NEW MULTIRING 
#ifdef __MULTIRING_FHC_RUNNING__
  cout << "MultiRing: Exp/SKMC >> " << "POT {MC; numu}       = " << fMCPOTMultiRingNuMu      << endl;
  cout << "MultiRing: Exp/SKMC >> " << "POT {MC; numubar}    = " << fMCPOTMultiRingNuMuBar   << endl;
  cout << "MultiRing: Exp/SKMC >> " << "POT {MC; nue}        = " << fMCPOTMultiRingNuE       << endl;
  cout << "MultiRing: Exp/SKMC >> " << "POT {MC; nuebar}     = " << fMCPOTMultiRingNuEBar    << endl;
  cout << "MultiRing: Exp/SKMC >> " << "POT {MC; osc_nue}    = " << fMCPOTMultiRingOscNuE    << endl;
  cout << "MultiRing: Exp/SKMC >> " << "POT {MC; osc_nuebar} = " << fMCPOTMultiRingOscNuEBar << endl;
#endif
  ////

  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if (s==kSKsUndefined) continue;
     
     cout << "Exp/SKMC >> " << "Nev {MC; numu;     "   <<AsStringCompact(s)<<"} = " << fMCTreeNuMu     [is]-> GetEntries()  << endl;
     cout << "Exp/SKMC >> " << "Nev {MC; numubar;  "   <<AsStringCompact(s)<<"} = " << fMCTreeNuMuBar  [is]-> GetEntries()  << endl;
     cout << "Exp/SKMC >> " << "Nev {MC; nue;      "   <<AsStringCompact(s)<<"} = " << fMCTreeNuE      [is]-> GetEntries()  << endl;
     cout << "Exp/SKMC >> " << "Nev {MC; nuebar;   "   <<AsStringCompact(s)<<"} = " << fMCTreeNuEBar   [is]-> GetEntries()  << endl;
     cout << "Exp/SKMC >> " << "Nev {MC; nue[osc]; "   <<AsStringCompact(s)<<"} = " << fMCTreeOscNuE   [is]-> GetEntries()  << endl;
     cout << "Exp/SKMC >> " << "Nev {MC; nuebar[osc]; "<<AsStringCompact(s)<<"} = " << fMCTreeOscNuEBar[is]-> GetEntries()  << endl;
     
  }
#endif//__show_init_mesg_info__

  // initialise physics hypothesis
  // need to initialise as 3 flavour oscillation, otherwise ExtTuning (NA61) and FinTuning (BANFF) is not applied to nuesig modes
  fPhysHypothesis.Init(kHypOsc3F);

  // fill-in MC templates from the input event tree
  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     for(int im = 0; im < knm; im++) {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       this->FillTemplate(s,m);
     }
  }

  // fill-in "special" MC templates (osc. numu CC)
  this->FillSpecialTemplates();

#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtOri, "nevents");
  this->NevDebug(kSKMCtNom, "mcnorm");
#endif

#ifdef __save_templates_debug__
    this->NevDebug(kSKMCtNom, "mcnorm_fluxtuned");
#endif

  this->CalculateMCStatError();
#ifdef __save_templates_debug_mcstat__
  this->MCStatErrDebug();
#endif

  // need to read in the xsec splines, before applying any BANFF tweaks
  fNuXSec -> BuildXsBinningMCTemplates(fMCTemplOri);
  if( !fNuXSec -> ReadSplines(xml,root_node) ) return false;

  // apply flux tuning (based on external data) to templates
  //this->ApplyFluxExtTuningToTemplates();



  // apply flux and cross-section tuning (based on nd280 data fit by BANFF) to templates
  this->ApplyBANFFTuningToTemplates();

  // fill the etrue spectra
  this->FillEtrueSpectra();


  //don't need Ori or Nom templates anymore
  for(int is = 0; is < kns; is++) {
    for(int im = 0; im < knm; im++) {
      if(fMCTemplOri[is][im]) delete fMCTemplOri[is][im];
      if(fMCTemplNom[is][im]) delete fMCTemplNom[is][im];
    }
  }

#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtTun, "mcnorm_fluxtuned_banfftuned");
#endif

#ifdef __show_init_mesg_info__      
  cout << "Exp/SKMC >> Done reading SuperK MC inputs!" << endl;
#endif
  
  return true;
}
//________________________________________________________________________________________
void SuperkMC::Book(void)
{
  //
  // Ereco histograms filled-in by SK spectrum prediction code
  // Histograms corresponding to different samples can have different binning
  //

  for(int is = 0; is < kns; is++) 
  {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;

     const int      n_ereco_bins    = kNErecoBins[is];
     const double * ereco_bin_edges = util::ErecoBinEdges(s);

     fMCErecoSpectrumTot[is] = new TH1D("","", n_ereco_bins, ereco_bin_edges);
     fMCErecoSpectrumTot[is] -> SetDirectory(0);
     fMCErecoSpectrumTot[is]->SetLineWidth(2);
     fMCErecoSpectrumTot[is]->SetLineColor(kBlack);
     
     fMCEtrueSpectrumTot[is] = new TH1D("","", kNEtrueBins[is], util::EtrueBinEdges(s));
     
     fMCEtrueSpectrumTot[is] -> SetDirectory(0);
     fMCEtrueSpectrumTot[is]->SetLineWidth(2);
     fMCEtrueSpectrumTot[is]->SetLineColor(kBlack);
     for(int im=0; im < knm; im++) 
     {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       fMCErecoSpectrum[is][im] = new TH1D("","", n_ereco_bins, ereco_bin_edges);
       fMCErecoSpectrum[is][im] -> SetDirectory(0);
       util::SetDefaultSKSpectrumStyle(s,m,fMCErecoSpectrum[is][im]);
     }//m
  }//s

  //
  // 2-D(Ereco,Etrue) MC templates
  // Histograms corresponding to different samples can have different binning
  //

  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;

     const int      n_ereco_bins    = kNErecoBins[is];
     const int      n_etrue_bins    = kNEtrueBins[is];
     const double * ereco_bin_edges = util::ErecoBinEdges(s);
     const int      var_n_ereco_bins    = kNXsErecoBins[is];
     const double * var_ereco_bin_edges = util::XsecErecoBinEdges(s);

     for(int im = 0; im < knm; im++) {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       const int      var_n_etrue_bins    = util::NumVariableMCEtrueBins (s,m);
       const double * var_etrue_bin_edges = util::VariableMCEtrueBinEdges(s,m);
       fMCTemplOri[is][im] = new TH2D(
          Form("nsk_ori_%d_%d",is,im), util::AsString(m).c_str(),
	  n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
       fMCTemplNom[is][im] = new TH2D(
	  Form("nsk_nom_%d_%d",is,im), util::AsString(m).c_str(),
	  n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
	  
       fMCTemplOri[is][im]->SetDirectory(0);
       fMCTemplOri[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fMCTemplOri[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
       fMCTemplNom[is][im]->SetDirectory(0);
       fMCTemplNom[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fMCTemplNom[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
       fMCVariableTempl[is][im] = new TH2D(
          Form("nsk_ori_var_%d_%d",is,im), util::AsString(m).c_str(),
	  var_n_ereco_bins, var_ereco_bin_edges, var_n_etrue_bins,  var_etrue_bin_edges);
       fMCVariableTempl[is][im]->SetDirectory(0);
       fMCVariableTempl[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fMCVariableTempl[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
     }
  }

}
//________________________________________________________________________________________
void SuperkMC::FillTemplate(SuperkSample_t s, SuperkMode_t m)
{
#ifdef __show_init_mesg_debug__      
  cout << "Exp/SKMC/CreateTemplate >> Creating MC template for: "
       << util::AsString(s) << " / " << util::AsString(m) << endl;
#endif

  //filled in FillSpecialTemplates()
  //~ if(util::IsOscNuMu(m) || util::IsOscNuMuBar(m) || util::IsOscNuEBar(m))
  if(util::IsOscNuMu(m) || util::IsOscNuMuBar(m))
    return;

  int sidx = arr_index->SuperkSampleIdx(s);
  int midx = arr_index->SuperkModeIdx(s,m);

  TTree * mctree = this->MCTree(s,m);
  double  mcpot  = this->MCPOT (s,m);

  if(!mctree) {
     cerr << "Exp/SKMC/CreateTemplate >> Couldn't find event tree for template: " 
          << util::AsString(s) << "/" << util::AsString(m) << endl;
     exit(1);
  }

  double pot_wght = this->AnalyzedSamplePOT(s) / mcpot;

#ifdef POTWEIGHTSAREONE
  pot_wght = 1;
#endif

#ifdef __show_init_mesg_debug__        
  cout << "Exp/SKMC/CreateTemplate >>      "
       << "Weight (data POT / MC POT) = " << pot_wght << endl;
#endif
  if(pot_wght < 0) {
     cerr << "Exp/SKMC/CreateTemplate >> Weight (data POT / MC POT) for " 
          << util::AsString(s) << "/" << util::AsString(m)
          << " is " << pot_wght << endl;
     exit(1);
  }
  
  const char * erec_branch      = "erec";
  const char * etrue_branch     = "pnu";
  const char * mode_branch      = "mode";
  const char * flux_wght_branch = "wgtflx";
  //~ const char * flux_wght_branch = "wgt";

  TBranch * branchEreco   = mctree->GetBranch(erec_branch      );
  TBranch * branchEtrue   = mctree->GetBranch(etrue_branch     );
  TBranch * branchGenCode = mctree->GetBranch(mode_branch      );
  TBranch * branchWght    = mctree->GetBranch(flux_wght_branch );
  
  branchEreco   -> SetAddress(&fBrEreco);
  branchEtrue   -> SetAddress(&fBrPtrue);
  branchGenCode -> SetAddress(&fBrGenCode);
  branchWght    -> SetAddress(&fBrWght);  

  mctree->SetBranchStatus("*",0);
  mctree->SetBranchStatus(erec_branch      , 1 );
  mctree->SetBranchStatus(etrue_branch     , 1 );
  mctree->SetBranchStatus(mode_branch      , 1 );
  mctree->SetBranchStatus(flux_wght_branch , 1 );

  for(int iev = 0; iev < mctree->GetEntries(); iev++) 
  {
     mctree->GetEntry(iev);

     if(fBrEreco<0) continue;
     
     // use current entry for mode m template?
     bool use_it = util::UseEventToFillInTemplate(m,util::SuperkMode2SKMCPdg(m),fBrGenCode);
     if(!use_it) continue;

     double ereco_fill = fBrEreco/1000.;
     
     double Etrue = fBrPtrue[0];
     
     //if(util::IsRHCSample(s))
#ifdef USE_FLUX_EXT_WEIGHT
     fBrWght = fNuFlux->SuperkFluxExtTuningWght(s,m,Etrue);
#endif
     
#ifdef INCLUDEEVENTSOVERTHIRTYGEV
     ereco_fill = TMath::Min(29.9, ereco_fill);
#endif
     this -> Template(s,m,kSKMCtOri) -> Fill (ereco_fill, Etrue);
     fMCVariableTempl[sidx][midx]    -> Fill (ereco_fill, Etrue);
#ifndef FLUXWEIGHTSAREONE
     this -> Template(s,m,kSKMCtNom) -> Fill (ereco_fill, Etrue, fBrWght * pot_wght);
     
#else
     this -> Template(s,m,kSKMCtNom) -> Fill (ereco_fill, Etrue, pot_wght);
#endif//FLUXWEIGHTSAREONE
     
#ifdef __show_init_mesg_debug__
     if(Etrue <= 0 || Etrue >= this->Template(s,m,kSKMCtOri)->GetYaxis()->GetBinUpEdge(this->Template(s,m,kSKMCtOri)->GetNbinsY()))
       cout << "Exp/SKMC/CreateTemplate >> Etrue is out of template range " << AsString(s) << " " << AsString(m) << " " << Etrue 
	    << " (Ereco = " << fBrEreco/1000 << ")" << endl;
#endif
#ifdef __show_init_mesg_debug2__
     if(fBrEreco/1000 <= 0)
       cout << "Exp/SKMC/CreateTemplate >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << fBrEreco/1000
	    << " (Etrue = " << Etrue << ")" << endl;
     //#endif __show_init_mesg_debug2__
     if(fBrEreco/1000 >= this->Template(s,m,kSKMCtOri)->GetXaxis()->GetBinUpEdge(this->Template(s,m,kSKMCtOri)->GetNbinsX()))
       cout << "Exp/SKMC/CreateTemplate >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << fBrEreco/1000
	    << " (Etrue = " << Etrue << ")" << endl;
#endif// __show_init_mesg_debug__

  }//iev

  if(sidx > -1 && midx > -1) {
    // loop over Etrue
    for(int it = 1; it <= this->Template(s,m,kSKMCtNom)->GetYaxis()->GetNbins(); it++) {
      // loop over Ereco
      for(int ir = 1; ir <= this->Template(s,m,kSKMCtNom)->GetXaxis()->GetNbins(); ir++) {

        
        fMCTemplArrayOri[sidx].GetEntryRef(midx, it-1 , ir-1) = this->Template(s,m,kSKMCtOri)->GetBinContent(ir,it);
        ////

        double nev = this->Template(s,m,kSKMCtNom)->GetBinContent(ir,it);
        fMCTemplArrayNom[sidx].GetEntryRef(midx, it-1 , ir-1) = nev;
        fMCTemplArrayTun[sidx].GetEntryRef(midx, it-1 , ir-1) = nev;
        fMCTemplArrayFin[sidx].GetEntryRef(midx, it-1 , ir-1) = nev;
        ////

      }//ir
    }//it
    
    int i_ne_var = 0;
    //now loop over the variable binning scheme (used to calculate weights only as many times as we need to)
    // This loop also keeps a record of which bins are non-empty
    //loop over variable Etrue
    for(int itv = 1; itv <= fMCVariableTempl[sidx][midx]->GetYaxis()->GetNbins(); itv++) {
      double EtrueVarMin = util::VariableMCEtrueBinEdges(s,m)[itv-1] + kOneMeV;
      double EtrueVarMax = util::VariableMCEtrueBinEdges(s,m)[itv]   - kOneMeV;

      int EtrueMinBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMin);
      int EtrueMaxBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMax);      
      /*
        cout <<
        "\nEmin: " << util::VariableMCEtrueBinEdges(s,m)[itv-1] <<
        " Emax: " << util::VariableMCEtrueBinEdges(s,m)[itv] <<
        "\nEtruemin Bin number: " << 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMin) <<
        " Etruemax Bin number: " << 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMax) << endl;
        */
      // loop over variable Ereco
      for(int irv = 1; irv <= fMCVariableTempl[sidx][midx]->GetXaxis()->GetNbins(); irv++) {
	double nev = fMCVariableTempl[sidx][midx]->GetBinContent(irv,itv);

	fMCVariableTemplArray[sidx].GetEntryRef(midx, itv-1, irv-1) = nev;
      
	// Fill nonempty bin numbers.
	if (nev < conv::kASmallNum) {
        //cout << irv << "exit" << endl;
	  continue;
	}
	
	fNonEmptyVariableBins[sidx][midx][i_ne_var] = conv::kBigIndex * irv + itv;
	
	double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
	double ErecoVarMax = util::XsecErecoBinEdges(s,irv)   - kOneMeV;

	int ErecoMinBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMin);
	int ErecoMaxBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMax);
/*
          cout <<
          "\n" <<
          irv <<
          "\nErmin: " << util::XsecErecoBinEdges(s,irv-1) <<
          " Ermax: " << util::XsecErecoBinEdges(s,irv) <<
          "\nErecomin Bin number: " << 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMin) <<
          " Erecomax Bin number: " << 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMax) << endl;
          
  */
	int i_ne_fixed = 0;
	// loop over fixed Etrue to count n nonempty
	for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
	  // loop over fixed Ereco to count
	  for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
	    double nevf = this->Template(s,m,kSKMCtOri)->GetBinContent(ir,it);
	    if (nevf < conv::kASmallNum) {
	      continue;
	    }
	    i_ne_fixed++;
	  }//ir
	}//it
	fNNonEmptyFixedBins[sidx][midx][i_ne_var] = i_ne_fixed;
	fNonEmptyFixedBins[sidx][midx][i_ne_var] = new int[i_ne_fixed];
	
	i_ne_fixed = 0;
	// loop over fixed Etrue to count n nonempty
	for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
	  // loop over fixed Ereco to count
	  for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
	    double nevf = this->Template(s,m,kSKMCtOri)->GetBinContent(ir,it);
	    if (nevf < conv::kASmallNum) {
	      continue;
	    }
	    fNonEmptyFixedBins[sidx][midx][i_ne_var][i_ne_fixed] = conv::kBigIndex * ir + it;
	    i_ne_fixed++;
	  }//ir
	}//it
          
          
	i_ne_var++;
      }//irv
    }//itv
    fNNonEmptyVariableBins[sidx][midx] = i_ne_var;
    assert(i_ne_var<kMaxVarBins);
    
  } //if valid index

#ifdef __show_init_mesg_debug__        
  cout << "Exp/SKMC/CreateTemplate >>      "
       << "Number of template entries = " 
       << this -> Template(s,m,kSKMCtNom) -> GetEntries() << endl;
#endif

}






//________________________________________________________________________________________
void SuperkMC::FillSpecialTemplates(void)
{
  //
  // Special templates for 1-ring mu-like sample:
  //

  //  SuperkSample_t s = kSKsFCFV1RingMuLike;
  //loop over samples to use it also for joint analysis
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    int sidx = arr_index->SuperkSampleIdx(s);
      
    EFluxType fNom = (IsRHCSample(s)) ? kRHCFlxNominal : kFlxNominal ;
    EFluxType fExt = (IsRHCSample(s)) ? kRHCFlxExtTuning : kFlxExtTuning ;
  
    //
    // Fill in `osc numu' templates from `numu' templates
    // and take into account the nue/numu flux weight
    // 
    
   const int n = 4;
    
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
  
    cout << "SKMC/FillSpecialTemplates >> Filling osc numu templates for " << util::AsString(s) << " sample" << endl;
    SuperkMode_t src_modes[n] = {  kSKmNuMuCCQE,    kSKmNuMuCC1pi,     kSKmNuMuCCcoh,    kSKmNuMuCCoth    };         
    SuperkMode_t tgt_modes[n] = {  kSKmOscNuMuCCQE, kSKmOscNuMuCC1pi,  kSKmOscNuMuCCcoh, kSKmOscNuMuCCoth };
    
    for(int i=0; i<n; i++) {
      
      SuperkMode_t tgt_m = tgt_modes[i];
      int tgt_midx = arr_index->SuperkModeIdx(s,tgt_modes[i]);
      int src_midx = arr_index->SuperkModeIdx(s,src_modes[i]);
      
      TH2D * src_ori = this->Template(s,src_modes[i],kSKMCtOri);
      TH2D * tgt_ori = this->Template(s,tgt_modes[i],kSKMCtOri);
      TH2D * src_nom = this->Template(s,src_modes[i],kSKMCtNom);
      TH2D * tgt_nom = this->Template(s,tgt_modes[i],kSKMCtNom);
      
      // loop over Etrue
      for(int it = 1; it <= src_ori->GetYaxis()->GetNbins(); it++) {
	double Etrue = src_ori->GetYaxis()->GetBinCenter(it);
	// calc weight for weighting numu flux to nue flux
	// (note: we have used `nominal' with different meaning in the flux and in the template enumerations.
	// perhaps we should rename kFlxNominal -> kFlxOriginal, for consistency).
	double flux_nue_ori  = fNuFlux->SuperkFlux (s, kPdgNuE,  Etrue, fNom   );
	double flux_nue_nom  = fNuFlux->SuperkFlux (s, kPdgNuE,  Etrue, fExt );
	double flux_numu_ori = fNuFlux->SuperkFlux (s, kPdgNuMu, Etrue, fNom   );
	double flux_numu_nom = fNuFlux->SuperkFlux (s, kPdgNuMu, Etrue, fExt ); 
	double wght_numu2nue_ori = (flux_numu_ori != 0) ? flux_nue_ori / flux_numu_ori : 0;
	double wght_numu2nue_nom = (flux_numu_nom != 0) ? flux_nue_nom / flux_numu_nom : 0;
	// loop over Ereco
	for(int ir = 1; ir <= src_ori->GetXaxis()->GetNbins(); ir++) {
	  double Ereco = src_ori->GetXaxis()->GetBinCenter(ir);
	  // get nevents
	  double nev_ori = src_ori-> GetBinContent(ir,it);
	  double nev_nom = src_nom-> GetBinContent(ir,it);
	  // set number of events
	  double weighted_nev_ori = nev_ori*wght_numu2nue_ori;
	  double weighted_nev_nom = nev_nom*wght_numu2nue_nom;
	  tgt_ori->Fill(Ereco,Etrue, weighted_nev_ori);
	  tgt_nom->Fill(Ereco,Etrue, weighted_nev_nom);
	  
	  double nev_tgt_ori = tgt_ori->GetBinContent(ir,it);
	  double nev_tgt_nom = tgt_nom->GetBinContent(ir,it);
            
            fMCTemplArrayOri[is].GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_ori;
            fMCTemplArrayNom[is].GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
            fMCTemplArrayTun[is].GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
            fMCTemplArrayFin[is].GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;

        }//ir
      }//it
      //one more loop, for the variable binning scheme.
      //just set the OscNuMu equal to the NuMu
      //loop over variable Etrue
      int i_ne_var = 0;
      for(int itv = 1; itv <= fMCVariableTempl[sidx][tgt_midx]->GetYaxis()->GetNbins(); itv++) {
	double EtrueVarMin = util::VariableMCEtrueBinEdges(s,tgt_m)[itv-1] + kOneMeV;
	double EtrueVarMax = util::VariableMCEtrueBinEdges(s,tgt_m)[itv]   - kOneMeV;

	int EtrueMinBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMin);
	int EtrueMaxBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMax);
	
	// loop over variable Ereco
	double nev = 0;
	for(int irv = 1; irv <= fMCVariableTempl[sidx][tgt_midx]->GetXaxis()->GetNbins(); irv++) {
	  
	  nev = fMCVariableTemplArray[is].GetEntry(src_midx, itv-1, irv-1);
	  fMCVariableTemplArray[is].GetEntryRef(tgt_midx, itv-1 , irv-1) = nev;
	    
	  // Fill nonempty bin numbers.
	  if (nev < conv::kASmallNum) {
	    continue;
	  }
	
	  fNonEmptyVariableBins[sidx][tgt_midx][i_ne_var] = conv::kBigIndex * irv + itv;

	  double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
	  double ErecoVarMax = util::XsecErecoBinEdges(s,irv)   - kOneMeV;
	  
	  int ErecoMinBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMin);
	  int ErecoMaxBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMax);

	  int i_ne_fixed = 0;
	  // loop over fixed Etrue to count n nonempty
	  for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
	    // loop over fixed Ereco to count
	    for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
	      double nevf = this->Template(s,tgt_m,kSKMCtOri)->GetBinContent(ir,it);
	      if (nevf < conv::kASmallNum) {
		continue;
	      }
	      i_ne_fixed++;
	    }//ir
	  }//it
	  fNNonEmptyFixedBins[sidx][tgt_midx][i_ne_var] = i_ne_fixed;
	  fNonEmptyFixedBins[sidx][tgt_midx][i_ne_var] = new int[i_ne_fixed];
	  
	  i_ne_fixed = 0;
	  // loop over fixed Etrue to count n nonempty
	  for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
	    // loop over fixed Ereco to count
	    for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
	      double nevf = this->Template(s,tgt_m,kSKMCtOri)->GetBinContent(ir,it);
	      if (nevf < conv::kASmallNum) {
		continue;
	      }
	      fNonEmptyFixedBins[sidx][tgt_midx][i_ne_var][i_ne_fixed] = conv::kBigIndex * ir + it;
	      i_ne_fixed++;
	    }//ir
	  }//it
	  i_ne_var++;
	}//irv
      }//itv
      fNNonEmptyVariableBins[sidx][tgt_midx] = i_ne_var;
      assert(i_ne_var<kMaxVarBins);

    }//i
    
#endif// __USE_OSC_NUMU_MC_TEMPLATES__



    //
    // Fill in `osc numubar/nuebar' templates from `numubar/nuebar' templates
    // and take into account the numubar/nuebar flux weight
    // 
    
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
  
    cout << "SKMC/FillSpecialTemplates >> Filling osc numubar templates for " << util::AsString(s) << " sample" << endl;
    SuperkMode_t src_modes_numubar[n] = {  kSKmNuMuBarCCQE,    kSKmNuMuBarCC1pi,     kSKmNuMuBarCCcoh,    kSKmNuMuBarCCoth    };    
    SuperkMode_t tgt_modes_numubar[n] = {  kSKmOscNuMuBarCCQE, kSKmOscNuMuBarCC1pi,  kSKmOscNuMuBarCCcoh, kSKmOscNuMuBarCCoth };
    
    for(int i=0; i<n; i++) {
      
      SuperkMode_t tgt_m = tgt_modes_numubar[i];
      int tgt_midx = arr_index->SuperkModeIdx(s,tgt_modes_numubar[i]);
      int src_midx = arr_index->SuperkModeIdx(s,src_modes_numubar[i]);
      
      TH2D * src_ori = this->Template(s,src_modes_numubar[i],kSKMCtOri);
      TH2D * tgt_ori = this->Template(s,tgt_modes_numubar[i],kSKMCtOri);
      TH2D * src_nom = this->Template(s,src_modes_numubar[i],kSKMCtNom);
      TH2D * tgt_nom = this->Template(s,tgt_modes_numubar[i],kSKMCtNom);
      
    // loop over Etrue
      for(int it = 1; it <= src_ori->GetYaxis()->GetNbins(); it++) {
	double Etrue = src_ori->GetYaxis()->GetBinCenter(it);
	// calc weight for weighting numu flux to nue flux
	// (note: we have used `nominal' with different meaning in the flux and in the template enumerations.
	// perhaps we should rename fNom -> kFlxOriginal, for consistency).
	double flux_nue_ori  = fNuFlux->SuperkFlux (s, kPdgNuEBar,  Etrue, fNom   );
	double flux_nue_nom  = fNuFlux->SuperkFlux (s, kPdgNuEBar,  Etrue, fExt ); 
	double flux_numu_ori = fNuFlux->SuperkFlux (s, kPdgNuMuBar, Etrue, fNom   );
	double flux_numu_nom = fNuFlux->SuperkFlux (s, kPdgNuMuBar, Etrue, fExt ); 
	double wght_numu2nue_ori = (flux_numu_ori != 0) ? flux_nue_ori / flux_numu_ori : 0;
	double wght_numu2nue_nom = (flux_numu_nom != 0) ? flux_nue_nom / flux_numu_nom : 0;
	// loop over Ereco
	for(int ir = 1; ir <= src_ori->GetXaxis()->GetNbins(); ir++) {
	  double Ereco = src_ori->GetXaxis()->GetBinCenter(ir);
	  // get nevents
	  double nev_ori = src_ori-> GetBinContent(ir,it);
	  double nev_nom = src_nom-> GetBinContent(ir,it);
	  // set number of events
	  double weighted_nev_ori = nev_ori*wght_numu2nue_ori;
	  double weighted_nev_nom = nev_nom*wght_numu2nue_nom;
	  tgt_ori->Fill(Ereco,Etrue, weighted_nev_ori);
	  tgt_nom->Fill(Ereco,Etrue, weighted_nev_nom);
	  
	  double nev_tgt_ori = tgt_ori->GetBinContent(ir,it);
	  double nev_tgt_nom = tgt_nom->GetBinContent(ir,it);
    
          fMCTemplArrayOri[is].GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_ori;
          fMCTemplArrayNom[is].GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
          fMCTemplArrayTun[is].GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
          fMCTemplArrayFin[is].GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
        }//ir
      }//it
      //one more loop, for the variable binning scheme.
      //just set the OscNuMu equal to the NuMu
      //loop over variable Etrue
      int i_ne_var = 0;

      for(int itv = 1; itv <= fMCVariableTempl[sidx][tgt_midx]->GetYaxis()->GetNbins(); itv++) {
	double EtrueVarMin = util::VariableMCEtrueBinEdges(s,tgt_m)[itv-1] + kOneMeV;
	double EtrueVarMax = util::VariableMCEtrueBinEdges(s,tgt_m)[itv]   - kOneMeV;
	
	int EtrueMinBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMin);
	int EtrueMaxBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMax);
	
	// loop over variable Ereco
	for(int irv = 1; irv <= fMCVariableTempl[sidx][tgt_midx]->GetXaxis()->GetNbins(); irv++) {
	  double nev = 0;
	  nev = fMCVariableTemplArray[is].GetEntry(src_midx, itv-1 , irv-1);
	  fMCVariableTemplArray[is].GetEntryRef(tgt_midx, itv-1, irv-1) = nev;
	  ////
	

	  // Fill nonempty bin numbers.
	  if (nev < conv::kASmallNum) {
	    continue;
	  }
	  
	  fNonEmptyVariableBins[sidx][tgt_midx][i_ne_var] = conv::kBigIndex * irv + itv;
	  
	  double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
	  double ErecoVarMax = util::XsecErecoBinEdges(s,irv)   - kOneMeV;
	  
	  int ErecoMinBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMin);
	  int ErecoMaxBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMax);

	  int i_ne_fixed = 0;
	  // loop over fixed Etrue to count n nonempty
	  for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
	    // loop over fixed Ereco to count
	    for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
	      double nevf = this->Template(s,tgt_m,kSKMCtOri)->GetBinContent(ir,it);
	      if (nevf < conv::kASmallNum) {
		continue;
	      }
	      i_ne_fixed++;
	    }//ir
	  }//it
	  fNNonEmptyFixedBins[sidx][tgt_midx][i_ne_var] = i_ne_fixed;
	  fNonEmptyFixedBins[sidx][tgt_midx][i_ne_var] = new int[i_ne_fixed];
	  
	  i_ne_fixed = 0;
	  // loop over fixed Etrue to count n nonempty
	  for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
	    // loop over fixed Ereco to count
	    for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
	      double nevf = this->Template(s,tgt_m,kSKMCtOri)->GetBinContent(ir,it);
	      if (nevf < conv::kASmallNum) {
		continue;
	      }
	      fNonEmptyFixedBins[sidx][tgt_midx][i_ne_var][i_ne_fixed] = conv::kBigIndex * ir + it;
	      i_ne_fixed++;
	    }//ir
	  }//it
	  i_ne_var++;
	}//irv
      }//itv
      fNNonEmptyVariableBins[sidx][tgt_midx] = i_ne_var;
      assert(i_ne_var<kMaxVarBins);
    }//i
#endif

#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    
    cout << "SKMC/FillSpecialTemplates >> Filling osc nuebar templates for " << util::AsString(s) << " sample" << endl;
  
    SuperkMode_t src_modes_nuebar[n] = {  kSKmNuEBarCCQE,    kSKmNuEBarCC1pi,     kSKmNuEBarCCcoh,    kSKmNuEBarCCoth    };   
    SuperkMode_t tgt_modes_nuebar[n] = {  kSKmOscNuEBarCCQE, kSKmOscNuEBarCC1pi,  kSKmOscNuEBarCCcoh, kSKmOscNuEBarCCoth };
    
    for(int i=0; i<n; i++) {
      
      SuperkMode_t tgt_m = tgt_modes_nuebar[i];
      int tgt_midx = arr_index->SuperkModeIdx(s,tgt_modes_nuebar[i]);
      int src_midx = arr_index->SuperkModeIdx(s,src_modes_nuebar[i]);
      
      TH2D * src_ori = this->Template(s,src_modes_nuebar[i],kSKMCtOri);
      TH2D * tgt_ori = this->Template(s,tgt_modes_nuebar[i],kSKMCtOri);
      TH2D * src_nom = this->Template(s,src_modes_nuebar[i],kSKMCtNom);
      TH2D * tgt_nom = this->Template(s,tgt_modes_nuebar[i],kSKMCtNom);
      
    // loop over Etrue
      for(int it = 1; it <= src_ori->GetYaxis()->GetNbins(); it++) {
	double Etrue = src_ori->GetYaxis()->GetBinCenter(it);
	// calc weight for weighting numu flux to nue flux
	// (note: we have used `nominal' with different meaning in the flux and in the template enumerations.
	// perhaps we should rename kFlxNominal -> kFlxOriginal, for consistency).
	double flux_nue_ori  = fNuFlux->SuperkFlux (s, kPdgNuEBar,  Etrue, fNom   );
	double flux_nue_nom  = fNuFlux->SuperkFlux (s, kPdgNuEBar,  Etrue, fExt ); 
	double flux_numu_ori = fNuFlux->SuperkFlux (s, kPdgNuMuBar, Etrue, fNom   );
	double flux_numu_nom = fNuFlux->SuperkFlux (s, kPdgNuMuBar, Etrue, fExt ); 
	double wght_nue2numu_ori = (flux_nue_ori != 0) ? flux_numu_ori / flux_nue_ori : 0;
	double wght_nue2numu_nom = (flux_nue_nom != 0) ? flux_numu_nom / flux_nue_nom : 0;
	// loop over Ereco
	for(int ir = 1; ir <= src_ori->GetXaxis()->GetNbins(); ir++) {
	  double Ereco = src_ori->GetXaxis()->GetBinCenter(ir);
	  // get nevents
	  double nev_ori = src_ori-> GetBinContent(ir,it);
	  double nev_nom = src_nom-> GetBinContent(ir,it);
	  // set number of events
	  double weighted_nev_ori = nev_ori*wght_nue2numu_ori;
	  double weighted_nev_nom = nev_nom*wght_nue2numu_nom;
	  tgt_ori->Fill(Ereco,Etrue, weighted_nev_ori);
	  tgt_nom->Fill(Ereco,Etrue, weighted_nev_nom);
	  
	  double nev_tgt_ori = tgt_ori->GetBinContent(ir,it);
	  double nev_tgt_nom = tgt_nom->GetBinContent(ir,it);
    
          fMCTemplArrayOri[is].GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_ori;
          fMCTemplArrayNom[is].GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
          fMCTemplArrayTun[is].GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
          fMCTemplArrayFin[is].GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
          
        }//ir
      }//it
      //one more loop, for the variable binning scheme.
      //just set the OscNuMu equal to the NuMu
      //loop over variable Etrue
      int i_ne_var = 0;

      for(int itv = 1; itv <= fMCVariableTempl[sidx][tgt_midx]->GetYaxis()->GetNbins(); itv++) {
	double EtrueVarMin = util::VariableMCEtrueBinEdges(s,tgt_m)[itv-1] + kOneMeV;
	double EtrueVarMax = util::VariableMCEtrueBinEdges(s,tgt_m)[itv]   - kOneMeV;

	int EtrueMinBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMin);
	int EtrueMaxBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMax);
	
	// loop over variable Ereco
	for(int irv = 1; irv <= fMCVariableTempl[sidx][tgt_midx]->GetXaxis()->GetNbins(); irv++) {
	  double nev = 0;

	  nev = fMCVariableTemplArray[sidx].GetEntry   (src_midx, itv-1, irv-1);
	  fMCVariableTemplArray[sidx]      .GetEntryRef(tgt_midx, itv-1, irv-1) = nev;
	  ////

	  
	  // Fill nonempty bin numbers.
	  if (nev < conv::kASmallNum) {
	    continue;
	  }
	  
	  fNonEmptyVariableBins[sidx][tgt_midx][i_ne_var] = conv::kBigIndex * irv + itv;
	  
	  double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
	  double ErecoVarMax = util::XsecErecoBinEdges(s,irv)   - kOneMeV;
	  
	  int ErecoMinBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMin);
	  int ErecoMaxBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMax);

	  int i_ne_fixed = 0;
	  // loop over fixed Etrue to count n nonempty
	  for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
	    // loop over fixed Ereco to count
	    for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
	      double nevf = this->Template(s,tgt_m,kSKMCtOri)->GetBinContent(ir,it);
	      if (nevf < conv::kASmallNum) {
		continue;
	      }
	      i_ne_fixed++;
	    }//ir
	  }//it
	  fNNonEmptyFixedBins[sidx][tgt_midx][i_ne_var] = i_ne_fixed;
	  fNonEmptyFixedBins[sidx][tgt_midx][i_ne_var] = new int[i_ne_fixed];
	  
	  i_ne_fixed = 0;
	  // loop over fixed Etrue to count n nonempty
	  for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
	    // loop over fixed Ereco to count
	    for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
	      double nevf = this->Template(s,tgt_m,kSKMCtOri)->GetBinContent(ir,it);
	      if (nevf < conv::kASmallNum) {
		continue;
	      }
	      fNonEmptyFixedBins[sidx][tgt_midx][i_ne_var][i_ne_fixed] = conv::kBigIndex * ir + it;
	      i_ne_fixed++;
	    }//ir
	  }//it
	  i_ne_var++;
	}//irv
      }//itv
      fNNonEmptyVariableBins[sidx][tgt_midx] = i_ne_var;
      assert(i_ne_var<kMaxVarBins);
      
    }//i
    
#endif// __USE_OSC_NUMU_MC_TEMPLATES__


  }
}
//________________________________________________________________________________________
void SuperkMC::CalculateMCStatError(void)
{
  // set MC frac err in each bin
  for(int is = 0; is < kns; is++) 
  {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    fMCStatFracErr[is].SetAll(0.);
    for(int im=0; im < knm; im++) 
    {
      SuperkMode_t m = kSKReactionModes[is][im];
      if(m == kSKmUndefined) continue;

      // loop over Etrue
      for(int it = 1; it <= this->Template(s,m,kSKMCtOri)->GetYaxis()->GetNbins(); it++) {
        // loop over Ereco
        for(int ir = 1; ir <= this->Template(s,m,kSKMCtOri)->GetXaxis()->GetNbins(); ir++) {
            double nev  = (double) this->Template(s,m,kSKMCtOri)->GetBinContent(ir,it);
            nev = TMath::Max(0.,nev);
            double dnev = TMath::Sqrt(nev);
            double fracerr = (nev > 0) ? dnev/nev : 0;

            fMCStatFracErr[is].GetEntryRef(im, it-1, ir-1) = fracerr;
        }//ir
      }//it
    }//m
  }//s

}
//________________________________________________________________________________________
void SuperkMC::ApplyMCStatFluctationToTemplates(bool shape_only)
{
  TRandom3 & rndgen = Rnd::Instance()->RndToyMCstat();

  for(int is = 0; is < kns; is++) 
  {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     for(int im=0; im < knm; im++) 
     {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       bool skip = this->SkipTemplate(m);
       if(skip) continue;
       double nev_old = 0;
       double nev_new = 0;

       // fluctuate bin contents
       int n_ereco_bins = kNErecoBins[is];
       int n_etrue_bins = kNEtrueBins[is];
       
       for(int ir = 1; ir <= n_ereco_bins; ir++) {
          for(int it = 1; it <= n_etrue_bins; it++) {
            double nev_old_bin  = fMCTemplArrayFin[is].GetEntryRef(im, it-1, ir-1); 
            double frac_error   = fMCStatFracErr  [is].GetEntryRef(im, it-1, ir-1);
            double error        = frac_error * nev_old_bin;
            double nev_new_bin  = rndgen.Gaus(nev_old_bin, error);
            nev_new_bin = TMath::Max(0.,nev_new_bin);
            nev_old += nev_old_bin;
            nev_new += nev_new_bin;
            fMCTemplArrayFin[is].GetEntryRef(im, it-1, ir-1) = nev_new_bin;
            
          }//it
       }//ir
       
       if(shape_only) {
          double scale = (nev_new > 0.) ? nev_old/nev_new : 0.;
          for(int ir = 1; ir <= n_ereco_bins; ir++) {
             for(int it = 1; it <= n_etrue_bins; it++) {
                double nev  = fMCTemplArrayFin[is].GetEntryRef(im, it-1, ir-1); 
                fMCTemplArrayFin[is].GetEntryRef(im, it-1, ir-1) = nev*scale;
             }//it
          }//ir
       }//shape_only?

     }//m
  }//s

}

//________________________________________________________________________________________
TH1D* SuperkMC::CurrErecoSpectrum(SuperkSample_t s) const                               
{ 
  const double * arr = CurrErecoSpectrumArray(s);
  if (!arr) return 0;
  
  string name  = Form("curr_ereco_%s"             ,AsStringCompact(s).c_str() );
  string title = Form("Current Ereco Spectrum: %s",       AsString(s).c_str() );

  int nbins = util::NumErecoBins(s);
  TH1D * h = new TH1D(name.c_str(), title.c_str(), nbins, util::ErecoBinEdges(s) );
  util::FormatErecoHist(h);
  for (int i = 0 ; i < nbins ; i++ ) {
    h->SetBinContent(i+1,arr[i]);
  }
  return h;
}
//________________________________________________________________________________________
TH1D* SuperkMC::CurrErecoSpectrum(SuperkSample_t s, SuperkMode_t m) const 
{ 
  const double * arr = CurrErecoSpectrumArray(s,m);
  if (!arr) return 0;
  
  string name  = Form("curr_ereco_%s_%s"             ,AsStringCompact(s).c_str(), AsStringCompact(m).c_str());
  string title = Form("Current Ereco Spectrum: %s %s",AsString(s).c_str()       , AsString(m).c_str()       );
  
  int nbins = util::NumErecoBins(s);
  TH1D * h = new TH1D(name.c_str(), title.c_str(), nbins, util::ErecoBinEdges(s) );
  util::FormatErecoHist(h);
  util::SetDefaultSKSpectrumStyle(s,m,h);
  for (int i = 0 ; i < nbins ; i++ ) {
    h->SetBinContent(i+1,arr[i]);
  }
  return h;
}
//________________________________________________________________________________________
const double* SuperkMC::CurrErecoSpectrumArray(SuperkSample_t s) const  
{ 
  int sidx = arr_index->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  
  return fMCErecoSpectrumTotArray[sidx].GetArray();
}
//________________________________________________________________________________________
double* SuperkMC::CurrErecoSpectrumArray(SuperkSample_t s)  
{ 
  int sidx = arr_index->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  
  return fMCErecoSpectrumTotArray[sidx].GetArray();
}
//________________________________________________________________________________________
const double* SuperkMC::CurrErecoSpectrumArray(SuperkSample_t s, SuperkMode_t m) const
{ 
  int sidx = arr_index->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = arr_index->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }
  
  return &(fMCErecoSpectrumArray[sidx].GetEntryRef(midx,0));
}
//________________________________________________________________________________________
double* SuperkMC::CurrErecoSpectrumArray(SuperkSample_t s, SuperkMode_t m)
{ 
  int sidx = arr_index->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = arr_index->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }

  return &(fMCErecoSpectrumArray[sidx].GetEntryRef(midx,0));
}
//________________________________________________________________________________________
TH2D * SuperkMC::Template(SuperkSample_t s, SuperkMode_t m, SKMCTemplateType_t type) const
{
  int sidx = arr_index->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = arr_index->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }            

  TH2D * h = 0;

  switch(type) {
   case (kSKMCtOri) : 
     h = fMCTemplOri[sidx][midx]; 
     break;
   case (kSKMCtNom) : 
     h = fMCTemplNom[sidx][midx]; 
     break;
   default:
     return 0;
  }
  assert(h);

  return h;
}
//________________________________________________________________________________________
TH2D * SuperkMC::VarTemplate(SuperkSample_t s, SuperkMode_t m) const
{
  int sidx = arr_index->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = arr_index->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }            

  TH2D * h = 0; 
  h = fMCVariableTempl[sidx][midx];
  assert(h);

  return h;
}
//________________________________________________________________________________________
void SuperkMC::ApplyBANFFTuningToTemplates(void)
{
// Take BANFF tuning into account (reweight MC to BANFF best-fit parameters)
// The weights returned by NuFlux::SuperkFluxFinTuningWght() are relative to the flux
// tuned based on external (NA61,...) data, so make sure that you run 
// SuperkMC::ApplyFluxExtTuningToTemplates() before running this method.
// Source template type: kSKMCtNom
// Target template type: kSKMCtTun
//

#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> Applying BANFF flux and cross-section tuning." << endl;
#endif

  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    // loop over sample modes
    for(int im = 0; im < knm; im++) { 
      SuperkMode_t m = kSKReactionModes[is][im];
      if(m == kSKmUndefined) continue;
      // consider template?
      bool skip = this->SkipTemplate(m);
      if(skip) continue;
      
#ifdef __show_banfftun_syst_var_wght_mesg_debug__ 
      cout << "Exp/SKMC >> Applying BANFF tuning to sample, mode: " << AsString(s) << "\t" << AsString(m) << endl;
#endif

      //
      // There are 600 50MeV Etrue energy bins in the templates, but
      // the flux/xsec systematics only have up to ~10 Etrue bins (same as flux binning), so
      // we can call NuXSec::FinTuningWght() and NuFlux::SuperkFluxFinTuningWght()  
      // a factor of 10 less often by looping over the kXsEtrueBinEdges.
      // Note that the Ereco binning is identical between XSec and SuperkMC templates
      //
      for (int ivb = 0 ; ivb < fNNonEmptyVariableBins[is][im] ; ivb++ ) {
	int itv =  fNonEmptyVariableBins[is][im][ivb] % conv::kBigIndex;
	int irv = (fNonEmptyVariableBins[is][im][ivb] - itv)/conv::kBigIndex;

	// Get the min/max energy of the current variable flux Etrue bin
	double EtrueVarMin = util::VariableMCEtrueBinEdges(s,m)[itv-1] + kOneMeV;

	// Plug-in BANFF flux systematics
        double flux_wght = fNuFlux->SuperkFluxFinTuningWght(s,m,EtrueVarMin);

	//skip if the big bin is empty
	if(fMCVariableTemplArray[is].GetEntryRef(im, itv-1, irv-1) <= 0) continue;
	
	double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
	
	// Plug-in BANFF & other xsec systematics 
	double xsec_wght = fNuXSec->FinTuningWght(s,m,ErecoVarMin,EtrueVarMin); 
      
	// Combined weight
	double wght = flux_wght * xsec_wght;

#ifdef BANFFWEIGHTSAREONE
	wght = 1;
#endif
	for (int ifb = 0 ; ifb < fNNonEmptyFixedBins[is][im][ivb] ; ifb++ ) {
	  int it =  fNonEmptyFixedBins[is][im][ivb][ifb] % conv::kBigIndex;
	  int ir = (fNonEmptyFixedBins[is][im][ivb][ifb] - it)/conv::kBigIndex;
    
	  //check for empty bins
	  if(fMCTemplArrayOri[is].GetEntryRef(im, it-1, ir-1) <= 0) continue;
	  ////
	  //double Etrue = util::EtrueBinEdges(s,it-1) + kOneMeV;
	  
	  // get number of events in current bin of source MC template
	  double nev = fMCTemplArrayNom[is].GetEntry(im, it-1, ir-1);
	  
	  nev = TMath::Max(0.,nev);
	  double nev_tuned = nev * wght;
	  
#ifdef NO_EVENTS_CAN_BE_NEGATIVE
	  if(nev_tuned < 0) {
#ifdef __show_negative_events_mesg_debug__
	    cout << "Exp/SKMC/ApplyBANFFTuningToTemplates >> set negative number of events " << nev_tuned
	      << " to " << 0 << "\t" << " Ereco: " << Ereco << " Etrue: " << util::EtrueBinEdges(s,it-1) + kOneMeV << endl;
#endif
	    nev_tuned = 0;
	  }
#endif //NO_EVENTS_CAN_BE_NEGATIVE
    
    fMCTemplArrayTun[is].GetEntryRef(im, it-1, ir-1) = nev_tuned;
    fMCTemplArrayFin[is].GetEntryRef(im, it-1, ir-1) = nev_tuned;
    
#ifdef __show_banfftun_syst_var_wght_mesg_debug__ 
	  // Get the energy of the current bin
	  double Etrue = util::EtrueBinEdges(s,it-1) + kOneMeV;
	  double Ereco = util::ErecoBinEdges(s,it-1) + kOneMeV;
	  if( flux_wght  != 1. ||
	      xsec_wght  != 1. )
	  {
	    cout << "Exp/SKMC/SystVarWght >> "
	      << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m)
	      << ", Ereco = " << Ereco << " GeV"
	      << ", Etrue = " << Etrue << " GeV => Weights: "
	      << ", Flux = "  << flux_wght
	      << ", XSec = "  << xsec_wght
	      << " (total = " << wght << ")"
	      << endl;
	  }
#endif //__show_banfftun_syst_var_wght_mesg_debug__ 

	}//ifb
      }//ivb
    }//im
  }//is

#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC/ApplyBANFFTuningToTemplates >> Complete" << endl;
#endif
}
//________________________________________________________________________________________
void SuperkMC::UpdateTemplates(void)
{
#ifdef __show_update_mesg_info__
  cout << "Exp/SKMC >> Updating MC templates...." << endl;
#endif
  
  this->ApplySystematicsToTemplates   ();
  this->ApplyPhysHypothesisToTemplates();

#ifdef __INCLUDE_MC_STAT_FLUCTUATIONS__
#ifdef __show_update_mesg_info__
  cout << "Exp/SKMC >> ** Note ** Fluctuating templates to take into account the effect of finite MC statistics...." << endl;
#endif
  this->ApplyMCStatFluctationToTemplates(true);
#endif
}
//________________________________________________________________________________________
void SuperkMC::ApplySystematicsToTemplates(void) 
{
// Apply systematic variations to MC templates
//
  this->ApplySystematicsToTemplates_Multiplicative   ();
}
//________________________________________________________________________________________
void SuperkMC::ApplySystematicsToTemplates_Multiplicative(void) 
{
  // All multiplicative systematics are handled in the following loop.
  //
#ifdef __show_update_mesg_info__
  cout << "Exp/SKMC >> Applying multiplicative systematics" << endl;
#endif

  // Loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    // Loop over sample modes
    for(int im = 0; im < knm; im++) { 
      SuperkMode_t m = kSKReactionModes[is][im];
      if(m == kSKmUndefined) continue;
      // Consider template?
      bool skip = this->SkipTemplate(m);
      if(skip) continue;
      
      const double * variableMCEtrueBinEdges = util::VariableMCEtrueBinEdges(s,m);
      
      // Recast VALORArray blocks so they can be used addressed like normal arrays (fast).
      // Hard to read, so should only be used in code for which speed is essential.
      int n_ereco_bins = kNErecoBins[is];
      int n_etrue_bins = kNEtrueBins[is];
      typedef double (templtype)[n_etrue_bins][n_ereco_bins];
      const templtype & templTun = *reinterpret_cast<const templtype * >(&(fMCTemplArrayTun[is].GetEntryRef(im,0,0)));
      templtype       & templFin = *reinterpret_cast<templtype       * >(&(fMCTemplArrayFin[is].GetEntryRef(im,0,0)));
      
      //~ const double(&templTun)[n_etrue_bins][n_ereco_bins] = *reinterpret_cast< const double (*) [n_etrue_bins][n_ereco_bins]> (&(fMCTemplArrayTun[is].GetEntryRef(im,0,0)));
      //~ double      (&templFin)[n_etrue_bins][n_ereco_bins] = *reinterpret_cast< double       (*) [n_etrue_bins][n_ereco_bins]> (&(fMCTemplArrayFin[is].GetEntryRef(im,0,0)));

      //
      // There are 600 50MeV Etrue energy bins in the templates, but
      // the flux/xsec systematics only have up to ~10 Etrue bins (same as flux binning), so
      // we can call NuXSec::FinTuningWght() and NuFlux::SuperkFluxFinTuningWght()  
      // a factor of 10 less often by looping over the kXsEtrueBinEdges.
      // Note that the Ereco binning is identical between XSec and SuperkMC templates
      //
      // Loop over variable flux Etrue bins
      int last_itv = -1;
      double flux_wght = -1;
      for (int ivb = 0 ; ivb < fNNonEmptyVariableBins[is][im] ; ivb++ ) {
	int itv =  fNonEmptyVariableBins[is][im][ivb] % conv::kBigIndex;
	int irv = (fNonEmptyVariableBins[is][im][ivb] - itv)/conv::kBigIndex;

	double EtrueVarMin = variableMCEtrueBinEdges[itv-1] + kOneMeV;

	// Plug-in BANFF flux systematics
	if (itv!=last_itv) {
	  flux_wght = fNuFlux->SuperkFluxSystWght(s,m,EtrueVarMin);
	  last_itv = itv;
	}
	
	//assert(fMCVariableTemplArray[is].GetEntryRef(im, itv-1, irv-1) > 0);
	
	double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;

	// Plug-in BANFF & other xsec systematics 
	double xsec_wght = fNuXSec->SystWght(s,m,ErecoVarMin,EtrueVarMin); 
  
	////////// HERE IT CRASHES //////////

	// Plug-in SK efficiency + FSI uncertainties
	double skfsi_wght = fSKDetFSI->SystWght(s,m,ErecoVarMin,EtrueVarMin);
  
	// if(s==kSKsFCFVMultiRingMuLike || s==kSKsFCFVMultiRingELike){
	//   exit(1);
	// }

	
	// Combined weight
	double wght = skfsi_wght * flux_wght * xsec_wght;
  
#ifdef SYSTWEIGHTSAREONE
	wght = 1;
#endif
	/*
	   if(util::IsMEC(m))
	   cout << "SK: " << skfsi_wght << "\tFlux: " << flux_wght << "\tXSec: " << xsec_wght << "\tTotal: " << wght << endl;
	 */
	// loop over 50MeV Ereco
	for (int ifb = 0 ; ifb < fNNonEmptyFixedBins[is][im][ivb] ; ifb++ ) {
	  int it =  fNonEmptyFixedBins[is][im][ivb][ifb] % conv::kBigIndex;
	  int ir = (fNonEmptyFixedBins[is][im][ivb][ifb] - it)/conv::kBigIndex;
	  //check for empty bins
	  
	  //assert(fMCTemplOri[is].GetEntryRef(im, itv-1, irv-1) > 0);
    
	  ////
	  //double Etrue = util::ErecoBinEdges(s)(s,it-1) + kOneMeV;
	  
	  // get number of events in current bin of source MC template
	  double nev = templTun[it-1][ir-1];
	  
	  nev = TMath::Max(0.,nev);
	  double nev_systm = nev * wght;

#ifdef __show_syst_var_wght_mesg_debug__ 
	    // Get the energy of the current bin
	  double Etrue = util::EtrueBinEdges(s,it-1) + kOneMeV;
	  double Ereco = util::ErecoBinEdges(s,it-1) + kOneMeV;
	  if( skfsi_wght != 1. ||
	      flux_wght  != 1. ||
	      xsec_wght  != 1. )
	  {
	    cout << "Exp/SKMC/SystVarWght >> "
	      << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m)
	      << ", Ereco = " << Ereco << " GeV"
	      << ", Etrue = " << Etrue << " GeV => Weights: "
	      << "SK+FSI = "  << skfsi_wght
	      << ", Flux = "  << flux_wght
	      << ", XSec = "  << xsec_wght
	      << " (total = " << wght << ")"
	      << endl;
	  }
#endif//__show_syst_var_wght_mesg_debug__ 

#ifdef NO_EVENTS_CAN_BE_NEGATIVE
	  if(nev_systm < 0) {
#ifdef __show_negative_events_mesg_debug__
	    cout << "Exp/SKMC/ApplySystematicsToTemplates_Multiplicative >> set negative number of events " << nev_systm
	      << " to " << 0 << "\t" << " Ereco: " << Ereco << endl;
#endif//__show_negative_events_mesg_debug__
	    nev_systm = 0;
	  }
#endif// NO_EVENTS_CAN_BE_NEGATIVE
	  // Update target MC template
	  
	  templFin[it-1][ir-1] = nev_systm;
	  
	}//ifb
      }//ivb
    }//im
      
  }//is

#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtSys, "mcnorm_fluxtuned_banfftuned_syst_multi");
#endif
}
//________________________________________________________________________________________
bool SuperkMC::ApplyEnergyScaleSystematicToTotalSpectrum(void) 
{
//
// Now that all other systematics and physics hypotheses (that depend on Etrue) have been 
// applied, apply the reconstructed energy scale systematic (that depends on Ereco).
// This systematic is applied by scaling the Ereco bin edges and assuming uniform 
// distribution of events within each bin.
#ifdef __show_update_mesg_info__
  cout << "Exp/SKMC >> Applying non-multiplicative systematics to total spectrum" << endl;
#endif

  bool applied_negevents_condition = false;
  
  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    
    double escale = fSKDetFSI->EnergyScale(s);
    
    const int n_ereco_bins = kNErecoBins[is];
    
    double * events_per_bin = this->CurrErecoSpectrumArray(s);
    
    const double * bin_edges = util::ErecoBinEdges(s);
    
    double new_events_per_bin [n_ereco_bins+1];// With overflow bin
    double scaled_bin_edges   [n_ereco_bins+1];
    double unscaled_bin_edges [n_ereco_bins+2];// With overflow bin
    
    // Calculate bin information
    double total_events = 0;
    for (int ibin = 0 ; ibin <= n_ereco_bins; ibin++) {
      if (ibin  != n_ereco_bins) {
        total_events += events_per_bin[ibin];
      }
      unscaled_bin_edges[ibin] = bin_edges[ibin];
      scaled_bin_edges  [ibin] = bin_edges[ibin] * escale;
      new_events_per_bin[ibin] = 0;
    }
    
    unscaled_bin_edges[n_ereco_bins+1] = 1E4; //Overflow bin (no underflow unless we start binning from something other than zero).
    
    if (total_events < kASmallNum) continue;
    
    double overlap = 0;
        
    for (int ibin_s = 0 ; ibin_s < n_ereco_bins ; ibin_s++) { // Loop over scaled bins
      double bin_min_s = scaled_bin_edges[ibin_s];
      double bin_max_s = scaled_bin_edges[ibin_s+1];
      for (int ibin_u = 0 ; ibin_u < n_ereco_bins+1 ; ibin_u++) {
        double bin_min_u = unscaled_bin_edges[ibin_u];
        double bin_max_u = unscaled_bin_edges[ibin_u+1];
        
        if     (bin_min_u >= bin_min_s && bin_max_u <= bin_max_s) overlap = bin_max_u - bin_min_u; // Entirely subsumed unscaled bin
        else if(bin_min_s >= bin_min_u && bin_max_s <= bin_max_u) overlap = bin_max_s - bin_min_s; // Entirely subsumed scaled bin
        else if(bin_max_u >  bin_max_s && bin_min_u >= bin_min_s) overlap = bin_max_s - bin_min_u; // Scaled bin extends low
        else if(bin_min_u <  bin_min_s && bin_max_u <= bin_max_s) overlap = bin_max_u - bin_min_s; // Scaled bin extends high
        else if(bin_min_u >  bin_max_s) break ; //Unscaled bin entirely larger than scaled bin, stop checking higher unscaled
        else overlap = -999; // No crossover
        
        if (overlap <= 0) continue;
        
        double new_events = events_per_bin[ibin_s] * overlap / (scaled_bin_edges[ibin_s+1]-scaled_bin_edges[ibin_s]);
        new_events_per_bin[ibin_u] += new_events;
      }
    }
    // Move events back into main array.
    for (int i = 0 ; i < n_ereco_bins ; i++) {
      events_per_bin[i] = new_events_per_bin[i];
    }
    events_per_bin[n_ereco_bins-1] += new_events_per_bin[n_ereco_bins]; //Put events from overflow bin back into top bin
    
#ifndef NDEBUG
    double new_total_events = 0;
    for (int i = 0 ; i < n_ereco_bins ; i++) {
      new_total_events+=events_per_bin[i];
    }
    if ( fabs((new_total_events/total_events)-1) > conv::kASmallNum ) {
      cerr<<"SuperkMC/ApplyEnergyScaleSystematicToTotalSpectrum - failed normalisation check for sample "<<AsString(s)<<endl;
      cerr<<"Events: "<<total_events<<" ==> "<<new_total_events<<"\t"<<new_total_events/total_events
          <<"\t"<<fabs((new_total_events/total_events)-1)<<"\t at scale "<<escale<<endl;
      assert(!"Broken");
    }
#endif //NDEBUG

  }//is
  return applied_negevents_condition;
}
//________________________________________________________________________________________
void SuperkMC::ApplyPhysHypothesisToTemplates(void)
{
// Apply physics hypothesis {oscillation (2-flavour / 3-flavour), decay or 
// decoherence} to MC templates
//

#ifdef __show_update_mesg_info__
  cout << "Exp/SKMC >> Applying physics hypothesis" << endl;
#endif

  PhysHypothesis_t hyp = fPhysHypothesis.Hypothesis();

  //
  // 3-flavour oscillations
  //
  if(hyp == kHypOsc3F)
    {
    // set oscillation parameters to the oscillation probability calculation code
    double mixing12 = fPhysHypothesis.CurrValue(kPhypMixing12);
    double mixing13 = fPhysHypothesis.CurrValue(kPhypMixing13);
    double mixing23 = fPhysHypothesis.CurrValue(kPhypMixing23);
    //    double dm12sq   = fPhysHypothesis.CurrValue(kPhypDm12sq);
    //double dm23sq   = fPhysHypothesis.CurrValue(kPhypDm23sq);
    double dmsq     = fPhysHypothesis.CurrValue(kPhypDm12sq);
    double Dmsq     = fPhysHypothesis.CurrValue(kPhypDm23sq);
    double delta_cp = fPhysHypothesis.CurrValue(kPhypCPV    );
#ifdef __show_hypothesis__
    cout << "**  Exp/SKMC >> " << fPhysHypothesis;
#endif
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    NuOscDmsqScheme_t  dmsqscheme = Opt::Instance()->gOscDmsqScheme;
    
    double p_osc[3][3];
    double p_osc_bar[3][3];

    // mixing params are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|                                                                                             
    if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsq2Theta)
      {
	fNuOsc3Flavour.SetOscParameters(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
      }
    // mixing params are: sin^2(theta)                                                                                                                                         
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|                                                                                              
    else if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsqTheta)
      {
	fNuOsc3Flavour.SetOscParametersSinSqTAll(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
      }
    // mixing params are: sin^2(2*theta)                                                                                                                                        
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|                                                                        
    else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsq2Theta)
      {
	fNuOsc3Flavour.SetOscParametersSSq2TDmAtmSol(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
      }
    // mixing params are: sin^2(theta)                                                                                                                                             
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|                                                                         
    else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsqTheta)
      {
	fNuOsc3Flavour.SetOscParametersDmsqAtmSol(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
      }
    // mixing paras are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
    else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsq2Theta)
      {
	fNuOsc3Flavour.SetOscParametersSSq2TDmsq13(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
      }
    // mixing paras are: sin^2(theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
    else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsqTheta)
      {
	fNuOsc3Flavour.SetOscParametersDmsq13(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
      }
    //unknown
    else {
      cerr << "**  Exp/SKMC >> Can't handle choice of "
	   << "\n oscillation mixing params: " << util::AsString(omp)
	   << "\n and independent squared mass-splittings: " << util::AsString(dmsqscheme)
	   << endl;
      exit(1);
    }

    // loop over samples
    for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;

     int n_ereco_bins = kNErecoBins[is];
     int n_etrue_bins = kNEtrueBins[is];

     //assert(this -> Template(s,kSKmNuMuCCQE,kSKMCtSys));
     //TAxis * ereco_axis = this->Template(s,kSKmNuMuCCQE,kSKMCtSys)->GetXaxis(); 
     //TAxis * etrue_axis = this->Template(s,kSKmNuMuCCQE,kSKMCtSys)->GetYaxis(); 

     // loop over Etrue
     for(int it = 1; it <= n_etrue_bins; it++) {

        // calculate oscillation probabilities (incl. matter effects) for the current Etrue value
       double Etrue = (util::EtrueBinEdges(s,it-1) + util::EtrueBinEdges(s,it) ) / 2;

	// skip unnecessary calculations
	//bool bin_is_not_empty = this->IsThereAnEvent(s, Etrue);
	//if(!bin_is_not_empty) continue;
	
	fNuOsc3Flavour.CalcMatterOscProbAll(false, Etrue, p_osc    );
	fNuOsc3Flavour.CalcMatterOscProbAll(true , Etrue, p_osc_bar);
	
	double p_numu_numu       = p_osc    [1][1];
	double p_numubar_numubar = p_osc_bar[1][1];
	double p_nue_nue         = p_osc    [0][0];
	double p_nuebar_nuebar   = p_osc_bar[0][0];
	double p_numu_nue        = p_osc    [1][0];
	double p_numubar_nuebar  = p_osc_bar[1][0];
	double p_nue_numu        = p_osc    [0][1];
	double p_nuebar_numubar  = p_osc_bar[0][1];
	
#ifdef __show_hypothesis_probabilities__
        cout << "Exp/SKMC >>  E = " << Etrue << " GeV" << endl;	
        cout << "Exp/SKMC >>  * P (numu     -> numu   ) = " << p_numu_numu       << endl;
        cout << "Exp/SKMC >>  * P (numubar  -> numubar) = " << p_numubar_numubar << endl;
        cout << "Exp/SKMC >>  * P (nue      -> nue    ) = " << p_nue_nue         << endl;
        cout << "Exp/SKMC >>  * P (nuebar   -> nuebar ) = " << p_nuebar_nuebar   << endl;
        cout << "Exp/SKMC >>  * P (numu     -> nue    ) = " << p_numu_nue        << endl;
        cout << "Exp/SKMC >>  * P (nue      -> numu   ) = " << p_nue_numu      << endl;
#endif


	// loop over modes
	for(int im=0; im < knm; im++) {
	  SuperkMode_t m = kSKReactionModes[is][im];
	  if(m == kSKmUndefined) continue;
	  double * array = &(fMCTemplArrayFin[is].GetEntryRef(im,it-1,0));
	  // osc. probability
	  double posc = 1.;
	  if(kIsCC[is][im])
	  {
	    if     ( util::IsNuMu(m)      ) { posc = p_numu_numu;       }
	    else if( util::IsNuMuBar(m)   ) { posc = p_numubar_numubar; }
	    else if( util::IsNuE(m)       ) { posc = p_nue_nue;         } 
	    else if( util::IsNuEBar(m)    ) { posc = p_nuebar_nuebar;   }
	    else if( util::IsOscNuMu(m)   ) { posc = p_nue_numu;        }
	    else if( util::IsOscNuMuBar(m)) { posc = p_nuebar_numubar;  }
	    else if( util::IsOscNuE(m)    ) { posc = p_numu_nue;        }
	    else if( util::IsOscNuEBar(m) ) { posc = p_numubar_nuebar;  }
	  }
	  // loop over Ereco
	  for(int ir = 1; ir <= n_ereco_bins ; ir++) {
	    // number of events in (Ereco,Etrue) bin (systematics applied)
	    double nev = 0;
	    
	    nev = array[ir-1];
	    ////
	    
#ifdef PHYSHYPWEIGHTSAREONE
	    posc = 1;
#ifdef PHYSHYPWEIGHTAREONEANDNUEISZERO
	    if(util::IsOscNuE(m) || util::IsOscNuEBar(m) || util::IsOscNuMuBar(m) || util::IsOscNuMu(m) ) posc = 0;
#endif
#endif
	    // number of events in (Ereco,Etrue) bin (oscillations + systematics applied)
	    double nev_osc = nev * posc;
	    
	    array[ir-1] = nev_osc;
	    ////
	  }//ir
	}//im
     }//it
   }//is
  }
  //
  // 3-flavour oscillations, separate nu and nubar
  //
  if(hyp == kHypOsc3FBar)
    {
    // set oscillation parameters to the oscillation probability calculation code
    double mixing12 = fPhysHypothesis.CurrValue(kPhypMixing12);
    double mixing13 = fPhysHypothesis.CurrValue(kPhypMixing13);
    double mixing23 = fPhysHypothesis.CurrValue(kPhypMixing23);
    //double dm12sq   = fPhysHypothesis.CurrValue(kPhypDm12sq);
    //double dm23sq   = fPhysHypothesis.CurrValue(kPhypDm23sq);
    double dmsq     = fPhysHypothesis.CurrValue(kPhypDm12sq);
    double Dmsq     = fPhysHypothesis.CurrValue(kPhypDm23sq);
    double delta_cp = fPhysHypothesis.CurrValue(kPhypCPV    );
    
    double mixing12bar = fPhysHypothesis.CurrValue(kPhypMixing12Bar);
    double mixing13bar = fPhysHypothesis.CurrValue(kPhypMixing13Bar);
    double mixing23bar = fPhysHypothesis.CurrValue(kPhypMixing23Bar);
    //double dm12sqbar   = fPhysHypothesis.CurrValue(kPhypDm12sqBar);
    //double dm23sqbar   = fPhysHypothesis.CurrValue(kPhypDm23sqBar);
    double dmsqbar     = fPhysHypothesis.CurrValue(kPhypDm12sqBar);
    double Dmsqbar     = fPhysHypothesis.CurrValue(kPhypDm23sqBar);
    double delta_cpbar = fPhysHypothesis.CurrValue(kPhypCPVBar );
#ifdef __show_hypothesis__
    cout << "**  Exp/SKMC >> " << fPhysHypothesis;
#endif
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    NuOscDmsqScheme_t  dmsqscheme = Opt::Instance()->gOscDmsqScheme;
    
    double p_osc[3][3];
    double p_osc_bar[3][3];

    // mixing params are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|                                                                                             
    if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsq2Theta)
    {
      fNuOsc3Flavour.SetOscParameters(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
      fNuOsc3FlavourBar.SetOscParameters(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing params are: sin^2(theta)                                                                                                                                         
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|                                                                                              
    else if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsqTheta)
    {
      fNuOsc3Flavour.SetOscParametersSinSqTAll(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
      fNuOsc3FlavourBar.SetOscParametersSinSqTAll(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing params are: sin^2(2*theta)                                                                                                                                        
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|                                                                        
    else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsq2Theta)
    {
      fNuOsc3Flavour.SetOscParametersSSq2TDmAtmSol(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
      fNuOsc3FlavourBar.SetOscParametersSSq2TDmAtmSol(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing params are: sin^2(theta)                                                                                                                                             
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|                                                                         
    else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsqTheta)
    {
      fNuOsc3Flavour.SetOscParametersDmsqAtmSol(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
      fNuOsc3FlavourBar.SetOscParametersDmsqAtmSol(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing paras are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
    else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsq2Theta)
    {
      fNuOsc3Flavour.SetOscParametersSSq2TDmsq13(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
      fNuOsc3FlavourBar.SetOscParametersSSq2TDmsq13(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing paras are: sin^2(theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
    else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsqTheta)
    {
      fNuOsc3Flavour.SetOscParametersDmsq13(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
      fNuOsc3FlavourBar.SetOscParametersDmsq13(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    //unknown
    else {
      cerr << "**  Exp/SKMC >> Can't handle choice of "
           << "\n oscillation mixing params: " << util::AsString(omp)
           << "\n and independent squared mass-splittings: " << util::AsString(dmsqscheme)
           << endl;
      exit(1);
    }

    // loop over samples
    for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;

     int n_ereco_bins = kNErecoBins[is];
     int n_etrue_bins = kNEtrueBins[is];

     //assert(this -> Template(s,kSKmNuMuCCQE,kSKMCtSys));
     //TAxis * ereco_axis = this->Template(s,kSKmNuMuCCQE,kSKMCtSys)->GetXaxis(); 
     //TAxis * etrue_axis = this->Template(s,kSKmNuMuCCQE,kSKMCtSys)->GetYaxis(); 

     // loop over Etrue
     for(int it = 1; it <= n_etrue_bins; it++) {

        // calculate oscillation probabilities (incl. matter effects) for the current Etrue value
        double Etrue = (util::EtrueBinEdges(s,it-1) + util::EtrueBinEdges(s,it) ) / 2;

        // skip unnecessary calculations
        //bool bin_is_not_empty = this->IsThereAnEvent(s, Etrue);
        //if(!bin_is_not_empty) continue;
        
        fNuOsc3Flavour   .CalcMatterOscProbAll(false, Etrue, p_osc    );
        fNuOsc3FlavourBar.CalcMatterOscProbAll(true , Etrue, p_osc_bar);
        
        double p_numu_numu       = p_osc    [1][1];
        double p_numubar_numubar = p_osc_bar[1][1];
        double p_nue_nue         = p_osc    [0][0];
        double p_nuebar_nuebar   = p_osc_bar[0][0];
        double p_numu_nue        = p_osc    [1][0];
        double p_numubar_nuebar  = p_osc_bar[1][0];
        double p_nue_numu        = p_osc    [0][1];
        double p_nuebar_numubar  = p_osc_bar[0][1];
        
#ifdef __show_hypothesis_probabilities__
        cout << "Exp/SKMC >>  E = " << Etrue << " GeV" << endl;        
        cout << "Exp/SKMC >>  * P (numu     -> numu   ) = " << p_numu_numu       << endl;
        cout << "Exp/SKMC >>  * P (numubar  -> numubar) = " << p_numubar_numubar << endl;
        cout << "Exp/SKMC >>  * P (nue      -> nue    ) = " << p_nue_nue         << endl;
        cout << "Exp/SKMC >>  * P (nuebar   -> nuebar ) = " << p_nuebar_nuebar   << endl;
        cout << "Exp/SKMC >>  * P (numu     -> nue    ) = " << p_numu_nue        << endl;
        cout << "Exp/SKMC >>  * P (nue      -> numu   ) = " << p_nue_numu      << endl;
#endif


        // loop over modes
        for(int im=0; im < knm; im++) {
          SuperkMode_t m = kSKReactionModes[is][im];
          if(m == kSKmUndefined) continue;
          double * array = &(fMCTemplArrayFin[is].GetEntryRef(im,it-1,0));
          // osc. probability
          double posc = 1.;
          if(kIsCC[is][im])
          {
            if     ( util::IsNuMu(m)      ) { posc = p_numu_numu;       }
            else if( util::IsNuMuBar(m)   ) { posc = p_numubar_numubar; }
            else if( util::IsNuE(m)       ) { posc = p_nue_nue;         } 
            else if( util::IsNuEBar(m)    ) { posc = p_nuebar_nuebar;   }
            else if( util::IsOscNuMu(m)   ) { posc = p_nue_numu;        }
            else if( util::IsOscNuMuBar(m)) { posc = p_nuebar_numubar;  }
            else if( util::IsOscNuE(m)    ) { posc = p_numu_nue;        }
            else if( util::IsOscNuEBar(m) ) { posc = p_numubar_nuebar;  }
          }
          // loop over Ereco
          for(int ir = 1; ir <= n_ereco_bins ; ir++) {
            // number of events in (Ereco,Etrue) bin (systematics applied)
            double nev = 0;
            
            nev = array[ir-1];
            ////
            
#ifdef PHYSHYPWEIGHTSAREONE
            posc = 1;
#ifdef PHYSHYPWEIGHTAREONEANDNUEISZERO
            if(util::IsOscNuE(m) || util::IsOscNuEBar(m) || util::IsOscNuMuBar(m) || util::IsOscNuMu(m) ) posc = 0;
#endif
#endif
            // number of events in (Ereco,Etrue) bin (oscillations + systematics applied)
            double nev_osc = nev * posc;
            
            array[ir-1] = nev_osc;
            ////
          }//ir
        }//im
     }//it
   }//is
  }
  //
  // 3+n flavour sterile oscillations
  //
  if(hyp == kHypOscSterile4F || hyp == kHypOscSterile5F) 
  {
    
    NuOscSterile osc;
    
    // set oscillation parameters to the oscillation probability calculation code
    double mixing12 = fPhysHypothesis.CurrValue(kPhypMixing12);
    double mixing13 = fPhysHypothesis.CurrValue(kPhypMixing13);
    double mixing14 = fPhysHypothesis.CurrValue(kPhypMixing14);
    double mixing23 = fPhysHypothesis.CurrValue(kPhypMixing23);
    double mixing24 = fPhysHypothesis.CurrValue(kPhypMixing24);
    double mixing34 = fPhysHypothesis.CurrValue(kPhypMixing34);
    double dmsq     = fPhysHypothesis.CurrValue(kPhypDm12sq);
    double Dmsq     = fPhysHypothesis.CurrValue(kPhypDm23sq);
    double dm41sq   = fPhysHypothesis.CurrValue(kPhypDm41sq);
    double delta_12 = fPhysHypothesis.CurrValue(kPhypDelta12);
    double delta_13 = fPhysHypothesis.CurrValue(kPhypDelta13);
    double delta_34 = fPhysHypothesis.CurrValue(kPhypDelta34);
    
    double dm51sq=dm41sq,mixing15=0,mixing25=0,mixing35=0,delta_15=0,delta_35=0;
    if (hyp == kHypOscSterile5F) {
      mixing15 = fPhysHypothesis.CurrValue(kPhypMixing15);
      mixing25 = fPhysHypothesis.CurrValue(kPhypMixing25);
      mixing35 = fPhysHypothesis.CurrValue(kPhypMixing35);
      dm51sq   = fPhysHypothesis.CurrValue(kPhypDm51sq);
      delta_15 = fPhysHypothesis.CurrValue(kPhypDelta15);
      delta_35 = fPhysHypothesis.CurrValue(kPhypDelta35);
    }
      
#ifdef __show_hypothesis__
    cout << "**  Exp/SKMC >> " << fPhysHypothesis;
#endif
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    
    if (omp == kNuMixSinsqTheta) {
      osc.SetMixingAngles(mixing12,mixing13,mixing14,mixing15,mixing23,mixing24,mixing25,mixing34,mixing35);
    }
    else if (omp == kNuMixSinsq2Theta) {
      if (Opt::Instance()->gOscOctantChoice == conv::kNuOctantOne) {
	osc.SetMixingDoubleAnglesOct1(mixing12,mixing13,mixing14,mixing15,mixing23,mixing24,mixing25,mixing34,mixing35);
      }
      else if (Opt::Instance()->gOscOctantChoice == conv::kNuOctantTwo) {
	osc.SetMixingDoubleAnglesOct2(mixing12,mixing13,mixing14,mixing15,mixing23,mixing24,mixing25,mixing34,mixing35);
      }
      else {
	  cerr << "**  Exp/SKMC >> Can't handle choice of "
	   << "\n octant: " << util::AsString(Opt::Instance()->gOscOctantChoice)
	   << endl;
	exit(1);
      }
    }
    else {
      cerr << "**  Exp/SKMC >> Can't handle choice of "
	   << "\n oscillation mixing params: " << util::AsString(omp)
	   << endl;
      exit(1);
    }
    
    // May need to put some decision on mass splitting convention here.
    //~ NuOscDmsqScheme_t  dmsqscheme = Opt::Instance()->gOscDmsqScheme;
    
    osc.SetDmsq(dmsq, Dmsq, dm41sq, dm51sq);
    osc.SetDeltaCP(delta_12, delta_13, delta_15, delta_34, delta_35);

    // loop over samples
    for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;

     int n_ereco_bins = kNErecoBins[is];
     int n_etrue_bins = kNEtrueBins[is];

     //assert(this -> Template(s,kSKmNuMuCCQE,kSKMCtSys));
     //TAxis * ereco_axis = this->Template(s,kSKmNuMuCCQE,kSKMCtSys)->GetXaxis(); 
     //TAxis * etrue_axis = this->Template(s,kSKmNuMuCCQE,kSKMCtSys)->GetYaxis(); 

     // loop over Etrue
     for(int it = 1; it <= n_etrue_bins; it++) {

        // calculate oscillation probabilities (incl. matter effects) for the current Etrue value
       double Etrue = (util::EtrueBinEdges(s,it-1) + util::EtrueBinEdges(s,it) ) / 2;

	// skip unnecessary calculations
	bool bin_is_not_empty = this->IsThereAnEvent(s, Etrue);
	if(!bin_is_not_empty) continue;
	
	double p_osc    [NuOscSterile::nFlav][NuOscSterile::nFlav] = {};
	double p_osc_bar[NuOscSterile::nFlav][NuOscSterile::nFlav] = {};
	
	osc.CalcMatterOscProbAll(false, Etrue, &p_osc    [0][0]);
	osc.CalcMatterOscProbAll(true , Etrue, &p_osc_bar[0][0]);
	
	double p_numu_numu       = p_osc    [1][1];
	double p_numubar_numubar = p_osc_bar[1][1];
	double p_nue_nue         = p_osc    [0][0];
	double p_nuebar_nuebar   = p_osc_bar[0][0];
	double p_numu_nue        = p_osc    [1][0];
	double p_numubar_nuebar  = p_osc_bar[1][0];
	double p_nue_numu        = p_osc    [0][1];
	double p_nuebar_numubar  = p_osc_bar[0][1];
	
	double p_numu_nc    = p_osc    [1][0]+p_osc    [1][1]+p_osc    [1][2];
	double p_numubar_nc = p_osc_bar[1][0]+p_osc_bar[1][1]+p_osc_bar[1][2];
	double p_nue_nc     = p_osc    [0][0]+p_osc    [0][1]+p_osc    [0][2];
	double p_nuebar_nc  = p_osc_bar[0][0]+p_osc_bar[0][1]+p_osc_bar[0][2];
	
#ifdef __show_hypothesis_probabilities__
	cout << "Exp/SKMC >>  E = " << Etrue << " GeV" << endl;
	cout << "Exp/SKMC >>  * P (numu     -> numu   ) = " << p_numu_numu       << endl;
	cout << "Exp/SKMC >>  * P (numubar  -> numubar) = " << p_numubar_numubar << endl;
	cout << "Exp/SKMC >>  * P (nue      -> nue    ) = " << p_nue_nue         << endl;
	cout << "Exp/SKMC >>  * P (nuebar   -> nuebar ) = " << p_nuebar_nuebar   << endl;
	cout << "Exp/SKMC >>  * P (numu     -> nue    ) = " << p_numu_nue        << endl;
	cout << "Exp/SKMC >>  * P (numubar  -> nuebar ) = " << p_numubar_nuebar  << endl;
	cout << "Exp/SKMC >>  * P (nue      -> numu   ) = " << p_nue_numu        << endl;
	cout << "Exp/SKMC >>  * P (nuebar   -> numubar) = " << p_nuebar_numubar  << endl;
	
	cout << "Exp/SKMC >>  * P (numu     -> active ) = " << p_numu_nc         << endl;
	cout << "Exp/SKMC >>  * P (numubar  -> active ) = " << p_numubar_nc      << endl;
	cout << "Exp/SKMC >>  * P (nue      -> active ) = " << p_nue_nc          << endl;
	cout << "Exp/SKMC >>  * P (nuebar   -> active ) = " << p_nuebar_nc       << endl;
	
#endif
        
	// loop over modes
	for(int im=0; im < knm; im++) {
	  SuperkMode_t m = kSKReactionModes[is][im];
	  if(m == kSKmUndefined) continue;
	  // osc. probability
	  double posc = 1.;
	  if(kIsCC[is][im])
	  {
	    if     ( util::IsNuMu(m)      ) { posc = p_numu_numu;       }
	    else if( util::IsNuMuBar(m)   ) { posc = p_numubar_numubar; }
	    else if( util::IsNuE(m)       ) { posc = p_nue_nue;         } 
	    else if( util::IsNuEBar(m)    ) { posc = p_nuebar_nuebar;   }
	    else if( util::IsOscNuMu(m)   ) { posc = p_nue_numu;        }
	    else if( util::IsOscNuMuBar(m)) { posc = p_nuebar_numubar;  }
	    else if( util::IsOscNuE(m)    ) { posc = p_numu_nue;        }
	    else if( util::IsOscNuEBar(m) ) { posc = p_numubar_nuebar;  }
	  }
	  else { // NC oscillations to sterile
	    if (util::IsNuMu(m)) {
	      posc = p_numu_nc;
	    }
	    else if (util::IsNuMuBar(m)) {
	      posc = p_numubar_nc;
	    }
	    else if (util::IsNuE(m)) {
	      posc = p_nue_nc;
	    }
	    else if (util::IsNuEBar(m)) {
	      posc = p_nuebar_nc;
	    }
	  }
	      
	  // loop over Ereco
	  for(int ir = 1; ir <= n_ereco_bins ; ir++) {
	    // number of events in (Ereco,Etrue) bin (systematics applied)
	    double nev = 0;
	    
	    nev = fMCTemplArrayFin[is].GetEntryRef(im, it-1, ir-1); 
	    ////
	    
#ifdef PHYSHYPWEIGHTSAREONE
	    posc = 1;
#ifdef PHYSHYPWEIGHTAREONEANDNUEISZERO
	    if(util::IsOscNuE(m)) posc = 0;
#endif
#endif
	    // number of events in (Ereco,Etrue) bin (oscillations + systematics applied)
	    double nev_osc = nev * posc;
	    
	    fMCTemplArrayFin[is].GetEntryRef(im, it-1, ir-1) = nev_osc;
	    ////
	  }//ir
	}//im
     }//it
   }//is
  } // end of angle/splitting/phase sterile hypotheses
/*
  //
  // 2-flavour oscillations / numu disappearance
  //
  else
  if(hyp == kHypNuMuDisapp2F) {
    double mixing23 = fPhysHypothesis.CurrValue(kPhypMixing23);
    double dm23sq   = fPhysHypothesis.CurrValue(kPhypDm23sq);
    // loop over samples
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];

      if(s!=kSKsFCFV1RingMuLike || s!=kSKsRHCFCFV1RingMuLike) continue; // MODIFIED ANTINU

      int n_ereco_bins = kNErecoBins[is];
      int n_etrue_bins = kNEtrueBins[is];

      // loop over sample modes
      for(int im = 0; im < knm; im++) { 
	 SuperkMode_t m = kSKReactionModes[is][im];
	 if(m == kSKmUndefined) continue;
         // consider template?
         bool skip = this->SkipTemplate(m);
         if(skip) continue;
         
         // loop over Etrue
         for(int it = 1; it <= n_etrue_bins; it++) {
           double Etrue = (util::EtrueBinEdges(s,it-1) + util::EtrueBinEdges(s,it)) / 2;
           double psurv = util::NuMuSurvProbOsc2Flavour(mixing23, dm23sq, Etrue, m, fUnphysicalOscProb);
           // loop over Ereco
	   for(int ir = 1; ir <= n_ereco_bins; ir++) {
	     
	     double nev = 0.0;
	     if     (s==kSKsFCFV1RingMuLike)    nev =  fMCTemplArrayFin1Rmu   [im][it-1][ir-1]; // MODIFIED
	     else if(s==kSKsRHCFCFV1RingMuLike) nev =  fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]; // ANTINU
	     
	     nev *= psurv;
	   
	     if     (s==kSKsFCFV1RingMuLike)    fMCTemplArrayFin1Rmu   [im][it-1][ir-1] = nev; // MODIFIED
	     else if(s==kSKsRHCFCFV1RingMuLike) fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1] = nev; // ANTINU
	     
	   }//ir
         }//it
      }//im
    }//is
  }

  //
  // decay
  //
  else
  if(hyp == kHypDecay) {
    double sinsq2theta = fPhysHypothesis.CurrValue(kPhypSinsq2ThetaDecay);
    double alpha       = fPhysHypothesis.CurrValue(kPhypAlphaDecay);
    // loop over samples
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];
     
      if(s!=kSKsFCFV1RingMuLike || s!=kSKsRHCFCFV1RingMuLike) continue;  // MODIFIED ANTINU
      
      int n_ereco_bins = kNErecoBins[is];
      int n_etrue_bins = kNEtrueBins[is];

      // loop over sample modes
      for(int im = 0; im < knm; im++) { 
         SuperkMode_t m = kSKReactionModes[is][im];
         if(m == kSKmUndefined) continue;
         // consider template?
         bool skip = this->SkipTemplate(m);
         if(skip) continue;
         // get source (systematics) and target (systematics + phys.hypothesis) template
         
         // loop over Etrue
         for(int it = 1; it <= n_etrue_bins; it++) {
           double Etrue = (util::EtrueBinEdges(s,it-1) + util::EtrueBinEdges(s,it)) / 2;
           double psurv = util::NuMuSurvProbDecay(sinsq2theta, alpha, Etrue, m);
           // loop over Ereco
           for(int ir = 1; ir <= n_ereco_bins; ir++) {
              
              double nev = 0.0;
              if     (s==kSKsFCFV1RingMuLike)    nev =  fMCTemplArrayFin1Rmu   [im][it-1][ir-1]; // MODIFIED
              else if(s==kSKsRHCFCFV1RingMuLike) nev =  fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]; // ANTINU
              
              nev *= psurv;
            
              if     (s==kSKsFCFV1RingMuLike)    fMCTemplArrayFin1Rmu   [im][it-1][ir-1] = nev; // MODIFIED
              else if(s==kSKsRHCFCFV1RingMuLike) fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1] = nev; // ANTINU
              
           }//ir
         }//it
      }//im
    }//is
  }


  //
  // decoherence
  //
  else
  if(hyp == kHypDecoh) {
    double sinsq2theta = fPhysHypothesis.CurrValue(kPhypSinsq2ThetaDecoh);
    double musq        = fPhysHypothesis.CurrValue(kPhypMusqDecoh);
    // loop over samples
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      
      if(s!=kSKsFCFV1RingMuLike || s!=kSKsRHCFCFV1RingMuLike) continue;  // MODIFIED ANTINU
      
      int n_ereco_bins = kNErecoBins[is];
      int n_etrue_bins = kNEtrueBins[is];

      // loop over sample modes
      for(int im = 0; im < knm; im++) { 
         SuperkMode_t m = kSKReactionModes[is][im];
         if(m == kSKmUndefined) continue;
         // consider template?
         bool skip = this->SkipTemplate(m);
         if(skip) continue;
         // get source (systematics) and target (systematics + phys.hypothesis) template
         TH2D * mctempl_tgt = this -> Template(s,m,kSKMCtFin);
         if(!mctempl_tgt) continue; 
         // loop over Etrue
         for(int it = 1; it <= n_etrue_bins; it++) {
           double Etrue = (util::EtrueBinEdges(s,it-1) + util::EtrueBinEdges(s,it)) / 2;
           double psurv = util::NuMuSurvProbDecoh(sinsq2theta, musq, Etrue, m);
           // loop over Ereco
           for(int ir = 1; ir <= n_ereco_bins; ir++) {
             
             double nev = 0.0;
             if     (s==kSKsFCFV1RingMuLike)    nev =  fMCTemplArrayFin1Rmu   [im][it-1][ir-1]; // MODIFIED
             else if(s==kSKsRHCFCFV1RingMuLike) nev =  fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]; // ANTINU
             
             nev *= psurv;
           
             if     (s==kSKsFCFV1RingMuLike)    fMCTemplArrayFin1Rmu   [im][it-1][ir-1] = nev; // MODIFIED
             else if(s==kSKsRHCFCFV1RingMuLike) fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1] = nev; // ANTINU
             
           }//ir
         }//it
      }//im
    }//is
  }
  * */
#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtFin, "mcnorm_fluxtuned_banfftuned_syst_phys");
#endif
}
//________________________________________________________________________________________
bool SuperkMC::UpdateErecoSpectra(void)
{
  bool applied_negevents_condition = false;
  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;

    int n_ereco_bins = kNErecoBins[is];
    int n_etrue_bins = kNEtrueBins[is];
    
    // Recast VALORArray blocks so they can be used addressed like normal arrays (fast).
    // Hard to read, so should only be used in code for which speed is essential.
    typedef double (templmtype)[knm][n_etrue_bins][n_ereco_bins];
    typedef double (ermtype)   [knm][n_ereco_bins];
    const templmtype & templFin = *reinterpret_cast< const templmtype * > (fMCTemplArrayFin[is].GetArray());
    ermtype & erecoMode         = *reinterpret_cast< ermtype          * > (fMCErecoSpectrumArray[is].GetArray());
    
    //~ const double(&templFin)[knm][n_etrue_bins][n_ereco_bins] = *reinterpret_cast< const double (*) [knm][n_etrue_bins][n_ereco_bins]> (fMCTemplArrayFin[is].GetArray());
    //~ double(&erecoMode)[knm][n_ereco_bins] = *reinterpret_cast<double (*) [knm][n_ereco_bins]> (fMCErecoSpectrumArray[is].GetArray());
    
    double * erecoTot = fMCErecoSpectrumTotArray[is].GetArray();

    // loop over Ereco
    for(int ir = 1; ir <= n_ereco_bins; ir++) {
      double nev_tot = 0;
      // loop over modes
      for(int im = 0; im < knm; im++) { 
        SuperkMode_t m = kSKReactionModes[is][im];
        if(m == kSKmUndefined) continue;
        // consider template?
        bool skip = this->SkipTemplate(m);
        if(skip) continue;
        
        // loop over Etrue
        double nev = 0;
        for(int it = 1; it <= n_etrue_bins; it++) {
          nev += templFin[im][it-1][ir-1];
        }//it
        nev_tot += nev;
#ifdef __apply_minevents_to_individual_modes__
        if(nev < 0) {
#ifdef __show_negative_events_mesg_debug__
        cout << "Exp/SKMC/UpdateErecoSpectra >> set negative number of events " << nev
             << " to zero " << 0 << "\t" << " Ereco: " << reco_spectrum_tot->GetXaxis()->GetBinCenter(ir) << endl;
#endif
          nev = 0;
        } // if nev < 0
#endif //__apply_minevents_to_individual_modes__

        erecoMode[im][ir-1] = nev;   
      }//im
      
      erecoTot[ir-1] = nev_tot;
      
    }//ir
  }//is
  
  applied_negevents_condition = this->ApplyEnergyScaleSystematicToTotalSpectrum();
  
  return applied_negevents_condition;
}
//________________________________________________________________________________________
void SuperkMC::FillEtrueSpectra(void)
{
  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    TH1D * true_spectrum_tot = fMCEtrueSpectrumTot[is];

    int n_ereco_bins = kNErecoBins[is];

    // loop over Etrue
    for(int it = 1; it <= true_spectrum_tot->GetNbinsX(); it++) {
       double nev_tot = 0;

       // loop over modes
       for(int im = 0; im < knm; im++) { 
          SuperkMode_t m = kSKReactionModes[is][im];
          if(m == kSKmUndefined) continue;
          // consider template?
          bool skip = this->SkipTemplate(m);
          if(skip) continue;

          TH2D * mctempl = this -> Template(s,m,kSKMCtOri);
          
          // loop over Ereco
          double nev = 0;
          for(int ir = 1; ir <= n_ereco_bins; ir++) {
            nev += mctempl->GetBinContent(ir,it);
          }//ir
          nev_tot += nev;
       }//im
      
       fMCEtrueSpectrumTotArray[is].GetEntryRef(it-1) = nev_tot;

    }//it
  }//is
}
//________________________________________________________________________________________
bool SuperkMC::IsThereAnEvent(SuperkSample_t s, double Etrue) const
{
  int is = arr_index->SuperkSampleIdx(s);
  assert(is != -1);

  TH1D * true_spectrum_tot = fMCEtrueSpectrumTot[is];
  assert(true_spectrum_tot);

  int it = true_spectrum_tot->GetXaxis()->FindBin(Etrue);
  bool ev = false;

  ev = fMCEtrueSpectrumTotArray[is].GetEntryRef(it-1) > kASmallNum;
  ////
  
#ifdef __show_isanevent_mesg_debug__
  if(ev){
    
    cout << "Exp/SuperkMC >> There are " << fMCEtrueSpectrumTotArray[is].GetEntryRef(it-1);
    cout << " events for sample = " << util::AsString(s) 
         << ", Etrue = " << Etrue
         << endl;
  }
#endif
  
  return ev;
} 
//________________________________________________________________________________________
double SuperkMC::AnalyzedSamplePOT(SuperkSample_t s) const
{
  if(util::IsRHCSample(s))
    return fPOTRHC;
  else if(util::IsFHCSample(s)){
    /// NEW MULTIRING
    if(util::IsMultiRingSample(s))
      return fPOTMultiRing;
    ///
    else
      return fPOT;  
  }
  else {
    cerr << "Sample not selected" << endl;
    exit(1);
  }
}
//________________________________________________________________________________________
TTree * SuperkMC::MCTree(SuperkSample_t s, SuperkMode_t m)
{
  if (s!=kSKsUndefined)
  {
    if ( util::IsNuMu(m)      ) return fMCTreeNuMu     [s];
    if ( util::IsNuMuBar(m)   ) return fMCTreeNuMuBar  [s];
    if ( util::IsNuE(m)       ) return fMCTreeNuE      [s];
    if ( util::IsNuEBar(m)    ) return fMCTreeNuEBar   [s];
    if ( util::IsOscNuE(m)    ) return fMCTreeOscNuE   [s];
    if ( util::IsOscNuEBar(m) ) return fMCTreeOscNuEBar[s];
  }

  return 0;
}
//________________________________________________________________________________________
double SuperkMC::MCPOT(SuperkSample_t s, SuperkMode_t m)
{
  if(util::IsFHCSample(s))
  {
    /// NEW MULTIRING 
    if(util::IsMultiRingSample(s)){
     if ( util::IsNuMu(m)      ) return fMCPOTMultiRingNuMu;
     if ( util::IsNuMuBar(m)   ) return fMCPOTMultiRingNuMuBar;
     if ( util::IsNuE(m)       ) return fMCPOTMultiRingNuE;
     if ( util::IsNuEBar(m)    ) return fMCPOTMultiRingNuEBar;
     if ( util::IsOscNuE(m)    ) return fMCPOTMultiRingOscNuE;
     if ( util::IsOscNuEBar(m) ) return fMCPOTMultiRingOscNuEBar;
    }
    ///
    else{
      if ( util::IsNuMu(m)      ) return fMCPOTNuMu;
      if ( util::IsNuMuBar(m)   ) return fMCPOTNuMuBar;
      if ( util::IsNuE(m)       ) return fMCPOTNuE;
      if ( util::IsNuEBar(m)    ) return fMCPOTNuEBar;
      if ( util::IsOscNuE(m)    ) return fMCPOTOscNuE;
      if ( util::IsOscNuEBar(m) ) return fMCPOTOscNuEBar;
    }
  }  
  else if(util::IsRHCSample(s))
  {
     if ( util::IsNuMu(m)      ) return fMCPOTRHCNuMu;
     if ( util::IsNuMuBar(m)   ) return fMCPOTRHCNuMuBar;
     if ( util::IsNuE(m)       ) return fMCPOTRHCNuE;
     if ( util::IsNuEBar(m)    ) return fMCPOTRHCNuEBar;
     if ( util::IsOscNuE(m)    ) return fMCPOTRHCOscNuE;
     if ( util::IsOscNuEBar(m) ) return fMCPOTRHCOscNuEBar;
  }
  
  ////
  
  return 0;
}
//________________________________________________________________________________________
PhysHypothesis & SuperkMC::Hypothesis(void) 
{ 
  return fPhysHypothesis; 
}
//________________________________________________________________________________________
void SuperkMC::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  if(opt.Contains("EtrueSpectra/NuMu")) { 
      fMCTreeNuMu[0]->Draw("Enu");
      fMCTreeNuMu[0]->Draw("Enu","Mode==1","same");
      fMCTreeNuMu[0]->Draw("Enu","Mode==11||Mode==12||Mode==13","same");
      c->Update();
  }
/*
  if(opt.Contains("1RingMuLike/ErecoEtrueTemplates")) { 
    double emin=0;
    double emax=6;
    const char * draw_opt = "COLZ";
    c->Divide(4,4);
    c->cd(1);
    c->GetPad(1)->SetLogz();
    c->GetPad(1)->SetBottomMargin(0.32);
    c->GetPad(1)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->Draw(draw_opt);
    c->cd(2);
    c->GetPad(2)->SetLogz();
    c->GetPad(2)->SetBottomMargin(0.32);
    c->GetPad(2)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->Draw(draw_opt);
    c->cd(3);
    c->GetPad(3)->SetLogz();
    c->GetPad(3)->SetBottomMargin(0.32);
    c->GetPad(3)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->Draw(draw_opt);
    c->cd(4);
    c->GetPad(4)->SetLogz();
    c->GetPad(4)->SetBottomMargin(0.32);
    c->GetPad(4)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->Draw(draw_opt);
    c->cd(5);
    c->GetPad(5)->SetLogz();
    c->GetPad(5)->SetBottomMargin(0.32);
    c->GetPad(5)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->Draw(draw_opt);
    c->cd(6);
    c->GetPad(6)->SetLogz();
    c->GetPad(6)->SetBottomMargin(0.32);
    c->GetPad(6)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->Draw(draw_opt);
    c->cd(7);
    c->GetPad(7)->SetLogz();
    c->GetPad(7)->SetBottomMargin(0.32);
    c->GetPad(7)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->Draw(draw_opt);
    c->cd(8);
    c->GetPad(8)->SetLogz();
    c->GetPad(8)->SetBottomMargin(0.32);
    c->GetPad(8)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->Draw(draw_opt);
    c->cd(9);
    c->GetPad(9)->SetLogz();
    c->GetPad(9)->SetBottomMargin(0.32);
    c->GetPad(9)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->Draw(draw_opt);
    c->cd(10);
    c->GetPad(10)->SetLogz();
    c->GetPad(10)->SetBottomMargin(0.32);
    c->GetPad(10)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->Draw(draw_opt);
    c->cd(11);
    c->GetPad(11)->SetLogz();
    c->GetPad(11)->SetBottomMargin(0.32);
    c->GetPad(11)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->Draw(draw_opt);
    c->Update();
    c->cd(12);
    c->GetPad(12)->SetLogz();
    c->GetPad(12)->SetBottomMargin(0.32);
    c->GetPad(12)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->Draw(draw_opt);
    c->cd(13);
    c->GetPad(13)->SetLogz();
    c->GetPad(13)->SetBottomMargin(0.32);
    c->GetPad(13)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->Draw(draw_opt);
    c->cd(14);
    c->GetPad(14)->SetLogz();
    c->GetPad(14)->SetBottomMargin(0.32);
    c->GetPad(14)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->Draw(draw_opt);
    c->cd(15);
    c->GetPad(15)->SetLogz();
    c->GetPad(15)->SetBottomMargin(0.32);
    c->GetPad(15)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth]->Draw(draw_opt);
    c->Update();
    c->cd(16);
    c->GetPad(16)->SetLogz();
    c->GetPad(16)->SetBottomMargin(0.32);
    c->GetPad(16)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->Draw(draw_opt);
    c->Update();
  }
*/
  c->Update();
}
//________________________________________________________________________________________
bool SuperkMC::SkipTemplate(SuperkMode_t m) 
{
// ignore nue_signal templates for all but the 3-flavour oscillation hypothesis

  bool mode_is_nuesig = util::IsOscNuE(m);
  bool mode_is_nuebarsig = util::IsOscNuEBar(m);


  PhysHypothesis_t hyp = fPhysHypothesis.Hypothesis();
  bool hyp_is_not_3flavour_osc = (hyp != kHypOsc3F);

  bool skip = ((mode_is_nuesig || mode_is_nuebarsig) && hyp_is_not_3flavour_osc );

#ifdef __show_skipped_templates__
  if(skip) {
    cout << "Exp/SKMC >> Skipping template: " << util::AsString(m) << endl;
  }
#endif

  return skip;
}
//________________________________________________________________________________________
void SuperkMC::NevDebug(SKMCTemplateType_t t, string name)
{
  
  for(int is = 0; is < kns; is++) {
      
    SuperkSample_t s = kSKSamples[is];
    
    if(s == kSKsUndefined) continue;
      
    TFile* oute = new TFile(Form("Template_%s%s.root", name.c_str(),util::AsStringCompact(s).c_str()),"RECREATE");
    static const int erecobins = util::NumErecoBins(s);//70;
    static const double * erecobinedges = util::ErecoBinEdges(s);
    
    //if(s == kSKsFCFV1RingELike) continue;
   // if(s == kSKsFCFV1RingMuLike) continue;
   // if(s == kSKsRHCFCFV1RingMuLike) continue;    // ANTINU
   // if(s == kSKsRHCFCFV1RingELike) continue;     // ANTINU

    double nsk = 0;
    for(int im = 0; im < knm; im++) {
      SuperkMode_t m = kSKReactionModes[is][im];
      if(m == kSKmUndefined) continue;
      TH2D h = this -> Get2DHist(s,m,t);
      TH2D * mctempl = &h;

      //first write the 2d template, and cout the mode integral
      mctempl->Write();
      //      cout << name << " " << util::AsString(s) << " " << util::AsString(m) << " " << mctempl->Integral() << endl;
      nsk += mctempl->Integral();

      //fill a 1D ereco spectrum
      TH1D * hrec = new TH1D(Form("erec_int_%s",AsStringCompact_KDI(m).c_str()), Form("erec_int_%s",AsStringCompact_KDI(m).c_str()), erecobins, erecobinedges);
      // loop over Ereco
      for(int ir = 1; ir <= hrec->GetNbinsX(); ir++) {
        //find low and high of ereco for current bin
        double rmin = hrec->GetXaxis()->GetBinLowEdge(ir) + kOneMeV;
        double rmax = hrec->GetXaxis()->GetBinUpEdge(ir)  - kOneMeV;
        //integrate the template
        int ir_min = mctempl->GetXaxis()->FindBin(rmin);
        int ir_max = mctempl->GetXaxis()->FindBin(rmax);
        int it_bins = mctempl->GetNbinsY();
        double nev = mctempl->Integral(ir_min, ir_max, 1, it_bins);
        hrec->SetBinContent(ir,nev);
      }
      hrec->Write();
      delete hrec;

      //fill a 1D ereco spectrum (fine binning)
      TH1D * hrecf = new TH1D(Form("erec_int_%s_fine",AsStringCompact_KDI(m).c_str()), Form("erec_int_%s",AsStringCompact_KDI(m).c_str()), mctempl->GetNbinsX(), mctempl->GetXaxis()->GetXbins()->GetArray());
      // loop over Ereco
      for(int ir = 1; ir <= hrecf->GetNbinsX(); ir++) {
        //find low and high of ereco for current bin
        double rmin = hrecf->GetXaxis()->GetBinLowEdge(ir) + kOneMeV;
        double rmax = hrecf->GetXaxis()->GetBinUpEdge(ir)  - kOneMeV;
        //integrate the template
        int ir_min = mctempl->GetXaxis()->FindBin(rmin);
        int ir_max = mctempl->GetXaxis()->FindBin(rmax);
        int it_bins = mctempl->GetNbinsY();
        double nev = mctempl->Integral(ir_min, ir_max, 1, it_bins);
        hrecf->SetBinContent(ir,nev);
      }
      hrecf->Write();
      delete hrecf;

      //fill a 1D etrue spectrum
      TH1D * htru = new TH1D(Form("enu_int_%s",AsStringCompact_KDI(m).c_str()), Form("enu_int_%s",AsStringCompact_KDI(m).c_str()), util::NumVariableMCEtrueBins(s,m), util::VariableMCEtrueBinEdges(s,m));
      // loop over Etrue
      for(int it = 1; it <= htru->GetNbinsX(); it++) {
        //find low and high of etrue for current bin
        double tmin = htru->GetXaxis()->GetBinLowEdge(it) + kOneMeV;
        double tmax = htru->GetXaxis()->GetBinUpEdge(it)  - kOneMeV;
        //integrate the template
        int it_min = mctempl->GetYaxis()->FindBin(tmin);
        int it_max = mctempl->GetYaxis()->FindBin(tmax);
        int ir_bins = mctempl->GetNbinsX();
        double nev = mctempl->Integral(1, ir_bins, it_min, it_max);
        htru->SetBinContent(it,nev);
      }
      htru->Write();
      delete htru;
    }//mode
    this->UpdateErecoSpectra();
    //get the total ereco spectrum
    TH1D * hfull = this->CurrErecoSpectrum(s);
    hfull->Write("ereco_spectrum");
    delete hfull;
    //    cout << "Total: " << nsk << endl << endl;
    oute->Close();
    delete oute;
  }//sample
 
}
//________________________________________________________________________________________
void SuperkMC::MCStatErrDebug()
{

  const int ntoys  = 10000;

  TFile out("Template_mcstaterr.root","RECREATE");
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];

    if(s == kSKsUndefined) continue;
    if(s == kSKsFCFV1RingELike) continue;

    if(s == kSKsRHCFCFV1RingMuLike) continue;   // RHC
    if(s == kSKsRHCFCFV1RingELike) continue;    // RHC

    /// NEW MULTIRING
    if(s == kSKsFCFVMultiRingMuLike) continue; // MultiRing FHC
    if(s == kSKsFCFVMultiRingELike) continue;  // MultiRing FHC
    ///

    int n_ereco_bins = kNErecoBins[is];

    //fill the nominal spectrum
    TH1D * spectrum           = new TH1D("spe",          "spe",          n_ereco_bins, util::ErecoBinEdges(kSKsFCFV1RingMuLike));
    for(int im = 0; im < knm; im++) {
      SuperkMode_t m = kSKReactionModes[is][im];
      if(m == kSKmUndefined) continue;
      TH2D * mctempl = this -> Template(s,m,kSKMCtOri);
      //and the 1d spectrum
      for(int ix = 1; ix <= mctempl->GetNbinsX(); ix++) {
	for(int iy = 1; iy <= mctempl->GetNbinsY(); iy++) {
	  double nev_new = mctempl->GetBinContent(ix, iy);
	  double nev_old = spectrum->GetBinContent(ix);
	  spectrum->SetBinContent(ix, nev_new + nev_old);
	}//iy
      }//ix
    }//im

    TH1D * spectrum_with2derr = new TH1D("spe_with2derr","spe_with2derr",n_ereco_bins, util::ErecoBinEdges(kSKsFCFV1RingMuLike));
    TH1D * spectrum_fracerr   = new TH1D("spe_fracerr",  "spe_fracerr",  n_ereco_bins, util::ErecoBinEdges(kSKsFCFV1RingMuLike));
    TH1D * spectrum_err       = new TH1D("spe_err",      "spe_err",      n_ereco_bins, util::ErecoBinEdges(kSKsFCFV1RingMuLike));
    TH1D * spectrum_rms       = new TH1D("spe_rms",      "spe_rms",      n_ereco_bins, util::ErecoBinEdges(kSKsFCFV1RingMuLike));
    for(int im = 0; im < knm; im++) {
      SuperkMode_t m = kSKReactionModes[is][im];
      if(m == kSKmUndefined) continue;
      TH2D * mctempl = this -> Template(s,m,kSKMCtOri);

      //first write the 2d template, and cout the mode integral
      mctempl->Write();

      //also write the 2d MC stat frac error template
      TH2D * mctemplfracerr = new TH2D(*mctempl);
      TH2D * mctemplerr     = new TH2D(*mctempl);
      for(int ix = 1; ix <= mctempl->GetNbinsX(); ix++) {
	for(int iy = 1; iy <= mctempl->GetNbinsY(); iy++) {
	  double nev = mctempl->GetBinContent(ix, iy);
	  nev = TMath::Max(0.,nev);
	  double dnev = TMath::Sqrt(nev);
	  double fracerr = (nev > 0) ? 1/dnev : 0;
	  mctemplfracerr->SetBinContent(ix, iy, fracerr);
	  mctemplerr    ->SetBinContent(ix, iy, (1 + fracerr) * nev);
	}//iy
      }//ix
      mctemplfracerr->Write(Form("%s_fracerr",mctempl->GetName()));
      mctemplerr    ->Write(Form("%s_err",    mctempl->GetName()));


      //and the 1d spectrum
      for(int ix = 1; ix <= mctempl->GetNbinsX(); ix++) {
	for(int iy = 1; iy <= mctempl->GetNbinsY(); iy++) {
	  double nev_new = mctemplerr->GetBinContent(ix, iy);
	  if(mctempl->GetBinContent(ix,iy) > kASmallNum)
	    if(mctempl->GetBinContent(ix,iy) > mctemplerr->GetBinContent(ix,iy))
	      cout << mctempl->GetBinContent(ix,iy) << "\t" << mctemplerr->GetBinContent(ix,iy) << endl;
	  double nev_old = spectrum_with2derr->GetBinContent(ix);
	  spectrum_with2derr->SetBinContent(ix, nev_new + nev_old);
	  for(int ie = 1; ie < ntoys; ie++) {
	    int nev_fluct = Rnd::Instance()->RndToyMCstat().Poisson((int)nev_new);
	    double current = spectrum_rms->GetBinContent(ix);
	    spectrum_rms->SetBinContent(ix, current + TMath::Power(nev_fluct-nev_new,2));
	  }//ie
	}//iy
      }//ix
      cout << "Exp/SKMC >> Finished MC stat error debug loops for mode: " << util::AsStringCompact(m) << endl;
    }//im

    //go from sum squares to RMS
    for(int ix = 1; ix <= spectrum_rms->GetNbinsX(); ix++) {
      double sumsq = spectrum_rms->GetBinContent(ix);
      spectrum_rms->SetBinContent(ix, TMath::Sqrt(sumsq/ntoys) + spectrum->GetBinContent(ix));
    }

    //look at what the error would of been on the spectrum, if the 1d bins didn't exist
    for(int ix = 1; ix <= spectrum->GetNbinsX(); ix++) {
      double nev = spectrum->GetBinContent(ix);
      nev = TMath::Max(0.,nev);
      double dnev = TMath::Sqrt(nev);
      double fracerr = (nev > 0) ? 1/dnev : 0;
      spectrum_fracerr->SetBinContent(ix, fracerr);
      spectrum_err    ->SetBinContent(ix, (1 + fracerr) * nev);
    }
      

    spectrum->Write();
    //spectrum_with2derr->Write();
    spectrum_fracerr->Write();
    spectrum_err->Write();
    spectrum_rms->Write("spe_with2derr");
  }//is
  out.Close();
}
//____________________________________________________________________________________
TH2D SuperkMC::Get2DHist(conv::SuperkSample_t s, conv::SKMCTemplateType_t type) const
{
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  const int      n_ereco_bins    = kNErecoBins[is];
  const int      n_etrue_bins    = kNEtrueBins[is];
  const double * ereco_bin_edges = util::ErecoBinEdges(s);
  
  TH2D h(Form("nsk_%d",is), util::AsString(s).c_str(),
    n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
    
  h.SetDirectory(0);
  h.GetXaxis()->SetTitle("E_{reco} (GeV)");
  h.GetYaxis()->SetTitle("E_{true} (GeV)");
  
  for (int im = 0 ; im < kNSuperkModes ; im++ ) {
    SuperkMode_t m = kSKReactionModes[is][im];
    if (m==kSKmUndefined) continue;
    for (int it = 0 ; it < n_etrue_bins ; it++) {
      for (int ir = 0 ; ir < n_ereco_bins ; ir++) {
        int binno = h.GetBin(ir+1,it+1);
        if (type == conv::kSKMCtOri) {
          h.SetBinContent(binno, fMCTemplArrayOri[is].GetEntryRef(im,it,ir));
        }
        if (type == conv::kSKMCtNom) {
          h.SetBinContent(binno, fMCTemplArrayNom[is].GetEntryRef(im,it,ir));
        }
        if (type == conv::kSKMCtTun) {
          h.SetBinContent(binno, fMCTemplArrayTun[is].GetEntryRef(im,it,ir));
        }
        if (type == conv::kSKMCtFin) {
          h.SetBinContent(binno, fMCTemplArrayFin[is].GetEntryRef(im,it,ir));
        }
      }
    }
  }
  return h;
}
//____________________________________________________________________________________
TH2D SuperkMC::Get2DHist(conv::SuperkSample_t s, conv::SuperkMode_t m , conv::SKMCTemplateType_t type) const
{
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  const int      n_ereco_bins    = kNErecoBins[is];
  const int      n_etrue_bins    = kNEtrueBins[is];
  const double * ereco_bin_edges = util::ErecoBinEdges(s);
  
  TH2D h(Form("nsk_%d_%d",is,im), util::AsString(m).c_str(),
    n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
    
  h.SetDirectory(0);
  h.GetXaxis()->SetTitle("E_{reco} (GeV)");
  h.GetYaxis()->SetTitle("E_{true} (GeV)");
  
  for (int it = 0 ; it < n_etrue_bins ; it++) {
    for (int ir = 0 ; ir < n_ereco_bins ; ir++) {
      int binno = h.GetBin(ir+1,it+1);
      if (type == conv::kSKMCtOri) {
        h.SetBinContent(binno, fMCTemplArrayOri[is].GetEntryRef(im,it,ir));
      }
      if (type == conv::kSKMCtNom) {
        h.SetBinContent(binno, fMCTemplArrayNom[is].GetEntryRef(im,it,ir));
      }
      if (type == conv::kSKMCtTun) {
        h.SetBinContent(binno, fMCTemplArrayTun[is].GetEntryRef(im,it,ir));
      }
      if (type == conv::kSKMCtFin) {
        h.SetBinContent(binno, fMCTemplArrayFin[is].GetEntryRef(im,it,ir));
      }
    }
  }
  
  return h;
}
