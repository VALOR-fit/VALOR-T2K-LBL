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
#define __show_init_mesg_debug__
//#define __ereco_scale_debug__
//#define __show_init_mesg_debug2__
//#define __show_update_mesg_info__
//#define __show_hypothesis__ 
//#define __show_hypothesis_probabilities__ 
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

#define APPLY_ESCALE_AT_SPECTRUM_STAGE
//#define POTWEIGHTSAREONE
//#define FLUXWEIGHTSAREONE
//#define BANFFWEIGHTSAREONE
//#define SYSTWEIGHTSAREONE
//#define PHYSHYPWEIGHTSAREONE
//#define PHYSHYPWEIGHTAREONEANDNUEISZERO
#define INCLUDEEVENTSOVERTHIRTYGEV // include Ereco overflow events in the last bin?
//#define RENORMALISE_MEC // NEEDS CHECKING renormalise non-MEC templates so that MEC isn't an arbitary addition of events, but respects the total normalisation

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

#include "SuperkMC.h"
#include "SuperkDetector.h"
#include "SuperkDetectorFSI.h"
#include "FSI.h"
#include "NuXSec.h"
#include "NuFlux.h"
#include "Utilities.h"
#include "Opt.h"
#include "Rnd.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;
using namespace t2kana::util;

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

  fMCTemplMECMartini = 0;
  fMCTemplMECNieves  = 0;
  for(int it=0; it<util::NumEtrueBins(kSKsFCFV1RingMuLike); it++) {
    for(int ir=0; ir<util::NumErecoBins(kSKsFCFV1RingMuLike); ir++) {
      fMCTemplArrayMECMartini[it][ir] = 0;
      fMCTemplArrayMECNieves [it][ir] = 0;
    }
  }
  
  //// ANTINU
  for(int it=0; it<util::NumEtrueBins(kSKsRHCFCFV1RingMuLike); it++) {
    for(int ir=0; ir<util::NumErecoBins(kSKsRHCFCFV1RingMuLike); ir++) {
      fMCTemplArrayRHCMECMartini[it][ir] = 0;
      fMCTemplArrayRHCMECNieves [it][ir] = 0;
    }
  }
  ////

  for(int is=0; is<kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    int n_etrue_bins = kNEtrueBins[is];
    int n_ereco_bins = kNErecoBins[is];
    int n_xs_ereco_bins = kNXsErecoBins[is];

    fMCErecoSpectrumTot[is] = 0;
    fMCEtrueSpectrumTot[is] = 0;
    for(int im=0; im<knm; im++) {
      fMCTemplOri      [is][im] = 0;
      fMCTemplNom      [is][im] = 0;
#ifndef PRODUCTIONMODE
      fMCTemplTun      [is][im] = 0;
      fMCTemplSys      [is][im] = 0;
      fMCTemplFin      [is][im] = 0;
#endif
      fMCErecoSpectrum [is][im] = 0;
      fMCVariableTempl [is][im] = 0;
      
      fNNonEmptyVariableBins [is][im] = 0;
      for (int ivb = 0; ivb < kMaxVarBins; ivb++) {
	fNNonEmptyFixedBins    [is][im][ivb] = 0;
	fNonEmptyFixedBins     [is][im][ivb] = NULL;
	fNonEmptyVariableBins  [is][im][ivb] = 0;
      }
      
      //////////////
      for(int it=0; it<n_etrue_bins; it++) {
	for(int ir=0; ir<n_ereco_bins; ir++) {
	  
	  // MODIFIED
	  if (s==kSKsFCFV1RingMuLike){
	    fMCTemplArrayOri1Rmu[im][it][ir] = 0;
	    fMCTemplArrayNom1Rmu[im][it][ir] = 0;
	    fMCTemplArrayTun1Rmu[im][it][ir] = 0;
	    fMCTemplArrayFin1Rmu[im][it][ir] = 0;
	  }
	  else if(s==kSKsFCFV1RingELike) {
	    fMCTemplArrayOri1Re[im][it][ir] = 0;
	    fMCTemplArrayNom1Re[im][it][ir] = 0;
	    fMCTemplArrayTun1Re[im][it][ir] = 0;
	    fMCTemplArrayFin1Re[im][it][ir] = 0;
	  }
	  ////
	  // ANTINU
	  else if (s==kSKsRHCFCFV1RingMuLike){
	    fMCTemplArrayOriRHC1Rmu[im][it][ir] = 0;
	    fMCTemplArrayNomRHC1Rmu[im][it][ir] = 0;
	    fMCTemplArrayTunRHC1Rmu[im][it][ir] = 0;
	    fMCTemplArrayFinRHC1Rmu[im][it][ir] = 0;
	  }
	  else if (s==kSKsRHCFCFV1RingELike) {
	    fMCTemplArrayOriRHC1Re[im][it][ir] = 0;
	    fMCTemplArrayNomRHC1Re[im][it][ir] = 0;
	    fMCTemplArrayTunRHC1Re[im][it][ir] = 0;
	    fMCTemplArrayFinRHC1Re[im][it][ir] = 0;
	  }
	  ////

	}//ir
      }//it
      
      int n_mc_etrue_bins = kMCNEtrueBins;

      for(int itv=0; itv<n_mc_etrue_bins; itv++) {   // MODIFIED
	for(int irv=0; irv<n_xs_ereco_bins; irv++) {

	  if (s==kSKsFCFV1RingMuLike) fMCVariableTemplArray1Rmu[im][itv][irv] = 0;  // MODIFIED
	  else if(s==kSKsFCFV1RingELike) fMCVariableTemplArray1Re[im][itv][irv] = 0; // MODIFIED
	  else if (s==kSKsRHCFCFV1RingMuLike) fMCVariableTemplArrayRHC1Rmu[im][itv][irv] = 0; // ANTINU
	  else if(s==kSKsRHCFCFV1RingELike) fMCVariableTemplArrayRHC1Re[im][itv][irv] = 0;    // ANTINU

	}//irv
      }//itv
    }//im
    for(int ir=0; ir<n_ereco_bins; ir++) {
      
      if (s==kSKsFCFV1RingMuLike) fMCErecoSpectrumTotArray1Rmu[ir] = 0;             // MODIFIED
      else if (s==kSKsFCFV1RingELike) fMCErecoSpectrumTotArray1Re[ir] = 0;          // MODIFIED
      else if (s==kSKsRHCFCFV1RingMuLike) fMCErecoSpectrumTotArrayRHC1Rmu[ir] = 0;  // ANTINU
      else if (s==kSKsRHCFCFV1RingELike) fMCErecoSpectrumTotArrayRHC1Re[ir] = 0;    // ANTINU
      
      for(int im=0; im<knm; im++) {

	if (s==kSKsFCFV1RingMuLike) fMCErecoSpectrumArray1Rmu[im][ir] = 0;             // MODIFIED
	else if (s==kSKsFCFV1RingELike) fMCErecoSpectrumArray1Re [im][ir] = 0;         // MODIFIED
	else if (s==kSKsRHCFCFV1RingMuLike) fMCErecoSpectrumArrayRHC1Rmu[im][ir] = 0;  // ANTINU
	else if (s==kSKsRHCFCFV1RingELike) fMCErecoSpectrumArrayRHC1Re [im][ir] = 0;   // ANTINU

      }//im
    }//ir
    for(int it=0; it<n_etrue_bins; it++) {

      if (s==kSKsFCFV1RingMuLike) fMCEtrueSpectrumTotArray1Rmu[it] = 0;            // MODIFIED
      else if (s==kSKsFCFV1RingELike) fMCEtrueSpectrumTotArray1Re[it] = 0;         // MODIFIED
      else if (s==kSKsRHCFCFV1RingMuLike) fMCEtrueSpectrumTotArrayRHC1Rmu[it] = 0; // ANTINU
      else if (s==kSKsRHCFCFV1RingELike) fMCEtrueSpectrumTotArrayRHC1Re[it] = 0;   // ANTINU

    }//it
  }//is
  
  fMCTreeNuMu1Rmu      = 0;
  fMCTreeNuMuBar1Rmu   = 0;
  fMCTreeNuE1Rmu       = 0;
  fMCTreeNuEBar1Rmu    = 0;
  fMCTreeOscNuE1Rmu    = 0;
  
  fMCTreeNuMu1Re       = 0;
  fMCTreeNuMuBar1Re    = 0;
  fMCTreeNuE1Re        = 0;
  fMCTreeNuEBar1Re     = 0;
  fMCTreeOscNuE1Re     = 0;

  // ANTINU
  fMCTreeNuMuRHC1Rmu      = 0;
  fMCTreeNuMuBarRHC1Rmu   = 0;
  fMCTreeNuERHC1Rmu       = 0;
  fMCTreeNuEBarRHC1Rmu    = 0;
  fMCTreeOscNuEBarRHC1Rmu = 0;
  
  fMCTreeNuMuRHC1Re       = 0;
  fMCTreeNuMuBarRHC1Re    = 0;
  fMCTreeNuERHC1Re        = 0;
  fMCTreeNuEBarRHC1Re     = 0;
  fMCTreeOscNuEBarRHC1Re  = 0;
  ////  

  fMCPOTNuMu           = 0;
  fMCPOTNuMuBar        = 0;
  fMCPOTNuE            = 0;
  fMCPOTNuEBar         = 0;
  fMCPOTOscNuE         = 0;

  // ANTINU
  fMCPOTRHCNuMu           = 0;
  fMCPOTRHCNuMuBar        = 0;
  fMCPOTRHCNuE            = 0;
  fMCPOTRHCNuEBar         = 0;
  fMCPOTRHCOscNuEBar      = 0;
  ////

  fBrEreco             = 0;
  fBrEtrue             = 0;
  fBrGenCode           = 0;
  
  fNuFlux              = 0;
  fNuXSec              = 0;  
  fFSI                 = 0;  
  fSKDet               = 0;
  fSKDetFSI            = 0;
  
  arr_index = util::ArrIdx::Instance();

  cout << "Exp/SKMC >> Init complete" << endl;

}
//________________________________________________________________________________________
void SuperkMC::CleanUp(void)
{
  /*started causing segfaults???
  if(fMCTreeNuMu1Rmu    ) delete fMCTreeNuMu1Rmu;
  if(fMCTreeNuMuBar1Rmu ) delete fMCTreeNuMuBar1Rmu;
  if(fMCTreeNuE1Rmu     ) delete fMCTreeNuE1Rmu;
  if(fMCTreeNuEBar1Rmu  ) delete fMCTreeNuEBar1Rmu;
  if(fMCTreeOscNuE1Rmu  ) delete fMCTreeOscNuE1Rmu;
  */
  /*
  if(fMCTreeNuMu1Re     ) delete fMCTreeNuMu1Re;
  if(fMCTreeNuMuBar1Re  ) delete fMCTreeNuMuBar1Re;
  if(fMCTreeNuE1Re      ) delete fMCTreeNuE1Re;
  if(fMCTreeNuEBar1Re   ) delete fMCTreeNuEBar1Re;
  if(fMCTreeOscNuE1Re   ) delete fMCTreeOscNuE1Re;
  */
  /*
  if(fMCTreeNuMuRHC1Rmu    ) delete fMCTreeNuMuRHC1Rmu;     // ANTINU
  if(fMCTreeNuMuBarRHC1Rmu ) delete fMCTreeNuMuBarRHC1Rmu;  // ANTINU
  if(fMCTreeNuERHC1Rmu     ) delete fMCTreeNuERHC1Rmu;      // ANTINU
  if(fMCTreeNuEBarRHC1Rmu  ) delete fMCTreeNuEBarRHC1Rmu;   // ANTINU
  if(fMCTreeOscNuERHC1Rmu  ) delete fMCTreeOscNuERHC1Rmu;   // ANTINU
  */
  /*
  if(fMCTreeNuMuRHC1Re     ) delete fMCTreeNuMuRHC1Re;      // ANTINU
  if(fMCTreeNuMuBarRHC1Re  ) delete fMCTreeNuMuBarRHC1Re;   // ANTINU
  if(fMCTreeNuERHC1Re      ) delete fMCTreeNuERHC1Re;       // ANTINU
  if(fMCTreeNuEBarRHC1Re   ) delete fMCTreeNuEBarRHC1Re;    // ANTINU
  if(fMCTreeOscNuERHC1Re   ) delete fMCTreeOscNuERHC1Re;    // ANTINU
  */

  if(fMCTemplMECMartini) delete fMCTemplMECMartini;
  if(fMCTemplMECNieves ) delete fMCTemplMECNieves;
  if(fMCTemplRHCMECMartini) delete fMCTemplRHCMECMartini;  // ANTINU
  if(fMCTemplRHCMECNieves ) delete fMCTemplRHCMECNieves;   // ANTINU
  for(int is = 0; is < kns; is++) {
    for(int im = 0; im < knm; im++) {
      //already deleted Ori & Nom (the if doesn't work!)
      //if( fMCTemplOri     [is][im] ) delete fMCTemplOri     [is][im];
      //if( fMCTemplNom     [is][im] ) delete fMCTemplNom     [is][im];
#ifndef PRODUCTIONMODE
      if( fMCTemplTun     [is][im] ) delete fMCTemplTun     [is][im];
      if( fMCTemplSys     [is][im] ) delete fMCTemplSys     [is][im];
      if( fMCTemplFin     [is][im] ) delete fMCTemplFin     [is][im];
#endif
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
  
  //
  // Get MC trees and MC POT
  //
  fMCPOTHistNuMu    = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numu/pot");
  fMCPOTHistNuMuBar = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numubar/pot");
  fMCPOTHistNuE     = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nue/pot");
  fMCPOTHistNuEBar  = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nuebar/pot");
  fMCPOTHistOscNuE  = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/signal_nue/pot");

  fMCPOTNuMu    = fMCPOTHistNuMu->GetBinContent(2)    / fMCPOTHistNuMu->GetBinContent(1);
  fMCPOTNuMuBar = fMCPOTHistNuMuBar->GetBinContent(2) / fMCPOTHistNuMuBar->GetBinContent(1);
  fMCPOTNuE     = fMCPOTHistNuE->GetBinContent(2)     / fMCPOTHistNuE->GetBinContent(1);
  fMCPOTNuEBar  = fMCPOTHistNuEBar->GetBinContent(2)  / fMCPOTHistNuEBar->GetBinContent(1);
  fMCPOTOscNuE  = fMCPOTHistOscNuE->GetBinContent(2)  / fMCPOTHistOscNuE->GetBinContent(1);


  // ANTINU
  fMCPOTHistRHCNuMu    = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numu/pot");
  fMCPOTHistRHCNuMuBar = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numubar/pot");
  fMCPOTHistRHCNuE     = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nue/pot");
  fMCPOTHistRHCNuEBar  = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nuebar/pot");
  fMCPOTHistRHCOscNuEBar  = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/signal_nue/pot");

  fMCPOTRHCNuMu    = fMCPOTHistRHCNuMu->GetBinContent(2)    / fMCPOTHistRHCNuMu->GetBinContent(1);
  fMCPOTRHCNuMuBar = fMCPOTHistRHCNuMuBar->GetBinContent(2) / fMCPOTHistRHCNuMuBar->GetBinContent(1);
  fMCPOTRHCNuE     = fMCPOTHistRHCNuE->GetBinContent(2)     / fMCPOTHistRHCNuE->GetBinContent(1);
  fMCPOTRHCNuEBar  = fMCPOTHistRHCNuEBar->GetBinContent(2)  / fMCPOTHistRHCNuEBar->GetBinContent(1);
  fMCPOTRHCOscNuEBar  = fMCPOTHistRHCOscNuEBar->GetBinContent(2)  / fMCPOTHistRHCOscNuEBar->GetBinContent(1);
  ////

  // MODIFIED
  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsFCFV1RingMuLike) {
       fMCTreeNuMu1Rmu    = util::XML2TTree (xml,root_node,"mc/superk/nominal/numu/one_mu_ring");
       fMCTreeNuMuBar1Rmu = util::XML2TTree (xml,root_node,"mc/superk/nominal/numubar/one_mu_ring");
       fMCTreeNuE1Rmu     = util::XML2TTree (xml,root_node,"mc/superk/nominal/nue/one_mu_ring");
       fMCTreeNuEBar1Rmu  = util::XML2TTree (xml,root_node,"mc/superk/nominal/nuebar/one_mu_ring");
       fMCTreeOscNuE1Rmu  = util::XML2TTree (xml,root_node,"mc/superk/nominal/signal_nue/one_mu_ring");
       assert (fMCTreeNuMu1Rmu   );
       assert (fMCTreeNuMuBar1Rmu);
       assert (fMCTreeNuE1Rmu    );
       assert (fMCTreeNuEBar1Rmu );
       assert (fMCTreeOscNuE1Rmu );
     }
     else if(s == kSKsFCFV1RingELike) {
       fMCTreeNuMu1Re     = util::XML2TTree (xml,root_node,"mc/superk/nominal/numu/one_e_ring");
       fMCTreeNuMuBar1Re  = util::XML2TTree (xml,root_node,"mc/superk/nominal/numubar/one_e_ring");
       fMCTreeNuE1Re      = util::XML2TTree (xml,root_node,"mc/superk/nominal/nue/one_e_ring");
       fMCTreeNuEBar1Re   = util::XML2TTree (xml,root_node,"mc/superk/nominal/nuebar/one_e_ring");
       fMCTreeOscNuE1Re   = util::XML2TTree (xml,root_node,"mc/superk/nominal/signal_nue/one_e_ring");
       assert (fMCTreeNuMu1Re    );
       assert (fMCTreeNuMuBar1Re );
       assert (fMCTreeNuE1Re     );
       assert (fMCTreeNuEBar1Re  );
       assert (fMCTreeOscNuE1Re  );
     }
     ////
	

					
     // ANTINU
     else if(s == kSKsRHCFCFV1RingMuLike) {
       fMCTreeNuMuRHC1Rmu    = util::XML2TTree (xml,root_node,"mc/superk/nominal/numu/one_mu_ring");
       fMCTreeNuMuBarRHC1Rmu = util::XML2TTree (xml,root_node,"mc/superk/nominal/numubar/one_mu_ring");
       fMCTreeNuERHC1Rmu     = util::XML2TTree (xml,root_node,"mc/superk/nominal/nue/one_mu_ring");
       fMCTreeNuEBarRHC1Rmu  = util::XML2TTree (xml,root_node,"mc/superk/nominal/nuebar/one_mu_ring");
       fMCTreeOscNuEBarRHC1Rmu  = util::XML2TTree (xml,root_node,"mc/superk/nominal/signal_nue/one_mu_ring");
       assert (fMCTreeNuMuRHC1Rmu   );
       assert (fMCTreeNuMuBarRHC1Rmu);
       assert (fMCTreeNuERHC1Rmu    );
       assert (fMCTreeNuEBarRHC1Rmu );
       assert (fMCTreeOscNuEBarRHC1Rmu );
     }
     else if(s == kSKsRHCFCFV1RingELike) {
       fMCTreeNuMuRHC1Re     = util::XML2TTree (xml,root_node,"mc/superk/nominal/numu/one_e_ring");
       fMCTreeNuMuBarRHC1Re  = util::XML2TTree (xml,root_node,"mc/superk/nominal/numubar/one_e_ring");
       fMCTreeNuERHC1Re      = util::XML2TTree (xml,root_node,"mc/superk/nominal/nue/one_e_ring");
       fMCTreeNuEBarRHC1Re   = util::XML2TTree (xml,root_node,"mc/superk/nominal/nuebar/one_e_ring");
       fMCTreeOscNuEBarRHC1Re   = util::XML2TTree (xml,root_node,"mc/superk/nominal/signal_nue/one_e_ring");
       assert (fMCTreeNuMuRHC1Re    );
       assert (fMCTreeNuMuBarRHC1Re );
       assert (fMCTreeNuERHC1Re     );
       assert (fMCTreeNuEBarRHC1Re  );
       assert (fMCTreeOscNuEBarRHC1Re  );
     }
     ////
  }
  

#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> " << "POT {MC; numu}    = " << fMCPOTNuMu    << endl;
  cout << "Exp/SKMC >> " << "POT {MC; numubar} = " << fMCPOTNuMuBar << endl;
  cout << "Exp/SKMC >> " << "POT {MC; nue}     = " << fMCPOTNuE     << endl;
  cout << "Exp/SKMC >> " << "POT {MC; nuebar}  = " << fMCPOTNuEBar  << endl;
  cout << "Exp/SKMC >> " << "POT {MC; osc_nue} = " << fMCPOTOscNuE  << endl;
  
  // ANTINU
  cout << "RHC: Exp/SKMC >> " << "POT {MC; numu}    = " << fMCPOTRHCNuMu    << endl;
  cout << "RHC: Exp/SKMC >> " << "POT {MC; numubar} = " << fMCPOTRHCNuMuBar << endl;
  cout << "RHC: Exp/SKMC >> " << "POT {MC; nue}     = " << fMCPOTRHCNuE     << endl;
  cout << "RHC: Exp/SKMC >> " << "POT {MC; nuebar}  = " << fMCPOTRHCNuEBar  << endl;
  cout << "RHC: Exp/SKMC >> " << "POT {MC; osc_nue} = " << fMCPOTRHCOscNuEBar  << endl;
  ////

  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsFCFV1RingMuLike) {
       cout << "Exp/SKMC >> " << "Nev {MC; numu;     1Rmu} = " << fMCTreeNuMu1Rmu    -> GetEntries()  << endl;
       cout << "Exp/SKMC >> " << "Nev {MC; numubar;  1Rmu} = " << fMCTreeNuMuBar1Rmu -> GetEntries()  << endl;
       cout << "Exp/SKMC >> " << "Nev {MC; nue;      1Rmu} = " << fMCTreeNuE1Rmu     -> GetEntries()  << endl;
       cout << "Exp/SKMC >> " << "Nev {MC; nuebar;   1Rmu} = " << fMCTreeNuEBar1Rmu  -> GetEntries()  << endl;
       cout << "Exp/SKMC >> " << "Nev {MC; nue[osc]; 1Rmu} = " << fMCTreeOscNuE1Rmu  -> GetEntries()  << endl;
     }
     else if(s == kSKsFCFV1RingELike) {
       cout << "Exp/SKMC >> " << "Nev {MC; numu;     1Re } = " << fMCTreeNuMu1Re     -> GetEntries()  << endl;
       cout << "Exp/SKMC >> " << "Nev {MC; numubar;  1Re } = " << fMCTreeNuMuBar1Re  -> GetEntries()  << endl;
       cout << "Exp/SKMC >> " << "Nev {MC; nue;      1Re } = " << fMCTreeNuE1Re      -> GetEntries()  << endl;
       cout << "Exp/SKMC >> " << "Nev {MC; nuebar;   1Re } = " << fMCTreeNuEBar1Re   -> GetEntries()  << endl;
       cout << "Exp/SKMC >> " << "Nev {MC; nue[osc]; 1Re } = " << fMCTreeOscNuE1Re   -> GetEntries()  << endl;
     }

     // ANTINU
     else if(s == kSKsRHCFCFV1RingMuLike) {
       cout << "RHC: Exp/SKMC >> " << "Nev {MC; numu;     1Rmu} = " << fMCTreeNuMuRHC1Rmu    -> GetEntries()  << endl;
       cout << "RHC: Exp/SKMC >> " << "Nev {MC; numubar;  1Rmu} = " << fMCTreeNuMuBarRHC1Rmu -> GetEntries()  << endl;
       cout << "RHC: Exp/SKMC >> " << "Nev {MC; nue;      1Rmu} = " << fMCTreeNuERHC1Rmu     -> GetEntries()  << endl;
       cout << "RHC: Exp/SKMC >> " << "Nev {MC; nuebar;   1Rmu} = " << fMCTreeNuEBarRHC1Rmu  -> GetEntries()  << endl;
       cout << "RHC: Exp/SKMC >> " << "Nev {MC; nue[osc]; 1Rmu} = " << fMCTreeOscNuEBarRHC1Rmu  -> GetEntries()  << endl;
     }
     else if(s == kSKsRHCFCFV1RingELike) {
       cout << "RHC: Exp/SKMC >> " << "Nev {MC; numu;     1Re } = " << fMCTreeNuMuRHC1Re     -> GetEntries()  << endl;
       cout << "RHC: Exp/SKMC >> " << "Nev {MC; numubar;  1Re } = " << fMCTreeNuMuBarRHC1Re  -> GetEntries()  << endl;
       cout << "RHC: Exp/SKMC >> " << "Nev {MC; nue;      1Re } = " << fMCTreeNuERHC1Re      -> GetEntries()  << endl;
       cout << "RHC: Exp/SKMC >> " << "Nev {MC; nuebar;   1Re } = " << fMCTreeNuEBarRHC1Re   -> GetEntries()  << endl;
       cout << "RHC: Exp/SKMC >> " << "Nev {MC; nue[osc]; 1Re } = " << fMCTreeOscNuEBarRHC1Re   -> GetEntries()  << endl;
     }
     ////

     
  }
#endif//__show_init_mesg_info__

  // get MEC template
  fMCTemplMECMartini    = util::XML2TH2D  (xml,root_node,"mec/templates/martini");
  fMCTemplMECNieves     = util::XML2TH2D  (xml,root_node,"mec/templates/nieves");

  // ANTINU
  fMCTemplRHCMECMartini    = util::XML2TH2D  (xml,root_node,"mec/templates/martini");
  fMCTemplRHCMECNieves     = util::XML2TH2D  (xml,root_node,"mec/templates/nieves");
  ////

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

  // fill-in "special" MC templates (MEC, osc. numu CC)
  this->FillSpecialTemplates();

#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtOri, "nevents");
  //this->NevDebug(kSKMCtNom, "mcnorm");
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

#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtNom, "mcnorm_fluxtuned");
#endif

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
#ifndef PRODUCTIONMODE
       fMCTemplTun[is][im] = new TH2D(
          Form("nsk_tun_%d_%d",is,im), util::AsString(m).c_str(),
	  n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
       fMCTemplSys[is][im] = new TH2D(
          Form("nsk_sys_%d_%d",is,im), util::AsString(m).c_str(),
	  n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
       fMCTemplFin[is][im] = new TH2D(
          Form("nsk_fin_%d_%d",is,im), util::AsString(m).c_str(),
	  n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
#endif
       fMCTemplOri[is][im]->SetDirectory(0);
       fMCTemplOri[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fMCTemplOri[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
       fMCTemplNom[is][im]->SetDirectory(0);
       fMCTemplNom[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fMCTemplNom[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
#ifndef PRODUCTIONMODE
       fMCTemplTun[is][im]->SetDirectory(0);
       fMCTemplTun[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fMCTemplTun[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
       fMCTemplSys[is][im]->SetDirectory(0);
       fMCTemplSys[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fMCTemplSys[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
       fMCTemplFin[is][im]->SetDirectory(0);
       fMCTemplFin[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fMCTemplFin[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
#endif
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
  if(util::IsMEC(m) || util::IsOscNuMu(m) || util::IsOscNuMuBar(m) || util::IsOscNuEBar(m))
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

  TBranch * branchEreco   = mctree->GetBranch("erec");
  TBranch * branchEtrue   = mctree->GetBranch("enu");
  TBranch * branchGenCode = mctree->GetBranch("mode");
  TBranch * branchWght    = mctree->GetBranch("wgt");

  branchEreco   -> SetAddress(&fBrEreco);
  branchEtrue   -> SetAddress(&fBrEtrue);
  branchGenCode -> SetAddress(&fBrGenCode);
  branchWght    -> SetAddress(&fBrWght);  

  mctree->SetBranchStatus("*",0);
  mctree->SetBranchStatus("erec",1);
  mctree->SetBranchStatus("enu",1);
  mctree->SetBranchStatus("mode",1);
  mctree->SetBranchStatus("wgt",1);

  for(int iev = 0; iev < mctree->GetEntries(); iev++) 
  {
     mctree->GetEntry(iev);
     
     if(fBrEreco<0) continue;

     // use current entry for mode m template?
     bool use_it = util::UseEventToFillInTemplate(m,util::SuperkMode2SKMCPdg(m),fBrGenCode);
     if(!use_it) continue;

     double ereco_fill = fBrEreco/1000.;

     if(util::IsRHCSample(s))
       fBrWght = fNuFlux->SuperkFluxExtTuningWght(s,m,fBrEtrue);
     
#ifdef INCLUDEEVENTSOVERTHIRTYGEV
     ereco_fill = TMath::Min(29.9, ereco_fill);
#endif
     this -> Template(s,m,kSKMCtOri) -> Fill (ereco_fill, fBrEtrue);
     fMCVariableTempl[sidx][midx]    -> Fill (ereco_fill, fBrEtrue);
#ifndef FLUXWEIGHTSAREONE
     this -> Template(s,m,kSKMCtNom) -> Fill (ereco_fill, fBrEtrue, fBrWght * pot_wght);
#ifndef PRODUCTIONMODE
     this -> Template(s,m,kSKMCtTun) -> Fill (ereco_fill, fBrEtrue, fBrWght * pot_wght);
     this -> Template(s,m,kSKMCtSys) -> Fill (ereco_fill, fBrEtrue, fBrWght * pot_wght);
     this -> Template(s,m,kSKMCtFin) -> Fill (ereco_fill, fBrEtrue, fBrWght * pot_wght);
#endif//PRODUCTIONMODE
     
#else
     this -> Template(s,m,kSKMCtNom) -> Fill (ereco_fill, fBrEtrue, pot_wght);
#ifndef PRODUCTIONMODE
     this -> Template(s,m,kSKMCtTun) -> Fill (ereco_fill, fBrEtrue, pot_wght);
     this -> Template(s,m,kSKMCtSys) -> Fill (ereco_fill, fBrEtrue, pot_wght);
     this -> Template(s,m,kSKMCtFin) -> Fill (ereco_fill, fBrEtrue, pot_wght);
#endif//PRODUCTIONMODE
#endif//FLUXWEIGHTSAREONE
     
#ifdef __show_init_mesg_debug__
     if(fBrEtrue <= 0 || fBrEtrue >= this->Template(s,m,kSKMCtOri)->GetYaxis()->GetBinUpEdge(this->Template(s,m,kSKMCtOri)->GetNbinsY()))
       cout << "Exp/SKMC/CreateTemplate >> Etrue is out of template range " << AsString(s) << " " << AsString(m) << " " << fBrEtrue 
	    << " (Ereco = " << fBrEreco/1000 << ")" << endl;
#endif
#ifdef __show_init_mesg_debug2__
     if(fBrEreco/1000 <= 0)
       cout << "Exp/SKMC/CreateTemplate >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << fBrEreco/1000
	    << " (Etrue = " << fBrEtrue << ")" << endl;
     //#endif __show_init_mesg_debug2__
     if(fBrEreco/1000 >= this->Template(s,m,kSKMCtOri)->GetXaxis()->GetBinUpEdge(this->Template(s,m,kSKMCtOri)->GetNbinsX()))
       cout << "Exp/SKMC/CreateTemplate >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << fBrEreco/1000
	    << " (Etrue = " << fBrEtrue << ")" << endl;
#endif// __show_init_mesg_debug__

  }//iev

  if(sidx > -1 && midx > -1) {
    // loop over Etrue
    for(int it = 1; it <= this->Template(s,m,kSKMCtNom)->GetYaxis()->GetNbins(); it++) {
      // loop over Ereco
      for(int ir = 1; ir <= this->Template(s,m,kSKMCtNom)->GetXaxis()->GetNbins(); ir++) {

	
	if (s==kSKsFCFV1RingMuLike) fMCTemplArrayOri1Rmu[midx][it-1][ir-1] = (int)this->Template(s,m,kSKMCtOri)->GetBinContent(ir,it);
	else if (s==kSKsFCFV1RingELike)
	  fMCTemplArrayOri1Re[midx][it-1][ir-1] = (int)this->Template(s,m,kSKMCtOri)->GetBinContent(ir,it);   // MODIFIED
	// ANTINU
	else if (s==kSKsRHCFCFV1RingMuLike)
	  fMCTemplArrayOriRHC1Rmu[midx][it-1][ir-1] = (int)this->Template(s,m,kSKMCtOri)->GetBinContent(ir,it);
	else if (s==kSKsRHCFCFV1RingELike)
	  fMCTemplArrayOriRHC1Re[midx][it-1][ir-1] = (int)this->Template(s,m,kSKMCtOri)->GetBinContent(ir,it);   
	////

	double nev = this->Template(s,m,kSKMCtNom)->GetBinContent(ir,it);
	if (s==kSKsFCFV1RingMuLike) {
	  fMCTemplArrayNom1Rmu[midx][it-1][ir-1] = nev;
	  fMCTemplArrayTun1Rmu[midx][it-1][ir-1] = nev;
	  fMCTemplArrayFin1Rmu[midx][it-1][ir-1] = nev;
	}
	// MODIFIED
	else if (s==kSKsFCFV1RingELike) {
	  fMCTemplArrayNom1Re[midx][it-1][ir-1] = nev;
	  fMCTemplArrayTun1Re[midx][it-1][ir-1] = nev;
	  fMCTemplArrayFin1Re[midx][it-1][ir-1] = nev;
	}
	////

	// ANTINU
	else if (s==kSKsRHCFCFV1RingMuLike) {
	  fMCTemplArrayNomRHC1Rmu[midx][it-1][ir-1] = nev;
	  fMCTemplArrayTunRHC1Rmu[midx][it-1][ir-1] = nev;
	  fMCTemplArrayFinRHC1Rmu[midx][it-1][ir-1] = nev;
	}
	else if (s==kSKsRHCFCFV1RingELike) {
	  fMCTemplArrayNomRHC1Re[midx][it-1][ir-1] = nev;
	  fMCTemplArrayTunRHC1Re[midx][it-1][ir-1] = nev;
	  fMCTemplArrayFinRHC1Re[midx][it-1][ir-1] = nev;
	}
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
      
      // loop over variable Ereco
      for(int irv = 1; irv <= fMCVariableTempl[sidx][midx]->GetXaxis()->GetNbins(); irv++) {
	int nev = (int)fMCVariableTempl[sidx][midx]->GetBinContent(irv,itv);

	if (s==kSKsFCFV1RingMuLike) fMCVariableTemplArray1Rmu[midx][itv-1][irv-1] = nev;
	else if (s==kSKsFCFV1RingELike) fMCVariableTemplArray1Re[midx][itv-1][irv-1] = nev;   // MODIFIED
	else if (s==kSKsRHCFCFV1RingMuLike) fMCVariableTemplArrayRHC1Rmu[midx][itv-1][irv-1] = nev; // ANTINU
	else if (s==kSKsRHCFCFV1RingELike) fMCVariableTemplArrayRHC1Re[midx][itv-1][irv-1] = nev;   // ANTINU
      
	// Fill nonempty bin numbers.
	if (nev < conv::kASmallNum) {
	  continue;
	}
	
	fNonEmptyVariableBins[sidx][midx][i_ne_var] = conv::kBigIndex * irv + itv;
	
	double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
	double ErecoVarMax = util::XsecErecoBinEdges(s,irv)   - kOneMeV;

	int ErecoMinBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMin);
	int ErecoMaxBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMax);

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
    
    //
    // Add MEC template, if included in the list of modes
    // Special template that isnt't in the tree! 
    // Don't need to apply POT weights or flux weights!
    // Initalise the MEC mode template as empty
    // Here, just put the histogram in arrays (for quicker access)
    // and also fill the variable binning array (for skipping calls to find weights)
    //
    SuperkMode_t m = kSKmMEC;                        
    int midx = arr_index->SuperkModeIdx(s,kSKmMEC);
    bool incl_mec = ( midx != -1 );
    if(incl_mec)
      {
	// loop over Etrue
	cout << "Exp/SKMC/FillSpecialTemplates >> Filling MEC templates for " << util::AsString(s) << " sample" <<endl;
	for(int it = 1; it <= this->Template(s,kSKmMEC,kSKMCtNom)->GetYaxis()->GetNbins(); it++) {
	  double Etrue = this->Template(s,kSKmMEC,kSKMCtNom)->GetYaxis()->GetBinCenter(it);
	  // loop over Ereco
	  for(int ir = 1; ir <= this->Template(s,kSKmMEC,kSKMCtNom)->GetXaxis()->GetNbins(); ir++) {
	    double Ereco = this->Template(s,kSKmMEC,kSKMCtNom)->GetXaxis()->GetBinCenter(ir);
	 
	    //get nevents


	    // MODIFIED
	    double nev_martini;
	    double nev_nieves;
	    if (s==kSKsFCFV1RingMuLike || s==kSKsFCFV1RingELike){
	      nev_martini = fMCTemplMECMartini->GetBinContent(ir,it);
	      nev_nieves  = fMCTemplMECNieves-> GetBinContent(ir,it);
	    }
	    ////

	    // ANTINU
	    else if (s==kSKsRHCFCFV1RingMuLike || s==kSKsRHCFCFV1RingELike){
	      nev_martini = fMCTemplRHCMECMartini->GetBinContent(ir,it);
	      nev_nieves  = fMCTemplRHCMECNieves-> GetBinContent(ir,it);
	    }
	    ////


	    // fill the variable template (with both)
	    fMCVariableTempl[sidx][midx]->Fill(Ereco,Etrue,nev_martini);
	    fMCVariableTempl[sidx][midx]->Fill(Ereco,Etrue,nev_nieves);
	    //fill the mec arrays
	  


	    // MODIFIED
	    if(s==kSKsFCFV1RingMuLike || s==kSKsFCFV1RingELike){
	      fMCTemplArrayMECMartini[it-1][ir-1] = nev_martini;
	      fMCTemplArrayMECNieves [it-1][ir-1] = nev_nieves;
	    }
	    ////
	    // ANTINU
	    else if(s==kSKsRHCFCFV1RingMuLike || s==kSKsRHCFCFV1RingELike){
	      fMCTemplArrayRHCMECMartini[it-1][ir-1] = nev_martini;
	      fMCTemplArrayRHCMECNieves [it-1][ir-1] = nev_nieves;
	    }
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
	  
	  // loop over variable Ereco
	  for(int irv = 1; irv <= fMCVariableTempl[sidx][midx]->GetXaxis()->GetNbins(); irv++) {
	    int nev = TMath::CeilNint(fMCVariableTempl[sidx][midx]->GetBinContent(irv,itv));

	    if (s==kSKsFCFV1RingMuLike) fMCVariableTemplArray1Rmu[midx][itv-1][irv-1] = nev;
	    else if (s==kSKsFCFV1RingELike) fMCVariableTemplArray1Re[midx][itv-1][irv-1] = nev;   // MODIFIED
	    // ANTINU
	    else if (s==kSKsRHCFCFV1RingMuLike) fMCVariableTemplArrayRHC1Rmu[midx][itv-1][irv-1] = nev;
	    else if (s==kSKsRHCFCFV1RingELike) fMCVariableTemplArrayRHC1Re[midx][itv-1][irv-1] = nev;   
	    ////


	    // Fill nonempty bin numbers.
	    if (nev < conv::kASmallNum) {
	      continue;
	    }
	    
	    fNonEmptyVariableBins[sidx][midx][i_ne_var] = conv::kBigIndex * irv + itv;
	    
	    double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
	    double ErecoVarMax = util::XsecErecoBinEdges(s,irv)   - kOneMeV;
	    
	    int ErecoMinBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMin);
	    int ErecoMaxBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMax);

	    int i_ne_fixed = 0;
	    // loop over fixed Etrue to count n nonempty
	    for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
	      // loop over fixed Ereco to count
	      for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
		double nevf = fMCVariableTempl[sidx][midx]->GetBinContent(ir,it);
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
		double nevf = fMCVariableTempl[sidx][midx]->GetBinContent(ir,it);
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
      }//incl_mec
    
  
    //
    // Fill in `osc numu' templates from `numu' templates
    // and take into account the nue/numu flux weight
    // 
    
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    
    cout << "SKMC/FillSpecialTemplates >> Filling osc numu templates for " << util::AsString(s) << " sample" << endl;
    const int n = 4;
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
	double flux_nue_ori  = fNuFlux->SuperkFlux (kPdgNuE,  Etrue, kFlxNominal   );
	double flux_nue_nom  = fNuFlux->SuperkFlux (kPdgNuE,  Etrue, kFlxExtTuning ); 
	double flux_numu_ori = fNuFlux->SuperkFlux (kPdgNuMu, Etrue, kFlxNominal   );
	double flux_numu_nom = fNuFlux->SuperkFlux (kPdgNuMu, Etrue, kFlxExtTuning ); 
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
#ifndef PRODUCTIONMODE
	  this -> Template(s,tgt_modes[i],kSKMCtTun) -> Fill (Ereco,Etrue, weighted_nev_nom);
	  this -> Template(s,tgt_modes[i],kSKMCtSys) -> Fill (Ereco,Etrue, weighted_nev_nom);
	  this -> Template(s,tgt_modes[i],kSKMCtFin) -> Fill (Ereco,Etrue, weighted_nev_nom);
#endif//PRODUCTIONMODE
	  int nev_tgt_ori = TMath::CeilNint(tgt_ori->GetBinContent(ir,it));
	  double nev_tgt_nom = tgt_nom->GetBinContent(ir,it);
	  if(s==kSKsFCFV1RingMuLike){
	    fMCTemplArrayOri1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_ori;
	    fMCTemplArrayNom1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayTun1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayFin1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	  }
	  
	  // MODIFIED
	  else if(s==kSKsFCFV1RingELike){
	    fMCTemplArrayOri1Re[tgt_midx][it-1][ir-1] = nev_tgt_ori;
	    fMCTemplArrayNom1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayTun1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayFin1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	  }
	  ////
	  // ANTINU
	  else if(s==kSKsRHCFCFV1RingMuLike){
	    fMCTemplArrayOriRHC1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_ori;
	    fMCTemplArrayNomRHC1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayTunRHC1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayFinRHC1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	  }
	  else if(s==kSKsRHCFCFV1RingELike){
	    fMCTemplArrayOriRHC1Re[tgt_midx][it-1][ir-1] = nev_tgt_ori;
	    fMCTemplArrayNomRHC1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayTunRHC1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayFinRHC1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	  }
	  ////
	  

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
	int nev = 0;
	for(int irv = 1; irv <= fMCVariableTempl[sidx][tgt_midx]->GetXaxis()->GetNbins(); irv++) {


	  if(s==kSKsFCFV1RingMuLike){
	    nev = fMCVariableTemplArray1Rmu[src_midx][itv-1][irv-1];
	    fMCVariableTemplArray1Rmu[tgt_midx][itv-1][irv-1] = nev;
	  }
	  // MODIFIED
	  else if(s==kSKsFCFV1RingELike){
	    nev = fMCVariableTemplArray1Re[src_midx][itv-1][irv-1];
	    fMCVariableTemplArray1Re[tgt_midx][itv-1][irv-1] = nev;
	  }
	  ////
	  // ANTINU
	  else if(s==kSKsRHCFCFV1RingMuLike){
	    nev = fMCVariableTemplArrayRHC1Rmu[src_midx][itv-1][irv-1];
	    fMCVariableTemplArrayRHC1Rmu[tgt_midx][itv-1][irv-1] = nev;
	  }
	  else if(s==kSKsRHCFCFV1RingELike){
	    nev = fMCVariableTemplArrayRHC1Re[src_midx][itv-1][irv-1];
	    fMCVariableTemplArrayRHC1Re[tgt_midx][itv-1][irv-1] = nev;
	  }
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
	// perhaps we should rename kFlxNominal -> kFlxOriginal, for consistency).
	double flux_nue_ori  = fNuFlux->SuperkFlux (kPdgNuEBar,  Etrue, kFlxNominal   );
	double flux_nue_nom  = fNuFlux->SuperkFlux (kPdgNuEBar,  Etrue, kFlxExtTuning ); 
	double flux_numu_ori = fNuFlux->SuperkFlux (kPdgNuMuBar, Etrue, kFlxNominal   );
	double flux_numu_nom = fNuFlux->SuperkFlux (kPdgNuMuBar, Etrue, kFlxExtTuning ); 
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
#ifndef PRODUCTIONMODE
	  this -> Template(s,tgt_modes_numubar[i],kSKMCtTun) -> Fill (Ereco,Etrue, weighted_nev_nom);
	  this -> Template(s,tgt_modes_numubar[i],kSKMCtSys) -> Fill (Ereco,Etrue, weighted_nev_nom);
	  this -> Template(s,tgt_modes_numubar[i],kSKMCtFin) -> Fill (Ereco,Etrue, weighted_nev_nom);
#endif//PRODUCTIONMODE
	  int    nev_tgt_ori = TMath::CeilNint(tgt_ori->GetBinContent(ir,it));
	  double nev_tgt_nom = tgt_nom->GetBinContent(ir,it);


	  if(s==kSKsFCFV1RingMuLike){
	    fMCTemplArrayOri1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_ori;
	    fMCTemplArrayNom1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayTun1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayFin1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	  }
	  // MODIFIED
	  else if(s==kSKsFCFV1RingELike){
	    fMCTemplArrayOri1Re[tgt_midx][it-1][ir-1] = nev_tgt_ori;
	    fMCTemplArrayNom1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayTun1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayFin1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	  }
	  ////
	  // ANTINU
	  else if(s==kSKsRHCFCFV1RingMuLike){
	    fMCTemplArrayOriRHC1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_ori;
	    fMCTemplArrayNomRHC1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayTunRHC1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayFinRHC1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	  }
	  else if(s==kSKsRHCFCFV1RingELike){
	    fMCTemplArrayOriRHC1Re[tgt_midx][it-1][ir-1] = nev_tgt_ori;
	    fMCTemplArrayNomRHC1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayTunRHC1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayFinRHC1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	  }
	  ////

  
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
	  int nev = 0;
	

	  if(s==kSKsFCFV1RingMuLike){
	    nev = fMCVariableTemplArray1Rmu[src_midx][itv-1][irv-1];
	    fMCVariableTemplArray1Rmu[tgt_midx][itv-1][irv-1] = nev;
	  }
	  // MODIFIED
	  else if(s==kSKsFCFV1RingELike){
	    nev = fMCVariableTemplArray1Re[src_midx][itv-1][irv-1];
	    fMCVariableTemplArray1Re[tgt_midx][itv-1][irv-1] = nev;
	  }
	  ////
	  // ANTINU
	  else if(s==kSKsRHCFCFV1RingMuLike){
	    nev = fMCVariableTemplArrayRHC1Rmu[src_midx][itv-1][irv-1];
	    fMCVariableTemplArrayRHC1Rmu[tgt_midx][itv-1][irv-1] = nev;
	  }
	  else if(s==kSKsRHCFCFV1RingELike){
	    nev = fMCVariableTemplArrayRHC1Re[src_midx][itv-1][irv-1];
	    fMCVariableTemplArrayRHC1Re[tgt_midx][itv-1][irv-1] = nev;
	  }
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
	double flux_nue_ori  = fNuFlux->SuperkFlux (kPdgNuEBar,  Etrue, kFlxNominal   );
	double flux_nue_nom  = fNuFlux->SuperkFlux (kPdgNuEBar,  Etrue, kFlxExtTuning ); 
	double flux_numu_ori = fNuFlux->SuperkFlux (kPdgNuMuBar, Etrue, kFlxNominal   );
	double flux_numu_nom = fNuFlux->SuperkFlux (kPdgNuMuBar, Etrue, kFlxExtTuning ); 
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
#ifndef PRODUCTIONMODE
	  this -> Template(s,tgt_modes_nuebar[i],kSKMCtTun) -> Fill (Ereco,Etrue, weighted_nev_nom);
	  this -> Template(s,tgt_modes_nuebar[i],kSKMCtSys) -> Fill (Ereco,Etrue, weighted_nev_nom);
	  this -> Template(s,tgt_modes_nuebar[i],kSKMCtFin) -> Fill (Ereco,Etrue, weighted_nev_nom);
#endif//PRODUCTIONMODE
	  int    nev_tgt_ori = TMath::CeilNint(tgt_ori->GetBinContent(ir,it));
	  double nev_tgt_nom = tgt_nom->GetBinContent(ir,it);
	  

	  
	  if(s==kSKsFCFV1RingMuLike){
	    fMCTemplArrayOri1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_ori;
	    fMCTemplArrayNom1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayTun1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayFin1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	  }
	  // MODIFIED
	  else if(s==kSKsFCFV1RingELike){
	    fMCTemplArrayOri1Re[tgt_midx][it-1][ir-1] = nev_tgt_ori;
	    fMCTemplArrayNom1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayTun1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayFin1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	  }
	  ////
	  // ANTINU
	  else if(s==kSKsRHCFCFV1RingMuLike){
	    fMCTemplArrayOriRHC1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_ori;
	    fMCTemplArrayNomRHC1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayTunRHC1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayFinRHC1Rmu[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	  }
	  else if(s==kSKsRHCFCFV1RingELike){
	    fMCTemplArrayOriRHC1Re[tgt_midx][it-1][ir-1] = nev_tgt_ori;
	    fMCTemplArrayNomRHC1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayTunRHC1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	    fMCTemplArrayFinRHC1Re[tgt_midx][it-1][ir-1] = nev_tgt_nom;
	  }
	  ////

	  
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
	  int nev = 0;

	  if(s==kSKsFCFV1RingMuLike){
	    nev = fMCVariableTemplArray1Rmu[src_midx][itv-1][irv-1];
	    fMCVariableTemplArray1Rmu[tgt_midx][itv-1][irv-1] = nev;
	  }
	  // MODIFIED
	  else if(s==kSKsFCFV1RingELike){
	    nev = fMCVariableTemplArray1Re[src_midx][itv-1][irv-1];
	    fMCVariableTemplArray1Re[tgt_midx][itv-1][irv-1] = nev;
	  }
	  ////
	  // ANTINU
	  else if(s==kSKsRHCFCFV1RingMuLike){
	    nev = fMCVariableTemplArrayRHC1Rmu[src_midx][itv-1][irv-1];
	    fMCVariableTemplArrayRHC1Rmu[tgt_midx][itv-1][irv-1] = nev;
	  }
	  else if(s==kSKsRHCFCFV1RingELike){
	    nev = fMCVariableTemplArrayRHC1Re[src_midx][itv-1][irv-1];
	    fMCVariableTemplArrayRHC1Re[tgt_midx][itv-1][irv-1] = nev;
	  }
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
  // init
  for(int im=0; im < knm; im++) {
    for(int it = 0; it < kNEtrueBins1Rmu; it++) {
      for(int ir = 0; ir < kNErecoBins1Rmu; ir++) {
	fMCStatFracErr1Rmu[im][it][ir] = 0;
      }
    }
    for(int it = 0; it < kNEtrueBins1Re; it++) {
      for(int ir = 0; ir < kNErecoBins1Re; ir++) {
	fMCStatFracErr1Re[im][it][ir] = 0;
      }
    }
    // ANTINU
    for(int it = 0; it < kRHCNEtrueBins1Rmu; it++) {
      for(int ir = 0; ir < kRHCNErecoBins1Rmu; ir++) {
	fMCStatFracErrRHC1Rmu[im][it][ir] = 0;
      }
    }
    for(int it = 0; it < kRHCNEtrueBins1Re; it++) {
      for(int ir = 0; ir < kRHCNErecoBins1Re; ir++) {
	fMCStatFracErrRHC1Re[im][it][ir] = 0;
      }
    }
    ////
     
  }
  
  // set MC frac err in each bin
  for(int is = 0; is < kns; is++) 
  {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
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

             if(s == kSKsFCFV1RingMuLike)         { fMCStatFracErr1Rmu[im][it-1][ir-1] = fracerr; }
             else if(s == kSKsFCFV1RingELike)     { fMCStatFracErr1Re [im][it-1][ir-1] = fracerr; } // MODIFIED
             else if(s == kSKsRHCFCFV1RingMuLike) { fMCStatFracErrRHC1Rmu[im][it-1][ir-1] = fracerr; } // ANTINU
             else if(s == kSKsRHCFCFV1RingELike)  { fMCStatFracErrRHC1Re [im][it-1][ir-1] = fracerr; } // ANTINU
	     
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

       // MODIFIED
       int n_ereco_bins = kNErecoBins[is];
       int n_etrue_bins = kNEtrueBins[is];

       for(int ir = 1; ir <= n_ereco_bins; ir++) {
          for(int it = 1; it <= n_etrue_bins; it++) {

            if (s==kSKsFCFV1RingMuLike)  
            { 
               double nev_old_bin  = fMCTemplArrayFin1Rmu[im][it-1][ir-1]; 
               double frac_error   = fMCStatFracErr1Rmu  [im][it-1][ir-1];
               double error        = frac_error * nev_old_bin;
               double nev_new_bin  = rndgen.Gaus(nev_old_bin, error);
               nev_new_bin = TMath::Max(0.,nev_new_bin);
               nev_old += nev_old_bin;
               nev_new += nev_new_bin;
               fMCTemplArrayFin1Rmu[im][it-1][ir-1] = nev_new_bin;
            }
	    // MODIFIED
            else if (s==kSKsFCFV1RingELike){
               double nev_old_bin  = fMCTemplArrayFin1Re[im][it-1][ir-1]; 
               double frac_error   = fMCStatFracErr1Re  [im][it-1][ir-1];
               double error        = frac_error * nev_old_bin;
               double nev_new_bin  = rndgen.Gaus(nev_old_bin, error);
               nev_new_bin = TMath::Max(0.,nev_new_bin);
               nev_old += nev_old_bin;
               nev_new += nev_new_bin;
               fMCTemplArrayFin1Re[im][it-1][ir-1] = nev_new_bin;
            }
	    ////
	    // ANTINU
            else if (s==kSKsRHCFCFV1RingMuLike)  
            { 
               double nev_old_bin  = fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]; 
               double frac_error   = fMCStatFracErrRHC1Rmu  [im][it-1][ir-1];
               double error        = frac_error * nev_old_bin;
               double nev_new_bin  = rndgen.Gaus(nev_old_bin, error);
               nev_new_bin = TMath::Max(0.,nev_new_bin);
               nev_old += nev_old_bin;
               nev_new += nev_new_bin;
               fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1] = nev_new_bin;
            }
	    else if (s==kSKsRHCFCFV1RingELike){
               double nev_old_bin  = fMCTemplArrayFinRHC1Re[im][it-1][ir-1]; 
               double frac_error   = fMCStatFracErrRHC1Re  [im][it-1][ir-1];
               double error        = frac_error * nev_old_bin;
               double nev_new_bin  = rndgen.Gaus(nev_old_bin, error);
               nev_new_bin = TMath::Max(0.,nev_new_bin);
               nev_old += nev_old_bin;
               nev_new += nev_new_bin;
               fMCTemplArrayFinRHC1Re[im][it-1][ir-1] = nev_new_bin;
            }
	    ////

          }//it
       }//ir

       if(shape_only) {
          double scale = (nev_new > 0.) ? nev_old/nev_new : 0.;
          for(int ir = 1; ir <= n_ereco_bins; ir++) {
             for(int it = 1; it <= n_etrue_bins; it++) {
              
	       

	       if (s==kSKsFCFV1RingMuLike) { 
                  double nev  = fMCTemplArrayFin1Rmu[im][it-1][ir-1]; 
                  fMCTemplArrayFin1Rmu[im][it-1][ir-1] = nev*scale;
               }
               // MODIFIED
	       else if (s==kSKsFCFV1RingELike){
                  double nev  = fMCTemplArrayFin1Re[im][it-1][ir-1]; 
                  fMCTemplArrayFin1Re[im][it-1][ir-1] = nev*scale;
               }
	       // ANTINU
	       else if (s==kSKsRHCFCFV1RingMuLike) { 
                  double nev  = fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]; 
                  fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1] = nev*scale;
               }
	       else if (s==kSKsRHCFCFV1RingELike){
                  double nev  = fMCTemplArrayFinRHC1Re[im][it-1][ir-1]; 
                  fMCTemplArrayFinRHC1Re[im][it-1][ir-1] = nev*scale;
               }
	       ////

	       
             }//it
          }//ir
       }//shape_only?

     }//m
  }//s

}

//________________________________________________________________________________________
TH1D* SuperkMC::CurrErecoSpectrum(SuperkSample_t s) const                               
{ 
  int sidx = arr_index->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  return fMCErecoSpectrumTot[sidx];
}
//________________________________________________________________________________________
TH1D* SuperkMC::CurrErecoSpectrum(SuperkSample_t s, SuperkMode_t m) const 
{ 
  int sidx = arr_index->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = arr_index->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }            
  return fMCErecoSpectrum[sidx][midx]; 
}
//________________________________________________________________________________________
double* SuperkMC::CurrErecoSpectrumArray(SuperkSample_t s)                               
{ 
  int sidx = arr_index->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  
  // ANTINU 
  if (s==kSKsRHCFCFV1RingMuLike) return fMCErecoSpectrumTotArrayRHC1Rmu;
  else if (s==kSKsRHCFCFV1RingELike) return fMCErecoSpectrumTotArrayRHC1Re; 
  else if (s==kSKsFCFV1RingELike) return fMCErecoSpectrumTotArray1Re;
  else return fMCErecoSpectrumTotArray1Rmu;
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

  // ANTINU 
  if (s==kSKsRHCFCFV1RingMuLike) return fMCErecoSpectrumArrayRHC1Rmu[midx];
  else if (s==kSKsRHCFCFV1RingELike) return fMCErecoSpectrumArrayRHC1Re[midx];
  else if (s==kSKsFCFV1RingELike) return fMCErecoSpectrumArray1Re[midx];
  else return fMCErecoSpectrumArray1Rmu[midx];
  
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
#ifndef PRODUCTIONMODE
   case (kSKMCtTun) : 
     h = fMCTemplTun[sidx][midx];
     break;
   case (kSKMCtSys) :
     h = fMCTemplSys[sidx][midx];
     break;
   case (kSKMCtFin) : 
     h = fMCTemplFin[sidx][midx];
     break;
#endif
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
void SuperkMC::ApplyFluxExtTuningToTemplates(void) 
{
// Take into account the flux tuning based on external (NA61,...) data.
// This method reweights/overwrites the nominal templates & it is run just once after
// filling-in the MC templates from the input mini-DSTs.

#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> Applying flux tuning based on external data (NA61,...)." << endl;
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
#ifndef PRODUCTIONMODE
       // get corresponding template
       TH2D * mctempl = this -> Template(s,m,kSKMCtNom);
       if(!mctempl) continue; 
#endif

       int n_ereco_bins = kNErecoBins[is];
       int n_etrue_bins = kNEtrueBins[is];
       
       // loop over Ereco
       for(int ir = 1; ir <= n_ereco_bins; ir++) {
         // loop over Etrue
         for(int it = 1; it <= n_etrue_bins; it++) {
	   //check for empty bins
	   

	   if (s==kSKsFCFV1RingMuLike) {
	     if(fMCTemplArrayOri1Rmu[im][it-1][ir-1] <= 0) continue;
	   }
	   // MODIFIED
	   else if (s==kSKsFCFV1RingELike){
	     if(fMCTemplArrayOri1Re[im][it-1][ir-1] <= 0) continue;
	   }
	   ////
	   // ANTINU
	   else if (s==kSKsRHCFCFV1RingMuLike) {
	     if(fMCTemplArrayOriRHC1Rmu[im][it-1][ir-1] <= 0) continue;
	   }
	   else if (s==kSKsRHCFCFV1RingELike){
	     if(fMCTemplArrayOriRHC1Re[im][it-1][ir-1] <= 0) continue;
	   }
	   ////



	   double Etrue = util::EtrueBinEdges(s,it-1) + kOneMeV;
           // get number of events
           double nev = 0;

	   if (s==kSKsFCFV1RingMuLike) nev = fMCTemplArrayNom1Rmu[im][it-1][ir-1];
	   else if (s==kSKsFCFV1RingELike) nev = fMCTemplArrayNom1Re[im][it-1][ir-1]; // MODIFIED
	   // ANTINU
	   else if (s==kSKsRHCFCFV1RingMuLike) nev = fMCTemplArrayNomRHC1Rmu[im][it-1][ir-1];
	   else if (s==kSKsRHCFCFV1RingELike) nev = fMCTemplArrayNomRHC1Re[im][it-1][ir-1]; // MODIFIED
	   ////

	   nev = TMath::Max(0.,nev);
           // get tuning weights
           double wght = fNuFlux->SuperkFluxExtTuningWght(s,m,Etrue);
	   // update nominal MC template
           nev *= wght;
#ifdef NO_EVENTS_CAN_BE_NEGATIVE
	     if(nev < 0) {
#ifdef __show_negative_events_mesg_debug__
	       cout << "Exp/SKMC/ApplyFluxExtTuningToTemplates >> set negative number of events " << nev
		    << " to " << 0 << "\t" 
		    << " Ereco: " << util::ErecoBinEdges(s,ir-1) + kOneMeV 
		    << " EtrueL " << Etrue << endl;
#endif
	       nev = 0;
	     }
#endif
#ifndef PRODUCTIONMODE
           mctempl->SetBinContent(ir,it,nev);
           // initialize "tun", "sys" & "fin" templates to "nom" values
           this->Template(s,m,kSKMCtTun)->SetBinContent(ir,it,nev);
           this->Template(s,m,kSKMCtSys)->SetBinContent(ir,it,nev);
           this->Template(s,m,kSKMCtFin)->SetBinContent(ir,it,nev);
#endif


	   if (s==kSKsFCFV1RingMuLike) {
	     fMCTemplArrayNom1Rmu[im][it-1][ir-1] = nev;
	     fMCTemplArrayTun1Rmu[im][it-1][ir-1] = nev;
	     fMCTemplArrayFin1Rmu[im][it-1][ir-1] = nev;
	   }
	   // MODIFIED
	   else if (s==kSKsFCFV1RingELike) {
	     fMCTemplArrayNom1Re[im][it-1][ir-1] = nev;
	     fMCTemplArrayTun1Re[im][it-1][ir-1] = nev;
	     fMCTemplArrayFin1Re[im][it-1][ir-1] = nev;
	   }
	   ////
	   // ANTINU
	   else if (s==kSKsRHCFCFV1RingMuLike) {
	     fMCTemplArrayNomRHC1Rmu[im][it-1][ir-1] = nev;
	     fMCTemplArrayTunRHC1Rmu[im][it-1][ir-1] = nev;
	     fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1] = nev;
	   }
	   else if (s==kSKsRHCFCFV1RingELike) {
	     fMCTemplArrayNomRHC1Re[im][it-1][ir-1] = nev;
	     fMCTemplArrayTunRHC1Re[im][it-1][ir-1] = nev;
	     fMCTemplArrayFinRHC1Re[im][it-1][ir-1] = nev;
	   }
	   ////


         }//it
       }//ir
    }//im
  }//is
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
#ifndef PRODUCTIONMODE
      // get corresponding template
      TH2D * mctempl_tgt = this -> Template(s,m,kSKMCtTun);
      if(!mctempl_tgt) continue; 
#endif

#if BANFF_BINS == 4
      double x_mec = fXMEC;
      if(util::IsMEC(m))
	x_mec = fNuXSec->MECBestFitValue();
#endif

#ifdef __show_syst_var_wght_mesg_debug__ 
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
	if(s == kSKsFCFV1RingMuLike) {
	  if(fMCVariableTemplArray1Rmu[im][itv-1][irv-1] <= 0) continue;
	}
	// MODIFIED
	else if(s == kSKsFCFV1RingELike) {
	  if(fMCVariableTemplArray1Re[im][itv-1][irv-1] <= 0) continue;
	}
	////
	// ANTINU
	else if(s == kSKsRHCFCFV1RingMuLike) {
	  if(fMCVariableTemplArrayRHC1Rmu[im][itv-1][irv-1] <= 0) continue;
	}
	else if(s == kSKsRHCFCFV1RingELike) {
	  if(fMCVariableTemplArrayRHC1Re[im][itv-1][irv-1] <= 0) continue;
	}
	////     
	
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



	  if (s==kSKsFCFV1RingMuLike) {
#if BANFF_BINS == 4
	    if(!util::IsMEC(m))
#endif
	      if(fMCTemplArrayOri1Rmu[im][it-1][ir-1] <= 0) continue;
	  }
	  // MODIFIED
	  else if (s==kSKsFCFV1RingELike){
	    if(fMCTemplArrayOri1Re[im][it-1][ir-1] <= 0) continue;
	  }
	  ////
	  // ANTINU
	  else if (s==kSKsRHCFCFV1RingMuLike) {
#if BANFF_BINS == 4
	    if(!util::IsMEC(m))
#endif
	      if(fMCTemplArrayOriRHC1Rmu[im][it-1][ir-1] <= 0) continue;
	  }
	  else if (s==kSKsRHCFCFV1RingELike){
	    if(fMCTemplArrayOriRHC1Re[im][it-1][ir-1] <= 0) continue;
	  }
	  ////



	  //double Etrue = util::EtrueBinEdges(s,it-1) + kOneMeV;

	  // get number of events in current bin of source MC template
	  double nev = 0;
#if BANFF_BINS == 4
	  if(util::IsMEC(m)) {
	    

	    
	    // MODIFIED
	    double nev_martini;
	    double nev_nieves;
	    if (s==kSKsFCFV1RingMuLike || s==kSKsFCFV1RingELike){
	      nev_martini = fMCTemplArrayMECMartini[it-1][ir-1];
	      nev_nieves  = fMCTemplArrayMECNieves [it-1][ir-1];
	    }
	    ////
	    // ANTINU
	    else if (s==kSKsRHCFCFV1RingMuLike || s==kSKsRHCFCFV1RingELike){
	      nev_martini = fMCTemplArrayRHCMECMartini[it-1][ir-1];
	      nev_nieves  = fMCTemplArrayRHCMECNieves [it-1][ir-1];
	    }
	    ////



	    //-1 nominal -> 0 nieves
	    if(x_mec < 0) {
	      nev  = nev_nieves * (1 + x_mec);
	    }
	    // 0 nieves -> +1 martini
	    else {
	      nev += TMath::Min(nev_martini, nev_martini * x_mec);
	      nev += nev_nieves * (1 - x_mec);
	    }
	  }//MEC
	  else 
#endif
	    
	    
	    if (s==kSKsFCFV1RingMuLike) nev = fMCTemplArrayNom1Rmu[im][it-1][ir-1];
	    else if (s==kSKsFCFV1RingELike) nev = fMCTemplArrayNom1Re[im][it-1][ir-1];  // MODIFIED
	    else if (s==kSKsRHCFCFV1RingMuLike) nev = fMCTemplArrayNomRHC1Rmu[im][it-1][ir-1]; // ANTINU
	    else if (s==kSKsRHCFCFV1RingELike) nev = fMCTemplArrayNomRHC1Re[im][it-1][ir-1];   // ANTINU
	  
	  
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

#ifndef PRODUCTIONMODE
	  mctempl_tgt->SetBinContent(ir,it,nev_tuned);
	  // initialize "sys" & "fin" templates to "tun" values
	  this->Template(s,m,kSKMCtSys)->SetBinContent(ir,it,nev_tuned);
	  this->Template(s,m,kSKMCtFin)->SetBinContent(ir,it,nev_tuned);
#endif

	  
	  if (s==kSKsFCFV1RingMuLike) {
	    fMCTemplArrayTun1Rmu[im][it-1][ir-1] = nev_tuned;
	    fMCTemplArrayFin1Rmu[im][it-1][ir-1] = nev_tuned;
	  }
	  // MODIFIED
	  else if (s==kSKsFCFV1RingELike){
	    fMCTemplArrayTun1Re[im][it-1][ir-1] = nev_tuned;
	    fMCTemplArrayFin1Re[im][it-1][ir-1] = nev_tuned;
	  }
	  ////
	  // ANTINU
	  else if (s==kSKsRHCFCFV1RingMuLike) {
	    fMCTemplArrayTunRHC1Rmu[im][it-1][ir-1] = nev_tuned;
	    fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1] = nev_tuned;
	  }
	  else if (s==kSKsRHCFCFV1RingELike){
	    fMCTemplArrayTunRHC1Re[im][it-1][ir-1] = nev_tuned;
	    fMCTemplArrayFinRHC1Re[im][it-1][ir-1] = nev_tuned;
	  }
	  ////

#ifdef __show_syst_var_wght_mesg_debug__ 
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
#endif//__show_syst_var_wght_mesg_debug__ 

	}//ifb
      }//ivb
    }//im
  }//is

  //fFSI    -> BuildFracMCTemplates(fMCTemplTun);

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
#ifdef RENORMALISE_MEC
  this->ApplyMECRenormalisationToTemplates();
#endif
#ifndef APPLY_ESCALE_AT_SPECTRUM_STAGE
  this->ApplySystematicsToTemplates_NonMultiplicative();
#endif
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
#ifndef PRODUCTIONMODE
      // Get corresponding templates
      TH2D * mctempl_tgt = this -> Template(s,m,kSKMCtSys);
      if(!mctempl_tgt) continue;
#endif

#if BANFF_BINS == 4
      double x_mec = fXMEC;
      if(util::IsMEC(m))
	x_mec = fNuXSec->CurrentValue(kNuisBANFFMECNorm);       // GUARDARE
#endif
      const double * variableMCEtrueBinEdges = util::VariableMCEtrueBinEdges(s,m);

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

	//skip if the big bin is empty

	


	if(s == kSKsFCFV1RingMuLike) assert(fMCVariableTemplArray1Rmu[im][itv-1][irv-1] > 0);
	else if(s == kSKsFCFV1RingELike) assert(fMCVariableTemplArray1Re[im][itv-1][irv-1] > 0); // MODIFIED
	else if(s == kSKsRHCFCFV1RingMuLike) assert(fMCVariableTemplArrayRHC1Rmu[im][itv-1][irv-1] > 0); // ANTINU
	else if(s == kSKsRHCFCFV1RingELike) assert(fMCVariableTemplArrayRHC1Re[im][itv-1][irv-1] > 0); // ANTINU




	double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;

	// Plug-in BANFF & other xsec systematics 
	double xsec_wght = fNuXSec->SystWght(s,m,ErecoVarMin,EtrueVarMin); 

	double wght = 1;

#ifdef SKFSITOGETHER
	// Plug-in SK efficiency + FSI uncertainties
	double skfsi_wght = fSKDetFSI->SystWght(s,m,ErecoVarMin,EtrueVarMin);
	// Combined weight					  
	wght = skfsi_wght * flux_wght * xsec_wght;
#else
	// Plug-in SK efficiency uncertainty due to fiducial volume, reduction chain, OD cut,
	// ring-counting, PID, momentum cut and decay e- cut.
	// The ring-counting error for CCQE has Ereco-dependence.
	double skeff_wght = fSKDet->EfficiencyWeight(s,m,ErecoVarMin);

	// Plug-in FSI weights
	double fsi_wght = fFSI->SystWght(s,m,ErecoVarMin,EtrueVarMin); 

	// Combined weight
	wght = skeff_wght * flux_wght * xsec_wght * fsi_wght;
#endif//!SKFSITOGETHER

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
	  
	  
	  if (s==kSKsFCFV1RingMuLike){
#if BANFF_BINS == 4
	    if(!util::IsMEC(m))
#endif
	      assert(fMCTemplArrayOri1Rmu[im][it-1][ir-1] > 0);
	  }
	  // MODIFIED
	  else if (s==kSKsFCFV1RingELike){
	    assert(fMCTemplArrayOri1Re[im][it-1][ir-1] > 0);
	  }
	  ////
	  // ANTINU
	  else if (s==kSKsRHCFCFV1RingMuLike){
#if BANFF_BINS == 4
	    if(!util::IsMEC(m))
#endif
	      assert(fMCTemplArrayOriRHC1Rmu[im][it-1][ir-1] > 0);
	  }
	  else if (s==kSKsRHCFCFV1RingELike){
	    assert(fMCTemplArrayOriRHC1Re[im][it-1][ir-1] > 0);
	  }
	  ////





	  //double Etrue = util::ErecoBinEdges(s)(s,it-1) + kOneMeV;

	  // get number of events in current bin of source MC template
	  double nev = 0;
#if BANFF_BINS == 4
	  if(util::IsMEC(m)) {

	  
	    // MODIFIED
	    double nev_martini;
	    double nev_nieves;
	    if (s==kSKsFCFV1RingMuLike || s==kSKsFCFV1RingELike){
	      nev_martini = fMCTemplArrayMECMartini[it-1][ir-1];
	      nev_nieves  = fMCTemplArrayMECNieves [it-1][ir-1];
	    }
	    ////
	    // ANTINU
	    else if (s==kSKsRHCFCFV1RingMuLike || s==kSKsRHCFCFV1RingELike){
	      nev_martini = fMCTemplArrayRHCMECMartini[it-1][ir-1];
	      nev_nieves  = fMCTemplArrayRHCMECNieves [it-1][ir-1];
	    }
	    ////

	    
	    //-1 nominal -> 0 nieves
	    if(x_mec < 0) {
	      nev  = nev_nieves * (1 + x_mec);
	    }
	    // 0 nieves -> +1 martini
	    else {
	      nev += TMath::Min(nev_martini, nev_martini * x_mec);
	      nev += nev_nieves * (1 - x_mec);
	    }
	  }//MEC
	  else
#endif


	    if (s==kSKsFCFV1RingMuLike) nev = fMCTemplArrayTun1Rmu[im][it-1][ir-1];
	    else if (s==kSKsFCFV1RingELike) nev = fMCTemplArrayTun1Re[im][it-1][ir-1]; // MODIFIED
	    else if (s==kSKsRHCFCFV1RingMuLike) nev = fMCTemplArrayTunRHC1Rmu[im][it-1][ir-1]; // ANTINU
	    else if (s==kSKsRHCFCFV1RingELike) nev = fMCTemplArrayTunRHC1Re[im][it-1][ir-1]; // ANTINU
	  
	  
	  
	  
	  nev = TMath::Max(0.,nev);
	  double nev_systm = nev * wght;

#ifdef __show_syst_var_wght_mesg_debug__ 
	    // Get the energy of the current bin
	    double Etrue = util::EtrueBinEdges(s,it-1) + kOneMeV;
	    double Ereco = util::ErecoBinEdges(s,it-1) + kOneMeV;
#ifndef SKFSITOGETHER	       
	  if( skeff_wght != 1. || 
	      flux_wght  != 1. || 
	      xsec_wght  != 1. ||
	      fsi_wght   != 1. )
	  {
	    cout << "Exp/SKMC/SystVarWght >> "
	      << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m) 
	      << ", Ereco = " << Ereco << " GeV" 
	      << ", Etrue = " << Etrue << " GeV => Weights: " 
	      << "SK eff = "  << skeff_wght
	      << ", Flux = "  << flux_wght
	      << ", XSec = "  << xsec_wght
	      << ", FSI  = "  << fsi_wght
	      << " (total = " << wght << ")"
	      << endl;
	  }
#else
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
#endif//SKFSITOGETHER
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
#ifndef PRODUCTIONMODE
	  mctempl_tgt->SetBinContent(ir,it,nev_systm);
#endif

	  
	  if (s==kSKsFCFV1RingMuLike) fMCTemplArrayFin1Rmu[im][it-1][ir-1] = nev_systm;
	  else if (s==kSKsFCFV1RingELike) fMCTemplArrayFin1Re[im][it-1][ir-1] = nev_systm;  // MODIFIED
	  else if (s==kSKsRHCFCFV1RingMuLike) fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1] = nev_systm; // ANTINU
	  else if (s==kSKsRHCFCFV1RingELike) fMCTemplArrayFinRHC1Re[im][it-1][ir-1] = nev_systm;  // ANTINU


	}//ifb
      }//ivb
    }//im
  }//is

#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtSys, "mcnorm_fluxtuned_banfftuned_syst_multi");
#endif
}
//________________________________________________________________________________________
void SuperkMC::ApplyMECRenormalisationToTemplates(void) 
{
// For use when the BANFF doesn't give us an updated MEC matrix
// MEC is just an addition of events, so we can renormalise non-MEC events
// to be more realistic with it's effect
//


  // MODIFIED
  // Declare variables used both for FHC and RHC modes. Avoid compilation errors when both RHC and FHC modes are used
  SuperkSample_t s=kSKsUndefined;
  int n_ereco_bins=0;
  int n_etrue_bins=0;
  double integrals[knm] = {0};
  double renorm_total = 0;
  double renorm_mec   = 0;
  double renorm_templ = 0;
  double renorm_factor = 0.;
  ////


#ifdef __show_update_mesg_info__
  cout << "Exp/SKMC >> Applying MEC renormalisation" << endl;
#endif

#ifdef __FHC_RUNNING__          
  
  //MEC only included in 1Rmu sample
  s = kSKsFCFV1RingMuLike;
  n_ereco_bins = kNErecoBins[0];
  n_etrue_bins = kNEtrueBins[0];
  
  //define which modes we're interested in normalising to.
  // to turn off, swap a mode with kSKmUndefined
  const SuperkMode_t RenormModes1Rmu[knm] = {
    kSKmNuMuCCQE,         
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    kSKmNuMuCCoth,                      
    kSKmUndefined, //kSKmNuMuNC1piC,
    kSKmUndefined, //kSKmNuMuNCoth_mu,          
    kSKmOscNuMuCCQE,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCcoh,
    kSKmOscNuMuCCoth,
    kSKmNuMuBarCCQE,      
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCcoh,
    kSKmNuMuBarCCoth,                  
    kSKmUndefined, //kSKmNuMuBarNC1piC,
    kSKmUndefined, //kSKmNuMuBarNCoth_mu,
    kSKmNuECCQE, 
    kSKmNuECC1pi,
    kSKmNuECCcoh,
    kSKmNuECCoth,     
                           
    kSKmUndefined, //kSKmNuENC1piC,
    kSKmUndefined, //kSKmNuENCoth_mu,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCCoth,
    kSKmUndefined, //kSKmNuEBarNC1piC,
    kSKmUndefined, //kSKmNuEBarNCoth_mu,
    kSKmOscNuECCQE,
    kSKmOscNuECC1pi,
    kSKmOscNuECCcoh,
    kSKmOscNuECCoth,
    kSKmUndefined,//kSKmMEC, //keep MEC templates here - they won't be renormalised!
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined
  };
  
  //FIRST integrate the templates
  // Loop over Renorm modes
  for(int im = 0; im < knm; im++) { 
    SuperkMode_t m = kSKReactionModes[0][im];//RenormModes1Rmu[im];
    if(m == kSKmUndefined) continue;
    // Consider template?
    bool skip = this->SkipTemplate(m);
    if(skip) continue;
    //initialise
    integrals[im] = 0;
    //loop over Ereco
    for(int ir = 1; ir <= n_ereco_bins; ir++) {
      // loop over Etrue
      for(int it = 1; it <= n_etrue_bins; it++) {
	double nev = fMCTemplArrayFin1Rmu[im][it-1][ir-1];
	nev = TMath::Max(0.,nev);
	integrals[im] += nev;
      }//it
    }//ir
    cout << integrals[im] << "\t" << util::AsStringCompact(m) << endl;
  }//im

  //NEXT get total number of events
  renorm_total = 0;
  renorm_mec   = 0;
  renorm_templ = 0;
  // Loop over Renorm modes
  for(int im = 0; im < knm; im++) { 
    SuperkMode_t m = kSKReactionModes[0][im];//RenormModes1Rmu[im];
    if(m == kSKmUndefined) continue;
    // Consider template?
    bool skip = this->SkipTemplate(m);
    if(skip) continue;
    if(util::IsMEC(m))
      renorm_mec   += integrals[im];
    else
      renorm_total += integrals[im];
  }
  for(int im = 0; im < knm; im++) { 
    SuperkMode_t m = RenormModes1Rmu[im];
    if(m == kSKmUndefined) continue;
    // Consider template?
    bool skip = this->SkipTemplate(m);
    if(skip) continue;
    renorm_templ += integrals[im];
  }
  //and calculate the renorm factor
  renorm_factor = (renorm_templ - renorm_mec)/renorm_templ;  
  cout << "FHC mode" << endl;    // MODIFIED
  cout << "Exp/SKMC >> " << renorm_mec   << " MEC events"           << endl;
  cout << "Exp/SKMC >> " << renorm_total << " total non-MEC events" << endl;
  cout << "Exp/SKMC >> " << renorm_templ << " events in templates we're renormalising" << endl;
  cout << "Exp/SKMC >> " << renorm_total - renorm_templ << " events in templates we're not renormalising" << endl;
  cout << "Exp/SKMC >> MEC renormalisation factor is " << renorm_factor << endl;
  
  //FINALLY apply the factor
  // Loop over Renormsample modes
  for(int im = 0; im < knm; im++) { 
    SuperkMode_t m = RenormModes1Rmu[im];
    if(m == kSKmUndefined) continue;
    // Consider template?
    bool skip = this->SkipTemplate(m);
    if(skip) continue;
    if(util::IsMEC(m)) continue;
    //loop over Ereco
    for(int ir = 1; ir <= n_ereco_bins; ir++) {
      // loop over Etrue
      for(int it = 1; it <= n_etrue_bins; it++) {
	fMCTemplArrayFin1Rmu[im][it-1][ir-1] *= renorm_factor;
#ifndef PRODUCTIONMODE
	TH2D * mctempl_tgt = this -> Template(s,m,kSKMCtSys);
	if(mctempl_tgt)
	  mctempl_tgt->SetBinContent(ir,it,fMCTemplArrayFin1Rmu[im][it-1][ir-1]);
#endif
      }//it
    }//ir
  }//im

#endif //__FHC_RUNNING__
    
#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtFin, "mcnorm_fluxtuned_banfftuned_syst_multi_mecrenorm");
#endif
  
  
  // ANTINU
#ifdef __RHC_RUNNING__          
  
  //MEC only included in 1Rmu sample
  s = kSKsRHCFCFV1RingMuLike;
  n_ereco_bins = kNErecoBins[2];
  n_etrue_bins = kNEtrueBins[2];
  
  //define which modes we're interested in normalising to.
  // to turn off, swap a mode with kSKmUndefined
  const SuperkMode_t RenormModesRHC1Rmu[knm] = {
    kSKmNuMuCCQE,         
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    kSKmNuMuCCoth,                      
    kSKmUndefined, //kSKmNuMuNC1piC,
    kSKmUndefined, //kSKmNuMuNCoth_mu,          
    kSKmOscNuMuCCQE,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCcoh,
    kSKmOscNuMuCCoth,
    kSKmNuMuBarCCQE,      
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCcoh,
    kSKmNuMuBarCCoth,                  
    kSKmUndefined, //kSKmNuMuBarNC1piC,
    kSKmUndefined, //kSKmNuMuBarNCoth_mu,
    kSKmNuECCQE, 
    kSKmNuECC1pi,
    kSKmNuECCcoh,
    kSKmNuECCoth,     
                           
    kSKmUndefined, //kSKmNuENC1piC,
    kSKmUndefined, //kSKmNuENCoth_mu,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCCoth,
    kSKmUndefined, //kSKmNuEBarNC1piC,
    kSKmUndefined, //kSKmNuEBarNCoth_mu,
    kSKmOscNuEBarCCQE,   // ANTINU --> OscNuEBar instead of OscNuE
    kSKmOscNuEBarCC1pi,  // ANTINU --> OscNuEBar instead of OscNuE
    kSKmOscNuEBarCCcoh,  // ANTINU --> OscNuEBar instead of OscNuE
    kSKmOscNuEBarCCoth,  // ANTINU --> OscNuEBar instead of OscNuE
    kSKmUndefined,//kSKmMEC, //keep MEC templates here - they won't be renormalised!
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined
  };
  
  //FIRST integrate the templates
  // Loop over Renorm modes
  for(int im = 0; im < knm; im++) { 
    SuperkMode_t m = kSKReactionModes[2][im];//RenormModesRHC1Rmu[im];  // is=2 --> RHC numu disapp // ANTINU
    if(m == kSKmUndefined) continue;
    // Consider template?
    bool skip = this->SkipTemplate(m);
    if(skip) continue;
    //initialise
    integrals[im] = 0;
    //loop over Ereco
    for(int ir = 1; ir <= n_ereco_bins; ir++) {
      // loop over Etrue
      for(int it = 1; it <= n_etrue_bins; it++) {
	double nev = fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1];
	nev = TMath::Max(0.,nev);
	integrals[im] += nev;
      }//it
    }//ir
    cout << integrals[im] << "\t" << util::AsStringCompact(m) << endl;
  }//im

  //NEXT get total number of events
  renorm_total = 0;
  renorm_mec   = 0;
  renorm_templ = 0;
  // Loop over Renorm modes
  for(int im = 0; im < knm; im++) { 
    SuperkMode_t m = kSKReactionModes[2][im];//RenormModesRHC1Rmu[im]; // is=2 --> RHC numu disapp // ANTINU
    if(m == kSKmUndefined) continue;
    // Consider template?
    bool skip = this->SkipTemplate(m);
    if(skip) continue;
    if(util::IsMEC(m))
      renorm_mec   += integrals[im];
    else
      renorm_total += integrals[im];
  }
  for(int im = 0; im < knm; im++) { 
    SuperkMode_t m = RenormModesRHC1Rmu[im];
    if(m == kSKmUndefined) continue;
    // Consider template?
    bool skip = this->SkipTemplate(m);
    if(skip) continue;
    renorm_templ += integrals[im];
  }
  //and calculate the renorm factor
  renorm_factor = (renorm_templ - renorm_mec)/renorm_templ;  
  cout << "RHC mode" << endl;    // MODIFIED
  cout << "Exp/SKMC >> " << renorm_mec   << " MEC events"           << endl;
  cout << "Exp/SKMC >> " << renorm_total << " total non-MEC events" << endl;
  cout << "Exp/SKMC >> " << renorm_templ << " events in templates we're renormalising" << endl;
  cout << "Exp/SKMC >> " << renorm_total - renorm_templ << " events in templates we're not renormalising" << endl;
  cout << "Exp/SKMC >> MEC renormalisation factor is " << renorm_factor << endl;
  
  //FINALLY apply the factor
  // Loop over Renormsample modes
  for(int im = 0; im < knm; im++) { 
    SuperkMode_t m = RenormModesRHC1Rmu[im];
    if(m == kSKmUndefined) continue;
    // Consider template?
    bool skip = this->SkipTemplate(m);
    if(skip) continue;
    if(util::IsMEC(m)) continue;
    //loop over Ereco
    for(int ir = 1; ir <= n_ereco_bins; ir++) {
      // loop over Etrue
      for(int it = 1; it <= n_etrue_bins; it++) {
	fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1] *= renorm_factor;
#ifndef PRODUCTIONMODE
	TH2D * mctempl_tgt = this -> Template(s,m,kSKMCtSys);
	if(mctempl_tgt)
	  mctempl_tgt->SetBinContent(ir,it,fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]);
#endif
      }//it
    }//ir
  }//im

#endif //__RHC_RUNNING__
    
#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtFin, "mcnorm_fluxtuned_banfftuned_syst_multi_mecrenorm");
#endif
 ////
  



}
//________________________________________________________________________________________
void SuperkMC::ApplySystematicsToTemplates_NonMultiplicative(void) 
{
//
// Now that all simple multiplicative systematics were applied, apply the reconstructed 
// energy scale systematic. This systematic is applied by scaling the Ereco bin edges
// and, by assuming uniform distribution of events within each bin, calculating how many
// events a lost to / gained from neighbouring bins.
//
  
#ifdef __show_update_mesg_info__
  cout << "Exp/SKMC >> Applying non-multiplicative systematics" << endl;
#endif
  
  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;

    int n_ereco_bins = kNErecoBins[is];
    int n_etrue_bins = kNEtrueBins[is];
    
    // loop over sample modes
    for(int im = 0; im < knm; im++) { 
      SuperkMode_t m = kSKReactionModes[is][im];
      if(m == kSKmUndefined) continue;
      // consider template?
      bool skip = this->SkipTemplate(m);
      if(skip) continue;
#ifndef PRODUCTIONMODE
      // get corresponding template
      TH2D * mctempl = this -> Template(s,m,kSKMCtSys);
      if(!mctempl) continue; 
      // clone
      TH2D * mctempl_no_ereco_scale = new TH2D(*mctempl);
#endif
      double mctemplarry_no_ereco_scale[kns][knm][n_etrue_bins][n_ereco_bins];

      if (s==kSKsFCFV1RingMuLike)
	memcpy (mctemplarry_no_ereco_scale, fMCTemplArrayFin1Rmu[im], sizeof (fMCTemplArrayFin1Rmu[im])) ;
      // MODIFIED
      else if (s==kSKsFCFV1RingELike)
	memcpy (mctemplarry_no_ereco_scale, fMCTemplArrayFin1Re[im], sizeof (fMCTemplArrayFin1Re[im]));
      ////
      // ANTINU
      else if (s==kSKsRHCFCFV1RingMuLike)
	memcpy (mctemplarry_no_ereco_scale, fMCTemplArrayFinRHC1Rmu[im], sizeof (fMCTemplArrayFinRHC1Rmu[im])) ;
      else if (s==kSKsRHCFCFV1RingELike)
	memcpy (mctemplarry_no_ereco_scale, fMCTemplArrayFinRHC1Re[im], sizeof (fMCTemplArrayFinRHC1Re[im]));
      ////





      // get energy scale
      double escale = 1;
#ifdef SKFSITOGETHER
      escale = fSKDetFSI->EnergyScale(s,m);
#endif 
#ifndef SKFSITOGETHER
      escale = fSKDet->EnergyScale(s,m);
#endif
      // energy scale > 1 : 
      // For a given bin, events are gained from bins on the left (lower Ereco bins) and 
      // are lost to bins on the right (higher Ereco bins)
      if(escale > 1. + kASmallNum) {
	// loop over Ereco
	 for(int ir = 1; ir <= n_ereco_bins; ir++) {
	   // nominal reco bin edges
	   double bin_max = util::ErecoBinEdges(s,ir-1);
	   double bin_min = util::ErecoBinEdges(s,ir);
	   // scaled reco bin edges
              bin_max *= escale;
              bin_min *= escale;
              // loop over Ereco' (> Ereco)
              for(int irp = ir+1; irp <= n_ereco_bins; irp++) {
                 // bin edges for reco bin ir'
		double bin_max_p = util::ErecoBinEdges(s,irp-1);
		double bin_min_p = util::ErecoBinEdges(s,irp);
                 // do bins ir and ir' overlap?
                 if(bin_max <= bin_min_p || bin_min >= bin_max_p) continue; // no overlap
                 // calculate overlap
                 double overlap = 0;
                 if(bin_min_p < bin_max && bin_min_p > bin_min) overlap = bin_max   - bin_min_p;
                 if(bin_max_p < bin_max && bin_max_p > bin_min) overlap = bin_max_p - bin_min;
                 if(bin_max_p > bin_max && bin_min_p < bin_min) overlap = bin_max   - bin_min;
                 if(bin_max_p < bin_max && bin_min_p > bin_min) overlap = bin_max_p - bin_min_p;                 
                 // loop over Etrue
                 for(int it = 1; it <= n_etrue_bins; it++) {
                      // calculate how many events from nominal bin ir get lost to bin ir'
 		      double n  = mctemplarry_no_ereco_scale[is][im][it-1][ir-1];
                      double dn = n * overlap / (bin_max-bin_min); ///< (overlap) * (nominal # of events) / (bin width)
                      // update the systematically varied templates
                      dn = TMath::Min(dn,n);
#ifndef PRODUCTIONMODE





		      if (s==kSKsFCFV1RingMuLike)  mctempl->SetBinContent(ir, it, fMCTemplArrayFin1Rmu[im][it-1][ir-1]  - dn);
		      // MODIFIED
		      else if (s==kSKsFCFV1RingELike) mctempl->SetBinContent(ir, it, fMCTemplArrayFin1Re[im][it-1][ir-1]  - dn);
		      ////
		      // ANTINU
		      else if (s==kSKsRHCFCFV1RingMuLike)  mctempl->SetBinContent(ir, it, fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]  - dn);
		      else if (s==kSKsRHCFCFV1RingELike) mctempl->SetBinContent(ir, it, fMCTemplArrayFinRHC1Re[im][it-1][ir-1]  - dn);
		      ////

                      if (s==kSKsFCFV1RingMuLike)  mctempl->SetBinContent(irp,it, fMCTemplArrayFin1Rmu[im][it-1][irp-1] + dn);
		      // MODIFIED
		      else if (s==kSKsFCFV1RingELike) mctempl->SetBinContent(irp,it, fMCTemplArrayFin1Re[im][it-1][irp-1] + dn);
		      ////
		      // ANTINU
                      else if (s==kSKsRHCFCFV1RingMuLike)  mctempl->SetBinContent(irp,it, fMCTemplArrayFinRHC1Rmu[im][it-1][irp-1] + dn);
		      else if (s==kSKsRHCFCFV1RingELike) mctempl->SetBinContent(irp,it, fMCTemplArrayFinRHC1Re[im][it-1][irp-1] + dn);
		      ////





#endif



		      if (s==kSKsFCFV1RingMuLike) fMCTemplArrayFin1Rmu[im][it-1][ir-1]  -= dn;
		      // MODIFIED
		      else if (s==kSKsFCFV1RingELike) fMCTemplArrayFin1Re[im][it-1][ir-1]  -= dn;
		      ////
		      // ANTINU
		      else if (s==kSKsRHCFCFV1RingMuLike) fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]  -= dn;
		      else if (s==kSKsRHCFCFV1RingELike) fMCTemplArrayFinRHC1Re[im][it-1][ir-1]  -= dn;
		      ////

		      
		      if (s==kSKsFCFV1RingMuLike) fMCTemplArrayFin1Rmu[im][it-1][irp-1] += dn;
		      // MODIFIED
		      else if (s==kSKsFCFV1RingELike) fMCTemplArrayFin1Re[im][it-1][irp-1] += dn;
		      ////
		      // ANTINU
		      else if (s==kSKsRHCFCFV1RingMuLike) fMCTemplArrayFinRHC1Rmu[im][it-1][irp-1] += dn;
		      else if (s==kSKsRHCFCFV1RingELike) fMCTemplArrayFinRHC1Re[im][it-1][irp-1] += dn;
		      ////




                 }//it
              }//irp
           }//ir
       }// scale > 1

      // energy scale < 1 : 
      // For a given bin, events are gained from bins on the right (higher Ereco bins) and 
      // are lost to bins on the left (lower Ereco bins)
      if(escale < 1. - kASmallNum) {
	// loop over Ereco
	for(int ir = 1; ir <= n_ereco_bins; ir++) {
	  // nominal reco bin edges
	  double bin_max = util::ErecoBinEdges(s,ir-1);
	  double bin_min = util::ErecoBinEdges(s,ir);
	  // scaled reco bin edges
	  bin_max *= escale;
	  bin_min *= escale;
	  // loop over Ereco' (< Ereco)
	  for(int irp = 1; irp <= ir-1; irp++) {
	    // bin edges for reco bin ir'
	    double bin_max_p = util::ErecoBinEdges(s,irp-1);
	    double bin_min_p = util::ErecoBinEdges(s,irp);
	    // do bins ir and ir' overlap?
	    if(bin_max <= bin_min_p || bin_min >= bin_max_p) continue; // no overlap
	    // calculate overlap
	    double overlap = 0;
	    if(bin_min_p < bin_max && bin_min_p > bin_min) overlap = bin_max   - bin_min_p;
	    if(bin_max_p < bin_max && bin_max_p > bin_min) overlap = bin_max_p - bin_min;
	    if(bin_max_p > bin_max && bin_min_p < bin_min) overlap = bin_max   - bin_min;
	    if(bin_max_p < bin_max && bin_min_p > bin_min) overlap = bin_max_p - bin_min_p;
	    // loop over Etrue
	    for(int it = 1; it <= n_etrue_bins; it++) {
	      // calculate how many events from nominal bin ir get lost to bin ir'
	      double n  = mctemplarry_no_ereco_scale[is][im][it-1][ir-1];
	      double dn = n * overlap / (bin_max-bin_min); ///< (overlap) * (nominal # of events) / (bin width)
	      // update the systematically varied templates
	      dn = TMath::Min(dn,n);

#ifndef PRODUCTIONMODE


	      if (s==kSKsFCFV1RingMuLike) mctempl->SetBinContent(ir, it, fMCTemplArrayFin1Rmu[im][it-1][ir-1]  - dn);
	      // MODIFIED
	      else if (s==kSKsFCFV1RingELike) mctempl->SetBinContent(ir, it, fMCTemplArrayFin1Re[im][it-1][ir-1]  - dn);
	      ////
	      // ANTINU
	      else if (s==kSKsRHCFCFV1RingMuLike) mctempl->SetBinContent(ir, it, fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]  - dn);
	      else if (s==kSKsRHCFCFV1RingELike) mctempl->SetBinContent(ir, it, fMCTemplArrayFinRHC1Re[im][it-1][ir-1]  - dn);
	      ////


	      if (s==kSKsFCFV1RingMuLike) mctempl->SetBinContent(irp,it, fMCTemplArrayFin1Rmu[im][it-1][irp-1] + dn);
	      // MODIFIED
	      else if (s==kSKsFCFV1RingELike) mctempl->SetBinContent(irp,it, fMCTemplArrayFin1Re[im][it-1][irp-1] + dn);
	      ////
	      // ANTINU
	      else if (s==kSKsRHCFCFV1RingMuLike) mctempl->SetBinContent(irp,it, fMCTemplArrayFinRHC1Rmu[im][it-1][irp-1] + dn);
	      else if (s==kSKsRHCFCFV1RingELike) mctempl->SetBinContent(irp,it, fMCTemplArrayFinRHC1Re[im][it-1][irp-1] + dn);
	      ////
	      

#endif

	      if (s==kSKsFCFV1RingMuLike) fMCTemplArrayFin1Rmu[im][it-1][ir-1]  -= dn;
	      // MODIFIED
	      else if (s==kSKsFCFV1RingELike) fMCTemplArrayFin1Re[im][it-1][ir-1]  -= dn; 
	      ////
	      // ANTINU
	      else if (s==kSKsRHCFCFV1RingMuLike) fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]  -= dn;
	      else if (s==kSKsRHCFCFV1RingELike) fMCTemplArrayFinRHC1Re[im][it-1][ir-1]  -= dn; 
	      ////

	      
	      if (s==kSKsFCFV1RingMuLike) fMCTemplArrayFin1Rmu[im][it-1][irp-1] += dn;
	      // MODIFIED
	      else if (s==kSKsFCFV1RingELike) fMCTemplArrayFin1Re[im][it-1][irp-1] += dn;
	      ////
	      // ANTINU 
	      else if (s==kSKsRHCFCFV1RingMuLike) fMCTemplArrayFinRHC1Rmu[im][it-1][irp-1] += dn;
	      else if (s==kSKsRHCFCFV1RingELike) fMCTemplArrayFinRHC1Re[im][it-1][irp-1] += dn;
	      ////


	    }//it
	  }//irp
	}//ir
      }// scale < 1
#ifndef PRODUCTIONMODE
      delete mctempl_no_ereco_scale;
#endif
      
#ifdef NO_EVENTS_CAN_BE_NEGATIVE
       // loop at the end, to get rid of any negatives
       // loop over Ereco
      for(int ir = 1; ir <= n_ereco_bins; ir++) {
	// loop over Etrue
	for(int it = 1; it <= n_etrue_bins; it++) {

	  
	  if (s==kSKsFCFV1RingMuLike){
	     if(fMCTemplArrayFin1Rmu[im][it-1][ir-1] < 0) {
#ifdef __show_negative_events_mesg_debug__
	     cout << "Exp/SKMC/ApplySystematicsToTemplates_NonMultiplicative " << fMCTemplArrayFin1Rmu[im][it-1][ir-1]
		  << " to " << 0 << "\t" << " Ereco: " << kErecoBinEdges1Rmu[ir-1]+kOneMeV
		  << " Etrue: " << kEtrueBinEdges1Rmu[it-1]+kOneMeV << endl;
#endif
#ifndef PRODUCTIONMODE
	     mctempl->SetBinContent(ir,  it,  0);
#endif
	     fMCTemplArrayFin1Rmu[im][it-1][ir-1] = 0;
	     }
	   }
	  else if (s==kSKsFCFV1RingELike) {                  // MODIFIED
	     if(fMCTemplArrayFin1Re[im][it-1][ir-1] < 0) {
#ifdef __show_negative_events_mesg_debug__
	     cout << "Exp/SKMC/ApplySystematicsToTemplates_NonMultiplicative " << fMCTemplArrayFin1Re[im][it-1][ir-1]
		  << " to " << 0 << "\t" << " Ereco: " << kErecoBinEdges1Re[ir-1]+kOneMeV
		  << " Etrue: " << kEtrueBinEdges1Re[it-1]+kOneMeV << endl;
#endif
#ifndef PRODUCTIONMODE
	     mctempl->SetBinContent(ir,  it,  0);
#endif
	     fMCTemplArrayFin1Re[im][it-1][ir-1] = 0;
	     }
	  }



	  // ANTINU
	  else if (s==kSKsRHCFCFV1RingMuLike){
	    if(fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1] < 0) {
#ifdef __show_negative_events_mesg_debug__
	      cout << "Exp/SKMC/ApplySystematicsToTemplates_NonMultiplicative " << fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]
		   << " to " << 0 << "\t" << " Ereco: " << kRHCErecoBinEdges1Rmu[ir-1]+kOneMeV
		   << " Etrue: " << kRHCEtrueBinEdges1Rmu[it-1]+kOneMeV << endl;
#endif
#ifndef PRODUCTIONMODE
	      mctempl->SetBinContent(ir,  it,  0);
#endif
	      fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1] = 0;
	    }
	  }
	  else if (s==kSKsRHCFCFV1RingELike) {                 
	    if(fMCTemplArrayFinRHC1Re[im][it-1][ir-1] < 0) {
#ifdef __show_negative_events_mesg_debug__
	      cout << "Exp/SKMC/ApplySystematicsToTemplates_NonMultiplicative " << fMCTemplArrayFinRHC1Re[im][it-1][ir-1]
		   << " to " << 0 << "\t" << " Ereco: " << kRHCErecoBinEdges1Re[ir-1]+kOneMeV
		   << " Etrue: " << kRHCEtrueBinEdges1Re[it-1]+kOneMeV << endl;
#endif
#ifndef PRODUCTIONMODE
	      mctempl->SetBinContent(ir,  it,  0);
#endif
	      fMCTemplArrayFinRHC1Re[im][it-1][ir-1] = 0;
	    }
	  }
	  ////
	  




	}
       }
#endif


    }//im
  }//is
#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtSys, "mcnorm_fluxtuned_banfftuned_syst_nonmulti");
#endif
}
//________________________________________________________________________________________
bool SuperkMC::ApplyEnergyScaleSystematicToTotalSpectrum(void) 
{
//
// Now that all other systematics and physics hypotheses (that depend on Etrue) have been 
// applied, apply the reconstructed energy scale systematic (that depends on Ereco).
// This systematic is applied by scaling the Ereco bin edges and, by assuming uniform 
//distribution of events within each bin, calculating how many events a lost to / gained
// from neighbouring bins.
//

#ifdef __show_update_mesg_info__
  cout << "Exp/SKMC >> Applying non-multiplicative systematics to total spectrum" << endl;
#endif

  bool applied_negevents_condition = false;
  
  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;

    int n_ereco_bins = kNErecoBins[is];

    //get corresponding spectrum
#ifndef PRODUCTIONMODE
    TH1D * reco_spectrum_tot = this->CurrErecoSpectrum(s);
    if(!reco_spectrum_tot) continue;
#endif
    //clone
    //TH1D * reco_spectrum_tot_no_ereco_scale = new TH1D(*reco_spectrum_tot);
    double reco_spectrum_tot_arry_no_ereco_scale[kns][n_ereco_bins];

    if (s==kSKsFCFV1RingMuLike) 
      memcpy(reco_spectrum_tot_arry_no_ereco_scale[is], fMCErecoSpectrumTotArray1Rmu, sizeof (fMCErecoSpectrumTotArray1Rmu)) ;
    // MODIFIED
    else if (s==kSKsFCFV1RingELike)
      memcpy(reco_spectrum_tot_arry_no_ereco_scale[is], fMCErecoSpectrumTotArray1Re, sizeof (fMCErecoSpectrumTotArray1Re)) ;
    ////
    // ANTINU
    else if (s==kSKsRHCFCFV1RingMuLike) 
      memcpy(reco_spectrum_tot_arry_no_ereco_scale[is], fMCErecoSpectrumTotArrayRHC1Rmu, sizeof (fMCErecoSpectrumTotArrayRHC1Rmu)) ;
    else if (s==kSKsRHCFCFV1RingELike)
      memcpy(reco_spectrum_tot_arry_no_ereco_scale[is], fMCErecoSpectrumTotArrayRHC1Re, sizeof (fMCErecoSpectrumTotArrayRHC1Re)) ;
    ////

     // get energy scale
     double escale = 1;
 #ifdef SKFSITOGETHER
     escale = fSKDetFSI->EnergyScale(s);
 #endif
 #ifndef SKFSITOGETHER
     escale = fSKDet->EnergyScale(s);
 #endif
     // energy scale > 1 : 
     // For a given bin, events are gained from bins on the left (lower Ereco bins) and 
     // are lost to bins on the right (higher Ereco bins)
     if(escale > 1. + kASmallNum) {
       // loop over Ereco
       for(int ir = 1; ir <= n_ereco_bins; ir++) {
	 // nominal reco bin edges
	 //double bin_max = reco_spectrum_tot->GetXaxis()->GetBinUpEdge (ir);
	 //double bin_min = reco_spectrum_tot->GetXaxis()->GetBinLowEdge(ir);
	 double bin_max = util::ErecoBinEdges(s, ir);
	 double bin_min = util::ErecoBinEdges(s, ir-1);
	 // scaled reco bin edges
	 bin_max *= escale;
	 bin_min *= escale;
	 // loop over Ereco' (> Ereco)
	 for(int irp = ir+1; irp <= n_ereco_bins; irp++) {
	   // bin edges for reco bin ir'
	   //double bin_max_p = reco_spectrum_tot->GetXaxis()->GetBinUpEdge (irp);
	   //double bin_min_p = reco_spectrum_tot->GetXaxis()->GetBinLowEdge(irp);
	   double bin_max_p = util::ErecoBinEdges(s, irp);
	   double bin_min_p = util::ErecoBinEdges(s, irp-1);
	   // do bins ir and ir' overlap?
	   if(bin_max <= bin_min_p || bin_min >= bin_max_p) continue; // no overlap
	   // calculate overlap
	   double overlap = 0;
	   if(bin_min_p < bin_max && bin_min_p > bin_min) overlap = bin_max   - bin_min_p;
	   if(bin_max_p < bin_max && bin_max_p > bin_min) overlap = bin_max_p - bin_min;
	   if(bin_max_p > bin_max && bin_min_p < bin_min) overlap = bin_max   - bin_min;
	   if(bin_max_p < bin_max && bin_min_p > bin_min) overlap = bin_max_p - bin_min_p;

	   // calculate how many events from nominal bin ir get lost to bin ir'
	   double n  = reco_spectrum_tot_arry_no_ereco_scale[is][ir-1];
	   double dn = n * overlap / (bin_max-bin_min); ///< (overlap) * (nominal # of events) / (bin width)
	   // update the systematically varied spectrum
	   dn = TMath::Min(dn,n);
 #ifndef PRODUCTIONMODE

	   

	   if (s==kSKsFCFV1RingMuLike) reco_spectrum_tot->SetBinContent(ir,  fMCErecoSpectrumTotArray1Rmu[ir-1]  - dn);
	   // MODIFIED
	   else if (s==kSKsFCFV1RingELike) reco_spectrum_tot->SetBinContent(ir,  fMCErecoSpectrumTotArray1Re[ir-1]  - dn);
	   ////
	   // ANTINU
	   else if (s==kSKsRHCFCFV1RingMuLike) reco_spectrum_tot->SetBinContent(ir,  fMCErecoSpectrumTotArrayRHC1Rmu[ir-1]  - dn);
	   else if (s==kSKsRHCFCFV1RingELike) reco_spectrum_tot->SetBinContent(ir,  fMCErecoSpectrumTotArrayRHC1Re[ir-1]  - dn);
	   ////

	   
	   if (s==kSKsFCFV1RingMuLike) reco_spectrum_tot->SetBinContent(irp, fMCErecoSpectrumTotArray1Rmu[irp-1] + dn);
	   // MODIFIED
	   else if (s==kSKsFCFV1RingELike) reco_spectrum_tot->SetBinContent(irp, fMCErecoSpectrumTotArray1Re[irp-1] + dn);
	   ////
	   // ANTINU
	   else if (s==kSKsRHCFCFV1RingMuLike) reco_spectrum_tot->SetBinContent(irp, fMCErecoSpectrumTotArrayRHC1Rmu[irp-1] + dn);
	   else if (s==kSKsRHCFCFV1RingELike) reco_spectrum_tot->SetBinContent(irp, fMCErecoSpectrumTotArrayRHC1Re[irp-1] + dn);
	   ////



#endif


	   if (s==kSKsFCFV1RingMuLike) fMCErecoSpectrumTotArray1Rmu[ir-1]  -= dn;
	   // MODIFIED
	   else if (s==kSKsFCFV1RingELike) fMCErecoSpectrumTotArray1Re[ir-1]  -= dn;
	   ////
	   // ANTINU
	   else if (s==kSKsRHCFCFV1RingMuLike) fMCErecoSpectrumTotArrayRHC1Rmu[ir-1]  -= dn;
	   else if (s==kSKsRHCFCFV1RingELike) fMCErecoSpectrumTotArrayRHC1Re[ir-1]  -= dn;
	   ////


	   if (s==kSKsFCFV1RingMuLike) fMCErecoSpectrumTotArray1Rmu[irp-1] += dn;
	   // MODIFIED
	   else if (s==kSKsFCFV1RingELike) fMCErecoSpectrumTotArray1Re[irp-1] += dn;
	   ////
	   // ANTINU
	   else if (s==kSKsRHCFCFV1RingMuLike) fMCErecoSpectrumTotArrayRHC1Rmu[irp-1] += dn;
	   else if (s==kSKsRHCFCFV1RingELike) fMCErecoSpectrumTotArrayRHC1Re[irp-1] += dn;
	   ////


	 }//irp
       }//ir
     }// scale > 1

     // energy scale < 1 : 
     // For a given bin, events are gained from bins on the right (higher Ereco bins) and 
     // are lost to bins on the left (lower Ereco bins)
     else if(escale < 1. - kASmallNum) {
       // loop over Ereco
       for(int ir = 1; ir <= n_ereco_bins; ir++) {
	 // nominal reco bin edges
	 //double bin_max = reco_spectrum_tot->GetXaxis()->GetBinUpEdge (ir);
	 //double bin_min = reco_spectrum_tot->GetXaxis()->GetBinLowEdge(ir);
	 double bin_max = util::ErecoBinEdges(s, ir);
	 double bin_min = util::ErecoBinEdges(s, ir-1);
	 // scaled reco bin edges
	 bin_max *= escale;
	 bin_min *= escale;
	 // loop over Ereco' (< Ereco)
	 for(int irp = 1; irp <= ir-1; irp++) {
	   // bin edges for reco bin ir'
	   //	   double bin_max_p = reco_spectrum_tot->GetXaxis()->GetBinUpEdge (irp);
	   //double bin_min_p = reco_spectrum_tot->GetXaxis()->GetBinLowEdge(irp);
	   double bin_max_p = util::ErecoBinEdges(s, irp);
	   double bin_min_p = util::ErecoBinEdges(s, irp-1);
	   // do bins ir and ir' overlap?
	   if(bin_max <= bin_min_p || bin_min >= bin_max_p) continue; // no overlap
	   // calculate overlap
	   double overlap = 0;
	   if(bin_min_p < bin_max && bin_min_p > bin_min) overlap = bin_max   - bin_min_p;
	   if(bin_max_p < bin_max && bin_max_p > bin_min) overlap = bin_max_p - bin_min;
	   if(bin_max_p > bin_max && bin_min_p < bin_min) overlap = bin_max   - bin_min;
	   if(bin_max_p < bin_max && bin_min_p > bin_min) overlap = bin_max_p - bin_min_p;

	   // calculate how many events from nominal bin ir get lost to bin ir'
	   double n  = reco_spectrum_tot_arry_no_ereco_scale[is][ir-1];
	   double dn = n * overlap / (bin_max-bin_min); ///< (overlap) * (nominal # of events) / (bin width)
	   // update the systematically varied templates
	   dn = TMath::Min(dn,n);
 #ifndef PRODUCTIONMODE


	   if (s==kSKsFCFV1RingMuLike) reco_spectrum_tot->SetBinContent(ir,  fMCErecoSpectrumTotArray1Rmu[ir-1]  - dn);
	   // MODIFIED
	   else if (s==kSKsFCFV1RingELike) reco_spectrum_tot->SetBinContent(ir,  fMCErecoSpectrumTotArray1Re[ir-1]  - dn);
	   ////
	   // ANTINU
	   else if (s==kSKsRHCFCFV1RingMuLike) reco_spectrum_tot->SetBinContent(ir,  fMCErecoSpectrumTotArrayRHC1Rmu[ir-1]  - dn);
	   else if (s==kSKsRHCFCFV1RingELike) reco_spectrum_tot->SetBinContent(ir,  fMCErecoSpectrumTotArrayRHC1Re[ir-1]  - dn);
	   ////

	   
	   if (s==kSKsFCFV1RingMuLike) reco_spectrum_tot->SetBinContent(irp, fMCErecoSpectrumTotArray1Rmu[irp-1] + dn);
	   // MODIFIED
	   else if (s==kSKsFCFV1RingELike) reco_spectrum_tot->SetBinContent(irp, fMCErecoSpectrumTotArray1Re[irp-1] + dn);
	   ////
	   // ANTINU
	   else if (s==kSKsRHCFCFV1RingMuLike) reco_spectrum_tot->SetBinContent(irp, fMCErecoSpectrumTotArrayRHC1Rmu[irp-1] + dn);
	   else if (s==kSKsRHCFCFV1RingELike) reco_spectrum_tot->SetBinContent(irp, fMCErecoSpectrumTotArrayRHC1Re[irp-1] + dn);
	   ////

 #endif


	   if (s==kSKsFCFV1RingMuLike) fMCErecoSpectrumTotArray1Rmu[ir-1]  -= dn;
	   // MODIFIED
	   else if (s==kSKsFCFV1RingELike) fMCErecoSpectrumTotArray1Re[ir-1]  -= dn;
	   ////
	   // ANTINU
	   else if (s==kSKsRHCFCFV1RingMuLike) fMCErecoSpectrumTotArrayRHC1Rmu[ir-1]  -= dn;
	   else if (s==kSKsRHCFCFV1RingELike) fMCErecoSpectrumTotArrayRHC1Re[ir-1]  -= dn;
	   ////

	   
	   if (s==kSKsFCFV1RingMuLike) fMCErecoSpectrumTotArray1Rmu[irp-1] += dn;
	   // MODIFIED
	   else if (s==kSKsFCFV1RingELike) fMCErecoSpectrumTotArray1Re[irp-1] += dn;
	   ////
	   // ANTINU
	   else if (s==kSKsRHCFCFV1RingMuLike) fMCErecoSpectrumTotArrayRHC1Rmu[irp-1] += dn;
	   // MODIFIED
	   else if (s==kSKsRHCFCFV1RingELike) fMCErecoSpectrumTotArrayRHC1Re[irp-1] += dn;
	   ////



	 }//irp
       }//ir
     }// scale < 1

    //delete reco_spectrum_tot_no_ereco_scale;

    // loop at the end, to get rid of any negatives
    // loop over Ereco
    for(int ir = 1; ir <= n_ereco_bins; ir++) {
      double nev_tot = 0;




      if (s==kSKsFCFV1RingMuLike) nev_tot = fMCErecoSpectrumTotArray1Rmu[ir-1];
      // MODIFIED 
      else if (s==kSKsFCFV1RingELike) nev_tot = fMCErecoSpectrumTotArray1Re[ir-1];
      ////
      // ANTINU
      else if (s==kSKsRHCFCFV1RingMuLike) nev_tot = fMCErecoSpectrumTotArrayRHC1Rmu[ir-1];
      else if (s==kSKsRHCFCFV1RingELike) nev_tot = fMCErecoSpectrumTotArrayRHC1Re[ir-1];
      ////


      bool MuLike = (s==kSKsFCFV1RingMuLike  && (ir >= 6 && nev_tot < kASmallNum));    // MODIFIED 
      bool ELike = (s==kSKsFCFV1RingELike  && (ir >= 3 && nev_tot < kASmallNum));      // MODIFIED
      bool RHCMuLike = (s==kSKsRHCFCFV1RingMuLike  && (ir >= 6 && nev_tot < kASmallNum));  // ANTINU
      bool RHCELike = (s==kSKsRHCFCFV1RingELike  && (ir >= 3 && nev_tot < kASmallNum));    // ANTINU
      
      //if(nev_tot < 0 || (s==kSKsFCFV1RingMuLike  && (ir >= 6 && nev_tot < kASmallNum)) || (s==kSKsFCFV1RingELike  && (ir >= 3 && nev_tot < kASmallNum))) {
      if(nev_tot < 0 || MuLike || ELike || RHCMuLike || RHCELike) { // MODIFIED ANTINU
	  // bin 5 is first bin expect events in (nominally empty, but can have events for negative energy scale) bin 4 is empty, except for more than 8 sigma tweak!
#ifdef __show_negative_events_mesg_debug__
	cout << "Exp/SKMC/ApplyEnergyScaleSystematicToTotalSpectrum >> set negative number of events " << nev_tot
	     << " to " << kSetNegEventsTo << "\t" << " for Ereco: " << ErecoBinEdges(s,ir-1)+kOneMeV << " in sample = " << util::AsString(s) << endl;
#endif
	applied_negevents_condition = true;
#ifndef PRODUCTIONMODE
	reco_spectrum_tot->SetBinContent(ir, kSetNegEventsTo);
#endif



	if (s==kSKsFCFV1RingMuLike) fMCErecoSpectrumTotArray1Rmu[ir-1] = kSetNegEventsTo;
	// MODIFIED
	else if (s==kSKsFCFV1RingELike) fMCErecoSpectrumTotArray1Re[ir-1] = kSetNegEventsTo;
	////
	// ANTINU
	else if (s==kSKsRHCFCFV1RingMuLike) fMCErecoSpectrumTotArrayRHC1Rmu[ir-1] = kSetNegEventsTo;
	else if (s==kSKsRHCFCFV1RingELike) fMCErecoSpectrumTotArrayRHC1Re[ir-1] = kSetNegEventsTo;
	////


      }
      /*      else if(ir == 5 && nev_tot < kASmallNum) {
#ifdef __show_negative_events_mesg_debug__
	cout << "Exp/SKMC/ApplyEnergyScaleSystematicToTotalSpectrum >> set negative number of events " << reco_spectrum_tot->GetBinContent(ir)
	     << " to " << kSetNegEventsTo << "\t" << " Ereco: " << reco_spectrum_tot->GetXaxis()->GetBinCenter(ir) << endl;
#endif
#ifndef PRODUCTIONMODE
	reco_spectrum_tot->SetBinContent(ir, kSetNegEventsTo);
#endif
        if (s==kSKsFCFV1RingMuLike) fMCErecoSpectrumTotArray1Rmu[ir-1] = kSetNegEventsTo;
	else fMCErecoSpectrumTotArray1Re[ir-1] = kSetNegEventsTo; 
	}*/
    }//ir

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
    double delta_cp = fPhysHypothesis.CurrValue(kPhypDeltaCP);
#ifdef __show_hypothesis__
    cout << "**  Exp/SKMC >> " << fPhysHypothesis;
#endif
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    NuOscDmsqScheme_t  dmsqscheme = Opt::Instance()->gOscDmsqScheme;

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
	bool bin_is_not_empty = this->IsThereAnEvent(s, Etrue);
	if(!bin_is_not_empty) continue;
	
	double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuMu,    kPdgNuMu,    Etrue);
	double p_numubar_numubar = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuMuBar, kPdgNuMuBar, Etrue);
	double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuE,     kPdgNuE,     Etrue);
	double p_nuebar_nuebar   = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuEBar,  kPdgNuEBar,  Etrue);
	double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuMu,    kPdgNuE,     Etrue);
	double p_numubar_nuebar  = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuMuBar, kPdgNuEBar,  Etrue);
	double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuE,     kPdgNuMu,    Etrue);
	double p_nuebar_numubar  = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuEBar,  kPdgNuMuBar, Etrue);
	
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

	    

	    if (s==kSKsFCFV1RingMuLike) nev = fMCTemplArrayFin1Rmu[im][it-1][ir-1];
	    // MODIFIED
	    else if (s==kSKsFCFV1RingELike) nev = fMCTemplArrayFin1Re[im][it-1][ir-1];
	    ////
	    // ANTINU
	    else if (s==kSKsRHCFCFV1RingMuLike) nev = fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1];
	    else if (s==kSKsRHCFCFV1RingELike) nev = fMCTemplArrayFinRHC1Re[im][it-1][ir-1];
	    ////


	    
#ifdef PHYSHYPWEIGHTSAREONE
              posc = 1;
#ifdef PHYSHYPWEIGHTAREONEANDNUEISZERO
              if(util::IsOscNuE(m)) posc = 0;
#endif
#endif
	    // number of events in (Ereco,Etrue) bin (oscillations + systematics applied)
	    double nev_osc = nev * posc;
#ifndef PRODUCTIONMODE
              this -> Template(s,m,kSKMCtFin) -> SetBinContent(ir,it, nev_osc);
#endif



	      if (s==kSKsFCFV1RingMuLike) fMCTemplArrayFin1Rmu[im][it-1][ir-1] = nev_osc;
	      // MODIFIED
	      else if (s==kSKsFCFV1RingELike) fMCTemplArrayFin1Re[im][it-1][ir-1] = nev_osc;
	      ////
	      // ANTINU
	      else if (s==kSKsRHCFCFV1RingMuLike) fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1] = nev_osc;
	      else if (s==kSKsRHCFCFV1RingELike) fMCTemplArrayFinRHC1Re[im][it-1][ir-1] = nev_osc;
	      ////



	  }//ir
	}//im
     }//it
   }//is
  }

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
#ifndef PRODUCTIONMODE
         // get source (systematics) and target (systematics + phys.hypothesis) template
         TH2D * mctempl_tgt = this -> Template(s,m,kSKMCtFin);
         if(!mctempl_tgt) continue; 
#endif
         // loop over Etrue
         for(int it = 1; it <= n_etrue_bins; it++) {
           double Etrue = (util::EtrueBinEdges(s,it-1) + util::EtrueBinEdges(s,it)) / 2;
           double psurv = util::NuMuSurvProbOsc2Flavour(mixing23, dm23sq, Etrue, m, fUnphysicalOscProb);
           // loop over Ereco
	   for(int ir = 1; ir <= n_ereco_bins; ir++) {

	     double nev;
	     if(s==kSKsFCFV1RingMuLike) nev =  fMCTemplArrayFin1Rmu[im][it-1][ir-1];            // MODIFIED
	     else if(s==kSKsRHCFCFV1RingMuLike) nev =  fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]; // ANTINU
	     
	     nev *= psurv;
#ifndef PRODUCTIONMODE
	     mctempl_tgt->SetBinContent(ir,it,nev);
#endif

	     if(s==kSKsFCFV1RingMuLike) fMCTemplArrayFin1Rmu[im][it-1][ir-1] = nev;  // MODIFIED
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
#ifndef PRODUCTIONMODE
         TH2D * mctempl_tgt = this -> Template(s,m,kSKMCtFin);
         if(!mctempl_tgt) continue; 
#endif
         // loop over Etrue
         for(int it = 1; it <= n_etrue_bins; it++) {
           double Etrue = (util::EtrueBinEdges(s,it-1) + util::EtrueBinEdges(s,it)) / 2;
           double psurv = util::NuMuSurvProbDecay(sinsq2theta, alpha, Etrue, m);
           // loop over Ereco
           for(int ir = 1; ir <= n_ereco_bins; ir++) {

	     double nev;
	     if(s==kSKsFCFV1RingMuLike) nev = fMCTemplArrayFin1Rmu[im][it-1][ir-1]; // MODIFIED
	     else if(s==kSKsRHCFCFV1RingMuLike) nev =  fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]; // ANTINU
	     
	     nev *= psurv;
#ifndef PRODUCTIONMODE
	     mctempl_tgt->SetBinContent(ir,it,nev);
#endif


	     if(s==kSKsFCFV1RingMuLike) fMCTemplArrayFin1Rmu[im][it-1][ir-1] = nev;  // MODIFIED
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
	     
	     double nev;
	     if(s==kSKsFCFV1RingMuLike) nev = fMCTemplArrayFin1Rmu[im][it-1][ir-1];  // MODIFIED
             else if(s==kSKsRHCFCFV1RingMuLike) nev =  fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]; // ANTINU
	     
	     nev *= psurv;
#ifndef PRODUCTIONMODE
	     mctempl_tgt->SetBinContent(ir,it,nev);
#endif
	   

	     if(s==kSKsFCFV1RingMuLike) fMCTemplArrayFin1Rmu[im][it-1][ir-1] = nev;   // MODIFIED
	     else if(s==kSKsRHCFCFV1RingMuLike) fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1] = nev; // ANTINU
	     

	   }//ir
         }//it
      }//im
    }//is
  }
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
#ifndef PRODUCTIONMODE
    TH1D * reco_spectrum_tot = this->CurrErecoSpectrum(s);
#endif

    int n_ereco_bins = kNErecoBins[is];
    int n_etrue_bins = kNEtrueBins[is];

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

          //TH2D * mctempl = this -> Template(s,m,kSKMCtFin);
#ifndef PRODUCTIONMODE
          TH1D * reco_spectrum = this->CurrErecoSpectrum(s,m);
#endif

          // loop over Etrue
          double nev = 0;
          for(int it = 1; it <= n_etrue_bins; it++) {
	    if (s==kSKsFCFV1RingMuLike)	nev += fMCTemplArrayFin1Rmu[im][it-1][ir-1];
	    else if (s==kSKsFCFV1RingELike) nev += fMCTemplArrayFin1Re[im][it-1][ir-1]; // MODIFIED
	    else if (s==kSKsRHCFCFV1RingMuLike)	nev += fMCTemplArrayFinRHC1Rmu[im][it-1][ir-1]; // ANTINU
	    else if (s==kSKsRHCFCFV1RingELike) nev += fMCTemplArrayFinRHC1Re[im][it-1][ir-1]; // ANTINU
          }//it
          nev_tot += nev;
#ifdef __apply_minevents_to_individual_modes__
	  if(nev < 0) {
#ifdef __show_negative_events_mesg_debug__
	    cout << "Exp/SKMC/UpdateErecoSpectra >> set negative number of events " << nev
		 << " to kSetNegEventsTo " << kSetNegEventsTo << "\t" << " Ereco: " << reco_spectrum_tot->GetXaxis()->GetBinCenter(ir) << endl;
#endif
	    nev = kSetNegEventsTo;
	  }
#endif //__apply_minevents_to_individual_modes__
#ifndef PRODUCTIONMODE
          reco_spectrum->SetBinContent(ir,nev);
#endif


	  if (s==kSKsFCFV1RingMuLike) 
	    fMCErecoSpectrumArray1Rmu[im][ir-1] = nev;   
	  else if (s==kSKsFCFV1RingELike)
	    fMCErecoSpectrumArray1Re [im][ir-1] = nev;  // MODIFIED
	  else if (s==kSKsRHCFCFV1RingMuLike) 
	    fMCErecoSpectrumArrayRHC1Rmu[im][ir-1] = nev;  // ANTINU 
	  else if (s==kSKsRHCFCFV1RingELike)
	    fMCErecoSpectrumArrayRHC1Re [im][ir-1] = nev;  // ANTINU
	 

       }//im
#ifndef APPLY_ESCALE_AT_SPECTRUM_STAGE

       bool RingMu = (s==kSKsFCFV1RingMuLike && (ir >= 5 && nev_tot < kASmallNum));  // MODIFIED
       bool RingE = (s==kSKsFCFV1RingELike && (ir >= 2 && nev_tot < kASmallNum));    // MODIFIED
       bool RHCRingMu = (s==kSKsRHCFCFV1RingMuLike && (ir >= 5 && nev_tot < kASmallNum));  // ANTINU
       bool RHCRingE = (s==kSKsRHCFCFV1RingELike && (ir >= 2 && nev_tot < kASmallNum));    // ANTINU
       
       //if(nev_tot < 0 || (s==kSKsFCFV1RingMuLike && (ir >= 5 && nev_tot < kASmallNum)) || (s==kSKsFCFV1RingELike && (ir >= 2 && nev_tot < kASmallNum))) {
       if(nev_tot < 0 || RingMu || RingE || RHCRingMu || RHCRingE) {  // ANTINU
	 // bin 5 is first bin expect events in (nominally empty, but can have events for negative energy scale) bin 4 is empty, except for more than 8 sigma tweak!
#ifdef __show_negative_events_mesg_debug__
	 cout << "Exp/SKMC/UpdateErecoSpectra >> set negative number of events " << nev_tot 
	      << " to kSetNegEventsTo " << kSetNegEventsTo << "\t" << " Ereco: " << reco_spectrum_tot->GetXaxis()->GetBinCenter(ir) << endl;
#endif
	 applied_negevents_condition = true;
	 nev_tot = kSetNegEventsTo;
       }
#endif //!APPLY_ESCALE_AT_SPECTRUM_STAGE
#ifndef PRODUCTIONMODE
       reco_spectrum_tot->SetBinContent(ir,nev_tot);
#endif



       if (s==kSKsFCFV1RingMuLike) fMCErecoSpectrumTotArray1Rmu[ir-1] = nev_tot;
       else if (s==kSKsFCFV1RingELike) fMCErecoSpectrumTotArray1Re[ir-1] = nev_tot;  // MODIFIED
       else if (s==kSKsRHCFCFV1RingMuLike) fMCErecoSpectrumTotArrayRHC1Rmu[ir-1] = nev_tot; // ANTINU
       else if (s==kSKsRHCFCFV1RingELike) fMCErecoSpectrumTotArrayRHC1Re[ir-1] = nev_tot;  // ANTINU



    }//ir
  }//is

#ifdef APPLY_ESCALE_AT_SPECTRUM_STAGE
  applied_negevents_condition = this->ApplyEnergyScaleSystematicToTotalSpectrum();
#endif

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
       
#ifndef PRODUCTIONMODE
       true_spectrum_tot->SetBinContent(it,nev_tot);
#endif

       if (s==kSKsFCFV1RingMuLike) fMCEtrueSpectrumTotArray1Rmu[it-1] = nev_tot;
       else if (s==kSKsFCFV1RingELike) fMCEtrueSpectrumTotArray1Re[it-1] = nev_tot;  // MODIFIED
       else if (s==kSKsRHCFCFV1RingMuLike) fMCEtrueSpectrumTotArrayRHC1Rmu[it-1] = nev_tot; // ANTINU
       else if (s==kSKsRHCFCFV1RingELike) fMCEtrueSpectrumTotArrayRHC1Re[it-1] = nev_tot;  // ANTINU

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




  if (s == kSKsFCFV1RingMuLike){
    ev = fMCEtrueSpectrumTotArray1Rmu[it-1] > kASmallNum;
  }
  // MODIFIED
  else if (s == kSKsFCFV1RingELike) {
    ev = fMCEtrueSpectrumTotArray1Re[it-1] > kASmallNum;
  }
  ////
  // ANTINU
  else if (s == kSKsRHCFCFV1RingMuLike){
    ev = fMCEtrueSpectrumTotArrayRHC1Rmu[it-1] > kASmallNum;
  }
  else if (s == kSKsRHCFCFV1RingELike) {
    ev = fMCEtrueSpectrumTotArrayRHC1Re[it-1] > kASmallNum;
  }
  ////




#ifdef __show_isanevent_mesg_debug__
  if(ev){


    if (s == kSKsFCFV1RingMuLike) cout << "Exp/SuperkMC >> There are " << fMCErecoSpectrumTotArray1Rmu[it-1];
    else if (s == kSKsFCFV1RingELike) cout << "Exp/SuperkMC >> There are " << fMCErecoSpectrumTotArray1Re[it-1]; // MODIFIED
    else if (s == kSKsRHCFCFV1RingMuLike) cout << "Exp/SuperkMC >> There are " << fMCErecoSpectrumTotArrayRHC1Rmu[it-1]; // ANTINU
    else if (s == kSKsRHCFCFV1RingELike) cout << "Exp/SuperkMC >> There are " << fMCErecoSpectrumTotArrayRHC1Re[it-1]; // ANTINU
 

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
  else
    return fPOT;  
}
//________________________________________________________________________________________
TTree * SuperkMC::MCTree(SuperkSample_t s, SuperkMode_t m)
{
  if (s==kSKsFCFV1RingMuLike)
  {
     if ( util::IsNuMu(m)      ) return fMCTreeNuMu1Rmu;
     if ( util::IsNuMuBar(m)   ) return fMCTreeNuMuBar1Rmu;
     if ( util::IsNuE(m)       ) return fMCTreeNuE1Rmu;
     if ( util::IsNuEBar(m)    ) return fMCTreeNuEBar1Rmu;
     if ( util::IsOscNuE(m)    ) return fMCTreeOscNuE1Rmu;
  }//1Rmu
  
  if (s==kSKsFCFV1RingELike)
  {
     if ( util::IsNuMu(m)      ) return fMCTreeNuMu1Re;
     if ( util::IsNuMuBar(m)   ) return fMCTreeNuMuBar1Re;
     if ( util::IsNuE(m)       ) return fMCTreeNuE1Re;
     if ( util::IsNuEBar(m)    ) return fMCTreeNuEBar1Re;
     if ( util::IsOscNuE(m)    ) return fMCTreeOscNuE1Re;
  }//1Re


  // ANTINU
  if (s==kSKsRHCFCFV1RingMuLike)
  {
     if ( util::IsNuMu(m)      ) return fMCTreeNuMuRHC1Rmu;
     if ( util::IsNuMuBar(m)   ) return fMCTreeNuMuBarRHC1Rmu;
     if ( util::IsNuE(m)       ) return fMCTreeNuERHC1Rmu;
     if ( util::IsNuEBar(m)    ) return fMCTreeNuEBarRHC1Rmu;
     if ( util::IsOscNuE(m)    ) return fMCTreeOscNuEBarRHC1Rmu;
  }//1Rmu
  if (s==kSKsRHCFCFV1RingELike)
  {
     if ( util::IsNuMu(m)      ) return fMCTreeNuMuRHC1Re;
     if ( util::IsNuMuBar(m)   ) return fMCTreeNuMuBarRHC1Re;
     if ( util::IsNuE(m)       ) return fMCTreeNuERHC1Re;
     if ( util::IsNuEBar(m)    ) return fMCTreeNuEBarRHC1Re;
     if ( util::IsOscNuE(m)    ) return fMCTreeOscNuEBarRHC1Re;
  }//1Re
  ////


  return 0;
}
//________________________________________________________________________________________
double SuperkMC::MCPOT(SuperkSample_t s, SuperkMode_t m)
{
  if(s == kSKsFCFV1RingMuLike || s == kSKsFCFV1RingELike)
  {
     if ( util::IsNuMu(m)      ) return fMCPOTNuMu;
     if ( util::IsNuMuBar(m)   ) return fMCPOTNuMuBar;
     if ( util::IsNuE(m)       ) return fMCPOTNuE;
     if ( util::IsNuEBar(m)    ) return fMCPOTNuEBar;
     if ( util::IsOscNuE(m)    ) return fMCPOTOscNuE;
  }

  
  // ANTINU
  else if(s == kSKsRHCFCFV1RingMuLike || s == kSKsRHCFCFV1RingELike)
  {
     if ( util::IsNuMu(m)      ) return fMCPOTRHCNuMu;
     if ( util::IsNuMuBar(m)   ) return fMCPOTRHCNuMuBar;
     if ( util::IsNuE(m)       ) return fMCPOTRHCNuE;
     if ( util::IsNuEBar(m)    ) return fMCPOTRHCNuEBar;
     if ( util::IsOscNuE(m)    ) return fMCPOTRHCOscNuEBar;
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
      fMCTreeNuMu1Rmu->Draw("Enu");
      fMCTreeNuMu1Rmu->Draw("Enu","Mode==1","same");
      fMCTreeNuMu1Rmu->Draw("Enu","Mode==11||Mode==12||Mode==13","same");
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

  PhysHypothesis_t hyp = fPhysHypothesis.Hypothesis();
  bool hyp_is_not_3flavour_osc = (hyp != kHypOsc3F);

  bool skip = (mode_is_nuesig && hyp_is_not_3flavour_osc);

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
  TFile out(Form("Template_%s.root", name.c_str()),"RECREATE");
  static const int erecobins = kNErecoBins1Rmu;//70;
  static const double * erecobinedges = kErecoBinEdges1Rmu;
  /*static const double erecobinedges[erecobins+1] = {
      0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450, 
      0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950, 
      1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  1.300,  1.350,  1.400,  1.450, 
      1.500,  1.550,  1.600,  1.650,  1.700,  1.750,  1.800,  1.850,  1.900,  1.950, 
      2.000,  2.050,  2.100,  2.150,  2.200,  2.250,  2.300,  2.350,  2.400,  2.450, 
      2.500,  2.550,  2.600,  2.650,  2.700,  2.750,  2.800,  2.850,  2.900,  2.950, 
      3.000,  3.100,  3.200,  3.300,  3.400,  3.500,  4.000,  5.000,  7.000, 10.000,
      30.000
      };*/
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    if(s == kSKsFCFV1RingELike) continue;

    //if(s == kSKsRHCFCFV1RingMuLike) continue;    // ANTINU
    if(s == kSKsRHCFCFV1RingELike) continue;     // ANTINU

    double nsk = 0;
    for(int im = 0; im < knm; im++) {
      SuperkMode_t m = kSKReactionModes[is][im];
      if(m == kSKmUndefined) continue;
      TH2D * mctempl = this -> Template(s,m,t);

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
      TH1D * htru = new TH1D(Form("enu_int_%s",AsStringCompact_KDI(m).c_str()), Form("enu_int_%s",AsStringCompact_KDI(m).c_str()), util::NumVariableXsEtrueBins(s,m), util::VariableXsEtrueBinEdges(s,m));
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

    //get the total ereco spectrum
    TH1D * hfull = this->CurrErecoSpectrum(s);
    hfull->Write("ereco_spectrum");

    //    cout << "Total: " << nsk << endl << endl;
  }//sample
  out.Close();
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

    if(s == kSKsRHCFCFV1RingMuLike) continue;   // ANTINU
    if(s == kSKsRHCFCFV1RingELike) continue;    // ANTINU

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
