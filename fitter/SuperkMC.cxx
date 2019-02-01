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
//#define __show_spectra_update_mesg_debug__c
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
#define __save_templates_debug_mcstat__ //do the extra MCStatErrDebug
#endif
//#define __apply_minevents_to_individual_modes__ //special debug mode that applies kSetNegEventsTo to all modes. Used for comparisons with the KDI 1 sigma spectra
//#define NO_EVENTS_CAN_BE_NEGATIVE // if defined, will stop events being negative at all stages (not just building the Ereco spectra)
// #define FakeDataWeightDebug //Print out fake data weights 
// #define POTWEIGHTSAREONE
// #define PHYSHYPWEIGHTSAREONE
// #define PHYSHYPWEIGHTAREONEANDNUEISZERO //needs PHYSHYPWEIGHTSAREONE to work
   
#define INCLUDEEVENTSOVERTHIRTYGEV // include Ereco overflow events in the last bin?

// NEW 2D TEMPLATES
#define INCLUDEEVENTSNEGATIVEQ2      // include Q2 underflow events in the first bin?
#define INCLUDEEVENTSOVERFLOWQ2      // include Q2 overflow events in the last bin?

#define USE_FLUX_EXT_WEIGHT // Instead of using 'wgt' branch in files, use '_ratio' flux plots.
//#define USE_RPA_WEIGHTS     // use the wgtRPA_nr and wgtRPA_r weights in the trees

//#define USE_OscBr_WEIGHTS   // use the ptheta branch with oscillation weights, for comparison purposes only

//#define __INCLUDE_MC_STAT_FLUCTUATIONS__ // fluctuate bin contents in final MC templates based on MC stat error in each bin

//#define __ALLOW_DIFFERING_HIERARCHIES__ // do we allow our oscillation calculators to have different mass hiearchies to the globals?

#define FULL_2017
//#define NO_WEIGHTS
//#define MACH3_NO_BANFF_TUNE
//#define TREE_WEIGHTS_BeRPA
//#define TREE_WEIGHTS_BeRPA_1

#if defined FULL_2017

   #define USE_SFRFG_WEIGHT    // use the wgtsfrfg weight in the trees - probably should not be combined with an RPA branch
   #define USE_BeRPA_WEIGHTS     // use the wgtBeRPA_nom weights in the trees
   #define USE_CCcoh_WEIGHTS    // use the pion CCcoh weight branch in the trees
   #define USE_NC1gamma_WEIGHTS // use the NC1gamma reweighting by a factor 2
   #define CANCEL_CCCoh_WEIGHTS // Override CCCoh weight for 2p2h and NC1Gamma events

#elif defined TREE_WEIGHTS_BeRPA

   #define SYSTWEIGHTSAREONE
   #define USE_SFRFG_WEIGHT    // use the wgtsfrfg weight in the trees - probably should not be combined with an RPA branch
   #define USE_BeRPA_WEIGHTS     // use the wgtBeRPA_nom weights in the trees
   #define USE_CCcoh_WEIGHTS    // use the pion CCcoh weight branch in the trees
   #define USE_NC1gamma_WEIGHTS // use the NC1gamma reweighting by a factor 2
   #define CANCEL_CCCoh_WEIGHTS

#elif defined TREE_WEIGHTS_BeRPA_1

   #define SYSTWEIGHTSAREONE
   #define USE_SFRFG_WEIGHT    // use the wgtsfrfg weight in the trees - probably should not be combined with an RPA branch
   #define USE_BeRPA_WEIGHTS     // use the wgtBeRPA_nom weights in the trees
   #define USE_CCcoh_WEIGHTS    // use the pion CCcoh weight branch in the trees
   #define USE_NC1gamma_WEIGHTS // use the NC1gamma reweighting by a factor 2
   #define FORCE_BeRPA_WEIGHTS_1
   #define CANCEL_CCCoh_WEIGHTS

#elif defined MACH3_NO_BANFF_TUNE

   #define FLUXWEIGHTSAREONE
   #define BANFFTUNEISONE
   #define XSECWEIGHTSAREONE
   #define SYSTWEIGHTSAREONE
   #define USE_SFRFG_WEIGHT    // use the wgtsfrfg weight in the trees - probably should not be combined with an RPA branch
   #define USE_BeRPA_WEIGHTS     // use the wgtBeRPA_nom weights in the trees
   #define FORCE_SFRFG_WEIGHT_1
   #define FORCE_BeRPA_WEIGHTS_1

#elif defined NO_WEIGHTS

   #define FLUXWEIGHTSAREONE
   #define BANFFTUNEISONE
   #define XSECWEIGHTSAREONE
   #define SYSTWEIGHTSAREONE
   #define USE_SFRFG_WEIGHT    // use the wgtsfrfg weight in the trees - probably should not be combined with an RPA branch
   #define USE_BeRPA_WEIGHTS     // use the wgtBeRPA_nom weights in the trees
   #define FORCE_SFRFG_WEIGHT_1
   #define FORCE_BeRPA_WEIGHTS_1

#endif

//#define FORCE_BeRPA_WEIGHTS_1

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <climits>
#include <iomanip>

#include <TTree.h>
#include <TH2D.h>

#include <TVirtualPad.h>
#include <TLegend.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLorentzVector.h> // NEW 2D TEMPLATES

#include "SuperkMC.h"
#include "SuperkDetectorFSI.h"
#include "NuXSec.h"
#include "NuFlux.h"
#include "Utilities.h"
#include "valor/Conventions/Opt.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Oscillation/NuOscSterile.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Messenger/Messenger.h"

using namespace valor;
using namespace valor::conv;
using namespace valor::util;

ClassImp(SuperkMC)

//________________________________________________________________________________________
SuperkMC::SuperkMC() :
   fMCPOTHistNuMu(0),
   fMCPOTHistNuMuBar(0),
   fMCPOTHistNuE(0),
   fMCPOTHistNuEBar(0),
   fMCPOTHistOscNuE(0),
   fMCPOTHistOscNuEBar(0),
   fMCPOTHistRHCNuMu(0),
   fMCPOTHistRHCNuMuBar(0),
   fMCPOTHistRHCNuE(0),
   fMCPOTHistRHCNuEBar(0),
   fMCPOTHistRHCOscNuE(0),
   fMCPOTHistRHCOscNuEBar(0),
   fMCPOTHistMultiRingNuMu(0),
   fMCPOTHistMultiRingNuMuBar(0),
   fMCPOTHistMultiRingNuE(0),
   fMCPOTHistMultiRingNuEBar(0),
   fMCPOTHistMultiRingOscNuE(0),
   fMCPOTHistMultiRingOscNuEBar(0)
{
    this->Init();
}
//________________________________________________________________________________________
SuperkMC::~SuperkMC()
{
    this->CleanUp();
    
    for (int is = 0 ; is < conv::kns ; is++) {
        if (fMCVariableTemplArray   [is]) { delete fMCVariableTemplArray   [is]; fMCVariableTemplArray   [is] = 0;}
        if (fMCTemplArrayOri        [is]) { delete fMCTemplArrayOri        [is]; fMCTemplArrayOri        [is] = 0;}
        if (fMCTemplArrayNom        [is]) { delete fMCTemplArrayNom        [is]; fMCTemplArrayNom        [is] = 0;}
        if (fMCTemplArrayTun        [is]) { delete fMCTemplArrayTun        [is]; fMCTemplArrayTun        [is] = 0;}
        if (fMCTemplArrayFin        [is]) { delete fMCTemplArrayFin        [is]; fMCTemplArrayFin        [is] = 0;}
        if (fMCTemplArrayRPAWght    [is]) { delete fMCTemplArrayRPAWght    [is]; fMCTemplArrayRPAWght    [is] = 0;}
        if (fMCStatFracErr          [is]) { delete fMCStatFracErr          [is]; fMCStatFracErr          [is] = 0;}
        if (fMCEtrueSpectrumTotArray[is]) { delete fMCEtrueSpectrumTotArray[is]; fMCEtrueSpectrumTotArray[is] = 0;}
        if (fMCErecoSpectrumTotArray[is]) { delete fMCErecoSpectrumTotArray[is]; fMCErecoSpectrumTotArray[is] = 0;}
        if (fMCErecoSpectrumArray   [is]) { delete fMCErecoSpectrumArray   [is]; fMCErecoSpectrumArray   [is] = 0;}
        
        // NEW 2D TEMPLATES
        if (fMCVariableTemplArray2D   [is])     { delete fMCVariableTemplArray2D   [is];     fMCVariableTemplArray2D   [is] = 0;}
        if (fMCTemplArrayOri2D        [is])     { delete fMCTemplArrayOri2D        [is];     fMCTemplArrayOri2D        [is] = 0;}
        if (fMCTemplArrayNom2D        [is])     { delete fMCTemplArrayNom2D        [is];     fMCTemplArrayNom2D        [is] = 0;}
        if (fMCTemplArrayTun2D        [is])     { delete fMCTemplArrayTun2D        [is];     fMCTemplArrayTun2D        [is] = 0;}
        if (fMCTemplArrayFin2D        [is])     { delete fMCTemplArrayFin2D        [is];     fMCTemplArrayFin2D        [is] = 0;}
        if (fMCTemplArrayRPAWght2D    [is])     { delete fMCTemplArrayRPAWght2D    [is];     fMCTemplArrayRPAWght2D    [is] = 0;}
        if (fMCStatFracErr2D          [is])     { delete fMCStatFracErr2D          [is];     fMCStatFracErr2D          [is] = 0;}
        //if (fMCEtrueSpectrumTotArray2D[is])     { delete fMCEtrueSpectrumTotArray2D[is];     fMCEtrueSpectrumTotArray2D[is] = 0;}
        if (fMCErecoSpectrumTotArray2D[is])     { delete fMCErecoSpectrumTotArray2D[is];     fMCErecoSpectrumTotArray2D[is] = 0;}
        if (fMCQ2SpectrumTotArray2D   [is])     { delete fMCQ2SpectrumTotArray2D   [is];     fMCQ2SpectrumTotArray2D   [is] = 0;}
        if (fMCErecoVsQ2SpectrumTotArray2D[is]) { delete fMCErecoVsQ2SpectrumTotArray2D[is]; fMCErecoVsQ2SpectrumTotArray2D[is] = 0;}
        if (fMCErecoVsQ2SpectrumArray [is])     { delete fMCErecoVsQ2SpectrumArray [is];     fMCErecoVsQ2SpectrumArray [is] = 0;}
        //
    }
    
    
}
//________________________________________________________________________________________
void SuperkMC::Init(void)
{
#if defined USE_RPA_WEIGHTS && defined USE_BeRPA_WEIGHTS
  LOG("SuperKMC", pERROR) << "Both RPA and BeRPA weights are being read from the MC tree!";
  util::Die(-1);
#endif

    fPOT = 0;
    fPOTRHC = 0;
    fX2p2h = -1;
    fUnphysicalOscProb = false;
    fTemplateArraysHaveNIWGTune = true;
    
    for(int is=0; is<kns; is++) {
        
        fMCErecoSpectrumTot[is] = 0;
        fMCEtrueSpectrumTot[is] = 0;
        
        fMCErecoVsQ2SpectrumTot[is] = 0; // NEW 2D TEMPLATES
        
        for(int im=0; im<knm; im++) {
            fMCTemplOri      [is][im] = 0;
            fMCTemplNom      [is][im] = 0;
            
            fMCErecoSpectrum [is][im] = 0;
            fMCVariableTempl [is][im] = 0;
            
            // NEW 2D TEMPLATES
            fMCTemplOri2D      [is][im] = 0;
            fMCTemplNom2D      [is][im] = 0;
            fMCTemplPOT2D      [is][im] = 0;
            fMCTemplRPA2D      [is][im] = 0;
            fMCErecoVsQ2Spectrum [is][im] = 0;
            fMCVariableTempl2D [is][im] = 0;
            //
            
            fNNonEmptyVariableBins [is][im] = 0;
            for (int ivb = 0; ivb < kMaxVarBins; ivb++) {
                fNNonEmptyFixedBins    [is][im][ivb] = 0;
                fNonEmptyFixedBins     [is][im][ivb] = NULL;
                fNonEmptyVariableBins  [is][im][ivb] = 0;
            }
            
            // NEW 2D TEMPLATES
            for (int ivb = 0; ivb < kMaxVarBins2D; ivb++) {
                fNNonEmptyFixedBins2D    [is][im][ivb] = 0;
                fNonEmptyFixedBins2D     [is][im][ivb] = NULL;
                fNonEmptyVariableBins2D  [is][im][ivb] = 0;
            }
            //
            
        }//im
        fMCTreeNuMu     .push_back(0);
        fMCTreeNuMuBar  .push_back(0);
        fMCTreeNuE      .push_back(0);
        fMCTreeNuEBar   .push_back(0);
        fMCTreeOscNuE   .push_back(0);
        fMCTreeOscNuEBar.push_back(0);
        
        SuperkSample_t s = kSKSamples[is];
        if (s==conv::kSKsUndefined) {
            fMCVariableTemplArray   .push_back( 0 );
            fMCTemplArrayOri        .push_back( 0 );
            fMCTemplArrayNom        .push_back( 0 );
            fMCTemplArrayTun        .push_back( 0 );
            fMCTemplArrayFin        .push_back( 0 );
            fMCTemplArrayRPAWght    .push_back( 0 );
            fMCStatFracErr          .push_back( 0 );
            fMCEtrueSpectrumTotArray.push_back( 0 );
            fMCErecoSpectrumTotArray.push_back( 0 );
            fMCErecoSpectrumArray   .push_back( 0 );
            
            
            // NEW 2D TEMPLATES
            fMCVariableTemplArray2D       .push_back( 0 );
            fMCTemplArrayOri2D            .push_back( 0 );
            fMCTemplArrayNom2D            .push_back( 0 );
            fMCTemplArrayTun2D            .push_back( 0 );
            fMCTemplArrayFin2D            .push_back( 0 );
            fMCTemplArrayRPAWght2D        .push_back( 0 );
            fMCStatFracErr2D              .push_back( 0 );
            //fMCEtrueSpectrumTotArray2D    .push_back( 0 );
            fMCErecoSpectrumTotArray2D    .push_back( 0 );
            fMCQ2SpectrumTotArray2D       .push_back( 0 );
            fMCErecoVsQ2SpectrumTotArray2D.push_back( 0 );
            fMCErecoVsQ2SpectrumArray     .push_back( 0 );
            //
        }
        else {
            int n_etrue_bins = kNEtrueBins[is];
            int n_ereco_bins = kNErecoBins[is];
            int n_mc_etrue_bins = kMCNEtrueBins;//not sample dependent
            int n_xs_ereco_bins = kNXsErecoBins[is];
            
            // NEW 2D TEMPLATES
            int n_q2_bins = kNQ2Bins[is];
            int n_xs_q2_bins = kNXsQ2Bins[is];
            //
            
            assert(n_etrue_bins < conv::kBigIndex);
            assert(n_ereco_bins < conv::kBigIndex);
            assert( (long) n_ereco_bins * (long)conv::kBigIndex + (long)n_etrue_bins < INT_MAX);
            
            // NEW 2D TEMPLATES
            assert(n_q2_bins < conv::kBigIndex);
            assert( (long) n_q2_bins * (long)conv::kBigIndex * (long)conv::kBigIndex + (long) n_ereco_bins * (long)conv::kBigIndex + (long)n_etrue_bins < INT_MAX);
            //
            
            
            fMCVariableTemplArray   .push_back( new VALORArray(knm, n_mc_etrue_bins, n_xs_ereco_bins) );
            fMCTemplArrayOri        .push_back( new VALORArray(knm, n_etrue_bins   , n_ereco_bins   ) );
            fMCTemplArrayNom        .push_back( new VALORArray(knm, n_etrue_bins   , n_ereco_bins   ) );
            fMCTemplArrayTun        .push_back( new VALORArray(knm, n_etrue_bins   , n_ereco_bins   ) );
            fMCTemplArrayFin        .push_back( new VALORArray(knm, n_etrue_bins   , n_ereco_bins   ) );
            fMCTemplArrayRPAWght    .push_back( new VALORArray(knm, n_etrue_bins   , n_ereco_bins   ) );
            fMCStatFracErr          .push_back( new VALORArray(knm, n_etrue_bins   , n_ereco_bins   ) );
            fMCEtrueSpectrumTotArray.push_back( new VALORArray(n_etrue_bins) );
            fMCErecoSpectrumTotArray.push_back( new VALORArray(n_ereco_bins) );
            fMCErecoSpectrumArray   .push_back( new VALORArray(knm, n_ereco_bins) );
            
            fMCVariableTemplArray   [is]->SetName( Form("mc_templ_variable_%s"    , util::AsStringCompact(s).c_str()));
            fMCTemplArrayOri        [is]->SetName( Form("mc_templ_ori_%s"         , util::AsStringCompact(s).c_str()));
            fMCTemplArrayNom        [is]->SetName( Form("mc_templ_nom_%s"         , util::AsStringCompact(s).c_str()));
            fMCTemplArrayTun        [is]->SetName( Form("mc_templ_tun_%s"         , util::AsStringCompact(s).c_str()));
            fMCTemplArrayFin        [is]->SetName( Form("mc_templ_fin_%s"         , util::AsStringCompact(s).c_str()));
            fMCStatFracErr          [is]->SetName( Form("mc_staterr_%s"           , util::AsStringCompact(s).c_str()));
            fMCEtrueSpectrumTotArray[is]->SetName( Form("mc_spectrum_etrue_tot_%s", util::AsStringCompact(s).c_str()));
            fMCErecoSpectrumTotArray[is]->SetName( Form("mc_spectrum_obs_tot_%s"  , util::AsStringCompact(s).c_str()));
            fMCErecoSpectrumArray   [is]->SetName( Form("mc_spectrum_obs_%s"      , util::AsStringCompact(s).c_str()));
            fMCTemplArrayRPAWght    [is]->SetName( Form("mc_rpa_weight_%s"        , util::AsStringCompact(s).c_str()));
            
            fMCVariableTemplArray   [is]->SetTitle( Form("MC variable template  %s"     , util::AsString(s).c_str()));
            fMCTemplArrayOri        [is]->SetTitle( Form("MC template ori %s"           , util::AsString(s).c_str()));
            fMCTemplArrayNom        [is]->SetTitle( Form("MC template nom %s"           , util::AsString(s).c_str()));
            fMCTemplArrayTun        [is]->SetTitle( Form("MC template tun %s"           , util::AsString(s).c_str()));
            fMCTemplArrayFin        [is]->SetTitle( Form("MC template fin %s"           , util::AsString(s).c_str()));
            fMCStatFracErr          [is]->SetTitle( Form("MC Fractional StatErr %s"     , util::AsString(s).c_str()));
            fMCEtrueSpectrumTotArray[is]->SetTitle( Form("MC Spectrum Etrue_Total %s"   , util::AsString(s).c_str()));
            fMCErecoSpectrumTotArray[is]->SetTitle( Form("MC Spectrum Ereco Total %s"   , util::AsString(s).c_str()));
            fMCErecoSpectrumArray   [is]->SetTitle( Form("MC Spectrum Ereco per mode %s", util::AsString(s).c_str()));
            fMCTemplArrayRPAWght    [is]->SetTitle( Form("MC RPA Weight %s" , util::AsString(s).c_str()));
            
            fMCTemplArrayRPAWght    [is]->SetAll(1);
            
            
            
            // NEW 2D TEMPLATES
            fMCVariableTemplArray2D       .push_back( new VALORArray(knm, n_mc_etrue_bins, n_xs_ereco_bins, n_xs_q2_bins) );
            fMCTemplArrayOri2D            .push_back( new VALORArray(knm, n_etrue_bins   , n_ereco_bins   , n_q2_bins   ) );
            fMCTemplArrayNom2D            .push_back( new VALORArray(knm, n_etrue_bins   , n_ereco_bins   , n_q2_bins  ) );
            fMCTemplArrayTun2D            .push_back( new VALORArray(knm, n_etrue_bins   , n_ereco_bins   , n_q2_bins  ) );
            fMCTemplArrayFin2D            .push_back( new VALORArray(knm, n_etrue_bins   , n_ereco_bins   , n_q2_bins   ) );
            fMCTemplArrayRPAWght2D        .push_back( new VALORArray(knm, n_etrue_bins   , n_ereco_bins   , n_q2_bins   ) );
            fMCStatFracErr2D              .push_back( new VALORArray(knm, n_etrue_bins   , n_ereco_bins   , n_q2_bins   ) );
            //fMCEtrueSpectrumTotArray2D    .push_back( new VALORArray(n_etrue_bins) );
            fMCErecoSpectrumTotArray2D    .push_back( new VALORArray(n_ereco_bins) );
            fMCQ2SpectrumTotArray2D       .push_back( new VALORArray(n_q2_bins) );
            fMCErecoVsQ2SpectrumTotArray2D.push_back( new VALORArray(n_ereco_bins, n_q2_bins) );
            fMCErecoVsQ2SpectrumArray     .push_back( new VALORArray(knm, n_ereco_bins, n_q2_bins) );
            
            fMCVariableTemplArray2D       [is]->SetName( Form("mc_templ_variable_2d_%s"        , util::AsStringCompact(s).c_str()));
            fMCTemplArrayOri2D            [is]->SetName( Form("mc_templ_ori_2d_%s"             , util::AsStringCompact(s).c_str()));
            fMCTemplArrayNom2D            [is]->SetName( Form("mc_templ_nom_2d_%s"             , util::AsStringCompact(s).c_str()));
            fMCTemplArrayTun2D            [is]->SetName( Form("mc_templ_tun_2d_%s"             , util::AsStringCompact(s).c_str()));
            fMCTemplArrayFin2D            [is]->SetName( Form("mc_templ_fin_2d_%s"             , util::AsStringCompact(s).c_str()));
            fMCStatFracErr2D              [is]->SetName( Form("mc_staterr_2d_%s"               , util::AsStringCompact(s).c_str()));
            //fMCEtrueSpectrumTotArray2D    [is]->SetName( Form("mc_spectrum_etrue_tot_2d_%s"    , util::AsStringCompact(s).c_str()));
            fMCErecoSpectrumTotArray2D    [is]->SetName( Form("mc_spectrum_obs_tot_ereco_2d_%s", util::AsStringCompact(s).c_str()));
            fMCErecoVsQ2SpectrumTotArray2D[is]->SetName( Form("mc_spectrum_obs_tot_erecoVsq2_2d_%s", util::AsStringCompact(s).c_str()));
            fMCQ2SpectrumTotArray2D       [is]->SetName( Form("mc_spectrum_obs_tot_q2_2d_%s"   , util::AsStringCompact(s).c_str()));
            fMCErecoVsQ2SpectrumArray     [is]->SetName( Form("mc_spectrum_obs_2d_%s"          , util::AsStringCompact(s).c_str()));
            fMCTemplArrayRPAWght2D        [is]->SetName( Form("mc_rpa_weight_2d_%s"            , util::AsStringCompact(s).c_str()));
            
            fMCVariableTemplArray2D       [is]->SetTitle( Form("MC variable template 2D  %s"        , util::AsString(s).c_str()));
            fMCTemplArrayOri2D            [is]->SetTitle( Form("MC template ori 2D %s"              , util::AsString(s).c_str()));
            fMCTemplArrayNom2D            [is]->SetTitle( Form("MC template nom 2D %s"              , util::AsString(s).c_str()));
            fMCTemplArrayTun2D            [is]->SetTitle( Form("MC template tun 2D %s"              , util::AsString(s).c_str()));
            fMCTemplArrayFin2D            [is]->SetTitle( Form("MC template fin 2D %s"              , util::AsString(s).c_str()));
            fMCStatFracErr2D              [is]->SetTitle( Form("MC Fractional StatErr 2D %s"        , util::AsString(s).c_str()));
            //fMCEtrueSpectrumTotArray2D    [is]->SetTitle( Form("MC Spectrum Etrue_Total 2D %s"      , util::AsString(s).c_str()));
            fMCErecoSpectrumTotArray2D    [is]->SetTitle( Form("MC Spectrum Ereco Total 2D %s"      , util::AsString(s).c_str()));
            fMCErecoVsQ2SpectrumTotArray2D[is]->SetTitle( Form("MC Spectrum Ereco Vs Q2 Total 2D %s"      , util::AsString(s).c_str()));
            fMCQ2SpectrumTotArray2D       [is]->SetTitle( Form("MC Spectrum Q2 Total 2D %s"         , util::AsString(s).c_str()));
            fMCErecoVsQ2SpectrumArray     [is]->SetTitle( Form("MC Spectrum Ereco Vs Q2 per mode %s", util::AsString(s).c_str()));
            fMCTemplArrayRPAWght2D        [is]->SetTitle( Form("MC RPA Weight %s"       , util::AsString(s).c_str()));
            
            fMCTemplArrayRPAWght2D        [is]->SetAll(1);
            //
            
        } // not kUndfined
        
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
    fBrWghtFlx           = 1;
    fBrWghtSFRFG         = 1;
    fBrWghtRPArel        = 1;
    fBrWghtRPAnonrel     = 1;
    fBrWghtCCcoh         = 1;
    fBrWghtOsc           = 1; // ptheta oscillations branch
    
    fNuFlux              = 0;
    fNuXSec              = 0;
    fSKDetFSI            = 0;
    
    arr_index = util::ArrIdx::Instance();
    
#ifdef USE_CCcoh_WEIGHTS
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC CCcoherent reweight is applied";
#endif   
#ifdef USE_NC1gamma_WEIGHTS 
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC NC1gamma reweight is applied";
#endif
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> Init complete";
    
}
//________________________________________________________________________________________
void SuperkMC::CleanUp(void)
{
    for(int is = 0; is < kns; is++) {
        for(int im = 0; im < knm; im++) {
            //already deleted Ori & Nom (the if doesn't work!)
            if( fMCTemplOri2D[is][im] )
            {
               delete fMCTemplOri2D[is][im];
               fMCTemplOri2D[is][im] = 0;
            }
            if( fMCTemplNom2D[is][im] )
            {
               delete fMCTemplNom2D[is][im];
               fMCTemplNom2D[is][im] = 0;
            }
            if( fMCTemplPOT2D[is][im] )
            {
               delete fMCTemplPOT2D[is][im];
               fMCTemplPOT2D[is][im] = 0;
            }
            if( fMCTemplRPA2D[is][im] )
            {
               delete fMCTemplRPA2D[is][im];
               fMCTemplRPA2D[is][im] = 0;
            }

            if( fMCErecoVsQ2Spectrum[is][im] )
            {
               delete fMCErecoVsQ2Spectrum[is][im];
               fMCErecoVsQ2Spectrum[is][im] = 0;
            }
            
            if( fMCErecoSpectrum[is][im] ) delete fMCErecoSpectrum[is][im];
            if( fMCVariableTempl[is][im] ) delete fMCVariableTempl[is][im];
            
            // NEW 2D TEMPLATES
            if( fMCErecoVsQ2Spectrum[is][im] ) delete fMCErecoVsQ2Spectrum[is][im];
            if( fMCVariableTempl2D  [is][im] ) delete fMCVariableTempl2D[is][im];
            //
            
            for (int ivb = 0; ivb < kMaxVarBins ; ivb++) {
                if (fNonEmptyFixedBins[is][im][ivb]) delete[] fNonEmptyFixedBins[is][im][ivb];
            }
            
            // NEW 2D TEMPLATES
            for (int ivb = 0; ivb < kMaxVarBins2D ; ivb++) {
                if (fNonEmptyFixedBins2D[is][im][ivb]) delete[] fNonEmptyFixedBins2D[is][im][ivb];
            }
            //
            
        } //im
        if( fMCErecoSpectrumTot [is] ) { delete fMCErecoSpectrumTot [is]; fMCErecoSpectrumTot [is] = 0; }
        if( fMCEtrueSpectrumTot [is] ) { delete fMCEtrueSpectrumTot [is]; fMCEtrueSpectrumTot [is] = 0; }
        
        // NEW 2D TEMPLATES
        if( fMCErecoVsQ2SpectrumTot [is] ) { delete fMCErecoVsQ2SpectrumTot [is]; fMCErecoVsQ2SpectrumTot [is] = 0; }
        //
    }

    if(fMCPOTHistNuMu) delete fMCPOTHistNuMu;
    if(fMCPOTHistNuMuBar) delete fMCPOTHistNuMuBar;
    if(fMCPOTHistNuE) delete fMCPOTHistNuE;
    if(fMCPOTHistNuEBar) delete fMCPOTHistNuEBar;
    if(fMCPOTHistOscNuE) delete fMCPOTHistOscNuE;
    if(fMCPOTHistOscNuEBar) delete fMCPOTHistOscNuEBar;

    if(fMCPOTHistRHCNuMu) delete fMCPOTHistRHCNuMu;
    if(fMCPOTHistRHCNuMuBar) delete fMCPOTHistRHCNuMuBar;
    if(fMCPOTHistRHCNuE) delete fMCPOTHistRHCNuE;
    if(fMCPOTHistRHCNuEBar) delete fMCPOTHistRHCNuEBar;
    if(fMCPOTHistRHCOscNuE) delete fMCPOTHistRHCOscNuE;
    if(fMCPOTHistRHCOscNuEBar) delete fMCPOTHistRHCOscNuEBar;
 
    if(fMCPOTHistMultiRingNuMu) delete fMCPOTHistMultiRingNuMu;
    if(fMCPOTHistMultiRingNuMuBar) delete fMCPOTHistMultiRingNuMuBar;
    if(fMCPOTHistMultiRingNuE) delete fMCPOTHistMultiRingNuE;
    if(fMCPOTHistMultiRingNuEBar) delete fMCPOTHistMultiRingNuEBar;
    if(fMCPOTHistMultiRingOscNuE) delete fMCPOTHistMultiRingOscNuE;
    if(fMCPOTHistMultiRingOscNuEBar) delete fMCPOTHistMultiRingOscNuEBar;
}
//________________________________________________________________________________________
bool SuperkMC::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
    this->CleanUp();
    
    // instantiate all necessary histograms
    this->Book();
  
  if(util::UtlMakeFakeData::Instance()->IsActive()){
    util::UtlMakeFakeData::Instance()->Read(xml,root_node);
  }
  
#ifdef __show_init_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Reading SuperK MC inputs from input XML file.";
#endif
    
    //
    // Get data POT
    //
    fPOT = util::XML2Dbl(xml,root_node,"data_pot/superk/fhc");
#ifdef __show_init_mesg_info__
    LOG("SuperkMC",pINFO) << "POT_{data} (FHC) = " << fPOT;
#endif
    //
    // Get data POT RHC
    //
    fPOTRHC = util::XML2Dbl(xml,root_node,"data_pot/superk/rhc");
#ifdef __show_init_mesg_info__
    LOG("SuperkMC",pINFO) << "POT_{data} (RHC) = " << fPOTRHC;
#endif
    
    // NEW MULTIRING
    
    //
    // Get data POT FHC MultiRing
    //
    
#ifdef __MULTIRING_FHC_RUNNING__
    fPOTMultiRing = util::XML2Dbl(xml,root_node,"data_pot/superk/multiring_fhc");
#ifdef __show_init_mesg_info__
    LOG("SuperkMC",pINFO) << "POT_{data} (MultiRing) = " << fPOTMultiRing;
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
    
    fMCPOTRHCNuMu      = fMCPOTHistRHCNuMu     ->GetBinContent(2) / fMCPOTHistRHCNuMu     ->GetBinContent(1);
    fMCPOTRHCNuMuBar   = fMCPOTHistRHCNuMuBar  ->GetBinContent(2) / fMCPOTHistRHCNuMuBar  ->GetBinContent(1);
    fMCPOTRHCNuE       = fMCPOTHistRHCNuE      ->GetBinContent(2) / fMCPOTHistRHCNuE      ->GetBinContent(1);
    fMCPOTRHCNuEBar    = fMCPOTHistRHCNuEBar   ->GetBinContent(2) / fMCPOTHistRHCNuEBar   ->GetBinContent(1);
    fMCPOTRHCOscNuE    = fMCPOTHistRHCOscNuE   ->GetBinContent(2) / fMCPOTHistRHCOscNuE   ->GetBinContent(1);
    fMCPOTRHCOscNuEBar = fMCPOTHistRHCOscNuEBar->GetBinContent(2) / fMCPOTHistRHCOscNuEBar->GetBinContent(1);
   
    // NEW MULTIRING
#ifdef __MULTIRING_FHC_RUNNING__
    fMCPOTHistMultiRingNuMu       = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numu/pot_fhc_cc1pi");
    fMCPOTHistMultiRingNuMuBar    = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numubar/pot_fhc_cc1pi");
    fMCPOTHistMultiRingNuE        = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nue/pot_fhc_cc1pi");
    fMCPOTHistMultiRingNuEBar     = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nuebar/pot_fhc_cc1pi");
    fMCPOTHistMultiRingOscNuE     = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/signal_nue/pot_fhc_cc1pi");
    fMCPOTHistMultiRingOscNuEBar  = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/signal_nuebar/pot_fhc_cc1pi");
    
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
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "FHC POT {MC; numu}    = "    << fMCPOTNuMu;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "FHC POT {MC; numubar} = "    << fMCPOTNuMuBar;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "FHC POT {MC; nue}     = "    << fMCPOTNuE;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "FHC POT {MC; nuebar}  = "    << fMCPOTNuEBar;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "FHC POT {MC; osc_nue} = "    << fMCPOTOscNuE;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "FHC POT {MC; osc_nuebar} = " << fMCPOTOscNuEBar;
    
    // RHC
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "RHC POT {MC; numu}    = "    << fMCPOTRHCNuMu;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "RHC POT {MC; numubar} = "    << fMCPOTRHCNuMuBar;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "RHC POT {MC; nue}     = "    << fMCPOTRHCNuE;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "RHC POT {MC; nuebar}  = "    << fMCPOTRHCNuEBar;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "RHC POT {MC; osc_nue} = "    << fMCPOTRHCOscNuE;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "RHC POT {MC; osc_nuebar} = " << fMCPOTRHCOscNuEBar;
    
    // NEW MULTIRING
#ifdef __MULTIRING_FHC_RUNNING__
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "MultiRing: POT {MC; numu}       = " << fMCPOTMultiRingNuMu;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "MultiRing: POT {MC; numubar}    = " << fMCPOTMultiRingNuMuBar;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "MultiRing: POT {MC; nue}        = " << fMCPOTMultiRingNuE;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "MultiRing: POT {MC; nuebar}     = " << fMCPOTMultiRingNuEBar;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "MultiRing: POT {MC; osc_nue}    = " << fMCPOTMultiRingOscNuE;
    LOG("SuperkMC",pINFO) << "Exp/SKMC >> " << "MultiRing: POT {MC; osc_nuebar} = " << fMCPOTMultiRingOscNuEBar;
#endif
    ////
    
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if (s==kSKsUndefined) continue;
        
        LOG("SuperkMC",pDEBUG) << "Nev {MC; numu;     "   <<AsStringCompact(s)<<"} = " << fMCTreeNuMu     [is]-> GetEntries();
        LOG("SuperkMC",pDEBUG) << "Nev {MC; numubar;  "   <<AsStringCompact(s)<<"} = " << fMCTreeNuMuBar  [is]-> GetEntries();
        LOG("SuperkMC",pDEBUG) << "Nev {MC; nue;      "   <<AsStringCompact(s)<<"} = " << fMCTreeNuE      [is]-> GetEntries();
        LOG("SuperkMC",pDEBUG) << "Nev {MC; nuebar;   "   <<AsStringCompact(s)<<"} = " << fMCTreeNuEBar   [is]-> GetEntries();
        LOG("SuperkMC",pDEBUG) << "Nev {MC; nue[osc]; "   <<AsStringCompact(s)<<"} = " << fMCTreeOscNuE   [is]-> GetEntries();
        LOG("SuperkMC",pDEBUG) << "Nev {MC; nuebar[osc]; "<<AsStringCompact(s)<<"} = " << fMCTreeOscNuEBar[is]-> GetEntries();
        
    }
#endif//__show_init_mesg_info__
  
  if(util::UtlMakeFakeData::Instance()->IsActive()){
    LOG("SuperkMC",pDEBUG) << "Reweighting events for " << util::UtlMakeFakeData::Instance()->GetReweightStyle();
  }
  
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
            this->FillTemplateHistogram(s,m);
        }
    }
    
    
    this->FillTemplateArrays(true); // Fill the arrays - with niwg weights by default
    
    // fill-in "special" MC templates (osc. numu CC)
    this->FillSpecialTemplates(); // NEW 2D TEMPLATES (not for 2D templates)
    
#ifdef __save_templates_debug__
#ifndef __USE_2D_TEMPLATES__
    this->NevDebug(kSKMCtOri, "nevents");
    this->NevDebug(kSKMCtNom, "mcnorm");
    this->NevDebug(kSKMCtRPA, "mcnorm_rpa");
#else
    // NEW 2D TEMPLATES
    this->NevDebug2D(kSKMCtOri, "nevents");
    this->NevDebug2D(kSKMCtNom, "mcnorm");
    this->NevDebug2D(kSKMCtRPA, "mcnorm_rpa");
    //
#endif
#endif
    
    
#ifndef __USE_2D_TEMPLATES__
    this->CalculateMCStatError(); // NEW 2D TEMPLATES (no update for 2D templates)
#ifdef __save_templates_debug_mcstat__
    this->MCStatErrDebug();
#endif
#endif
    
    
    
    // need to read in the xsec splines, before applying any BANFF tweaks
    
#ifndef __USE_2D_TEMPLATES__ // 1D templates
    
    fNuXSec -> BuildXsBinningMCTemplates(fMCTemplOri);
    
#else // 2D templates
    
    // Initialize fMCTemplOri 1D template, otherwise it crashes
    // Only needed to build the fXsBinningTempl in NuXSec.cxx, which are 1D
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        const int      n_ereco_bins    = kNErecoBins[is];
        const int      n_etrue_bins    = kNEtrueBins[is];
        const double * ereco_bin_edges = util::ErecoBinEdges(s);
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            fMCTemplOri[is][im] = new TH2D(Form("nsk_ori_%d_%d",is,im), util::AsString(m).c_str(),
                                           n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
            fMCTemplOri[is][im]->SetDirectory(0);
            fMCTemplOri[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
            fMCTemplOri[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
        }
    }
    //this->ProjectTemplate3Dto2D(fMCTemplOri2D,fMCTemplOri);
    
    fNuXSec -> BuildXsBinningMCTemplates(fMCTemplOri2D);
    
#endif
    
#ifndef __USE_2D_TEMPLATES__
    if( !fNuXSec -> ReadSplines(xml,root_node) ) return false; // old 1d
#else
    if( !fNuXSec -> ReadSplines2D(xml,root_node) ) return false; // new 2d
#endif
    
    // apply flux tuning (based on external data) to templates
    //this->ApplyFluxExtTuningToTemplates();
    
    
    // apply flux and cross-section tuning (based on nd280 data fit by BANFF) to templates
#ifndef __USE_2D_TEMPLATES__
    this->ApplyBANFFTuningToTemplates();
#else
    this->ApplyBANFFTuningToTemplates2D(); // NEW 2D TEMPLATES
#endif
    
    // fill the etrue spectra
#ifndef __USE_2D_TEMPLATES__
    this->FillEtrueSpectra();
#else
    this->FillEtrueSpectra2D(); // NEW 2D TEMPLATES
#endif
    
#ifdef __save_templates_debug__
#ifndef __USE_2D_TEMPLATES__
    this->NevDebug(kSKMCtTun, "mcnorm_fluxtuned_banfftuned");
#else
    this->NevDebug2D(kSKMCtTun, "mcnorm_fluxtuned_banfftuned"); // NEW 2D TEMPLATES
#endif
#endif
    
#ifdef __show_init_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> Done reading SuperK MC inputs!";
#endif
    
    return true;
}
//________________________________________________________________________________________
bool SuperkMC::ReadSmearedMomentum(TXMLEngine* xml, XMLNodePointer_t root_node)
{
    this->CleanUp();
    
    // instantiate all necessary histograms
    this->Book();
    

#ifdef __show_init_mesg_info__
    // FHC
    LOG("SuperkMC",pDEBUG) << "FHC POT {MC; numu}    = "    << fMCPOTNuMu;
    LOG("SuperkMC",pDEBUG) << "FHC POT {MC; numubar} = "    << fMCPOTNuMuBar;
    LOG("SuperkMC",pDEBUG) << "FHC POT {MC; nue}     = "    << fMCPOTNuE;
    LOG("SuperkMC",pDEBUG) << "FHC POT {MC; nuebar}  = "    << fMCPOTNuEBar;
    LOG("SuperkMC",pDEBUG) << "FHC POT {MC; osc_nue} = "    << fMCPOTOscNuE;
    LOG("SuperkMC",pDEBUG) << "FHC POT {MC; osc_nuebar} = " << fMCPOTOscNuEBar;
    
    // RHC
    LOG("SuperkMC",pDEBUG) << "RHC POT {MC; numu}    = "    << fMCPOTRHCNuMu;
    LOG("SuperkMC",pDEBUG) << "RHC POT {MC; numubar} = "    << fMCPOTRHCNuMuBar;
    LOG("SuperkMC",pDEBUG) << "RHC POT {MC; nue}     = "    << fMCPOTRHCNuE;
    LOG("SuperkMC",pDEBUG) << "RHC POT {MC; nuebar}  = "    << fMCPOTRHCNuEBar;
    LOG("SuperkMC",pDEBUG) << "RHC POT {MC; osc_nue} = "    << fMCPOTRHCOscNuE;
    LOG("SuperkMC",pDEBUG) << "RHC POT {MC; osc_nuebar} = " << fMCPOTRHCOscNuEBar;
    
    // NEW MULTIRING
#ifdef __MULTIRING_FHC_RUNNING__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> " << "MultiRing: POT {MC; numu}       = " << fMCPOTMultiRingNuMu;
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> " << "MultiRing: POT {MC; numubar}    = " << fMCPOTMultiRingNuMuBar;
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> " << "MultiRing: POT {MC; nue}        = " << fMCPOTMultiRingNuE;
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> " << "MultiRing: POT {MC; nuebar}     = " << fMCPOTMultiRingNuEBar;
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> " << "MultiRing: POT {MC; osc_nue}    = " << fMCPOTMultiRingOscNuE;
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> " << "MultiRing: POT {MC; osc_nuebar} = " << fMCPOTMultiRingOscNuEBar;
#endif
    ////
    
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if (s==kSKsUndefined) continue;
        
        LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> " << "Nev {MC; numu;     "   <<AsStringCompact(s)<<"} = " << fMCTreeNuMu     [is]-> GetEntries()  ;
        LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> " << "Nev {MC; numubar;  "   <<AsStringCompact(s)<<"} = " << fMCTreeNuMuBar  [is]-> GetEntries()  ;
        LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> " << "Nev {MC; nue;      "   <<AsStringCompact(s)<<"} = " << fMCTreeNuE      [is]-> GetEntries()  ;
        LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> " << "Nev {MC; nuebar;   "   <<AsStringCompact(s)<<"} = " << fMCTreeNuEBar   [is]-> GetEntries()  ;
        LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> " << "Nev {MC; nue[osc]; "   <<AsStringCompact(s)<<"} = " << fMCTreeOscNuE   [is]-> GetEntries()  ;
        LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> " << "Nev {MC; nuebar[osc]; "<<AsStringCompact(s)<<"} = " << fMCTreeOscNuEBar[is]-> GetEntries()  ;
        
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
            this->FillTemplateHistogramSmeared(s,m);
        }
    }
    
    this->FillTemplateArrays(true); // Fill the arrays - with niwg weights by default
    
    // fill-in "special" MC templates (osc. numu CC)
    this->FillSpecialTemplates(); // NEW 2D TEMPLATES (not for 2D templates)
    
#ifdef __save_templates_debug__
#ifndef __USE_2D_TEMPLATES__
    this->NevDebug(kSKMCtOri, "nevents");
    this->NevDebug(kSKMCtNom, "mcnorm");
    this->NevDebug(kSKMCtRPA, "mcnorm_rpa");
#else
    // NEW 2D TEMPLATES
    this->NevDebug2D(kSKMCtOri, "nevents");
    this->NevDebug2D(kSKMCtNom, "mcnorm");
    this->NevDebug2D(kSKMCtRPA, "mcnorm_rpa");
    //
#endif
#endif
    
    
#ifndef __USE_2D_TEMPLATES__
    this->CalculateMCStatError(); // NEW 2D TEMPLATES (no update for 2D templates)
#ifdef __save_templates_debug_mcstat__
    this->MCStatErrDebug();
#endif
#endif
    
    
    // need to read in the xsec splines, before applying any BANFF tweaks
    
#ifndef __USE_2D_TEMPLATES__ // 1D templates
    
    fNuXSec -> BuildXsBinningMCTemplates(fMCTemplOri);
    
#else // 2D templates
    
    // Initialize fMCTemplOri 1D template, otherwise it crashes
    // Only needed to build the fXsBinningTempl in NuXSec.cxx, which are 1D
    
    // Don't need this any more as we have 2D splines
    
    /*
     for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     const int      n_ereco_bins    = kNErecoBins[is];
     const int      n_etrue_bins    = kNEtrueBins[is];
     const double * ereco_bin_edges = util::ErecoBinEdges(s);
     for(int im = 0; im < knm; im++) {
     SuperkMode_t m = kSKReactionModes[is][im];
     if(m == kSKmUndefined) continue;
     fMCTemplOri[is][im] = new TH2D(Form("nsk_ori_%d_%d",is,im), util::AsString(m).c_str(),
     n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
     fMCTemplOri[is][im]->SetDirectory(0);
     fMCTemplOri[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
     fMCTemplOri[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
     }
     }
     */
    //this->ProjectTemplate3Dto2D(fMCTemplOri2D,fMCTemplOri);
    
    fNuXSec -> BuildXsBinningMCTemplates(fMCTemplOri2D);
    
#endif
    
    //#ifndef __USE_2D_TEMPLATES__
    //    if( !fNuXSec -> ReadSplines(xml,root_node) ) return false; // old 1d
    //#else
    //    if( !fNuXSec -> ReadSplines2D(xml,root_node) ) return false; // new 2d
    //#endif
    
    // apply flux tuning (based on external data) to templates
    //this->ApplyFluxExtTuningToTemplates();
    
    
    // apply flux and cross-section tuning (based on nd280 data fit by BANFF) to templates
    //#ifndef __USE_2D_TEMPLATES__
    //    this->ApplyBANFFTuningToTemplates();
    //#else
    //    this->ApplyBANFFTuningToTemplates2D(); // NEW 2D TEMPLATES
    //#endif
    
    // fill the etrue spectra
#ifndef __USE_2D_TEMPLATES__
    this->FillEtrueSpectra();
#else
    this->FillEtrueSpectra2D(); // NEW 2D TEMPLATES
#endif
    
#ifdef __save_templates_debug__
#ifndef __USE_2D_TEMPLATES__
    this->NevDebug(kSKMCtTun, "mcnorm_fluxtuned_banfftuned");
#else
    this->NevDebug2D(kSKMCtTun, "mcnorm_fluxtuned_banfftuned"); // NEW 2D TEMPLATES
#endif
#endif
    
#ifdef __show_init_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> Done reading SuperK MC inputs!" ;
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
    
#ifndef __USE_2D_TEMPLATES__
    
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
    
    
    
#else  // __USE_2D_TEMPLATES__
    
    // NEW 2D TEMPLATES
    
    for(int is = 0; is < kns; is++)
    {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        const int      n_ereco_bins    = kNErecoBins[is];
        const double * ereco_bin_edges = util::ErecoBinEdges(s);
        const int      n_q2_bins    = kNQ2Bins[is];
        const double * q2_bin_edges = util::Q2BinEdges(s);
        
        fMCErecoVsQ2SpectrumTot[is] = new TH2D("","", n_q2_bins, q2_bin_edges, n_ereco_bins, ereco_bin_edges);
        fMCErecoVsQ2SpectrumTot[is] -> SetDirectory(0);
        fMCErecoVsQ2SpectrumTot[is]->SetLineWidth(2);
        fMCErecoVsQ2SpectrumTot[is]->SetLineColor(kBlack);
        
        // Use the same Etrue template as for 1D templates mode
        fMCEtrueSpectrumTot[is] = new TH1D("","", kNEtrueBins[is], util::EtrueBinEdges(s));
        fMCEtrueSpectrumTot[is] -> SetDirectory(0);
        fMCEtrueSpectrumTot[is]->SetLineWidth(2);
        fMCEtrueSpectrumTot[is]->SetLineColor(kBlack);
        
        for(int im=0; im < knm; im++)
        {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            fMCErecoVsQ2Spectrum[is][im] = new TH2D("","", n_q2_bins, q2_bin_edges, n_ereco_bins, ereco_bin_edges);
            fMCErecoVsQ2Spectrum[is][im] -> SetDirectory(0);
            //util::SetDefaultSKSpectrumStyle(s,m,fMCErecoVsQ2Spectrum[is][im]);
        }//m
    }//s
    
    //
    
#endif // __USE_2D_TEMPLATES__
    
    
    
    //
    // 2-D(Ereco,Etrue) MC templates
    // Histograms corresponding to different samples can have different binning
    //
    
#ifndef __USE_2D_TEMPLATES__
    
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
            fMCTemplOri[is][im] = new TH2D(Form("nsk_ori_%d_%d",is,im), util::AsString(m).c_str(),
                                           n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
            fMCTemplNom[is][im] = new TH2D(Form("nsk_nom_%d_%d",is,im), util::AsString(m).c_str(),
                                           n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
            fMCTemplRPA[is][im] = new TH2D(Form("nsk_rpa_%d_%d",is,im), util::AsString(m).c_str(),
                                           n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
            fMCTemplPOT[is][im] = new TH2D(Form("nsk_pot_%d_%d",is,im), util::AsString(m).c_str(),
                                           n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
            
            fMCTemplOri[is][im]->SetDirectory(0);
            fMCTemplOri[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
            fMCTemplOri[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
            fMCTemplNom[is][im]->SetDirectory(0);
            fMCTemplNom[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
            fMCTemplNom[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
            fMCTemplRPA[is][im]->SetDirectory(0);
            fMCTemplRPA[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
            fMCTemplRPA[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
            fMCTemplPOT[is][im]->SetDirectory(0);
            fMCTemplPOT[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
            fMCTemplPOT[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
            
            fMCVariableTempl[is][im] = new TH2D(Form("nsk_ori_var_%d_%d",is,im), util::AsString(m).c_str(),
                                                var_n_ereco_bins, var_ereco_bin_edges, var_n_etrue_bins,  var_etrue_bin_edges);
            fMCVariableTempl[is][im]->SetDirectory(0);
            fMCVariableTempl[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
            fMCVariableTempl[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
        }
    }
    
#else  // __USE_2D_TEMPLATES__
    
    
    // NEW 2D TEMPLATES
    
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        const int      n_q2_bins       = kNQ2Bins[is];
        const int      n_ereco_bins    = kNErecoBins[is];
        const int      n_etrue_bins    = kNEtrueBins[is];
        const double * q2_bin_edges    = util::Q2BinEdges(s);
        const double * ereco_bin_edges = util::ErecoBinEdges(s);
        const int      var_n_q2_bins       = kNXsQ2Bins[is];
        const double * var_q2_bin_edges    = util::XsecQ2BinEdges(s);
        const int      var_n_ereco_bins    = kNXsErecoBins[is];
        const double * var_ereco_bin_edges = util::XsecErecoBinEdges(s);
        
        //const double * etrue_bin_edges = util::EtrueBinEdges(s);
        
        
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            const int      var_n_etrue_bins    = util::NumVariableMCEtrueBins (s,m);
            const double * var_etrue_bin_edges = util::VariableMCEtrueBinEdges(s,m);
            fMCTemplOri2D[is][im] = new TH3D(Form("nsk_ori_%d_%d_2d",is,im), util::AsString(m).c_str(),
                                             n_q2_bins, q2_bin_edges, n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
            fMCTemplNom2D[is][im] = new TH3D(Form("nsk_nom_%d_%d",is,im), util::AsString(m).c_str(),
                                             n_q2_bins, q2_bin_edges, n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
            fMCTemplRPA2D[is][im] = new TH3D(Form("nsk_rpa_%d_%d",is,im), util::AsString(m).c_str(),
                                             n_q2_bins, q2_bin_edges, n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
            fMCTemplPOT2D[is][im] = new TH3D(Form("nsk_pot_%d_%d",is,im), util::AsString(m).c_str(),
                                             n_q2_bins, q2_bin_edges, n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
            
            fMCTemplOri2D[is][im]->SetDirectory(0);
            fMCTemplOri2D[is][im]->GetXaxis()->SetTitle("Q^2 (GeV)");
            fMCTemplOri2D[is][im]->GetYaxis()->SetTitle("E_{reco} (GeV)");
            fMCTemplOri2D[is][im]->GetZaxis()->SetTitle("E_{true} (GeV)");
            fMCTemplNom2D[is][im]->SetDirectory(0);
            fMCTemplNom2D[is][im]->GetXaxis()->SetTitle("Q^2 (GeV)");
            fMCTemplNom2D[is][im]->GetYaxis()->SetTitle("E_{reco} (GeV)");
            fMCTemplNom2D[is][im]->GetZaxis()->SetTitle("E_{true} (GeV)");
            fMCTemplRPA2D[is][im]->SetDirectory(0);
            fMCTemplRPA2D[is][im]->GetXaxis()->SetTitle("Q^2 (GeV)");
            fMCTemplRPA2D[is][im]->GetYaxis()->SetTitle("E_{reco} (GeV)");
            fMCTemplRPA2D[is][im]->GetZaxis()->SetTitle("E_{true} (GeV)");
            fMCTemplPOT2D[is][im]->SetDirectory(0);
            fMCTemplPOT2D[is][im]->GetXaxis()->SetTitle("Q^2 (GeV)");
            fMCTemplPOT2D[is][im]->GetYaxis()->SetTitle("E_{reco} (GeV)");
            fMCTemplPOT2D[is][im]->GetZaxis()->SetTitle("E_{true} (GeV)");
            
            fMCVariableTempl2D[is][im] = new TH3D(Form("nsk_ori_var_%d_%d",is,im), util::AsString(m).c_str(),
                                                  var_n_q2_bins, var_q2_bin_edges, var_n_ereco_bins, var_ereco_bin_edges, var_n_etrue_bins,  var_etrue_bin_edges);
            fMCVariableTempl2D[is][im]->SetDirectory(0);
            fMCVariableTempl2D[is][im]->GetXaxis()->SetTitle("Q^2 (GeV)");
            fMCVariableTempl2D[is][im]->GetYaxis()->SetTitle("E_{reco} (GeV)");
            fMCVariableTempl2D[is][im]->GetZaxis()->SetTitle("E_{true} (GeV)");
        }
    }
    //
#endif // __USE_2D_TEMPLATES__
    
}
//________________________________________________________________________________________
void SuperkMC::FillTemplateHistogram(SuperkSample_t s, SuperkMode_t m)
{
#ifdef __show_init_mesg_debug__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate >> Creating MC template for: "
    << util::AsString(s) << " / " << util::AsString(m);
#endif
    
    //filled in FillSpecialTemplates()
    if(util::IsOscNuMu(m) || util::IsOscNuMuBar(m))
        return;
    
#ifdef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
    if (util::IsOscNuEBar(m)) return;// with that flag, these are filled in FillSpecial
#endif
    
    int sidx = arr_index->SuperkSampleIdx(s);
    int midx = arr_index->SuperkModeIdx(s,m);
    
    TTree * mctree = this->MCTree(s,m);
    double  mcpot  = this->MCPOT (s,m);
    
    if(!mctree) {
        LOG("SuperkMC",pERROR) << "Exp/SKMC/FillTemplate >> Couldn't find event tree for template: "
        << util::AsString(s) << "/" << util::AsString(m);
        util::Die(1); return;
    }
    
    double pot_wght = this->AnalyzedSamplePOT(s) / mcpot;
    
#ifdef POTWEIGHTSAREONE
    pot_wght = 1;
#endif
    
#ifdef __show_init_mesg_debug__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate >>      "
    << "Weight (data POT / MC POT) = " << pot_wght;
#endif
    if(pot_wght < 0) {
        LOG("SuperkMC",pERROR) << "Exp/SKMC/FillTemplate >> Weight (data POT / MC POT) for "
        << util::AsString(s) << "/" << util::AsString(m)
        << " is " << pot_wght;
        util::Die(1); return;
    }
    
    const char * erec_branch          = "erec"; // MeV/c
    const char * etrue_branch         = "pnu";  // GeV/c
    const char * mode_branch          = "mode";
    const char * flux_wght_branch     = "wgtflx";
    const char * berpa_wght_branch    = "wgtBeRPA_nom";
    const char * sf_rfg_wght_branch   = "wgtsfrfg";
    const char * cc_coh_wght_branch   = "wght_CCcoh";
    const char * wall_branch = "fqwall";
    //const char * towall_branch = "fqtowall";

    
    //~ const char * flux_wght_branch = "wgt";
    // // NEW 2D TEMPLATES
    // const char * numnu_branch = "numnu";
    // const char * nring_branch = "nring";
    // const char * amomm_branch = "amomm"; // MeV/c
    // const char * amome_branch = "amome"; // MeV/c
    // const char * dir_branch   = "dir";
    // const char * dirnu_branch = "dirnu";

    // fiTQun branches have different names 
    const char * amomm_branch = "fqmomm"; // MeV/c
    const char * amome_branch = "fqmome"; // MeV/c
    const char * dir_branch   = "fq1rdir";
    const char * dirnu_branch = "dirnu";

    // //
    
    // SetBranchAddress?
    TBranch * branchWall = mctree->GetBranch(wall_branch);

    TBranch * branchEreco   = mctree->GetBranch(erec_branch      );
    TBranch * branchEtrue   = mctree->GetBranch(etrue_branch     );
    TBranch * branchGenCode = mctree->GetBranch(mode_branch      );
    TBranch * branchWghtFlx = mctree->GetBranch(flux_wght_branch );
    // // NEW 2D TEMPLATES
    // TBranch * branchNumnu   = mctree->GetBranch(numnu_branch      );
    // TBranch * branchNring   = mctree->GetBranch(nring_branch      );
    TBranch * branchAmomm   = mctree->GetBranch(amomm_branch      );
    TBranch * branchAmome   = mctree->GetBranch(amome_branch      );
    TBranch * branchDir     = mctree->GetBranch(dir_branch        );
    TBranch * branchDirnu   = mctree->GetBranch(dirnu_branch      );
    // //

    branchEreco    -> SetAddress(&fBrEreco);
    branchEtrue    -> SetAddress(&fBrPtrue);
    branchGenCode  -> SetAddress(&fBrGenCode);
    branchWghtFlx  -> SetAddress(&fBrWghtFlx);

    branchWall->SetAddress(&fBrWall);

    // // NEW 2D TEMPLATES
    // branchNumnu    -> SetAddress(&fBrNumnu);
    // branchNring    -> SetAddress(&fBrNring);
    branchAmomm    -> SetAddress(&fBrAmomm);
    branchAmome    -> SetAddress(&fBrAmome);
    branchDir      -> SetAddress(&fBrDir);
    branchDirnu    -> SetAddress(&fBrDirNu);
    // //

    mctree->SetBranchStatus("*",0);
    mctree->SetBranchStatus(wall_branch, 1);
    mctree->SetBranchStatus(erec_branch      , 1 );
    mctree->SetBranchStatus(etrue_branch     , 1 );
    mctree->SetBranchStatus(mode_branch      , 1 );
    mctree->SetBranchStatus(flux_wght_branch , 1 );
    // // NEW 2D TEMPLATES
    // mctree->SetBranchStatus(numnu_branch     , 1 );
    // mctree->SetBranchStatus(nring_branch     , 1 );
    mctree->SetBranchStatus(amomm_branch     , 1 );
    mctree->SetBranchStatus(amome_branch     , 1 );
    mctree->SetBranchStatus(dir_branch       , 1 );
    mctree->SetBranchStatus(dirnu_branch     , 1 );
    // //


    
#ifdef USE_SFRFG_WEIGHT
    TBranch * branchSfRfg   = mctree->GetBranch(sf_rfg_wght_branch);
    mctree->SetBranchStatus(sf_rfg_wght_branch , 1 );
    branchSfRfg->SetAddress(&fBrWghtSFRFG);
#endif
#ifdef USE_RPA_WEIGHTS
    // Don't have the rpa_nrel_wght_branch in the MultiRing sample
    TBranch * branchRPAnr   = mctree->GetBranch(rpa_nrel_wght_branch);
    if (branchRPAnr){
      mctree->SetBranchStatus(rpa_nrel_wght_branch , 1 );
      branchRPAnr->SetAddress(&fBrWghtRPAnonrel);
    }    
    TBranch * branchRPAr   = mctree->GetBranch(rpa_rel_wght_branch);
    mctree->SetBranchStatus(rpa_rel_wght_branch , 1 );
    branchRPAr->SetAddress(&fBrWghtRPArel);
#endif
#ifdef USE_BeRPA_WEIGHTS
    TBranch * branchBeRPA = mctree->GetBranch(berpa_wght_branch);
    if (branchBeRPA){
      mctree->SetBranchStatus(berpa_wght_branch , 1 );
      branchBeRPA->SetAddress(&fBrWghtBeRPA);
    } 
#endif    
#ifdef USE_CCcoh_WEIGHTS
    TBranch * branchCCcoh   = mctree->GetBranch(cc_coh_wght_branch);
    mctree->SetBranchStatus(cc_coh_wght_branch , 1 );
    branchCCcoh->SetAddress(&fBrWghtCCcoh);
#endif
#ifdef USE_OscBr_WEIGHTS
    if (s == kSKsFCFVMultiRingELike ){
    TBranch * branchOscWgt   = mctree->GetBranch(osc_wgt_branch);
    mctree->SetBranchStatus(osc_wgt_branch, 1 );
    branchOscWgt->SetAddress(&fBrWghtOsc);
    }
#endif
  


    bool is_fake_data_active = util::UtlMakeFakeData::Instance()->IsActive();
  
    double ntotevents = 0.;
    
    for(int iev = 0; iev < mctree->GetEntries(); iev++)
    {
        mctree->GetEntry(iev);
        
        if(fBrEreco<0) continue;
	if(fBrWall<200) continue; //use old cut
        
        // use current entry for mode m template?
        bool use_it = util::UseEventToFillInTemplate(m,util::SuperkMode2SKMCPdg(m),fBrGenCode);
        if(!use_it) continue;
        
        // Get reconstructed energy
        double ereco_fill = fBrEreco/1000.; // GeV
        
        
        // NEW 2D TEMPLATES
        
#if defined(__ERECO_VS_Q2_TEMPLATES__) // Use Ereco Vs Q2 parametrization
        
        // Get momentum and energy of most energetic lepton (only 1Ring samples)
        
        double pmax_lep = 0.; // momentum of most energetic lepton
        double mass_lep = 0.;
        if( (s == kSKsFCFV1RingELike) || (s == kSKsRHCFCFV1RingELike) ){ // e-like ring
            mass_lep = kElectronMassGeV; // GeV
            pmax_lep = fBrAmome[0] / 1000.; // GeV
        }
        else if( (s == kSKsFCFV1RingMuLike) || (s == kSKsRHCFCFV1RingMuLike) ){ // mu-like ring
            mass_lep = kMuonMassGeV; // GeV
            pmax_lep = fBrAmomm[0] / 1000.; // GeV
        }
        else{
            LOG("SuperkMC",pERROR) ;
            LOG("SuperkMC",pERROR) << "SuperkMC::FillTemplateHistogram" ;
            LOG("SuperkMC",pERROR) << "Using sample " << util::AsString(s) ;
            LOG("SuperkMC",pERROR) << "Samples other than: " ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsFCFV1RingMuLike) ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsFCFV1RingELike) ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsRHCFCFV1RingMuLike);
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsRHCFCFV1RingMuLike);
            LOG("SuperkMC",pERROR) << "are not implemented" ;
            util::Die(1); return;
        }
        
        double mass_nucl = 0.;
        if( (s == kSKsRHCFCFV1RingELike) || (s == kSKsRHCFCFV1RingMuLike) ){ // RHC
            mass_nucl = kProtonMassGeV; // nubar interacts with proton
        }
        else if( (s == kSKsFCFV1RingELike) || (s == kSKsFCFV1RingMuLike) ){ // FHC
            mass_nucl = kNeutronMassGeV; // nubar interacts with neutron
        }
        else{
            LOG("SuperkMC",pERROR) ;
            LOG("SuperkMC",pERROR) << "SuperkMC::FillTemplateHistogram" ;
            LOG("SuperkMC",pERROR) << "Using sample " << util::AsString(s);
            LOG("SuperkMC",pERROR) << "Samples other than: ";
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsFCFV1RingMuLike);
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsFCFV1RingELike);
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsRHCFCFV1RingMuLike);
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsRHCFCFV1RingMuLike);
            LOG("SuperkMC",pERROR) << "are not implemented";
            util::Die(1); return;
        }
        
        // Get lepton energy
        
        double emax_lep = sqrt( pmax_lep*pmax_lep + mass_lep*mass_lep );
        
        // Get lepton direction
        // Costheta is the scalar product between the neutrino beam and lepton directions.
        double costheta = -999;
        if (util::Is1ReSample(s) || util::IsMultiReSample(s)) {
          double costheta = kBeamDir[0] * fBrDir[0][1][0] + kBeamDir[1] * fBrDir[0][1][1] + kBeamDir[2] * fBrDir[0][1][2];
        }
        if util::Is1RMuSample(s) {
            double costheta = kBeamDir[0] * fBrDir[0][2][0] + kBeamDir[1] * fBrDir[0][2][1] + kBeamDir[2] * fBrDir[0][2][2];    
        }
        if( costheta<-1.00001 || costheta>+1.00001 ){
            LOG("SuperkMC",pERROR);
            LOG("SuperkMC",pERROR) << "SuperkMC::FillTemplateHistogram";
            LOG("SuperkMC",pERROR) << "costheta = " << costheta << " --> outside allowed range";
            LOG("SuperkMC",pERROR);
            util::Die(1); return;
        }
        
        double sintheta = 1. - costheta*costheta;
        if( sintheta<-1.00001 || sintheta>+1.00001 ){
            LOG("SuperkMC",pERROR);
            LOG("SuperkMC",pERROR) << "SuperkMC::FillTemplateHistogram";
            LOG("SuperkMC",pERROR) << "sintheta = " << sintheta << " --> outside allowed range";
            LOG("SuperkMC",pERROR) << "costheta = " << costheta;
            LOG("SuperkMC",pERROR);
            util::Die(1); return;
        }
        
        double LepTheta        = TMath::ACos(costheta); // radians
        double LepThetaDegrees = LepTheta / TMath::Pi() * 180.; // degrees
        if(LepTheta<0.){
            LOG("SuperkMC",pERROR);
            LOG("SuperkMC",pERROR) << "Theta lepton is negative: " << LepThetaDegrees << " degrees";
            LOG("SuperkMC",pERROR);
            util::Die(1); return;
        }
        
        
        double ereco_nu = ereco_fill; // neutrino reconstructed energy (same for momentum)
        
        //double Q2 = this->CalcQ2(ereco_nu,emax_lep,pmax_lep,costheta,mass_lep); // Q^2
        //double Q2 = this->CalcQd2(ereco_nu,emax_lep,pmax_lep,costheta,mass_lep); // Qd^2
        //double Q2 = this->CalcEhad(ereco_nu,emax_lep);
        //double Q2 = this->CalcInelasticity(ereco_nu,emax_lep);
        //double Q2 = this->CalcBjorkenScale(ereco_nu,emax_lep,pmax_lep,costheta,mass_lep,mass_nucl);
        //double Q2 = this->CalcWhadr(ereco_nu,emax_lep,pmax_lep,costheta,mass_lep,mass_nucl);
        double Q2 = this->CalcXSecTerm(ereco_nu,emax_lep,pmax_lep,costheta,mass_lep,mass_nucl);
        
        double q2_fill = Q2;
        
#ifdef INCLUDEEVENTSNEGATIVEQ2
        if(q2_fill < 0. && q2_fill> -0.01){
            q2_fill = 0.01;
        }
#endif
#ifdef INCLUDEEVENTSOVERFLOWQ2
        int nbins = util::NumQ2Bins(s);
        double MaxQ2_up   = util::Q2BinEdges(s,nbins);
        double MaxQ2_low  = util::Q2BinEdges(s,nbins-1);
        double MaxQ2_cent = (MaxQ2_up + MaxQ2_low)/2.;
        q2_fill = TMath::Min(MaxQ2_cent, q2_fill);
        
        if(q2_fill > MaxQ2_up){
            LOG("SuperkMC",pDEBUG) << q2_fill;
            util::Die(1); return;
        }
#endif
        // Will need to be updated if we use the fitqun 
        if( q2_fill < util::Q2BinEdges(s,0) ){
            LOG("SuperkMC",pDEBUG) << "Q2 < Q2BinEdge(s,0)";
            LOG("SuperkMC",pDEBUG) << "nu direction: " << kBeamDir[0] << " " << kBeamDir[1] << " " << kBeamDir[2] << " ";
            LOG("SuperkMC",pDEBUG) << "lepton direction: " << fBrDir[0][0] << " " << fBrDir[0][1] << " " << fBrDir[0][2] << " ";
            LOG("SuperkMC",pDEBUG) << "theta = " << LepThetaDegrees << " ";
            LOG("SuperkMC",pDEBUG) << "costheta = " << costheta << " ";
            //LOG("SuperkMC",pDEBUG) << "sintheta = " << sintheta << " ";
            LOG("SuperkMC",pDEBUG) << "emax_lep = " << emax_lep << " ";
            LOG("SuperkMC",pDEBUG) << "ereco_nu = " << ereco_nu << " ";
            LOG("SuperkMC",pDEBUG) << "Q2 = " << Q2 << " GeV";
            util::Die(1); return;
        }
        if( q2_fill > util::Q2BinEdges(s,nbins) ){
            LOG("SuperkMC",pDEBUG) << "Q2 > Q2BinEdge(s,0)";
            LOG("SuperkMC",pDEBUG) << "costheta = " << costheta << " ";
            LOG("SuperkMC",pDEBUG) << "sintheta = " << sintheta << " ";
            LOG("SuperkMC",pDEBUG) << "emax_lep = " << emax_lep << " ";
            LOG("SuperkMC",pDEBUG) << "ereco_nu = " << ereco_nu << " ";
            LOG("SuperkMC",pDEBUG) << "Q2 = " << Q2 << " GeV";
            util::Die(1); return;
        }
        
        
#elif defined(__ERECO_VS_THETA_TEMPLATES__) // Use Ereco Vs Theta parametrization

	double q2_fill = CalcThetaDegrees(fBrDir,s);

#ifdef INCLUDEEVENTSNEGATIVEQ2
        if(q2_fill < 0.){
            q2_fill = util::Q2BinEdges(s,0)+conv::kASmallNum;
        }
#endif
#ifdef INCLUDEEVENTSOVERFLOWQ2
        int nbins = util::NumQ2Bins(s);
        double MaxQ2_up   = util::Q2BinEdges(s,nbins);
        double MaxQ2_low  = util::Q2BinEdges(s,nbins-1);
        double MaxQ2_cent = (MaxQ2_up + MaxQ2_low)/2.;
        q2_fill = TMath::Min(MaxQ2_cent, q2_fill);
        
        if(q2_fill > MaxQ2_up){
            LOG("SuperkMC",pDEBUG) << q2_fill;
            util::Die(1); return;
        }
#endif
        
        if( q2_fill < util::Q2BinEdges(s,0) ){
          LOG("SuperkMC",pDEBUG) << "q2_fill = " << q2_fill;
	  util::Die(1); return;
        }
        
        if( q2_fill > util::Q2BinEdges(s,nbins) ){
          LOG("SuperkMC",pDEBUG) << "q2_fill = " << q2_fill;
	  util::Die(1); return;
        }
        
#else
        LOG("SuperkMC",pERROR) << "SuperkMC::FillTemplateHistogram()";
        LOG("SuperkMC",pERROR) << "Need to choose between Ereco Vs Q2 or Ereco Vs Theta parametrization";
        util::Die(1); return;
#endif //__ERECO_VS_Q2_TEMPLATES__ or __ERECO_VS_THETA_TEMPLATES__
        
        double Etrue = fBrPtrue[0];
        
        //if(util::IsRHCSample(s))
#ifdef USE_FLUX_EXT_WEIGHT
        fBrWghtFlx = fNuFlux->SuperkFluxExtTuningWght(s,m,Etrue);
#endif
        
#ifdef FLUXWEIGHTSAREONE
        fBrWghtFlx = 1.;
#endif
        
#ifdef INCLUDEEVENTSOVERTHIRTYGEV
        ereco_fill = TMath::Min(29.9, ereco_fill);
#endif
        
        
        //q2_fill = ereco_fill; // NEW 2D TEMPLATES (TODO) // be careful to overflows!!!
        
        
        
#ifndef __USE_2D_TEMPLATES__
        this -> Template(s,m,kSKMCtOri) -> Fill (ereco_fill, Etrue);
        fMCVariableTempl[sidx][midx]    -> Fill (ereco_fill, Etrue);
        ntotevents += 1.;
#else
        // NEW 2D TEMPLATES
        this -> Template2D(s,m,kSKMCtOri) -> Fill (q2_fill, ereco_fill, Etrue);
        fMCVariableTempl2D[sidx][midx]    -> Fill (q2_fill, ereco_fill, Etrue);
        ntotevents += 1.;
        //
#endif
        
        double pot_flx_weight = pot_wght * fBrWghtFlx;
        double NC1GammaWght = 1.;       
#ifdef USE_NC1gamma_WEIGHTS
        if(util::IsNC1Gamma(TMath::Abs(fBrGenCode))){
            NC1GammaWght = 2.;
        }
#endif

#ifdef CANCEL_CCCoh_WEIGHTS
        // In 2017 analysis SuperK mtuples included NC1Gamma and 2p2h weights in CCCoh weight
        // creating problems with double-weighting this events. As CCCoh does not apply to these
        // events, we override the mtuple weight and apply the weighting manually
        if(util::IsNC1Gamma(TMath::Abs(fBrGenCode)) || util::Is2p2h(TMath::Abs(fBrGenCode))){
            fBrWghtCCcoh = 1.;
        }
#endif

        double fake_data_weight = 1.;
      
      if(is_fake_data_active){
        util::UtlMakeFakeData::Instance()->SetEventInfo(fBrPtrue,fBrDirNu,fBrDir,Etrue,m,ereco_fill,s,fBrAmome,fBrAmomm);
        fake_data_weight = util::UtlMakeFakeData::Instance()->GetWeight();
        ereco_fill = util::UtlMakeFakeData::Instance()->GetEreco();
#ifdef FakeDataWeightDebug
        LOG("SuperkMC",pDEBUG) << "Fake data reweight: " << util::UtlMakeFakeData::Instance()->GetReweightStyle() <<  " " << fake_data_weight
                               << " Ereco: " << ereco_fill << " Sample: " << util::AsStringCompact(s)
                               << " mode: " << util::AsStringCompact(m);
#endif
      }
      
#ifndef __USE_2D_TEMPLATES__
        this -> Template(s,m,kSKMCtPOT) -> Fill (ereco_fill, Etrue, pot_flx_weight * fake_data_weight);
#else
        //NEW 2D TEMPLATES
        this -> Template2D(s,m,kSKMCtPOT) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fake_data_weight);
        //
#endif
        
        // Actually fill the template
        
#ifndef __USE_2D_TEMPLATES__
#ifdef USE_RPA_WEIGHTS || USE_OscBr_WEIGHTS 
#  ifdef USE_RPA_WEIGHTS && !defined(USE_OscBr_WEIGHTS)        
	if(branchRPAnr){
	  this -> Template(s,m,kSKMCtNom) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtRPAnonrel * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fake_data_weight);
	} else {
	  this -> Template(s,m,kSKMCtNom) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fake_data_weight);
	}
        this -> Template(s,m,kSKMCtRPA) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtRPArel * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fake_data_weight);

#  elif !defined(USE_RPA_WEIGHTS) && defined(USE_OscBr_WEIGHTS)
        if (s == kSKsFCFVMultiRingELike ){
        this -> Template(s,m,kSKMCtNom) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);
        this -> Template(s,m,kSKMCtRPA) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);
        }
#  elif defined(USE_RPA_WEIGHTS) && defined(USE_OscBr_WEIGHTS)
        if (s == kSKsFCFVMultiRingELike ){
        if(branchRPAnr){
	  this -> Template(s,m,kSKMCtNom) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtRPAnonrel * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);
	} else {
	  this -> Template(s,m,kSKMCtNom) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);
	}
        this -> Template(s,m,kSKMCtRPA) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtRPArel * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);
        }
#  endif
#else
        this -> Template(s,m,kSKMCtNom) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fake_data_weight);
        this -> Template(s,m,kSKMCtRPA) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fake_data_weight);
#endif
#else
        // NEW 2D TEMPLATES
#if defined(USE_RPA_WEIGHTS) || defined(USE_OscBr_WEIGHTS) || defined(USE_BeRPA_WEIGHTS) 
#  if defined(USE_RPA_WEIGHTS) && !defined(USE_OscBr_WEIGHTS)       
	if(branchRPAnr){
	  this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtRPAnonrel * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fake_data_weight);
	} else {
          // Adding the line below gives exact agreement with the output from the p-theta code, wonder about how this is treated
          // further on in the code see that the ratio is used?
	  //this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtRPArel * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght);
          this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fake_data_weight);
	}
#  elif !defined(USE_RPA_WEIGHTS) && defined(USE_OscBr_WEIGHTS)
          this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);
          this -> Template2D(s,m,kSKMCtRPA) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);

#  elif defined(USE_RPA_WEIGHTS) && defined(USE_OscBr_WEIGHTS)
        LOG("SuperkMC",pINFO) << "Applying RPA and Oscillation branch weights" ;
        if(branchRPAnr){
	  this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtRPAnonrel * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);
	} else {
	  this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);
	}
        this -> Template2D(s,m,kSKMCtRPA) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtRPArel * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);
        }

#  elif defined(USE_BeRPA_WEIGHTS) && !defined(USE_OscBr_WEIGHTS)       
   #ifdef FORCE_SFRFG_WEIGHT_1
      fBrWghtSFRFG = 1.0; // Can't just turn off USE_SFRFG_WEIGHT, so force to 1
   #endif
   #ifdef FORCE_BeRPA_WEIGHTS_1
      fBrWghtBeRPA = 1.0; // Can't just turn off USE_BeRPA_WEIGHTS, so force to 1
   #endif
        this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtBeRPA * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fake_data_weight);
        this -> Template2D(s,m,kSKMCtRPA) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtBeRPA * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fake_data_weight);
        
#  elif !defined(USE_BeRPA_WEIGHTS) && defined(USE_OscBr_WEIGHTS)
        this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);
        this -> Template2D(s,m,kSKMCtRPA) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);
        
#  elif defined(USE_BeRPA_WEIGHTS) && defined(USE_OscBr_WEIGHTS)
   #ifdef FORCE_SFRFG_WEIGHT_1
      fBrWghtSFRFG = 1.0; // Can't just turn off USE_SFRFG_WEIGHT, so force to 1
   #endif
   #ifdef FORCE_BeRPA_WEIGHTS_1
      fBrWghtBeRPA = 1.0; // Can't just turn off USE_BeRPA_WEIGHTS, so force to 1
   #endif
        this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtBeRPA * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);
        this -> Template2D(s,m,kSKMCtRPA) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtBeRPA * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtOsc * fake_data_weight);

#  else        
        this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fake_data_weight);
        this -> Template2D(s,m,kSKMCtRPA) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fake_data_weight);
#  endif
#endif //

        // if (fBrWghtBeRPA > conv::kASmallNum){
        //   LOG("SuperkMC", pINFO) << "Value of BeRPA weight"  << fBrWghtBeRPA;
        // }
        
#endif // __USE_2D_TEMPLATES__
        
        
        
#ifndef __USE_2D_TEMPLATES__
        
#ifdef __show_init_mesg_debug__
        if(Etrue <= 0 || Etrue >= this->Template(s,m,kSKMCtOri)->GetYaxis()->GetBinUpEdge(this->Template(s,m,kSKMCtOri)->GetNbinsY()))
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate >> Etrue is out of template range " << AsString(s) << " " << AsString(m) << " " << Etrue
            << " (Ereco = " << fBrEreco/1000 << ")";
#endif
#ifdef __show_init_mesg_debug2__
        if(fBrEreco/1000 <= 0)
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << fBrEreco/1000
            << " (Etrue = " << Etrue << ")";
        //#endif __show_init_mesg_debug2__
        if(fBrEreco/1000 >= this->Template(s,m,kSKMCtOri)->GetXaxis()->GetBinUpEdge(this->Template(s,m,kSKMCtOri)->GetNbinsX()))
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << fBrEreco/1000
            << " (Etrue = " << Etrue << ")";
#endif
        
#else
        
        // NEW 2D TEMPLATES
        
#ifdef __show_init_mesg_debug__
        if(Etrue <= 0 || Etrue >= this->Template2D(s,m,kSKMCtOri)->GetZaxis()->GetBinUpEdge(this->Template2D(s,m,kSKMCtOri)->GetNbinsZ())){
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate >> Etrue is out of template range " << AsString(s) << " " << AsString(m) << " " << Etrue
            << " (Ereco = " << fBrEreco/1000 << ")"
            << " (Q2 = " << q2_fill << ")";
        }
#endif
#ifdef __show_init_mesg_debug2__
        if(fBrEreco/1000 <= 0)
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << fBrEreco/1000
            << " (Etrue = " << Etrue << ")";
        //#endif __show_init_mesg_debug2__
        if(fBrEreco/1000 >= this->Template2D(s,m,kSKMCtOri)->GetYaxis()->GetBinUpEdge(this->Template2D(s,m,kSKMCtOri)->GetNbinsY()))
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << fBrEreco/1000
            << " (Etrue = " << Etrue << ")"
            << " (Q2 = " << q2_fill << ")";
        
        if(q2_fill <= 0)
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate >> Q2 is out of template range " << AsString(s) << " " << AsString(m) << " " << q2_fill
            << " (Ereco = " << fBrEreco/1000 << ")"
            << " (Etrue = " << Etrue << ")";
        //#endif __show_init_mesg_debug2__
        if(q2_fill >= this->Template2D(s,m,kSKMCtOri)->GetXaxis()->GetBinUpEdge(this->Template2D(s,m,kSKMCtOri)->GetNbinsX()))
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate >> Q2 is out of template range " << AsString(s) << " " << AsString(m) << " " << q2_fill
            << " (Ereco = " << fBrEreco/1000 << ")"
            << " (Etrue = " << Etrue << ")";
        //
#endif// __show_init_mesg_debug__
        
#endif
        
    }//iev loop over events
    
#ifdef __show_init_mesg_debug__
    
#ifndef __USE_2D_TEMPLATES__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate >>      "
    << "Number of template entries = "
    << this -> Template(s,m,kSKMCtNom) -> GetEntries();
#else
    // NEW 2D TEMPLATES
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate2D >>      "
    << "Number of template entries = "
    << this -> Template2D(s,m,kSKMCtNom) -> GetEntries();
    //
#endif

    LOG("SuperKMC", pDEBUG) << util::AsString(s) <<  " " << util::AsString(m)  << "Total events in MC: " << ntotevents;
    
#endif
    
}

//________________________________________________________________________________________
void SuperkMC::FillTemplateHistogramSmeared(SuperkSample_t s, SuperkMode_t m)
{
#ifdef __show_init_mesg_debug__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate >> Creating MC template for: "
    << util::AsString(s) << " / " << util::AsString(m);
#endif
    
    //filled in FillSpecialTemplates()
    if(util::IsOscNuMu(m) || util::IsOscNuMuBar(m))
        return;
    
#ifdef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
    if (util::IsOscNuEBar(m)) return;// with that flag, these are filled in FillSpecial
#endif
    
    int sidx = arr_index->SuperkSampleIdx(s);
    int midx = arr_index->SuperkModeIdx(s,m);
    
    TTree * mctree = this->MCTree(s,m);
    double  mcpot  = this->MCPOT (s,m);
    
    if(!mctree) {
        LOG("SuperkMC",pERROR) << "Exp/SKMC/FillTemplate >> Couldn't find event tree for template: "
        << util::AsString(s) << "/" << util::AsString(m);
        util::Die(1); return;
    }
    
    double pot_wght = this->AnalyzedSamplePOT(s) / mcpot;
    
#ifdef POTWEIGHTSAREONE
    pot_wght = 1;
#endif
    
#ifdef __show_init_mesg_debug__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate >>      "
    << "Weight (data POT / MC POT) = " << pot_wght;
#endif
    if(pot_wght < 0) {
        LOG("SuperkMC",pERROR) << "Exp/SKMC/FillTemplate >> Weight (data POT / MC POT) for "
        << util::AsString(s) << "/" << util::AsString(m)
        << " is " << pot_wght;
        util::Die(1); return;
    }
    
    const char * erec_branch          = "erec"; // MeV/c
    const char * etrue_branch         = "pnu";  // GeV/c
    const char * mode_branch          = "mode";
    const char * flux_wght_branch     = "wgtflx";
    const char * berpa_wght_branch    = "wgtBeRPA_nom";
    const char * sf_rfg_wght_branch   = "wgtsfrfg";
    const char * cc_coh_wght_branch   = "wght_CCcoh";
    
    //~ const char * flux_wght_branch = "wgt";
    // // NEW 2D TEMPLATES
    // const char * numnu_branch = "numnu";
    // const char * nring_branch = "nring";
    const char * amomm_branch = "amomm"; // MeV/c
    const char * amome_branch = "amome"; // MeV/c
    const char * dir_branch   = "dir";
    const char * dirnu_branch = "dirnu";
    // //
    
    TBranch * branchEreco   = mctree->GetBranch(erec_branch      );
    TBranch * branchEtrue   = mctree->GetBranch(etrue_branch     );
    TBranch * branchGenCode = mctree->GetBranch(mode_branch      );
    TBranch * branchWghtFlx = mctree->GetBranch(flux_wght_branch );
    // // NEW 2D TEMPLATES
    // TBranch * branchNumnu   = mctree->GetBranch(numnu_branch      );
    // TBranch * branchNring   = mctree->GetBranch(nring_branch      );
    TBranch * branchAmomm   = mctree->GetBranch(amomm_branch      );
    TBranch * branchAmome   = mctree->GetBranch(amome_branch      );
    TBranch * branchDir     = mctree->GetBranch(dir_branch        );
    TBranch * branchDirnu   = mctree->GetBranch(dirnu_branch      );
    // //
    
    branchEreco    -> SetAddress(&fBrEreco);
    branchEtrue    -> SetAddress(&fBrPtrue);
    branchGenCode  -> SetAddress(&fBrGenCode);
    branchWghtFlx  -> SetAddress(&fBrWghtFlx);
    // // NEW 2D TEMPLATES
    // branchNumnu    -> SetAddress(&fBrNumnu);
    // branchNring    -> SetAddress(&fBrNring);
    branchAmomm    -> SetAddress(&fBrAmomm);
    branchAmome    -> SetAddress(&fBrAmome);
    branchDir      -> SetAddress(&fBrDir);
    branchDirnu    -> SetAddress(&fBrDirNu);
    // //
    
    mctree->SetBranchStatus("*",0);
    mctree->SetBranchStatus(erec_branch      , 1 );
    mctree->SetBranchStatus(etrue_branch     , 1 );
    mctree->SetBranchStatus(mode_branch      , 1 );
    mctree->SetBranchStatus(flux_wght_branch , 1 );
    // // NEW 2D TEMPLATES
    // mctree->SetBranchStatus(numnu_branch     , 1 );
    // mctree->SetBranchStatus(nring_branch     , 1 );
    mctree->SetBranchStatus(amomm_branch     , 1 );
    mctree->SetBranchStatus(amome_branch     , 1 );
    mctree->SetBranchStatus(dir_branch       , 1 );
    mctree->SetBranchStatus(dirnu_branch     , 1 );
    // //
    
#ifdef USE_SFRFG_WEIGHT
    TBranch * branchSfRfg   = mctree->GetBranch(sf_rfg_wght_branch);
    mctree->SetBranchStatus(sf_rfg_wght_branch , 1 );
    branchSfRfg->SetAddress(&fBrWghtSFRFG);
#endif
#ifdef USE_RPA_WEIGHTS
    // Branch may not (doesn't currently) exist for MultiRing samples
    TBranch * branchRPAnr   = mctree->GetBranch(rpa_nrel_wght_branch);
    if (branchRPAnr){
      mctree->SetBranchStatus(rpa_nrel_wght_branch , 1 );
      branchRPAnr->SetAddress(&fBrWghtRPAnonrel);
    }    
    TBranch * branchRPAr   = mctree->GetBranch(rpa_rel_wght_branch);
    mctree->SetBranchStatus(rpa_rel_wght_branch , 1 );
    branchRPAr->SetAddress(&fBrWghtRPArel);
#endif
#ifdef USE_BeRPA_WEIGHTS
    TBranch * branchBeRPA = mctree->GetBranch(berpa_wght_branch);
    if (branchBeRPA){
      mctree->SetBranchStatus(berpa_wght_branch , 1 );
      branchBeRPA->SetAddress(&fBrWghtBeRPA);
    } 
#endif    
#ifdef USE_CCcoh_WEIGHTS
    TBranch * branchCCcoh   = mctree->GetBranch(cc_coh_wght_branch);
    mctree->SetBranchStatus(cc_coh_wght_branch , 1 );
    branchCCcoh->SetAddress(&fBrWghtCCcoh);
#endif
    

    
    double ntotevents = 0.;
    
    for(int iev = 0; iev < mctree->GetEntries(); iev++)
    {
        mctree->GetEntry(iev);
        
        if(fBrEreco<0) continue;
        
        // use current entry for mode m template?
        bool use_it = util::UseEventToFillInTemplate(m,util::SuperkMode2SKMCPdg(m),fBrGenCode);
        if(!use_it) continue;
        
        //Get True energy
        double Etrue = fBrPtrue[0];

        // Get reconstructed energy
        double ereco_fill = fBrEreco/1000.; // GeV
        
        // NEW 2D TEMPLATES
        
        // Get lepton direction
        // Costheta is the scalar product between the neutrino beam and lepton directions.
        double costheta = -999;
        double dir_nu[3] = {0.669764, -0.742179, 0.024223}; // T2K beam direction
        if (util::Is1ReSample(s) || util::IsMultiReSample(s)) {
          costheta = dir_nu[0] * fBrDir[0][1][0] + dir_nu[1] * fBrDir[0][1][1] + dir_nu[2] * fBrDir[0][1][2];
        }
        if (util::Is1RMuSample(s)) {
          costheta = dir_nu[0] * fBrDir[0][2][0] + dir_nu[1] * fBrDir[0][2][1] + dir_nu[2] * fBrDir[0][2][2];
        }
        if( costheta<-1.00001 || costheta>+1.00001 ){
            LOG("SuperkMC",pERROR) << "SuperkMC::FillTemplateHistogram" ;
            LOG("SuperkMC",pERROR) << "costheta = " << costheta << " --> outside allowed range" ;
            util::Die(1); return;
        }
        
        double sintheta = 1. - costheta*costheta;
        if( sintheta<-1.00001 || sintheta>+1.00001 ){
            LOG("SuperkMC",pERROR) << "SuperkMC::FillTemplateHistogram" ;
            LOG("SuperkMC",pERROR) << "sintheta = " << sintheta << " --> outside allowed range" ;
            LOG("SuperkMC",pERROR) << "costheta = " << costheta ;
            util::Die(1); return;
        }
        
        double LepTheta        = TMath::ACos(costheta); // radiant
        double LepThetaDegrees = LepTheta / TMath::Pi() * 180.; // degrees
        if(LepTheta<0.){
            LOG("SuperkMC",pERROR) << "Theta lepton is negative: " << LepThetaDegrees << " degrees" ;
            util::Die(1); return;
        }
        
        
        //
        //Smear theta based on gaussian around the true value and error from http://arxiv.org/pdf/1502.01550v2.pdf
        //Different muon and electron theta resolution, I took that to mean muon sample and electron sample.
        //
        double theta_error_degrees = 0.;
        
        if( (s == kSKsFCFV1RingELike) || (s == kSKsRHCFCFV1RingELike) ){ // e-like ring
            theta_error_degrees = 1.8;
        }
        else if( (s == kSKsFCFV1RingMuLike) || (s == kSKsRHCFCFV1RingMuLike) ){ // mu-like ring
            theta_error_degrees = 3.;
        }
        
        double LepThetaDegreesSmeared = Rnd::Instance()->RndOther().Gaus(LepThetaDegrees,theta_error_degrees);
        
        //Set all values based on smeared quantitiy
        
        double costhetaSmeared = TMath::Cos(LepTheta);
        
        // Get the true lepton momentum from the tree, position in array depends on int mode.

        int lepton_int =0;
        if(util::Is2p2h(m)){
            lepton_int = 3;
        }
        else{
            lepton_int = 2;
        }
        
        double ptrue_lep = fBrPtrue[lepton_int];
        double pmax_lepSmeared = ptrue_lep;
        
        if( (s == kSKsFCFV1RingELike) || (s == kSKsRHCFCFV1RingELike) ){ // e-like ring
            pmax_lepSmeared = SmearMomentum(ptrue_lep,"electron");
        }
        else if( (s == kSKsFCFV1RingMuLike) || (s == kSKsRHCFCFV1RingMuLike) ){ // mu-like ring
            pmax_lepSmeared = SmearMomentum(ptrue_lep,"muon");
        }


        // Set the correct lepton mass
        
        double mass_lep = 0.;
        if( (s == kSKsFCFV1RingELike) || (s == kSKsRHCFCFV1RingELike) ){ // e-like ring
            mass_lep    = kElectronMassGeV; // GeV
        }
        else if( (s == kSKsFCFV1RingMuLike) || (s == kSKsRHCFCFV1RingMuLike) ){ // mu-like ring
            mass_lep = kMuonMassGeV; // GeV
        }
        else{
            LOG("SuperkMC",pERROR) << "SuperkMC::FillTemplateHistogram" ;
            LOG("SuperkMC",pERROR) << "Using sample " << util::AsString(s);
            LOG("SuperkMC",pERROR) << "Samples other than: " ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsFCFV1RingMuLike) ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsFCFV1RingELike) ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsRHCFCFV1RingMuLike) ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsRHCFCFV1RingMuLike);
            LOG("SuperkMC",pERROR) << "are not implemented" ;
            util::Die(1); return;
        }
        
        /*
         
         Recalculate ERec using new plep from distribution!
         
         */
        
        double ereco_nuSmeared = 0;
        
        if(util::IsRHCSample(s)){
            ereco_nuSmeared = ErecoFunc_NuBar(pmax_lepSmeared,costhetaSmeared,mass_lep);
        }
        else{
            ereco_nuSmeared = ErecoFunc_Nu(pmax_lepSmeared,costhetaSmeared,mass_lep);
        }
        
        if(ereco_fill < 0) ereco_fill = 0;
        if(ereco_nuSmeared < 0) ereco_nuSmeared = 0;

#if defined(__ERECO_VS_Q2_TEMPLATES__) // Use Ereco Vs Q2 parametrization
        
        // Get momentum and energy of most energetic lepton (only 1Ring samples)
        
        double pmax_lep = 0.; // momentum of most energetic lepton
        if( (s == kSKsFCFV1RingELike) || (s == kSKsRHCFCFV1RingELike) ){ // e-like ring
            pmax_lep    = fBrAmome[0] /1000.; // GeV
        }
        else if( (s == kSKsFCFV1RingMuLike) || (s == kSKsRHCFCFV1RingMuLike) ){ // mu-like ring
            pmax_lep = fBrAmomm[0] / 1000.; // GeV
        }
        else{
            LOG("SuperkMC",pERROR) << "SuperkMC::FillTemplateHistogram";
            LOG("SuperkMC",pERROR) << "Using sample " << util::AsString(s) ;
            LOG("SuperkMC",pERROR) << "Samples other than: " ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsFCFV1RingMuLike) ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsFCFV1RingELike) ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsRHCFCFV1RingMuLike) ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsRHCFCFV1RingMuLike) ;
            LOG("SuperkMC",pERROR) << "are not implemented" ;
            util::Die(1); return;
        }
        
        double mass_nucl = 0.;
        if( (s == kSKsRHCFCFV1RingELike) || (s == kSKsRHCFCFV1RingMuLike) ){ // RHC
            mass_nucl = kProtonMassGeV; // nubar interacts with proton
        }
        else if( (s == kSKsFCFV1RingELike) || (s == kSKsFCFV1RingMuLike) ){ // FHC
            mass_nucl = kNeutronMassGeV; // nubar interacts with neutron
        }
        else{
            LOG("SuperkMC",pERROR) << "SuperkMC::FillTemplateHistogram" ;
            LOG("SuperkMC",pERROR) << "Using sample " << util::AsString(s) ;
            LOG("SuperkMC",pERROR) << "Samples other than: " ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsFCFV1RingMuLike) ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsFCFV1RingELike) ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsRHCFCFV1RingMuLike) ;
            LOG("SuperkMC",pERROR) << " - " << util::AsString(kSKsRHCFCFV1RingMuLike) ;
            LOG("SuperkMC",pERROR) << "are not implemented" ;
            util::Die(1); return;
        }
        
        double ereco_nu = ereco_fill; // neutrino reconstructed energy (same for momentum)
        
        //double Q2 = this->CalcQ2(ereco_nuSmeared,emax_lepSmeared,pmax_lepSmeared,costhetaSmeared,mass_lep); // Q^2
        //double Q2 = this->CalcQd2(ereco_nuSmeared,emax_lepSmeared,pmax_lepSmeared,costhetaSmeared,mass_lep); // Qd^2
        //double Q2 = this->CalcEhad(ereco_nuSmeared,emax_lepSmeared);
        //double Q2 = this->CalcInelasticity(ereco_nuSmeared,emax_lepSmeared);
        //double Q2 = this->CalcBjorkenScale(ereco_nuSmeared,emax_lepSmeared,pmax_lepSmeared,costhetaSmeared,mass_lep,mass_nucl);
        //double Q2 = this->CalcWhadr(ereco_nuSmeared,emax_lepSmeared,pmax_lepSmeared,costhetaSmeared,mass_lep,mass_nucl);
        double Q2 = this->CalcXSecTerm(ereco_nuSmeared,emax_lepSmeared,pmax_lepSmeared,costhetaSmeared,mass_lep,mass_nucl);
        
        double q2_fill = Q2;
        
#ifdef INCLUDEEVENTSNEGATIVEQ2
        if(q2_fill < 0. && q2_fill> -0.01){
            q2_fill = 0.01;
        }
#endif
#ifdef INCLUDEEVENTSOVERFLOWQ2
        int nbins = util::NumQ2Bins(s);
        double MaxQ2_up   = util::Q2BinEdges(s,nbins);
        double MaxQ2_low  = util::Q2BinEdges(s,nbins-1);
        double MaxQ2_cent = (MaxQ2_up + MaxQ2_low)/2.;
        q2_fill = TMath::Min(MaxQ2_cent, q2_fill);
        
        if(q2_fill > MaxQ2_up){
            LOG("SuperkMC",pDEBUG) << q2_fill;
            util::Die(1); return;
        }
#endif
        
        if( q2_fill < util::Q2BinEdges(s,0) ){
            LOG("SuperkMC",pDEBUG) << "Q2 < Q2BinEdge(s,0)";
            LOG("SuperkMC",pDEBUG) << "nu direction: " << kBeamDir[0] << " " << kBeamDir[1] << " " << kBeamDir[2] << " ";
            LOG("SuperkMC",pDEBUG) << "lepton direction: " << fBrDir[0][0] << " " << fBrDir[0][1] << " " << fBrDir[0][2] << " ";
            LOG("SuperkMC",pDEBUG) << "theta = " << LepThetaDegreesSmeared << " ";
            LOG("SuperkMC",pDEBUG) << "costheta = " << costhetaSmeared << " ";
            //LOG("SuperkMC",pDEBUG) << "sintheta = " << sintheta << " ";
            LOG("SuperkMC",pDEBUG) << "emax_lep = " << emax_lepSmeared << " ";
            LOG("SuperkMC",pDEBUG) << "ereco_nu = " << ereco_nuSmeared << " ";
            LOG("SuperkMC",pDEBUG) << "Q2 = " << Q2 << " GeV";
            util::Die(1); return;
        }
        if( q2_fill > util::Q2BinEdges(s,nbins) ){
            LOG("SuperkMC",pDEBUG) << "Q2 > Q2BinEdge(s,0)";
            LOG("SuperkMC",pDEBUG) << "costheta = " << costhetaSmeared << " ";
            LOG("SuperkMC",pDEBUG) << "sintheta = " << TMath::Sin(LepThetaSmeared) << " ";
            LOG("SuperkMC",pDEBUG) << "emax_lep = " << emax_lepSmeared << " ";
            LOG("SuperkMC",pDEBUG) << "ereco_nu = " << ereco_nuSmeared << " ";
            LOG("SuperkMC",pDEBUG) << "Q2 = " << Q2 << " GeV";
            util::Die(1); return;
        }
        
        
#elif defined(__ERECO_VS_THETA_TEMPLATES__) // Use Ereco Vs Theta parametrization

        double q2_fill = LepThetaDegreesSmeared;
        
#ifdef INCLUDEEVENTSNEGATIVEQ2
        if(q2_fill < 0.){
            q2_fill = util::Q2BinEdges(s,0)+conv::kASmallNum;
        }
#endif
#ifdef INCLUDEEVENTSOVERFLOWQ2
        int nbins = util::NumQ2Bins(s);
        double MaxQ2_up   = util::Q2BinEdges(s,nbins);
        double MaxQ2_low  = util::Q2BinEdges(s,nbins-1);
        double MaxQ2_cent = (MaxQ2_up + MaxQ2_low)/2.;
        q2_fill = TMath::Min(MaxQ2_cent, q2_fill);
        
        if(q2_fill > MaxQ2_up){
            LOG("SuperkMC",pDEBUG) << q2_fill;
            util::Die(1); return;
        }
#endif
        
        if( q2_fill < util::Q2BinEdges(s,0) ){
            LOG("SuperkMC",pDEBUG) << "costheta = " << costhetaSmeared;
            LOG("SuperkMC",pDEBUG) << "LepThetaDegrees = " << LepThetaDegreesSmeared;
            util::Die(1); return;
        }
        
        if( q2_fill > util::Q2BinEdges(s,nbins) ){
            LOG("SuperkMC",pDEBUG) << "costheta = " << costhetaSmeared;
            LOG("SuperkMC",pDEBUG) << "LepThetaDegrees = " << LepThetaDegreesSmeared;
            util::Die(1); return;
        }
        
#else
        LOG("SuperkMC",pERROR) << "SuperkMC::FillTemplateHistogram()" ;
        LOG("SuperkMC",pERROR) << "Need to choose between Ereco Vs Q2 or Ereco Vs Theta parametrization" ;
        util::Die(1); return;
#endif //__ERECO_VS_Q2_TEMPLATES__ or __ERECO_VS_THETA_TEMPLATES__
        
        //if(util::IsRHCSample(s))
#ifdef USE_FLUX_EXT_WEIGHT
        fBrWghtFlx = fNuFlux->SuperkFluxExtTuningWght(s,m,Etrue);
#endif
        
#ifdef FLUXWEIGHTSAREONE
        fBrWghtFlx = 1.;
#endif
        
#ifdef INCLUDEEVENTSOVERTHIRTYGEV
        ereco_fill = TMath::Min(29.9, ereco_fill);
#endif
        
        
        //q2_fill = ereco_fill; // NEW 2D TEMPLATES (TODO) // be careful to overflows!!!
        
        
        
#ifndef __USE_2D_TEMPLATES__
        this -> Template(s,m,kSKMCtOri) -> Fill (ereco_nuSmeared, Etrue);
        fMCVariableTempl[sidx][midx]    -> Fill (ereco_nuSmeared, Etrue);
        ntotevents += 1.;
#else
        // NEW 2D TEMPLATES
        this -> Template2D(s,m,kSKMCtOri) -> Fill (q2_fill, ereco_nuSmeared, Etrue);
        fMCVariableTempl2D[sidx][midx]    -> Fill (q2_fill, ereco_nuSmeared, Etrue);
        ntotevents += 1.;
        //
#endif
        
        double pot_flx_weight = pot_wght * fBrWghtFlx;
        double NC1GammaWght = 1.;
#ifdef USE_NC1gamma_WEIGHTS
        if(util::IsNC1Gamma(TMath::Abs(fBrGenCode))){
            NC1GammaWght = 2.;
        }
#endif

#ifndef __USE_2D_TEMPLATES__
        this -> Template(s,m,kSKMCtPOT) -> Fill (ereco_nuSmeared, Etrue, pot_flx_weight);
#else
        //NEW 2D TEMPLATES
        this -> Template2D(s,m,kSKMCtPOT) -> Fill (q2_fill, ereco_nuSmeared, Etrue, pot_flx_weight);
        //
#endif
        
        // Actually fill the template
        
#ifndef __USE_2D_TEMPLATES__
#ifdef USE_RPA_WEIGHTS
	if(branchRPAnr){
	  this -> Template(s,m,kSKMCtNom) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtRPAnonrel * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght);
	} else {
	  this -> Template(s,m,kSKMCtNom) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght);
	}
        this -> Template(s,m,kSKMCtRPA) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtRPArel * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght);
#else
        this -> Template(s,m,kSKMCtNom) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght);
        this -> Template(s,m,kSKMCtRPA) -> Fill (ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght);
#endif
#else
        // NEW 2D TEMPLATES
#ifdef USE_RPA_WEIGHTS
	if(branchRPAnr){
	  this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtRPAnonrel * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght);
	} else {
	  this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght);
	}
        this -> Template2D(s,m,kSKMCtRPA) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtRPArel    * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght);
#elif defined(USE_BeRPA_WEIGHTS)
        this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtBeRPA );
        this -> Template2D(s,m,kSKMCtRPA) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght * fBrWghtBeRPA );        
#else
        this -> Template2D(s,m,kSKMCtNom) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght);
        this -> Template2D(s,m,kSKMCtRPA) -> Fill (q2_fill, ereco_fill, Etrue, pot_flx_weight * fBrWghtSFRFG * fBrWghtCCcoh * NC1GammaWght);
#endif
        //
#endif

#ifndef __USE_2D_TEMPLATES__
        
#ifdef __show_init_mesg_debug__
        if(Etrue <= 0 || Etrue >= this->Template(s,m,kSKMCtOri)->GetYaxis()->GetBinUpEdge(this->Template(s,m,kSKMCtOri)->GetNbinsY()))
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplateSmeared >> Etrue is out of template range " << AsString(s) << " " << AsString(m) << " " << Etrue
            << " (Ereco = " << ereco_nuSmeared << ")";
#endif
#ifdef __show_init_mesg_debug2__
        if(fBrEreco/1000 <= 0)
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplateSmeared >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << ereco_nuSmeared
            << " (Etrue = " << Etrue << ")";
        //#endif __show_init_mesg_debug2__
        if(fBrEreco/1000 >= this->Template(s,m,kSKMCtOri)->GetXaxis()->GetBinUpEdge(this->Template(s,m,kSKMCtOri)->GetNbinsX()))
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplateSmeared >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << ereco_nuSmeared
            << " (Etrue = " << Etrue << ")";
#endif
        
#else
        
        // NEW 2D TEMPLATES
        
#ifdef __show_init_mesg_debug__
        if(Etrue <= 0 || Etrue >= this->Template2D(s,m,kSKMCtOri)->GetZaxis()->GetBinUpEdge(this->Template2D(s,m,kSKMCtOri)->GetNbinsZ())){
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplateSmeared >> Etrue is out of template range " << AsString(s) << " " << AsString(m) << " " << Etrue
            << " (Ereco = " << ereco_nuSmeared << ")"
            << " (Q2 = " << q2_fill << ")";
        }
#endif
#ifdef __show_init_mesg_debug2__
        if(fBrEreco/1000 <= 0)
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplateSmeared >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << ereco_nuSmeared
            << " (Etrue = " << Etrue << ")";
        //#endif __show_init_mesg_debug2__
        if(fBrEreco/1000 >= this->Template2D(s,m,kSKMCtOri)->GetYaxis()->GetBinUpEdge(this->Template2D(s,m,kSKMCtOri)->GetNbinsY()))
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplateSmeared >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << ereco_nuSmeared
            << " (Etrue = " << Etrue << ")"
            << " (Q2 = " << q2_fill << ")";
        
        if(q2_fill <= 0)
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplateSmeared >> Q2 is out of template range " << AsString(s) << " " << AsString(m) << " " << q2_fill
            << " (Ereco = " << ereco_nuSmeared << ")"
            << " (Etrue = " << Etrue << ")";
        //#endif __show_init_mesg_debug2__
        if(q2_fill >= this->Template2D(s,m,kSKMCtOri)->GetXaxis()->GetBinUpEdge(this->Template2D(s,m,kSKMCtOri)->GetNbinsX()))
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplateSmeared >> Q2 is out of template range " << AsString(s) << " " << AsString(m) << " " << q2_fill
            << " (Ereco = " << ereco_nuSmeared << ")"
            << " (Etrue = " << Etrue << ")";
        //
#endif// __show_init_mesg_debug__
        
#endif
        
    }//iev
    
#ifdef __show_init_mesg_debug__
    
#ifndef __USE_2D_TEMPLATES__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplateSmeared >>      "
    << "Number of template entries = "
    << this -> Template(s,m,kSKMCtNom) -> GetEntries();
#else
    // NEW 2D TEMPLATES
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC/FillTemplate2DSmeared >>      "
    << "Number of template entries = "
    << this -> Template2D(s,m,kSKMCtNom) -> GetEntries();
    //
#endif
    
#endif
    
}
//________________________________________________________________________________________
void SuperkMC::FillTemplateArrays(bool with_niwg_weights) // refill the template arrays, with or without NIWG weights
{
    SKMCTemplateType_t template_to_use = (with_niwg_weights ? kSKMCtNom : kSKMCtPOT);
    
    fTemplateArraysHaveNIWGTune = with_niwg_weights;
    
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        const int      n_ereco_bins    = kNErecoBins[is];
        const int      n_etrue_bins    = kNEtrueBins[is];
        const double * ereco_bin_edges = util::ErecoBinEdges(s);
        const double * etrue_bin_edges = util::EtrueBinEdges(s);
        const double * var_ereco_bin_edges = util::XsecErecoBinEdges(s);
#ifdef __USE_2D_TEMPLATES__ //2D templates
        const int      n_q2_bins        = kNQ2Bins[is]; // NEW 2D TEMPLATES
        const double * q2_bin_edges     = util::Q2BinEdges(s); // NEW 2D TEMPLATES
        const double * var_q2_bin_edges = util::XsecQ2BinEdges(s); // NEW 2D TEMPLATES
#endif
        
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            
            const double * var_etrue_bin_edges = util::VariableMCEtrueBinEdges(s,m);
            
            
#ifndef __USE_2D_TEMPLATES__
            
            // loop over Etrue
            for(int it = 1; it <= this->Template(s,m,template_to_use)->GetYaxis()->GetNbins(); it++) {
                // loop over Ereco
                for(int ir = 1; ir <= this->Template(s,m,template_to_use)->GetXaxis()->GetNbins(); ir++) {
                    
                    fMCTemplArrayOri[is]->GetEntryRef(im, it-1 , ir-1) = this->Template(s,m,kSKMCtOri)->GetBinContent(ir,it);
                    ////
                    
                    double nev = this->Template(s,m,template_to_use)->GetBinContent(ir,it);
                    
                    fMCTemplArrayNom[is]->GetEntryRef(im, it-1 , ir-1) = nev;
                    fMCTemplArrayTun[is]->GetEntryRef(im, it-1 , ir-1) = nev;
                    fMCTemplArrayFin[is]->GetEntryRef(im, it-1 , ir-1) = nev;
                    
#ifdef USE_RPA_WEIGHTS || USE_BeRPA_WEIGHTS
                    double nev_relrpa = this->Template(s,m,kSKMCtRPA)->GetBinContent(ir,it);
                    if ( (!with_niwg_weights) || (nev_relrpa < conv::kASmallNum && nev < conv::kASmallNum)) {
                        fMCTemplArrayRPAWght[is]->GetEntryRef(im, it-1 , ir-1) = 1.;
                    }
                    else if  (nev < conv::kASmallNum) {
                        LOG("SuperkMC",pERROR)<<"Warning: divide by zero in RPA rel weighting, setting weight of 1";
                        fMCTemplArrayRPAWght[is]->GetEntryRef(im, it-1 , ir-1) = 1.;
                    }
                    else{
                        fMCTemplArrayRPAWght[is]->GetEntryRef(im, it-1 , ir-1) = nev_relrpa/nev;
                    }
#endif
                    ////
                    
                }//ir
            }//it
            
            int i_ne_var = 0;
            //now loop over the variable binning scheme (used to calculate weights only as many times as we need to)
            // This loop also keeps a record of which bins are non-empty
            //loop over variable Etrue
            for(int itv = 1; itv <= fMCVariableTempl[is][im]->GetYaxis()->GetNbins(); itv++) {
                double EtrueVarMin = var_etrue_bin_edges[itv-1] + kOneMeV;
                double EtrueVarMax = var_etrue_bin_edges[itv]   - kOneMeV;
                
                int EtrueMinBin = 1 + TMath::BinarySearch(n_etrue_bins,etrue_bin_edges,EtrueVarMin);
                int EtrueMaxBin = 1 + TMath::BinarySearch(n_etrue_bins,etrue_bin_edges,EtrueVarMax);
                
                // loop over variable Ereco
                for(int irv = 1; irv <= fMCVariableTempl[is][im]->GetXaxis()->GetNbins(); irv++) {
                    double nev = fMCVariableTempl[is][im]->GetBinContent(irv,itv);
                    
                    fMCVariableTemplArray[is]->GetEntryRef(im, itv-1, irv-1) = nev;
                    
                    // Fill nonempty bin numbers.
                    if (nev < conv::kASmallNum) {
                        continue;
                    }
                    
                    fNonEmptyVariableBins[is][im][i_ne_var] = conv::kBigIndex * irv + itv;
                    
                    double ErecoVarMin = var_ereco_bin_edges[irv-1] + kOneMeV;
                    double ErecoVarMax = var_ereco_bin_edges[irv]   - kOneMeV;
                    
                    int ErecoMinBin = 1 + TMath::BinarySearch(n_ereco_bins,ereco_bin_edges,ErecoVarMin);
                    int ErecoMaxBin = 1 + TMath::BinarySearch(n_ereco_bins,ereco_bin_edges,ErecoVarMax);
                    
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
                    fNNonEmptyFixedBins[is][im][i_ne_var] = i_ne_fixed;
                    fNonEmptyFixedBins [is][im][i_ne_var] = new int[i_ne_fixed];
                    
                    i_ne_fixed = 0;
                    // loop over fixed Etrue to count n nonempty
                    for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
                        // loop over fixed Ereco to count
                        for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
                            double nevf = this->Template(s,m,kSKMCtOri)->GetBinContent(ir,it);
                            if (nevf < conv::kASmallNum) {
                                continue;
                            }
                            fNonEmptyFixedBins[is][im][i_ne_var][i_ne_fixed] = conv::kBigIndex * ir + it;
                            i_ne_fixed++;
                        }//ir
                    }//it
                    
                    
                    i_ne_var++;
                }//irv
            }//itv
            fNNonEmptyVariableBins[is][im] = i_ne_var;
            assert(i_ne_var<kMaxVarBins);
            
            
            
#else  //__USE_2D_TEMPLATES__
            
            // NEW 2D TEMPLATES
            
            // loop over Etrue
            for(int it = 1; it <= this->Template2D(s,m,template_to_use)->GetZaxis()->GetNbins(); it++) {
                // loop over Ereco
                for(int ir = 1; ir <= this->Template2D(s,m,template_to_use)->GetYaxis()->GetNbins(); ir++) {
                    // loop over Q2
                    for(int iq = 1; iq <= this->Template2D(s,m,template_to_use)->GetXaxis()->GetNbins(); iq++) {
                        
                        fMCTemplArrayOri2D[is]->GetEntryRef(im, it-1 , ir-1, iq-1) = this->Template2D(s,m,kSKMCtOri)->GetBinContent(iq,ir,it);
                        
                        double nev = this->Template2D(s,m,template_to_use)->GetBinContent(iq,ir,it);
                        
                        fMCTemplArrayNom2D[is]->GetEntryRef(im, it-1 , ir-1, iq-1) = nev;
                        fMCTemplArrayTun2D[is]->GetEntryRef(im, it-1 , ir-1, iq-1) = nev;
                        fMCTemplArrayFin2D[is]->GetEntryRef(im, it-1 , ir-1, iq-1) = nev;
                        
#if defined(USE_RPA_WEIGHTS) || defined(USE_BeRPA_WEIGHTS)
                        double nev_relrpa = this->Template2D(s,m,kSKMCtRPA)->GetBinContent(iq,ir,it);
                        if ( (!with_niwg_weights) || (nev_relrpa < conv::kASmallNum && nev < conv::kASmallNum)) {
                            fMCTemplArrayRPAWght2D[is]->GetEntryRef(im, it-1 , ir-1, iq-1) = 1.;
                        }
                        else if  (nev < conv::kASmallNum) {
                            LOG("SuperkMC",pERROR)<<"Warning: divide by zero in RPA rel weighting, setting weight of 1";
                            fMCTemplArrayRPAWght2D[is]->GetEntryRef(im, it-1 , ir-1, iq-1) = 1.;
                        }
                        else{
                            fMCTemplArrayRPAWght2D[is]->GetEntryRef(im, it-1 , ir-1, iq-1) = nev_relrpa/nev;
                        }
#endif
                        ////
                    }//iq
                }//ir
            }//it
            
            int i_ne_var_2d = 0;
            //now loop over the variable binning scheme (used to calculate weights only as many times as we need to)
            //This loop also keeps a record of which bins are non-empty
            
            //loop over variable Etrue
            for(int itv = 1; itv <= fMCVariableTempl2D[is][im]->GetZaxis()->GetNbins(); itv++) {
                
                double EtrueVarMin = var_etrue_bin_edges[itv-1] + kOneMeV;
                double EtrueVarMax = var_etrue_bin_edges[itv]   - kOneMeV;
                int EtrueMinBin = 1 + TMath::BinarySearch(n_etrue_bins,etrue_bin_edges,EtrueVarMin);
                int EtrueMaxBin = 1 + TMath::BinarySearch(n_etrue_bins,etrue_bin_edges,EtrueVarMax);
                
                // loop over variable Ereco
                for(int irv = 1; irv <= fMCVariableTempl2D[is][im]->GetYaxis()->GetNbins(); irv++) {
                    
                    double ErecoVarMin = var_ereco_bin_edges[irv-1] + kOneMeV;
                    double ErecoVarMax = var_ereco_bin_edges[irv]   - kOneMeV;
                    int ErecoMinBin = 1 + TMath::BinarySearch(n_ereco_bins,ereco_bin_edges,ErecoVarMin);
                    int ErecoMaxBin = 1 + TMath::BinarySearch(n_ereco_bins,ereco_bin_edges,ErecoVarMax);
                    
                    // loop over variable Q2
                    for(int iqv = 1; iqv <= fMCVariableTempl2D[is][im]->GetXaxis()->GetNbins(); iqv++) {
                        
                        double Q2VarMin = var_q2_bin_edges[iqv-1] + kOneMeV;
                        double Q2VarMax = var_q2_bin_edges[iqv]   - kOneMeV;
                        int Q2MinBin = 1 + TMath::BinarySearch(n_q2_bins,q2_bin_edges,Q2VarMin);
                        int Q2MaxBin = 1 + TMath::BinarySearch(n_q2_bins,q2_bin_edges,Q2VarMax);
                        
                        double nev = fMCVariableTempl2D[is][im]->GetBinContent(iqv,irv,itv);
                        
                        fMCVariableTemplArray2D[is]->GetEntryRef(im, itv-1, irv-1, iqv-1) = nev;
                        
                        // Fill nonempty bin numbers.
                        if (nev < conv::kASmallNum) {
                            continue;
                        }
                        
                        fNonEmptyVariableBins2D[is][im][i_ne_var_2d] = conv::kBigIndex*conv::kBigIndex * iqv + conv::kBigIndex * irv + itv;
                        
                        int i_ne_fixed = 0;
                        // loop over fixed Etrue to count n nonempty
                        for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
                            // loop over fixed Ereco to count
                            for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
                                // loop over fixed Q2 to count
                                for(int iq = Q2MinBin; iq <= Q2MaxBin; iq++) {
                                    
                                    double nevf = this->Template2D(s,m,kSKMCtOri)->GetBinContent(iq,ir,it);
                                    if (nevf < conv::kASmallNum) {
                                        continue;
                                    }
                                    i_ne_fixed++;
                                }//iq
                            }//ir
                        }//it
                        fNNonEmptyFixedBins2D[is][im][i_ne_var_2d] = i_ne_fixed;
                        fNonEmptyFixedBins2D [is][im][i_ne_var_2d] = new int[i_ne_fixed];
                        
                        i_ne_fixed = 0;
                        // loop over fixed Etrue to count n nonempty
                        for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
                            // loop over fixed Ereco to count
                            for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
                                // loop over fixed Q2 to count
                                for(int iq = Q2MinBin; iq <= Q2MaxBin; iq++) {
                                    
                                    double nevf = this->Template2D(s,m,kSKMCtOri)->GetBinContent(iq,ir,it);
                                    if (nevf < conv::kASmallNum) {
                                        continue;
                                    }
                                    fNonEmptyFixedBins2D[is][im][i_ne_var_2d][i_ne_fixed] = conv::kBigIndex*conv::kBigIndex * iq + conv::kBigIndex * ir + it;
                                    i_ne_fixed++;
                                }//iq
                            }//ir
                        }//it
                        
                        i_ne_var_2d++;
                    }//iqv
                }//irv
            }//itv
            fNNonEmptyVariableBins2D[is][im] = i_ne_var_2d;
            assert(i_ne_var_2d<kMaxVarBins2D);
            
            //
            
#endif // __USE_2D_TEMPLATES__
            
            
        } //modes
    } //samples
    
    return;
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
        
#if defined(__USE_OSC_NUMU_MC_TEMPLATES__) || defined(__USE_OSC_NUMUBAR_MC_TEMPLATES__) || defined(__USE_OLD_OSC_NUEBAR_MC_TEMPLATES__)
        int sidx = arr_index->SuperkSampleIdx(s);
        
        EFluxType fNom = (IsRHCSample(s)) ? kRHCFlxNominal : kFlxNominal ;
        EFluxType fExt = (IsRHCSample(s)) ? kRHCFlxExtTuning : kFlxExtTuning ;
        
        const int n = 4;
#endif
        
        //
        // Fill in `osc numu' templates from `numu' templates
        // and take into account the nue/numu flux weight
        //
        
        
        
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        
#ifndef __USE_2D_TEMPLATES__ // 1D templates
        
        LOG("SuperkMC",pDEBUG) << "SKMC/FillSpecialTemplates >> Filling osc numu templates for " << util::AsString(s) << " sample";
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
                    
                    fMCTemplArrayOri[is]->GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_ori;
                    fMCTemplArrayNom[is]->GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
                    fMCTemplArrayTun[is]->GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
                    fMCTemplArrayFin[is]->GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
                    
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
                    
                    nev = fMCVariableTemplArray[is]->GetEntry(src_midx, itv-1, irv-1);
                    fMCVariableTemplArray[is]->GetEntryRef(tgt_midx, itv-1 , irv-1) = nev;
                    
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
        
        
#else // 2D templates
        
        // NEW 2D TEMPLATES
        
        LOG("SuperkMC",pDEBUG) << "SKMC/FillSpecialTemplates2D >> Filling osc numu templates for " << util::AsString(s) << " sample";
        SuperkMode_t src_modes[n] = {  kSKmNuMuCCQE,    kSKmNuMuCC1pi,     kSKmNuMuCCcoh,    kSKmNuMuCCoth    };
        SuperkMode_t tgt_modes[n] = {  kSKmOscNuMuCCQE, kSKmOscNuMuCC1pi,  kSKmOscNuMuCCcoh, kSKmOscNuMuCCoth };
        
        for(int i=0; i<n; i++) {
            
            SuperkMode_t tgt_m = tgt_modes[i];
            int tgt_midx = arr_index->SuperkModeIdx(s,tgt_modes[i]);
            int src_midx = arr_index->SuperkModeIdx(s,src_modes[i]);
            
            TH3D * src_ori = this->Template2D(s,src_modes[i],kSKMCtOri);
            TH3D * tgt_ori = this->Template2D(s,tgt_modes[i],kSKMCtOri);
            TH3D * src_nom = this->Template2D(s,src_modes[i],kSKMCtNom);
            TH3D * tgt_nom = this->Template2D(s,tgt_modes[i],kSKMCtNom);
            
            // loop over Etrue
            for(int it = 1; it <= src_ori->GetZaxis()->GetNbins(); it++) {
                double Etrue = src_ori->GetZaxis()->GetBinCenter(it);
                
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
                for(int ir = 1; ir <= src_ori->GetYaxis()->GetNbins(); ir++) {
                    double Ereco = src_ori->GetYaxis()->GetBinCenter(ir);
                    
                    // loop over Q2
                    for(int iq = 1; iq <= src_ori->GetXaxis()->GetNbins(); iq++) {
                        double Q2 = src_ori->GetXaxis()->GetBinCenter(iq);
                        
                        // get nevents
                        double nev_ori = src_ori-> GetBinContent(iq,ir,it);
                        double nev_nom = src_nom-> GetBinContent(iq,ir,it);
                        // set number of events
                        double weighted_nev_ori = nev_ori*wght_numu2nue_ori;
                        double weighted_nev_nom = nev_nom*wght_numu2nue_nom;
                        tgt_ori->Fill(Q2,Ereco,Etrue, weighted_nev_ori);
                        tgt_nom->Fill(Q2,Ereco,Etrue, weighted_nev_nom);
                        
                        double nev_tgt_ori = tgt_ori->GetBinContent(iq,ir,it);
                        double nev_tgt_nom = tgt_nom->GetBinContent(iq,ir,it);
                        
                        fMCTemplArrayOri2D[is]->GetEntryRef(tgt_midx, it-1, ir-1, iq-1) = nev_tgt_ori;
                        fMCTemplArrayNom2D[is]->GetEntryRef(tgt_midx, it-1, ir-1, iq-1) = nev_tgt_nom;
                        fMCTemplArrayTun2D[is]->GetEntryRef(tgt_midx, it-1, ir-1, iq-1) = nev_tgt_nom;
                        fMCTemplArrayFin2D[is]->GetEntryRef(tgt_midx, it-1, ir-1, iq-1) = nev_tgt_nom;
                        
                    }//iq
                }//ir
            }//it
            
            //one more loop, for the variable binning scheme.
            //just set the OscNuMu equal to the NuMu
            //loop over variable Etrue
            
            int i_ne_var = 0;
            
            for(int itv = 1; itv <= fMCVariableTempl2D[sidx][tgt_midx]->GetZaxis()->GetNbins(); itv++) {
                double EtrueVarMin = util::VariableMCEtrueBinEdges(s,tgt_m)[itv-1] + kOneMeV;
                double EtrueVarMax = util::VariableMCEtrueBinEdges(s,tgt_m)[itv]   - kOneMeV;
                
                int EtrueMinBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMin);
                int EtrueMaxBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMax);
                
                double nev = 0;
                
                // loop over variable Ereco
                for(int irv = 1; irv <= fMCVariableTempl2D[sidx][tgt_midx]->GetYaxis()->GetNbins(); irv++) {
                    
                    // loop over variable Q2
                    for(int iqv = 1; iqv <= fMCVariableTempl2D[sidx][tgt_midx]->GetXaxis()->GetNbins(); iqv++) {
                        
                        nev = fMCVariableTemplArray2D[is]->GetEntry(src_midx, itv-1, irv-1, iqv-1);
                        fMCVariableTemplArray2D[is]->GetEntryRef(tgt_midx, itv-1 , irv-1, iqv-1) = nev;
                        
                        // Fill nonempty bin numbers.
                        if (nev < conv::kASmallNum) {
                            continue;
                        }
                        
                        fNonEmptyVariableBins2D[sidx][tgt_midx][i_ne_var] = conv::kBigIndex*conv::kBigIndex * iqv + conv::kBigIndex * irv + itv;
                        
                        double Q2VarMin = util::XsecQ2BinEdges(s,iqv-1) + kOneMeV;
                        double Q2VarMax = util::XsecQ2BinEdges(s,iqv)   - kOneMeV;
                        int Q2MinBin = 1 + TMath::BinarySearch(kNQ2Bins[sidx],util::Q2BinEdges(s),Q2VarMin);
                        int Q2MaxBin = 1 + TMath::BinarySearch(kNQ2Bins[sidx],util::Q2BinEdges(s),Q2VarMax);
                        
                        double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
                        double ErecoVarMax = util::XsecErecoBinEdges(s,irv)   - kOneMeV;
                        int ErecoMinBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMin);
                        int ErecoMaxBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMax);
                        
                        int i_ne_fixed = 0;
                        // loop over fixed Etrue to count n nonempty
                        for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
                            // loop over fixed Ereco to count
                            for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
                                // loop over fixed Q2 to count
                                for(int iq = Q2MinBin; iq <= Q2MaxBin; iq++) {
                                    
                                    double nevf = this->Template2D(s,tgt_m,kSKMCtOri)->GetBinContent(iq,ir,it);
                                    if (nevf < conv::kASmallNum) {
                                        continue;
                                    }
                                    i_ne_fixed++;
                                }//iq
                            }//ir
                        }//it
                        fNNonEmptyFixedBins2D[sidx][tgt_midx][i_ne_var] = i_ne_fixed;
                        fNonEmptyFixedBins2D[sidx][tgt_midx][i_ne_var] = new int[i_ne_fixed];
                        
                        i_ne_fixed = 0;
                        // loop over fixed Etrue to count n nonempty
                        for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
                            // loop over fixed Ereco to count
                            for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
                                // loop over fixed Q2 to count
                                for(int iq = Q2MinBin; iq <= Q2MaxBin; iq++) {
                                    
                                    double nevf = this->Template2D(s,tgt_m,kSKMCtOri)->GetBinContent(iq,ir,it);
                                    if (nevf < conv::kASmallNum) {
                                        continue;
                                    }
                                    fNonEmptyFixedBins2D[sidx][tgt_midx][i_ne_var][i_ne_fixed] = conv::kBigIndex*conv::kBigIndex * iq + conv::kBigIndex * ir + it;
                                    i_ne_fixed++;
                                }//iq
                            }//ir
                        }//it
                        i_ne_var++;
                    }//iqv
                }//irv
            }//itv
            fNNonEmptyVariableBins2D[sidx][tgt_midx] = i_ne_var;
            assert(i_ne_var<kMaxVarBins2D);
            
        }//i
        
#endif // __USE_2D_TEMPLATES__
        
#endif// __USE_OSC_NUMU_MC_TEMPLATES__
        
        
        //
        // Fill in `osc numubar/nuebar' templates from `numubar/nuebar' templates
        // and take into account the numubar/nuebar flux weight
        //
        
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        
#ifndef __USE_2D_TEMPLATES__ // 1D templates
        
        LOG("SuperkMC",pDEBUG) << "SKMC/FillSpecialTemplates >> Filling osc numubar templates for " << util::AsString(s) << " sample";
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
                    
                    fMCTemplArrayOri[is]->GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_ori;
                    fMCTemplArrayNom[is]->GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
                    fMCTemplArrayTun[is]->GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
                    fMCTemplArrayFin[is]->GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
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
                    nev = fMCVariableTemplArray[is]->GetEntry(src_midx, itv-1 , irv-1);
                    fMCVariableTemplArray[is]->GetEntryRef(tgt_midx, itv-1, irv-1) = nev;
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
        
#else // 2D templates
        
        // NEW 2D TEMPLATES
        
        LOG("SuperkMC",pDEBUG) << "SKMC/FillSpecialTemplates2D >> Filling osc numubar templates for " << util::AsString(s) << " sample";
        SuperkMode_t src_modes_numubar[n] = {  kSKmNuMuBarCCQE,    kSKmNuMuBarCC1pi,     kSKmNuMuBarCCcoh,    kSKmNuMuBarCCoth    };
        SuperkMode_t tgt_modes_numubar[n] = {  kSKmOscNuMuBarCCQE, kSKmOscNuMuBarCC1pi,  kSKmOscNuMuBarCCcoh, kSKmOscNuMuBarCCoth };
        
        for(int i=0; i<n; i++) {
            
            SuperkMode_t tgt_m = tgt_modes_numubar[i];
            int tgt_midx = arr_index->SuperkModeIdx(s,tgt_modes_numubar[i]);
            int src_midx = arr_index->SuperkModeIdx(s,src_modes_numubar[i]);
            
            TH3D * src_ori = this->Template2D(s,src_modes_numubar[i],kSKMCtOri);
            TH3D * tgt_ori = this->Template2D(s,tgt_modes_numubar[i],kSKMCtOri);
            TH3D * src_nom = this->Template2D(s,src_modes_numubar[i],kSKMCtNom);
            TH3D * tgt_nom = this->Template2D(s,tgt_modes_numubar[i],kSKMCtNom);
            
            // loop over Etrue
            for(int it = 1; it <= src_ori->GetZaxis()->GetNbins(); it++) {
                double Etrue = src_ori->GetZaxis()->GetBinCenter(it);
                
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
                for(int ir = 1; ir <= src_ori->GetYaxis()->GetNbins(); ir++) {
                    double Ereco = src_ori->GetYaxis()->GetBinCenter(ir);
                    
                    // loop over Q2
                    for(int iq = 1; iq <= src_ori->GetXaxis()->GetNbins(); iq++) {
                        double Q2 = src_ori->GetXaxis()->GetBinCenter(iq);
                        
                        // get nevents
                        double nev_ori = src_ori-> GetBinContent(iq,ir,it);
                        double nev_nom = src_nom-> GetBinContent(iq,ir,it);
                        // set number of events
                        double weighted_nev_ori = nev_ori*wght_numu2nue_ori;
                        double weighted_nev_nom = nev_nom*wght_numu2nue_nom;
                        tgt_ori->Fill(Ereco,Etrue, weighted_nev_ori);
                        tgt_nom->Fill(Ereco,Etrue, weighted_nev_nom);
                        
                        double nev_tgt_ori = tgt_ori->GetBinContent(iq,ir,it);
                        double nev_tgt_nom = tgt_nom->GetBinContent(iq,ir,it);
                        
                        fMCTemplArrayOri2D[is]->GetEntryRef(tgt_midx, it-1, ir-1, iq-1) = nev_tgt_ori;
                        fMCTemplArrayNom2D[is]->GetEntryRef(tgt_midx, it-1, ir-1, iq-1) = nev_tgt_nom;
                        fMCTemplArrayTun2D[is]->GetEntryRef(tgt_midx, it-1, ir-1, iq-1) = nev_tgt_nom;
                        fMCTemplArrayFin2D[is]->GetEntryRef(tgt_midx, it-1, ir-1, iq-1) = nev_tgt_nom;
                    }//iq
                }//ir
            }//it
            //one more loop, for the variable binning scheme.
            //just set the OscNuMu equal to the NuMu
            //loop over variable Etrue
            int i_ne_var = 0;
            
            for(int itv = 1; itv <= fMCVariableTempl2D[sidx][tgt_midx]->GetZaxis()->GetNbins(); itv++) {
                double EtrueVarMin = util::VariableMCEtrueBinEdges(s,tgt_m)[itv-1] + kOneMeV;
                double EtrueVarMax = util::VariableMCEtrueBinEdges(s,tgt_m)[itv]   - kOneMeV;
                
                int EtrueMinBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMin);
                int EtrueMaxBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMax);
                
                // loop over variable Ereco
                for(int irv = 1; irv <= fMCVariableTempl2D[sidx][tgt_midx]->GetYaxis()->GetNbins(); irv++) {
                    
                    // loop over variable Q2
                    for(int iqv = 1; iqv <= fMCVariableTempl2D[sidx][tgt_midx]->GetXaxis()->GetNbins(); iqv++) {
                        
                        double nev = 0;
                        nev = fMCVariableTemplArray2D[is]->GetEntry(src_midx, itv-1, irv-1, iqv-1);
                        fMCVariableTemplArray2D[is]->GetEntryRef(tgt_midx, itv-1, irv-1, iqv-1) = nev;
                        ////
                        
                        // Fill nonempty bin numbers.
                        if (nev < conv::kASmallNum) {
                            continue;
                        }
                        
                        fNonEmptyVariableBins2D[sidx][tgt_midx][i_ne_var] = conv::kBigIndex*conv::kBigIndex * iqv + conv::kBigIndex * irv + itv;
                        
                        double Q2VarMin = util::XsecQ2BinEdges(s,iqv-1) + kOneMeV;
                        double Q2VarMax = util::XsecQ2BinEdges(s,iqv)   - kOneMeV;
                        int Q2MinBin = 1 + TMath::BinarySearch(kNQ2Bins[sidx],util::Q2BinEdges(s),Q2VarMin);
                        int Q2MaxBin = 1 + TMath::BinarySearch(kNQ2Bins[sidx],util::Q2BinEdges(s),Q2VarMax);
                        
                        double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
                        double ErecoVarMax = util::XsecErecoBinEdges(s,irv)   - kOneMeV;
                        int ErecoMinBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMin);
                        int ErecoMaxBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMax);
                        
                        int i_ne_fixed = 0;
                        // loop over fixed Etrue to count n nonempty
                        for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
                            // loop over fixed Ereco to count
                            for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
                                // loop over fixed Q2 to count
                                for(int iq = Q2MinBin; iq <= Q2MaxBin; iq++) {
                                    
                                    double nevf = this->Template2D(s,tgt_m,kSKMCtOri)->GetBinContent(iq,ir,it);
                                    if (nevf < conv::kASmallNum) {
                                        continue;
                                    }
                                    i_ne_fixed++;
                                }//iq
                            }//ir
                        }//it
                        fNNonEmptyFixedBins2D[sidx][tgt_midx][i_ne_var] = i_ne_fixed;
                        fNonEmptyFixedBins2D[sidx][tgt_midx][i_ne_var] = new int[i_ne_fixed];
                        
                        i_ne_fixed = 0;
                        // loop over fixed Etrue to count n nonempty
                        for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
                            // loop over fixed Ereco to count
                            for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
                                // loop over fixed Q2 to count
                                for(int iq = Q2MinBin; iq <= Q2MaxBin; iq++) {
                                    
                                    double nevf = this->Template2D(s,tgt_m,kSKMCtOri)->GetBinContent(iq,ir,it);
                                    if (nevf < conv::kASmallNum) {
                                        continue;
                                    }
                                    fNonEmptyFixedBins2D[sidx][tgt_midx][i_ne_var][i_ne_fixed] = conv::kBigIndex*conv::kBigIndex * iq + conv::kBigIndex * ir + it;
                                    i_ne_fixed++;
                                }//iq
                            }//ir
                        }//it
                        i_ne_var++;
                    }//iqv
                }//irv
            }//itv
            fNNonEmptyVariableBins2D[sidx][tgt_midx] = i_ne_var;
            assert(i_ne_var<kMaxVarBins2D);
        }//i
        
        
#endif // __USE_2D_TEMPLATES__
        
#endif // __USE_OSC_NUMUBAR_MC_TEMPLATES__
        
        
        
        
        // Old Osc NueBar MC templates
        
#ifdef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
        
#ifndef __USE_2D_TEMPLATES__ // 1D templates
        
        LOG("SuperkMC",pDEBUG) << "SKMC/FillSpecialTemplates >> Filling osc nuebar templates for " << util::AsString(s) << " sample";
        
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
                    
                    fMCTemplArrayOri[is]->GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_ori;
                    fMCTemplArrayNom[is]->GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
                    fMCTemplArrayTun[is]->GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
                    fMCTemplArrayFin[is]->GetEntryRef(tgt_midx, it-1, ir-1) = nev_tgt_nom;
                    
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
                    
                    nev = fMCVariableTemplArray[sidx]->GetEntry   (src_midx, itv-1, irv-1);
                    fMCVariableTemplArray[sidx]      ->GetEntryRef(tgt_midx, itv-1, irv-1) = nev;
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
        
        
#else // 2D templates
        
        // NEW 2D TEMPLATES
        
        LOG("SuperkMC",pDEBUG) << "SKMC/FillSpecialTemplates2D >> Filling osc nuebar templates for " << util::AsString(s) << " sample";
        
        SuperkMode_t src_modes_nuebar[n] = {  kSKmNuEBarCCQE,    kSKmNuEBarCC1pi,     kSKmNuEBarCCcoh,    kSKmNuEBarCCoth    };
        SuperkMode_t tgt_modes_nuebar[n] = {  kSKmOscNuEBarCCQE, kSKmOscNuEBarCC1pi,  kSKmOscNuEBarCCcoh, kSKmOscNuEBarCCoth };
        
        for(int i=0; i<n; i++) {
            
            SuperkMode_t tgt_m = tgt_modes_nuebar[i];
            int tgt_midx = arr_index->SuperkModeIdx(s,tgt_modes_nuebar[i]);
            int src_midx = arr_index->SuperkModeIdx(s,src_modes_nuebar[i]);
            
            TH3D * src_ori = this->Template2D(s,src_modes_nuebar[i],kSKMCtOri);
            TH3D * tgt_ori = this->Template2D(s,tgt_modes_nuebar[i],kSKMCtOri);
            TH3D * src_nom = this->Template2D(s,src_modes_nuebar[i],kSKMCtNom);
            TH3D * tgt_nom = this->Template2D(s,tgt_modes_nuebar[i],kSKMCtNom);
            
            // loop over Etrue
            for(int it = 1; it <= src_ori->GetZaxis()->GetNbins(); it++) {
                double Etrue = src_ori->GetZaxis()->GetBinCenter(it);
                
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
                for(int ir = 1; ir <= src_ori->GetYaxis()->GetNbins(); ir++) {
                    double Ereco = src_ori->GetYaxis()->GetBinCenter(ir);
                    
                    // loop over Q2
                    for(int iq = 1; iq <= src_ori->GetXaxis()->GetNbins(); iq++) {
                        double Q2 = src_ori->GetXaxis()->GetBinCenter(iq);
                        
                        // get nevents
                        double nev_ori = src_ori-> GetBinContent(iq,ir,it);
                        double nev_nom = src_nom-> GetBinContent(iq,ir,it);
                        // set number of events
                        double weighted_nev_ori = nev_ori*wght_nue2numu_ori;
                        double weighted_nev_nom = nev_nom*wght_nue2numu_nom;
                        tgt_ori->Fill(Q2,Ereco,Etrue, weighted_nev_ori);
                        tgt_nom->Fill(Q2,Ereco,Etrue, weighted_nev_nom);
                        
                        double nev_tgt_ori = tgt_ori->GetBinContent(iq,ir,it);
                        double nev_tgt_nom = tgt_nom->GetBinContent(iq,ir,it);
                        
                        fMCTemplArrayOri2D[is]->GetEntryRef(tgt_midx, it-1, ir-1, iq-1) = nev_tgt_ori;
                        fMCTemplArrayNom2D[is]->GetEntryRef(tgt_midx, it-1, ir-1, iq-1) = nev_tgt_nom;
                        fMCTemplArrayTun2D[is]->GetEntryRef(tgt_midx, it-1, ir-1, iq-1) = nev_tgt_nom;
                        fMCTemplArrayFin2D[is]->GetEntryRef(tgt_midx, it-1, ir-1, iq-1) = nev_tgt_nom;
                        
                    }//iq
                }//ir
            }//it
            //one more loop, for the variable binning scheme.
            //just set the OscNuMu equal to the NuMu
            //loop over variable Etrue
            int i_ne_var = 0;
            
            for(int itv = 1; itv <= fMCVariableTempl2D[sidx][tgt_midx]->GetZaxis()->GetNbins(); itv++) {
                double EtrueVarMin = util::VariableMCEtrueBinEdges(s,tgt_m)[itv-1] + kOneMeV;
                double EtrueVarMax = util::VariableMCEtrueBinEdges(s,tgt_m)[itv]   - kOneMeV;
                
                int EtrueMinBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMin);
                int EtrueMaxBin = 1 + TMath::BinarySearch(kNEtrueBins[sidx],util::EtrueBinEdges(s),EtrueVarMax);
                
                // loop over variable Ereco
                for(int irv = 1; irv <= fMCVariableTempl2D[sidx][tgt_midx]->GetYaxis()->GetNbins(); irv++) {
                    
                    // loop over variable Q2
                    for(int iqv = 1; iqv <= fMCVariableTempl2D[sidx][tgt_midx]->GetXaxis()->GetNbins(); iqv++) {
                        
                        double nev = 0;
                        nev = fMCVariableTemplArray2D[sidx]->GetEntry(src_midx, itv-1, irv-1, iqv-1);
                        fMCVariableTemplArray2D[sidx]->GetEntryRef(tgt_midx, itv-1, irv-1, iqv-1) = nev;
                        ////
                        
                        // Fill nonempty bin numbers.
                        if (nev < conv::kASmallNum) {
                            continue;
                        }
                        
                        fNonEmptyVariableBins2D[sidx][tgt_midx][i_ne_var] = conv::kBigIndex*conv::kBigIndex * iqv + conv::kBigIndex * irv + itv;
                        
                        double Q2VarMin = util::XsecQ2BinEdges(s,iqv-1) + kOneMeV;
                        double Q2VarMax = util::XsecQ2BinEdges(s,iqv)   - kOneMeV;
                        int Q2MinBin = 1 + TMath::BinarySearch(kNQ2Bins[sidx],util::Q2BinEdges(s),Q2VarMin);
                        int Q2MaxBin = 1 + TMath::BinarySearch(kNQ2Bins[sidx],util::Q2BinEdges(s),Q2VarMax);
                        
                        double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
                        double ErecoVarMax = util::XsecErecoBinEdges(s,irv)   - kOneMeV;
                        int ErecoMinBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMin);
                        int ErecoMaxBin = 1 + TMath::BinarySearch(kNErecoBins[sidx],util::ErecoBinEdges(s),ErecoVarMax);
                        
                        int i_ne_fixed = 0;
                        // loop over fixed Etrue to count n nonempty
                        for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
                            // loop over fixed Ereco to count
                            for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
                                // loop over fixed Q2 to count
                                for(int iq = Q2MinBin; iq <= Q2MaxBin; iq++) {
                                    
                                    double nevf = this->Template2D(s,tgt_m,kSKMCtOri)->GetBinContent(iq,ir,it);
                                    if (nevf < conv::kASmallNum) {
                                        continue;
                                    }
                                    i_ne_fixed++;
                                }//iq
                            }//ir
                        }//it
                        fNNonEmptyFixedBins2D[sidx][tgt_midx][i_ne_var] = i_ne_fixed;
                        fNonEmptyFixedBins2D[sidx][tgt_midx][i_ne_var] = new int[i_ne_fixed];
                        
                        i_ne_fixed = 0;
                        // loop over fixed Etrue to count n nonempty
                        for(int it = EtrueMinBin; it <= EtrueMaxBin; it++) {
                            // loop over fixed Ereco to count
                            for(int ir = ErecoMinBin; ir <= ErecoMaxBin; ir++) {
                                // loop over fixed Q2 to count
                                for(int iq = Q2MinBin; iq <= Q2MaxBin; iq++) {
                                    
                                    double nevf = this->Template2D(s,tgt_m,kSKMCtOri)->GetBinContent(iq,ir,it);
                                    if (nevf < conv::kASmallNum) {
                                        continue;
                                    }
                                    fNonEmptyFixedBins2D[sidx][tgt_midx][i_ne_var][i_ne_fixed] = conv::kBigIndex*conv::kBigIndex * ir + conv::kBigIndex * ir + it;
                                    i_ne_fixed++;
                                }//iq
                            }//ir
                        }//it
                        i_ne_var++;
                    }//iqv
                }//irv
            }//itv
            fNNonEmptyVariableBins2D[sidx][tgt_midx] = i_ne_var;
            assert(i_ne_var<kMaxVarBins2D);
            
        }//i
        
#endif // __USE_2D_TEMPLATES__
        
#endif// __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
        
        
    } // end loop kns
}
//________________________________________________________________________________________
void SuperkMC::CalculateMCStatError(void)
{
    // set MC frac err in each bin
    for(int is = 0; is < kns; is++)
    {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        fMCStatFracErr[is]->SetAll(0.);
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
                    
                    fMCStatFracErr[is]->GetEntryRef(im, it-1, ir-1) = fracerr;
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
                    double nev_old_bin  = fMCTemplArrayFin[is]->GetEntryRef(im, it-1, ir-1);
                    double frac_error   = fMCStatFracErr  [is]->GetEntryRef(im, it-1, ir-1);
                    double error        = frac_error * nev_old_bin;
                    double nev_new_bin  = rndgen.Gaus(nev_old_bin, error);
                    nev_new_bin = TMath::Max(0.,nev_new_bin);
                    nev_old += nev_old_bin;
                    nev_new += nev_new_bin;
                    fMCTemplArrayFin[is]->GetEntryRef(im, it-1, ir-1) = nev_new_bin;
                    
                }//it
            }//ir
            
            if(shape_only) {
                double scale = (nev_new > 0.) ? nev_old/nev_new : 0.;
                for(int ir = 1; ir <= n_ereco_bins; ir++) {
                    for(int it = 1; it <= n_etrue_bins; it++) {
                        double nev  = fMCTemplArrayFin[is]->GetEntryRef(im, it-1, ir-1);
                        fMCTemplArrayFin[is]->GetEntryRef(im, it-1, ir-1) = nev*scale;
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
    
    return fMCErecoSpectrumTotArray[sidx]->GetArray();
}
//________________________________________________________________________________________
double* SuperkMC::CurrErecoSpectrumArray(SuperkSample_t s)
{
    int sidx = arr_index->SuperkSampleIdx(s);
    if(sidx == -1) {
        return 0;
    }
    
    return fMCErecoSpectrumTotArray[sidx]->GetArray();
}
//________________________________________________________________________________________
const VALORArray * SuperkMC::CurrSpectrumVArray(SuperkSample_t s) const
{
    int sidx = arr_index->SuperkSampleIdx(s);
    if(sidx == -1) {
        return 0;
    }
    
    return fMCTemplArrayFin[sidx];
}
//________________________________________________________________________________________
VALORArray * SuperkMC::CurrSpectrumVArray(SuperkSample_t s)
{
    int sidx = arr_index->SuperkSampleIdx(s);
    if(sidx == -1) {
        return 0;
    }
    
    return fMCTemplArrayFin[sidx];
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
    
    return &(fMCErecoSpectrumArray[sidx]->GetEntryRef(midx,0));
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
    
    return &(fMCErecoSpectrumArray[sidx]->GetEntryRef(midx,0));
}
//________________________________________________________________________________________









// NEW 2D TEMPLATES
TH2D* SuperkMC::CurrErecoVsQ2Spectrum(SuperkSample_t s) const                               
{ 
  int is = arr_index->SuperkSampleIdx(s);
  int n_q2_bins    = kNQ2Bins[is];
  int n_ereco_bins = kNErecoBins[is];
  //int n_etrue_bins = kNEtrueBins[is];
  typedef double (erq2type)  [n_ereco_bins][n_q2_bins];
  const erq2type & arr     = *reinterpret_cast< const erq2type         * > ( CurrErecoVsQ2SpectrumArray(s) );
  if (!arr) return 0;
  //const double * arr = CurrErecoSpectrumArray(s,m);
  //if (!arr) return 0;

  string name  = Form("curr_erecoVsq2_%s"               , AsStringCompact(s).c_str() );
  string title = Form("Current Ereco Vs Q2 Spectrum: %s", AsString(s).c_str() );

  int nbins_q2    = util::NumQ2Bins(s);
  int nbins_ereco = util::NumErecoBins(s);
  TH2D * h = new TH2D(name.c_str(), title.c_str(), nbins_q2, util::Q2BinEdges(s), nbins_ereco, util::ErecoBinEdges(s) );
  util::FormatErecoVsQ2Hist(h);
  for (int iq = 0 ; iq < nbins_q2 ; iq++ ) {
    for (int ir = 0 ; ir < nbins_ereco ; ir++ ) {
      h->SetBinContent(iq+1,ir+1,arr[ir][iq]);
    }
  }
  return h;
}
//________________________________________________________________________________________
TH2D* SuperkMC::CurrErecoVsQ2Spectrum(SuperkSample_t s, SuperkMode_t m) const
{
    int is = arr_index->SuperkSampleIdx(s);
    int n_q2_bins    = kNQ2Bins[is];
    int n_ereco_bins = kNErecoBins[is];
    //int n_etrue_bins = kNEtrueBins[is];
    typedef double (erq2type)  [n_ereco_bins][n_q2_bins];
    const erq2type & arr     = *reinterpret_cast< const erq2type         * > ( CurrErecoVsQ2SpectrumArray(s) );
    if (!arr) return 0;
    //const double * arr = CurrErecoSpectrumArray(s,m);
    //if (!arr) return 0;
    
    string name  = Form("curr_erecoVsq2_%s_%s"               ,AsStringCompact(s).c_str(), AsStringCompact(m).c_str());
    string title = Form("Current Ereco Vs Q2 Spectrum: %s %s",AsString(s).c_str()       , AsString(m).c_str()       );
    
    int nbins_q2    = util::NumQ2Bins(s);
    int nbins_ereco = util::NumErecoBins(s);
    TH2D * h = new TH2D(name.c_str(), title.c_str(), nbins_q2, util::Q2BinEdges(s), nbins_ereco, util::ErecoBinEdges(s) );
    util::FormatErecoVsQ2Hist(h);
    util::SetDefaultSKSpectrumStyle(s,m,h);
    for (int iq = 0 ; iq < nbins_q2 ; iq++ ) {
        for (int ir = 0 ; ir < nbins_ereco ; ir++ ) {
            h->SetBinContent(iq+1,ir+1,arr[ir][iq]);
        }
    }
    return h;
}
//________________________________________________________________________________________
const double* SuperkMC::CurrErecoVsQ2SpectrumArray(SuperkSample_t s) const
{
    int sidx = arr_index->SuperkSampleIdx(s);
    if(sidx == -1) {
        return 0;
    }
    
    return fMCErecoVsQ2SpectrumTotArray2D[sidx]->GetArray();
}
//________________________________________________________________________________________
double* SuperkMC::CurrErecoVsQ2SpectrumArray(SuperkSample_t s)
{
    int sidx = arr_index->SuperkSampleIdx(s);
    if(sidx == -1) {
        return 0;
    }
    
    return fMCErecoVsQ2SpectrumTotArray2D[sidx]->GetArray();
}
//________________________________________________________________________________________
const VALORArray * SuperkMC::CurrSpectrumVArray2D(SuperkSample_t s) const
{
    int sidx = arr_index->SuperkSampleIdx(s);
    if(sidx == -1) {
        return 0;
    }
    
    return fMCTemplArrayFin2D[sidx];
}
//________________________________________________________________________________________
VALORArray * SuperkMC::CurrSpectrumVArray2D(SuperkSample_t s)
{
    int sidx = arr_index->SuperkSampleIdx(s);
    if(sidx == -1) {
        return 0;
    }
    
    return fMCTemplArrayFin2D[sidx];
}
//________________________________________________________________________________________
const double* SuperkMC::CurrErecoVsQ2SpectrumArray(SuperkSample_t s, SuperkMode_t m) const
{
    int sidx = arr_index->SuperkSampleIdx(s);
    if(sidx == -1) {
        return 0;
    }
    int midx = arr_index->SuperkModeIdx(s,m);
    if(midx == -1) {
        return 0;
    }
    
    return &(fMCErecoVsQ2SpectrumArray[sidx]->GetEntryRef(midx,0,0));
}
//________________________________________________________________________________________
double* SuperkMC::CurrErecoVsQ2SpectrumArray(SuperkSample_t s, SuperkMode_t m)
{
    int sidx = arr_index->SuperkSampleIdx(s);
    if(sidx == -1) {
        return 0;
    }
    int midx = arr_index->SuperkModeIdx(s,m);
    if(midx == -1) {
        return 0;
    }
    
    return &(fMCErecoVsQ2SpectrumArray[sidx]->GetEntryRef(midx,0,0));
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
        case (kSKMCtPOT) :
            h = fMCTemplPOT[sidx][midx];
            break;
        case (kSKMCtNom) :
            h = fMCTemplNom[sidx][midx];
            break;
        case (kSKMCtRPA) :
            h = fMCTemplRPA[sidx][midx];
            break;
        default:
            return 0;
    }
    assert(h);
    
    return h;
}
//________________________________________________________________________________________
// NEW 2D TEMPLATES
TH3D * SuperkMC::Template2D(SuperkSample_t s, SuperkMode_t m, SKMCTemplateType_t type) const
{
    int sidx = arr_index->SuperkSampleIdx(s);
    if(sidx == -1) {
        return 0;
    }
    int midx = arr_index->SuperkModeIdx(s,m);
    if(midx == -1) {
        return 0;
    }
    
    TH3D * h = 0;
    
    switch(type) {
        case (kSKMCtOri) :
            h = fMCTemplOri2D[sidx][midx];
            break;
        case (kSKMCtPOT) :
            h = fMCTemplPOT2D[sidx][midx];
            break;
        case (kSKMCtNom) :
            h = fMCTemplNom2D[sidx][midx];
            break;
        case (kSKMCtRPA) :
            h = fMCTemplRPA2D[sidx][midx];
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




// NEW 2D TEMPLATES

TH3D * SuperkMC::VarTemplate2D(SuperkSample_t s, SuperkMode_t m) const
{
    int sidx = arr_index->SuperkSampleIdx(s);
    if(sidx == -1) {
        return 0;
    }
    int midx = arr_index->SuperkModeIdx(s,m);
    if(midx == -1) {
        return 0;
    }
    
    TH3D * h = 0;
    h = fMCVariableTempl2D[sidx][midx];
    assert(h);
    
    return h;
}
//





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
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> Applying BANFF flux and cross-section tuning.";
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
            LOG("SuperkMC", pDEBUG) << "Applying BANFF tuning to sample, mode: " << AsString(s) << "\t" << AsString(m);
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
                if(fMCVariableTemplArray[is]->GetEntryRef(im, itv-1, irv-1) <= 0) continue;
                
                double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
                
                // Plug-in BANFF & other xsec systematics
                double xsec_wght = fNuXSec->FinTuningWght(s,m,ErecoVarMin,EtrueVarMin);
                
                // Combined weight
                double wght = flux_wght * xsec_wght;
                
#ifdef BANFFTUNEISONE
                wght = 1;
#endif
                for (int ifb = 0 ; ifb < fNNonEmptyFixedBins[is][im][ivb] ; ifb++ ) {
                    int it =  fNonEmptyFixedBins[is][im][ivb][ifb] % conv::kBigIndex;
                    int ir = (fNonEmptyFixedBins[is][im][ivb][ifb] - it)/conv::kBigIndex;
                    
                    //check for empty bins
                    if(fMCTemplArrayOri[is]->GetEntryRef(im, it-1, ir-1) <= 0) continue;
                    ////
                    //double Etrue = util::EtrueBinEdges(s,it-1) + kOneMeV;
                    
                    // get number of events in current bin of source MC template
                    double nev = fMCTemplArrayNom[is]->GetEntry(im, it-1, ir-1);
                    
                    nev = TMath::Max(0.,nev);
                    double nev_tuned = nev * wght;
                    
#ifdef NO_EVENTS_CAN_BE_NEGATIVE
                    if(nev_tuned < 0) {
#ifdef __show_negative_events_mesg_debug__
                        LOG("SuperkMC",pDEBUG) << "Exp/SKMC/ApplyBANFFTuningToTemplates >> set negative number of events " << nev_tuned
                        << " to " << 0 << "\t" << " Ereco: " << Ereco << " Etrue: " << util::EtrueBinEdges(s,it-1) + kOneMeV;
#endif
                        nev_tuned = 0;
                    }
#endif //NO_EVENTS_CAN_BE_NEGATIVE
                    
                    fMCTemplArrayTun[is]->GetEntryRef(im, it-1, ir-1) = nev_tuned;
                    fMCTemplArrayFin[is]->GetEntryRef(im, it-1, ir-1) = nev_tuned;
                    
#ifdef __show_banfftun_syst_var_wght_mesg_debug__
                    // Get the energy of the current bin
                    double Etrue = util::EtrueBinEdges(s,it-1) + kOneMeV;
                    double Ereco = util::ErecoBinEdges(s,it-1) + kOneMeV;
                    if( flux_wght  != 1. ||
                       xsec_wght  != 1. )
                    {
                      LOG("SuperkMC", pDEBUG) << "SystVarWght >> "
                        << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m)
                        << ", Ereco = " << Ereco << " GeV"
                        << ", Etrue = " << Etrue << " GeV => Weights: "
                        << ", Flux = "  << flux_wght
                        << ", XSec = "  << xsec_wght
                        << " (total = " << wght << ")";
                    }
#endif //__show_banfftun_syst_var_wght_mesg_debug__
                    
                }//ifb
            }//ivb
        }//im
    }//is
    
#ifdef __show_init_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC/ApplyBANFFTuningToTemplates >> Complete";
#endif
}
//________________________________________________________________________________________










// NEW 2D TEMPLATES

void SuperkMC::ApplyBANFFTuningToTemplates2D(void)
{
    // Take BANFF tuning into account (reweight MC to BANFF best-fit parameters)
    // The weights returned by NuFlux::SuperkFluxFinTuningWght() are relative to the flux
    // tuned based on external (NA61,...) data, so make sure that you run
    // SuperkMC::ApplyFluxExtTuningToTemplates() before running this method.
    // Source template type: kSKMCtNom
    // Target template type: kSKMCtTun
    //
    
#ifdef __show_init_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> Applying BANFF 2D flux and cross-section tuning.";
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
             LOG("SuperkMC", pDEBUG) << "Applying BANFF 2D tuning to sample, mode: " << AsString(s) << "\t" << AsString(m);
#endif
            
            //
            // There are 600 50MeV Etrue energy bins in the templates, but
            // the flux/xsec systematics only have up to ~10 Etrue bins (same as flux binning), so
            // we can call NuXSec::FinTuningWght() and NuFlux::SuperkFluxFinTuningWght()
            // a factor of 10 less often by looping over the kXsEtrueBinEdges.
            // Note that the Ereco binning is identical between XSec and SuperkMC templates
            //
            for (int ivb = 0 ; ivb < fNNonEmptyVariableBins2D[is][im] ; ivb++ ) {
                
                // The array is built as following: fNonEmptyVariableBins2D = conv::kBigIndex*conv::kBigIndex * iqv + conv::kBigIndex * irv + itv;
                int itv = ( fNonEmptyVariableBins2D[is][im][ivb] % (conv::kBigIndex * conv::kBigIndex) ) % conv::kBigIndex;
                int irv = ( ( fNonEmptyVariableBins2D[is][im][ivb] % (conv::kBigIndex * conv::kBigIndex) ) - itv ) / conv::kBigIndex;
                int iqv = ( fNonEmptyVariableBins2D[is][im][ivb] - (conv::kBigIndex * irv + itv) ) / (conv::kBigIndex * conv::kBigIndex);
                
                // Get the min/max energy of the current variable flux Etrue bin
                double EtrueVarMin = util::VariableMCEtrueBinEdges(s,m)[itv-1] + kOneMeV;
                
                // Plug-in BANFF flux systematics
                double flux_wght = fNuFlux->SuperkFluxFinTuningWght(s,m,EtrueVarMin);
                
                //skip if the big bin is empty
                if(fMCVariableTemplArray2D[is]->GetEntryRef(im, itv-1, irv-1, iqv-1) <= 0) continue;
                
                double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
                double Q2VarMin    = util::XsecQ2BinEdges(s,iqv-1) +kOneMeV;
                // Plug-in BANFF & other xsec systematics
                //double xsec_wght = fNuXSec->FinTuningWght(s,m,ErecoVarMin,EtrueVarMin,Q2VarMin); // Apply systematics in Ereco and Q2
                // Don't apply the xsec weight here, this is now handled in ApplySystematicsToTemplates2D
                double xsec_wght = 1.0;
                
                // Combined weight
                double wght = flux_wght * xsec_wght;
                
#ifdef BANFFTUNEISONE
                wght = 1;
#endif
                for (int ifb = 0 ; ifb < fNNonEmptyFixedBins2D[is][im][ivb] ; ifb++ ) {
                    
                    // The array is built as following: fNonEmptyVariableBins2D[is][im][i_ne_var_2d] = conv::kBigIndex*conv::kBigIndex * iqv + conv::kBigIndex * irv + itv;
                    int it = ( fNonEmptyFixedBins2D[is][im][ivb][ifb] % (conv::kBigIndex * conv::kBigIndex) ) % conv::kBigIndex;
                    int ir = ( ( fNonEmptyFixedBins2D[is][im][ivb][ifb] % (conv::kBigIndex * conv::kBigIndex) ) - it ) / conv::kBigIndex;
                    int iq = ( fNonEmptyFixedBins2D[is][im][ivb][ifb] - (conv::kBigIndex * ir + it) ) / (conv::kBigIndex * conv::kBigIndex);
                    
                    //check for empty bins
                    if(fMCTemplArrayOri2D[is]->GetEntryRef(im, it-1, ir-1, iq-1) <= 0) continue;
                    ////
                    //double Etrue = util::EtrueBinEdges(s,it-1) + kOneMeV;
                    
                    // get number of events in current bin of source MC template
                    double nev = fMCTemplArrayNom2D[is]->GetEntry(im, it-1, ir-1, iq-1);
                    
                    nev = TMath::Max(0.,nev);
                    double nev_tuned = nev * wght;
                    
#ifdef NO_EVENTS_CAN_BE_NEGATIVE
                    if(nev_tuned < 0) {
#ifdef __show_negative_events_mesg_debug__
                        LOG("SuperkMC",pDEBUG) << "Exp/SKMC/ApplyBANFFTuningToTemplates2D >> set negative number of events " << nev_tuned
                        << " to " << 0 << "\t" << " Ereco: " << Ereco << " Etrue: " << util::EtrueBinEdges(s,it-1) + kOneMeV;
#endif
                        nev_tuned = 0;
                    }
#endif //NO_EVENTS_CAN_BE_NEGATIVE
                    
                    fMCTemplArrayTun2D[is]->GetEntryRef(im, it-1, ir-1, iq-1) = nev_tuned;
                    fMCTemplArrayFin2D[is]->GetEntryRef(im, it-1, ir-1, iq-1) = nev_tuned;
                    
#ifdef __show_banfftun_syst_var_wght_mesg_debug__
                    // Get the energy of the current bin
                    double Etrue = util::EtrueBinEdges(s,it-1) + kOneMeV;
                    double Ereco = util::ErecoBinEdges(s,ir-1) + kOneMeV;
                    double Q2    = util::Q2BinEdges(s,iq-1)    + kOneMeV;
                    if( flux_wght  != 1. ||
                       xsec_wght  != 1. )
                    {
                      LOG("SuperkMC", pDEBUG) <<  "\nSystVarWght >> "
                        << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m)
                        << ", Q2 = " << Q2 << " GeV"
                        << ", Ereco = " << Ereco << " GeV"
                        << ", Etrue = " << Etrue << " GeV => Weights: "
                        << ", Flux = "  << flux_wght
                        << ", XSec = "  << xsec_wght
                        << " (total = " << wght << ")"
                        << " Nom events = " << nev
                        << " New events = " << nev_tuned;
                    }
#endif //__show_banfftun_syst_var_wght_mesg_debug__
                    
                }//ifb
            }//ivb
        }//im
    }//is
    
#ifdef __show_init_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC/ApplyBANFFTuningToTemplates2D >> Complete";
#endif
}

//



//________________________________________________________________________________________
void SuperkMC::UpdateTemplates(void)
{
#ifdef __show_update_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> Updating MC templates....";
#endif
    
#ifndef __USE_2D_TEMPLATES__
    this->ApplySystematicsToTemplates   ();
    this->ApplyPhysHypothesisToTemplates();
#else
    this->ApplySystematicsToTemplates2D   ();
    this->ApplyPhysHypothesisToTemplates2D();
#endif
    
#ifdef __INCLUDE_MC_STAT_FLUCTUATIONS__
#ifdef __show_update_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> ** Note ** Fluctuating templates to take into account the effect of finite MC statistics....";
#endif
#ifndef __USE_2D_TEMPLATES__
    this->ApplyMCStatFluctationToTemplates(true); // Implemented only for 1D templates
#endif
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




// NEW 2D TEMPLATES
void SuperkMC::ApplySystematicsToTemplates2D(void)
{
    // Apply systematic variations to MC templates
    //
    this->ApplySystematicsToTemplates_Multiplicative_2D   ();
}
//



//________________________________________________________________________________________
void SuperkMC::ApplySystematicsToTemplates_Multiplicative(void)
{
    // All multiplicative systematics are handled in the following loop.
    //
#ifdef __show_update_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> Applying multiplicative systematics";
#endif
    
    //~ fNuXSec->SetSKMCArraysHaveNIWGTune(fTemplateArraysHaveNIWGTune);
    
    // Loop over samples
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        bool is_rhc = util::IsRHCSample(s);
        
        double rpa_tweakval = fNuXSec->CurrentValue( is_rhc ? conv::kRHCNuisXSecRPA : conv::kNuisXSecRPA);
        
        bool relativisticRPA = (bool) rpa_tweakval;
        assert ( fabs ( (double)relativisticRPA - rpa_tweakval) < conv::kASmallNum); // Check RPA was exactly 1 or 0
        
        // Loop over sample modes
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            // Consider template?
            bool skip = this->SkipTemplate(m);
            if(skip) continue;
            
            bool applyRPAtoThisMode = false;
            if (util::IsRHCSample(s)) applyRPAtoThisMode = fNuXSec->ApplySystToSKMode(conv::kRHCNuisXSecRPA, s, m) ;
            else                      applyRPAtoThisMode = fNuXSec->ApplySystToSKMode(conv::kNuisXSecRPA   , s, m) ;
            
            const double * variableMCEtrueBinEdges = util::VariableMCEtrueBinEdges(s,m);
            
            // Recast VALORArray blocks so they can be used addressed like normal arrays (fast).
            // Hard to read, so should only be used in code for which speed is essential.
            int n_ereco_bins = kNErecoBins[is];
            int n_etrue_bins = kNEtrueBins[is];
            typedef double (templtype)[n_etrue_bins][n_ereco_bins];
            const templtype & templTun = *reinterpret_cast<const templtype * >(&(fMCTemplArrayTun    [is]->GetEntryRef(im,0,0)));
            const templtype & templRPA = *reinterpret_cast<const templtype * >(&(fMCTemplArrayRPAWght[is]->GetEntryRef(im,0,0)));
            templtype       & templFin = *reinterpret_cast<templtype       * >(&(fMCTemplArrayFin    [is]->GetEntryRef(im,0,0)));
            
            //~ const double(&templTun)[n_etrue_bins][n_ereco_bins] = *reinterpret_cast< const double (*) [n_etrue_bins][n_ereco_bins]> (&(fMCTemplArrayTun[is]->GetEntryRef(im,0,0)));
            //~ double      (&templFin)[n_etrue_bins][n_ereco_bins] = *reinterpret_cast< double       (*) [n_etrue_bins][n_ereco_bins]> (&(fMCTemplArrayFin[is]->GetEntryRef(im,0,0)));
            
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
                
                //assert(fMCVariableTemplArray[is]->GetEntryRef(im, itv-1, irv-1) > 0);
                
                double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
                
                // Plug-in BANFF & other xsec systematics
                double xsec_wght = fNuXSec->SystWght(s,m,ErecoVarMin,EtrueVarMin);
                
#ifdef XSECWEIGHTSAREONE
                xsec_wght = 1;
#endif
                // Plug-in SK efficiency + FSI uncertainties
                double skfsi_wght = fSKDetFSI->SystWght(s,m,ErecoVarMin,EtrueVarMin);
                
                // if(s==kSKsFCFVMultiRingMuLike || s==kSKsFCFVMultiRingELike){
                //   util::Die(1); return;
                // }
                
                
                // Combined weight
                double wght = skfsi_wght * flux_wght * xsec_wght;
                
#ifdef SYSTWEIGHTSAREONE
                wght = 1;
#endif
                // loop over 50MeV Ereco
                for (int ifb = 0 ; ifb < fNNonEmptyFixedBins[is][im][ivb] ; ifb++ ) {
                    int it =  fNonEmptyFixedBins[is][im][ivb][ifb] % conv::kBigIndex;
                    int ir = (fNonEmptyFixedBins[is][im][ivb][ifb] - it)/conv::kBigIndex;
                    //check for empty bins
                    
                    //assert(fMCTemplOri[is]->GetEntryRef(im, itv-1, irv-1) > 0);
                    
                    ////
                    //double Etrue = util::ErecoBinEdges(s)(s,it-1) + kOneMeV;
                    
                    // get number of events in current bin of source MC template
                    double nev = templTun[it-1][ir-1];
                    
                    nev = TMath::Max(0.,nev);
                    double nev_systm = nev * wght;
#ifndef SYSTWEIGHTSAREONE
                    if (relativisticRPA && applyRPAtoThisMode) {
                        nev_systm *= templRPA[it-1][ir-1]; // apply nonrel to rel RPA tune
                    }
#endif
                    
#ifdef __show_syst_var_wght_mesg_debug__
                    // Get the energy of the current bin
                    double Etrue = util::EtrueBinEdges(s,it-1) + kOneMeV;
                    double Ereco = util::ErecoBinEdges(s,it-1) + kOneMeV;
		    if( skfsi_wght > 0){
                    if( skfsi_wght != 1. || 
                        flux_wght  != 1. ||
                        xsec_wght  != 1. )
                    {
                        LOG("SuperkMC",pDEBUG) << "Exp/SKMC/SystVarWght >> "
                        << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m)
                        << ", Ereco = " << Ereco << " GeV"
                        << ", Etrue = " << Etrue << " GeV \n=> Weights: "
                        << "SK+FSI = "  << skfsi_wght 
                        << ", Flux = "  << flux_wght
                        << ", XSec = "  << xsec_wght
                        << " (total = " << wght << ")";
                    }
		    }
#endif//__show_syst_var_wght_mesg_debug__
                    
#ifdef NO_EVENTS_CAN_BE_NEGATIVE
                    if(nev_systm < 0) {
#ifdef __show_negative_events_mesg_debug__
                        LOG("SuperkMC",pDEBUG) << "Exp/SKMC/ApplySystematicsToTemplates_Multiplicative >> set negative number of events " << nev_systm
                        << " to " << 0 << "\t" << " Ereco: " << Ereco;
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





// NEW 2D TEMPLATES

void SuperkMC::ApplySystematicsToTemplates_Multiplicative_2D(void)
{
    // All multiplicative systematics are handled in the following loop.
    //
#ifdef __show_update_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> Applying multiplicative systematics";
#endif
    
    // Loop over samples
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        bool is_rhc = util::IsRHCSample(s);
        
        double rpa_tweakval = fNuXSec->CurrentValue( is_rhc ? conv::kRHCNuisXSecRPA : conv::kNuisXSecRPA);
        
        bool relativisticRPA = (bool) rpa_tweakval;
        assert ( fabs ( (double)relativisticRPA - rpa_tweakval) < conv::kASmallNum); // Check RPA was exactly 1 or 0
        
        // Loop over sample modes
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            // Consider template?
            bool skip = this->SkipTemplate(m);
            if(skip) continue;
            
            bool applyRPAtoThisMode = false;
            if (util::IsRHCSample(s)) applyRPAtoThisMode = fNuXSec->ApplySystToSKMode(conv::kRHCNuisXSecRPA, s, m) ;
            else                      applyRPAtoThisMode = fNuXSec->ApplySystToSKMode(conv::kNuisXSecRPA   , s, m) ;
            
            const double * variableMCEtrueBinEdges = util::VariableMCEtrueBinEdges(s,m);
            
            // Recast VALORArray blocks so they can be used addressed like normal arrays (fast).
            // Hard to read, so should only be used in code for which speed is essential.
            int n_q2_bins    = kNQ2Bins[is];
            int n_ereco_bins = kNErecoBins[is];
            int n_etrue_bins = kNEtrueBins[is];
            typedef double (templtype)[n_etrue_bins][n_ereco_bins][n_q2_bins];
            const templtype & templTun = *reinterpret_cast<const templtype * >(&(fMCTemplArrayTun2D    [is]->GetEntryRef(im,0,0,0)));
            const templtype & templRPA = *reinterpret_cast<const templtype * >(&(fMCTemplArrayRPAWght2D[is]->GetEntryRef(im,0,0,0)));
            templtype       & templFin = *reinterpret_cast<templtype       * >(&(fMCTemplArrayFin2D    [is]->GetEntryRef(im,0,0,0)));
            
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
            for (int ivb = 0 ; ivb < fNNonEmptyVariableBins2D[is][im] ; ivb++ ) {
                
                // The array is built as following: fNonEmptyVariableBins2D[is][im][i_ne_var_2d] = conv::kBigIndex*conv::kBigIndex * iqv + conv::kBigIndex * irv + itv;
                int itv = ( fNonEmptyVariableBins2D[is][im][ivb] % (conv::kBigIndex * conv::kBigIndex) ) % conv::kBigIndex;
                int irv = ( ( fNonEmptyVariableBins2D[is][im][ivb] % (conv::kBigIndex * conv::kBigIndex) ) - itv ) / conv::kBigIndex;
                int iqv = ( fNonEmptyVariableBins2D[is][im][ivb] - (conv::kBigIndex * irv + itv) ) / (conv::kBigIndex * conv::kBigIndex);
                
                double EtrueVarMin = variableMCEtrueBinEdges[itv-1] + kOneMeV;
                
                // Plug-in BANFF flux systematics
                if (itv!=last_itv) {
                    flux_wght = fNuFlux->SuperkFluxSystWght(s,m,EtrueVarMin);
                    last_itv = itv;
                }
                
                //assert(fMCVariableTemplArray[is]->GetEntryRef(im, itv-1, irv-1) > 0);
                
                double ErecoVarMin = util::XsecErecoBinEdges(s,irv-1) + kOneMeV;
                double Q2VarMin = util::XsecQ2BinEdges(s,iqv-1) + kOneMeV;
                
                // Plug-in BANFF & other xsec systematics
                double xsec_wght = fNuXSec->SystWght(s,m,ErecoVarMin,EtrueVarMin,Q2VarMin);
                
#ifdef XSECWEIGHTSAREONE
                xsec_wght = 1;
#endif
                // Plug-in SK efficiency + FSI uncertainties
                double skfsi_wght = fSKDetFSI->SystWght(s,m,ErecoVarMin,EtrueVarMin);
		
		// assert(skfsi_wght > 0 && "skfsi_wght is  <= 0");
		
                // if(s==kSKsFCFVMultiRingMuLike || s==kSKsFCFVMultiRingELike){
                //   util::Die(1); return;;
                // }
                
                
                // Combined weight
                double wght = skfsi_wght * flux_wght * xsec_wght;
                
#ifdef SYSTWEIGHTSAREONE
                wght = 1;
#endif
                // loop over 50MeV Ereco
                for (int ifb = 0 ; ifb < fNNonEmptyFixedBins2D[is][im][ivb] ; ifb++ ) {
                    
                    // The array is built as following: fNonEmptyVariableBins2D[is][im][i_ne_var_2d] = conv::kBigIndex*conv::kBigIndex * iqv + conv::kBigIndex * irv + itv;
                    int it = ( fNonEmptyFixedBins2D[is][im][ivb][ifb] % (conv::kBigIndex * conv::kBigIndex) ) % conv::kBigIndex;
                    int ir = ( ( fNonEmptyFixedBins2D[is][im][ivb][ifb] % (conv::kBigIndex * conv::kBigIndex) ) - it ) / conv::kBigIndex;
                    int iq = ( fNonEmptyFixedBins2D[is][im][ivb][ifb] - (conv::kBigIndex * ir + it) ) / (conv::kBigIndex * conv::kBigIndex);
                    
                    //check for empty bins
                    
                    //assert(fMCTemplOri[is]->GetEntryRef(im, itv-1, irv-1) > 0);
                    
                    ////
                    //double Etrue = util::ErecoBinEdges(s)(s,it-1) + kOneMeV;
                    
                    // get number of events in current bin of source MC template
                    double nev = templTun[it-1][ir-1][iq-1];
                    
                    nev = TMath::Max(0.,nev);
                    double nev_systm = nev * wght;
#ifndef SYSTWEIGHTSAREONE
                    if (relativisticRPA && applyRPAtoThisMode) {
                        nev_systm *= templRPA[it-1][ir-1][iq-1]; // apply nonrel to rel RPA tune
                    }
#endif
                    
#ifdef __show_syst_var_wght_mesg_debug__
                    // Get the energy of the current bin
                    double Etrue = util::EtrueBinEdges(s,it-1) + kOneMeV;
                    double Ereco = util::ErecoBinEdges(s,ir-1) + kOneMeV;
                    double Q2    = util::Q2BinEdges(s,iq-1)    + kOneMeV;
		    if( skfsi_wght > 0){
                    if( skfsi_wght != 1. ||
                        flux_wght  != 1. ||
                        xsec_wght  != 1. )
                    {
                        LOG("SuperkMC",pDEBUG) << "Exp/SKMC/SystVarWght >> "
                        << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m)
                        << ", Q2 = " << Q2 << " GeV"
                        << ", Ereco = " << Ereco << " GeV"
                        << ", Etrue = " << Etrue << " GeV \n=> Weights: "
                        << "SK+FSI = "  << skfsi_wght
                        << ", Flux = "  << flux_wght
                        << ", XSec = "  << xsec_wght
                        << ", RPA  = "  << templRPA[it-1][ir-1][iq-1]
                        << " (total = " << wght << ")"
                        << " Nom events = " << nev
                        << " New events = " << nev_systm;
                    }
		    }
#endif//__show_syst_var_wght_mesg_debug__
                    
#ifdef NO_EVENTS_CAN_BE_NEGATIVE
                    if(nev_systm < 0) {
#ifdef __show_negative_events_mesg_debug__
                        LOG("SuperkMC",pDEBUG) << "Exp/SKMC/ApplySystematicsToTemplates_Multiplicative >> set negative number of events " << nev_systm
                        << " to " << 0 << "\t" << " Ereco: " << Ereco << ", Q2: " << Q2;
#endif//__show_negative_events_mesg_debug__
                        nev_systm = 0;
                    }
#endif// NO_EVENTS_CAN_BE_NEGATIVE
                    // Update target MC template
                    
                    templFin[it-1][ir-1][iq-1] = nev_systm;
                    
                }//ifb
            }//ivb
        }//im
        
    }//is
    
#ifdef __save_templates_debug__
    this->NevDebug2D(kSKMCtSys, "mcnorm_fluxtuned_banfftuned_syst_multi");
#endif
}

//













//________________________________________________________________________________________
bool SuperkMC::ApplyEnergyScaleSystematicToTotalSpectrum(void)
{
    //
    // Now that all other systematics and physics hypotheses (that depend on Etrue) have been
    // applied, apply the reconstructed energy scale systematic (that depends on Ereco).
    // This systematic is applied by scaling the Ereco bin edges and assuming uniform
    // distribution of events within each bin.
#ifdef __show_update_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> Applying non-multiplicative systematics to total spectrum";
#endif
    
    bool applied_negevents_condition = false;
    
    // loop over samples
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        double escale = fSKDetFSI->EnergyScale(s);
        
        double * events_per_bin = this->CurrErecoSpectrumArray(s);
        
        util::ApplyEnergyScaleSystematic(s,events_per_bin,escale);
    }//is
    return applied_negevents_condition;
}
//________________________________________________________________________________________










// NEW 2D TEMPLATES

bool SuperkMC::ApplyEnergyScaleSystematicToTotalSpectrum2D(void)
{
    //
    // Now that all other systematics and physics hypotheses (that depend on Etrue) have been
    // applied, apply the reconstructed energy scale systematic (that depends on Ereco).
    // This systematic is applied by scaling the Ereco bin edges and assuming uniform
    // distribution of events within each bin.
#ifdef __show_update_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> Applying non-multiplicative systematics to total spectrum";
#endif
    
    bool applied_negevents_condition = false;
    
    // loop over samples
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        double escale = fSKDetFSI->EnergyScale(s);
        //double escale = fSKDetFSI->CurrentValue(conv::kNuisSuperkEnergyScale);
        
        double * events_per_bin = this->CurrErecoVsQ2SpectrumArray(s);
        
        util::ApplyEnergyScaleSystematic2D(s,events_per_bin,escale);
    }//is
    return applied_negevents_condition;
}
//









//________________________________________________________________________________________
void SuperkMC::ApplyPhysHypothesisToTemplates(void)
{
    // Apply physics hypothesis {oscillation (2-flavour / 3-flavour), decay or
    // decoherence} to MC templates
    //
    
#ifdef __show_update_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> Applying physics hypothesis";
#endif
    
    PhysHypothesis_t hyp = fPhysHypothesis.Hypothesis();
    
    //
    // 3-flavour oscillations
    //
    if(hyp == kHypOsc3F)
    {
        this->ApplyPhysHypothesisToTemplates_3fOsc(false); // NEW 2D TEMPLATES
    }
    else if(hyp == kHypOsc3FBar)
    {
        this->ApplyPhysHypothesisToTemplates_3fOscBar(false); // NEW 2D TEMPLATES
    }
    else if(hyp == kHypOsc3FBarDissap)
    {
      this->ApplyPhysHypothesisToTemplates_3fOscBarDissap(false); // NEW 2D TEMPLATES
    }
    else if(hyp == kHypOsc3FBarNormBar)
    {
        this->ApplyPhysHypothesisToTemplates_3fOscBarNormBar(false); // NEW 2D TEMPLATES
    }
    else if(hyp == kHypOsc3FNormBar)
    {
        this->ApplyPhysHypothesisToTemplates_3fOscNormBar(false); // NEW 2D TEMPLATES
    }
    else if(hyp == kHypOscSterile4F || hyp == kHypOscSterile5F)
    {
        this->ApplyPhysHypothesisToTemplates_Sterile(false); // NEW 2D TEMPLATES
    }
    else if (hyp == kHypVoid) {
        // Do not apply PhysHypothesis
    }
    else
    {
        LOG("SuperkMC",pERROR)<<"SKMC/ApplyPhysHypothesisToTemplates >> Unimplemented hypothesis: "<<util::AsString(hyp);
        util::Die(-1); return;
    }
    
#ifdef __save_templates_debug__
    this->NevDebug(kSKMCtFin, "mcnorm_fluxtuned_banfftuned_syst_phys");
#endif
    
    return;
}
//________________________________________________________________________________________







// NEW 2D TEMPLATES

void SuperkMC::ApplyPhysHypothesisToTemplates2D(void)
{
    // Apply physics hypothesis {oscillation (2-flavour / 3-flavour), decay or
    // decoherence} to MC templates
    //
    
#ifdef __show_update_mesg_info__
    LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> Applying physics hypothesis";
#endif
    
    PhysHypothesis_t hyp = fPhysHypothesis.Hypothesis();
    
    //
    // 3-flavour oscillations
    //
    if(hyp == kHypOsc3F)
    {
        this->ApplyPhysHypothesisToTemplates_3fOsc(true);
    }
    else if(hyp == kHypOsc3FBar)
    {
        this->ApplyPhysHypothesisToTemplates_3fOscBar(true);
    }
    else if(hyp == kHypOsc3FBarDissap)
    {
      this->ApplyPhysHypothesisToTemplates_3fOscBarDissap(true);
    }
    else if(hyp == kHypOsc3FBarNormBar)
    {
        this->ApplyPhysHypothesisToTemplates_3fOscBarNormBar(true);
    }
    else if(hyp == kHypOsc3FNormBar)
    {
        this->ApplyPhysHypothesisToTemplates_3fOscNormBar(true);
    }
    else if(hyp == kHypOscSterile4F || hyp == kHypOscSterile5F)
    {
        this->ApplyPhysHypothesisToTemplates_Sterile(true);
    }
    else if (hyp == kHypVoid) {
        // Do not apply PhysHypothesis
    }
    else
    {
        LOG("SuperkMC",pERROR)<<"SKMC/ApplyPhysHypothesisToTemplates2D >> Unimplemented hypothesis: "<<util::AsString(hyp);
        util::Die(-1); return;
    }
    
#ifdef __save_templates_debug__
    this->NevDebug2D(kSKMCtFin, "mcnorm_fluxtuned_banfftuned_syst_phys");
#endif
    
    return;
}
//



//________________________________________________________________________________________
//void SuperkMC::ApplyPhysHypothesisToTemplates_3fOsc(void)
void SuperkMC::ApplyPhysHypothesisToTemplates_3fOsc(const bool Use2Dtemplates) // NEW 2D TEMPLATES
{
#ifndef __ALLOW_DIFFERING_HIERARCHIES__
    assert(Opt::Instance()->gMassHierarchyAssumed==fNuOsc3Flavour.GetMassHierarchy());
#endif
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
    LOG("SuperkMC",pDEBUG) << "**  SKMC/ApplyPhysHypothesisToTemplates_3fOsc >> " << fPhysHypothesis;
#endif
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    NuOscDmsqScheme_t  dmsqscheme = Opt::Instance()->GetEffectiveDmsqScheme();
    
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
        LOG("SuperkMC",pERROR) << "**  Exp/SKMC >> Can't handle choice of "
        << "\n oscillation mixing params: " << util::AsString(omp)
        << "\n and independent squared mass-splittings: " << util::AsString(dmsqscheme);
        util::Die(1); return;
    }
    
    // loop over samples
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        int n_q2_bins    = kNQ2Bins[is]; // NEW 2D TEMPLATES
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
            
            // double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuMu,    conv::kPdgNuMu,    Etrue);
            // double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuE,     conv::kPdgNuE,     Etrue);
            // double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuMu,    conv::kPdgNuE,     Etrue);
            // double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuE,    conv::kPdgNuMu,     Etrue);
            // double p_numubar_numubar = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuMuBar, Etrue);
            // double p_nuebar_nuebar   = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuEBar,  conv::kPdgNuEBar,  Etrue);
            // double p_numubar_nuebar  = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuEBar,  Etrue);
            // double p_nuebar_numubar  = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuEBar, conv::kPdgNuMuBar,  Etrue);
            
            
#ifdef __show_hypothesis_probabilities__
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  E = " << Etrue << " GeV" << endl;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numu     -> numu   ) = " << p_numu_numu      ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numubar  -> numubar) = " << p_numubar_numubar;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nue      -> nue    ) = " << p_nue_nue        ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nuebar   -> nuebar ) = " << p_nuebar_nuebar  ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numu     -> nue    ) = " << p_numu_nue       ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numubar  -> nuebar ) = " << p_numubar_nuebar ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nue      -> numu   ) = " << p_nue_numu       ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nuebar   -> numubar) = " << p_nuebar_numubar ;
#endif
            
            
            
            
            // NEW 2D TEMPLATES
            
            if(Use2Dtemplates){ // 2D templates
                
                // loop over modes
                for(int im=0; im < knm; im++) {
                    SuperkMode_t m = kSKReactionModes[is][im];
                    if(m == kSKmUndefined) continue;
                    //double * array = &(fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,0,0));
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
                        // loop over Q2
                        for(int iq = 1; iq <= n_q2_bins ; iq++) {
                            
                            // number of events in (Q2,Ereco,Etrue) bin (systematics applied)
                            double nev = 0;
                            //nev = array[ir-1][iq-1];
                            nev = fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,ir-1,iq-1);
                            ////
                            
#ifdef PHYSHYPWEIGHTSAREONE
                            posc = 1;
#ifdef PHYSHYPWEIGHTAREONEANDNUEISZERO
                            if(util::IsOscNuE(m) || util::IsOscNuEBar(m) || util::IsOscNuMuBar(m) || util::IsOscNuMu(m) ) posc = 0;
#endif
#endif
                            // number of events in (Ereco,Etrue) bin (oscillations + systematics applied)
                            double nev_osc = nev * posc;
                            //array[ir-1][iq-1] = nev_osc;
                            fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,ir-1,iq-1) = nev_osc;
                            ////
                        }//iq
                    }//ir
                }//im
                
            }
            else{ // 1D templates
                
                // loop over modes
                for(int im=0; im < knm; im++) {
                    SuperkMode_t m = kSKReactionModes[is][im];
                    if(m == kSKmUndefined) continue;
                    double * array = &(fMCTemplArrayFin[is]->GetEntryRef(im,it-1,0));
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
                
            } // end if Use2Dtemplates)
            //
            
        }//it
    }//is
    
    return;
}

//________________________________________________________________________________________
//void SuperkMC::ApplyPhysHypothesisToTemplates_3fOscBar(void)
void SuperkMC::ApplyPhysHypothesisToTemplates_3fOscBar(const bool Use2Dtemplates) // NEW 2D TEMPLATES
{
    //
    // 3-flavour oscillations, separate nu and nubar
    //
    
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
    LOG("SuperkMC",pDEBUG) << "**  Exp/SKMC >> " << fPhysHypothesis;
#endif
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    NuOscDmsqScheme_t  dmsqscheme = Opt::Instance()->GetEffectiveDmsqScheme();
    
    double p_osc[3][3];
    double p_osc_bar[3][3];
    
    // mixing params are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|
    if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsq2Theta)
    {
        fNuOsc3Flavour   .SetOscParameters(mixing12,   mixing13,   mixing23,   dmsq   ,Dmsq   ,delta_cp);
        fNuOsc3FlavourBar.SetOscParameters(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing params are: sin^2(theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|
    else if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsqTheta)
    {
        fNuOsc3Flavour   .SetOscParametersSinSqTAll(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
        fNuOsc3FlavourBar.SetOscParametersSinSqTAll(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing params are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|
    else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsq2Theta)
    {
        fNuOsc3Flavour   .SetOscParametersSSq2TDmAtmSol(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
        fNuOsc3FlavourBar.SetOscParametersSSq2TDmAtmSol(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing params are: sin^2(theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|
    else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsqTheta)
    {
        fNuOsc3Flavour   .SetOscParametersDmsqAtmSol(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
        fNuOsc3FlavourBar.SetOscParametersDmsqAtmSol(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing paras are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
    else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsq2Theta)
    {
        fNuOsc3Flavour   .SetOscParametersSSq2TDmsq13(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
        fNuOsc3FlavourBar.SetOscParametersSSq2TDmsq13(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing paras are: sin^2(theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
    else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsqTheta)
    {
        fNuOsc3Flavour   .SetOscParametersDmsq13(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
        fNuOsc3FlavourBar.SetOscParametersDmsq13(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    //unknown
    else {
        LOG("SuperkMC",pERROR) << "**  Exp/SKMC >> Can't handle choice of "
        << "\n oscillation mixing params: " << util::AsString(omp)
        << "\n and independent squared mass-splittings: " << util::AsString(dmsqscheme);
        util::Die(1); return;
    }
    
    // loop over samples
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        int n_q2_bins    = kNQ2Bins[is]; // NEW 2D TEMPLATES
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
            
            // double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuMu,    conv::kPdgNuMu,    Etrue);
            // double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuE,     conv::kPdgNuE,     Etrue);
            // double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuMu,    conv::kPdgNuE,     Etrue);
            // double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuE,    conv::kPdgNuMu,     Etrue);
            // double p_numubar_numubar = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuMuBar, Etrue);
            // double p_nuebar_nuebar   = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuEBar,  conv::kPdgNuEBar,  Etrue);
            // double p_numubar_nuebar  = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuEBar,  Etrue);
            // double p_nuebar_numubar  = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuEBar, conv::kPdgNuMuBar,  Etrue);
            
            
#ifdef __show_hypothesis_probabilities__
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  E = " << Etrue << " GeV" << endl;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numu     -> numu   ) = " << p_numu_numu      ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numubar  -> numubar) = " << p_numubar_numubar;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nue      -> nue    ) = " << p_nue_nue        ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nuebar   -> nuebar ) = " << p_nuebar_nuebar  ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numu     -> nue    ) = " << p_numu_nue       ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numubar  -> nuebar ) = " << p_numubar_nuebar ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nue      -> numu   ) = " << p_nue_numu       ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nuebar   -> numubar) = " << p_nuebar_numubar ;
#endif
            
            
            
            // NEW 2D TEMPLATES
            
            if(Use2Dtemplates){ // 2D templates
                
                // loop over modes
                for(int im=0; im < knm; im++) {
                    SuperkMode_t m = kSKReactionModes[is][im];
                    if(m == kSKmUndefined) continue;
                    //double * array = &(fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,0,0));
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
                        // loop over Q2
                        for(int iq = 1; iq <= n_q2_bins ; iq++) {
                            // number of events in (Q2,Ereco,Etrue) bin (systematics applied)
                            double nev = 0;
                            //nev = array[ir-1][iq-1];
                            nev = fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,ir-1,iq-1);
                            ////
                            
#ifdef PHYSHYPWEIGHTSAREONE
                            posc = 1;
#ifdef PHYSHYPWEIGHTAREONEANDNUEISZERO
                            if(util::IsOscNuE(m) || util::IsOscNuEBar(m) || util::IsOscNuMuBar(m) || util::IsOscNuMu(m) ) posc = 0;
#endif
#endif
                            // number of events in (Q2,Ereco,Etrue) bin (oscillations + systematics applied)
                            double nev_osc = nev * posc;
                            //array[ir-1][iq-1] = nev_osc;
                            fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,ir-1,iq-1) = nev_osc;
                            ////
                        }//iq
                    }//ir
                }//im
                
            }
            else{ // 1D templates
                
                // loop over modes
                for(int im=0; im < knm; im++) {
                    SuperkMode_t m = kSKReactionModes[is][im];
                    if(m == kSKmUndefined) continue;
                    double * array = &(fMCTemplArrayFin[is]->GetEntryRef(im,it-1,0));
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
                
            } // end if Use2DTemplates
            //
            
        }//it
    }//is
    return;
}
//________________________________________________________________________________________
void SuperkMC::ApplyPhysHypothesisToTemplates_3fOscBarDissap(const bool Use2Dtemplates) // NEW 2D TEMPLATES
{
  //
  // 3-flavour oscillations, separate nu and nubar
  //
  
  // set oscillation parameters to the oscillation probability calculation code
  double mixing12 = fPhysHypothesis.CurrValue(kPhypMixing12);
  double mixing13 = fPhysHypothesis.CurrValue(kPhypMixing13);
  double mixing23 = fPhysHypothesis.CurrValue(kPhypMixing23);
  //double dm12sq   = fPhysHypothesis.CurrValue(kPhypDm12sq);
  //double dm23sq   = fPhysHypothesis.CurrValue(kPhypDm23sq);
  double dmsq     = fPhysHypothesis.CurrValue(kPhypDm12sq);
  double Dmsq     = fPhysHypothesis.CurrValue(kPhypDm23sq);
  double delta_cp = fPhysHypothesis.CurrValue(kPhypCPV    );
  
  double mixing23bar = fPhysHypothesis.CurrValue(kPhypMixing23Bar);
  double Dmsqbar     = fPhysHypothesis.CurrValue(kPhypDm23sqBar);
  
#ifdef __show_hypothesis__
  LOG("SuperkMC",pDEBUG) << "**  Exp/SKMC >> " << fPhysHypothesis;
#endif
  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  NuOscDmsqScheme_t  dmsqscheme = Opt::Instance()->GetEffectiveDmsqScheme();
  
  double p_osc[3][3];
  double p_osc_bar[3][3];
  
  // mixing params are: sin^2(2*theta)
  // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|
  if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsq2Theta)
  {
    fNuOsc3Flavour   .SetOscParameters(mixing12,   mixing13,   mixing23,   dmsq   ,Dmsq   ,delta_cp);
    fNuOsc3FlavourBar.SetOscParameters(mixing12,mixing13,mixing23bar,dmsq,Dmsqbar,delta_cp);
  }
  // mixing params are: sin^2(theta)
  // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|
  else if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsqTheta)
  {
    fNuOsc3Flavour   .SetOscParametersSinSqTAll(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
    fNuOsc3FlavourBar.SetOscParametersSinSqTAll(mixing12,mixing13,mixing23bar,dmsq,Dmsqbar,delta_cp);
  }
  // mixing params are: sin^2(2*theta)
  // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|
  else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsq2Theta)
  {
    fNuOsc3Flavour   .SetOscParametersSSq2TDmAtmSol(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
    fNuOsc3FlavourBar.SetOscParametersSSq2TDmAtmSol(mixing12,mixing13,mixing23bar,dmsq,Dmsqbar,delta_cp);
  }
  // mixing params are: sin^2(theta)
  // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|
  else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsqTheta)
  {
    fNuOsc3Flavour   .SetOscParametersDmsqAtmSol(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
    fNuOsc3FlavourBar.SetOscParametersDmsqAtmSol(mixing12,mixing13,mixing23bar,dmsq,Dmsqbar,delta_cp);
  }
  // mixing paras are: sin^2(2*theta)
  // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
  else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsq2Theta)
  {
    fNuOsc3Flavour   .SetOscParametersSSq2TDmsq13(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
    fNuOsc3FlavourBar.SetOscParametersSSq2TDmsq13(mixing12,mixing13,mixing23bar,dmsq,Dmsqbar,delta_cp);
  }
  // mixing paras are: sin^2(theta)
  // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
  else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsqTheta)
  {
    fNuOsc3Flavour   .SetOscParametersDmsq13(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
    fNuOsc3FlavourBar.SetOscParametersDmsq13(mixing12,mixing13,mixing23bar,dmsq,Dmsqbar,delta_cp);
  }
  //unknown
  else {
    LOG("SuperkMC",pERROR) << "**  Exp/SKMC >> Can't handle choice of "
    << "\n oscillation mixing params: " << util::AsString(omp)
    << "\n and independent squared mass-splittings: " << util::AsString(dmsqscheme);
    util::Die(1); return;
  }
  
  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    
    int n_q2_bins    = kNQ2Bins[is]; // NEW 2D TEMPLATES
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
      
      // double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuMu,    conv::kPdgNuMu,    Etrue);
      // double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuE,     conv::kPdgNuE,     Etrue);
      // double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuMu,    conv::kPdgNuE,     Etrue);
      // double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuE,    conv::kPdgNuMu,     Etrue);
      // double p_numubar_numubar = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuMuBar, Etrue);
      // double p_nuebar_nuebar   = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuEBar,  conv::kPdgNuEBar,  Etrue);
      // double p_numubar_nuebar  = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuEBar,  Etrue);
      // double p_nuebar_numubar  = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuEBar, conv::kPdgNuMuBar,  Etrue);
      
      
#ifdef __show_hypothesis_probabilities__
      LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  E = " << Etrue << " GeV" ;
      LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numu     -> numu   ) = " << p_numu_numu       ;
      LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numubar  -> numubar) = " << p_numubar_numubar ;
      LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nue      -> nue    ) = " << p_nue_nue         ;
      LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nuebar   -> nuebar ) = " << p_nuebar_nuebar   ;
      LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numu     -> nue    ) = " << p_numu_nue        ;
      LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numubar  -> nuebar ) = " << p_numubar_nuebar  ;
      LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nue      -> numu   ) = " << p_nue_numu        ;
      LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nuebar   -> numubar) = " << p_nuebar_numubar  ;
#endif
      
      
      
      // NEW 2D TEMPLATES
      
      if(Use2Dtemplates){ // 2D templates
        
        // loop over modes
        for(int im=0; im < knm; im++) {
          SuperkMode_t m = kSKReactionModes[is][im];
          if(m == kSKmUndefined) continue;
          //double * array = &(fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,0,0));
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
            // loop over Q2
            for(int iq = 1; iq <= n_q2_bins ; iq++) {
              // number of events in (Q2,Ereco,Etrue) bin (systematics applied)
              double nev = 0;
              //nev = array[ir-1][iq-1];
              nev = fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,ir-1,iq-1);
              ////
              
#ifdef PHYSHYPWEIGHTSAREONE
              posc = 1;
#ifdef PHYSHYPWEIGHTAREONEANDNUEISZERO
              if(util::IsOscNuE(m) || util::IsOscNuEBar(m) || util::IsOscNuMuBar(m) || util::IsOscNuMu(m) ) posc = 0;
#endif
#endif
              // number of events in (Q2,Ereco,Etrue) bin (oscillations + systematics applied)
              double nev_osc = nev * posc;
              //array[ir-1][iq-1] = nev_osc;
              fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,ir-1,iq-1) = nev_osc;
              ////
            }//iq
          }//ir
        }//im
        
      }
      else{ // 1D templates
        
        // loop over modes
        for(int im=0; im < knm; im++) {
          SuperkMode_t m = kSKReactionModes[is][im];
          if(m == kSKmUndefined) continue;
          double * array = &(fMCTemplArrayFin[is]->GetEntryRef(im,it-1,0));
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
        
      } // end if Use2DTemplates
        //
      
    }//it
  }//is
  return;
}
//________________________________________________________________________________________
//void SuperkMC::ApplyPhysHypothesisToTemplates_3fOscBarNormBar(void)
void SuperkMC::ApplyPhysHypothesisToTemplates_3fOscBarNormBar(const bool Use2Dtemplates) // NEW 2D TEMPLATES
{
    //
    // 3-flavour oscillations, separate nu and nubar
    //
  
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
    
    double betabar     = fPhysHypothesis.CurrValue(kPhypBetaBar);  // add betabar
    
#ifdef __show_hypothesis__
    LOG("SuperkMC",pDEBUG) << "**  Exp/SKMC >> " << fPhysHypothesis;
#endif
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    NuOscDmsqScheme_t  dmsqscheme = Opt::Instance()->GetEffectiveDmsqScheme();
    
    double p_osc[3][3];
    double p_osc_bar[3][3];
    
    // mixing params are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|
    if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsq2Theta)
    {
        fNuOsc3Flavour   .SetOscParameters(mixing12,   mixing13,   mixing23,   dmsq   ,Dmsq   ,delta_cp);
        fNuOsc3FlavourBar.SetOscParameters(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing params are: sin^2(theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|
    else if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsqTheta)
    {
        fNuOsc3Flavour   .SetOscParametersSinSqTAll(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
        fNuOsc3FlavourBar.SetOscParametersSinSqTAll(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing params are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|
    else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsq2Theta)
    {
        fNuOsc3Flavour   .SetOscParametersSSq2TDmAtmSol(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
        fNuOsc3FlavourBar.SetOscParametersSSq2TDmAtmSol(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing params are: sin^2(theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|
    else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsqTheta)
    {
        fNuOsc3Flavour   .SetOscParametersDmsqAtmSol(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
        fNuOsc3FlavourBar.SetOscParametersDmsqAtmSol(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing paras are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
    else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsq2Theta)
    {
        fNuOsc3Flavour   .SetOscParametersSSq2TDmsq13(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
        fNuOsc3FlavourBar.SetOscParametersSSq2TDmsq13(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    // mixing paras are: sin^2(theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
    else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsqTheta)
    {
        fNuOsc3Flavour   .SetOscParametersDmsq13(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
        fNuOsc3FlavourBar.SetOscParametersDmsq13(mixing12bar,mixing13bar,mixing23bar,dmsqbar,Dmsqbar,delta_cpbar);
    }
    //unknown
    else {
        LOG("SuperkMC",pERROR) << "**  Exp/SKMC >> Can't handle choice of "
        << "\n oscillation mixing params: " << util::AsString(omp)
        << "\n and independent squared mass-splittings: " << util::AsString(dmsqscheme);
        util::Die(1); return;
    }
    
    // loop over samples
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        int n_q2_bins    = kNQ2Bins[is]; // NEW 2D TEMPLATES
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
            
            // double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuMu,    conv::kPdgNuMu,    Etrue);
            // double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuE,     conv::kPdgNuE,     Etrue);
            // double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuMu,    conv::kPdgNuE,     Etrue);
            // double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuE,    conv::kPdgNuMu,     Etrue);
            // double p_numubar_numubar = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuMuBar, Etrue);
            // double p_nuebar_nuebar   = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuEBar,  conv::kPdgNuEBar,  Etrue);
            // double p_numubar_nuebar  = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuEBar,  Etrue);
            // double p_nuebar_numubar  = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuEBar, conv::kPdgNuMuBar,  Etrue);
            
            
#ifdef __show_hypothesis_probabilities__
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  E = " << Etrue << " GeV" ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numu     -> numu   ) = " << p_numu_numu       ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numubar  -> numubar) = " << p_numubar_numubar ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nue      -> nue    ) = " << p_nue_nue         ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nuebar   -> nuebar ) = " << p_nuebar_nuebar   ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numu     -> nue    ) = " << p_numu_nue        ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numubar  -> nuebar ) = " << p_numubar_nuebar  ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nue      -> numu   ) = " << p_nue_numu        ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nuebar   -> numubar) = " << p_nuebar_numubar  ;
#endif
            
            
            
            
            
            // NEW 2D TEMPLATES
            
            if(Use2Dtemplates){ // 2D templates
                // loop over modes
                for(int im=0; im < knm; im++) {
                    SuperkMode_t m = kSKReactionModes[is][im];
                    if(m == kSKmUndefined) continue;
                    //double * array = &(fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,0,0));
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
                    
                    if(util::IsOscNuEBar(m)) posc *= betabar;   // Enable or disable nuebar appearance
                    
                    // loop over Ereco
                    for(int ir = 1; ir <= n_ereco_bins ; ir++) {
                        // loop over Q2
                        for(int iq = 1; iq <= n_q2_bins ; iq++) {
                            // number of events in (Q2,Ereco,Etrue) bin (systematics applied)
                            double nev = 0;
                            //nev = array[ir-1][iq-1];
                            nev = fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,ir-1,iq-1);
                            ////
                            
#ifdef PHYSHYPWEIGHTSAREONE
                            posc = 1;
#ifdef PHYSHYPWEIGHTAREONEANDNUEISZERO
                            if(util::IsOscNuE(m) || util::IsOscNuEBar(m) || util::IsOscNuMuBar(m) || util::IsOscNuMu(m) ) posc = 0;
#endif
#endif
                            // number of events in (Q2,Ereco,Etrue) bin (oscillations + systematics applied)
                            double nev_osc = nev * posc;
                            //array[ir-1][iq-1] = nev_osc;
                            fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,ir-1,iq-1) = nev_osc;
                            ////
                        }//iq
                    }//ir
                }//im
                
            }
            else{ // 1D templates
                
                // loop over modes
                for(int im=0; im < knm; im++) {
                    SuperkMode_t m = kSKReactionModes[is][im];
                    if(m == kSKmUndefined) continue;
                    double * array = &(fMCTemplArrayFin[is]->GetEntryRef(im,it-1,0));
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
                    
                    if(util::IsOscNuEBar(m)) posc *= betabar;   // Enable or disable nuebar appearance
                    
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
            } // end if Use2Dtemplates
            
            //
            
            
        }//it
    }//is
    return;
}
//________________________________________________________________________________________
//void SuperkMC::ApplyPhysHypothesisToTemplates_3fOscNormBar(void)
void SuperkMC::ApplyPhysHypothesisToTemplates_3fOscNormBar(const bool Use2Dtemplates) // NEW 2D TEMPLATES
{
    //
    // 3-flavour oscillations + norm nuebar
    //
    
    // set oscillation parameters to the oscillation probability calculation code
    double mixing12 = fPhysHypothesis.CurrValue(kPhypMixing12);
    double mixing13 = fPhysHypothesis.CurrValue(kPhypMixing13);
    double mixing23 = fPhysHypothesis.CurrValue(kPhypMixing23);
    //double dm12sq   = fPhysHypothesis.CurrValue(kPhypDm12sq);
    //double dm23sq   = fPhysHypothesis.CurrValue(kPhypDm23sq);
    double dmsq     = fPhysHypothesis.CurrValue(kPhypDm12sq);
    double Dmsq     = fPhysHypothesis.CurrValue(kPhypDm23sq);
    double delta_cp = fPhysHypothesis.CurrValue(kPhypCPV    );
    
    double betabar     = fPhysHypothesis.CurrValue(kPhypBetaBar);  // add betabar
    
#ifdef __show_hypothesis__
    LOG("SuperkMC",pDEBUG) << "**  Exp/SKMC >> " << fPhysHypothesis;
#endif
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    NuOscDmsqScheme_t  dmsqscheme = Opt::Instance()->GetEffectiveDmsqScheme();
    
    double p_osc[3][3];
    double p_osc_bar[3][3];
    
    // mixing params are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|
    if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsq2Theta)
    {
        fNuOsc3Flavour   .SetOscParameters(mixing12,   mixing13,   mixing23,   dmsq   ,Dmsq   ,delta_cp);
    }
    // mixing params are: sin^2(theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{23}|
    else if(dmsqscheme == kNuDmsq12andAbs23 && omp == kNuMixSinsqTheta)
    {
        fNuOsc3Flavour   .SetOscParametersSinSqTAll(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
    }
    // mixing params are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|
    else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsq2Theta)
    {
        fNuOsc3Flavour   .SetOscParametersSSq2TDmAtmSol(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
    }
    // mixing params are: sin^2(theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|
    else if(dmsqscheme == kNuDmsqFogliLisi && omp == kNuMixSinsqTheta)
    {
        fNuOsc3Flavour   .SetOscParametersDmsqAtmSol(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
    }
    // mixing paras are: sin^2(2*theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
    else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsq2Theta)
    {
        fNuOsc3Flavour   .SetOscParametersSSq2TDmsq13(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
    }
    // mixing paras are: sin^2(theta)
    // independent squared mass splittings are: dm^2 = Dm^2_{12}, Dm^2 = |Dm^2_{31}|
    else if(dmsqscheme == kNuDmsq12andAbs31 && omp == kNuMixSinsqTheta)
    {
        fNuOsc3Flavour   .SetOscParametersDmsq13(mixing12,   mixing13,   mixing23,   dmsq,   Dmsq,   delta_cp);
    }
    //unknown
    else {
        LOG("SuperkMC",pERROR) << "**  Exp/SKMC >> Can't handle choice of "
        << "\n oscillation mixing params: " << util::AsString(omp)
        << "\n and independent squared mass-splittings: " << util::AsString(dmsqscheme);
        util::Die(1); return;
    }
    
    // loop over samples
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        int n_q2_bins    = kNQ2Bins[is]; // NEW 2D TEMPLATES
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
            
            // double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuMu,    conv::kPdgNuMu,    Etrue);
            // double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuE,     conv::kPdgNuE,     Etrue);
            // double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuMu,    conv::kPdgNuE,     Etrue);
            // double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb   (conv::kPdgNuE,    conv::kPdgNuMu,     Etrue);
            // double p_numubar_numubar = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuMuBar, Etrue);
            // double p_nuebar_nuebar   = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuEBar,  conv::kPdgNuEBar,  Etrue);
            // double p_numubar_nuebar  = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuEBar,  Etrue);
            // double p_nuebar_numubar  = fNuOsc3FlavourBar.CalcMatterOscProb(conv::kPdgNuEBar, conv::kPdgNuMuBar,  Etrue);
            
            
#ifdef __show_hypothesis_probabilities__
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  E = " << Etrue << " GeV" ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numu     -> numu   ) = " << p_numu_numu       ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numubar  -> numubar) = " << p_numubar_numubar ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nue      -> nue    ) = " << p_nue_nue         ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nuebar   -> nuebar ) = " << p_nuebar_nuebar   ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numu     -> nue    ) = " << p_numu_nue        ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numubar  -> nuebar ) = " << p_numubar_nuebar  ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nue      -> numu   ) = " << p_nue_numu        ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nuebar   -> numubar) = " << p_nuebar_numubar  ;
#endif
            
            
            
            // NEW 2D TEMPLATES
            if(Use2Dtemplates){ // 2D templates
                
                // loop over modes
                for(int im=0; im < knm; im++) {
                    SuperkMode_t m = kSKReactionModes[is][im];
                    if(m == kSKmUndefined) continue;
                    //double * array = &(fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,0,0));
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
                    
                    if(util::IsOscNuEBar(m)) posc *= betabar;   // Enable or disable nuebar appearance
                    
                    // loop over Ereco
                    for(int ir = 1; ir <= n_ereco_bins ; ir++) {
                        // loop over Q2
                        for(int iq = 1; iq <= n_q2_bins ; iq++) {
                            // number of events in (Q2,Ereco,Etrue) bin (systematics applied)
                            double nev = 0;
                            //nev = array[ir-1][iq-1];
                            nev = fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,ir-1,iq-1);
                            ////
                            
#ifdef PHYSHYPWEIGHTSAREONE
                            posc = 1;
#ifdef PHYSHYPWEIGHTAREONEANDNUEISZERO
                            if(util::IsOscNuE(m) || util::IsOscNuEBar(m) || util::IsOscNuMuBar(m) || util::IsOscNuMu(m) ) posc = 0;
#endif
#endif
                            // number of events in (Q2,Ereco,Etrue) bin (oscillations + systematics applied)
                            double nev_osc = nev * posc;
                            //array[ir-1][iq-1] = nev_osc;
                            fMCTemplArrayFin2D[is]->GetEntryRef(im,it-1,ir-1,iq-1) = nev_osc;
                            ////
                        }//iq
                    }//ir
                }//im
                
            }
            else{ // 1D templates
                
                // loop over modes
                for(int im=0; im < knm; im++) {
                    SuperkMode_t m = kSKReactionModes[is][im];
                    if(m == kSKmUndefined) continue;
                    double * array = &(fMCTemplArrayFin[is]->GetEntryRef(im,it-1,0));
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
                    
                    if(util::IsOscNuEBar(m)) posc *= betabar;   // Enable or disable nuebar appearance
                    
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
            } // end if Use2Dtemplates
            
            //
            
        }//it
    }//is
    return;
}
//________________________________________________________________________________________
//void SuperkMC::ApplyPhysHypothesisToTemplates_Sterile(void)
void SuperkMC::ApplyPhysHypothesisToTemplates_Sterile(const bool Use2Dtemplates) // NEW 2D TEMPLATES
{
#ifndef __ALLOW_DIFFERING_HIERARCHIES__
    assert(Opt::Instance()->gMassHierarchyAssumed==fNuOscSterile.GetMassHierarchy());
#endif
    //
    // 3+n flavour sterile oscillations
    //
    
    PhysHypothesis_t hyp = fPhysHypothesis.Hypothesis();
    
    //NuOscSterile osc;
    // Use fNuOscSterile
    
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
    LOG("SuperkMC",pDEBUG) << "**  Exp/SKMC >> " << fPhysHypothesis;
#endif
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    
    if (omp == kNuMixSinsqTheta) {
        fNuOscSterile.SetMixingAngles(mixing12,mixing13,mixing14,mixing15,mixing23,mixing24,mixing25,mixing34,mixing35);
    }
    else if (omp == kNuMixSinsq2Theta) {
        if (Opt::Instance()->gOscOctantChoice == conv::kNuOctantOne) {
            fNuOscSterile.SetMixingDoubleAnglesOct1(mixing12,mixing13,mixing14,mixing15,mixing23,mixing24,mixing25,mixing34,mixing35);
        }
        else if (Opt::Instance()->gOscOctantChoice == conv::kNuOctantTwo) {
            fNuOscSterile.SetMixingDoubleAnglesOct2(mixing12,mixing13,mixing14,mixing15,mixing23,mixing24,mixing25,mixing34,mixing35);
        }
        else {
            LOG("SuperkMC",pERROR) << "**  Exp/SKMC >> Can't handle choice of "
            << "\n octant: " << util::AsString(Opt::Instance()->gOscOctantChoice);
            util::Die(1); return;
        }
    }
    else {
        LOG("SuperkMC",pERROR) << "**  Exp/SKMC >> Can't handle choice of "
        << "\n oscillation mixing params: " << util::AsString(omp);
        util::Die(1); return;
    }
    
    // May need to put some decision on mass splitting convention here.
    //~ NuOscDmsqScheme_t  dmsqscheme = Opt::Instance()->GetEffectiveDmsqScheme();
    
    fNuOscSterile.SetDmsq(dmsq, Dmsq, dm41sq, dm51sq);
    fNuOscSterile.SetDeltaCP(delta_12, delta_13, delta_15, delta_34, delta_35);
    
    // loop over samples
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        int n_q2_bins    = kNQ2Bins[is]; // NEW 2D TEMPLATES
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
            
            fNuOscSterile.CalcMatterOscProbAll(false, Etrue, &p_osc    [0][0]);
            fNuOscSterile.CalcMatterOscProbAll(true , Etrue, &p_osc_bar[0][0]);
            
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
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  E = " << Etrue << " GeV" << endl;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numu     -> numu   ) = " << p_numu_numu       ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numubar  -> numubar) = " << p_numubar_numubar ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nue      -> nue    ) = " << p_nue_nue         ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nuebar   -> nuebar ) = " << p_nuebar_nuebar   ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numu     -> nue    ) = " << p_numu_nue        ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numubar  -> nuebar ) = " << p_numubar_nuebar  ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nue      -> numu   ) = " << p_nue_numu        ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nuebar   -> numubar) = " << p_nuebar_numubar  ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numu     -> active ) = " << p_numu_nc         ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (numubar  -> active ) = " << p_numubar_nc      ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nue      -> active ) = " << p_nue_nc          ;
            LOG("SuperkMC",pDEBUG) << "Exp/SKMC >>  * P (nuebar   -> active ) = " << p_nuebar_nc       ;
            
#endif
            
            
            
            // NEW 2D TEMPLATES
            
            if(Use2Dtemplates){ // 2D templates
                
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
                        // loop over Q2
                        for(int iq = 1; iq <= n_q2_bins ; iq++) {
                            // number of events in (Q2,Ereco,Etrue) bin (systematics applied)
                            double nev = 0;
                            nev = fMCTemplArrayFin2D[is]->GetEntryRef(im, it-1, ir-1, iq-1);
                            ////
                            
#ifdef PHYSHYPWEIGHTSAREONE
                            posc = 1;
#ifdef PHYSHYPWEIGHTAREONEANDNUEISZERO
                            if(util::IsOscNuE(m)) posc = 0;
#endif
#endif
                            // number of events in (Q2,Ereco,Etrue) bin (oscillations + systematics applied)
                            double nev_osc = nev * posc;
                            fMCTemplArrayFin2D[is]->GetEntryRef(im, it-1, ir-1, iq-1) = nev_osc;
                            ////
                        }//iq
                    }//ir
                }//im
                
            }
            else{ // 1D templates
                
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
                        
                        nev = fMCTemplArrayFin[is]->GetEntryRef(im, it-1, ir-1);
                        ////
                        
#ifdef PHYSHYPWEIGHTSAREONE
                        posc = 1;
#ifdef PHYSHYPWEIGHTAREONEANDNUEISZERO
                        if(util::IsOscNuE(m)) posc = 0;
#endif
#endif
                        // number of events in (Ereco,Etrue) bin (oscillations + systematics applied)
                        double nev_osc = nev * posc;
                        
                        fMCTemplArrayFin[is]->GetEntryRef(im, it-1, ir-1) = nev_osc;
                        ////
                    }//ir
                }//im
            } // end if Use2Dtemplates
            
            //
            
        }//it
    }//is
    return;
}
//________________________________________________________________________________________
/*
 void SuperkMC::ApplyPhysHypothesisToTemplates_Other(void)
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
        const templmtype & templFin = *reinterpret_cast< const templmtype * > (fMCTemplArrayFin[is]->GetArray());
        ermtype & erecoMode         = *reinterpret_cast< ermtype          * > (fMCErecoSpectrumArray[is]->GetArray());
        
        //~ const double(&templFin)[knm][n_etrue_bins][n_ereco_bins] = *reinterpret_cast< const double (*) [knm][n_etrue_bins][n_ereco_bins]> (fMCTemplArrayFin[is]->GetArray());
        //~ double(&erecoMode)[knm][n_ereco_bins] = *reinterpret_cast<double (*) [knm][n_ereco_bins]> (fMCErecoSpectrumArray[is]->GetArray());
        
        double * erecoTot = fMCErecoSpectrumTotArray[is]->GetArray();
        
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
                    LOG("SuperkMC",pDEBUG) << "Exp/SKMC/UpdateErecoSpectra >> set negative number of events " << nev
                    << " to zero " << 0 << "\t" << " Ereco: " << reco_spectrum_tot->GetXaxis()->GetBinCenter(ir) ;
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






// NEW 2D TEMPLATES

bool SuperkMC::UpdateErecoVsQ2Spectra(void)
{
    bool applied_negevents_condition = false;
    // loop over samples
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        int n_q2_bins    = kNQ2Bins[is];
        int n_ereco_bins = kNErecoBins[is];
        int n_etrue_bins = kNEtrueBins[is];
        
        // Recast VALORArray blocks so they can be used addressed like normal arrays (fast).
        // Hard to read, so should only be used in code for which speed is essential.
        typedef double (templmtype)[knm][n_etrue_bins][n_ereco_bins][n_q2_bins];
        typedef double (erq2mtype) [knm][n_ereco_bins][n_q2_bins];
        typedef double (erq2type)  [n_ereco_bins][n_q2_bins];
        const templmtype & templFin = *reinterpret_cast< const templmtype * > (fMCTemplArrayFin2D[is]->GetArray());
        erq2mtype & erecoVsq2Mode   = *reinterpret_cast< erq2mtype        * > (fMCErecoVsQ2SpectrumArray[is]->GetArray());
        erq2type & erecoVsq2Tot     = *reinterpret_cast< erq2type         * > (fMCErecoVsQ2SpectrumTotArray2D[is]->GetArray());
        
        // loop over Q2
        for(int iq = 1; iq <= n_q2_bins; iq++) {
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
                        nev += templFin[im][it-1][ir-1][iq-1];
                    }//it
                    nev_tot += nev;
#ifdef __apply_minevents_to_individual_modes__
                    if(nev < 0) {
#ifdef __show_negative_events_mesg_debug__
                        LOG("SuperkMC",pDEBUG) << "Exp/SKMC/UpdateErecoSpectra >> set negative number of events " << nev
                        << " to zero " << 0 << "\t" << " Ereco: " << reco_spectrum_tot->GetYaxis()->GetBinCenter(ir)
                        << " to zero " << 0 << "\t" << " Q2: "    << q2_spectrum_tot->GetXaxis()->GetBinCenter(iq) ;
#endif
                        nev = 0;
                    } // if nev < 0
#endif //__apply_minevents_to_individual_modes__
                    
                    erecoVsq2Mode[im][ir-1][iq-1] = nev;
                    
                }//im
                
                erecoVsq2Tot[ir-1][iq-1] = nev_tot;
                
            }//iq
        }//ir
    }//is
    
    applied_negevents_condition = this->ApplyEnergyScaleSystematicToTotalSpectrum2D();
    
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
            
            fMCEtrueSpectrumTotArray[is]->GetEntryRef(it-1) = nev_tot;
            
        }//it
    }//is
}
//________________________________________________________________________________________
// NEW 2D TEMPLATES

void SuperkMC::FillEtrueSpectra2D(void)
{
    // loop over samples
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        TH1D * true_spectrum_tot = fMCEtrueSpectrumTot[is];
        
        int n_q2_bins = kNQ2Bins[is];
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
                
                TH3D * mctempl = this -> Template2D(s,m,kSKMCtOri);
                
                double nev = 0;
                // loop over Ereco
                for(int ir = 1; ir <= n_ereco_bins; ir++) {
                    // loop over Q2
                    for(int iq = 1; iq <= n_q2_bins; iq++) {
                        nev += mctempl->GetBinContent(iq,ir,it);
                    }//ir
                }//iq
                
                nev_tot += nev;
            }//im
            
            fMCEtrueSpectrumTotArray[is]->GetEntryRef(it-1) = nev_tot; // Use the same for 1D and 2D
            
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
    
    ev = fMCEtrueSpectrumTotArray[is]->GetEntryRef(it-1) > kASmallNum;
    ////
    
#ifdef __show_isanevent_mesg_debug__
    if(ev){
        
        LOG("SuperkMC",pDEBUG) << "Exp/SuperkMC >> There are " << fMCEtrueSpectrumTotArray[is]->GetEntryRef(it-1);
        LOG("SuperkMC",pDEBUG) << " events for sample = " << util::AsString(s)
        << ", Etrue = " << Etrue;
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
        LOG("SuperkMC",pERROR) << "Sample not selected";
        util::Die(1); return 0;
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
    
    //~ bool mode_is_nuesig    = util::IsOscNuE(m);
    //~ bool mode_is_nuebarsig = util::IsOscNuEBar(m);
    //~
    //~
    //~ PhysHypothesis_t hyp = fPhysHypothesis.Hypothesis();
    //~ bool hyp_is_not_3flavour_osc = (hyp != kHypOsc3F);
    //~
    //~ bool skip = ((mode_is_nuesig || mode_is_nuebarsig) && hyp_is_not_3flavour_osc );
    
    bool skip = false;
    
#ifdef __show_skipped_templates__
    if(skip) {
        LOG("SuperkMC",pDEBUG) << "Exp/SKMC >> Skipping template: " << util::AsString(m);
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
        const int erecobins = util::NumErecoBins(s);//70;
        const double * erecobinedges = util::ErecoBinEdges(s);
        
        //if(s == kSKsFCFV1RingELike) continue;
        // if(s == kSKsFCFV1RingMuLike) continue;
        // if(s == kSKsRHCFCFV1RingMuLike) continue;    // ANTINU
        // if(s == kSKsRHCFCFV1RingELike) continue;     // ANTINU

        if(s == kSKsFCFV1RingELike || s == kSKsFCFV1RingMuLike || s == kSKsRHCFCFV1RingMuLike) continue; //Skip everything except RHC1Re and Multiring

	
        double nsk = 0;
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            TH2D h = this -> Get2DHist(s,m,t);
            TH2D * mctempl = &h;
            
            //first write the 2d template, and cout the mode integral
            mctempl->Write();
            nsk += mctempl->Integral();
            
            //fill a 1D ereco spectrum
            TH1D * hrec = new TH1D(Form("erec_int_%s",AsStringCompact(m).c_str()), Form("erec_int_%s",AsStringCompact(m).c_str()), erecobins, erecobinedges);
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
            TH1D * hrecf = new TH1D(Form("erec_int_%s_fine",AsStringCompact(m).c_str()), Form("erec_int_%s",AsStringCompact(m).c_str()), mctempl->GetNbinsX(), mctempl->GetXaxis()->GetXbins()->GetArray());
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
            TH1D * htru = new TH1D(Form("enu_int_%s",AsStringCompact(m).c_str()), Form("enu_int_%s",AsStringCompact(m).c_str()), util::NumVariableMCEtrueBins(s,m), util::VariableMCEtrueBinEdges(s,m));
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
        oute->Close();
        delete oute;
    }//sample
    
}
//________________________________________________________________________________________
void SuperkMC::NevDebug2D(SKMCTemplateType_t t, string name)
{
    
    for(int is = 0; is < kns; is++) {
        
        SuperkSample_t s = kSKSamples[is];
        
        if(s == kSKsUndefined) continue;

	// if(s == kSKsFCFV1RingELike || s == kSKsFCFV1RingMuLike || s == kSKsRHCFCFV1RingMuLike) {
	//   LOG("SuperkMC", pDEBUG) << "Skipping templates other than RHC1Re and Multiring";
	//   continue; //Skip everything except RHC1Re and Multiring
	// }
        TFile* oute = new TFile(Form("Template2D_%s%s.root", name.c_str(),util::AsStringCompact(s).c_str()),"RECREATE");
        
        const int q2bins = util::NumQ2Bins(s);
        const double * q2binedges = util::Q2BinEdges(s);
        
        const int erecobins = util::NumErecoBins(s);
        const double * erecobinedges = util::ErecoBinEdges(s);
        
        double nsk = 0;
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            TH3D h = this -> Get3DHist(s,m,t);
            TH3D * mctempl = &h;
            
            
            //this->DrawProjectionXY( mctempl ); // OK!!!
            //util::Die(1); return;);
            
            
            //first write the 3d template, and cout the mode integral
            mctempl->Write();
            LOG("SuperkMC",pINFO) << name << " " << util::AsString(s) << " " << util::AsString(m) << " " << mctempl->Integral();
            nsk += mctempl->Integral();
            
            //fill a 1D q2 spectrum
            TH1D * hq2 = new TH1D(Form("q2_int_2d_%s",AsStringCompact(m).c_str()), Form("q2_int_2d_%s",AsStringCompact(m).c_str()), q2bins, q2binedges);
            // loop over Q2
            for(int iq = 1; iq <= hq2->GetNbinsX(); iq++) {
                //find low and high of q2 for current bin
                double qmin = hq2->GetXaxis()->GetBinLowEdge(iq) + kOneMeV;
                double qmax = hq2->GetXaxis()->GetBinUpEdge(iq)  - kOneMeV;
                //integrate the template
                int iq_min = mctempl->GetXaxis()->FindBin(qmin);
                int iq_max = mctempl->GetXaxis()->FindBin(qmax);
                int ir_bins = mctempl->GetNbinsY();
                int it_bins = mctempl->GetNbinsZ();
                double nev = mctempl->Integral(iq_min, iq_max, 1, ir_bins, 1, it_bins);
                hq2->SetBinContent(iq,nev);
            }
            hq2->Write();
            delete hq2;
            
            //fill a 1D ereco spectrum (fine binning)
            TH1D * hq2f = new TH1D(Form("q2_int_2d_%s_fine",AsStringCompact(m).c_str()), Form("q2_int_2d_%s",AsStringCompact(m).c_str()), mctempl->GetNbinsX(), mctempl->GetXaxis()->GetXbins()->GetArray());
            // loop over Q2
            for(int iq = 1; iq <= hq2f->GetNbinsX(); iq++) {
                //find low and high of q2 for current bin
                double qmin = hq2f->GetXaxis()->GetBinLowEdge(iq) + kOneMeV;
                double qmax = hq2f->GetXaxis()->GetBinUpEdge(iq)  - kOneMeV;
                //integrate the template
                int iq_min = mctempl->GetXaxis()->FindBin(qmin);
                int iq_max = mctempl->GetXaxis()->FindBin(qmax);
                int ir_bins = mctempl->GetNbinsY();
                int it_bins = mctempl->GetNbinsZ();
                double nev = mctempl->Integral(iq_min, iq_max, 1, ir_bins, 1, it_bins);
                hq2f->SetBinContent(iq,nev);
            }
            hq2f->Write();
            delete hq2f;
            
            //fill a 1D ereco spectrum
            TH1D * hrec = new TH1D(Form("erec_int_2d_%s",AsStringCompact(m).c_str()), Form("erec_int_2d_%s",AsStringCompact(m).c_str()), erecobins, erecobinedges);
            // loop over Ereco
            for(int ir = 1; ir <= hrec->GetNbinsX(); ir++) {
                //find low and high of ereco for current bin
                double rmin = hrec->GetXaxis()->GetBinLowEdge(ir) + kOneMeV;
                double rmax = hrec->GetXaxis()->GetBinUpEdge(ir)  - kOneMeV;
                //integrate the template
                int iq_bins = mctempl->GetNbinsX();
                int ir_min = mctempl->GetYaxis()->FindBin(rmin);
                int ir_max = mctempl->GetYaxis()->FindBin(rmax);
                int it_bins = mctempl->GetNbinsZ();
                double nev = mctempl->Integral(1, iq_bins, ir_min, ir_max, 1, it_bins);
                hrec->SetBinContent(ir,nev);
            }
            hrec->Write();
            delete hrec;
            
            //fill a 1D ereco spectrum (fine binning)
            TH1D * hrecf = new TH1D(Form("erec_int_2d_%s_fine",AsStringCompact(m).c_str()), Form("erec_int_2d_%s",AsStringCompact(m).c_str()), mctempl->GetNbinsY(), mctempl->GetYaxis()->GetXbins()->GetArray());
            // loop over Ereco
            for(int ir = 1; ir <= hrecf->GetNbinsX(); ir++) {
                //find low and high of ereco for current bin
                double rmin = hrecf->GetXaxis()->GetBinLowEdge(ir) + kOneMeV;
                double rmax = hrecf->GetXaxis()->GetBinUpEdge(ir)  - kOneMeV;
                //integrate the template
                int iq_bins = mctempl->GetNbinsX();
                int ir_min = mctempl->GetYaxis()->FindBin(rmin);
                int ir_max = mctempl->GetYaxis()->FindBin(rmax);
                int it_bins = mctempl->GetNbinsZ();
                double nev = mctempl->Integral(1, iq_bins, ir_min, ir_max, 1, it_bins);
                hrecf->SetBinContent(ir,nev);
            }
            hrecf->Write();
            delete hrecf;
            
            //fill a 1D etrue spectrum
            TH1D * htru = new TH1D(Form("enu_int_2d_%s",AsStringCompact(m).c_str()), Form("enu_int_2d_%s",AsStringCompact(m).c_str()), util::NumVariableMCEtrueBins(s,m), util::VariableMCEtrueBinEdges(s,m));
            // loop over Etrue
            for(int it = 1; it <= htru->GetNbinsX(); it++) {
                //find low and high of etrue for current bin
                double tmin = htru->GetXaxis()->GetBinLowEdge(it) + kOneMeV;
                double tmax = htru->GetXaxis()->GetBinUpEdge(it)  - kOneMeV;
                //integrate the template
                int it_min = mctempl->GetZaxis()->FindBin(tmin);
                int it_max = mctempl->GetZaxis()->FindBin(tmax);
                int iq_bins = mctempl->GetNbinsX();
                int ir_bins = mctempl->GetNbinsY();
                double nev = mctempl->Integral(1, iq_bins, 1, ir_bins, it_min, it_max);
                htru->SetBinContent(it,nev);
            }
            htru->Write();
            delete htru;
            
        }//mode
        
        
        //this->UpdateErecoSpectra();
        this->UpdateErecoVsQ2Spectra();
        
        //get the total spectrum
        
        //TH1D * hfull = this->CurrErecoSpectrum(s);
        //hfull->Write("ereco_spectrum");
        //delete hfull;
        
        TH2D * hfull2D = this->CurrErecoVsQ2Spectrum(s);
        hfull2D->Write("erecoVsq2_spectrum");
        delete hfull2D;
        
        LOG("SuperkMC",pINFO) << "Total: " << nsk;
        oute->Close();
        delete oute;
    }//sample
}

//








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
        TH1D * spectrum = new TH1D("spe","spe",n_ereco_bins, util::ErecoBinEdges(kSKsFCFV1RingMuLike));
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
                            LOG("SuperkMC",pINFO) << mctempl->GetBinContent(ix,iy) << "\t" << mctemplerr->GetBinContent(ix,iy);
                    double nev_old = spectrum_with2derr->GetBinContent(ix);
                    spectrum_with2derr->SetBinContent(ix, nev_new + nev_old);
                    for(int ie = 1; ie < ntoys; ie++) {
                        int nev_fluct = Rnd::Instance()->RndToyMCstat().Poisson((int)nev_new);
                        double current = spectrum_rms->GetBinContent(ix);
                        spectrum_rms->SetBinContent(ix, current + TMath::Power(nev_fluct-nev_new,2));
                    }//ie
                }//iy
            }//ix
            LOG("SuperkMC",pDEBUG) << "Finished MC stat error debug loops for mode: " << util::AsStringCompact(m);
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
                    h.SetBinContent(binno, fMCTemplArrayOri[is]->GetEntryRef(im,it,ir));
                }
                if (type == conv::kSKMCtNom) {
                    h.SetBinContent(binno, fMCTemplArrayNom[is]->GetEntryRef(im,it,ir));
                }
                if (type == conv::kSKMCtTun) {
                    h.SetBinContent(binno, fMCTemplArrayTun[is]->GetEntryRef(im,it,ir));
                }
                if (type == conv::kSKMCtFin) {
                    h.SetBinContent(binno, fMCTemplArrayFin[is]->GetEntryRef(im,it,ir));
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
                h.SetBinContent(binno, fMCTemplArrayOri[is]->GetEntryRef(im,it,ir));
            }
            if (type == conv::kSKMCtNom) {
                h.SetBinContent(binno, fMCTemplArrayNom[is]->GetEntryRef(im,it,ir));
            }
            if (type == conv::kSKMCtTun) {
                h.SetBinContent(binno, fMCTemplArrayTun[is]->GetEntryRef(im,it,ir));
            }
            if (type == conv::kSKMCtFin) {
                h.SetBinContent(binno, fMCTemplArrayFin[is]->GetEntryRef(im,it,ir));
            }
        }
    }
    
    return h;
}
//____________________________________________________________________________________








// NEW 2D TEMPLATES

TH3D SuperkMC::Get3DHist(conv::SuperkSample_t s, conv::SKMCTemplateType_t type) const
{
    int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
    
    const int      n_q2_bins       = kNQ2Bins[is];
    const int      n_ereco_bins    = kNErecoBins[is];
    const int      n_etrue_bins    = kNEtrueBins[is];
    const double * q2_bin_edges    = util::Q2BinEdges(s);
    const double * ereco_bin_edges = util::ErecoBinEdges(s);
    
    TH3D h(Form("nsk_2d_%d",is), util::AsString(s).c_str(),
           n_q2_bins, q2_bin_edges, n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
    
    h.SetDirectory(0);
    h.GetXaxis()->SetTitle("Q^{2} (GeV)");
    h.GetYaxis()->SetTitle("E_{reco} (GeV)");
    h.GetZaxis()->SetTitle("E_{true} (GeV)");
    
    for (int im = 0 ; im < kNSuperkModes ; im++ ) {
        SuperkMode_t m = kSKReactionModes[is][im];
        if (m==kSKmUndefined) continue;
        for (int it = 0 ; it < n_etrue_bins ; it++) {
            for (int ir = 0 ; ir < n_ereco_bins ; ir++) {
                for (int iq = 0 ; iq < n_q2_bins ; iq++) {
                    int binno = h.GetBin(iq+1,ir+1,it+1);
                    if (type == conv::kSKMCtOri) {
                        h.SetBinContent(binno, fMCTemplArrayOri2D[is]->GetEntryRef(im,it,ir,iq));
                    }
                    if (type == conv::kSKMCtNom) {
                        h.SetBinContent(binno, fMCTemplArrayNom2D[is]->GetEntryRef(im,it,ir,iq));
                    }
                    if (type == conv::kSKMCtTun) {
                        h.SetBinContent(binno, fMCTemplArrayTun2D[is]->GetEntryRef(im,it,ir,iq));
                    }
                    if (type == conv::kSKMCtFin) {
                        h.SetBinContent(binno, fMCTemplArrayFin2D[is]->GetEntryRef(im,it,ir,iq));
                    }
                }
            }
        }
    }
    return h;
}
//____________________________________________________________________________________
TH3D SuperkMC::Get3DHist(conv::SuperkSample_t s, conv::SuperkMode_t m , conv::SKMCTemplateType_t type) const
{
    int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
    int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
    const int      n_q2_bins       = kNQ2Bins[is];
    const int      n_ereco_bins    = kNErecoBins[is];
    const int      n_etrue_bins    = kNEtrueBins[is];
    const double * q2_bin_edges    = util::Q2BinEdges(s);
    const double * ereco_bin_edges = util::ErecoBinEdges(s);
    
    TH3D h(Form("nsk_2d_%d_%d",is,im), util::AsString(m).c_str(),
           n_q2_bins, q2_bin_edges, n_ereco_bins, ereco_bin_edges, n_etrue_bins, util::EtrueBinEdges(s));
    
    h.SetDirectory(0);
    h.GetXaxis()->SetTitle("Q^{2} (GeV)");
    h.GetYaxis()->SetTitle("E_{reco} (GeV)");
    h.GetZaxis()->SetTitle("E_{true} (GeV)");
    
    for (int it = 0 ; it < n_etrue_bins ; it++) {
        for (int ir = 0 ; ir < n_ereco_bins ; ir++) {
            for (int iq = 0 ; iq < n_q2_bins ; iq++) {
                int binno = h.GetBin(iq+1,ir+1,it+1);
                if (type == conv::kSKMCtOri) {
                    h.SetBinContent(binno, fMCTemplArrayOri2D[is]->GetEntryRef(im,it,ir,iq));
                }
                if (type == conv::kSKMCtNom) {
                    h.SetBinContent(binno, fMCTemplArrayNom2D[is]->GetEntryRef(im,it,ir,iq));
                }
                if (type == conv::kSKMCtTun) {
                    h.SetBinContent(binno, fMCTemplArrayTun2D[is]->GetEntryRef(im,it,ir,iq));
                }
                if (type == conv::kSKMCtFin) {
                    h.SetBinContent(binno, fMCTemplArrayFin2D[is]->GetEntryRef(im,it,ir,iq));
                }
            }
        }
    }
    return h;
}
//



//____________________________________________________________________________________
bool SuperkMC::DoTemplateArraysHaveNIWGTune(void) const
{
    return fTemplateArraysHaveNIWGTune;
}
//____________________________________________________________________________________







// NEW 2D TEMPLATES
void SuperkMC::DrawProjectionXY(TH3D *h){ // for debugging
    
    //Template2D(kSKsFCFV1RingELike,kSKmOscNuECCQE,kSKMCtNom);
    
    //this->DrawProjectionXY( mctempl ); // OK!!!
    //util::Die(1); return;;
    
    TCanvas *c = new TCanvas("c","c");
    h->DrawCopy("colz");
    
    TCanvas *c1 = new TCanvas("c1","c1");
    h->Project3D("xy")->DrawCopy("colz");
    
    c->Print("ciao.pdf");
    c1->Print("proj_xy.pdf");
    
    return;
}
//____________________________________________________________________________________

void SuperkMC::ProjectTemplate3Dto2D(TH3D *h3d[conv::kns][conv::knm], TH2D *h2d[conv::kns][conv::knm]){
    
    // Project from (Q2,Ereco,Etrue) to (Ereco,Etrue)
    
    for(int is=0;is<conv::kns;is++){
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        const int      n_q2_bins       = kNQ2Bins[is];
        const int      n_ereco_bins    = kNErecoBins[is];
        const int      n_etrue_bins    = kNEtrueBins[is];
        
        for(int im=0;im<conv::knm;im++){
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            
            for (int it = 0 ; it < n_etrue_bins ; it++) {
                for (int ir = 0 ; ir < n_ereco_bins ; ir++) {
                    
                    double cont = 0.;
                    for (int iq = 0 ; iq < n_q2_bins ; iq++) {
                        cont += h3d[is][im]->GetBinContent(iq+1,ir+1,it+1);
                    }
                    
                    h2d[is][im]->SetBinContent(ir+1,it+1,cont);
                }
            }
            
        }
    }
    
    return;
}
//____________________________________________________________________________________
double SuperkMC::GetTemplateIntegral(TH3D *h3d[conv::kns][conv::knm]){
    
    double integral = 0.;
    
    for(int is=0;is<conv::kns;is++){
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        const int      n_q2_bins       = kNQ2Bins[is];
        const int      n_ereco_bins    = kNErecoBins[is];
        const int      n_etrue_bins    = kNEtrueBins[is];
        
        for(int im=0;im<conv::knm;im++){
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            
            for (int it = 0 ; it < n_etrue_bins ; it++) {
                for (int ir = 0 ; ir < n_ereco_bins ; ir++) {
                    for (int iq = 0 ; iq < n_q2_bins ; iq++) {
                        integral += h3d[is][im]->GetBinContent(iq+1,ir+1,it+1);
                    }
                }
            }
            
        }
    }
    
    return integral;
}
//____________________________________________________________________________________
double SuperkMC::GetTemplateIntegral(TH2D *h2d[conv::kns][conv::knm]){
    
    double integral = 0.;
    
    for(int is=0;is<conv::kns;is++){
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
        //const int      n_q2_bins       = kNQ2Bins[is];
        const int      n_ereco_bins    = kNErecoBins[is];
        const int      n_etrue_bins    = kNEtrueBins[is];
        
        for(int im=0;im<conv::knm;im++){
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            
            for (int it = 0 ; it < n_etrue_bins ; it++) {
                for (int ir = 0 ; ir < n_ereco_bins ; ir++) {
                    integral += h2d[is][im]->GetBinContent(ir+1,it+1);
                }
            }
            
        }
    }
    
    return integral;
}
//____________________________________________________________________________________
double SuperkMC::GetTemplateIntegral(TH3D *h3d){
    
    double integral = 0.;
    
    for (int it = 0 ; it < h3d->GetZaxis()->GetNbins() ; it++) {
        for (int ir = 0 ; ir < h3d->GetYaxis()->GetNbins() ; ir++) {
            for (int iq = 0 ; iq < h3d->GetXaxis()->GetNbins() ; iq++) {
                integral += h3d->GetBinContent(iq+1,ir+1,it+1);
            }
        }
    }
    
    return integral;
}
//____________________________________________________________________________________
double SuperkMC::GetTemplateIntegral(TH2D *h2d){
    
    double integral = 0.;
    
    for (int it = 0 ; it < h2d->GetZaxis()->GetNbins() ; it++) {
        for (int ir = 0 ; ir < h2d->GetYaxis()->GetNbins() ; ir++) {
            integral += h2d->GetBinContent(ir+1,it+1);
        }
    }
    
    return integral;
}

//____________________________________________________________________________________
double util::CalcThetaDegrees(float dir[][7][3], conv::SuperkSample_t s )
{
  //double dir_nu[3] = {0.669764, -0.742179, 0.024223}; // T2K beam direction
  //double costheta = dir_nu[0] * dir[0][0] + dir_nu[1] * dir[0][1] + dir_nu[2] * dir[0][2];
  // LOG("util", pERROR) << "Filling theta for sample: " << util::SampleAsString(s)
  //                         << " Is1RESample?: "  << util::Is1ReSample(s)
  //                         << " Is1RMuSample?: " << util::Is1RMuSample(s);
  
  double costheta = -999;
  if (util::Is1ReSample(s) || util::IsMultiReSample(s)) {
    costheta = kBeamDir[0] * dir[0][1][0] + kBeamDir[1] * dir[0][1][1] + kBeamDir[2] * dir[0][1][2];
  }
  if (util::Is1RMuSample(s)) {
    costheta = kBeamDir[0] * dir[0][2][0] + kBeamDir[1] * dir[0][2][1] + kBeamDir[2] * dir[0][2][2];    
  }
  if( costheta<-1.00001 || costheta>+1.00001 ){
    LOG("util",pERROR) << "util::FillTemplateHistogram" ;
    LOG("util",pERROR) << "costheta = " << costheta << " --> outside allowed range" ;
    util::Die(1); return 0;
    exit(1);
  }    
  double LepTheta        = TMath::ACos(costheta); // radiant    
  double LepThetaDegrees = LepTheta / TMath::Pi() * 180.; // degrees   
  if(LepTheta<0.){
    LOG("util",pERROR) << "Theta lepton is negative: " << LepThetaDegrees << " degrees" ;
    exit(1);
  }
  return LepThetaDegrees;
}
//____________________________________________________________________________________
double util::CalcCosTheta(float dir[][7][3], conv::SuperkSample_t s )
{
  //double dir_nu[3] = {0.669764, -0.742179, 0.024223}; // T2K beam direction
  //double costheta = dir_nu[0] * dir[0][0] + dir_nu[1] * dir[0][1] + dir_nu[2] * dir[0][2];
  // LOG("util", pERROR) << "Filling theta for sample: " << util::SampleAsString(s)
  //                         << " Is1RESample?: "  << util::Is1ReSample(s)
  //                         << " Is1RMuSample?: " << util::Is1RMuSample(s);
  
  double costheta = -999;
  if (util::Is1ReSample(s) || util::IsMultiReSample(s)) {
    costheta = kBeamDir[0] * dir[0][1][0] + kBeamDir[1] * dir[0][1][1] + kBeamDir[2] * dir[0][1][2];
  }
  if (util::Is1RMuSample(s)) {
    costheta = kBeamDir[0] * dir[0][2][0] + kBeamDir[1] * dir[0][2][1] + kBeamDir[2] * dir[0][2][2];    
  }
  if( costheta<-1.00001 || costheta>+1.00001 ){
    LOG("util",pERROR) << "util::FillTemplateHistogram" ;
    LOG("util",pERROR) << "costheta = " << costheta << " --> outside allowed range" ;
    util::Die(1); return 0;
    exit(1);
  }    
  return costheta;
}
//____________________________________________________________________________________
double util::CalcQ2(double Enu,double Elep,double Plep,double costheta,double Mlep)
{
    // Enu      --> neutrino reco energy
    // Elep     --> lepton reco energy
    // Plep     --> lepton reco momentum
    // costheta --> angle between neutrino and lepton
    
    // Q2 in the laboratory reference system (Giunti section 5.3.1)
    double Q2 = 2*Enu*(Elep-Plep*costheta) - Mlep*Mlep;
    
    return Q2;
}
//____________________________________________________________________________________
double util::CalcQd2(double Enu,double Elep,double Plep,double costheta,double Mlep)
{
    // Qd^2 (variation of Q^2)
    
    // Enu      --> neutrino reco energy
    // Elep     --> lepton reco energy
    // Plep     --> lepton reco momentum
    // costheta --> angle between neutrino and lepton
    
    double Q2 = CalcQ2(Enu,Elep,Plep,costheta,Mlep);
    
    double Qd2 = TMath::Power(1+Q2/0.7, -1);
    
    return Qd2;
}
//____________________________________________________________________________________
double util::CalcEhad(double Enu,double Elep)
{
    // Energy of the hadronic recoil system
    
    // Enu      --> neutrino reco energy
    // Elep     --> lepton reco energy
    
    double nu = Enu - Elep;
    
    return nu;
}
//____________________________________________________________________________________
double util::CalcInelasticity(double Enu,double Elep)
{
    // Inelasticity = (Enu - Elep) / Enu
    
    // Enu      --> neutrino reco energy
    // Elep     --> lepton reco energy
    
    double Inel = (Enu - Elep) / Enu;
    
    return Inel;
}
//____________________________________________________________________________________
double util::CalcBjorkenScale(double Enu,double Elep,double Plep,double costheta,double Mlep,double Mn)
{
    // x = Q2 / (2*Mn*Enu*Inelasticity)
    
    // Enu      --> neutrino reco energy
    // Elep     --> lepton reco energy
    // Plep     --> lepton reco momentum
    // costheta --> angle between neutrino and lepton
    // Mn       --> nucleon mass
    
    double Q2 = CalcQ2(Enu,Elep,Plep,costheta,Mlep);
    double Inel = CalcInelasticity(Enu,Elep);
    double x = Q2 / (2*Mn*Enu*Inel);
    
    return x;
}
//____________________________________________________________________________________
double util::CalcWhadr(double Enu,double Elep,double Plep,double costheta,double Mlep,double Mn)
{
  // x = Q2 / (2*Mn*Enu*Inelasticity)
    
  // Enu      --> neutrino reco energy
  // Elep     --> lepton reco energy
  // Plep     --> lepton reco momentum
  // costheta --> angle between neutrino and lepton
  // Mn       --> nucleon mass
    
  double Q2 = CalcQ2(Enu,Elep,Plep,costheta,Mlep);
  double nu = CalcEhad(Enu,Elep);
  
  double W = Mn*Mn + 2*Mn*nu - Q2;
  
  return W;
}
//____________________________________________________________________________________
double util::CalcXSecTerm(double Enu,double Elep,double Plep,double costheta,double Mlep,double Mn){
  // Term of CCQE dsigma/dQ^2 which sign is different for nu/nubar (arXiv:1305.7513v1, eq. 57-61)

  // Enu      --> neutrino reco energy 
  // Elep     --> lepton reco energy 
  // Plep     --> lepton reco momentum
  // costheta --> angle between neutrino and lepton 
  // Mn       --> nucleon mass

  double MA = 1.2257; // Axial Mass, BANFF post-fit (GeV) 
  double gA = 1.2694; // FA(Q2=0) 
  
  double Q2 = CalcQ2(Enu,Elep,Plep,costheta,Mlep);

  double F1 = 0.5;
  double F2 = 0.5;
  double FA = gA / TMath::Power( 1 + Q2/(MA*MA), 2. );
  double B = Q2/(Mn*Mn)*FA*(F1+F2);

  double XSecTerm = (4*Mn*Enu - Q2 - Mlep*Mlep) / (Mn*Mn) * B;
 
  return XSecTerm;
}
//____________________________________________________________________________________
double util::ErecoFunc_Nu(double Mom,double costh,double mass_lep)
{
    const double Eb = 27. / 1000.;     // oxygen
    const double mp = 938.272 / 1000.; // proton mass (MeV)
    const double mn = 939.565 / 1000.; // neutron mass (MeV)
    
    double Ereco=0.;
    double E_lep=0.;
    
    E_lep = sqrt( Mom*Mom + mass_lep*mass_lep );    // Reco energy muon or electron
    Ereco = mp*mp - (mn - Eb)*(mn - Eb) - mass_lep*mass_lep + 2*(mn - Eb)*E_lep;   // right formula
    Ereco /= 2*(mn - Eb - E_lep + Mom*costh);
    
    return Ereco;  // GeV
}
//____________________________________________________________________________________
double util::ErecoFunc_NuBar(double Mom,double costh,double mass_lep)
{    
    const double Eb = 27. / 1000.;     // oxygen
    const double mp = 938.272 / 1000.; // proton mass (MeV)
    const double mn = 939.565 / 1000.; // neutron mass (MeV)
    
    double Ereco=0.;
    double E_lep=0.;
    
    E_lep = sqrt( Mom*Mom + mass_lep*mass_lep );    // Reco energy lepton
    Ereco = mn*mn - (mp - Eb)*(mp - Eb) - mass_lep*mass_lep + 2*(mp - Eb)*E_lep;
    Ereco /= 2*(mp - Eb - E_lep + Mom*costh);
    
    return Ereco;  // GeV
}
//____________________________________________________________________________________
double util::ErecoFuncCCOnePi_Nu(double Mom,double costh,double mass_lep)
{
  const double m_start = 938.272 / 1000.; //initial nucleon mass
  const double mn_end = 1232.0 / 1000.;   //final nucleon mass
  
  double Ereco=0.;
  double E_lep=0.;
  
  E_lep = sqrt( Mom*Mom + mass_lep*mass_lep );    // Reco energy muon or electron
  Ereco = mn_end*mn_end - (m_start)*(m_start) - mass_lep*mass_lep + 2*(m_start)*E_lep;   // right formula
  Ereco /= 2*(m_start - E_lep + Mom*costh);
  
  return Ereco;  // GeV
}
//____________________________________________________________________________________
double util::SmearMomentum(double fit_momentum,Option_t * option ) 
{
  //Taken from T2K paper: http://arxiv.org/pdf/1502.01550v2.pdf (bottom of p53)
    TRandom3 & rndgen = Rnd::Instance()->RndOther();
    
    double error = 1;
    TString opt = option;
    if(opt.Contains("electron")){
        error = fit_momentum * (0.006 + 0.026 / TMath::Sqrt(fit_momentum) );
    }
    else if (opt.Contains("muon"))
    {
        error = fit_momentum * (0.017 + 0.007 / TMath::Sqrt(fit_momentum));
    }
    
    return rndgen.Gaus(fit_momentum,error);
}
//____________________________________________________________________________________




