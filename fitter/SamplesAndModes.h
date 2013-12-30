//____________________________________________________________________________
/*!

\class    t2kana::conv

\brief    Conventions

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

Changes since 2011_v0.9.0:

17-12-2013  Nick  Split t2kana::conv namespace into several smaller files

*/
//____________________________________________________________________________

#ifndef samples_modes_h
#define samples_modes_h

#include "Enumerations.h"
#include "AnalysisOptions.h"

namespace t2kana {
namespace conv   {

//
// ------------------------------------------------------------------------
// ** NOTE: IMPORTANT CHOICES MADE HERE **
//
// Choice of samples included in the analysis 
// Choice of modes (per sample) we are going to keep track of.
//
// Note that some of the arrays defined further below (such as formatting 
// choices or binning choices for systematic parameter weight splines) 
// depend on kSKReactionModes. Make sure changes are made consistently.
// ------------------------------------------------------------------------
//

//
// Define which SK samples are considered
//

static const SuperkSample_t kSKSamples[kns] = {

#if defined __FHC_RUNNING__ && defined __RHC_RUNNING__
#ifdef __NUMUONLY__
  kSKsFCFV1RingMuLike,
  kSKsUndefined,
  kSKsRHCFCFV1RingMuLike,
  kSKsUndefined
#else
#ifdef __NUEONLY__
  kSKsUndefined,
  kSKsFCFV1RingELike,
  kSKsUndefined,
  kSKsRHCFCFV1RingELike
#else
  kSKsFCFV1RingMuLike,
  kSKsFCFV1RingELike,
  kSKsRHCFCFV1RingMuLike,
  kSKsRHCFCFV1RingELike
#endif
#endif

#else
#ifdef __FHC_RUNNING__
#ifdef __NUMUONLY__
  kSKsFCFV1RingMuLike,
  kSKsUndefined,
  kSKsUndefined,
  kSKsUndefined
#else
#ifdef __NUEONLY__
  kSKsUndefined,
  kSKsFCFV1RingELike,
  kSKsUndefined,
  kSKsUndefined
#else
  kSKsFCFV1RingMuLike,
  kSKsFCFV1RingELike,
  kSKsUndefined,
  kSKsUndefined
#endif
#endif

#else
#ifdef __RHC_RUNNING__
#ifdef __NUMUONLY__
  kSKsUndefined,
  kSKsUndefined,
  kSKsRHCFCFV1RingMuLike,
  kSKsUndefined
#else
#ifdef __NUEONLY__
  kSKsUndefined,
  kSKsUndefined,
  kSKsUndefined,
  kSKsRHCFCFV1RingELike
#else
  kSKsUndefined,
  kSKsUndefined,
  kSKsRHCFCFV1RingMuLike,
  kSKsRHCFCFV1RingELike
#endif
#endif
#endif
#endif
#endif
};

//
// Define which SK modes are considered for each SK sample included in the analysis
//

static const SuperkMode_t kSKReactionModes[kns][knm] = {
#if defined __FHC_RUNNING__ && defined __RHC_RUNNING__
  // > 1 muon-like ring (numu disappearance)
  {
    kSKmNuMuCCQE,         
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    //    kSKmNuMuCCoth_mu,                      
    kSKmNuMuCCoth,                      
    kSKmNuMuNC1piC,
    kSKmNuMuNCoth_mu,          
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmOscNuMuCCQE,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCcoh,
    //    kSKmOscNuMuCCoth_mu,
    kSKmOscNuMuCCoth,
#endif
    kSKmNuMuBarCCQE,      
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCcoh,
    //    kSKmNuMuBarCCoth_mu,                  
    kSKmNuMuBarCCoth,                  
    kSKmNuMuBarNC1piC,
    kSKmNuMuBarNCoth_mu,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCQE,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCcoh,
    //    kSKmOscNuMuCCoth_mu,
    kSKmOscNuMuBarCCoth,
#endif
    kSKmNuECCQE, 
    kSKmNuECC1pi,
    kSKmNuECCcoh,
    //    kSKmNuECCoth_mu,                                
    kSKmNuECCoth,                                
    kSKmNuENC1piC,
    kSKmNuENCoth_mu,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    //    kSKmNuEBarCCoth_mu,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1piC,
    kSKmNuEBarNCoth_mu,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmOscNuEBarCCQE,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCcoh,
    //    kSKmOscNuMuCCoth_mu,
    kSKmOscNuEBarCCoth,
#endif
    kSKmSigNuECCQE,
    kSKmSigNuECC1pi,
    kSKmSigNuECCcoh,
    //    kSKmSigNuECCoth_mu,
    kSKmSigNuECCoth,
#ifdef USE_MEC
    kSKmMEC,
#else
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
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
  }

  ,

  // > 1 electron-like ring (nue appearance)
  {
    kSKmNuMuCCQE,
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    //    kSKmNuMuCCoth_e, 
    kSKmNuMuCCoth, 
    kSKmNuMuNC1pi0, 
    kSKmNuMuNCcoh, 
    kSKmNuMuNCoth_e, 
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__    
    kSKmOscNuMuCCQE,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCcoh,
    //    kSKmOscNuMuCCoth_e,                            
    kSKmOscNuMuCCoth,               
#endif             
    kSKmNuMuBarCCQE, 
    kSKmNuMuBarCC1pi, 
    kSKmNuMuBarCCcoh, 
    //    kSKmNuMuBarCCoth_e,  
    kSKmNuMuBarCCoth,  
    kSKmNuMuBarNC1pi0, 
    kSKmNuMuBarNCcoh, 
    kSKmNuMuBarNCoth_e,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCQE,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCcoh,
    //    kSKmOscNuMuCCoth_e,                            
    kSKmOscNuMuBarCCoth,               
#endif             
    kSKmNuECCQE, 
    kSKmNuECC1pi, 
    kSKmNuECCcoh, 
    //    kSKmNuECCoth_e,
    kSKmNuECCoth,
    kSKmNuENC1pi0,  
    kSKmNuENCcoh, 
    kSKmNuENCoth_e,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    //    kSKmNuEBarCCoth_e,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1pi0,
    kSKmNuEBarNCcoh,
    kSKmNuEBarNCoth_e,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmOscNuEBarCCQE,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCcoh,
    //    kSKmOscNuMuCCoth_e,                            
    kSKmOscNuEBarCCoth,               
#endif          
    kSKmSigNuECCQE,  
    kSKmSigNuECC1pi,  
    kSKmSigNuECCcoh,
    //    kSKmSigNuECCoth_e,
    kSKmSigNuECCoth,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined
  }

  ,

  // > 1 muon-like ring (numubar disappearance) with reverse horn current
  {
    kSKmNuMuCCQE,         
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    //    kSKmNuMuCCoth_mu,                      
    kSKmNuMuCCoth,                      
    kSKmNuMuNC1piC,
    kSKmNuMuNCoth_mu,          
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmOscNuMuCCQE,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCcoh,
    //    kSKmOscNuMuCCoth_mu,
    kSKmOscNuMuCCoth,
#endif
    kSKmNuMuBarCCQE,      
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCcoh,
    //    kSKmNuMuBarCCoth_mu,                  
    kSKmNuMuBarCCoth,                  
    kSKmNuMuBarNC1piC,
    kSKmNuMuBarNCoth_mu,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCQE,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCcoh,
    //    kSKmOscNuMuCCoth_mu,
    kSKmOscNuMuBarCCoth,
#endif
    kSKmNuECCQE, 
    kSKmNuECC1pi,
    kSKmNuECCcoh,
    //    kSKmNuECCoth_mu,                                
    kSKmNuECCoth,                                
    kSKmNuENC1piC,
    kSKmNuENCoth_mu,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    //    kSKmNuEBarCCoth_mu,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1piC,
    kSKmNuEBarNCoth_mu,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmOscNuEBarCCQE,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCcoh,
    //    kSKmOscNuMuCCoth_mu,
    kSKmOscNuEBarCCoth,
#endif
    kSKmSigNuECCQE,
    kSKmSigNuECC1pi,
    kSKmSigNuECCcoh,
    //    kSKmSigNuECCoth_mu,
    kSKmSigNuECCoth,
#ifdef USE_MEC
    kSKmMEC,
#else
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
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
  }

  ,

  // > 1 electron-like ring (nuebar appearance) with reverse horn current
  {
    kSKmNuMuCCQE,
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    //    kSKmNuMuCCoth_e, 
    kSKmNuMuCCoth,
    kSKmNuMuNC1pi0, 
    kSKmNuMuNCcoh, 
    kSKmNuMuNCoth_e, 
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__    
    kSKmOscNuMuCCQE,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCcoh,
    //    kSKmOscNuMuCCoth_e,                            
    kSKmOscNuMuCCoth,               
#endif             
    kSKmNuMuBarCCQE, 
    kSKmNuMuBarCC1pi, 
    kSKmNuMuBarCCcoh, 
    //    kSKmNuMuBarCCoth_e,  
    kSKmNuMuBarCCoth,  
    kSKmNuMuBarNC1pi0, 
    kSKmNuMuBarNCcoh, 
    kSKmNuMuBarNCoth_e,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCQE,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCcoh,
    //    kSKmOscNuMuCCoth_e,                            
    kSKmOscNuMuBarCCoth,               
#endif             
    kSKmNuECCQE, 
    kSKmNuECC1pi, 
    kSKmNuECCcoh, 
    //    kSKmNuECCoth_e,
    kSKmNuECCoth,
    kSKmNuENC1pi0,  
    kSKmNuENCcoh, 
    kSKmNuENCoth_e,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    //    kSKmNuEBarCCoth_e,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1pi0,
    kSKmNuEBarNCcoh,
    kSKmNuEBarNCoth_e,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmOscNuEBarCCQE,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCcoh,
    //    kSKmOscNuMuCCoth_e,                            
    kSKmOscNuEBarCCoth,               
#endif     
    kSKmSigNuECCQE,  
    kSKmSigNuECC1pi,  
    kSKmSigNuECCcoh,
    //    kSKmSigNuECCoth_e,
    kSKmSigNuECCoth,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined
  }

#else

#ifdef __FHC_RUNNING__
  // > 1 muon-like ring (numu disappearance)
  {
    kSKmNuMuCCQE,         
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    //    kSKmNuMuCCoth_mu,                      
    kSKmNuMuCCoth,                      
    kSKmNuMuNC1piC,
    kSKmNuMuNCoth_mu,          
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmOscNuMuCCQE,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCcoh,
    //    kSKmOscNuMuCCoth_mu,
    kSKmOscNuMuCCoth,
#endif
    kSKmNuMuBarCCQE,      
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCcoh,
    //    kSKmNuMuBarCCoth_mu,                  
    kSKmNuMuBarCCoth,                  
    kSKmNuMuBarNC1piC,
    kSKmNuMuBarNCoth_mu,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCQE,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCcoh,
    //    kSKmOscNuMuCCoth_mu,
    kSKmOscNuMuBarCCoth,
#endif
    kSKmNuECCQE, 
    kSKmNuECC1pi,
    kSKmNuECCcoh,
    //    kSKmNuECCoth_mu,                                
    kSKmNuECCoth,                                
    kSKmNuENC1piC,
    kSKmNuENCoth_mu,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    //    kSKmNuEBarCCoth_mu,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1piC,
    kSKmNuEBarNCoth_mu,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmOscNuEBarCCQE,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCcoh,
    //    kSKmOscNuMuCCoth_mu,
    kSKmOscNuEBarCCoth,
#endif
    kSKmSigNuECCQE,
    kSKmSigNuECC1pi,
    kSKmSigNuECCcoh,
    //    kSKmSigNuECCoth_mu,
    kSKmSigNuECCoth,
#ifdef USE_MEC
    kSKmMEC,
#else
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
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
  }

  ,

  // > 1 electron-like ring (nue appearance)
  {
    kSKmNuMuCCQE,
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    //    kSKmNuMuCCoth_e, 
    kSKmNuMuCCoth, 
    kSKmNuMuNC1pi0, 
    kSKmNuMuNCcoh, 
    kSKmNuMuNCoth_e, 
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__    
    kSKmOscNuMuCCQE,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCcoh,
    //    kSKmOscNuMuCCoth_e,                            
    kSKmOscNuMuCCoth,               
#endif             
    kSKmNuMuBarCCQE, 
    kSKmNuMuBarCC1pi, 
    kSKmNuMuBarCCcoh, 
    //    kSKmNuMuBarCCoth_e,  
    kSKmNuMuBarCCoth,  
    kSKmNuMuBarNC1pi0, 
    kSKmNuMuBarNCcoh, 
    kSKmNuMuBarNCoth_e,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCQE,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCcoh,
    //    kSKmOscNuMuCCoth_e,                            
    kSKmOscNuMuBarCCoth,               
#endif             
    kSKmNuECCQE, 
    kSKmNuECC1pi, 
    kSKmNuECCcoh, 
    //    kSKmNuECCoth_e,
    kSKmNuECCoth,
    kSKmNuENC1pi0,  
    kSKmNuENCcoh, 
    kSKmNuENCoth_e,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    //    kSKmNuEBarCCoth_e,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1pi0,
    kSKmNuEBarNCcoh,
    kSKmNuEBarNCoth_e,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmOscNuEBarCCQE,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCcoh,
    //    kSKmOscNuMuCCoth_e,                            
    kSKmOscNuEBarCCoth,               
#endif             
    kSKmSigNuECCQE,  
    kSKmSigNuECC1pi,  
    kSKmSigNuECCcoh,
    //    kSKmSigNuECCoth_e,
    kSKmSigNuECCoth,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined
  }

    ,

  // > 1 muon-like ring (numubar disappearance) with reverse horn current
  {
    kSKmUndefined,         
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,                      
    kSKmUndefined,                      
    kSKmUndefined,
    kSKmUndefined,          
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,
    kSKmUndefined,
#endif
    kSKmUndefined,      
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,                  
    kSKmUndefined,                  
    kSKmUndefined,
    kSKmUndefined,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,
    kSKmUndefined,
#endif
    kSKmUndefined, 
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,                                
    kSKmUndefined,                                
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,
    kSKmUndefined,
#endif
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,
    kSKmUndefined,
#ifdef USE_MEC
    kSKmUndefined,
#else
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
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
  }

  ,

  // > 1 electron-like ring (nuebar appearance) with reverse horn current
  {
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //   kSKmUndefined, 
    kSKmUndefined, 
    kSKmUndefined, 
    kSKmUndefined, 
    kSKmUndefined, 
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__    
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,                            
    kSKmUndefined,               
#endif             
    kSKmUndefined, 
    kSKmUndefined, 
    kSKmUndefined, 
    //   kSKmUndefined,  
    kSKmUndefined,  
    kSKmUndefined, 
    kSKmUndefined, 
    kSKmUndefined,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,                            
    kSKmUndefined,               
#endif             
    kSKmUndefined, 
    kSKmUndefined, 
    kSKmUndefined, 
    //    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,  
    kSKmUndefined, 
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //   kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,                            
    kSKmUndefined,               
#endif     
    kSKmUndefined,  
    kSKmUndefined,  
    kSKmUndefined,
    //    kSKmUndefined,
    kSKmUndefined,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined
  }

#else

#ifdef __RHC_RUNNING__
  // > 1 muon-like ring (numu disappearance)
  {
    kSKmUndefined,         
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,                      
    kSKmUndefined,                      
    kSKmUndefined,
    kSKmUndefined,          
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,
    kSKmUndefined,
#endif
    kSKmUndefined,      
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,                  
    kSKmUndefined,                  
    kSKmUndefined,
    kSKmUndefined,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,
    kSKmUndefined,
#endif
    kSKmUndefined, 
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,                                
    kSKmUndefined,                                
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,
    kSKmUndefined,
#endif
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,
    kSKmUndefined,
#ifdef USE_MEC
    kSKmUndefined,
#else
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
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
  }

  ,

  // > 1 electron-like ring (nue appearance) 
  {
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //   kSKmUndefined, 
    kSKmUndefined, 
    kSKmUndefined, 
    kSKmUndefined, 
    kSKmUndefined, 
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__    
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,                            
    kSKmUndefined,               
#endif             
    kSKmUndefined, 
    kSKmUndefined, 
    kSKmUndefined, 
    //   kSKmUndefined,  
    kSKmUndefined,  
    kSKmUndefined, 
    kSKmUndefined, 
    kSKmUndefined,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,                            
    kSKmUndefined,               
#endif             
    kSKmUndefined, 
    kSKmUndefined, 
    kSKmUndefined, 
    //    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,  
    kSKmUndefined, 
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //   kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    //    kSKmUndefined,                            
    kSKmUndefined,               
#endif  
    kSKmUndefined,  
    kSKmUndefined,  
    kSKmUndefined,
    //    kSKmUndefined,
    kSKmUndefined,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined
  }

  ,

  // 1 muon-like ring (numubar disappearance) with reverse horn current                                                                                                                         
  {
    kSKmRHCNuMuCCQE,
    kSKmRHCNuMuCC1pi,
    kSKmRHCNuMuCCcoh,
    kSKmRHCNuMuCCoth,
    kSKmRHCNuMuNC1piC,
    kSKmRHCNuMuNCoth_mu,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmRHCOscNuMuCCQE,
    kSKmRHCOscNuMuCC1pi,
    kSKmRHCOscNuMuCCcoh,
    kSKmRHCOscNuMuCCoth,
#endif
    kSKmRHCNuMuBarCCQE,
    kSKmRHCNuMuBarCC1pi,
    kSKmRHCNuMuBarCCcoh,
    kSKmRHCNuMuBarCCoth,
    kSKmRHCNuMuBarNC1piC,
    kSKmRHCNuMuBarNCoth_mu,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmRHCOscNuMuBarCCQE,
    kSKmRHCOscNuMuBarCC1pi,
    kSKmRHCOscNuMuBarCCcoh,
    kSKmRHCOscNuMuBarCCoth,
#endif
    kSKmRHCNuECCQE,
    kSKmRHCNuECC1pi,
    kSKmRHCNuECCcoh,
    kSKmRHCNuECCoth,
    kSKmRHCNuENC1piC,
    kSKmRHCNuENCoth_mu,
    kSKmRHCNuEBarCCQE,
    kSKmRHCNuEBarCC1pi,
    kSKmRHCNuEBarCCcoh,
    kSKmRHCNuEBarCCoth,
    kSKmRHCNuEBarNC1piC,
    kSKmRHCNuEBarNCoth_mu,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmRHCOscNuECCQE,
    kSKmRHCOscNuECC1pi,
    kSKmRHCOscNuECCcoh,
    kSKmRHCOscNuECCoth,
#endif
    kSKmRHCSigNuEBarCCQE,
    kSKmRHCSigNuEBarCC1pi,
    kSKmRHCSigNuEBarCCcoh,
    kSKmRHCSigNuEBarCCoth,
#ifdef USE_MEC
    kSKmRHCMEC,
#else
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
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
  }

  ,

  // 1 electron-like ring (nuebar appearance) with reverse horn current                                                                                                                         
  {
    kSKmRHCNuMuCCQE,
    kSKmRHCNuMuCC1pi,
    kSKmRHCNuMuCCcoh,
    kSKmRHCNuMuCCoth,
    kSKmRHCNuMuNC1pi0,
    kSKmRHCNuMuNCcoh,
    kSKmRHCNuMuNCoth_e,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmRHCOscNuMuCCQE,
    kSKmRHCOscNuMuCC1pi,
    kSKmRHCOscNuMuCCcoh,
    kSKmRHCOscNuMuCCoth,
#endif
    kSKmRHCNuMuBarCCQE,
    kSKmRHCNuMuBarCC1pi,
    kSKmRHCNuMuBarCCcoh,
    kSKmRHCNuMuBarCCoth,
    kSKmRHCNuMuBarNC1pi0,
    kSKmRHCNuMuBarNCcoh,
    kSKmRHCNuMuBarNCoth_e,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmRHCOscNuMuBarCCQE,
    kSKmRHCOscNuMuBarCC1pi,
    kSKmRHCOscNuMuBarCCcoh,
    kSKmRHCOscNuMuBarCCoth,
#endif
    kSKmRHCNuECCQE,
    kSKmRHCNuECC1pi,
    kSKmRHCNuECCcoh,
    kSKmRHCNuECCoth,
    kSKmRHCNuENC1pi0,
    kSKmRHCNuENCcoh,
    kSKmRHCNuENCoth_e,
    kSKmRHCNuEBarCCQE,
    kSKmRHCNuEBarCC1pi,
    kSKmRHCNuEBarCCcoh,
    kSKmRHCNuEBarCCoth,
    kSKmRHCNuEBarNC1pi0,
    kSKmRHCNuEBarNCcoh,
    kSKmRHCNuEBarNCoth_e,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined
  }
  
#endif  
#endif
#endif
};



//
// -----------------------------------------------------------------------
//
// Info on flavours and CC/NC modes
//
// -----------------------------------------------------------------------
//

 static const bool kIsCC[kns][knm]= {
#if defined __FHC_RUNNING__ && defined __RHC_RUNNING__
   // > 1 muon-like ring
   { 1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
#ifdef USE_MEC
     1,
#else
     -1,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
    -1
   }
   ,
   // > 1 e-like ring
   {
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     0,
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
    -1

   }

   ,

   // > 1 muon-like ring with reverse horn current
   { 1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
#ifdef USE_MEC
     1,
#else
     -1,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
    -1
   }
   ,
   // > 1 e-like ring with reverse horn current
   {
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     0,
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
    -1

   }

#else

#ifdef __FHC_RUNNING__

   // > 1 muon-like ring
   { 1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
#ifdef USE_MEC
     1,
#else
     -1,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
    -1
   }
   ,
   // > 1 e-like ring
   {
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     0,
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
    -1

   }

   ,

   // > 1 muon-like ring with reverse horn current
   { 1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
#ifdef USE_MEC
     1,
#else
     -1,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
    -1
   }
   ,
   // > 1 e-like ring with reverse horn current
   {
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     0,
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
    -1

   }

#else

#ifdef __RHC_RUNNING__

   // > 1 muon-like ring
   { 1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
#ifdef USE_MEC
     1,
#else
     -1,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
    -1
   }
   ,
   // > 1 e-like ring
   {
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     0,
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
    -1
   }
   ,

   // > 1 muon-like ring with reverse horn current
   { 1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     1,
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
#ifdef USE_MEC
     1,
#else
     -1,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
    -1
   }
   ,
   // > 1 e-like ring with reverse horn current
   {
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
     0,
     0,
     0,
     1,
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1,
     1,
     1,
     1,
#endif
     1,
     1,
     1,
     1,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     -1,
     -1,
     -1,
     -1,
#endif
     -1,
     -1,
     -1,
     -1,
     -1,
     -1,
    -1
   }

#endif
#endif
#endif 
};



} // conv
} // t2kana
  
#endif 

