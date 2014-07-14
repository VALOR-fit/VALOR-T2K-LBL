//____________________________________________________________________________
/*!

\class    valor::conv

\brief    Conventions

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

Changes since 2011_v0.9.0:

1702-2013  Nick  
            Split valor::conv namespace into several smaller files

3102-2013  Nick
Add samples and modes for reverse horn current running

*/
//____________________________________________________________________________

#ifndef samples_modes_h
#define samples_modes_h

#include "Enumerations.h"
#include "AnalysisOptions.h"

namespace valor {
  
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

#ifdef __FHC_RUNNING__
#ifdef __NUMUONLY__
  kSKsFCFV1RingMuLike,
  kSKsUndefined,
#else
#ifdef __NUEONLY__
  kSKsUndefined,
  kSKsFCFV1RingELike,
#else
  kSKsFCFV1RingMuLike,
  kSKsFCFV1RingELike,
#endif
#endif
#else
  kSKsUndefined,
  kSKsUndefined,
#endif

#ifdef __RHC_RUNNING__
#ifdef __NUMUONLY__
  kSKsRHCFCFV1RingMuLike,
  kSKsUndefined,
#else
#ifdef __NUEONLY__
  kSKsUndefined,
  kSKsRHCFCFV1RingELike,
#else
  kSKsRHCFCFV1RingMuLike,
  kSKsRHCFCFV1RingELike,
#endif
#endif
#else
  kSKsUndefined,
  kSKsUndefined,
#endif

#ifdef __MULTIRING_FHC_RUNNING__
#ifdef __NUMUONLY__
  kSKsFCFVMultiRingMuLike,
  kSKsUndefined,
#else
#ifdef __NUEONLY__
  kSKsUndefined,
  kSKsFCFVMultiRingELike,
#else
  kSKsFCFVMultiRingMuLike,
  kSKsFCFVMultiRingELike,
#endif
#endif
#else
  kSKsUndefined,
  kSKsUndefined,
#endif

};

//
// Define which SK modes are considered for each SK sample included in the analysis
//

static const SuperkMode_t kSKReactionModes[kns][knm] = {

#ifdef __FHC_RUNNING__
  // > 1 muon-like ring (numu disappearance)
  {
    kSKmNuMuCCQE,         
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    kSKmNuMuCCoth,                      
    kSKmNuMuNC1piC,
    kSKmNuMuNCoth_mu,          
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmOscNuMuCCQE,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCcoh,
    kSKmOscNuMuCCoth,
#endif
    kSKmNuMuBarCCQE,      
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCcoh,
    kSKmNuMuBarCCoth,                  
    kSKmNuMuBarNC1piC,
    kSKmNuMuBarNCoth_mu,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCQE,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCcoh,
    kSKmOscNuMuBarCCoth,
#endif
    kSKmNuECCQE, 
    kSKmNuECC1pi,
    kSKmNuECCcoh,
    kSKmNuECCoth,                                
    kSKmNuENC1piC,
    kSKmNuENCoth_mu,
    kSKmOscNuECCQE,
    kSKmOscNuECC1pi,
    kSKmOscNuECCcoh,
    kSKmOscNuECCoth,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1piC,
    kSKmNuEBarNCoth_mu,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmOscNuEBarCCQE,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCcoh,
    kSKmOscNuEBarCCoth,
#endif
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
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
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
  }
  ,
  // > 1 electron-like ring (nue appearance)
  {
    kSKmNuMuCCQE,
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    kSKmNuMuCCoth, 
    kSKmNuMuNC1pi0, 
    kSKmNuMuNCcoh, 
    kSKmNuMuNCoth_e, 
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__    
    kSKmOscNuMuCCQE,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCcoh,
    kSKmOscNuMuCCoth,               
#endif
    kSKmNuMuBarCCQE, 
    kSKmNuMuBarCC1pi, 
    kSKmNuMuBarCCcoh, 
    kSKmNuMuBarCCoth,  
    kSKmNuMuBarNC1pi0, 
    kSKmNuMuBarNCcoh, 
    kSKmNuMuBarNCoth_e,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCQE,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCcoh,
    kSKmOscNuMuBarCCoth,               
#endif             
    kSKmNuECCQE, 
    kSKmNuECC1pi, 
    kSKmNuECCcoh, 
    kSKmNuECCoth,
    kSKmNuENC1pi0,  
    kSKmNuENCcoh, 
    kSKmNuENCoth_e,
    kSKmOscNuECCQE,  
    kSKmOscNuECC1pi,  
    kSKmOscNuECCcoh,
    kSKmOscNuECCoth,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1pi0,
    kSKmNuEBarNCcoh,
    kSKmNuEBarNCoth_e,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmOscNuEBarCCQE,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCcoh,
    kSKmOscNuEBarCCoth,               
#endif          
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
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
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
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
    kSKmUndefined
  }

#else //not using FHC
  // > 1 muon-like ring (numu disappearance)
  {
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined
  }
  ,
  // > 1 electron-like ring (nue appearance)
  {
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined
  }
#endif

  //FHC arrays now complete
  ,

#ifdef __RHC_RUNNING__
  // > 1 muon-like ring (numu disappearance) (RHC)
  {
    kSKmNuMuCCQE,         
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    kSKmNuMuCCoth,                      
    kSKmNuMuNC1piC,
    kSKmNuMuNCoth_mu,          
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmOscNuMuCCQE,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCcoh,
    kSKmOscNuMuCCoth,
#endif
    kSKmNuMuBarCCQE,      
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCcoh,
    kSKmNuMuBarCCoth,                  
    kSKmNuMuBarNC1piC,
    kSKmNuMuBarNCoth_mu,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCQE,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCcoh,
    kSKmOscNuMuBarCCoth,
#endif
    kSKmNuECCQE, 
    kSKmNuECC1pi,
    kSKmNuECCcoh,
    kSKmNuECCoth,                                
    kSKmNuENC1piC,
    kSKmNuENCoth_mu,
    kSKmOscNuECCQE,
    kSKmOscNuECC1pi,
    kSKmOscNuECCcoh,
    kSKmOscNuECCoth,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1piC,
    kSKmNuEBarNCoth_mu,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmOscNuEBarCCQE,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCcoh,
    kSKmOscNuEBarCCoth,
#endif
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
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
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
  }
  ,
  // > 1 electron-like ring (nue appearance) (RHC)
  {
        kSKmNuMuCCQE,
        kSKmNuMuCC1pi,
        kSKmNuMuCCcoh,
        kSKmNuMuCCoth,
        kSKmNuMuNC1pi0,
        kSKmNuMuNCcoh,
        kSKmNuMuNCoth_e,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        kSKmOscNuMuCCQE,
        kSKmOscNuMuCC1pi,
        kSKmOscNuMuCCcoh,
        kSKmOscNuMuCCoth,
#endif
        kSKmNuMuBarCCQE,
        kSKmNuMuBarCC1pi,
        kSKmNuMuBarCCcoh,
        kSKmNuMuBarCCoth,
        kSKmNuMuBarNC1pi0,
        kSKmNuMuBarNCcoh,
        kSKmNuMuBarNCoth_e,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        kSKmOscNuMuBarCCQE,
        kSKmOscNuMuBarCC1pi,
        kSKmOscNuMuBarCCcoh,
        kSKmOscNuMuBarCCoth,
#endif
        kSKmNuECCQE,
        kSKmNuECC1pi,
        kSKmNuECCcoh,
        kSKmNuECCoth,
        kSKmNuENC1pi0,
        kSKmNuENCcoh,
        kSKmNuENCoth_e,
        kSKmOscNuECCQE,
        kSKmOscNuECC1pi,
        kSKmOscNuECCcoh,
        kSKmOscNuECCoth,
        kSKmNuEBarCCQE,
        kSKmNuEBarCC1pi,
        kSKmNuEBarCCcoh,
        kSKmNuEBarCCoth,
        kSKmNuEBarNC1pi0,
        kSKmNuEBarNCcoh,
        kSKmNuEBarNCoth_e,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
        kSKmOscNuEBarCCQE,
        kSKmOscNuEBarCC1pi,
        kSKmOscNuEBarCCcoh,
        kSKmOscNuEBarCCoth,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
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
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
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
        kSKmUndefined
    }

#else //not using RHC
  // > 1 muon-like ring (numu disappearance) (RHC)
  {
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined
  }
  ,
  // > 1 electron-like ring (nue appearance) (RHC)
  {
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined
  }
#endif  
  
  //RHC arrays now complete


  ,


#ifdef __MULTIRING_FHC_RUNNING__
  // > 1 muon-like ring (numu disappearance)
  {
    kSKmNuMuCCQE,         
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    kSKmNuMuCCoth,                      
    kSKmNuMuNC1piC,
    kSKmNuMuNCoth_mu,          
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmOscNuMuCCQE,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCcoh,
    kSKmOscNuMuCCoth,
#endif
    kSKmNuMuBarCCQE,      
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCcoh,
    kSKmNuMuBarCCoth,                  
    kSKmNuMuBarNC1piC,
    kSKmNuMuBarNCoth_mu,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCQE,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCcoh,
    kSKmOscNuMuBarCCoth,
#endif
    kSKmNuECCQE, 
    kSKmNuECC1pi,
    kSKmNuECCcoh,
    kSKmNuECCoth,                                
    kSKmNuENC1piC,
    kSKmNuENCoth_mu,
    kSKmOscNuECCQE,
    kSKmOscNuECC1pi,
    kSKmOscNuECCcoh,
    kSKmOscNuECCoth,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1piC,
    kSKmNuEBarNCoth_mu,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmOscNuEBarCCQE,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCcoh,
    kSKmOscNuEBarCCoth,
#endif
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
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
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
  }
  ,
  // > 1 electron-like ring (nue appearance)
  {
    kSKmNuMuCCQE,
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    kSKmNuMuCCoth, 
    kSKmNuMuNC1pi0, 
    kSKmNuMuNCcoh, 
    kSKmNuMuNCoth_e, 
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__    
    kSKmOscNuMuCCQE,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCcoh,
    kSKmOscNuMuCCoth,               
#endif
    kSKmNuMuBarCCQE, 
    kSKmNuMuBarCC1pi, 
    kSKmNuMuBarCCcoh, 
    kSKmNuMuBarCCoth,  
    kSKmNuMuBarNC1pi0, 
    kSKmNuMuBarNCcoh, 
    kSKmNuMuBarNCoth_e,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCQE,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCcoh,
    kSKmOscNuMuBarCCoth,               
#endif             
    kSKmNuECCQE, 
    kSKmNuECC1pi, 
    kSKmNuECCcoh, 
    kSKmNuECCoth,
    kSKmNuENC1pi0,  
    kSKmNuENCcoh, 
    kSKmNuENCoth_e,
    kSKmOscNuECCQE,  
    kSKmOscNuECC1pi,  
    kSKmOscNuECCcoh,
    kSKmOscNuECCoth,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1pi0,
    kSKmNuEBarNCcoh,
    kSKmNuEBarNCoth_e,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    kSKmOscNuEBarCCQE,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCcoh,
    kSKmOscNuEBarCCoth,               
#endif          
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
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
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
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
    kSKmUndefined
  }

#else //not using MultiRing FHC
  // > 1 muon-like ring (numu disappearance)
  {
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined
  }
  ,
  // > 1 electron-like ring (nue appearance)
  {
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
    kSKmUndefined, kSKmUndefined, kSKmUndefined
  }
#endif

  //MultiRing FHC arrays now complete
  
};

//
// -----------------------------------------------------------------------
//
// Info on flavours and CC/NC modes
//
// -----------------------------------------------------------------------
//

static const bool kIsCC[kns][knm]= {

  // > 1 muon-like ring
  {  1, //numu
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1, //oscnumu
     1,
     1,
     1,
#endif
     1, //numubar
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1, //oscnumubar
     1,
     1,
     1,
#endif
     1, //nue
     1,
     1,
     1,
     0,
     0,
     1, //oscnue
     1,
     1,
     1,
     1, //nuebar
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1, //oscnuebar
     1,
     1,
     1,
#endif
#ifdef USE_MEC
     1, //mec
#else
     0,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     0,
     0,
     0,
     0,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     0,
     0,
     0,
     0,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     0,
     0,
     0,
     0,
#endif
     0, //always undefined (12)
     0,
     0,
     0,
     0,
     0,
     0,
     0,
     0,
     0,
     0,
     0
  }
  ,
  // > 1 e-like ring
  {  1, //numu
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1, //oscnumu
     1,
     1,
     1,
#endif
     1, //numubar
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1, //oscnumubar
     1,
     1,
     1,
#endif
     1, //nue
     1,
     1,
     1,
     0,
     0,
     0,
     1, //oscnue
     1,
     1,
     1,
     1, //nuebar
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1, //oscnuebar
     1,
     1,
     1,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
#endif
      0, //always undefined (9)
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
      0  //undefined
  }
  ,
  // > 1 muon-like ring (RHC)
  {  1, //numu
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1, //oscnumu
     1,
     1,
     1,
#endif
     1, //numubar
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1, //oscnumubar
     1,
     1,
     1,
#endif
     1, //nue
     1,
     1,
     1,
     0,
     0,
     1, //oscnue
     1,
     1,
     1,
     1, //nuebar
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1, //oscnuebar
     1,
     1,
     1,
#endif
#ifdef USE_MEC
     1, //mec
#else
     0,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     0,
     0,
     0,
     0,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     0,
     0,
     0,
     0,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     0,
     0,
     0,
     0,
#endif
     0, //always undefined (12)
     0,
     0,
     0,
     0,
     0,
     0,
     0,
     0,
     0,
     0,
     0
  }
  ,
  // > 1 e-like ring (RHC)
    {  1, //numu
        1,
        1,
        1,
        0,
        0,
        0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        1, //oscnumu
        1,
        1,
        1,
#endif
        1, //numubar
        1,
        1,
        1,
        0,
        0,
        0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        1, //oscnumubar
        1,
        1,
        1,
#endif
        1, //nue
        1,
        1,
        1,
        0,
        0,
        0,
        1, //oscnue
        1,
        1,
        1,
        1, //nuebar
        1,
        1,
        1,
        0,
        0,
        0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
        1, //oscnuebar
        1,
        1,
        1,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        0, //undefined
        0, //undefined
        0, //undefined
        0, //undefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        0, //undefined
        0, //undefined
        0, //undefined
        0, //undefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
        0, //undefined
        0, //undefined
        0, //undefined
        0, //undefined
#endif
        0, //always undefined (9)
        0, //undefined
        0, //undefined
        0, //undefined
        0, //undefined
        0, //undefined
        0, //undefined
        0, //undefined
        0  //undefined
    }
  ,
  // > multi muon-like ring (FHC)
  {  1, //numu
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1, //oscnumu
     1,
     1,
     1,
#endif
     1, //numubar
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1, //oscnumubar
     1,
     1,
     1,
#endif
     1, //nue
     1,
     1,
     1,
     0,
     0,
     1, //oscnue
     1,
     1,
     1,
     1, //nuebar
     1,
     1,
     1,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1, //oscnuebar
     1,
     1,
     1,
#endif
#ifdef USE_MEC
     1, //mec
#else
     0,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
     0,
     0,
     0,
     0,
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     0,
     0,
     0,
     0,
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
     0,
     0,
     0,
     0,
#endif
     0, //always undefined (12)
     0,
     0,
     0,
     0,
     0,
     0,
     0,
     0,
     0,
     0,
     0
  }
  ,
  // > multi e-like ring
  {  1, //numu
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
     1, //oscnumu
     1,
     1,
     1,
#endif
     1, //numubar
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
     1, //oscnumubar
     1,
     1,
     1,
#endif
     1, //nue
     1,
     1,
     1,
     0,
     0,
     0,
     1, //oscnue
     1,
     1,
     1,
     1, //nuebar
     1,
     1,
     1,
     0,
     0,
     0,
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
     1, //oscnuebar
     1,
     1,
     1,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
#endif
      0, //always undefined (9)
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
      0, //undefined
      0  //undefined
  }

};

} // conv
} // valor
  
#endif 

