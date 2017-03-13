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

  21-7-2014 Davide
  Added samples for multiring forward horn current

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
#  if defined (__NUMUONLY__)  || defined(__MULTIRING_FHC_NUMU__) && !defined(__MULTIRING_FHC_NUE__)
      kSKsFCFVMultiRingMuLike,
      kSKsUndefined,
#  elif defined (__NUEONLY__) || defined(__MULTIRING_FHC_NUE__) && !defined(__MULTIRING_FHC_NUMU__)
      kSKsUndefined,
      kSKsFCFVMultiRingELike,
#  else
      kSKsFCFVMultiRingMuLike,
      kSKsFCFVMultiRingELike,
#  endif
      
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
        kSKmNuMuCC2p2h,
        kSKmNuMuCCoth,
        kSKmNuMuNC1piC,
        kSKmNuMuNC1pi0,
        kSKmNuMuNCcoh,
        kSKmNuMuNCoth,
        kSKmNuMuNC1Gamma,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        kSKmOscNuMuCCQE,
        kSKmOscNuMuCC1pi,
        kSKmOscNuMuCCcoh,
        kSKmOscNuMuCC2p2h,
        kSKmOscNuMuCCoth,
#endif
        kSKmNuMuBarCCQE,
        kSKmNuMuBarCC1pi,
        kSKmNuMuBarCCcoh,
        kSKmNuMuBarCC2p2h,
        kSKmNuMuBarCCoth,
        kSKmNuMuBarNC1piC,
        kSKmNuMuBarNC1pi0,
        kSKmNuMuBarNCcoh,
        kSKmNuMuBarNCoth,
        kSKmNuMuBarNC1Gamma,

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        kSKmOscNuMuBarCCQE,
        kSKmOscNuMuBarCC1pi,
        kSKmOscNuMuBarCCcoh,
        kSKmOscNuMuBarCC2p2h,
        kSKmOscNuMuBarCCoth,
#endif
        kSKmNuECCQE,
        kSKmNuECC1pi,
        kSKmNuECCcoh,
        kSKmNuECC2p2h,
        kSKmNuECCoth,
        kSKmNuENC1piC,
        kSKmNuENC1pi0,
        kSKmNuENCcoh,
        kSKmNuENCoth,
        kSKmNuENC1Gamma,

        kSKmOscNuECCQE,
        kSKmOscNuECC1pi,
        kSKmOscNuECCcoh,
        kSKmOscNuECC2p2h,
        kSKmOscNuECCoth,
        kSKmNuEBarCCQE,
        kSKmNuEBarCC1pi,
        kSKmNuEBarCCcoh,
        kSKmNuEBarCC2p2h,
        kSKmNuEBarCCoth,
        kSKmNuEBarNC1piC,
        kSKmNuEBarNC1pi0,
        // NEW FHC+RHC JOINT
        kSKmNuEBarNCcoh,
        kSKmNuEBarNCoth,
        kSKmNuEBarNC1Gamma,

        //
        kSKmOscNuEBarCCQE,
        kSKmOscNuEBarCC1pi,
        kSKmOscNuEBarCCcoh,
        kSKmOscNuEBarCC2p2h,
        kSKmOscNuEBarCCoth,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        kSKmUndefined,
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
        kSKmUndefined,
#endif
      }
      ,
      // > 1 electron-like ring (nue appearance)
      {
        kSKmNuMuCCQE,
        kSKmNuMuCC1pi,
        kSKmNuMuCCcoh,
        kSKmNuMuCC2p2h,
        kSKmNuMuCCoth,
        kSKmNuMuNC1piC,
        kSKmNuMuNC1pi0,
        kSKmNuMuNCcoh,
        kSKmNuMuNCoth,
        kSKmNuMuNC1Gamma,

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        kSKmOscNuMuCCQE,
        kSKmOscNuMuCC1pi,
        kSKmOscNuMuCCcoh,
        kSKmOscNuMuCC2p2h,
        kSKmOscNuMuCCoth,
#endif
        kSKmNuMuBarCCQE,
        kSKmNuMuBarCC1pi,
        kSKmNuMuBarCCcoh,
        kSKmNuMuBarCC2p2h,
        kSKmNuMuBarCCoth,
        kSKmNuMuBarNC1piC,
        kSKmNuMuBarNC1pi0,
        kSKmNuMuBarNCcoh,
        kSKmNuMuBarNCoth,
        kSKmNuMuBarNC1Gamma,

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        kSKmOscNuMuBarCCQE,
        kSKmOscNuMuBarCC1pi,
        kSKmOscNuMuBarCCcoh,
        kSKmOscNuMuBarCC2p2h,
        kSKmOscNuMuBarCCoth,
#endif
        kSKmNuECCQE,
        kSKmNuECC1pi,
        kSKmNuECCcoh,
        kSKmNuECC2p2h,
        kSKmNuECCoth,
        kSKmNuENC1piC,
        kSKmNuENC1pi0,
        kSKmNuENCcoh,
        kSKmNuENCoth,
        kSKmNuENC1Gamma,

        kSKmOscNuECCQE,
        kSKmOscNuECC1pi,
        kSKmOscNuECCcoh,
        kSKmOscNuECC2p2h,
        kSKmOscNuECCoth,
        kSKmNuEBarCCQE,
        kSKmNuEBarCC1pi,
        kSKmNuEBarCCcoh,
        kSKmNuEBarCC2p2h,
        kSKmNuEBarCCoth,
        kSKmNuEBarNC1piC,
        kSKmNuEBarNC1pi0,
        kSKmNuEBarNCcoh,
        kSKmNuEBarNCoth,
        kSKmNuEBarNC1Gamma,

        kSKmOscNuEBarCCQE,
        kSKmOscNuEBarCC1pi,
        kSKmOscNuEBarCCcoh,
        kSKmOscNuEBarCC2p2h,
        kSKmOscNuEBarCCoth,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        kSKmUndefined,
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
        kSKmUndefined,
#endif
      }

#else //not using FHC
      // > 1 muon-like ring (numu disappearance)
      {
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined
      }
      ,
      // > 1 electron-like ring (nue appearance)
      {
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined
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
        kSKmNuMuCC2p2h,
        kSKmNuMuCCoth,
        kSKmNuMuNC1piC,
        kSKmNuMuNC1pi0,
        kSKmNuMuNCcoh,
        kSKmNuMuNCoth,
        kSKmNuMuNC1Gamma,

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        kSKmOscNuMuCCQE,
        kSKmOscNuMuCC1pi,
        kSKmOscNuMuCCcoh,
        kSKmOscNuMuCC2p2h,
        kSKmOscNuMuCCoth,
#endif
        kSKmNuMuBarCCQE,
        kSKmNuMuBarCC1pi,
        kSKmNuMuBarCCcoh,
        kSKmNuMuBarCC2p2h,
        kSKmNuMuBarCCoth,
        kSKmNuMuBarNC1piC,
        kSKmNuMuBarNC1pi0,
        kSKmNuMuBarNCcoh,
        kSKmNuMuBarNCoth,
        kSKmNuMuBarNC1Gamma,

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        kSKmOscNuMuBarCCQE,
        kSKmOscNuMuBarCC1pi,
        kSKmOscNuMuBarCCcoh,
        kSKmOscNuMuBarCC2p2h,
        kSKmOscNuMuBarCCoth,
#endif
        kSKmNuECCQE,
        kSKmNuECC1pi,
        kSKmNuECCcoh,
        kSKmNuECC2p2h,
        kSKmNuECCoth,
        kSKmNuENC1piC,
        kSKmNuENC1pi0,
        kSKmNuENCcoh,
        kSKmNuENCoth,
        kSKmNuENC1Gamma,

        kSKmOscNuECCQE,
        kSKmOscNuECC1pi,
        kSKmOscNuECCcoh,
        kSKmOscNuECC2p2h,
        kSKmOscNuECCoth,
        kSKmNuEBarCCQE,
        kSKmNuEBarCC1pi,
        kSKmNuEBarCCcoh,
        kSKmNuEBarCC2p2h,
        kSKmNuEBarCCoth,
        kSKmNuEBarNC1piC,
        kSKmNuEBarNC1pi0,
        kSKmNuEBarNCcoh,
        kSKmNuEBarNCoth,
        kSKmNuEBarNC1Gamma,

        kSKmOscNuEBarCCQE,
        kSKmOscNuEBarCC1pi,
        kSKmOscNuEBarCCcoh,
        kSKmOscNuEBarCC2p2h,
        kSKmOscNuEBarCCoth,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        kSKmUndefined,
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
        kSKmUndefined,
#endif
      }
      ,
      // > 1 electron-like ring (nue appearance) (RHC)
      {
        kSKmNuMuCCQE,
        kSKmNuMuCC1pi,
        kSKmNuMuCCcoh,
        kSKmNuMuCC2p2h,
        kSKmNuMuCCoth,
        kSKmNuMuNC1piC,
        kSKmNuMuNC1pi0,
        kSKmNuMuNCcoh,
        kSKmNuMuNCoth,
        kSKmNuMuNC1Gamma,

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        kSKmOscNuMuCCQE,
        kSKmOscNuMuCC1pi,
        kSKmOscNuMuCCcoh,
        kSKmOscNuMuCC2p2h,
        kSKmOscNuMuCCoth,
#endif
        kSKmNuMuBarCCQE,
        kSKmNuMuBarCC1pi,
        kSKmNuMuBarCCcoh,
        kSKmNuMuBarCC2p2h,
        kSKmNuMuBarCCoth,
        kSKmNuMuBarNC1piC,
        kSKmNuMuBarNC1pi0,
        kSKmNuMuBarNCcoh,
        kSKmNuMuBarNCoth,
        kSKmNuMuBarNC1Gamma,

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        kSKmOscNuMuBarCCQE,
        kSKmOscNuMuBarCC1pi,
        kSKmOscNuMuBarCCcoh,
        kSKmOscNuMuBarCC2p2h,
        kSKmOscNuMuBarCCoth,
#endif
        kSKmNuECCQE,
        kSKmNuECC1pi,
        kSKmNuECCcoh,
        kSKmNuECC2p2h,
        kSKmNuECCoth,
        kSKmNuENC1piC,
        kSKmNuENC1pi0,
        kSKmNuENCcoh,
        kSKmNuENCoth,
        kSKmNuENC1Gamma,

        kSKmOscNuECCQE,
        kSKmOscNuECC1pi,
        kSKmOscNuECCcoh,
        kSKmOscNuECC2p2h,
        kSKmOscNuECCoth,
        kSKmNuEBarCCQE,
        kSKmNuEBarCC1pi,
        kSKmNuEBarCCcoh,
        kSKmNuEBarCC2p2h,
        kSKmNuEBarCCoth,
        kSKmNuEBarNC1piC,
        kSKmNuEBarNC1pi0,
        kSKmNuEBarNCcoh,
        kSKmNuEBarNCoth,
        kSKmNuEBarNC1Gamma,

        kSKmOscNuEBarCCQE,
        kSKmOscNuEBarCC1pi,
        kSKmOscNuEBarCCcoh,
        kSKmOscNuEBarCC2p2h,
        kSKmOscNuEBarCCoth,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        kSKmUndefined,
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
        kSKmUndefined,
#endif
      }

#else //not using RHC
      // > 1 muon-like ring (numu disappearance) (RHC)
      {
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined
      }
      ,
      // > 1 electron-like ring (nue appearance) (RHC)
      {
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined
      }
#endif

      //RHC arrays now complete


      ,


#ifdef __MULTIRING_FHC_RUNNING__
#ifdef __MULTIRING_FHC_NUMU__
      // > MultiRmu (1 muon-like ring (numu disappearance))
      {
        kSKmNuMuCCQE,
        kSKmNuMuCC1pi,
        kSKmNuMuCCcoh,
        kSKmNuMuCC2p2h,
        kSKmNuMuCCoth,
        kSKmNuMuNC1piC,
        kSKmNuMuNC1pi0,
        kSKmNuMuNCcoh,
        kSKmNuMuNCoth,
        kSKmNuMuNC1Gamma,

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        kSKmOscNuMuCCQE,
        kSKmOscNuMuCC1pi,
        kSKmOscNuMuCCcoh,
        kSKmOscNuMuCC2p2h,
        kSKmOscNuMuCCoth,
#endif
        kSKmNuMuBarCCQE,
        kSKmNuMuBarCC1pi,
        kSKmNuMuBarCCcoh,
        kSKmNuMuBarCC2p2h,
        kSKmNuMuBarCCoth,
        kSKmNuMuBarNC1piC,
        kSKmNuMuBarNC1pi0,
        kSKmNuMuBarNCcoh,
        kSKmNuMuBarNCoth,
        kSKmNuMuBarNC1Gamma,

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        kSKmOscNuMuBarCCQE,
        kSKmOscNuMuBarCC1pi,
        kSKmOscNuMuBarCCcoh,
        kSKmOscNuMuBarCC2p2h,
        kSKmOscNuMuBarCCoth,
#endif
        kSKmNuECCQE,
        kSKmNuECC1pi,
        kSKmNuECCcoh,
        kSKmNuECC2p2h,
        kSKmNuECCoth,
        kSKmNuENC1piC,
        kSKmNuENC1pi0,
        kSKmNuENCcoh,
        kSKmNuENCoth,
        kSKmNuENC1Gamma,

        kSKmOscNuECCQE,
        kSKmOscNuECC1pi,
        kSKmOscNuECCcoh,
        kSKmOscNuECC2p2h,
        kSKmOscNuECCoth,
        kSKmNuEBarCCQE,
        kSKmNuEBarCC1pi,
        kSKmNuEBarCCcoh,
        kSKmNuEBarCC2p2h,
        kSKmNuEBarCCoth,
        kSKmNuEBarNC1piC,
        kSKmNuEBarNC1pi0,
        kSKmNuEBarNCcoh,
        kSKmNuEBarNCoth,
        kSKmNuEBarNC1Gamma,

        kSKmOscNuEBarCCQE,
        kSKmOscNuEBarCC1pi,
        kSKmOscNuEBarCCcoh,
        kSKmOscNuEBarCC2p2h,
        kSKmOscNuEBarCCoth,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        kSKmUndefined,
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
        kSKmUndefined,
#endif
      }
      ,
#endif // MultiRmu
#if !defined(__MULTIRING_FHC_NUMU__)
      {
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined
      }
      ,      
#endif      
#ifdef __MULTIRING_FHC_NUE__      
      // > MultiRe (1 electron-like ring (nue appearance))
      {
        kSKmNuMuCCQE,
        kSKmNuMuCC1pi,
        kSKmNuMuCCcoh,
        kSKmNuMuCC2p2h,
        kSKmNuMuCCoth,
        kSKmNuMuNC1piC,
        kSKmNuMuNC1pi0,
        kSKmNuMuNCcoh,
        kSKmNuMuNCoth,
        kSKmNuMuNC1Gamma,

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
        kSKmOscNuMuCCQE,
        kSKmOscNuMuCC1pi,
        kSKmOscNuMuCCcoh,
        kSKmOscNuMuCC2p2h,
        kSKmOscNuMuCCoth,
#endif
        kSKmNuMuBarCCQE,
        kSKmNuMuBarCC1pi,
        kSKmNuMuBarCCcoh,
        kSKmNuMuBarCC2p2h,
        kSKmNuMuBarCCoth,
        kSKmNuMuBarNC1piC,
        kSKmNuMuBarNC1pi0,
        kSKmNuMuBarNCcoh,
        kSKmNuMuBarNCoth,
        kSKmNuMuBarNC1Gamma,

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
        kSKmOscNuMuBarCCQE,
        kSKmOscNuMuBarCC1pi,
        kSKmOscNuMuBarCCcoh,
        kSKmOscNuMuBarCC2p2h,
        kSKmOscNuMuBarCCoth,
#endif
        kSKmNuECCQE,
        kSKmNuECC1pi,
        kSKmNuECCcoh,
        kSKmNuECC2p2h,
        kSKmNuECCoth,
        kSKmNuENC1piC,
        kSKmNuENC1pi0,
        kSKmNuENCcoh,
        kSKmNuENCoth,
        kSKmNuENC1Gamma,

        kSKmOscNuECCQE,
        kSKmOscNuECC1pi,
        kSKmOscNuECCcoh,
        kSKmOscNuECC2p2h,
        kSKmOscNuECCoth,
        kSKmNuEBarCCQE,
        kSKmNuEBarCC1pi,
        kSKmNuEBarCCcoh,
        kSKmNuEBarCC2p2h,
        kSKmNuEBarCCoth,
        kSKmNuEBarNC1piC,
        kSKmNuEBarNC1pi0,
        kSKmNuEBarNCcoh,
        kSKmNuEBarNCoth,
        kSKmNuEBarNC1Gamma,

        kSKmOscNuEBarCCQE,
        kSKmOscNuEBarCC1pi,
        kSKmOscNuEBarCCcoh,
        kSKmOscNuEBarCC2p2h,
        kSKmOscNuEBarCCoth,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
        kSKmUndefined,
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
        kSKmUndefined,
#endif
      }
#endif // MultiRe FHC
#if !defined(__MULTIRING_FHC_NUE__)
            {
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined
      }
#endif
      
#else //not using MR FHC
      // > 1 muon-like ring (numu disappearance)
      {
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined
      }
      ,
      // > 1 electron-like ring (nue appearance)
      {
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined,
        kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined, kSKmUndefined
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
         1,
         0,
         0,
         0,
         0,
         0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
         1, //oscnumu
         1,
         1,
         1,
         1,
#endif
         1, //numubar
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
         1, //oscnumubar
         1,
         1,
         1,
         1,
#endif
         1, //nue
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
         1, //oscnue
         1,
         1,
         1,
         1,
         1, //nuebar
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
         1, //oscnuebar
         1,
         1,
         1,
         1,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
         0,
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
         0,
#endif
      }
      ,
      // > 1 e-like ring
      {  1, //numu
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
         1, //oscnumu
         1,
         1,
         1,
         1,
#endif
         1, //numubar
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
         1, //oscnumubar
         1,
         1,
         1,
         1,
#endif
         1, //nue
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
         1, //oscnue
         1,
         1,
         1,
         1,
         1, //nuebar
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
         1, //oscnuebar
         1,
         1,
         1,
         1,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
         0, //undefined
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
         0, //undefined
#endif
      }
      ,
      // > 1 muon-like ring (RHC)
      {  1, //numu
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
         1, //oscnumu
         1,
         1,
         1,
         1,
#endif
         1, //numubar
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
         1, //oscnumubar
         1,
         1,
         1,
         1,
#endif
         1, //nue
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
         1, //oscnue
         1,
         1,
         1,
         1,
         1, //nuebar
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
         1, //oscnuebar
         1,
         1,
         1,
         1,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
         0,
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
         0,
#endif
      }
      ,
      // > 1 e-like ring (RHC)
      {   1, //numu
          1,
          1,
          1,
          1,
          0,
          0,
          0,
          0,
          0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
          1, //oscnumu
          1,
          1,
          1,
          1,
#endif
          1, //numubar
          1,
          1,
          1,
          1,
          0,
          0,
          0,
          0,
          0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
          1, //oscnumubar
          1,
          1,
          1,
          1,
#endif
          1, //nue
          1,
          1,
          1,
          1,
          0,
          0,
          0,
          0,
          0,
          1, //oscnue
          1,
          1,
          1,
          1,
          1, //nuebar
          1,
          1,
          1,
          1,
          0,
          0,
          0,
          0,
          0,
          1, //oscnuebar
          1,
          1,
          1,
          1,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
          0, //undefined
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
          0, //undefined
#endif
      }
      ,
      // > multi muon-like ring (FHC)
      {  1, //numu
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
         1, //oscnumu
         1,
         1,
         1,
         1,
#endif
         1, //numubar
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
         1, //oscnumubar
         1,
         1,
         1,
         1,
#endif
         1, //nue
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
         1, //oscnue
         1,
         1,
         1,
         1,
         1, //nuebar
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
         1, //oscnuebar
         1,
         1,
         1,
         1,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
         0,
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
         0,
#endif
      }
      ,
      // > multi e-like ring
      {  1, //numu
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
         1, //oscnumu
         1,
         1,
         1,
         1,
#endif
         1, //numubar
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
         1, //oscnumubar
         1,
         1,
         1,
         1,
#endif
         1, //nue
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
         1, //oscnue
         1,
         1,
         1,
         1,
         1, //nuebar
         1,
         1,
         1,
         1,
         0,
         0,
         0,
         0,
         0,
         1, //oscnuebar
         1,
         1,
         1,
         1,
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
         0, //undefined
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
         0, //undefined
#endif
      }

    };

  } // conv
namespace util {

  
} // util
} // valor
#endif

