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

#ifndef formatting_h
#define formatting_h

#include "Binning.h"

namespace t2kana {
namespace conv   {

// Indices for array with default graphic attributes (see below)
static const int kNGraphicAttributes = 5; 
typedef enum EGraphicAttributeIdx
{   
   kGrAttUndefined  = -1, 
   kGrAttFillStyle  =  0, 
   kGrAttFillColour =  1, 
   kGrAttLineStyle  =  2, 
   kGrAttLineWidth  =  3, 
   kGrAttLineColour =  4 
} 
GraphicAttributeIdx_t;

//
// ------------------------------------------------------------------------
// Formatting etc
// ------------------------------------------------------------------------
//

//
// Default graphic attributes for output spectra. Each line refers to a mode of a given sample.
// List of samples and modes should match the one given in kSKReactionModes[kns][knm] above.
// The 5 listed attributes for each sample are: Fill style, Fill color, Line style, Line width, Line color.
//
static const int kSKHistoAttributes[kns][knm][kNGraphicAttributes] = {

  // > 1 muon-like ring 
  {

    {   -1,  kWhite,    1,  1,  kBlack}, // kSKmNuMuCCQE
    { 3544,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCC1pi
    { 3244,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCCcoh
    { 3006,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCCoth_mu
    { 3554,  kRed,      1,  1,  kBlack}, // kSKmNuMuNC1piC
    { 3006,  kRed,      1,  1,  kBlack}, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    { 3001,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCQE
    { 3544,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC1pi
    { 3244,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCcoh
    { 3006,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCoth_mu
#endif
    { 3001,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCQE
    { 3345,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCC1pi
    { 3244,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCcoh
    { 3354,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCoth_mu
    { 3544,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNC1piC
    { 3006,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    { 3001,  kAzure-5,  1,  1,  kBlack}, // kSKmNuMuBarCCQE
    { 3345,  kAzure-5,  1,  1,  kBlack}, // kSKmNuMuBarCC1pi
    { 3244,  kAzure-5,  1,  1,  kBlack}, // kSKmNuMuBarCCcoh
    { 3354,  kAzure-5,  1,  1,  kBlack}, // kSKmNuMuBarCCoth_mu
#endif

    { 3001,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCQE
    { 3345,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECC1pi
    { 3244,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCcoh
    { 3354,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCoth_mu
    { 3544,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENC1pic
    { 3006,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENCoth_mu

    { 3001,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCQE
    { 3345,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCC1pi
    { 3244,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCcoh
    { 3354,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCoth_mu
    { 3544,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNC1pic
    { 3006,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    { 3001,  kBlue-3,     1,  1,  kBlack}, // kSKmNuEBarCCQE
    { 3345,  kBlue-3,     1,  1,  kBlack}, // kSKmNuEBarCC1pi
    { 3244,  kBlue-3,     1,  1,  kBlack}, // kSKmNuEBarCCcoh
    { 3354,  kBlue-3,     1,  1,  kBlack}, // kSKmNuEBarCCoth_mu
#endif

    { 3001,  kGreen-3,  1,  1,  kBlack}, // kSKmSigNuECCQE
    { 3345,  kGreen-3,  1,  1,  kBlack}, // kSKmSigNuECC1pi
    { 3244,  kGreen-3,  1,  1,  kBlack}, // kSKmSigNuECCcoh
    { 3354,  kGreen-3,  1,  1,  kBlack}, // kSKmSigNuECCoth_mu

    { 3244, kOrange+1,  1,  1,  kBlack}, // kSKmMEC

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
#endif
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }  // kSKmUndefined
  }
  ,

  // > 1 electron-like ring 
  {
    {   -1,    kWhite,  1,  1, kBlack }, // kSKmNuMuCCQE
    { 3544,    kBlack,  1,  1, kBlack }, // kSKmNuMuCC1pi
    { 3244,    kGreen,  1,  1, kBlack }, // kSKmNuMuCCcoh
    { 3244,     kBlue,  1,  1, kBlack }, // kSKmNuMuCCoth_e
    { 3006,      kRed,  1,  1, kBlack }, // kSKmNuMuNC1pi0
    { 3006,    kRed-1,  1,  1, kBlack }, // kSKmNuMuNCcoh
    { 3006,    kRed-2,  1,  1, kBlack }, // kSKmNuMuNCoth_e

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    { 3001,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCQE     
    { 3544,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC1pi                                                                                                    
    { 3244,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCcoh 
    { 3006,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCoth_mu     
#endif

    { 3001,  kAzure-3,  1,  1, kBlack }, // kSKmNuMuBarCCQE
    { 3345,  kAzure-4,  1,  1, kBlack }, // kSKmNuMuBarCC1pi
    { 3354,  kAzure-6,  1,  1, kBlack }, // kSKmNuMuBarCCcoh
    { 3354,  kMagenta,  1,  1, kBlack }, // kSKmNuMuBarCCoth_e
    { 3006,  kAzure-7,  1,  1, kBlack }, // kSKmNuMuBarNC1pi0
    { 3006,  kAzure-9,  1,  1, kBlack }, // kSKmNuMuBarNCcoh
    { 3006,  kMagenta,  1,  1, kBlack }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    { 3001,  kAzure-7,  1,  1, kBlack }, // kSKmNuMuBarCCQE
    { 3345,  kAzure-8,  1,  1, kBlack }, // kSKmNuMuBarCC1pi
    { 3354,  kAzure-9,  1,  1, kBlack }, // kSKmNuMuBarCCcoh
    { 3354,  kMagenta-3,  1,  1, kBlack }, // kSKmNuMuBarCCoth_e
#endif

    { 3001,  kGreen+1,  1,  1, kBlack }, // kSKmNuECCQE
    { 3345,  kGreen+2,  1,  1, kBlack }, // kSKmNuECC1pi
    { 3354,  kGreen+2,  1,  1, kBlack }, // kSKmNuECCcoh
    { 3354,    kBlack,  1,  1, kBlack }, // kSKmNuECCoth_e
    { 3006,  kGreen+5,  1,  1, kBlack }, // kSKmNuENC1pi0
    { 3006,  kGreen+7,  1,  1, kBlack }, // kSKmNuENCcoh
    { 3006,    kBlack,  1,  1, kBlack }, // kSKmNuENCoth_e

    { 3001,  kBlue+1,   1,  1, kBlack }, // kSKmNuEBarCCQE
    { 3345,  kBlue+2,   1,  1, kBlack }, // kSKmNuEBarCC1pi
    { 3354,  kBlue+2,   1,  1, kBlack }, // kSKmNuEBarCCcoh
    { 3354,    kBlack,  1,  1, kBlack }, // kSKmNuEBarCCoth_e
    { 3006,  kBlue+5,   1,  1, kBlack }, // kSKmNuEBarNC1pi0
    { 3006,  kBlue+7,   1,  1, kBlack }, // kSKmNuEBarNCcoh
    { 3006,    kBlack,  1,  1, kBlack }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    { 3001,  kBlue+3,   1,  1, kBlack }, // kSKmNuEBarCCQE
    { 3345,  kBlue+5,   1,  1, kBlack }, // kSKmNuEBarCC1pi
    { 3354,  kBlue+5,   1,  1, kBlack }, // kSKmNuEBarCCcoh
    { 3354,    kBlack-3,  1,  1, kBlack }, // kSKmNuEBarCCoth_e
#endif

    { 3001,  kGreen+8,  1,  1, kBlack }, // kSKmSigNuECCQE
    { 3345, kGreen+10,  1,  1, kBlack }, // kSKmSigNuECC1pi
    { 3354,      kRed,  1,  1, kBlack }, // kSKmSigNuECCcoh
    { 3354,    kRed-2,  1,  1, kBlack }, // kSKmSigNuECCoth_e

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
#endif
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }  // kSKmUndefined

  }
  
};

//
// Define the order in which TH1D spectra are added in a THStack.
// The order is chosen to increase visibility of smaller components.
//

static const SuperkMode_t kSKReactionModesOrderInStack[kns][knm] = {

  // > 1 muon-like ring (numu disappearance)
  {
#ifdef USE_MEC
    kSKmMEC,
#else
    kSKmUndefined,
#endif
    kSKmSigNuECCcoh, 
    kSKmSigNuECC1pi,
    //    kSKmSigNuECCoth_mu,
    kSKmSigNuECCoth,
    kSKmSigNuECCQE, 
    kSKmNuEBarNC1piC,
    kSKmNuEBarNCoth_mu,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCC1pi,
    //    kSKmNuEBarCCoth_mu,                                
    kSKmNuEBarCCoth,                                
    kSKmNuEBarCCQE, 
    kSKmNuENC1piC,
    kSKmNuENCoth_mu,
    kSKmNuECCcoh,
    kSKmNuECC1pi,
    //    kSKmNuECCoth_mu,                                
    kSKmNuECCoth,                                
    kSKmNuECCQE, 
    kSKmNuMuBarNC1piC, 
    kSKmNuMuBarNCoth_mu,
    kSKmNuMuBarCCcoh, 
    kSKmNuMuBarCC1pi, 
    //    kSKmNuMuBarCCoth_mu,                  
    kSKmNuMuBarCCoth,                  
    kSKmNuMuBarCCQE,      
    kSKmNuMuNC1piC,
    kSKmNuMuNCoth_mu,
    //    kSKmNuMuCCoth_mu,                      
    kSKmNuMuCCoth,                      
    kSKmNuMuCCcoh, 
    kSKmNuMuCC1pi, 
    kSKmNuMuCCQE,         
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
    kSKmNuMuBarCCQE, 
    kSKmNuMuBarCC1pi, 
    kSKmNuMuBarCCcoh, 
    //    kSKmNuMuBarCCoth_e,  
    kSKmNuMuBarCCoth,  
    kSKmNuMuBarNC1pi0, 
    kSKmNuMuBarNCcoh, 
    kSKmNuMuBarNCoth_e,
    kSKmNuEBarCCQE, 
    kSKmNuEBarCC1pi, 
    kSKmNuEBarCCcoh, 
    //    kSKmNuEBarCCoth_e,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1pi0,  
    kSKmNuEBarNCcoh, 
    kSKmNuEBarNCoth_e,
    kSKmNuECCQE, 
    kSKmNuECC1pi, 
    kSKmNuECCcoh, 
    //    kSKmNuECCoth_e,
    kSKmNuECCoth,
    kSKmNuENC1pi0,  
    kSKmNuENCcoh, 
    kSKmNuENCoth_e,
    kSKmSigNuECCQE,  
    kSKmSigNuECC1pi,  
    kSKmSigNuECCcoh,
    //    kSKmSigNuECCoth_e,
    kSKmSigNuECCoth,
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
    kSKmUndefined,
    kSKmUndefined
  }
  
};
  
} // conv
} // t2kana
  
#endif 

