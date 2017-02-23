//____________________________________________________________________________
/*!

\class    valor::conv

\brief    Conventions

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  Oct 20, 2010

Changes since 2011_v0.9.0:

17-12-2013  Nick
Split valor::conv namespace into several smaller files

*/
//____________________________________________________________________________

#ifndef formatting_h
#define formatting_h

#include "AnalysisOptions.h"
#include "SamplesAndModes.h"

namespace valor {
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
// The 5 listed attributes for each sample are fill style, fill colour, line style, line width, line colour.
//
static const int kSKHistoAttributes[kns][knm][kNGraphicAttributes] = {
  // > 1 muon-like ring
  {
    {   -1,  kWhite,    1,  1,  kBlack}, // kSKmNuMuCCQE
    { 3544,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCC1pi
    { 3244,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCCcoh
    { 3016,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCC2p2h
    { 3006,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCCoth
    { 3554,  kRed,      1,  1,  kBlack}, // kSKmNuMuNC1piC
    { 3006,  kRed,      1,  1,  kBlack}, // kSKmNuMuNC1pi0
    { 3006,  kRed,      1,  1,  kBlack}, // kSKmNuMuNCcoh
    { 3006,  kRed,      1,  1,  kBlack}, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    { 3001,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCQE
    { 3544,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC1pi
    { 3244,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCcoh
    { 3016,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC2p2h
    { 3006,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCoth
#endif
    { 3001,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCQE
    { 3345,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCC1pi
    { 3244,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCcoh
    { 3016,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCC2p2h
    { 3354,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCoth
    { 3544,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNC1piC
    { 3006,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNC1pi0
    { 3006,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNCcoh
    { 3006,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    { 3001,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCCQE
    { 3345,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCC1pi
    { 3244,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCCcoh
    { 3016,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCC2p2h
    { 3354,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCCoth
#endif
    { 3001,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCQE
    { 3345,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECC1pi
    { 3244,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCcoh
    { 3016,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECC2p2h
    { 3354,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCoth
    { 3544,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENC1piC
    { 3006,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENC1pi0
    { 3006,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENCcoh
    { 3006,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENCoth

    { 3001,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECCQE
    { 3345,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECC1pi
    { 3244,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECCcoh
    { 3016,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECC2p2h
    { 3354,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECCoth

    { 3001,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCQE
    { 3345,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCC1pi
    { 3244,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCcoh
    { 3016,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCC2p2h
    { 3354,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCoth
    { 3544,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNC1piC
    { 3006,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNC1pi0
    { 3006,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNCcoh
    { 3006,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNCoth

    { 3001,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCCQE
    { 3345,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCC1pi
    { 3244,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCCcoh
    { 3016,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCC2p2h
    { 3354,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
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
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
#endif

  }//1Rmu
  ,
  // > 1 electron-like ring
  {
    {   -1,    kWhite,  1,  1, kBlack }, // kSKmNuMuCCQE
    { 3544,    kBlack,  1,  1, kBlack }, // kSKmNuMuCC1pi
    { 3244,    kGreen,  1,  1, kBlack }, // kSKmNuMuCCcoh
    { 3016,    kGreen,  1,  1, kBlack }, // kSKmNuMuCC2p2h
    { 3244,     kBlue,  1,  1, kBlack }, // kSKmNuMuCCoth
    { 3554,      kRed,  1,  1, kBlack }, // kSKmNuMuNC1piC
    { 3006,      kRed,  1,  1, kBlack }, // kSKmNuMuNC1pi0
    { 3006,    kRed-1,  1,  1, kBlack }, // kSKmNuMuNCcoh
    { 3006,    kRed-2,  1,  1, kBlack }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    { 3001,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCQE
    { 3544,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC1pi
    { 3244,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCcoh
    { 3016,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC2p2h
    { 3006,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCoth
#endif
    { 3001,  kAzure-3,  1,  1, kBlack }, // kSKmNuMuBarCCQE
    { 3345,  kAzure-4,  1,  1, kBlack }, // kSKmNuMuBarCC1pi
    { 3354,  kAzure-6,  1,  1, kBlack }, // kSKmNuMuBarCCcoh
    { 3016,  kAzure-6,  1,  1, kBlack }, // kSKmNuMuBarCC2p2h
    { 3354,  kMagenta,  1,  1, kBlack }, // kSKmNuMuBarCCoth
    { 3554,  kAzure-7,  1,  1, kBlack }, // kSKmNuMuBarNC1piC
    { 3006,  kAzure-7,  1,  1, kBlack }, // kSKmNuMuBarNC1pi0
    { 3006,  kAzure-9,  1,  1, kBlack }, // kSKmNuMuBarNCcoh
    { 3006,  kMagenta,  1,  1, kBlack }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    { 3001,  kAzure-7,  1,  1, kBlack }, // kSKmOscNuMuBarCCQE
    { 3345,  kAzure-8,  1,  1, kBlack }, // kSKmOscNuMuBarCC1pi
    { 3354,  kAzure-9,  1,  1, kBlack }, // kSKmOscNuMuBarCCcoh
    { 3016,  kAzure-9,  1,  1, kBlack }, // kSKmOscNuMuBarCC2p2h
    { 3354,  kMagenta-3,  1,  1, kBlack }, // kSKmOscNuMuBarCCoth
#endif
    { 3001,  kGreen+1,  1,  1, kBlack }, // kSKmNuECCQE
    { 3345,  kGreen+2,  1,  1, kBlack }, // kSKmNuECC1pi
    { 3354,  kGreen+2,  1,  1, kBlack }, // kSKmNuECCcoh
    { 3016,  kGreen+2,  1,  1, kBlack }, // kSKmNuECC2p2h
    { 3354,    kBlack,  1,  1, kBlack }, // kSKmNuECCoth
    { 3554,  kGreen+5,  1,  1, kBlack }, // kSKmNuENC1piC
    { 3006,  kGreen+5,  1,  1, kBlack }, // kSKmNuENC1pi0
    { 3006,  kGreen+7,  1,  1, kBlack }, // kSKmNuENCcoh
    { 3006,    kBlack,  1,  1, kBlack }, // kSKmNuENCoth

    { 3001,  kGreen+8,  1,  1, kBlack }, // kSKmOscNuECCQE
    { 3345, kGreen+10,  1,  1, kBlack }, // kSKmOscNuECC1pi
    { 3354,      kRed,  1,  1, kBlack }, // kSKmOscNuECCcoh
    { 3016,      kRed,  1,  1, kBlack }, // kSKmOscNuECC2p2h
    { 3354,    kRed-2,  1,  1, kBlack }, // kSKmOscNuECCoth

    { 3001,  kBlue+1,   1,  1, kBlack }, // kSKmNuEBarCCQE
    { 3345,  kBlue+2,   1,  1, kBlack }, // kSKmNuEBarCC1pi
    { 3354,  kBlue+2,   1,  1, kBlack }, // kSKmNuEBarCCcoh
    { 3016,  kBlue+2,   1,  1, kBlack }, // kSKmNuEBarCC2p2h
    { 3354,    kBlack,  1,  1, kBlack }, // kSKmNuEBarCCoth
    { 3006,  kBlue+5,   1,  1, kBlack }, // kSKmNuEBarNC1piC
    { 3006,  kBlue+5,   1,  1, kBlack }, // kSKmNuEBarNC1pi0
    { 3006,  kBlue+7,   1,  1, kBlack }, // kSKmNuEBarNCcoh
    { 3006,    kBlack,  1,  1, kBlack }, // kSKmNuEBarNCoth

    { 3001,  kBlue+3,   1,  1, kBlack }, // kSKmOscNuEBarCCQE
    { 3345,  kBlue+5,   1,  1, kBlack }, // kSKmOscNuEBarCC1pi
    { 3354,  kBlue+5,   1,  1, kBlack }, // kSKmOscNuEBarCCcoh
    { 3016,  kBlue+5,   1,  1, kBlack }, // kSKmOscNuEBarCC2p2h
    { 3354,  kBlack-3,  1,  1, kBlack }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
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
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
#endif
  }//1Re
  ,
  // > 1 muon-like ring RHC
  {
    {   -1,  kWhite,    1,  1,  kBlack}, // kSKmNuMuCCQE
    { 3544,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCC1pi
    { 3244,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCCcoh
    { 3016,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCC2p2h
    { 3006,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCCoth
    { 3554,  kRed,      1,  1,  kBlack}, // kSKmNuMuNC1piC
    { 3006,  kRed,      1,  1,  kBlack}, // kSKmNuMuNC1pi0
    { 3006,  kRed,      1,  1,  kBlack}, // kSKmNuMuNCcoh
    { 3006,  kRed,      1,  1,  kBlack}, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    { 3001,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCQE
    { 3544,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC1pi
    { 3244,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCcoh
    { 3016,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC2p2h
    { 3006,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCoth
#endif
    { 3001,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCQE
    { 3345,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCC1pi
    { 3244,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCcoh
    { 3016,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCC2p2h
    { 3354,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCoth
    { 3544,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNC1piC
    { 3006,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNC1pi0
    { 3006,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNCcoh
    { 3006,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    { 3001,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCCQE
    { 3345,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCC1pi
    { 3244,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCCcoh
    { 3016,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCC2p2h
    { 3354,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCCoth
#endif
    { 3001,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCQE
    { 3345,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECC1pi
    { 3244,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCcoh
    { 3016,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECC2p2h
    { 3354,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCoth
    { 3544,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENC1piC
    { 3006,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENC1pi0
    { 3006,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENCcoh
    { 3006,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENCoth

    { 3001,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECCQE
    { 3345,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECC1pi
    { 3244,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECCcoh
    { 3016,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECC2p2h
    { 3354,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECCoth

    { 3001,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCQE
    { 3345,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCC1pi
    { 3244,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCcoh
    { 3016,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCC2p2h
    { 3354,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCoth
    { 3544,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNC1piC
    { 3006,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNC1pi0
    { 3006,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNCcoh
    { 3006,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNCoth

    { 3001,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCCQE
    { 3345,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCC1pi
    { 3244,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCCcoh
    { 3016,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCC2p2h
    { 3354,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
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
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
#endif
  }//1Rmu RHC
  ,
  // > 1 electron-like ring RHC
  {
    {   -1,    kWhite,  1,  1, kBlack }, // kSKmNuMuCCQE
    { 3544,    kBlack,  1,  1, kBlack }, // kSKmNuMuCC1pi
    { 3244,    kGreen,  1,  1, kBlack }, // kSKmNuMuCCcoh
    { 3016,    kGreen,  1,  1, kBlack }, // kSKmNuMuCC2p2h
    { 3244,     kBlue,  1,  1, kBlack }, // kSKmNuMuCCoth
    { 3554,      kRed,  1,  1, kBlack }, // kSKmNuMuNC1piC
    { 3006,      kRed,  1,  1, kBlack }, // kSKmNuMuNC1pi0
    { 3006,    kRed-1,  1,  1, kBlack }, // kSKmNuMuNCcoh
    { 3006,    kRed-2,  1,  1, kBlack }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    { 3001,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCQE
    { 3544,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC1pi
    { 3244,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCcoh
    { 3016,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC2p2h
    { 3006,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCoth
#endif
    { 3001,  kAzure-3,  1,  1, kBlack }, // kSKmNuMuBarCCQE
    { 3345,  kAzure-4,  1,  1, kBlack }, // kSKmNuMuBarCC1pi
    { 3354,  kAzure-6,  1,  1, kBlack }, // kSKmNuMuBarCCcoh
    { 3016,  kAzure-6,  1,  1, kBlack }, // kSKmNuMuBarCC2p2h
    { 3354,  kMagenta,  1,  1, kBlack }, // kSKmNuMuBarCCoth
    { 3554,  kAzure-7,  1,  1, kBlack }, // kSKmNuMuBarNC1piC
    { 3006,  kAzure-7,  1,  1, kBlack }, // kSKmNuMuBarNC1pi0
    { 3006,  kAzure-9,  1,  1, kBlack }, // kSKmNuMuBarNCcoh
    { 3006,  kMagenta,  1,  1, kBlack }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    { 3001,  kAzure-7,  1,  1, kBlack }, // kSKmOscNuMuBarCCQE
    { 3345,  kAzure-8,  1,  1, kBlack }, // kSKmOscNuMuBarCC1pi
    { 3354,  kAzure-9,  1,  1, kBlack }, // kSKmOscNuMuBarCCcoh
    { 3016,  kAzure-9,  1,  1, kBlack }, // kSKmOscNuMuBarCC2p2h
    { 3354,  kMagenta-3,  1,  1, kBlack }, // kSKmOscNuMuBarCCoth
#endif
    { 3001,  kGreen+1,  1,  1, kBlack }, // kSKmNuECCQE
    { 3345,  kGreen+2,  1,  1, kBlack }, // kSKmNuECC1pi
    { 3354,  kGreen+2,  1,  1, kBlack }, // kSKmNuECCcoh
    { 3016,  kGreen+2,  1,  1, kBlack }, // kSKmNuECC2p2h
    { 3354,    kBlack,  1,  1, kBlack }, // kSKmNuECCoth
    { 3554,  kGreen+5,  1,  1, kBlack }, // kSKmNuENC1piC
    { 3006,  kGreen+5,  1,  1, kBlack }, // kSKmNuENC1pi0
    { 3006,  kGreen+7,  1,  1, kBlack }, // kSKmNuENCcoh
    { 3006,    kBlack,  1,  1, kBlack }, // kSKmNuENCoth

    { 3001,  kGreen+8,  1,  1, kBlack }, // kSKmOscNuECCQE
    { 3345, kGreen+10,  1,  1, kBlack }, // kSKmOscNuECC1pi
    { 3354,      kRed,  1,  1, kBlack }, // kSKmOscNuECCcoh
    { 3016,      kRed,  1,  1, kBlack }, // kSKmOscNuECC2p2h
    { 3354,    kRed-2,  1,  1, kBlack }, // kSKmOscNuECCoth

    { 3001,  kBlue+1,   1,  1, kBlack }, // kSKmNuEBarCCQE
    { 3345,  kBlue+2,   1,  1, kBlack }, // kSKmNuEBarCC1pi
    { 3354,  kBlue+2,   1,  1, kBlack }, // kSKmNuEBarCCcoh
    { 3016,  kBlue+2,   1,  1, kBlack }, // kSKmNuEBarCC2p2h
    { 3354,    kBlack,  1,  1, kBlack }, // kSKmNuEBarCCoth
    { 3554,  kBlue+5,   1,  1, kBlack }, // kSKmNuEBarNC1piC
    { 3006,  kBlue+5,   1,  1, kBlack }, // kSKmNuEBarNC1pi0
    { 3006,  kBlue+7,   1,  1, kBlack }, // kSKmNuEBarNCcoh
    { 3006,    kBlack,  1,  1, kBlack }, // kSKmNuEBarNCoth

    { 3001,  kBlue+1,   1,  1, kBlack }, // kSKmOscNuEBarCCQE
    //{ 3345,  kBlue+5,   1,  1, kBlack }, // kSKmOscNuEBarCC1pi
    //{ 3354,  kBlue+5,   1,  1, kBlack }, // kSKmOscNuEBarCCcoh
    //{ 3016,  kBlue+5,   1,  1, kBlack }, // kSKmOscNuEBarCC2p2h
    { 3345,  kBlue-5,   1,  1, kBlack }, // kSKmOscNuEBarCC1pi
    { 3354,  kBlue-3,   1,  1, kBlack }, // kSKmOscNuEBarCCcoh
    { 3016,  kBlue+3,   1,  1, kBlack }, // kSKmOscNuEBarCC2p2h

    { 3354,  kBlack-3,  1,  1, kBlack }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
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
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
#endif
  }//1Re RHC
  ,
  // > Multi muon-like ring
  {
    {   -1,  kWhite,    1,  1,  kBlack}, // kSKmNuMuCCQE
    { 3544,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCC1pi
    { 3244,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCCcoh
    { 3016,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCC2p2h
    { 3006,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCCoth
    { 3554,  kRed,      1,  1,  kBlack}, // kSKmNuMuNC1piC
    { 3006,  kRed,      1,  1,  kBlack}, // kSKmNuMuNC1pi0
    { 3006,  kRed,      1,  1,  kBlack}, // kSKmNuMuNCcoh
    { 3006,  kRed,      1,  1,  kBlack}, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    { 3001,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCQE
    { 3544,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC1pi
    { 3244,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCcoh
    { 3016,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC2p2h
    { 3006,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCoth
#endif
    { 3001,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCQE
    { 3345,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCC1pi
    { 3244,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCcoh
    { 3016,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCC2p2h
    { 3354,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCoth
    { 3544,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNC1piC
    { 3006,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNC1pi0
    { 3006,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNCcoh
    { 3006,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    { 3001,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCCQE
    { 3345,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCC1pi
    { 3244,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCCcoh
    { 3016,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCC2p2h
    { 3354,  kAzure-5,  1,  1,  kBlack}, // kSKmOscNuMuBarCCoth
#endif
    { 3001,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCQE
    { 3345,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECC1pi
    { 3244,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCcoh
    { 3016,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECC2p2h
    { 3354,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCoth
    { 3544,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENC1piC
    { 3006,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENC1pi0
    { 3006,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENCcoh
    { 3006,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENCoth

    { 3001,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECCQE
    { 3345,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECC1pi
    { 3244,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECCcoh
    { 3016,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECC2p2h
    { 3354,  kGreen-3,  1,  1,  kBlack}, // kSKmOscNuECCoth

    { 3001,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCQE
    { 3345,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCC1pi
    { 3244,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCcoh
    { 3016,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCC2p2h
    { 3354,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCoth
    { 3544,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNC1piC
    { 3006,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNC1pi0
    { 3006,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNCcoh
    { 3006,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNCoth

    { 3001,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCCQE
    { 3345,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCC1pi
    { 3244,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCCcoh
    { 3016,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCC2p2h
    { 3354,  kBlue-3,     1,  1,  kBlack}, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
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
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
#endif

  }//MultiRmu
  ,
  // > Multi electron-like ring
  {
    {   -1,    kWhite,  1,  1, kBlack }, // kSKmNuMuCCQE
    { 3544,    kBlack,  1,  1, kBlack }, // kSKmNuMuCC1pi
    { 3244,    kGreen,  1,  1, kBlack }, // kSKmNuMuCCcoh
    { 3016,    kGreen,  1,  1, kBlack }, // kSKmNuMuCC2p2h
    { 3244,     kBlue,  1,  1, kBlack }, // kSKmNuMuCCoth
    { 3554,      kRed,  1,  1, kBlack }, // kSKmNuMuNC1piC
    { 3006,      kRed,  1,  1, kBlack }, // kSKmNuMuNC1pi0
    { 3006,    kRed-1,  1,  1, kBlack }, // kSKmNuMuNCcoh
    { 3006,    kRed-2,  1,  1, kBlack }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    { 3001,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCQE
    { 3544,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC1pi
    { 3244,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCcoh
    { 3016,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCC2p2h
    { 3006,  kGray,     1,  1,  kBlack}, // kSKmOscNuMuCCoth
#endif
    { 3001,  kAzure-3,  1,  1, kBlack }, // kSKmNuMuBarCCQE
    { 3345,  kAzure-4,  1,  1, kBlack }, // kSKmNuMuBarCC1pi
    { 3354,  kAzure-6,  1,  1, kBlack }, // kSKmNuMuBarCCcoh
    { 3016,  kAzure-6,  1,  1, kBlack }, // kSKmNuMuBarCC2p2h
    { 3354,  kMagenta,  1,  1, kBlack }, // kSKmNuMuBarCCoth
    { 3554,  kAzure-7,  1,  1, kBlack }, // kSKmNuMuBarNC1piC
    { 3006,  kAzure-7,  1,  1, kBlack }, // kSKmNuMuBarNC1pi0
    { 3006,  kAzure-9,  1,  1, kBlack }, // kSKmNuMuBarNCcoh
    { 3006,  kMagenta,  1,  1, kBlack }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    { 3001,  kAzure-7,  1,  1, kBlack }, // kSKmOscNuMuBarCCQE
    { 3345,  kAzure-8,  1,  1, kBlack }, // kSKmOscNuMuBarCC1pi
    { 3354,  kAzure-9,  1,  1, kBlack }, // kSKmOscNuMuBarCCcoh
    { 3016,  kAzure-9,  1,  1, kBlack }, // kSKmOscNuMuBarCC2p2h
    { 3354,  kMagenta-3,  1,  1, kBlack }, // kSKmOscNuMuBarCCoth
#endif
    { 3001,  kGreen+1,  1,  1, kBlack }, // kSKmNuECCQE
    { 3345,  kGreen+2,  1,  1, kBlack }, // kSKmNuECC1pi
    { 3354,  kGreen+2,  1,  1, kBlack }, // kSKmNuECCcoh
    { 3016,  kGreen+2,  1,  1, kBlack }, // kSKmNuECC2p2h
    { 3354,    kBlack,  1,  1, kBlack }, // kSKmNuECCoth
    { 3554,  kGreen+5,  1,  1, kBlack }, // kSKmNuENC1piC
    { 3006,  kGreen+5,  1,  1, kBlack }, // kSKmNuENC1pi0
    { 3006,  kGreen+7,  1,  1, kBlack }, // kSKmNuENCcoh
    { 3006,    kBlack,  1,  1, kBlack }, // kSKmNuENCoth

    { 3001,  kGreen+8,  1,  1, kBlack }, // kSKmOscNuECCQE
    { 3345, kGreen+10,  1,  1, kBlack }, // kSKmOscNuECC1pi
    { 3354,      kRed,  1,  1, kBlack }, // kSKmOscNuECCcoh
    { 3016,      kRed,  1,  1, kBlack }, // kSKmOscNuECC2p2h
    { 3354,    kRed-2,  1,  1, kBlack }, // kSKmOscNuECCoth

    { 3001,  kBlue+1,   1,  1, kBlack }, // kSKmNuEBarCCQE
    { 3345,  kBlue+2,   1,  1, kBlack }, // kSKmNuEBarCC1pi
    { 3354,  kBlue+2,   1,  1, kBlack }, // kSKmNuEBarCCcoh
    { 3016,  kBlue+2,   1,  1, kBlack }, // kSKmNuEBarCC2p2h
    { 3354,    kBlack,  1,  1, kBlack }, // kSKmNuEBarCCoth
    { 3554,  kBlue+5,   1,  1, kBlack }, // kSKmNuEBarNC1piC
    { 3006,  kBlue+5,   1,  1, kBlack }, // kSKmNuEBarNC1pi0
    { 3006,  kBlue+7,   1,  1, kBlack }, // kSKmNuEBarNCcoh
    { 3006,    kBlack,  1,  1, kBlack }, // kSKmNuEBarNCoth

    { 3001,  kBlue+3,   1,  1, kBlack }, // kSKmOscNuEBarCCQE
    { 3345,  kBlue+5,   1,  1, kBlack }, // kSKmOscNuEBarCC1pi
    { 3354,  kBlue+5,   1,  1, kBlack }, // kSKmOscNuEBarCCcoh
    { 3016,  kBlue+5,   1,  1, kBlack }, // kSKmOscNuEBarCC2p2h
    { 3354,  kBlack-3,  1,  1, kBlack }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
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
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
#endif
  }//MultiRe

};

//
// Define the order in which TH1D spectra are added in a THStack.
// The order is chosen to increase visibility of smaller components.
//

static const SuperkMode_t kSKReactionModesOrderInStack[kns][knm] = {

  // > 1 muon-like ring (numu disappearance)
  {
    kSKmOscNuECCcoh,
    kSKmOscNuECC2p2h,
    kSKmOscNuECC1pi,
    kSKmOscNuECCoth,
    kSKmOscNuECCQE,
    kSKmOscNuEBarCCcoh,
    kSKmOscNuEBarCC2p2h,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCoth,
    kSKmOscNuEBarCCQE,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCcoh,
    kSKmOscNuMuBarCC2p2h,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCoth,
    kSKmOscNuMuBarCCQE,
#endif
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmOscNuMuCCcoh,
    kSKmOscNuMuCC2p2h,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCoth,
    kSKmOscNuMuCCQE,
#endif
    kSKmNuEBarNC1piC,
    kSKmNuEBarNC1pi0,
    kSKmNuEBarNCcoh,
    kSKmNuEBarNCoth,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCC2p2h,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCoth,
    kSKmNuEBarCCQE,
    kSKmNuENC1piC,
    kSKmNuENC1pi0,
    kSKmNuENCcoh,
    kSKmNuENCoth,
    kSKmNuECCcoh,
    kSKmNuECC2p2h,
    kSKmNuECC1pi,
    kSKmNuECCoth,
    kSKmNuECCQE,
    kSKmNuMuBarNC1piC,
    kSKmNuMuBarNC1pi0,
    kSKmNuMuBarNCcoh,
    kSKmNuMuBarNCoth,
    kSKmNuMuBarCCcoh,
    kSKmNuMuBarCC2p2h,
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCoth,
    kSKmNuMuBarCCQE,
    kSKmNuMuNC1piC,
    kSKmNuMuNC1pi0,
    kSKmNuMuNCcoh,
    kSKmNuMuNCoth,
    kSKmNuMuCCoth,
    kSKmNuMuCCcoh,
    kSKmNuMuCC2p2h,
    kSKmNuMuCC1pi,
    kSKmNuMuCCQE,
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
  }//1Rmu
  ,
  // > 1 electron-like ring (nue appearance)
  {
    kSKmOscNuEBarCCcoh,
    kSKmOscNuEBarCC2p2h,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCoth,
    kSKmOscNuEBarCCQE,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCcoh,
    kSKmOscNuMuBarCC2p2h,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCoth,
    kSKmOscNuMuBarCCQE,
#endif
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmOscNuMuCCcoh,
    kSKmOscNuMuCC2p2h,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCoth,
    kSKmOscNuMuCCQE,
#endif
    kSKmNuMuCCQE,
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    kSKmNuMuCC2p2h,
    kSKmNuMuCCoth,
    kSKmNuMuNC1piC,
    kSKmNuMuNC1pi0,
    kSKmNuMuNCcoh,
    kSKmNuMuNCoth,
    kSKmNuMuBarCCQE,
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCcoh,
    kSKmNuMuBarCC2p2h,
    kSKmNuMuBarCCoth,
    kSKmNuMuBarNC1piC,
    kSKmNuMuBarNC1pi0,
    kSKmNuMuBarNCcoh,
    kSKmNuMuBarNCoth,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCC2p2h,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1piC,
    kSKmNuEBarNC1pi0,
    kSKmNuEBarNCcoh,
    kSKmNuEBarNCoth,
    kSKmNuECCQE,
    kSKmNuECC1pi,
    kSKmNuECCcoh,
    kSKmNuECC2p2h,
    kSKmNuECCoth,
    kSKmNuENC1piC,
    kSKmNuENC1pi0,
    kSKmNuENCcoh,
    kSKmNuENCoth,
    kSKmOscNuECCQE,
    kSKmOscNuECC1pi,
    kSKmOscNuECCcoh,
    kSKmOscNuECC2p2h,
    kSKmOscNuECCoth,
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
  }//1Re
  ,
  // > 1 muon-like ring (numu disappearance) RHC
  {
    kSKmOscNuECCcoh,
    kSKmOscNuECC2p2h,
    kSKmOscNuECC1pi,
    kSKmOscNuECCoth,
    kSKmOscNuECCQE,
    kSKmOscNuEBarCCcoh,
    kSKmOscNuEBarCC2p2h,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCoth,
    kSKmOscNuEBarCCQE,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCcoh,
    kSKmOscNuMuBarCC2p2h,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCoth,
    kSKmOscNuMuBarCCQE,
#endif
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmOscNuMuCCcoh,
    kSKmOscNuMuCC2p2h,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCoth,
    kSKmOscNuMuCCQE,
#endif
    kSKmNuEBarNC1piC,
    kSKmNuEBarNC1pi0,
    kSKmNuEBarNCcoh,
    kSKmNuEBarNCoth,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCC2p2h,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCoth,
    kSKmNuEBarCCQE,
    kSKmNuENC1piC,
    kSKmNuENC1pi0,
    kSKmNuENCcoh,
    kSKmNuENCoth,
    kSKmNuECCcoh,
    kSKmNuECC2p2h,
    kSKmNuECC1pi,
    kSKmNuECCoth,
    kSKmNuECCQE,
    kSKmNuMuBarNC1piC,
    kSKmNuMuBarNC1pi0,
    kSKmNuMuBarNCcoh,
    kSKmNuMuBarNCoth,
    kSKmNuMuBarCCcoh,
    kSKmNuMuBarCC2p2h,
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCoth,
    kSKmNuMuBarCCQE,
    kSKmNuMuNC1piC,
    kSKmNuMuNC1pi0,
    kSKmNuMuNCcoh,
    kSKmNuMuNCoth,
    kSKmNuMuCCoth,
    kSKmNuMuCCcoh,
    kSKmNuMuCC2p2h,
    kSKmNuMuCC1pi,
    kSKmNuMuCCQE,
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
  }//1Rmu RHC
  ,
  // > 1 electron-like ring (nue appearance) RHC
  {
    kSKmOscNuEBarCCcoh,
    kSKmOscNuEBarCC2p2h,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCoth,
    kSKmOscNuEBarCCQE,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCcoh,
    kSKmOscNuMuBarCC2p2h,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCoth,
    kSKmOscNuMuBarCCQE,
#endif
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmOscNuMuCCcoh,
    kSKmOscNuMuCC2p2h,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCoth,
    kSKmOscNuMuCCQE,
#endif
    kSKmNuMuCCQE,
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    kSKmNuMuCC2p2h,
    kSKmNuMuCCoth,
    kSKmNuMuNC1piC,
    kSKmNuMuNC1pi0,
    kSKmNuMuNCcoh,
    kSKmNuMuNCoth,
    kSKmNuMuBarCCQE,
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCcoh,
    kSKmNuMuBarCC2p2h,
    kSKmNuMuBarCCoth,
    kSKmNuMuBarNC1piC,
    kSKmNuMuBarNC1pi0,
    kSKmNuMuBarNCcoh,
    kSKmNuMuBarNCoth,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCC2p2h,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1piC,
    kSKmNuEBarNC1pi0,
    kSKmNuEBarNCcoh,
    kSKmNuEBarNCoth,
    kSKmNuECCQE,
    kSKmNuECC1pi,
    kSKmNuECCcoh,
    kSKmNuECC2p2h,
    kSKmNuECCoth,
    kSKmNuENC1piC,
    kSKmNuENC1pi0,
    kSKmNuENCcoh,
    kSKmNuENCoth,
    kSKmOscNuECCQE,
    kSKmOscNuECC1pi,
    kSKmOscNuECCcoh,
    kSKmOscNuECC2p2h,
    kSKmOscNuECCoth,
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
  }//1Re RHC
  ,
  // > Multi muon-like ring (numu disappearance)
  {
    kSKmOscNuECCcoh,
    kSKmOscNuECC2p2h,
    kSKmOscNuECC1pi,
    kSKmOscNuECCoth,
    kSKmOscNuECCQE,
    kSKmOscNuEBarCCcoh,
    kSKmOscNuEBarCC2p2h,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCoth,
    kSKmOscNuEBarCCQE,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCcoh,
    kSKmOscNuMuBarCC2p2h,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCoth,
    kSKmOscNuMuBarCCQE,
#endif
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmOscNuMuCCcoh,
    kSKmOscNuMuCC2p2h,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCoth,
    kSKmOscNuMuCCQE,
#endif
    kSKmNuEBarNC1piC,
    kSKmNuEBarNC1pi0,
    kSKmNuEBarNCcoh,
    kSKmNuEBarNCoth,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCC2p2h,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCoth,
    kSKmNuEBarCCQE,
    kSKmNuENC1piC,
    kSKmNuENC1pi0,
    kSKmNuENCcoh,
    kSKmNuENCoth,
    kSKmNuECCcoh,
    kSKmNuECC2p2h,
    kSKmNuECC1pi,
    kSKmNuECCoth,
    kSKmNuECCQE,
    kSKmNuMuBarNC1piC,
    kSKmNuMuBarNC1pi0,
    kSKmNuMuBarNCcoh,
    kSKmNuMuBarNCoth,
    kSKmNuMuBarCCcoh,
    kSKmNuMuBarCC2p2h,
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCoth,
    kSKmNuMuBarCCQE,
    kSKmNuMuNC1piC,
    kSKmNuMuNC1pi0,
    kSKmNuMuNCcoh,
    kSKmNuMuNCoth,
    kSKmNuMuCCoth,
    kSKmNuMuCCcoh,
    kSKmNuMuCC2p2h,
    kSKmNuMuCC1pi,
    kSKmNuMuCCQE,
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
  }//MultiRmu
  ,
  // > Multi electron-like ring (nue appearance)
  {
    kSKmOscNuEBarCCcoh,
    kSKmOscNuEBarCC2p2h,
    kSKmOscNuEBarCC1pi,
    kSKmOscNuEBarCCoth,
    kSKmOscNuEBarCCQE,
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmOscNuMuBarCCcoh,
    kSKmOscNuMuBarCC2p2h,
    kSKmOscNuMuBarCC1pi,
    kSKmOscNuMuBarCCoth,
    kSKmOscNuMuBarCCQE,
#endif
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmOscNuMuCCcoh,
    kSKmOscNuMuCC2p2h,
    kSKmOscNuMuCC1pi,
    kSKmOscNuMuCCoth,
    kSKmOscNuMuCCQE,
#endif
    kSKmNuMuCCQE,
    kSKmNuMuCC1pi,
    kSKmNuMuCCcoh,
    kSKmNuMuCC2p2h,
    kSKmNuMuCCoth,
    kSKmNuMuNC1piC,
    kSKmNuMuNC1pi0,
    kSKmNuMuNCcoh,
    kSKmNuMuNCoth,
    kSKmNuMuBarCCQE,
    kSKmNuMuBarCC1pi,
    kSKmNuMuBarCCcoh,
    kSKmNuMuBarCC2p2h,
    kSKmNuMuBarCCoth,
    kSKmNuMuBarNC1piC,
    kSKmNuMuBarNC1pi0,
    kSKmNuMuBarNCcoh,
    kSKmNuMuBarNCoth,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCC2p2h,
    kSKmNuEBarCCoth,
    kSKmNuEBarNC1piC,
    kSKmNuEBarNC1pi0,
    kSKmNuEBarNCcoh,
    kSKmNuEBarNCoth,
    kSKmNuECCQE,
    kSKmNuECC1pi,
    kSKmNuECCcoh,
    kSKmNuECC2p2h,
    kSKmNuECCoth,
    kSKmNuENC1piC,
    kSKmNuENC1pi0,
    kSKmNuENCcoh,
    kSKmNuENCoth,
    kSKmOscNuECCQE,
    kSKmOscNuECC1pi,
    kSKmOscNuECCcoh,
    kSKmOscNuECC2p2h,
    kSKmOscNuECCoth,
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
    kSKmUndefined,
#endif
  }//MultiRe


};

} // conv
} // valor

#endif

