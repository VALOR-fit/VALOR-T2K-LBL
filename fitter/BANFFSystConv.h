//____________________________________________________________________________
/*!

\class    t2kana::conv

\brief    Conventions

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

Changes since 2011_v0.9.0:

17-12-2013  Nick  
Split t2kana::conv namespace into several smaller files

22-12-2013  Nick
Added extra arrays for reverse horn current systematic parameters
*/
//____________________________________________________________________________

#ifndef banff_syst_conv_h
#define banff_syst_conv_h

#include "AnalysisOptions.h"
#include "Enumerations.h"
#include "SamplesAndModes.h"

namespace t2kana {
namespace conv   {

//
// Number of BANFF fit parameters
//
#if BANFF_BINS == 1 //joint
static const int kNumOfParamsInBANFFfit = 33;
static const int kNumOfFluxParamsInBANFFfit = kNumOfParamsInBANFFfit - 8;
static const int kNumOfRHCParamsInBANFFfit = 33;
static const int kNumOfRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - 8;
#elif BANFF_BINS == 2 //disap
static const int kNumOfParamsInBANFFfit = 23;
static const int kNumOfFluxParamsInBANFFfit = kNumOfParamsInBANFFfit - 7;
static const int kNumOfRHCParamsInBANFFfit = 23;
static const int kNumOfRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - 7;
#elif BANFF_BINS == 3 //app
static const int kNumOfParamsInBANFFfit = 27;
static const int kNumOfFluxParamsInBANFFfit = kNumOfParamsInBANFFfit - 5;
static const int kNumOfRHCParamsInBANFFfit = 27;
static const int kNumOfRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - 5;
#elif BANFF_BINS == 4 //disap + MEC
static const int kNumOfParamsInBANFFfit = 25;
static const int kNumOfFluxParamsInBANFFfit = kNumOfParamsInBANFFfit - 7;
static const int kNumOfRHCParamsInBANFFfit = 25;
static const int kNumOfRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - 7;
#endif

//
// ** Flux systematics
//

// Mark which BANFF parameters are SuperK flux params
static const bool kParamInBANFFfitIsSuperkFluxParam[kNumOfParamsInBANFFfit] =
{
   true,                 // kNuisBANFF00
   true,                 // kNuisBANFF01
   true,                 // kNuisBANFF02
   true,                 // kNuisBANFF03
   true,                 // kNuisBANFF04
   true,                 // kNuisBANFF05
   true,                 // kNuisBANFF06
   true,                 // kNuisBANFF07
   true,                 // kNuisBANFF08
   true,                 // kNuisBANFF09
   true,                 // kNuisBANFF10
   true,                 // kNuisBANFF11
   true,                 // kNuisBANFF12
   true,                 // kNuisBANFF13
   true,                 // kNuisBANFF14
   true,                 // kNuisBANFF15
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
   true,                 // kNuisBANFF16
   true,                 // kNuisBANFF17
   true,                 // kNuisBANFF18
   true,                 // kNuisBANFF19
   true,                 // kNuisBANFF20
   true,                 // kNuisBANFF21
#endif
#if BANFF_BINS == 1
   true,                 // kNuisBANFF22
   true,                 // kNuisBANFF23
   true,                 // kNuisBANFF24
#endif
   false,                // kNuisBANFFMAQE
   false,                // kNuisBANFFMARES
   false,                // kNuisBANFFCCQEE0
#if BANFF_BINS != 3
   false,                // kNuisBANFFCCQEE1
   false,                // kNuisBANFFCCQEE2
#endif
   false,                // kNuisBANFFCC1piE0
#if BANFF_BINS != 3
   false,                // kNuisBANFFCC1piE1
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
   false                 // kNuisBANFFNC1pi0
#endif
#if BANFF_BINS == 4
   false,                // kNuisBANFFMECNorm
   false                 // kNuisBANFFPDDDoubleCount
#endif
};

static const bool kRHCParamInBANFFfitIsSuperkFluxParam[kNumOfRHCParamsInBANFFfit] =
  {
    true,                 // kRHCNuisBANFF00
    true,                 // kRHCNuisBANFF01
    true,                 // kRHCNuisBANFF02
    true,                 // kRHCNuisBANFF03
    true,                 // kRHCNuisBANFF04
    true,                 // kRHCNuisBANFF05
    true,                 // kRHCNuisBANFF06
    true,                 // kRHCNuisBANFF07
    true,                 // kRHCNuisBANFF08
    true,                 // kRHCNuisBANFF09
    true,                 // kRHCNuisBANFF10
    true,                 // kRHCNuisBANFF11
    true,                 // kRHCNuisBANFF12
    true,                 // kRHCNuisBANFF13
    true,                 // kRHCNuisBANFF14
    true,                 // kRHCNuisBANFF15
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    true,                 // kRHCNuisBANFF16
    true,                 // kRHCNuisBANFF17
    true,                 // kRHCNuisBANFF18
    true,                 // kRHCNuisBANFF19
    true,                 // kRHCNuisBANFF20
    true,                 // kRHCNuisBANFF21
#endif
#if BANFF_BINS == 1
    true,                 // kRHCNuisBANFF22
    true,                 // kRHCNuisBANFF23
    true,                 // kRHCNuisBANFF24
#endif
    false,                // kRHCNuisBANFFMAQE
    false,                // kRHCNuisBANFFMARES
    false,                // kRHCNuisBANFFCCQEE0
#if BANFF_BINS != 3
    false,                // kRHCNuisBANFFCCQEE1
    false,                // kRHCNuisBANFFCCQEE2
#endif
    false,                // kRHCNuisBANFFCC1piE0
#if BANFF_BINS != 3
    false,                // kRHCNuisBANFFCC1piE1
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    false                 // kRHCNuisBANFFNC1pi0
#endif
#if BANFF_BINS == 4
    false,                // kRHCNuisBANFFMECNorm
    false                 // kRHCNuisBANFFPDDDoubleCount
#endif
  };



// flux-related BANFF parameters --> flux neutrino PDG code
static const int kNuPdgCodeOfFluxParamInBANFFfit[kNumOfParamsInBANFFfit] =
{
#if (BANFF_BINS == 2) || (BANFF_BINS == 4)
   kPdgNuMu,             // kNuisBANFF00
   kPdgNuMu,             // kNuisBANFF01
   kPdgNuMu,             // kNuisBANFF02
   kPdgNuMu,             // kNuisBANFF03
   kPdgNuMu,             // kNuisBANFF04
   kPdgNuMu,             // kNuisBANFF05
   kPdgNuMu,             // kNuisBANFF06
   kPdgNuMu,             // kNuisBANFF07
   kPdgNuMu,             // kNuisBANFF08
   kPdgNuMu,             // kNuisBANFF09
   kPdgNuMu,             // kNuisBANFF10
   kPdgNuMuBar,          // kNuisBANFF11
   kPdgNuMuBar,          // kNuisBANFF12
   kPdgNuMuBar,          // kNuisBANFF13
   kPdgNuMuBar,          // kNuisBANFF14
   kPdgNuMuBar,          // kNuisBANFF15
#endif
#if BANFF_BINS == 3
   kPdgNuMu,             // kNuisBANFF00        
   kPdgNuMu,             // kNuisBANFF01                                                                                                                                          
   kPdgNuMu,             // kNuisBANFF02                                                                                                                                            
   kPdgNuMu,             // kNuisBANFF03                                                                                                                                           
   kPdgNuMu,             // kNuisBANFF04                                                                                                                                           
   kPdgNuMu,             // kNuisBANFF05                                                                                                                                          
   kPdgNuMu,             // kNuisBANFF06                                                                                                                                           
   kPdgNuMu,             // kNuisBANFF07                                                                                                                                          
   kPdgNuMu,             // kNuisBANFF08                                                                                                                                          
   kPdgNuMu,             // kNuisBANFF09                                                                                                                                          
   kPdgNuMu,             // kNuisBANFF10                                                                                                                                          
   kPdgNuMuBar,          // kNuisBANFF11                                                                                                                                           
   kPdgNuMuBar,          // kNuisBANFF12                                                                                                                                          
   kPdgNuE,              // kNuisBANFF13
   kPdgNuE,              // kNuisBANFF14
   kPdgNuE,              // kNuisBANFF15
   kPdgNuE,              // kNuisBANFF16
   kPdgNuE,              // kNuisBANFF17
   kPdgNuE,              // kNuisBANFF18
   kPdgNuE,              // kNuisBANFF19
   kPdgNuEBar,           // kNuisBANFF20
   kPdgNuEBar,           // kNuisBANFF21
#endif
#if BANFF_BINS == 1
   kPdgNuMu,             // kNuisBANFF00
   kPdgNuMu,             // kNuisBANFF01
   kPdgNuMu,             // kNuisBANFF02
   kPdgNuMu,             // kNuisBANFF03
   kPdgNuMu,             // kNuisBANFF04
   kPdgNuMu,             // kNuisBANFF05
   kPdgNuMu,             // kNuisBANFF06
   kPdgNuMu,             // kNuisBANFF07
   kPdgNuMu,             // kNuisBANFF08
   kPdgNuMu,             // kNuisBANFF09
   kPdgNuMu,             // kNuisBANFF10
   kPdgNuMuBar,          // kNuisBANFF11
   kPdgNuMuBar,          // kNuisBANFF12
   kPdgNuMuBar,          // kNuisBANFF13
   kPdgNuMuBar,          // kNuisBANFF14
   kPdgNuMuBar,          // kNuisBANFF15
   kPdgNuE,              // kNuisBANFF16
   kPdgNuE,              // kNuisBANFF17
   kPdgNuE,              // kNuisBANFF18
   kPdgNuE,              // kNuisBANFF19
   kPdgNuE,              // kNuisBANFF20
   kPdgNuE,              // kNuisBANFF21
   kPdgNuE,              // kNuisBANFF22
   kPdgNuEBar,           // kNuisBANFF23
   kPdgNuEBar,           // kNuisBANFF24
#endif
   0,                    // kNuisBANFFMAQE, not a flux param
   0,                    // kNuisBANFFMARES, not a flux param
   0,                    // kNuisBANFFCCQEE0, not a flux param
#if BANFF_BINS != 3
   0,                    // kNuisBANFFCCQEE1, not a flux param
   0,                    // kNuisBANFFCCQEE2, not a flux param
#endif
   0,                    // kNuisBANFFCC1piE0, not a flux param
#if BANFF_BINS != 3
   0,                    // kNuisBANFFCC1piE1, not a flux param
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
   0                     // kNuisBANFFNC1pi0, not a flux param
#endif
#if BANFF_BINS == 4
   0,                    // kNuisBANFFMECNorm
   0                     // kNuisBANFFPDDDoubleCount
#endif
};

static const int kNuPdgCodeOfRHCFluxParamInBANFFfit[kNumOfRHCParamsInBANFFfit] =
  {
#if (BANFF_BINS == 2) || (BANFF_BINS == 4)
    kPdgNuMu,             // kRHCNuisBANFF00
    kPdgNuMu,             // kRHCNuisBANFF01
    kPdgNuMu,             // kRHCNuisBANFF02
    kPdgNuMu,             // kRHCNuisBANFF03
    kPdgNuMu,             // kRHCNuisBANFF04
    kPdgNuMu,             // kRHCNuisBANFF05
    kPdgNuMu,             // kRHCNuisBANFF06
    kPdgNuMu,             // kRHCNuisBANFF07
    kPdgNuMu,             // kRHCNuisBANFF08
    kPdgNuMu,             // kRHCNuisBANFF09
    kPdgNuMu,             // kRHCNuisBANFF10
    kPdgNuMuBar,          // kRHCNuisBANFF11
    kPdgNuMuBar,          // kRHCNuisBANFF12
    kPdgNuMuBar,          // kRHCNuisBANFF13
    kPdgNuMuBar,          // kRHCNuisBANFF14
    kPdgNuMuBar,          // kRHCNuisBANFF15
#endif
#if BANFF_BINS == 3
    kPdgNuMu,             // kRHCNuisBANFF00        
    kPdgNuMu,             // kRHCNuisBANFF01                                                                                                                                          
    kPdgNuMu,             // kRHCNuisBANFF02                                                                                                                                            
    kPdgNuMu,             // kRHCNuisBANFF03                                                                                                                                           
    kPdgNuMu,             // kRHCNuisBANFF04                                                                                                                                           
    kPdgNuMu,             // kRHCNuisBANFF05                                                                                                                                          
    kPdgNuMu,             // kRHCNuisBANFF06                                                                                                                                           
    kPdgNuMu,             // kRHCNuisBANFF07                                                                                                                                          
    kPdgNuMu,             // kRHCNuisBANFF08                                                                                                                                          
    kPdgNuMu,             // kRHCNuisBANFF09                                                                                                                                          
    kPdgNuMu,             // kRHCNuisBANFF10                                                                                                                                          
    kPdgNuMuBar,          // kRHCNuisBANFF11                                                                                                                                           
    kPdgNuMuBar,          // kRHCNuisBANFF12                                                                                                                                          
    kPdgNuE,              // kRHCNuisBANFF13
    kPdgNuE,              // kRHCNuisBANFF14
    kPdgNuE,              // kRHCNuisBANFF15
    kPdgNuE,              // kRHCNuisBANFF16
    kPdgNuE,              // kRHCNuisBANFF17
    kPdgNuE,              // kRHCNuisBANFF18
    kPdgNuE,              // kRHCNuisBANFF19
    kPdgNuEBar,           // kRHCNuisBANFF20
    kPdgNuEBar,           // kRHCNuisBANFF21
#endif
#if BANFF_BINS == 1
    kPdgNuMu,             // kRHCNuisBANFF00
    kPdgNuMu,             // kRHCNuisBANFF01
    kPdgNuMu,             // kRHCNuisBANFF02
    kPdgNuMu,             // kRHCNuisBANFF03
    kPdgNuMu,             // kRHCNuisBANFF04
    kPdgNuMu,             // kRHCNuisBANFF05
    kPdgNuMu,             // kRHCNuisBANFF06
    kPdgNuMu,             // kRHCNuisBANFF07
    kPdgNuMu,             // kRHCNuisBANFF08
    kPdgNuMu,             // kRHCNuisBANFF09
    kPdgNuMu,             // kRHCNuisBANFF10
    kPdgNuMuBar,          // kRHCNuisBANFF11
    kPdgNuMuBar,          // kRHCNuisBANFF12
    kPdgNuMuBar,          // kRHCNuisBANFF13
    kPdgNuMuBar,          // kRHCNuisBANFF14
    kPdgNuMuBar,          // kRHCNuisBANFF15
    kPdgNuE,              // kRHCNuisBANFF16
    kPdgNuE,              // kRHCNuisBANFF17
    kPdgNuE,              // kRHCNuisBANFF18
    kPdgNuE,              // kRHCNuisBANFF19
    kPdgNuE,              // kRHCNuisBANFF20
    kPdgNuE,              // kRHCNuisBANFF21
    kPdgNuE,              // kRHCNuisBANFF22
    kPdgNuEBar,           // kRHCNuisBANFF23
    kPdgNuEBar,           // kRHCNuisBANFF24
#endif
    0,                    // kRHCNuisBANFFMAQE, not a flux param
    0,                    // kRHCNuisBANFFMARES, not a flux param
    0,                    // kRHCNuisBANFFCCQEE0, not a flux param
#if BANFF_BINS != 3
    0,                    // kRHCNuisBANFFCCQEE1, not a flux param
    0,                    // kRHCNuisBANFFCCQEE2, not a flux param
#endif
    0,                    // kRHCNuisBANFFCC1piE0, not a flux param
#if BANFF_BINS != 3
    0,                    // kRHCNuisBANFFCC1piE1, not a flux param
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    0                     // kRHCNuisBANFFNC1pi0, not a flux param
#endif
#if BANFF_BINS == 4
    0,                    // kRHCNuisBANFFMECNorm
    0                     // kRHCNuisBANFFPDDDoubleCount
#endif
  };


// Etrue energy range for flux-related BANFF parameters
static const double kEtrueRangeOfFluxParamInBANFFfit[kNumOfParamsInBANFFfit][2] =
{
#if BANFF_BINS == 2
   { 0.0,  0.4 },        // kNuisBANFF00
   { 0.4,  0.5 },        // kNuisBANFF01
   { 0.5,  0.6 },        // kNuisBANFF02
   { 0.6,  0.7 },        // kNuisBANFF03
   { 0.7,  1.0 },        // kNuisBANFF04
   { 1.0,  1.5 },        // kNuisBANFF05
   { 1.5,  2.5 },        // kNuisBANFF06
   { 2.5,  3.5 },        // kNuisBANFF07
   { 3.5,  5.0 },        // kNuisBANFF08
   { 5.0,  7.0 },        // kNuisBANFF09
   { 7.0, 30.0 },        // kNuisBANFF10
   { 0.0,  0.7 },        // kNuisBANFF11
   { 0.7,  1.0 },        // kNuisBANFF12
   { 1.0,  1.5 },        // kNuisBANFF13
   { 1.5,  2.5 },        // kNuisBANFF14
   { 2.5, 30.0 },        // kNuisBANFF15
#endif
#if BANFF_BINS == 3
   { 0.0,  0.4 },        // kNuisBANFF00
   { 0.4,  0.5 },        // kNuisBANFF01
   { 0.5,  0.6 },        // kNuisBANFF02
   { 0.6,  0.7 },        // kNuisBANFF03
   { 0.7,  1.0 },        // kNuisBANFF04
   { 1.0,  1.5 },        // kNuisBANFF05
   { 1.5,  2.5 },        // kNuisBANFF06
   { 2.5,  3.5 },        // kNuisBANFF07
   { 3.5,  5.0 },        // kNuisBANFF08
   { 5.0,  7.0 },        // kNuisBANFF09
   { 7.0, 30.0 },        // kNuisBANFF10
   { 0.0,  1.5 },        // kNuisBANFF11
   { 1.5, 30.0 },        // kNuisBANFF12
   { 0.0,  0.5 },        // kNuisBANFF13
   { 0.5,  0.7 },        // kNuisBANFF14
   { 0.7,  0.8 },        // kNuisBANFF15
   { 0.8,  1.5 },        // kNuisBANFF16
   { 1.5,  2.5 },        // kNuisBANFF17
   { 2.5,  4.0 },        // kNuisBANFF18
   { 4.0, 30.0 },        // kNuisBANFF19
   { 0.0,  2.5 },        // kNuisBANFF20
   { 2.5, 30.0 },        // kNuisBANFF21
#endif
#if BANFF_BINS == 1
   { 0.0,  0.4 },        // kNuisBANFF00
   { 0.4,  0.5 },        // kNuisBANFF01
   { 0.5,  0.6 },        // kNuisBANFF02
   { 0.6,  0.7 },        // kNuisBANFF03
   { 0.7,  1.0 },        // kNuisBANFF04
   { 1.0,  1.5 },        // kNuisBANFF05
   { 1.5,  2.5 },        // kNuisBANFF06
   { 2.5,  3.5 },        // kNuisBANFF07
   { 3.5,  5.0 },        // kNuisBANFF08
   { 5.0,  7.0 },        // kNuisBANFF09
   { 7.0, 30.0 },        // kNuisBANFF10
   { 0.0,  0.7 },        // kNuisBANFF11
   { 0.7,  1.0 },        // kNuisBANFF12
   { 1.0,  1.5 },        // kNuisBANFF13
   { 1.5,  2.5 },        // kNuisBANFF14
   { 2.5, 30.0 },        // kNuisBANFF15
   { 0.0,  0.5 },        // kNuisBANFF16
   { 0.5,  0.7 },        // kNuisBANFF17
   { 0.7,  0.8 },        // kNuisBANFF18
   { 0.8,  1.5 },        // kNuisBANFF19
   { 1.5,  2.5 },        // kNuisBANFF20
   { 2.5,  4.0 },        // kNuisBANFF21
   { 4.0, 30.0 },        // kNuisBANFF22
   { 0.0,  2.5 },        // kNuisBANFF23
   { 2.5, 30.0 },        // kNuisBANFF24
#endif
   {  -1,   -1 },        // kNuisBANFFMAQE, not a flux param
   {  -1,   -1 },        // kNuisBANFFMARES, not a flux param
   {  -1,   -1 },        // kNuisBANFFCCQEE0, not a flux param
#if BANFF_BINS != 3
   {  -1,   -1 },        // kNuisBANFFCCQEE1, not a flux param
   {  -1,   -1 },        // kNuisBANFFCCQEE2, not a flux param
#endif
   {  -1,   -1 },        // kNuisBANFFCC1piE0, not a flux param
#if BANFF_BINS != 3
   {  -1,   -1 },        // kNuisBANFFCC1piE1, not a flux param
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
   {  -1,   -1 }         // kNuisBANFFNC1pi0, not a flux param
#endif
#if BANFF_BINS == 4
   {  -1,   -1 },        // kNuisBANFFMECNorm
   {  -1,   -1 }         // kNuisBANFFPDDDoubleCount
#endif
};

static const double kEtrueRangeOfRHCFluxParamInBANFFfit[kNumOfRHCParamsInBANFFfit][2] =
  {
#if BANFF_BINS == 2
    { 0.0,  0.4 },        // kRHCNuisBANFF00
    { 0.4,  0.5 },        // kRHCNuisBANFF01
    { 0.5,  0.6 },        // kRHCNuisBANFF02
    { 0.6,  0.7 },        // kRHCNuisBANFF03
    { 0.7,  1.0 },        // kRHCNuisBANFF04
    { 1.0,  1.5 },        // kRHCNuisBANFF05
    { 1.5,  2.5 },        // kRHCNuisBANFF06
    { 2.5,  3.5 },        // kRHCNuisBANFF07
    { 3.5,  5.0 },        // kRHCNuisBANFF08
    { 5.0,  7.0 },        // kRHCNuisBANFF09
    { 7.0, 30.0 },        // kRHCNuisBANFF10
    { 0.0,  0.7 },        // kRHCNuisBANFF11
    { 0.7,  1.0 },        // kRHCNuisBANFF12
    { 1.0,  1.5 },        // kRHCNuisBANFF13
    { 1.5,  2.5 },        // kRHCNuisBANFF14
    { 2.5, 30.0 },        // kRHCNuisBANFF15
#endif
#if BANFF_BINS == 3
    { 0.0,  0.4 },        // kRHCNuisBANFF00
    { 0.4,  0.5 },        // kRHCNuisBANFF01
    { 0.5,  0.6 },        // kRHCNuisBANFF02
    { 0.6,  0.7 },        // kRHCNuisBANFF03
    { 0.7,  1.0 },        // kRHCNuisBANFF04
    { 1.0,  1.5 },        // kRHCNuisBANFF05
    { 1.5,  2.5 },        // kRHCNuisBANFF06
    { 2.5,  3.5 },        // kRHCNuisBANFF07
    { 3.5,  5.0 },        // kRHCNuisBANFF08
    { 5.0,  7.0 },        // kRHCNuisBANFF09
    { 7.0, 30.0 },        // kRHCNuisBANFF10
    { 0.0,  1.5 },        // kRHCNuisBANFF11
    { 1.5, 30.0 },        // kRHCNuisBANFF12
    { 0.0,  0.5 },        // kRHCNuisBANFF13
    { 0.5,  0.7 },        // kRHCNuisBANFF14
    { 0.7,  0.8 },        // kRHCNuisBANFF15
    { 0.8,  1.5 },        // kRHCNuisBANFF16
    { 1.5,  2.5 },        // kRHCNuisBANFF17
    { 2.5,  4.0 },        // kRHCNuisBANFF18
    { 4.0, 30.0 },        // kRHCNuisBANFF19
    { 0.0,  2.5 },        // kRHCNuisBANFF20
    { 2.5, 30.0 },        // kRHCNuisBANFF21
#endif
#if BANFF_BINS == 1
    { 0.0,  0.4 },        // kRHCNuisBANFF00
    { 0.4,  0.5 },        // kRHCNuisBANFF01
    { 0.5,  0.6 },        // kRHCNuisBANFF02
    { 0.6,  0.7 },        // kRHCNuisBANFF03
    { 0.7,  1.0 },        // kRHCNuisBANFF04
    { 1.0,  1.5 },        // kRHCNuisBANFF05
    { 1.5,  2.5 },        // kRHCNuisBANFF06
    { 2.5,  3.5 },        // kRHCNuisBANFF07
    { 3.5,  5.0 },        // kRHCNuisBANFF08
    { 5.0,  7.0 },        // kRHCNuisBANFF09
    { 7.0, 30.0 },        // kRHCNuisBANFF10
    { 0.0,  0.7 },        // kRHCNuisBANFF11
    { 0.7,  1.0 },        // kRHCNuisBANFF12
    { 1.0,  1.5 },        // kRHCNuisBANFF13
    { 1.5,  2.5 },        // kRHCNuisBANFF14
    { 2.5, 30.0 },        // kRHCNuisBANFF15
    { 0.0,  0.5 },        // kRHCNuisBANFF16
    { 0.5,  0.7 },        // kRHCNuisBANFF17
    { 0.7,  0.8 },        // kRHCNuisBANFF18
    { 0.8,  1.5 },        // kRHCNuisBANFF19
    { 1.5,  2.5 },        // kRHCNuisBANFF20
    { 2.5,  4.0 },        // kRHCNuisBANFF21
    { 4.0, 30.0 },        // kRHCNuisBANFF22
    { 0.0,  2.5 },        // kRHCNuisBANFF23
    { 2.5, 30.0 },        // kRHCNuisBANFF24
#endif
    {  -1,   -1 },        // kRHCNuisBANFFMAQE, not a flux param
    {  -1,   -1 },        // kRHCNuisBANFFMARES, not a flux param
    {  -1,   -1 },        // kRHCNuisBANFFCCQEE0, not a flux param
#if BANFF_BINS != 3
    {  -1,   -1 },        // kRHCNuisBANFFCCQEE1, not a flux param
    {  -1,   -1 },        // kRHCNuisBANFFCCQEE2, not a flux param
#endif
    {  -1,   -1 },        // kRHCNuisBANFFCC1piE0, not a flux param
#if BANFF_BINS != 3
    {  -1,   -1 },        // kRHCNuisBANFFCC1piE1, not a flux param
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    {  -1,   -1 }         // kRHCNuisBANFFNC1pi0, not a flux param
#endif
#if BANFF_BINS == 4
    {  -1,   -1 },        // kRHCNuisBANFFMECNorm
    {  -1,   -1 }         // kRHCNuisBANFFPDDDoubleCount
#endif
  };


} // conv
} // t2kana
  
#endif 

