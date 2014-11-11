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

31-12-2013  Nick
Added extra arrays for reverse horn current modes

21-7-2014 Davide
Added samples for multiring forward horn current

*/
//____________________________________________________________________________

#ifndef sk_fsi_syst_conv_h
#define sk_fsi_syst_conv_h

#include "AnalysisOptions.h"
#include "Enumerations.h"
#include "SamplesAndModes.h"

namespace valor {
namespace conv   {

//
// ** SKDET+FSI systematics
//

// Number of SKDET+FSI systematic params

#if FSI_BINS == 1
static const int kNumOfSKDetFSIParams     = 18;
static const int kNumOfSKDetFSIParams1Rmu = 6;
static const int kNumOfSKDetFSIParams1Re  = 12;
static const int kNumOfRHCSKDetFSIParams     = 18;
static const int kNumOfRHCSKDetFSIParams1Rmu = 6;
static const int kNumOfRHCSKDetFSIParams1Re  = 12;
/// NEW MULTIRING 2
static const int kNumOfSKDetFSIParamsMultiRing = 18;
static const int kNumOfSKDetFSIParamsMultiRmu = 6;
static const int kNumOfSKDetFSIParamsMultiRe  = 12;
////
#elif FSI_BINS == 2
static const int kNumOfSKDetFSIParams     = 6;
static const int kNumOfSKDetFSIParams1Rmu = 6;
static const int kNumOfSKDetFSIParams1Re  = 0;
static const int kNumOfRHCSKDetFSIParams     = 6;
static const int kNumOfRHCSKDetFSIParams1Rmu = 6;
static const int kNumOfRHCSKDetFSIParams1Re  = 0;
/// NEW MULTIRING 2
static const int kNumOfSKDetFSIParamsMultiRing     = 6;
static const int kNumOfSKDetFSIParamsMultiRmu = 6;
static const int kNumOfSKDetFSIParamsMultiRe  = 0;
/// 
#elif FSI_BINS == 3
static const int kNumOfSKDetFSIParams     = 12;
static const int kNumOfSKDetFSIParams1Rmu = 0;
static const int kNumOfSKDetFSIParams1Re  = 12;
static const int kNumOfRHCSKDetFSIParams     = 12;
static const int kNumOfRHCSKDetFSIParams1Rmu = 0;
static const int kNumOfRHCSKDetFSIParams1Re  = 12;
/// NEW MULTIRING 2
static const int kNumOfSKDetFSIParamsMultiRing     = 12;
static const int kNumOfSKDetFSIParamsMultiRmu = 0;
static const int kNumOfSKDetFSIParamsMultiRe  = 12;
///
#endif

// The following is the list of the above kNumOfFSIParams
static const PhysParamId_t kListOfSKDetFSIParams[kNumOfSKDetFSIParams] = 
{
#if FSI_BINS != 3
    kNuisSKDetFSI000,
    kNuisSKDetFSI001,
    kNuisSKDetFSI002,
    kNuisSKDetFSI003,
    kNuisSKDetFSI004,
    kNuisSKDetFSI005,
#endif
#if FSI_BINS != 2
    kNuisSKDetFSI006,
    kNuisSKDetFSI007,
    kNuisSKDetFSI008,
    kNuisSKDetFSI009,
    kNuisSKDetFSI010,
    kNuisSKDetFSI011,
    kNuisSKDetFSI012,
    kNuisSKDetFSI013,
    kNuisSKDetFSI014,
    kNuisSKDetFSI015,
    kNuisSKDetFSI016,
    kNuisSKDetFSI017
#endif
};

static const PhysParamId_t kListOfRHCSKDetFSIParams[kNumOfRHCSKDetFSIParams] = 
  {
#if FSI_BINS != 3
    kRHCNuisSKDetFSI000,
    kRHCNuisSKDetFSI001,
    kRHCNuisSKDetFSI002,
    kRHCNuisSKDetFSI003,
    kRHCNuisSKDetFSI004,
    kRHCNuisSKDetFSI005,
#endif
#if FSI_BINS != 2
    kRHCNuisSKDetFSI006,
    kRHCNuisSKDetFSI007,
    kRHCNuisSKDetFSI008,
    kRHCNuisSKDetFSI009,
    kRHCNuisSKDetFSI010,
    kRHCNuisSKDetFSI011,
    kRHCNuisSKDetFSI012,
    kRHCNuisSKDetFSI013,
    kRHCNuisSKDetFSI014,
    kRHCNuisSKDetFSI015,
    kRHCNuisSKDetFSI016,
    kRHCNuisSKDetFSI017
#endif
  };

/// NEW MULTIRING 2
static const PhysParamId_t kListOfMultiRingSKDetFSIParams[kNumOfSKDetFSIParamsMultiRing] = 
  {
#if FSI_BINS != 3
    kMultiRingNuisSKDetFSI000,
    kMultiRingNuisSKDetFSI001,
    kMultiRingNuisSKDetFSI002,
    kMultiRingNuisSKDetFSI003,
    kMultiRingNuisSKDetFSI004,
    kMultiRingNuisSKDetFSI005,
#endif
#if FSI_BINS != 2
    kMultiRingNuisSKDetFSI006,
    kMultiRingNuisSKDetFSI007,
    kMultiRingNuisSKDetFSI008,
    kMultiRingNuisSKDetFSI009,
    kMultiRingNuisSKDetFSI010,
    kMultiRingNuisSKDetFSI011,
    kMultiRingNuisSKDetFSI012,
    kMultiRingNuisSKDetFSI013,
    kMultiRingNuisSKDetFSI014,
    kMultiRingNuisSKDetFSI015,
    kMultiRingNuisSKDetFSI016,
    kMultiRingNuisSKDetFSI017
#endif
  };
////

static const double kSKDetFSIParamsErecoCut[kNumOfSKDetFSIParams][2] =
{
#if FSI_BINS != 3
  { 0.0,   0.4 }, // FSI000
  { 0.4,   1.1 }, // FSI001
  { 1.1,  30.0 }, // FSI002
  { 0.0,  30.0 }, // FSI003
  { 0.0,  30.0 }, // FSI004
  { 0.0,  30.0 }, // FSI005
#endif
#if FSI_BINS != 2
  { 0.0,  0.35 }, // FSI1Re000
  { 0.35,  0.8 }, // FSI1Re001
  { 0.8,  1.25 }, // FSI1Re002
  { 0.0,  0.35 }, // FSI1Re003
  { 0.35,  0.8 }, // FSI1Re004
  { 0.8,  1.25 }, // FSI1Re005
  { 0.0,  0.35 }, // FSI1Re006
  { 0.35,  0.8 }, // FSI1Re007
  { 0.8,  1.25 }, // FSI1Re008
  { 0.0,  0.35 }, // FSI1Re009
  { 0.35,  0.8 }, // FSI1Re010
  { 0.8,  1.25 }  // FSI1Re011
#endif
};

static const double kRHCSKDetFSIParamsErecoCut[kNumOfRHCSKDetFSIParams][2] =
  {
#if FSI_BINS != 3
    { 0.0,   0.4 }, // FSI000
    { 0.4,   1.1 }, // FSI001
    { 1.1,  30.0 }, // FSI002
    { 0.0,  30.0 }, // FSI003
    { 0.0,  30.0 }, // FSI004
    { 0.0,  30.0 }, // FSI005
#endif
#if FSI_BINS != 2
    { 0.0,  0.35 }, // FSI1Re000
    { 0.35,  0.8 }, // FSI1Re001
    { 0.8,  1.25 }, // FSI1Re002
    { 0.0,  0.35 }, // FSI1Re003
    { 0.35,  0.8 }, // FSI1Re004
    { 0.8,  1.25 }, // FSI1Re005
    { 0.0,  0.35 }, // FSI1Re006
    { 0.35,  0.8 }, // FSI1Re007
    { 0.8,  1.25 }, // FSI1Re008
    { 0.0,  0.35 }, // FSI1Re009
    { 0.35,  0.8 }, // FSI1Re010
    { 0.8,  1.25 }  // FSI1Re011
#endif
  };

/// NEW MULTIRING 2
static const double kMultiRingSKDetFSIParamsErecoCut[kNumOfSKDetFSIParamsMultiRing][2] =
  {
#if FSI_BINS != 3
    { 0.0,   0.4 }, // FSI000
    { 0.4,   1.1 }, // FSI001
    { 1.1,  30.0 }, // FSI002
    { 0.0,  30.0 }, // FSI003
    { 0.0,  30.0 }, // FSI004
    { 0.0,  30.0 }, // FSI005
#endif
#if FSI_BINS != 2
    { 0.0,  0.35 }, // FSI1Re000
    { 0.35,  0.8 }, // FSI1Re001
    { 0.8,  1.25 }, // FSI1Re002
    { 0.0,  0.35 }, // FSI1Re003
    { 0.35,  0.8 }, // FSI1Re004
    { 0.8,  1.25 }, // FSI1Re005
    { 0.0,  0.35 }, // FSI1Re006
    { 0.35,  0.8 }, // FSI1Re007
    { 0.8,  1.25 }, // FSI1Re008
    { 0.0,  0.35 }, // FSI1Re009
    { 0.35,  0.8 }, // FSI1Re010
    { 0.8,  1.25 }  // FSI1Re011
#endif
  };
///

/// NEW MULTIRING 2

#if FSI_BINS == 1

//Todo: make this work if kNumOfSKDetFSIParams is different from kNumOfRHCSKDetFSIParams and kNumOfSKDetFSIParamsMultiRing

static const bool kSKDetFSIParamToSKMode[2][knm][kNumOfSKDetFSIParams] = {
  // > 1 muon-like ring 
  {
    // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
    //
    // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1Rmu
  ,
  // > 1 electron-like ring 
  {
    // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
    //
    // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth_e
    
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1Re
};

 static const bool kRHCSKDetFSIParamToSKMode[2][knm][kNumOfRHCSKDetFSIParams] = { 
  // > 1 muon-like ring RHC 
  {
    // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
    //
    // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1RmuRHC
  ,
  // > 1 electron-like ring RHC
  {
    // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
    //
    // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth_e
    
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1ReRHC
};

static const bool kMultiRingSKDetFSIParamToSKMode[2][knm][kNumOfSKDetFSIParamsMultiRing] = {
  // > Multi muon-like ring 
  {
    // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
    //
    // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//MultiRmu
  ,
  // > Multi electron-like ring 
  {
    // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
    //
    // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth_e
    
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//MultiRe
};


#elif FSI_BINS == 2

//Todo: make this work if kNumOfSKDetFSIParams is different from kNumOfRHCSKDetFSIParams
static const bool kSKDetFSIParamToSKMode[2][knm][kNumOfSKDetFSIParams] = {
  // > 1 muon-like ring 
  {
    // FSI000 FSI001 FSI002  FSI003 FSI004 FSI005
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1Rmu
  ,
  // > 1 electron-like ring 
  {
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e
    
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1Re
};


static const bool kRHCSKDetFSIParamToSKMode[2][knm][kNumOfRHCSKDetFSIParams] = {
  // > 1 muon-like ring RHC
  {
    // FSI000 FSI001 FSI002  FSI003 FSI004 FSI005
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1RmuRHC
  ,
  // > 1 electron-like ring RHC 
  {
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e
    
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1ReRHC
};


static const bool kMultiRingSKDetFSIParamToSKMode[2][knm][kNumOfSKDetFSIParamsMultiRing] = {
  // > Multi muon-like ring 
  {
    // FSI000 FSI001 FSI002  FSI003 FSI004 FSI005
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCcoh
    //{  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//MultiRmu
  ,
  // > Multi electron-like ring 
  {
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e
    
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//MultiRe
};



#endif


///


















// OLD

/* 

#if FSI_BINS == 1

//Todo: make this work if kNumOfSKDetFSIParams is different from kNumOfRHCSKDetFSIParams and kNumOfSKDetFSIParamsMultiRing

//static const bool kSKDetFSIApplySKDetFSIParamToSKMode[kns][knm][kNumOfSKDetFSIParams] = {
static const bool kSKDetFSIParamToSKMode[kns][knm][kNumOfSKDetFSIParams] = {

  // > 1 muon-like ring 
  {
    // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
    //
    // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1Rmu
  ,
  // > 1 electron-like ring 
  {
    // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
    //
    // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth_e
    
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1Re
  ,
  // > 1 muon-like ring RHC
  {
    // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
    //
    // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1RmuRHC
  ,
  // > 1 electron-like ring RHC
  {
    // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
    //
    // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth_e
    
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1ReRHC

  /// NEW MULTRING
   ,
  // > Multi muon-like ring
  {
    // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
    //
    // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_mu
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//MultiRmu
  ,
  // > Multi electron-like ring
  {
    // CCQEmu(bar) x 3  / CCnQEmu(bar) / CCe / NC
    //
    // oscnueCC x 3 / numuCC x 3 / nueCC x 3 / NC x 3
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth_e
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth_e
    
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//MultiRe

  ////
};
#elif FSI_BINS == 2

//Todo: make this work if kNumOfSKDetFSIParams is different from kNumOfRHCSKDetFSIParams

//static const bool kSKDetFSIApplySKDetFSIParamToSKMode[kns][knm][kNumOfSKDetFSIParams] = {
static const bool kSKDetFSIParamToSKMode[kns][knm][kNumOfSKDetFSIParams] = {
  // > 1 muon-like ring 
  {
    // FSI000 FSI001 FSI002  FSI003 FSI004 FSI005
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCoth
    
    {  1, 1, 1, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1Rmu
  ,
  // > 1 electron-like ring 
  {
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e
    
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1Re
  ,
  // > 1 muon-like ring RHC
  {
    // FSI000 FSI001 FSI002  FSI003 FSI004 FSI005
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCoth
    
    {  1, 1, 1, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1RmuRHC
  ,
  // > 1 electron-like ring RHC
  {
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e
    
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//1ReRHC  

  /// NEW MULTIRING 
  ,
  // > Multi muon-like ring 
  {
    // FSI000 FSI001 FSI002  FSI003 FSI004 FSI005
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCoth
    
    {  1, 1, 1, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//MultiRmu
  ,
  // > Multi electron-like ring 
  {
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e
     
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }//MultiRe  

///
};
#endif
*/






















//
// ** FSI systematics
//

// Number of FSI systematic params

#if FSI_BINS == 1
static const int kNumOfFSIParams     = 18;
static const int kNumOfFSIParams1Rmu = 6;
static const int kNumOfFSIParams1Re  = 12;
#elif FSI_BINS == 2
static const int kNumOfFSIParams     = 6;
static const int kNumOfFSIParams1Rmu = 6;
static const int kNumOfFSIParams1Re  = 0;
#elif FSI_BINS == 3
static const int kNumOfFSIParams     = 12;
static const int kNumOfFSIParams1Rmu = 0;
static const int kNumOfFSIParams1Re  = 12;
#endif

static const double kFSIParamsErecoCut[kNumOfFSIParams][2] =
{
#if FSI_BINS != 3
  { 0.0,   0.4 }, // FSI000
  { 0.4,   1.1 }, // FSI001
  { 1.1,  30.0 }, // FSI002
  { 0.0,  30.0 }, // FSI003
  { 0.0,  30.0 }, // FSI004
  { 0.0,  30.0 }, // FSI005
#endif
#if FSI_BINS != 2
  { 0.0,  0.35 }, // FSI1Re000
  { 0.35,  0.8 }, // FSI1Re001
  { 0.8,  1.25 }, // FSI1Re002
  { 0.0,  0.35 }, // FSI1Re003
  { 0.35,  0.8 }, // FSI1Re004
  { 0.8,  1.25 }, // FSI1Re005
  { 0.0,  0.35 }, // FSI1Re006
  { 0.35,  0.8 }, // FSI1Re007
  { 0.8,  1.25 }, // FSI1Re008
  { 0.0,  0.35 }, // FSI1Re009
  { 0.35,  0.8 }, // FSI1Re010
  { 0.8,  1.25 }  // FSI1Re011
#endif
};


#if FSI_BINS == 1
static const bool kFSIApplyFSIParamToSKMode[kns][knm][kNumOfFSIParams] = {
  // > 1 muon-like ring 
  {
    // FSI000 FSI001 FSI002  FSI003 FSI004 FSI005
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
#endif

    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth

    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  1, 1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmMEC

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring 
  {
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth_e

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
#endif

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth_e
    
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth

    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 0, 0,   1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  
};
#elif FSI_BINS == 2
static const bool kFSIApplyFSIParamToSKMode[kns][knm][kNumOfFSIParams] = {
  
  // Todo: add MEC for 1 e-like ring 

  // > 1 muon-like ring 
  {
    // FSI000 FSI001 FSI002  FSI003 FSI004 FSI005
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCoth
#endif
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 1, 1, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth
#endif

    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNCoth_mu
    
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth

    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCoth

    {  1, 1, 1, 0, 0, 0 }, // kSKmMEC

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring 
  {
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e

    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e

    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e

    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  
};
#elif FSI_BINS == 3
static const bool kFSIApplyFSIParamToSKMode[kns][knm][kNumOfFSIParams] = {
  // > 1 muon-like ring 
  {

    // Todo: add MEC for 1 e-lke ring

    // FSI000 FSI001 FSI002  FSI003 FSI004 FSI005
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmMEC

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring 
  {
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth_e

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif

    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
#endif

    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth_e
    
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth

    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  
};
#endif

#if FSI_BINS == 1
static const int kNumOfSKDetParams1Rmu = 8;
static const int kNumOfSKDetParams1Re = 12; //only the specific for 1Re to read the matrix, Ereco scale not included then
static const int kNumOfSKDetParams    = 20;
#elif FSI_BINS == 2
static const int kNumOfSKDetParams1Rmu = 8;
static const int kNumOfSKDetParams1Re  = 0; //only the specific for 1Re to read the matrix, Ereco scale not included then
static const int kNumOfSKDetParams     = 8;
#else
static const int kNumOfSKDetParams1Rmu = 0;
static const int kNumOfSKDetParams1Re  = 13; //only the specific for 1Re to read the matrix, Ereco scale not included then
static const int kNumOfSKDetParams     = 13; //12 + 1 energy scale
#endif

// The following is the list of the efficiency parameters in 1Re
#if FSI_BINS != 2
static const PhysParamId_t kListOfSKDet1ReParams[kNumOfSKDetParams1Re] = 
{
  kNuisSuperkEff1ReNueSigCC0,    // SK 1Re CC eff for osc nue Ereco bin 0
  kNuisSuperkEff1ReNueSigCC1,    // SK 1Re CC eff for osc nue Ereco bin 1
  kNuisSuperkEff1ReNueSigCC2,    // SK 1Re CC eff for osc nue Ereco bin 2
  kNuisSuperkEff1ReNumuCC0,      // SK 1Re CC eff for numu Ereco bin 0
  kNuisSuperkEff1ReNumuCC1,      // SK 1Re CC eff for numu Ereco bin 1
  kNuisSuperkEff1ReNumuCC2,      // SK 1Re CC eff for numu Ereco bin 2
  kNuisSuperkEff1ReNueCC0,       // SK 1Re CC eff for nue Ereco bin 0
  kNuisSuperkEff1ReNueCC1,       // SK 1Re CC eff for nue Ereco bin 1
  kNuisSuperkEff1ReNueCC2,       // SK 1Re CC eff for nue Ereco bin 2
  kNuisSuperkEff1ReNC0,          // SK 1Re NC eff Ereco bin 0
  kNuisSuperkEff1ReNC1,          // SK 1Re NC eff Ereco bin 1
  kNuisSuperkEff1ReNC2,          // SK 1Re NC eff Ereco bin 2
};



static const double kSKDet1ReParamsErecoCut[kNumOfSKDetParams1Re][2] =
{
  { 0.0,  0.35 }, // kNuisSuperkEff1ReNueSigCC0
  { 0.35,  0.8 }, // kNuisSuperkEff1ReNueSigCC1
  { 0.8,  1.25 }, // kNuisSuperkEff1ReNueSigCC2
  { 0.0,  0.35 }, // kNuisSuperkEff1ReNumuCC0
  { 0.35,  0.8 }, // kNuisSuperkEff1ReNumuCC1
  { 0.8,  1.25 }, // kNuisSuperkEff1ReNumuCC2
  { 0.0,  0.35 }, // kNuisSuperkEff1ReNueCC0
  { 0.35,  0.8 }, // kNuisSuperkEff1ReNueCC1
  { 0.8,  1.25 }, // kNuisSuperkEff1ReNueCC2
  { 0.0,  0.35 }, // kNuisSuperkEff1ReNC0
  { 0.35,  0.8 }, // kNuisSuperkEff1ReNC1
  { 0.8,  1.25 }  // kNuisSuperkEff1ReNC2
};

static const bool kSKDet1ReApplyParamToSKMode[kns][knm][kNumOfSKDetParams1Re] = {
  // > 1 muon-like ring 
  {
    // 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmMEC

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined

  }
  ,
  // > 1 electron-like ring 
  {
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuNCoth_e

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif

    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
#endif

    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuECCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0 }, // kSKmNuEBarCCoth
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1 }, // kSKmNuEBarNCoth_e
    
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
    

    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
    {  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined

  }
  
};
#endif

//                                
// Number of nuisance params parameters for 1Rmu, 1Re, 1Rmu+1Re analysis 
// 1Rmu : 57 nuisance params     
//   - 1 SK E scale              
//   - 7 SK efficiencies (3 CCQE-RC, CCQE-Oth, CCnonQE, NC, nue CC) 
//   - 31 BANFF flux             
//   - 8 BANFF XSec              
//   - 6 XSec parameters treated with splines
//   - 5 XSec normalisations      
//   - 15 FSI                     
// 1Re :                          
//   - 25 BANFF flux              
//   -  8 BANFF xsec              
// Joint :                       
//   - 25 BANFF flux            
//   -  8 BANFF xsec            
//                                 

// NOT USED AT ALL
//static const int kNumOfNuisanceParams1Rmu  = kNumOfParamsInBANFFfit + kNumOfSKDetFSIParams1Rmu + kNumOfUncXSecSyst;                                                                                
//static const int kNumOfNuisanceParams1Re   = kNumOfParamsInBANFFfit + kNumOfSKDetFSIParams1Re  + kNumOfUncXSecSyst - 1 + 1 ;                                                                       
//static const int kNumOfNuisanceParamsJoint = kNumOfParamsInBANFFfit + kNumOfSKDetFSIParams     + kNumOfXSecSyst - (kNumOfParamsInBANFFfit - kNumOfFluxParamsInBANFFfit) + 1;                      
///

 
// TODO
// It doesn't support RHC systematic parameters if not fully correlated w/ FHC sample --> to fix it 
#if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)
 
    #if defined(__MULTIRING_FHC_RUNNING__) && defined(__FHC_RUNNING__) // fhc all + skdetfsi multiring (no escale)
    static const int kNumOfNuisanceParams = kNumOfParamsInBANFFfit + kNumOfSKDetFSIParams + 1/*fhc energy scale*/  + kNumOfUncXSecSyst + kNumOfSKDetFSIParamsMultiRing + 1/*multiring escale*/;  // Tot # par FHC+MultiRing
    
    // NEW FLUX
    #elif defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) // fhc+rhc
    static const int kNumOfNuisanceParams = kNumOfParamsInBANFFfit + kNumOfSKDetFSIParams + 1/*energy scale*/  + kNumOfUncXSecSyst ; // Todo FHC+RHC with correlations
    //

    #elif defined(__FHC_RUNNING__)  // fhc only
    static const int kNumOfNuisanceParams = kNumOfParamsInBANFFfit + kNumOfSKDetFSIParams + 1/*energy scale*/  + kNumOfUncXSecSyst ;
    
    // NEW FLUX
    #elif defined(__RHC_RUNNING__) // rhc only
    static const int kNumOfNuisanceParams = kNumOfRHCParamsInBANFFfit + kNumOfRHCSKDetFSIParams + 1/*energy scale*/  + kNumOfUncXSecSyst ; // uncorr xsec is identical to fhc-only
    //
    
    #endif

#else
    #ifdef __MULTIRING_FHC_RUNNING__ // multiring fhc only 
    static const int kNumOfNuisanceParams = kNumOfParamsInBANFFfit + kNumOfSKDetFSIParamsMultiRing + 1/*energy scale*/  + kNumOfUncXSecSyst;
#else
    cerr << "No samples have been chosen for the analysis!!!" << endl;
    exit(1);    
    #endif
#endif // sample analysis option



} // conv
} // valor
  
#endif 

