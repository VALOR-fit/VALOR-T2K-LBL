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

31-12-2013  Nick
Add samples, modes and systematic parameters for reverse horn current running

*/
//____________________________________________________________________________

#ifndef enumerations_h
#define enumerations_h

#include "AnalysisOptions.h"
#include "PhysParamId.h"    
#include "SuperkModes.h"    
#include "SuperkSamples.h"  
#include "ConfLevel.h"      
#include "PDG.h"            

namespace t2kana {
namespace conv   {

//
// ------------------------------------------------------------------------
// Enumerations
// ------------------------------------------------------------------------
//

// Physics hypotheses
typedef enum EPhysHypothesis 
{   
   kHypUndefined = -1,
   kHypDecay,             // neutrino decay 
   kHypDecoh,             // neutrino decoherence 
   kHypNuMuDisapp2F,      // muon-neutrino disappearance, simple 2-flavour oscillation framework 
   kHypOsc3F,             // neutrino oscillations, 3-flavour framework
   kHypOscSterile4F,      // neutrino oscillations with steriles, 3+1
   kHypOscSterile5F,      // neutrino oscillations with steriles, 3+2
   kHypNSI3F              // Nonstandard matter interactions
} 
PhysHypothesis_t;

// Combined or single analysis 
typedef enum EOscAnaType 
{   
   kAnaUndefined = -1, 
   kAnaCombined,
   kAnaNuMuDisappOnly,
   kAnaNuEAppOnly,
   kAnaRHCCombined,
   kAnaRHCNuMuDisappOnly,
   kAnaRHCNuEAppOnly,
   kAnaFHCAndRHCNuMuDisappOnly,
   kAnaFHCAndRHCNuEAppOnly,
   kAnaFHCAndRHCCombined
} 
OscAnaType_t;

// Neutrino mass hierarchy
typedef enum ENuMassHierarchy 
{   
   kNmhUndefined = -1, 
   kNmhNormal,
   kNmhInverted
} 
NuMassHierarchy_t;

// Choices of oscillation mixing parameter
typedef enum ENuOscMixParam
{   
   kNuMixUndefined = -1, 
   kNuMixSinsq2Theta,      // sin^2(2*theta)
   kNuMixSinsqTheta        // sin^2(theta), preferred choice in global nue appearance + numu disappearance fit
} 
NuOscMixParam_t;

// Choices of oscillation mixing octant
typedef enum ENuOscOctant
{   
   kNuOctantUndefined = -1, 
   kNuOctantOne,      // sin^2(theta_{23}) <= 0.5
   kNuOctantTwo       // sin^2(theta_{23}) >= 0.5
} 
NuOscOctant_t;

// Choices of squared mass splittings                               
typedef enum ENuOscDmsqScheme
{
  kNuDmsqUndefined = -1,
  kNuDmsq12andAbs23,  // use |Dm^2_{23}| and Dm^2_{12} as the 2 independent squared mass splittings                                                                            
  kNuDmsqFogliLisi,   // scheme used by Fogli and Lisi, dm^2 = Dm^2_{12} and Dm^2 = |m^2_{3} - (m^2_{1} + m^2_{1})/2|                                                            
  kNuDmsq12andAbs31   // use |Dm^2_{31}| and Dm^2_{12} as the 2 independent squared mass splittings                                                                            
}
 NuOscDmsqScheme_t;

// Choices of deltaCP                              
typedef enum ENuOscDeltaCP
{
  kNuDeltaCPUndefined = -1,
  kNuDeltaCPPhase,  // use deltaCP 
  kNuDeltaCPSin   // use sin(deltaCP)
}
 NuOscDeltaCP_t;

// Oscillation channels
typedef enum EOscChannelId 
{   
   kOscUndefined = -1, 
   kOscNuE2NuE,
   kOscNuE2NuMu,
   kOscNuE2NuTau,
   kOscNuMu2NuE,
   kOscNuMu2NuMu,
   kOscNuMu2NuTau,
   kOscNuTau2NuE,
   kOscNuTau2NuMu,
   kOscNuTau2NuTau,
   kOscNuEBar2NuEBar,
   kOscNuEBar2NuMuBar,
   kOscNuEBar2NuTauBar,
   kOscNuMuBar2NuEBar,
   kOscNuMuBar2NuMuBar,
   kOscNuMuBar2NuTauBar,
   kOscNuTauBar2NuEBar,
   kOscNuTauBar2NuMuBar,
   kOscNuTauBar2NuTauBar,
   //oscillations to or from sterile neutrinos NuS4 or NuS5
   kOscNuE2NuS4,
   kOscNuE2NuS5,
   kOscNuMu2NuS4,
   kOscNuMu2NuS5,
   kOscNuTau2NuS4,
   kOscNuTau2NuS5,
   kOscNuS42NuE,
   kOscNuS52NuE,
   kOscNuS42NuMu,
   kOscNuS52NuMu,
   kOscNuS42NuTau,
   kOscNuS52NuTau,
   kOscNuEBar2NuS4Bar,
   kOscNuEBar2NuS5Bar,
   kOscNuMuBar2NuS4Bar,
   kOscNuMuBar2NuS5Bar,
   kOscNuTauBar2NuS4Bar,
   kOscNuTauBar2NuS5Bar,
   kOscNuS4Bar2NuEBar,
   kOscNuS5Bar2NuEBar,
   kOscNuS4Bar2NuMuBar,
   kOscNuS5Bar2NuMuBar,
   kOscNuS4Bar2NuTauBar,
   kOscNuS5Bar2NuTauBar
} 
OscChannelId_t;

// Type of SK MC template
typedef enum ESKMCTemplateType
{   
   kSKMCtUndefined = -1, 
   kSKMCtOri, // nominal MC
   kSKMCtNom, // nominal MC + tuning (ext) weights 
   kSKMCtTun, // nominal MC + tuning (ext) weights + tuning (T2K) weights
   kSKMCtSys, // nominal MC + tuning (ext) weights + tuning (T2K) weights + systematic variations
   kSKMCtFin  // nominal MC + tuning (ext) weights + tuning (T2K) weights + systematic variations + physics hypothesis 
} 
SKMCTemplateType_t;

// Enumeration of flux histograms types held by NuFlux
typedef enum EFluxType
{   
   kFlxUndefined      = -1, 
   kFlxNominal        =  0,  // flux used for MC production
   kFlxExtTuning      =  1,  // flux used for MC production + tuning wghts based on external data (NA61,...)
   kFlxFinalTuning    =  2,  // flux used for MC production + tuning wghts based on external data (NA61,...) + tuning wghts based on T2K data (ND280 data fit by BANFF,...)
   kFlxCurrent        =  3,  // current flux from systematic param variation of kFlxFinalTuning flux
   kRHCFlxNominal     =  4,  // RHC flux used for MC production
   kRHCFlxExtTuning   =  5,  // RHC flux used for MC production + tuning wghts based on external data (NA61,...)
   kRHCFlxFinalTuning =  6,  // RHC flux used for MC production + tuning wghts based on external data (NA61,...) + tuning wghts based on T2K data (ND280 data fit by BANFF,...)
   kRHCFlxCurrent     =  7   // RHC current flux from systematic param variation of kFlxFinalTuning flux
} 
FluxType_t;

} // conv
} // t2kana
  
#endif 

