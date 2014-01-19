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

namespace t2kana {
namespace conv   {

//
// ------------------------------------------------------------------------
// Enumerations
// ------------------------------------------------------------------------
//

// Confidence level
typedef enum EConfLevel {
   kCLUndefined = -1, 
   kCL68, 
   kCL90, 
   kCL95, 
   kCL99
} 
ConfLevel_t; 

// Physics hypotheses
typedef enum EPhysHypothesis 
{   
   kHypUndefined = -1,
   kHypDecay,             // neutrino decay 
   kHypDecoh,             // neutrino decoherence 
   kHypNuMuDisapp2F,      // muon-neutrino disappearance, simple 2-flavour oscillation framework 
   kHypOsc3F              // neutrino oscillations, 3-flavour framework 
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
   kOscNuTauBar2NuTauBar
} 
OscChannelId_t;

// Enumeration of event samples that may be included in the analysis
// Actual choice for current analysis is given in `Configuration' section.
static const int kNSuperkSamples  = 4;  // number of samples 
static const int kns = kNSuperkSamples; // as above - shorter name
typedef enum ESuperkSample 
{
   kSKsUndefined = -1, 
   kSKsFCFV1RingMuLike, // fully-contained fiducial-volume single-ring mu-like
   kSKsFCFV1RingELike,  // fully-contained fiducial-volume single-ring e-like
   kSKsRHCFCFV1RingMuLike, //fully-contained fiducial-volume single-ring mu-like with reverse horn current 
   kSKsRHCFCFV1RingELike  //fully-contained fiducial-volume single-ring e-like with reverse horn current 
} 
SuperkSample_t; 

// Enumeration of (true) reaction modes that may be considered separately within each SuperK sample
// Actual choice for current analysis is given in `Configuration' section.
static const int kNSuperkModes = 53;  // number of reaction modes kept track of within each sample
static const int knm = kNSuperkModes; // as above - shorter name
typedef enum ESuperkMode 
{
   kSKmUndefined = -1, 
   kSKmNuMuCCQE,         // numu CC QE
   kSKmNuMuCC1pi,        // numu CC 1pi
   kSKmNuMuCCcoh,        // numu CC coherent
   kSKmNuMuCCoth,        // numu CC other
   kSKmNuMuNC1piC,       // numu NC 1 charged pi
   kSKmNuMuNC1pi0,       // numu NC 1pi0
   kSKmNuMuNCcoh,        // numu NC coherent
   kSKmNuMuNCoth_mu,     // numu NC other, nu_mu appearance analysis definition
   kSKmNuMuNCoth_e,      // numu NC other, nu_e appearance analysis definition
   kSKmOscNuMuCCQE,      // oscillated numu CCQE
   kSKmOscNuMuCC1pi,     // oscillated numu CC 1pi
   kSKmOscNuMuCCcoh,     // oscillated numu coh
   kSKmOscNuMuCCoth,     // oscillated numu CC other
   kSKmNuMuBarCCQE,      // numubar CCQE
   kSKmNuMuBarCC1pi,     // numubar CC 1pi
   kSKmNuMuBarCCcoh,     // numubar CC coh
   kSKmNuMuBarCCoth,     // numubar CC other
   kSKmNuMuBarNC1piC,    // numubar NC 1 charged pi
   kSKmNuMuBarNC1pi0,    // numubar NC 1pi0
   kSKmNuMuBarNCcoh,     // numubar NC coherent
   kSKmNuMuBarNCoth_mu,  // numubar NC other, nu_mu appearance analysis definition
   kSKmNuMuBarNCoth_e,   // numubar NC other, nu_e appearance analysis definition
   kSKmOscNuMuBarCCQE,   // oscillated numubar CCQE
   kSKmOscNuMuBarCC1pi,  // oscillated numubar CC 1pi
   kSKmOscNuMuBarCCcoh,  // oscillated numubar CC coh
   kSKmOscNuMuBarCCoth,  // oscillated numubar CC other
   kSKmNuECCQE,          // nue CC QE
   kSKmNuECC1pi,         // nue CC 1pi
   kSKmNuECCcoh,         // nue CC coherent
   kSKmNuECCoth,         // nue CC other
   kSKmNuENC1piC,        // nue NC 1 charged pi
   kSKmNuENC1pi0,        // nue NC 1pi0
   kSKmNuENCcoh,         // nue NC coherent
   kSKmNuENCoth_mu,      // nue NC other, nu_mu appearance analysis definition
   kSKmNuENCoth_e,       // nue NC other, nu_e appearance analysis definition
   kSKmOscNuEBarCCQE,    // oscillated nuebar CC QE
   kSKmOscNuEBarCC1pi,   // oscillated nuebar CC 1pi
   kSKmOscNuEBarCCcoh,   // oscillated nuebar CC coherent
   kSKmOscNuEBarCCoth,   // oscillated nuebar CC other
   kSKmNuEBarCCQE,       // nuebar CC QE
   kSKmNuEBarCC1pi,      // nuebar CC 1pi
   kSKmNuEBarCCcoh,      // nuebar CC coherent
   kSKmNuEBarCCoth,      // nuebar CC other
   kSKmNuEBarNC1piC,     // nuebar NC 1 charged pi
   kSKmNuEBarNC1pi0,     // nuebar NC 1pi0
   kSKmNuEBarNCcoh,      // nuebar NC coherent
   kSKmNuEBarNCoth_mu,   // nuebar NC other, nu_mu appearance analysis definition
   kSKmNuEBarNCoth_e,    // nuebar NC other, nu_e appearance analysis definition
   kSKmOscNuECCQE,       // oscillated nue CCQE
   kSKmOscNuECC1pi,      // oscillated nue CC 1pi
   kSKmOscNuECCcoh,      // oscillated nue CC coh
   kSKmOscNuECCoth,      // oscillated nue CC other
   kSKmMEC               // numu CC MEC
} 
SuperkMode_t; 

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

// Fit parameters
static const int knppmax = 200; // maximum number of physics params listed below (for booking large enough param arrays)

typedef enum EPhysParamId 
{   
   kPhypUndefined = -1, 
   kPhypMixing23,              // 23-mixing parameter (either sin^2(2*\theta_{23}) or sin^2(\theta_{23}) depending on config choices)
   kPhypMixing13,              // 13-mixing parameter (either sin^2(2*\theta_{13}) or sin^2(\theta_{13}) depending on config choices)
   kPhypMixing12,              // 12-mixing parameter (either sin^2(2*\theta_{12}) or sin^2(\theta_{12}) depending on config choices)
   kPhypDeltaCP,               // CP-invariance violating phase 
   kPhypDm23sq,                // modulus of squared-mass splitting \Deltam_{23}^{2} (always positive independent of mass hierarchy) 
   kPhypDm12sq,                // squared-mass splitting \Deltam_{12}^{2} 
   kPhypAlphaDecay,            //
   kPhypSinsq2ThetaDecay,      // 
   kPhypMusqDecoh,             //
   kPhypSinsq2ThetaDecoh,      //

   /*
    * List of 2011 analysis nuisance parameters 
    * -----------------------------------------
   */

   //********* FORWARD HORN CURRENT ************************ 

   /*
    * SuperK detector errors (energy scale and efficiencies)
   */   
   kNuisSuperkEnergyScale,     // SK energy scale
#if FSI_BINS != 3
   kNuisSuperkEff1RmuCCQEoth,  // SK 1Rmu CCQE efficiency (other)
   kNuisSuperkEff1RmuCCQErc0u, // SK 1Rmu CCQE efficiency (ring counting, Ereco bin 1)
   kNuisSuperkEff1RmuCCQErc1u, // SK 1Rmu CCQE efficiency (ring counting, Ereco bin 2)
   kNuisSuperkEff1RmuCCQErc2u, // SK 1Rmu CCQE efficiency (ring counting, Ereco bin 3)
   kNuisSuperkEff1RmuCCnonQE,  // SK 1Rmu CCnonQE efficiency
   kNuisSuperkEff1RmuNC,       // SK 1Rmu NC efficiency
   kNuisSuperkEff1RmuNueCC,    // SK 1Rmu nue CC efficiency
#endif

#if FSI_BINS != 2
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
#endif

   /*
    * (flux and cross-section) correlated parameters used in the BANFF fit
   */   
   kNuisBANFF00,               // # 00, SK     numu    flux,  E =  0.0 -  0.4 GeV
   kNuisBANFF01,               // # 01, SK     numu    flux,  E =  0.4 -  0.5 GeV
   kNuisBANFF02,               // # 02, SK     numu    flux,  E =  0.5 -  0.6 GeV
   kNuisBANFF03,               // # 03, SK     numu    flux,  E =  0.6 -  0.7 GeV
   kNuisBANFF04,               // # 04, SK     numu    flux,  E =  0.7 -  1.0 GeV
   kNuisBANFF05,               // # 05, SK     numu    flux,  E =  1.0 -  1.5 GeV
   kNuisBANFF06,               // # 06, SK     numu    flux,  E =  1.5 -  2.5 GeV
   kNuisBANFF07,               // # 07, SK     numu    flux,  E =  2.5 -  3.5 GeV
   kNuisBANFF08,               // # 08, SK     numu    flux,  E =  3.5 -  5.0 GeV
   kNuisBANFF09,               // # 09, SK     numu    flux,  E =  5.0 -  7.0 GeV
   kNuisBANFF10,               // # 10, SK     numu    flux,  E =  7.0 - 30.0 GeV
   kNuisBANFF11,               // # 11, different meanings depending on samples defined
   kNuisBANFF12,               // # 12, different meanings depending on samples defined
   kNuisBANFF13,               // # 13, different meanings depending on samples defined
   kNuisBANFF14,               // # 14, different meanings depending on samples defined
   kNuisBANFF15,               // # 15, different meanings depending on samples defined
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
   kNuisBANFF16,               // # 16, different meanings depending on samples defined
   kNuisBANFF17,               // # 17, different meanings depending on samples defined
   kNuisBANFF18,               // # 18, different meanings depending on samples defined
   kNuisBANFF19,               // # 19, different meanings depending on samples defined
   kNuisBANFF20,               // # 20, different meanings depending on samples defined
   kNuisBANFF21,               // # 21, different meanings depending on samples defined
#endif
#if BANFF_BINS == 1
   kNuisBANFF22,               // # 22, different meanings depending on samples defined
   kNuisBANFF23,               // # 23, different meanings depending on samples defined
   kNuisBANFF24,               // # 24, different meanings depending on samples defined
#endif
   kNuisBANFFMAQE,             // # M_{A}  QE
   kNuisBANFFMARES,            // # M_{A}  RES
   kNuisBANFFCCQEE0,           // # CCQE   normalization, 0 - 1.5 GeV
#if BANFF_BINS != 3
   kNuisBANFFCCQEE1,           // # CCQE   normalization, 1.5 -  3.5 GeV
   kNuisBANFFCCQEE2,           // # CCQE   normalization, 3.5 - 30.0 GeV
#endif
   kNuisBANFFCC1piE0,          // # CC1pi  normalization, 0 - 2.5 GeV
#if BANFF_BINS != 3
   kNuisBANFFCC1piE1,          // # CC1pi  normalization, 2.5 - 30.0 GeV
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4) 
   kNuisBANFFNC1pi0,           // # NC1pi0 normalization, 0 - 30.0 GeV
#endif
#if BANFF_BINS == 4
   kNuisBANFFMECNorm,          // MEC template normalisation
   kNuisBANFFPDDDoubleCount,   // PDD double counting
#endif

   /*
    * Cross-section parameters (excluding cross-section params which were included in BANFF fit and are correlated with the flux params)
    * These other paramers are taken to be uncorrelated and are applied to SK only.
   */   
   //systematics which use splines
   kNuisXSecFermiMom,          // Fermi momentum
   kNuisXSecSpecFunc,          // spectral function on/off
   kNuisXSecWshape,            // W-shape
   kNuisXSecShapeCCoth,        // CC other shape
   kNuisXSecPiLessDelta,       // Pion-less delta decay
   kNuisXSecBindE,             // Binding energy
   //systematics that are normalisations
   kNuisXSecNormCCcoh,         // CC coherent norm
   kNuisXSecNormNC1piC,        // NC 1piC norm
   kNuisXSecNormNCcoh,         // NC coh norm
   kNuisXSecNormNCoth,         // NC other norm
   kNuisXSecNormNuEToNuMu,     // nue / numu cross section
   kNuisXSecNormAntiNuToNu,    // antinu / nu cross section

   /*
    * FSI parameters.
    * Assuming here a 10x10 covariance matrix will be made available.
    * Each parameter is a scaling factor to apply at some Ereco range (to be defined) for some
    * SK reaction mode (or the total? to be defined)
    */

   // FSI matrix 6x6 for 1Rmu
#if FSI_BINS != 3
   kNuisFSI000,
   kNuisFSI001,
   kNuisFSI002,
   kNuisFSI003,
   kNuisFSI004,
   kNuisFSI005,
#endif
   // FSI matrix 12x12 for 1Re
#if FSI_BINS != 2
   kNuisFSI006,
   kNuisFSI007,
   kNuisFSI008,
   kNuisFSI009,
   kNuisFSI010,
   kNuisFSI011,
   kNuisFSI012,
   kNuisFSI013,
   kNuisFSI014,
   kNuisFSI015,
   kNuisFSI016,
   kNuisFSI017,
#endif

   // SKDET+FSI matrix 18x18
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
   kNuisSKDetFSI017,
#endif

   //********* REVERSE HORN CURRENT **************************

   /*
    * Special RHC to FHC errors
    */

   kNuisRHCToFHC,

   /*
    * SuperK detector errors (energy scale and efficiencies)
    */
   
   kRHCNuisSuperkEnergyScale,     // SK energy scale
   
   /*                                               
    * (flux and cross-section) correlated parameters used in the BANFF fit               
    */
   
   kRHCNuisBANFF00,               // # 00, SK     numu    flux,  E =  0.0 -  0.4 GeV
   kRHCNuisBANFF01,               // # 01, SK     numu    flux,  E =  0.4 -  0.5 GeV
   kRHCNuisBANFF02,               // # 02, SK     numu    flux,  E =  0.5 -  0.6 GeV   
   kRHCNuisBANFF03,               // # 03, SK     numu    flux,  E =  0.6 -  0.7 GeV
   kRHCNuisBANFF04,               // # 04, SK     numu    flux,  E =  0.7 -  1.0 GeV  
   kRHCNuisBANFF05,               // # 05, SK     numu    flux,  E =  1.0 -  1.5 GeV
   kRHCNuisBANFF06,               // # 06, SK     numu    flux,  E =  1.5 -  2.5 GeV 
   kRHCNuisBANFF07,               // # 07, SK     numu    flux,  E =  2.5 -  3.5 GeV
   kRHCNuisBANFF08,               // # 08, SK     numu    flux,  E =  3.5 -  5.0 GeV
   kRHCNuisBANFF09,               // # 09, SK     numu    flux,  E =  5.0 -  7.0 GeV
   kRHCNuisBANFF10,               // # 10, SK     numu    flux,  E =  7.0 - 30.0 GeV 
   kRHCNuisBANFF11,               // # 11, different meanings depending on samples defined
   kRHCNuisBANFF12,               // # 12, different meanings depending on samples defined
   kRHCNuisBANFF13,               // # 13, different meanings depending on samples defined       
   kRHCNuisBANFF14,               // # 14, different meanings depending on samples defined      
   kRHCNuisBANFF15,               // # 15, different meanings depending on samples defined       
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
   kRHCNuisBANFF16,               // # 16, different meanings depending on samples defined
   kRHCNuisBANFF17,               // # 17, different meanings depending on samples defined
   kRHCNuisBANFF18,               // # 18, different meanings depending on samples defined
   kRHCNuisBANFF19,               // # 19, different meanings depending on samples defined       
   kRHCNuisBANFF20,               // # 20, different meanings depending on samples defined
   kRHCNuisBANFF21,               // # 21, different meanings depending on samples defined      
#endif
#if BANFF_BINS == 1
   kRHCNuisBANFF22,               // # 22, different meanings depending on samples defined
   kRHCNuisBANFF23,               // # 23, different meanings depending on samples defined       
   kRHCNuisBANFF24,               // # 24, different meanings depending on samples defined
#endif
   kRHCNuisBANFFMAQE,             // # M_{A}  QE                                     
   kRHCNuisBANFFMARES,            // # M_{A}  RES                                   
   kRHCNuisBANFFCCQEE0,           // # CCQE   normalization, 0 - 1.5 GeV   
#if BANFF_BINS != 3
   kRHCNuisBANFFCCQEE1,           // # CCQE   normalization, 1.5 -  3.5 GeV        
   kRHCNuisBANFFCCQEE2,           // # CCQE   normalization, 3.5 - 30.0 GeV        
#endif
   kRHCNuisBANFFCC1piE0,          // # CC1pi  normalization, 0 - 2.5 GeV           
#if BANFF_BINS != 3
   kRHCNuisBANFFCC1piE1,          // # CC1pi  normalization, 2.5 - 30.0 GeV       
#endif
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
   kRHCNuisBANFFNC1pi0,           // # NC1pi0 normalization, 0 - 30.0 GeV         
#endif
#if BANFF_BINS == 4
   kRHCNuisBANFFMECNorm,          // MEC template normalisation                     
   kRHCNuisBANFFPDDDoubleCount,   // PDD double counting                          
#endif
   
   /*                                                                             
    * Cross-section parameters (excluding cross-section params which were included in BANFF fit and are correlated with the flux params)
    * These other paramers are taken to be uncorrelated and are applied to SK only.     
    */
   
   //systematics which use splines                                               
   kRHCNuisXSecFermiMom,          // Fermi momentum                                 
   kRHCNuisXSecSpecFunc,          // spectral function on/off         
   kRHCNuisXSecWshape,            // W-shape                          
   kRHCNuisXSecShapeCCoth,        // CC other shape               
   kRHCNuisXSecPiLessDelta,       // Pion-less delta decay                
   kRHCNuisXSecBindE,             // Binding energy               

   //systematics that are normalisations                                         
   kRHCNuisXSecNormCCcoh,         // CC coherent norm                              
   kRHCNuisXSecNormNC1piC,        // NC 1piC norm                                   
   kRHCNuisXSecNormNCcoh,         // NC coh norm                                   
   kRHCNuisXSecNormNCoth,         // NC other norm                       
   kRHCNuisXSecNormNuEToNuMu,     // nue / numu cross section                  
   kRHCNuisXSecNormAntiNuToNu,    // antinu / nu cross section        
   
   /*                                                                          
    * SK and FSI parameters.                                 
    * Assuming here a 10x10 covariance matrix will be made available.            
    * Each parameter is a scaling factor to apply at some Ereco range (to be defined) for some
    * SK reaction mode (or the total? to be defined)                              
    */
   
   // SKDET+FSI matrix 18x18                                                  
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
   
} 
PhysParamId_t;

} // conv
} // t2kana
  
#endif 

