//____________________________________________________________________________
/*!

\class    t2kana::conv

\brief    Conventions

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory
          
\created  Oct 20, 2010

Changes since 2011_v0.9.0:

*/
//____________________________________________________________________________

#ifndef _t2kana_conv_h_
#define _t2kana_conv_h_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

namespace t2kana {
namespace conv   {

//
// ------------------------------------------------------------------------
// Constants & commonly used codes
// ------------------------------------------------------------------------
//

// PDG Codes
static const int kPdgNuE      =  12;
static const int kPdgNuMu     =  14;
static const int kPdgNuTau    =  16;
static const int kPdgNuEBar   = -12;
static const int kPdgNuMuBar  = -14;
static const int kPdgNuTauBar = -16;

//
// Oscillation parameters, recent global fit values & 1 sigma errors and
// typical initialization and range values for fitter
//

// Dm^2 global fit results
static const double kGlobFitDm23sq           = 2.35E-3;
static const double kGlobFitDm12sq           = 7.58E-5;
static const double kGlobFitErrDm23sq        = 0.05; // ~5%
static const double kGlobFitErrDm12sq        = 0.03; // ~3%
// sinsq2theta global fit results
static const double kGlobFitSinsq2Theta23    = 1.00;
static const double kGlobFitSinsq2Theta13    = 0.10;
static const double kGlobFitSinsq2Theta12    = 0.86;
static const double kGlobFitErrSinsq2Theta23 = 0.15; // ~15%
static const double kGlobFitErrSinsq2Theta12 = 0.06; // ~6%
// sinsqtheta global fit results
static const double kGlobFitSinsqTheta23    = 0.420;
static const double kGlobFitSinsqTheta13    = 0.025;
static const double kGlobFitSinsqTheta12    = 0.312;
static const double kGlobFitErrSinsqTheta23 = 0.261;  // ~26% (1 sigma range)
static const double kGlobFitErrSinsqTheta12 = 0.110;  // ~11% (1 sigma range)
static const double kGlobFitErrSinsqTheta13 = 0.560;  // ~56% (1 sigma range)

// Dm^2 fitter init values and range
static const double kInitFitDm23sq           = kGlobFitDm23sq;
static const double kMinFitDm23sq            = 1.0E-3;
static const double kMaxFitDm23sq            = 4.0E-3;
static const double kInitFitDm12sq           = kGlobFitDm12sq;
static const double kMinFitDm12sq            = 6.0E-5;
static const double kMaxFitDm12sq            = 9.0E-5;
// sinsq2theta fitter init values and range
static const double kInitFitSinsq2Theta23    = 0.90;
static const double kMinFitSinsq2Theta23     = 0.60;
static const double kMaxFitSinsq2Theta23     = 1.00;
static const double kInitFitSinsq2Theta13    = kGlobFitSinsq2Theta13;
static const double kMinFitSinsq2Theta13     = 0.00;
static const double kMaxFitSinsq2Theta13     = 0.40;
static const double kInitFitSinsq2Theta12    = kGlobFitSinsq2Theta12;
static const double kMinFitSinsq2Theta12     = 0.70;
static const double kMaxFitSinsq2Theta12     = 1.00;
// sinsqtheta fitter init values and range
static const double kInitFitSinsqTheta23     = 0.90;
static const double kMinFitSinsqTheta23      = 0.60;
static const double kMaxFitSinsqTheta23      = 1.00;
static const double kInitFitSinsqTheta13     = kGlobFitSinsqTheta13;
static const double kMinFitSinsqTheta13      = 0.00;
static const double kMaxFitSinsqTheta13      = 0.40;
static const double kInitFitSinsqTheta12     = kGlobFitSinsqTheta12;
static const double kMinFitSinsqTheta12      = 0.70;
static const double kMaxFitSinsqTheta12      = 1.00;

//
// Other constants
//
static const double kASmallNum      = 1E-10;                         // 
static const double kOneMeV         = 0.001;                         // 
static const double kBaseline       = 295.335;                       // Tokai->Kamioka (km)
static const double kOscDimensFactor= 1.26693281;                    // Dimensional factor (GeV/eV^2/km) appearing in netrino oscillation formula
static const double kCrustDensity   = 2.7;                           // Average density of Earth's continental crust (gr / cm^3)
static const double kpi             = 3.141592653589793;             // pi
static const double khbar_SI        = 6.582119e-16;                  // hbar (eV s)
static const double kc_SI           = 2.9979246e10;                  // c (cm/s)
static const double khbarc_SI       = khbar_SI*kc_SI;                // hbar*c (eV * cm)
static const double khbarc3_SI      = khbarc_SI*khbarc_SI*khbarc_SI; // (hbar*c)^3 (eV^3 * cm^3)
static const double kGF             = 1.16637e-5;                    // Fermi coupling constant, (GeV^-2, natural units)
static const double kGF_SI          = kGF*khbarc3_SI*1.0e-18;        // Fermi coupling constant, (eV cm^3, SI units), 1e-18 is used for GeV^-2 -> eV^-2
static const double kNA             = 6.022142e23;                   // Avogadro number (mol^-1)

//
// ------------------------------------------------------------------------
// Enumerations
// ------------------------------------------------------------------------
//

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
   kAnaNuEAppOnly
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
   kNuMixSinsqTheta        // sin^2(theta), preferred choice in globa nue appearance + numu disappearance fit
} 
NuOscMixParam_t;

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
static const int kNSuperkSamples  = 2;  // number of samples 
static const int kns = kNSuperkSamples; // as above - shorter name
typedef enum ESuperkSample 
{
   kSKsUndefined = -1, 
   kSKsFCFV1RingMuLike, // fully-contained fiducial-volume single-ring mu-like
   kSKsFCFV1RingELike,  // fully-contained fiducial-volume single-ring e-like
} 
SuperkSample_t; 

// Enumeration of (true) reaction modes that may be considered separately within each SuperK sample
// Actual choice for current analysis is given in `Configuration' section.
static const int kNSuperkModes = 39;  // number of reaction modes kept track of within each sample
static const int knm = kNSuperkModes; // as above - shorter name
typedef enum ESuperkMode 
{
   kSKmUndefined = -1, 
   kSKmNuMuCCQE,         // numu CC QE
   kSKmNuMuCC1pi,        // numu CC 1pi
   kSKmNuMuCCcoh,        // numu CC coherent
   kSKmNuMuCCoth_mu,     // numu CC other, nu_mu disappearance analysis definition
   kSKmNuMuCCoth_e,      // numu CC other, nu_e appearance analysis definition
   kSKmNuMuNC,           // numu NC (all)
   kSKmNuMuNC1pi0,       // numu NC 1pi0
   kSKmNuMuNCcoh,        // numu NC coherent
   kSKmNuMuNCoth_e,      // numu NC other, nu_e appearance analysis definition
   kSKmNuMuBarCCQE,      // numubar CCQE
   kSKmNuMuBarCC1pi,     // numubar CC 1pi
   kSKmNuMuBarCCcoh,     // numubar CC coh
   kSKmNuMuBarCCoth_mu,  // numubar CC other, nu_mu disappearance analysis definition
   kSKmNuMuBarCCoth_e,   // numubar CC other, nu_e appearance analysis definition
   kSKmNuMuBarNC,        // numubar NC (all)
   kSKmNuMuBarNC1pi0,    // numubar NC 1pi0
   kSKmNuMuBarNCcoh,     // numubar NC coherent
   kSKmNuMuBarNCoth_e,   // numubar NC other, nu_e appearance analysis definition
   kSKmNuECCQE,          // nue CC QE
   kSKmNuECC1pi,         // nue CC 1pi
   kSKmNuECCcoh,         // nue CC coherent
   kSKmNuECCoth_mu,      // nue CC other, nu_mu disappearance analysis definition
   kSKmNuECCoth_e,       // nue CC other, nu_e appearance analysis definition
   kSKmNuENC,            // nue NC (all)
   kSKmNuENC1pi0,        // nue NC 1pi0
   kSKmNuENCcoh,         // nue NC coherent
   kSKmNuENCoth_e,       // nue NC other, nu_e appearance analysis definition
   kSKmNuEBarCCQE,       // nuebar CC QE
   kSKmNuEBarCC1pi,      // nuebar CC 1pi
   kSKmNuEBarCCcoh,      // nuebar CC coherent
   kSKmNuEBarCCoth_mu,   // nuebar CC other, nu_mu disappearance analysis definition
   kSKmNuEBarCCoth_e,    // nuebar CC other, nu_e appearance analysis definition
   kSKmNuEBarNC,         // nuebar NC (all)
   kSKmNuEBarNC1pi0,     // nuebar NC 1pi0
   kSKmNuEBarNCcoh,      // nuebar NC coherent
   kSKmNuEBarNCoth_e,    // nuebar NC other, nu_e appearance analysis definition
   kSKmSigNuECCQE,       // sig nue
   kSKmSigNuECC1pi,      // sig nue
   kSKmSigNuECCoth       // sig nue
} 
SuperkMode_t; 

// Type of SK MC template
typedef enum ESKMCTemplateType
{   
   kSKMCtUndefined = -1, 
   kSKMCtNom,   // nominal MC + tuning (ext) weights 
   kSKMCtTun,   // nominal MC + tuning (ext) weights + tuning (T2K) weights
   kSKMCtSys,   // nominal MC + tuning (ext) weights + tuning (T2K) weights + systematic variations
   kSKMCtFin,   // nominal MC + tuning (ext) weights + tuning (T2K) weights + systematic variations + physics hypothesis 
} 
SKMCTemplateType_t;

// Enumeration of flux histograms types held by NuFlux
typedef enum EFluxType
{   
   kFlxUndefined      = -1, 
   kFlxNominal        =  0,  // flux used for MC production
   kFlxExtTuning      =  1,  // flux used for MC production + tuning wghts based on external data (NA61,...)
   kFlxFinalTuning    =  2,  // flux used for MC production + tuning wghts based on external data (NA61,...) + tuning wghts based on T2K data (ND280 data fit by BANFF,...)
   kFlxCurrent        =  3   // current flux from systematic param variation of kFlxFinalTuning flux
} 
FluxType_t;

// Fit parameters
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

   /*
    * SuperK detector errors (energy scale and efficiencies)
   */   
   kNuisSuperkEnergyScale,     // SK energy scale
   kNuisSuperkEff1RmuCCQEoth,  // SK 1Rmu CCQE efficiency (other)
   kNuisSuperkEff1RmuCCQErc0u, // SK 1Rmu CCQE efficiency (ring counting, Ereco bin 1)
   kNuisSuperkEff1RmuCCQErc1u, // SK 1Rmu CCQE efficiency (ring counting, Ereco bin 2)
   kNuisSuperkEff1RmuCCQErc2u, // SK 1Rmu CCQE efficiency (ring counting, Ereco bin 3)
   kNuisSuperkEff1RmuCCnonQE,  // SK 1Rmu CCnonQE efficiency
   kNuisSuperkEff1RmuNC,       // SK 1Rmu NC efficiency
   kNuisSuperkEff1RmuNueCC,    // SK 1Rmu nue CC efficiency

   /*
    * 26 (flux and cross-section) corellated parameters used in the BANFF fit
   */   
   kNuisBANFF00,               // # 00, SK    numu    flux,  E =  0.0 -  0.4 GeV
   kNuisBANFF01,               // # 01, SK    numu    flux,  E =  0.4 -  0.5 GeV
   kNuisBANFF02,               // # 02, SK    numu    flux,  E =  0.5 -  0.6 GeV
   kNuisBANFF03,               // # 03, SK    numu    flux,  E =  0.6 -  0.7 GeV
   kNuisBANFF04,               // # 04, SK    numu    flux,  E =  0.7 -  1.0 GeV
   kNuisBANFF05,               // # 05, SK    numu    flux,  E =  1.0 -  1.5 GeV
   kNuisBANFF06,               // # 06, SK    numu    flux,  E =  1.5 -  2.5 GeV
   kNuisBANFF07,               // # 07, SK    numu    flux,  E =  2.5 -  3.5 GeV
   kNuisBANFF08,               // # 08, SK    numu    flux,  E =  3.5 -  5.0 GeV
   kNuisBANFF09,               // # 09, SK    numu    flux,  E =  5.0 -  7.0 GeV
   kNuisBANFF10,               // # 10, SK    numu    flux,  E =  7.0 - 30.0 GeV
   kNuisBANFF11,               // # 11, SK    numubar flux,  E =  0.0 -  1.5 GeV
   kNuisBANFF12,               // # 12, SK    numubar flux,  E =  1.5 - 30.0 GeV
   kNuisBANFF13,               // # 13, SK    nue     flux,  E =  0.0 -  0.5 GeV
   kNuisBANFF14,               // # 14, SK    nue     flux,  E =  0.5 -  0.7 GeV
   kNuisBANFF15,               // # 15, SK    nue     flux,  E =  0.7 -  0.8 GeV
   kNuisBANFF16,               // # 16, SK    nue     flux,  E =  0.8 -  1.5 GeV
   kNuisBANFF17,               // # 17, SK    nue     flux,  E =  1.5 -  2.5 GeV
   kNuisBANFF18,               // # 18, SK    nue     flux,  E =  2.5 -  4.0 GeV
   kNuisBANFF19,               // # 19, SK    nue     flux,  E =  4.0 - 30.0 GeV
   kNuisBANFF20,               // # 20, SK    nuebar  flux,  E =  0.0 -  2.5 GeV
   kNuisBANFF21,               // # 21, SK    nuebar  flux,  E =  2.5 - 30.0 GeV
   kNuisBANFF22,               // # 22, M_{A} QE
   kNuisBANFF23,               // # 23, M_{A} RES
   kNuisBANFF24,               // # 24, CCQE  normalization, 0 - 1.5 GeV
   kNuisBANFF25,               // # 25, CCRES normalization, 0 - 2.5 GeV

   /*
    * Cross-section parameters (excluding cross-section params which were included in BANFF fit and are correlated with the flux params)
    * These other paramers are taken to be uncorrelated and are applied to SK only.
   */   
   kNuisXSecFermiMom,          // Fermi momentum
   kNuisXSecBindE,             // binding energy
   kNuisXSecSpecFunc,          // spectral function on/off
   kNuisXSecNormCCQE1,         // CCQE norm 1.5 - 3.5 GeV
   kNuisXSecNormCCQE2,         // CCQE norm > 3.5 GeV
   kNuisXSecDeltaMass,         // Delta mean mass 
   kNuisXSecDeltaWidth,        // Delta width
   kNuisXSecNormCC1pi1,        // CC1pi norm > 2.5 GeV   
   kNuisXSecNormCCcoh,         // CC coherent norm
   kNuisXSecNormCCoth,         // CC other norm

   /*
    * FSI parameters.
    * Assuming here a 10x10 covariance matrix will be made available.
    * Each parameter is a scaling factor to apply at some Ereco range (to be defined) for some
    * SK reaction mode (or the total? to be defined)
    */
    kNuisFSI00,                // # 00, ??
    kNuisFSI01,                // # 01, ??
    kNuisFSI02,                // # 02, ??
    kNuisFSI03,                // # 03, ??
    kNuisFSI04,                // # 04, ??
    kNuisFSI05,                // # 05, ??
    kNuisFSI06,                // # 06, ??
    kNuisFSI07,                // # 07, ??
    kNuisFSI08,                // # 08, ??
    kNuisFSI09                 // # 09, ??
} 
PhysParamId_t;


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
    kSKsFCFV1RingMuLike,
    kSKsUndefined            //<-- add kSKsFCFV1RingELike for the combined fit
};

//
// Define which SK modes are considered for each SK sample included in the analysis
//

static const SuperkMode_t kSKReactionModes[kns][knm] = {

  // > 1 muon-like ring (numu disappearance)
  {
    kSKmNuMuCCQE,         
    kSKmNuMuCC1pi, 
    kSKmNuMuCCoth_mu,                      
    kSKmNuMuNC,          
    kSKmNuMuBarCCQE,      
    kSKmNuMuBarCC1pi, 
    kSKmNuMuBarCCoth_mu,                  
    kSKmNuMuBarNC,       
    kSKmNuECCQE, 
    kSKmNuECC1pi, 
    kSKmNuECCoth_mu,                                
    kSKmNuENC,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCoth_mu,
    kSKmNuEBarNC,
    kSKmSigNuECCQE, 
    kSKmSigNuECC1pi, 
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
    kSKmUndefined,
    kSKmUndefined
  },

  // > 1 electron-like ring (nue appearance)
  {
    kSKmNuMuCCQE, 
    kSKmNuMuCC1pi, 
    kSKmNuMuCCcoh, 
    kSKmNuMuCCoth_e, 
    kSKmNuMuNC1pi0, 
    kSKmNuMuNCcoh, 
    kSKmNuMuNCoth_e,                            
    kSKmNuMuBarCCQE, 
    kSKmNuMuBarCC1pi, 
    kSKmNuMuBarCCcoh, 
    kSKmNuMuBarCCoth_e,  
    kSKmNuMuBarNC1pi0, 
    kSKmNuMuBarNCcoh, 
    kSKmNuMuBarNCoth_e,
    kSKmNuECCQE, 
    kSKmNuECC1pi, 
    kSKmNuECCcoh, 
    kSKmNuECCoth_e,
    kSKmNuENC1pi0,  
    kSKmNuENCcoh, 
    kSKmNuENCoth_e,
    kSKmNuEBarCCQE,
    kSKmNuEBarCC1pi,
    kSKmNuEBarCCcoh,
    kSKmNuEBarCCoth_e,
    kSKmNuEBarNC1pi0,
    kSKmNuEBarNCcoh,
    kSKmNuEBarNCoth_e,
    kSKmSigNuECCQE,  
    kSKmSigNuECC1pi,  
    kSKmSigNuECCoth,
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

//
// ------------------------------------------------------------------------
//
// Info on systematic parameters 
//
// ------------------------------------------------------------------------
//

//
// Number of BANFF fit parameters
//
static const int kNumOfParamsInBANFFfit = 26;  

//
// ** Flux systematics
//

// Etrue energy range for flux-related BANFF parameters
static const double kEtrueRangeOfFluxParamInBANFFfit[kNumOfParamsInBANFFfit][2] =
{
   { 0.0,  0.4 },              // # 00, SK    numu    flux,  E =  0.0 -  0.4 GeV
   { 0.4,  0.5 },              // # 01, SK    numu    flux,  E =  0.4 -  0.5 GeV
   { 0.5,  0.6 },              // # 02, SK    numu    flux,  E =  0.5 -  0.6 GeV
   { 0.6,  0.7 },              // # 03, SK    numu    flux,  E =  0.6 -  0.7 GeV
   { 0.7,  1.0 },              // # 04, SK    numu    flux,  E =  0.7 -  1.0 GeV
   { 1.0,  1.5 },              // # 05, SK    numu    flux,  E =  1.0 -  1.5 GeV
   { 1.5,  2.5 },              // # 06, SK    numu    flux,  E =  1.5 -  2.5 GeV
   { 2.5,  3.5 },              // # 07, SK    numu    flux,  E =  2.5 -  3.5 GeV
   { 3.5,  5.0 },              // # 08, SK    numu    flux,  E =  3.5 -  5.0 GeV
   { 5.0,  7.0 },              // # 09, SK    numu    flux,  E =  5.0 -  7.0 GeV
   { 7.0, 30.0 },              // # 10, SK    numu    flux,  E =  7.0 - 30.0 GeV
   { 0.0,  1.5 },              // # 11, SK    numubar flux,  E =  0.0 -  1.5 GeV
   { 1.5, 30.0 },              // # 12, SK    numubar flux,  E =  1.5 - 30.0 GeV
   { 0.0,  0.5 },              // # 13, SK    nue     flux,  E =  0.0 -  0.5 GeV
   { 0.5,  0.7 },              // # 14, SK    nue     flux,  E =  0.5 -  0.7 GeV
   { 0.7,  0.8 },              // # 15, SK    nue     flux,  E =  0.7 -  0.8 GeV
   { 0.8,  1.5 },              // # 16, SK    nue     flux,  E =  0.8 -  1.5 GeV
   { 1.5,  2.5 },              // # 17, SK    nue     flux,  E =  1.5 -  2.5 GeV
   { 2.5,  4.0 },              // # 18, SK    nue     flux,  E =  2.5 -  4.0 GeV
   { 4.0, 30.0 },              // # 19, SK    nue     flux,  E =  4.0 - 30.0 GeV
   { 0.0,  2.5 },              // # 20, SK    nuebar  flux,  E =  0.0 -  2.5 GeV
   { 2.5, 30.0 },              // # 21, SK    nuebar  flux,  E =  2.5 - 30.0 GeV
   {  -1,   -1 },              // # 22, not a flux param
   {  -1,   -1 },              // # 23, not a flux param
   {  -1,   -1 },              // # 24, not a flux param
   {  -1,   -1 }               // # 25, not a flux param
};

//
// ** Cross-section systematics
//

// Nominal values of cross-section related params 
static const double kNomMAQE  = 1.21;
static const double kNomMARES = 1.21;

// Number of cross-section related systematics 
static const int kNumOfXSecSyst = 4;

// The following is the list of the above kNumOfXSecSyst parameters
static const PhysParamId_t kListOfXSecSyst[kNumOfXSecSyst] = 
{
   kNuisBANFF22,               // # 22, M_{A} QE
   kNuisBANFF23,               // # 23, M_{A} RES
   kNuisBANFF24,               // # 24, CCQE  normalization, 0 - 1.5 GeV
   kNuisBANFF25                // # 25, CCRES normalization, 0 - 2.5 GeV
};

// The following list specifies which cross-section related systematic has its
// effect on MC templates pre-computed and stored in cubic splines.
// The order of systematic parameters should match the one in kListOfXSecSyst.
static const bool kXSecSystHasWghtSplines[kNumOfXSecSyst] = 
{
   true,                      // kNuisBANFF22
   true,                      // kNuisBANFF23
   false,                     // kNuisBANFF24
   false                      // kNuisBANFF25
};

// Ereco and Etrue cuts for the application of cross-section systematics.
// The order of systematic parameters should match the one in kListOfXSecSyst.
static const double kXSecSystErecoCut[kNumOfXSecSyst][2] =
{
   { 0.0, 100.0 },            // kNuisBANFF22
   { 0.0, 100.0 },            // kNuisBANFF23
   { 0.0, 100.0 },            // kNuisBANFF24
   { 0.0, 100.0 }             // kNuisBANFF25
};
static const double kXSecSystEtrueCut[kNumOfXSecSyst][2] =
{
   { 0.0, 100.0 },            // kNuisBANFF22
   { 0.0, 100.0 },            // kNuisBANFF23
   { 0.0,   1.5 },            // kNuisBANFF24
   { 0.0,   2.5 }             // kNuisBANFF25
};

// The following list describles in which MC templates we are supposed to apply the
// weights as a result of tweaking any of the above systematic parameters.
// The order of systematic parameters should match the one in kListOfXSecSyst.
static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {

  // > 1 muon-like ring 
  {
    // banff22  banff23  banff24  banff25
    {  true,    false,   true,    false }, // kSKmNuMuCCQE         
    {  false,   true,    false,   true  }, // kSKmNuMuCC1pi 
    {  false,   false,   false,   false }, // kSKmNuMuCCoth_mu                      
    {  false,   false,   false,   false }, // kSKmNuMuNC          
    {  true,    false,   true,    false }, // kSKmNuMuBarCCQE      
    {  false,   true,    false,   true  }, // kSKmNuMuBarCC1pi 
    {  false,   false,   false,   false }, // kSKmNuMuBarCCoth_mu                  
    {  false,   false,   false,   false }, // kSKmNuMuBarNC       
    {  true,    false,   true,    false }, // kSKmNuECCQE 
    {  false,   true,    false,   true  }, // kSKmNuECC1pi 
    {  false,   false,   false,   false }, // kSKmNuECCoth_mu   
    {  false,   false,   false,   false }, // kSKmNuENC
    {  true,    false,   true,    false }, // kSKmNuEBarCCQE
    {  false,   true,    false,   true  }, // kSKmNuEBarCC1pi
    {  false,   false,   false,   false }, // kSKmNuEBarCCoth_mu
    {  false,   false,   false,   false }, // kSKmNuEBarNC
    {  true,    false,   true,    false }, // kSKmSigNuECCQE 
    {  false,   true,    false,   true  }, // kSKmSigNuECC1pi 
    {  false,   false,   false,   false }, // kSKmSigNuECCoth
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }  // kSKmUndefined
  },

  // > 1 electron-like ring 
  {
    // banff22  banff23  banff24  banff25
    {  false,   false,   false,   false }, // kSKmNuMuCCQE 
    {  false,   false,   false,   false }, // kSKmNuMuCC1pi 
    {  false,   false,   false,   false }, // kSKmNuMuCCcoh 
    {  false,   false,   false,   false }, // kSKmNuMuCCoth_e 
    {  false,   false,   false,   false }, // kSKmNuMuNC1pi0 
    {  false,   false,   false,   false }, // kSKmNuMuNCcoh 
    {  false,   false,   false,   false }, // kSKmNuMuNCoth_e   
    {  false,   false,   false,   false }, // kSKmNuMuBarCCQE 
    {  false,   false,   false,   false }, // kSKmNuMuBarCC1pi 
    {  false,   false,   false,   false }, // kSKmNuMuBarCCcoh 
    {  false,   false,   false,   false }, // kSKmNuMuBarCCoth_e  
    {  false,   false,   false,   false }, // kSKmNuMuBarNC1pi0 
    {  false,   false,   false,   false }, // kSKmNuMuBarNCcoh 
    {  false,   false,   false,   false }, // kSKmNuMuBarNCoth_e
    {  false,   false,   false,   false }, // kSKmNuECCQE 
    {  false,   false,   false,   false }, // kSKmNuECC1pi 
    {  false,   false,   false,   false }, // kSKmNuECCcoh 
    {  false,   false,   false,   false }, // kSKmNuECCoth_e
    {  false,   false,   false,   false }, // kSKmNuENC1pi0  
    {  false,   false,   false,   false }, // kSKmNuENCcoh 
    {  false,   false,   false,   false }, // kSKmNuENCoth_e
    {  false,   false,   false,   false }, // kSKmNuEBarCCQE
    {  false,   false,   false,   false }, // kSKmNuEBarCC1pi
    {  false,   false,   false,   false }, // kSKmNuEBarCCcoh
    {  false,   false,   false,   false }, // kSKmNuEBarCCoth_e
    {  false,   false,   false,   false }, // kSKmNuEBarNC1pi0
    {  false,   false,   false,   false }, // kSKmNuEBarNCcoh
    {  false,   false,   false,   false }, // kSKmNuEBarNCoth_e
    {  false,   false,   false,   false }, // kSKmSigNuECCQE  
    {  false,   false,   false,   false }, // kSKmSigNuECC1pi  
    {  false,   false,   false,   false }, // kSKmSigNuECCoth
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }, // kSKmUndefined
    {  false,   false,   false,   false }  // kSKmUndefined
  }
};

//
// The weight splines are pre-computed for the following Ereco bins.
// For the time being, the same Ereco binning is used for all samples and all modes.
//

static const int kXsNErecoBins = 28;

static const double kXsErecoBinEdges[kXsNErecoBins+1] = 
{
   0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450, 
   0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950, 
   1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  
   3.000, 10.000, 30.000
};

//
// The weight splines are pre-computed in the following Etrue bins (for each Ereco bin)
// The Etrue bins can, in principle, be different for each sample/mode combination.
// For the time-being, the Etrue bins are the BANFF ones (depend on the neutrino code).
// Note that the order of modes below should match the one defined in kSKReactionModes array.
//

static const int kXsNEtrueBins = 11;

static const double kXsEtrueBinEdges[kns][knm][kXsNEtrueBins+1] = {

  // > 1 muon-like ring 
  {
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCCQE         
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCC1pi 
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCCoth_mu                      
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuNC          
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCCQE      
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCC1pi 
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCCoth_mu                  
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarNC       
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECCQE 
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECC1pi 
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECCoth_mu                                
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuENC
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCCQE
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCC1pi
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCCoth_mu
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarNC
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmSigNuECCQE 
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmSigNuECC1pi 
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmSigNuECCoth
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }  // kSKmUndefined
  },

  // > 1 electron-like ring 
  {
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCCQE 
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCC1pi 
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCCcoh 
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuCCoth_e 
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuNC1pi0 
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuNCcoh 
    {  0.000,  0.400,  0.500,  0.600,  0.700,  1.000,  1.500,  2.500,  3.500,  5.000,  7.000, 30.000 }, // kSKmNuMuNCoth_e                            
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCCQE 
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCC1pi 
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCCcoh 
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarCCoth_e  
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarNC1pi0 
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarNCcoh 
    {  0.000,  1.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuMuBarNCoth_e
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECCQE 
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECC1pi 
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECCcoh 
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuECCoth_e
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuENC1pi0  
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuENCcoh 
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuENCoth_e
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCCQE
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCC1pi
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCCcoh
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarCCoth_e
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarNC1pi0
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarNCcoh
    {  0.000,  2.500, 30.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmNuEBarNCoth_e
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmSigNuECCQE  
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmSigNuECC1pi  
    {  0.000,  0.500,  0.700,  0.800,  1.500,  2.500,  4.000, 30.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmSigNuECCoth
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }, // kSKmUndefined
    { -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000, -1.000 }  // kSKmUndefined
  }
};

// Number of cross-section weight splines
static const int kNumOfXSecSystWghtSpl = kNSuperkSamples * kNSuperkModes * kNumOfXSecSyst * kXsNErecoBins * kXsNEtrueBins;
static const int knwspl = kNumOfXSecSystWghtSpl; // as above - shorter name

//
// ** FSI systematics
//

static const int kNumOfFSIParams = 10;


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
    { 3244,  kBlack,    1,  1,  kBlack}, // kSKmNuMuCCoth_mu
    { 3006,  kRed,      1,  1,  kBlack}, // kSKmNuMuNC
    { 3001,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCQE
    { 3345,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCC1pi
    { 3354,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarCCoth_mu
    { 3006,  kAzure-3,  1,  1,  kBlack}, // kSKmNuMuBarNC
    { 3001,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCQE
    { 3345,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECC1pi
    { 3354,  kGreen+1,  1,  1,  kBlack}, // kSKmNuECCoth_mu
    { 3006,  kGreen+1,  1,  1,  kBlack}, // kSKmNuENC
    { 3001,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCQE
    { 3345,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCC1pi
    { 3354,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarCCoth_mu
    { 3006,  kBlue,     1,  1,  kBlack}, // kSKmNuEBarNC
    { 3001,  kGreen-3,  1,  1,  kBlack}, // kSKmSigNuECCQE
    { 3345,  kGreen-3,  1,  1,  kBlack}, // kSKmSigNuECC1pi
    { 3354,  kGreen-3,  1,  1,  kBlack}, // kSKmSigNuECCoth
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
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
  },

  // > 1 electron-like ring 
  {
    {   -1,    kWhite,  1,  1, kBlack }, // kSKmNuMuCCQE
    { 3544,    kBlack,  1,  1, kBlack }, // kSKmNuMuCC1pi
    { 3244,    kGreen,  1,  1, kBlack }, // kSKmNuMuCCcoh
    { 3244,     kBlue,  1,  1, kBlack }, // kSKmNuMuCCoth_e
    { 3006,      kRed,  1,  1, kBlack }, // kSKmNuMuNC1pi0
    { 3006,    kRed-1,  1,  1, kBlack }, // kSKmNuMuNCcoh
    { 3006,    kRed-2,  1,  1, kBlack }, // kSKmNuMuNCoth_e
    { 3001,  kAzure-3,  1,  1, kBlack }, // kSKmNuMuBarCCQE
    { 3345,  kAzure-4,  1,  1, kBlack }, // kSKmNuMuBarCC1pi
    { 3354,  kAzure-6,  1,  1, kBlack }, // kSKmNuMuBarCCcoh
    { 3354,  kMagenta,  1,  1, kBlack }, // kSKmNuMuBarCCoth_e
    { 3006,  kAzure-7,  1,  1, kBlack }, // kSKmNuMuBarNC1pi0
    { 3006,  kAzure-9,  1,  1, kBlack }, // kSKmNuMuBarNCcoh
    { 3006,  kMagenta,  1,  1, kBlack }, // kSKmNuMuBarNCoth_e
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
    { 3001,  kGreen+8,  1,  1, kBlack }, // kSKmSigNuECCQE
    { 3345, kGreen+10,  1,  1, kBlack }, // kSKmSigNuECC1pi
    { 3354,      kRed,  1,  1, kBlack }, // kSKmSigNuECCoth
    {   -1,        -1, -1, -1,     -1 }, // kSKmUndefined
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
    kSKmSigNuECC1pi, 
    kSKmSigNuECCoth,
    kSKmSigNuECCQE, 
    kSKmNuEBarNC,                                                               
    kSKmNuEBarCC1pi, 
    kSKmNuEBarCCoth_mu,                                
    kSKmNuEBarCCQE, 
    kSKmNuENC,                                                               
    kSKmNuECC1pi, 
    kSKmNuECCoth_mu,                                
    kSKmNuECCQE, 
    kSKmNuMuBarNC,       
    kSKmNuMuBarCC1pi, 
    kSKmNuMuBarCCoth_mu,                  
    kSKmNuMuBarCCQE,      
    kSKmNuMuNC,          
    kSKmNuMuCCoth_mu,                      
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
    kSKmUndefined
  },

  // > 1 electron-like ring (nue appearance)
  {
    kSKmNuMuCCQE, 
    kSKmNuMuCC1pi, 
    kSKmNuMuCCcoh, 
    kSKmNuMuCCoth_e, 
    kSKmNuMuNC1pi0, 
    kSKmNuMuNCcoh, 
    kSKmNuMuNCoth_e,                            
    kSKmNuMuBarCCQE, 
    kSKmNuMuBarCC1pi, 
    kSKmNuMuBarCCcoh, 
    kSKmNuMuBarCCoth_e,  
    kSKmNuMuBarNC1pi0, 
    kSKmNuMuBarNCcoh, 
    kSKmNuMuBarNCoth_e,
    kSKmNuEBarCCQE, 
    kSKmNuEBarCC1pi, 
    kSKmNuEBarCCcoh, 
    kSKmNuEBarCCoth_e,
    kSKmNuEBarNC1pi0,  
    kSKmNuEBarNCcoh, 
    kSKmNuEBarNCoth_e,
    kSKmNuECCQE, 
    kSKmNuECC1pi, 
    kSKmNuECCcoh, 
    kSKmNuECCoth_e,
    kSKmNuENC1pi0,  
    kSKmNuENCcoh, 
    kSKmNuENCoth_e,
    kSKmSigNuECCQE,  
    kSKmSigNuECC1pi,  
    kSKmSigNuECCoth,
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

//
// ------------------------------------------------------------------------
// Binning 
// ------------------------------------------------------------------------
//

//
// Etrue binning for SK MC templates
//

static const double kEtrueMin    = 0;
static const double kEtrueMax    = 12;
static const double kEtrueBinSz  = 0.05;
static const int    kNEtrueBins  = (int) ((kEtrueMax - kEtrueMin)/kEtrueBinSz);

//
// Ereco binning for SK 1Rmu sample
//

// 200 50-MeV bins from 0 to 10 GeV
static const int kNErecoBins1Rmu = 200;
static const double kErecoBinEdges1Rmu[kNErecoBins1Rmu+1] = 
{
   0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450, 
   0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950, 
   1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  1.300,  1.350,  1.400,  1.450, 
   1.500,  1.550,  1.600,  1.650,  1.700,  1.750,  1.800,  1.850,  1.900,  1.950, 
   2.000,  2.050,  2.100,  2.150,  2.200,  2.250,  2.300,  2.350,  2.400,  2.450, 
   2.500,  2.550,  2.600,  2.650,  2.700,  2.750,  2.800,  2.850,  2.900,  2.950, 
   3.000,  3.050,  3.100,  3.150,  3.200,  3.250,  3.300,  3.350,  3.400,  3.450, 
   3.500,  3.550,  3.600,  3.650,  3.700,  3.750,  3.800,  3.850,  3.900,  3.950, 
   4.000,  4.050,  4.100,  4.150,  4.200,  4.250,  4.300,  4.350,  4.400,  4.450, 
   4.500,  4.550,  4.600,  4.650,  4.700,  4.750,  4.800,  4.850,  4.900,  4.950, 
   5.000,  5.050,  5.100,  5.150,  5.200,  5.250,  5.300,  5.350,  5.400,  5.450, 
   5.500,  5.550,  5.600,  5.650,  5.700,  5.750,  5.800,  5.850,  5.900,  5.950, 
   6.000,  6.050,  6.100,  6.150,  6.200,  6.250,  6.300,  6.350,  6.400,  6.450, 
   6.500,  6.550,  6.600,  6.650,  6.700,  6.750,  6.800,  6.850,  6.900,  6.950, 
   7.000,  7.050,  7.100,  7.150,  7.200,  7.250,  7.300,  7.350,  7.400,  7.450, 
   7.500,  7.550,  7.600,  7.650,  7.700,  7.750,  7.800,  7.850,  7.900,  7.950, 
   8.000,  8.050,  8.100,  8.150,  8.200,  8.250,  8.300,  8.350,  8.400,  8.450, 
   8.500,  8.550,  8.600,  8.650,  8.700,  8.750,  8.800,  8.850,  8.900,  8.950, 
   9.000,  9.050,  9.100,  9.150,  9.200,  9.250,  9.300,  9.350,  9.400,  9.450, 
   9.500,  9.550,  9.600,  9.650,  9.700,  9.750,  9.800,  9.850,  9.900,  9.950, 
  10.000
};
/*
// Fewer bins at higher energies for improved speed
// -   50 MeV bins in 0-3 GeV
// -  200 MeV bins in 3-4 GeV
// -  500 MeV bins in 4-6 GeV
// - 1000 MeV bins > 6 GeV
static const int kNErecoBins1Rmu = 73;
static const double kErecoBinEdges1Rmu[kNErecoBins1Rmu+1] = 
{
   0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450, 
   0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950, 
   1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  1.300,  1.350,  1.400,  1.450, 
   1.500,  1.550,  1.600,  1.650,  1.700,  1.750,  1.800,  1.850,  1.900,  1.950, 
   2.000,  2.050,  2.100,  2.150,  2.200,  2.250,  2.300,  2.350,  2.400,  2.450, 
   2.500,  2.550,  2.600,  2.650,  2.700,  2.750,  2.800,  2.850,  2.900,  2.950, 
   3.000,  3.200,  3.400,  3.600,  3.800,
   4.000,  4.500,  5.000,  5.500,  
   6.000,  7.000,  8.000,  9.000, 10.000
};
*/

//
// Ereco binning for SK 1Re sample
//
// 50 bins from 0 to 10 GeV (same as the 2010a NIWG input for SKFSI)
static const int kNErecoBins1Re = 50;
static const double kErecoBinEdges1Re[kNErecoBins1Re+1] = 
{
  0.000,  0.200,  0.400, 0.600, 0.800,  1.000,  1.200,  1.400, 1.600, 1.800,
  2.000,  2.200,  2.400, 2.600, 2.800,  3.000,  3.200,  3.400, 3.600, 3.800,
  4.000,  4.200,  4.400, 4.600, 4.800,  5.000,  5.200,  5.400, 5.600, 5.800,
  6.000,  6.200,  6.400, 6.600, 6.800,  7.000,  7.200,  7.400, 7.600, 7.800,
  8.000,  8.200,  8.400, 8.600, 8.800,  9.000,  9.200,  9.400, 9.600, 9.800,
  10.000
};
/*
// 1 bin from 0 to 10 GeV (counting expt)
//WARNING this gives wrong errors for anything that is Ereco dependant (SKFSI)
static const int kNErecoBins1Re = 1;
static const double kErecoBinEdges1Re[kNErecoBins1Re+1] = 
{ 
   0., 10. 
};
*/
/*
// Same binning as the standard 1Rmu one: 200 50-MeV bins from 0 to 10 GeV
static const int kNErecoBins1Re = 200;
static const double kErecoBinEdges1Re[kNErecoBins1Re+1] = 
{
   0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450, 
   0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950, 
   1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  1.300,  1.350,  1.400,  1.450, 
   1.500,  1.550,  1.600,  1.650,  1.700,  1.750,  1.800,  1.850,  1.900,  1.950, 
   2.000,  2.050,  2.100,  2.150,  2.200,  2.250,  2.300,  2.350,  2.400,  2.450, 
   2.500,  2.550,  2.600,  2.650,  2.700,  2.750,  2.800,  2.850,  2.900,  2.950, 
   3.000,  3.050,  3.100,  3.150,  3.200,  3.250,  3.300,  3.350,  3.400,  3.450, 
   3.500,  3.550,  3.600,  3.650,  3.700,  3.750,  3.800,  3.850,  3.900,  3.950, 
   4.000,  4.050,  4.100,  4.150,  4.200,  4.250,  4.300,  4.350,  4.400,  4.450, 
   4.500,  4.550,  4.600,  4.650,  4.700,  4.750,  4.800,  4.850,  4.900,  4.950, 
   5.000,  5.050,  5.100,  5.150,  5.200,  5.250,  5.300,  5.350,  5.400,  5.450, 
   5.500,  5.550,  5.600,  5.650,  5.700,  5.750,  5.800,  5.850,  5.900,  5.950, 
   6.000,  6.050,  6.100,  6.150,  6.200,  6.250,  6.300,  6.350,  6.400,  6.450, 
   6.500,  6.550,  6.600,  6.650,  6.700,  6.750,  6.800,  6.850,  6.900,  6.950, 
   7.000,  7.050,  7.100,  7.150,  7.200,  7.250,  7.300,  7.350,  7.400,  7.450, 
   7.500,  7.550,  7.600,  7.650,  7.700,  7.750,  7.800,  7.850,  7.900,  7.950, 
   8.000,  8.050,  8.100,  8.150,  8.200,  8.250,  8.300,  8.350,  8.400,  8.450, 
   8.500,  8.550,  8.600,  8.650,  8.700,  8.750,  8.800,  8.850,  8.900,  8.950, 
   9.000,  9.050,  9.100,  9.150,  9.200,  9.250,  9.300,  9.350,  9.400,  9.450, 
   9.500,  9.550,  9.600,  9.650,  9.700,  9.750,  9.800,  9.850,  9.900,  9.950, 
  10.000
};
*/


} // conv
} // t2kana
  
#endif 


