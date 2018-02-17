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
 
 22-12-2013  Nick
 Added extra arrays for reverse horn current systematic parameters
 */
//____________________________________________________________________________

#ifndef banff_syst_conv_h
#define banff_syst_conv_h

#include "AnalysisOptions.h"
#include "Enumerations.h"
#include "SamplesAndModes.h"

namespace valor {
  namespace conv   {
    
    static const int kNBANFFluxPDGs = 4;
    static const int kOrderOfBANFFFluxPDG[kNBANFFluxPDGs] = { kPdgNuMu, kPdgNuMuBar, kPdgNuE, kPdgNuEBar };
      
    //
    // Number of BANFF fit parameters
    //
        
    // NEW XSEC NIWG
#if defined(__USE_NIWG2015__)
#  if BANFF_BINS == 1 //joint
#  elif BANFF_BINS == 2 //disap (no scca and sccv)
#  elif BANFF_BINS == 3 //app
#  endif
#endif
#if defined(__USE_NIWG2017__)
#  if BANFF_BINS == 1 //joint
    static const int kNumOfParamsInBANFFfit = 56;
    static const int kNumOfXSecParamsInBANFFfit = 31;
    static const int kNumOfFluxParamsInBANFFfit = kNumOfParamsInBANFFfit - kNumOfXSecParamsInBANFFfit;
    static const int kNumOfRHCParamsInBANFFfit = 56;
    static const int kNumOfRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - kNumOfXSecParamsInBANFFfit;
    // NEW FLUX FHC+RHC
    static const int kNumOfFHCpRHCParamsInBANFFfit = 81;
    static const int kNumOfFHCpRHCFluxParamsInBANFFfit = kNumOfFHCpRHCParamsInBANFFfit - kNumOfXSecParamsInBANFFfit ;     
#  endif
#endif
        
        
    static const PhysParamId_t kListOfBANFFParams[kNumOfParamsInBANFFfit] =
    {
        // NEW XSEC NIWG
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )
      // NEW XSEC
      kNuisBANFF00 ,
      kNuisBANFF01 ,
      kNuisBANFF02 ,
      kNuisBANFF03 ,
      kNuisBANFF04 ,
      kNuisBANFF05 ,
      kNuisBANFF06 ,
      kNuisBANFF07 ,
      kNuisBANFF08 ,
      kNuisBANFF09 ,
      kNuisBANFF10 ,
      kNuisBANFF11 ,
      kNuisBANFF12 ,
      kNuisBANFF13 ,
      kNuisBANFF14 ,
      kNuisBANFF15 ,
      kNuisBANFF16 ,
      kNuisBANFF17 ,
      kNuisBANFF18 ,
      kNuisBANFF19 ,
      kNuisBANFF20 ,
      kNuisBANFF21 ,
      kNuisBANFF22 ,
      kNuisBANFF23 ,
      kNuisBANFF24 ,
            
      kNuisBANFFNorm2p2h     , //
      kNuisBANFFCA5          , //
      kNuisBANFFBgRES        , //
      kNuisBANFFMAQE         , //
      kNuisBANFFMARES        , //
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      kNuisBANFFSCCA         ,
      kNuisBANFFSCCV         ,
#endif
      //
      kNuisBANFFFermiMom      , //
      kNuisBANFFShapeCCoth    , //
      kNuisBANFFNormCCcoh     , //
      kNuisBANFFNormNCcoh     , //
      kNuisBANFFNormNCoth     ,  //
      kNuisBANFFNormNuEToNuMu ,  //
      kNuisBANFFNormNC1Gamma      , // NC 1 gamma
      kNuisBANFFNormNuEBarToNuMuBar , //
      kNuisBANFFNorm2p2hBar          , //
#  if defined(__USE_NIWG2017__)      
      kNuisBANFFBeRPA_A     , //
      kNuisBANFFBeRPA_B     , //
      kNuisBANFFBeRPA_D     , //
      kNuisBANFFBeRPA_E     , //
      kNuisBANFFBeRPA_U     , //
      kNuisBANFFShape2p2h   , //
      kNuisBANFFNorm2p2hCtoO, //
      kNuisBANFFNormCCNu,
      kNuisBANFFNormCCNuBar,
      kNuisBANFFSKFSIAbsorption,
      kNuisBANFFSKFSIQELow,
      kNuisBANFFSKFSICXLow,
      kNuisBANFFSKFSIQEHigh,
      kNuisBANFFSKFSIInelastic,
      kNuisBANFFEb
#  endif      
#endif
    };
        
    static const PhysParamId_t kListOfRHCBANFFParams[kNumOfRHCParamsInBANFFfit] =
    {            
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )     
      // NEW XSEC
      kRHCNuisBANFF00 ,
      kRHCNuisBANFF01 ,
      kRHCNuisBANFF02 ,
      kRHCNuisBANFF03 ,
      kRHCNuisBANFF04 ,
      kRHCNuisBANFF05 ,
      kRHCNuisBANFF06 ,
      kRHCNuisBANFF07 ,
      kRHCNuisBANFF08 ,
      kRHCNuisBANFF09 ,
      kRHCNuisBANFF10 ,
      kRHCNuisBANFF11 ,
      kRHCNuisBANFF12 ,
      kRHCNuisBANFF13 ,
      kRHCNuisBANFF14 ,
      kRHCNuisBANFF15 ,
      kRHCNuisBANFF16 ,
      kRHCNuisBANFF17 ,
      kRHCNuisBANFF18 ,
      kRHCNuisBANFF19 ,
      kRHCNuisBANFF20 ,
      kRHCNuisBANFF21 ,
      kRHCNuisBANFF22 ,
      kRHCNuisBANFF23 ,
      kRHCNuisBANFF24 ,
            
      kRHCNuisBANFF2p2h   ,
      kRHCNuisBANFFCA5   ,
      kRHCNuisBANFFBgRES ,
      kRHCNuisBANFFMAQE   ,
      kRHCNuisBANFFMARES  ,
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      kRHCNuisBANFFSCCA   ,
      kRHCNuisBANFFSCCV   ,
#endif
      //
      kRHCNuisBANFFFermiMom      ,
      kRHCNuisBANFFShapeCCoth    ,
      kRHCNuisBANFFNormCCcoh     ,
      kRHCNuisBANFFNormNCcoh     ,
      kRHCNuisBANFFNormNCoth     ,
      kRHCNuisBANFFNormNuEToNuMu ,
      kRHCNuisBANFFNormNC1Gamma      , // NC 1 gamma
      kRHCNuisBANFFNormNuEBarToNuMuBar ,
      kRHCNuisBANFF2p2hBar   ,
#  if defined(__USE_NIWG2017__)      
      kRHCNuisBANFFBeRPA_A ,
      kRHCNuisBANFFBeRPA_B     ,
      kRHCNuisBANFFBeRPA_D     ,
      kRHCNuisBANFFBeRPA_E     ,
      kRHCNuisBANFFBeRPA_U     ,
      kRHCNuisBANFFShape2p2h   ,
      kRHCNuisBANFFNorm2p2hCtoO,
      kRHCNuisBANFFNormCCNu,
      kRHCNuisBANFFNormCCNuBar,
      kRHCNuisBANFFSKFSIAbsorption,
      kRHCNuisBANFFSKFSIQELow,
      kRHCNuisBANFFSKFSICXLow,
      kRHCNuisBANFFSKFSIQEHigh,
      kRHCNuisBANFFSKFSIInelastic,
      kRHCNuisBANFFEb
#  endif
      
#endif
    };
        
    static const PhysParamId_t kListOfFHCpRHCBANFFParams[kNumOfFHCpRHCParamsInBANFFfit] =
    {            
#if ( defined(__USE_NIWG2015__)  || defined(__USE_NIWG2017__) )
      // NEW XSEC
      kNuisBANFF00 ,
      kNuisBANFF01 ,
      kNuisBANFF02 ,
      kNuisBANFF03 ,
      kNuisBANFF04 ,
      kNuisBANFF05 ,
      kNuisBANFF06 ,
      kNuisBANFF07 ,
      kNuisBANFF08 ,
      kNuisBANFF09 ,
      kNuisBANFF10 ,
      kNuisBANFF11 ,
      kNuisBANFF12 ,
      kNuisBANFF13 ,
      kNuisBANFF14 ,
      kNuisBANFF15 ,
      kNuisBANFF16 ,
      kNuisBANFF17 ,
      kNuisBANFF18 ,
      kNuisBANFF19 ,
      kNuisBANFF20 ,
      kNuisBANFF21 ,
      kNuisBANFF22 ,
      kNuisBANFF23 ,
      kNuisBANFF24 ,
            
      //
      kRHCNuisBANFF00 ,
      kRHCNuisBANFF01 ,
      kRHCNuisBANFF02 ,
      kRHCNuisBANFF03 ,
      kRHCNuisBANFF04 ,
      kRHCNuisBANFF05 ,
      kRHCNuisBANFF06 ,
      kRHCNuisBANFF07 ,
      kRHCNuisBANFF08 ,
      kRHCNuisBANFF09 ,
      kRHCNuisBANFF10 ,
      kRHCNuisBANFF11 ,
      kRHCNuisBANFF12 ,
      kRHCNuisBANFF13 ,
      kRHCNuisBANFF14 ,
      kRHCNuisBANFF15 ,
      kRHCNuisBANFF16 ,
      kRHCNuisBANFF17 ,
      kRHCNuisBANFF18 ,
      kRHCNuisBANFF19 ,
      kRHCNuisBANFF20 ,
      kRHCNuisBANFF21 ,
      kRHCNuisBANFF22 ,
      kRHCNuisBANFF23 ,
      kRHCNuisBANFF24 ,
            
      //
      kNuisBANFFNorm2p2h,
      kNuisBANFFCA5,
      kNuisBANFFBgRES,
      kNuisBANFFMAQE   ,
      kNuisBANFFMARES  ,
      // NEW SPLIT ANA
#  if (BANFF_BINS != 2)
      kNuisBANFFSCCA   ,
      kNuisBANFFSCCV   ,
#  endif
      kNuisBANFFFermiMom      ,
      kNuisBANFFShapeCCoth    ,
      kNuisBANFFNormCCcoh     ,
      kNuisBANFFNormNCcoh     ,
      kNuisBANFFNormNCoth     ,
      kNuisBANFFNormNuEToNuMu ,
      kNuisBANFFNormNC1Gamma      , // NC 1 gamma
      kNuisBANFFNormNuEBarToNuMuBar ,
      kNuisBANFFNorm2p2hBar,
#  if defined(__USE_NIWG2017__)      
      kNuisBANFFBeRPA_A ,
      kNuisBANFFBeRPA_B     ,
      kNuisBANFFBeRPA_D     ,
      kNuisBANFFBeRPA_E     ,
      kNuisBANFFBeRPA_U     ,
      kNuisBANFFShape2p2h   ,
      kNuisBANFFNorm2p2hCtoO,
      kNuisBANFFNormCCNu,
      kNuisBANFFNormCCNuBar,
      kNuisBANFFSKFSIAbsorption,
      kNuisBANFFSKFSIQELow,
      kNuisBANFFSKFSICXLow,
      kNuisBANFFSKFSIQEHigh,
      kNuisBANFFSKFSIInelastic,
      kNuisBANFFEb
#  endif
#endif
    };
        
        
    //static const PhysParamId_t kFirstBANFFXsecParam        = kListOfBANFFParams       [kNumOfFluxParamsInBANFFfit];
    //static const PhysParamId_t kFirstRHCBANFFXsecParam     = kListOfRHCBANFFParams    [kNumOfRHCFluxParamsInBANFFfit];
    //static const PhysParamId_t kFirstFHCpRHCBANFFXsecParam = kListOfFHCpRHCBANFFParams[kNumOfFHCpRHCFluxParamsInBANFFfit];
        
    //
    // ** Flux systematics
    //
        
    // Mark which BANFF parameters are SuperK flux params
    static const bool kParamInBANFFfitIsSuperkFluxParam[kNumOfParamsInBANFFfit] =
    {
        // NEW XSEC NIWG            
#if ( defined(__USE_NIWG2015__)  || defined(__USE_NIWG2017__) )
      // NEW XSEC
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
      true,                 // kNuisBANFF16
      true,                 // kNuisBANFF17
      true,                 // kNuisBANFF18
      true,                 // kNuisBANFF19
      true,                 // kNuisBANFF20
      true,                 // kNuisBANFF21
      true,                 // kNuisBANFF22
      true,                 // kNuisBANFF23
      true,                 // kNuisBANFF24
      //
      false,                // kNuisBANFFNorm2p2h
      false,                // kNuisBANFFCA5
      false,                // kNuisBANFFBgRES
      false,                // kNuisBANFFMAQE
      false,                // kNuisBANFFMARES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      false,                // kNuisBANFFSCCA
      false,                // kNuisBANFFSCCV
#endif
      false,                // kNuisBANFFFermiMom
      false,                // kNuisBANFFShapeCCoth
      false,                // kNuisBANFFNormCCcoh
      false,                // kNuisBANFFNormNCcoh
      false,                // kNuisBANFFNormNCoth
      false,                // kNuisBANFFNormNuEToNuMu
      false,                // kNuisBANFFNormNC1Gamma
      false,                // kNuisBANFFNormNuEBarToNuMuBar
      false,                // kNuisBANFFNorm2p2hBar
#  if defined(__USE_NIWG2017__)      
      false,                // kNuisBANFFBeRPA_A     
      false,                // kNuisBANFFBeRPA_B     
      false,                // kNuisBANFFBeRPA_D     
      false,                // kNuisBANFFBeRPA_E     
      false,                // kNuisBANFFBeRPA_U     
      false,                // kNuisBANFFShape2p2h   
      false,                // kNuisBANFFNorm2p2hCtoO
      false,                // kNuisBANFFNormCCNu,
      false,                // kNuisBANFFNormCCNuBar,
      false,                // kNuisBANFFSKFSIAbsorption,
      false,                // kNuisBANFFSKFSIQELow,
      false,                // kNuisBANFFSKFSICXLow,
      false,                // kNuisBANFFSKFSIQEHigh,
      false,                // kNuisBANFFSKFSIInelastic
      false                 // kNuisBANFFEb
#  endif

#endif
    };
        
    static const bool kRHCParamInBANFFfitIsSuperkFluxParam[kNumOfRHCParamsInBANFFfit] =
    {
        // NEW XSEC NIWG            
#if ( defined(__USE_NIWG2015__)  || defined(__USE_NIWG2017__) )
      // NEW XSEC
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
      true,                 // kRHCNuisBANFF16
      true,                 // kRHCNuisBANFF17
      true,                 // kRHCNuisBANFF18
      true,                 // kRHCNuisBANFF19
      true,                 // kRHCNuisBANFF20
      true,                 // kRHCNuisBANFF21
      true,                 // kRHCNuisBANFF22
      true,                 // kRHCNuisBANFF23
      true,                 // kRHCNuisBANFF24
      //
      false,                // kRHCNuisBANFF2p2h
      false,                // kRHCNuisBANFFCA5
      false,                // kRHCNuisBANFFBgRES
      false,                // kRHCNuisBANFFMAQE
      false,                // kRHCNuisBANFFMARES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      false,                // kRHCNuisBANFFSCCA
      false,                // kRHCNuisBANFFSCCV
#endif
      false,                // kRHCNuisBANFFFermiMom
      false,                // kRHCNuisBANFFShapeCCoth
      false,                // kRHCNuisBANFFNormCCcoh
      false,                // kRHCNuisBANFFNormNCcoh
      false,                // kRHCNuisBANFFNormNCoth
      false,                // kRHCNuisBANFFNormNuEToNuMu
      false,                // kNuisBANFFNormNC1Gamma
      false,                // kRHCNuisBANFFNormNuEBarToNuMuBar
      false,                // kRHCNuisBANFF2p2hBar
#  if defined(__USE_NIWG2017__)      
      false,                // kRHCNuisBANFFBeRPA_A     
      false,                // kRHCNuisBANFFBeRPA_B     
      false,                // kRHCNuisBANFFBeRPA_D     
      false,                // kRHCNuisBANFFBeRPA_E     
      false,                // kRHCNuisBANFFBeRPA_U     
      false,                // kRHCNuisBANFFShape2p2h   
      false,                // kRHCNuisBANFFNorm2p2hCtoO
      false,                // kRHCNuisBANFFNormCCNu,
      false,                // kRHCNuisBANFFNormCCNuBar,
      false,                // kRHCNuisBANFFSKFSIAbsorption,
      false,                // kRHCNuisBANFFSKFSIQELow,
      false,                // kRHCNuisBANFFSKFSICXLow,
      false,                // kRHCNuisBANFFSKFSIQEHigh,
      false,                // kRHCNuisBANFFSKFSIInelastic
      false                 // kRHCNuisBANFFEb
#  endif

#endif
    };
        
        
        
    // NEW FLUX FHC+RHC
    static const bool kFHCpRHCParamInBANFFfitIsSuperkFluxParam[kNumOfFHCpRHCParamsInBANFFfit] =
    {
      // NEW XSEC NIWG
#if ( defined(__USE_NIWG2015__)  || defined(__USE_NIWG2017__) )
      // NEW XSEC
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
      true,                 // kNuisBANFF16
      true,                 // kNuisBANFF17
      true,                 // kNuisBANFF18
      true,                 // kNuisBANFF19
      true,                 // kNuisBANFF20
      true,                 // kNuisBANFF21
      true,                 // kNuisBANFF22
      true,                 // kNuisBANFF23
      true,                 // kNuisBANFF24
            
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
      true,                 // kRHCNuisBANFF16
      true,                 // kRHCNuisBANFF17
      true,                 // kRHCNuisBANFF18
      true,                 // kRHCNuisBANFF19
      true,                 // kRHCNuisBANFF20
      true,                 // kRHCNuisBANFF21
      true,                 // kRHCNuisBANFF22
      true,                 // kRHCNuisBANFF23
      true,                 // kRHCNuisBANFF24
            
      //
      false,                // kRHCNuisBANFF2p2h
      false,                // kRHCNuisBANFFCA5
      false,                // kRHCNuisBANFFBgRES
      false,                // kRHCNuisBANFFMAQE
      false,                // kRHCNuisBANFFMARES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      false,                // kRHCNuisBANFFSCCA
      false,                // kRHCNuisBANFFSCCV
#endif
      false,                // kRHCNuisBANFFFermiMom
      false,                // kRHCNuisBANFFShapeCCoth
      false,                // kRHCNuisBANFFNormCCcoh
      false,                // kRHCNuisBANFFNormNCcoh
      false,                // kRHCNuisBANFFNormNCoth
      false,                // kRHCNuisBANFFNormNuEToNuMu
      false,                // kRHCNuisBANFFNormNC1Gamma
      false,                // kRHCNuisBANFFNormNuEBarToNuMuBar
      false,                // kRHCNuisBANFF2p2hBar
#  if defined(__USE_NIWG2017__)      
      false,                // kRHCNuisBANFFBeRPA_A     
      false,                // kRHCNuisBANFFBeRPA_B     
      false,                // kRHCNuisBANFFBeRPA_D     
      false,                // kRHCNuisBANFFBeRPA_E     
      false,                // kRHCNuisBANFFBeRPA_U     
      false,                // kRHCNuisBANFFShape2p2h   
      false,                // kRHCNuisBANFFNorm2p2hCtoO
      false,                // kRHCNuisBANFFNormCCNu,
      false,                // kRHCNuisBANFFNormCCNuBar,
      false,                // kRHCNuisBANFFSKFSIAbsorption,
      false,                // kRHCNuisBANFFSKFSIQELow,
      false,                // kRHCNuisBANFFSKFSICXLow,
      false,                // kRHCNuisBANFFSKFSIQEHigh,
      false,                // kRHCNuisBANFFSKFSIInelastic
      false                 // kRHCNuisBANFFEb

#  endif            
                  
#endif
    };
        
        
    // flux-related BANFF parameters --> flux neutrino PDG code
    static const int kNuPdgCodeOfFluxParamInBANFFfit[kNumOfParamsInBANFFfit] =
    {
      // NEW XSEC NIWG            
#if ( defined(__USE_NIWG2015__)  || defined(__USE_NIWG2017__) )
      // NEW XSEC
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
            
      //
      0,                    // kNuisBANFFNorm2p2h
      0,                    // kNuisBANFFCA5
      0,                    // kNuisBANFFBgRES
      0,                    // kNuisBANFFMAQE
      0,                    // kNuisBANFFMARES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      0,                    // kNuisBANFFSCCA
      0,                    // kNuisBANFFSCCV
#endif
      //
      0,                    // kNuisBANFFFermiMom
      0,                    // kNuisBANFFShapeCCoth
      0,                    // kNuisBANFFNormCCcoh
      0,                    // kNuisBANFFNormNCcoh
      0,                    // kNuisBANFFNormNCoth
      0,                    // kNuisBANFFNormNuEToNuMu
      0,                    // kNuisBANFFNormNC1Gamma
      0,                    // kNuisBANFFNormNuEBarToNuMuBar
      0,                    // kNuisBANFFNorm2p2hBar
#  if defined(__USE_NIWG2017__)      
      0,                    // kNuisBANFFBeRPA_A     
      0,                    // kNuisBANFFBeRPA_B     
      0,                    // kNuisBANFFBeRPA_D     
      0,                    // kNuisBANFFBeRPA_E     
      0,                    // kNuisBANFFBeRPA_U     
      0,                    // kNuisBANFFShape2p2h   
      0,                    // kNuisBANFFNorm2p2hCtoO
      0,                    // kNuisBANFFNormCCNu,
      0,                    // kNuisBANFFNormCCNuBar,
      0,                    // kNuisBANFFSKFSIAbsorption,
      0,                    // kNuisBANFFSKFSIQELow,
      0,                    // kNuisBANFFSKFSICXLow,
      0,                    // kNuisBANFFSKFSIQEHigh,
      0,                    // kNuisBANFFSKFSIInelastic
      0                     // kNuisBANFFEb

#  endif
            
#endif
    };
        
    // NEW FLUX (invert *Nu* w/ *Nu*Bar for RHC mode)
        
    static const int kNuPdgCodeOfRHCFluxParamInBANFFfit[kNumOfRHCParamsInBANFFfit] =
    {            
#if ( defined(__USE_NIWG2015__)  || defined(__USE_NIWG2017__) )
      // NEW XSEC
      kPdgNuMu,             // kRHCNuisBANFF00
      kPdgNuMu,             // kRHCNuisBANFF01
      kPdgNuMu,             // kRHCNuisBANFF02
      kPdgNuMu,             // kRHCNuisBANFF03
      kPdgNuMu,             // kRHCNuisBANFF04
      kPdgNuMuBar,          // kRHCNuisBANFF05
      kPdgNuMuBar,          // kRHCNuisBANFF06
      kPdgNuMuBar,          // kRHCNuisBANFF07
      kPdgNuMuBar,          // kRHCNuisBANFF08
      kPdgNuMuBar,          // kRHCNuisBANFF09
      kPdgNuMuBar,          // kRHCNuisBANFF10
      kPdgNuMuBar,          // kRHCNuisBANFF11
      kPdgNuMuBar,          // kRHCNuisBANFF12
      kPdgNuMuBar,          // kRHCNuisBANFF13
      kPdgNuMuBar,          // kRHCNuisBANFF14
      kPdgNuMuBar,          // kRHCNuisBANFF15
      kPdgNuE,              // kRHCNuisBANFF16
      kPdgNuE,              // kRHCNuisBANFF17
      kPdgNuEBar,           // kRHCNuisBANFF18
      kPdgNuEBar,           // kRHCNuisBANFF19
      kPdgNuEBar,           // kRHCNuisBANFF20
      kPdgNuEBar,           // kRHCNuisBANFF21
      kPdgNuEBar,           // kRHCNuisBANFF22
      kPdgNuEBar,           // kRHCNuisBANFF23
      kPdgNuEBar,           // kRHCNuisBANFF24
      //
      0,                    // kRHCNuisBANFF2p2h
      0,                    // kRHCNuisBANFFCA5
      0,                    // kRHCNuisBANFFBgRES
      0,                    // kRHCNuisBANFFMAQE
      0,                    // kRHCNuisBANFFMARES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      0,                    // kRHCNuisBANFFSCCA
      0,                    // kRHCNuisBANFFSCCV
#endif
      0,                    // kRHCNuisBANFFFermiMom
      0,                    // kRHCNuisBANFFShapeCCoth
      0,                    // kRHCNuisBANFFNormCCcoh
      0,                    // kRHCNuisBANFFNormNCcoh
      0,                    // kRHCNuisBANFFNormNCoth
      0,                    // kRHCNuisBANFFNormNuEToNuMu
      0,                    // kRHCNuisBANFFNormNC1Gamma
      0,                    // kRHCNuisBANFFNormNuEBarToNuMuBar
      0,                    // kRHCNuisBANFF2p2hBar
#  if defined(__USE_NIWG2017__)      
      0,                    // kRHCNuisBANFFBeRPA_A     
      0,                    // kRHCNuisBANFFBeRPA_B     
      0,                    // kRHCNuisBANFFBeRPA_D     
      0,                    // kRHCNuisBANFFBeRPA_E     
      0,                    // kRHCNuisBANFFBeRPA_U     
      0,                    // kRHCNuisBANFFShape2p2h   
      0,                    // kRHCNuisBANFFNorm2p2hCtoO
      0,                    // kRHCNuisBANFFNormCCNu,
      0,                    // kRHCNuisBANFFNormCCNuBar,
      0,                    // kRHCNuisBANFFSKFSIAbsorption,
      0,                    // kRHCNuisBANFFSKFSIQELow,
      0,                    // kRHCNuisBANFFSKFSICXLow,
      0,                    // kRHCNuisBANFFSKFSIQEHigh,
      0,                    // kRHCNuisBANFFSKFSIInelastic
      0                     // kRHCNuisBANFFEb

#  endif
      
#endif
    };
    ///
        
        
    // NEW FLUX FHC+RHC (invert *Nu* w/ *Nu*Bar for RHC mode)
        
    static const int kNuPdgCodeOfFHCpRHCFluxParamInBANFFfit[kNumOfFHCpRHCParamsInBANFFfit] =
    {
      // NEW XSEC NIWG            
#if ( defined(__USE_NIWG2015__)  || defined(__USE_NIWG2017__) )
      // NEW XSEC
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
      //
      kPdgNuMu,             // kRHCNuisBANFF00
      kPdgNuMu,             // kRHCNuisBANFF01
      kPdgNuMu,             // kRHCNuisBANFF02
      kPdgNuMu,             // kRHCNuisBANFF03
      kPdgNuMu,             // kRHCNuisBANFF04
      kPdgNuMuBar,          // kRHCNuisBANFF05
      kPdgNuMuBar,          // kRHCNuisBANFF06
      kPdgNuMuBar,          // kRHCNuisBANFF07
      kPdgNuMuBar,          // kRHCNuisBANFF08
      kPdgNuMuBar,          // kRHCNuisBANFF09
      kPdgNuMuBar,          // kRHCNuisBANFF10
      kPdgNuMuBar,          // kRHCNuisBANFF11
      kPdgNuMuBar,          // kRHCNuisBANFF12
      kPdgNuMuBar,          // kRHCNuisBANFF13
      kPdgNuMuBar,          // kRHCNuisBANFF14
      kPdgNuMuBar,          // kRHCNuisBANFF15
            
      kPdgNuE,              // kRHCNuisBANFF16
      kPdgNuE,              // kRHCNuisBANFF17
      kPdgNuEBar,           // kRHCNuisBANFF18
      kPdgNuEBar,           // kRHCNuisBANFF19
      kPdgNuEBar,           // kRHCNuisBANFF20
      kPdgNuEBar,           // kRHCNuisBANFF21
      kPdgNuEBar,           // kRHCNuisBANFF22
      kPdgNuEBar,           // kRHCNuisBANFF23
      kPdgNuEBar,           // kRHCNuisBANFF24
            
      //
      0,                    // kRHCNuisBANFF2p2h
      0,                    // kRHCNuisBANFFCA5
      0,                    // kRHCNuisBANFFBgRES
      0,                    // kRHCNuisBANFFMAQE
      0,                    // kRHCNuisBANFFMARES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      0,                    // kRHCNuisBANFFSCCA
      0,                    // kRHCNuisBANFFSCCV
#endif
      0,                    // kRHCNuisBANFFFermiMom
      0,                    // kRHCNuisBANFFShapeCCoth
      0,                    // kRHCNuisBANFFNormCCcoh
      0,                    // kRHCNuisBANFFNormNCcoh
      0,                    // kRHCNuisBANFFNormNCoth
      0,                    // kRHCNuisBANFFNormNuEToNuMu
      0,                    // kRHCNuisBANFFNormNC1Gamma
      0,                    // kNuisBANFFNormNuEBarToNuMuBar
      0,                    // kRHCNuisBANFF2p2hBar
#  if defined(__USE_NIWG2017__)      
      0,                    // kRHCNuisBANFFBeRPA_A     
      0,                    // kRHCNuisBANFFBeRPA_B     
      0,                    // kRHCNuisBANFFBeRPA_D     
      0,                    // kRHCNuisBANFFBeRPA_E     
      0,                    // kRHCNuisBANFFBeRPA_U     
      0,                    // kRHCNuisBANFFShape2p2h   
      0,                    // kRHCNuisBANFFNorm2p2hCtoO
      0,                    // kRHCNuisBANFFNormCCNu,
      0,                    // kRHCNuisBANFFNormCCNuBar,
      0,                    // kRHCNuisBANFFSKFSIAbsorption,
      0,                    // kRHCNuisBANFFSKFSIQELow,
      0,                    // kRHCNuisBANFFSKFSICXLow,
      0,                    // kRHCNuisBANFFSKFSIQEHigh,
      0,                    // kRHCNuisBANFFSKFSIInelastic
      0                     // kRHCNuisBANFFEb

#  endif            
#endif
    };
    ///
        
        
        
    // Etrue energy range for flux-related BANFF parameters
    static const double kEtrueRangeOfFluxParamInBANFFfit[kNumOfParamsInBANFFfit][2] =
    {
      // NEW XSEC NIWG            
#if ( defined(__USE_NIWG2015__)  || defined(__USE_NIWG2017__) )
      // NEW XSEC
      { 0.0,  0.4 },        // kNuisBANFF00 //numu
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
            
      { 0.0,  0.7 },        // kNuisBANFF11 //numub
      { 0.7,  1.0 },        // kNuisBANFF12
      { 1.0,  1.5 },        // kNuisBANFF13
      { 1.5,  2.5 },        // kNuisBANFF14
      { 2.5, 30.0 },        // kNuisBANFF15
            
      { 0.0,  0.5 },        // kNuisBANFF16 //nue
      { 0.5,  0.7 },        // kNuisBANFF17
      { 0.7,  0.8 },        // kNuisBANFF18
      { 0.8,  1.5 },        // kNuisBANFF19
      { 1.5,  2.5 },        // kNuisBANFF20
      { 2.5,  4.0 },        // kNuisBANFF21
      { 4.0, 30.0 },        // kNuisBANFF22
            
      { 0.0,  2.5 },        // kNuisBANFF23 //nueb
      { 2.5, 30.0 },        // kNuisBANFF24
            
      //
      {  -1,   -1 },        // kNuisBANFFNorm2p2h not a flux param
      {  -1,   -1 },        // kNuisBANFFCA5 not a flux param
      {  -1,   -1 },        // kNuisBANFFBgRES not a flux param
      {  -1,   -1 },        // kNuisBANFFMAQE, not a flux param
      {  -1,   -1 },        // kNuisBANFFMARES, not a flux param
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      {  -1,   -1 },        // kNuisBANFFSCCA, not a flux param
      {  -1,   -1 },        // kNuisBANFFSCCV, not a flux param
#endif
      {  -1,   -1 },        // kNuisBANFFFermiMom, not a flux param
      {  -1,   -1 },        // kNuisBANFFShapeCCoth, not a flux param
      {  -1,   -1 },        // kNuisBANFFNormCCcoh, not a flux param
      {  -1,   -1 },        // kNuisBANFFNormNCcoh, not a flux param
      {  -1,   -1 },        // kNuisBANFFNormNCoth, not a flux param
      {  -1,   -1 },        // kNuisBANFFNormNuEToNuMu, not a flux param
      {  -1,   -1 },        // kNuisBANFFNormNC1Gamma, not a flux param
      {  -1,   -1 },        // kNuisBANFFNormNuEBarToNuMuBar, not a flux param
      {  -1,   -1 },        // kNuisBANFFNorm2p2hBar not a flux param
#  if defined(__USE_NIWG2017__)      
      {  -1,   -1 },        // kNuisBANFFBeRPA_A     
      {  -1,   -1 },        // kNuisBANFFBeRPA_B     
      {  -1,   -1 },        // kNuisBANFFBeRPA_D     
      {  -1,   -1 },        // kNuisBANFFBeRPA_E     
      {  -1,   -1 },        // kNuisBANFFBeRPA_U     
      {  -1,   -1 },        // kNuisBANFFShape2p2h   
      {  -1,   -1 },        // kNuisBANFFNorm2p2hCtoO
      {  -1,   -1 },        // kNuisBANFFNormCCNu,
      {  -1,   -1 },        // kNuisBANFFNormCCNuBar,
      {  -1,   -1 },        // kNuisBANFFSKFSIAbsorption,
      {  -1,   -1 },        // kNuisBANFFSKFSIQELow,
      {  -1,   -1 },        // kNuisBANFFSKFSICXLow,
      {  -1,   -1 },        // kNuisBANFFSKFSIQEHigh,
      {  -1,   -1 },        // kNuisBANFFSKFSIInelastic
      {  -1,   -1 }         // kNuisBANFFEb

#  endif            
#endif
    };
        
        
    // NEW FLUX
        
    static const double kEtrueRangeOfRHCFluxParamInBANFFfit[kNumOfRHCParamsInBANFFfit][2] =
    {
      // NEW XSEC NIWG            
#if ( defined(__USE_NIWG2015__)  || defined(__USE_NIWG2017__) )
      // NEW XSEC
      // numu
      { 0.0,  0.7 },        // kRHCNuisBANFF00
      { 0.7,  1.0 },        // kRHCNuisBANFF01
      { 1.0,  1.5 },        // kRHCNuisBANFF02
      { 1.5,  2.5 },        // kRHCNuisBANFF03
      { 2.5, 30.0 },        // kRHCNuisBANFF04
      // anti-numu
      { 0.0,  0.4 },        // kRHCNuisBANFF05
      { 0.4,  0.5 },        // kRHCNuisBANFF06
      { 0.5,  0.6 },        // kRHCNuisBANFF07
      { 0.6,  0.7 },        // kRHCNuisBANFF08
      { 0.7,  1.0 },        // kRHCNuisBANFF09
      { 1.0,  1.5 },        // kRHCNuisBANFF10
      { 1.5,  2.5 },        // kRHCNuisBANFF11
      { 2.5,  3.5 },        // kRHCNuisBANFF12
      { 3.5,  5.0 },        // kRHCNuisBANFF13
      { 5.0,  7.0 },        // kRHCNuisBANFF14
      { 7.0, 30.0 },        // kRHCNuisBANFF15
      // nue
      { 0.0,  2.5 },        // kRHCNuisBANFF16
      { 2.5, 30.0 },        // kRHCNuisBANFF17
      // anti-nue 
      { 0.0,  0.5 },        // kRHCNuisBANFF18
      { 0.5,  0.7 },        // kRHCNuisBANFF19
      { 0.7,  0.8 },        // kRHCNuisBANFF20
      { 0.8,  1.5 },        // kRHCNuisBANFF21
      { 1.5,  2.5 },        // kRHCNuisBANFF22
      { 2.5,  4.0 },        // kRHCNuisBANFF23
      { 4.0, 30.0 },        // kRHCNuisBANFF24
            
      //
      {  -1,   -1 },        // kRHCNuisBANFF2p2h not a flux param
      {  -1,   -1 },        // kRHCNuisBANFFCA5 not a flux param  
      {  -1,   -1 },        // kRHCNuisBANFFBgRES not a flux param
      {  -1,   -1 },        // kRHCNuisBANFFMAQE, not a flux param
      {  -1,   -1 },        // kRHCNuisBANFFMARES, not a flux param
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      {  -1,   -1 },        // kRHCNuisBANFFSCCA, not a flux param
      {  -1,   -1 },        // kRHCNuisBANFFSCCV, not a flux param
#endif
      {  -1,   -1 },        // kRHCNuisBANFFFermiMom, not a flux param
      {  -1,   -1 },        // kRHCNuisBANFFShapeCCoth, not a flux param
      {  -1,   -1 },        // kRHCNuisBANFFNormCCcoh, not a flux param
      {  -1,   -1 },        // kRHCNuisBANFFNormNCcoh, not a flux param
      {  -1,   -1 },        // kRHCNuisBANFFNormNCoth, not a flux param
      {  -1,   -1 },        // kRHCNuisBANFFNormNuEToNuMu, not a flux param
      {  -1,   -1 },        // kRHCNuisBANFFNormNC1Gamma, not a flux param
      {  -1,   -1 },        // kRHCNuisBANFFNormNuEBarToNuMuBar, not a flux param
      {  -1,   -1 },        // kRHCNuisBANFF2p2hBar not a flux param
#  if defined(__USE_NIWG2017__)      
      {  -1,   -1 },        // kRHCNuisBANFFBeRPA_A     
      {  -1,   -1 },        // kRHCNuisBANFFBeRPA_B     
      {  -1,   -1 },        // kRHCNuisBANFFBeRPA_D     
      {  -1,   -1 },        // kRHCNuisBANFFBeRPA_E     
      {  -1,   -1 },        // kRHCNuisBANFFBeRPA_U     
      {  -1,   -1 },        // kRHCNuisBANFFShape2p2h   
      {  -1,   -1 },        // kRHCNuisBANFFNorm2p2hCtoO
      {  -1,   -1 },        // kRHCNuisBANFFNormCCNu,
      {  -1,   -1 },        // kRHCNuisBANFFNormCCNuBar,
      {  -1,   -1 },        // kRHCNuisBANFFSKFSIAbsorption,
      {  -1,   -1 },        // kRHCNuisBANFFSKFSIQELow,
      {  -1,   -1 },        // kRHCNuisBANFFSKFSICXLow,
      {  -1,   -1 },        // kRHCNuisBANFFSKFSIQEHigh,
      {  -1,   -1 },        // kRHCNuisBANFFSKFSIInelastic
      {  -1,   -1 }         // kRHCNuisBANFFEb

#  endif                  
#endif
    };
        
        
    // NEW FLUX FHC+RHC
    // NEW FLUX
        
    static const double kEtrueRangeOfFHCpRHCFluxParamInBANFFfit[kNumOfFHCpRHCParamsInBANFFfit][2] =
    {
        // NEW XSEC NIWG            
#if ( defined(__USE_NIWG2015__)  || defined(__USE_NIWG2017__) )
      // NEW XSEC
      // FHC
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
            
      // RHC
      // numu
      { 0.0,  0.7 },        // kRHCNuisBANFF00
      { 0.7,  1.0 },        // kRHCNuisBANFF01
      { 1.0,  1.5 },        // kRHCNuisBANFF02
      { 1.5,  2.5 },        // kRHCNuisBANFF03
      { 2.5, 30.0 },        // kRHCNuisBANFF04
      // anti-numu
      { 0.0,  0.4 },        // kRHCNuisBANFF05
      { 0.4,  0.5 },        // kRHCNuisBANFF06
      { 0.5,  0.6 },        // kRHCNuisBANFF07
      { 0.6,  0.7 },        // kRHCNuisBANFF08
      { 0.7,  1.0 },        // kRHCNuisBANFF09
      { 1.0,  1.5 },        // kRHCNuisBANFF10
      { 1.5,  2.5 },        // kRHCNuisBANFF11
      { 2.5,  3.5 },        // kRHCNuisBANFF12
      { 3.5,  5.0 },        // kRHCNuisBANFF13
      { 5.0,  7.0 },        // kRHCNuisBANFF14
      { 7.0, 30.0 },        // kRHCNuisBANFF15
            
      // nue
      { 0.0,  2.5 },        // kRHCNuisBANFF16
      { 2.5, 30.0 },        // kRHCNuisBANFF17
      // anti-nue 
      { 0.0,  0.5 },        // kRHCNuisBANFF18
      { 0.5,  0.7 },        // kRHCNuisBANFF19
      { 0.7,  0.8 },        // kRHCNuisBANFF20
      { 0.8,  1.5 },        // kRHCNuisBANFF21
      { 1.5,  2.5 },        // kRHCNuisBANFF22
      { 2.5,  4.0 },        // kRHCNuisBANFF23
      { 4.0, 30.0 },        // kRHCNuisBANFF24
            
      //
      {  -1,   -1 },        // kNuisBANFFNorm2p2h not a flux param
      {  -1,   -1 },        // kNuisBANFFCA5 not a flux param  
      {  -1,   -1 },        // kNuisBANFFBgRES not a flux param
      {  -1,   -1 },        // kNuisBANFFMAQE, not a flux param
      {  -1,   -1 },        // kNuisBANFFMARES, not a flux param
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      {  -1,   -1 },        // kNuisBANFFSCCA, not a flux param
      {  -1,   -1 },        // kNuisBANFFSCCV, not a flux param
#endif
      {  -1,   -1 },        // kNuisBANFFFermiMom, not a flux param
      {  -1,   -1 },        // kNuisBANFFShapeCCoth, not a flux param
      {  -1,   -1 },        // kNuisBANFFNormCCcoh, not a flux param
      {  -1,   -1 },        // kNuisBANFFNormNCcoh, not a flux param
      {  -1,   -1 },        // kNuisBANFFNormNCoth, not a flux param
      {  -1,   -1 },        // kNuisBANFFNormNuEToNuMu, not a flux param
      {  -1,   -1 },        // kNuisBANFFNormNC1Gamma, not a flux param
      {  -1,   -1 },        // kNuisBANFFNormNuEBarToNuMuBar, not a flux param
      {  -1,   -1 },        // kNuisBANFFNorm2p2hBar not a flux param
#  if defined(__USE_NIWG2017__)      
      {  -1,   -1 },        // kNuisBANFFBeRPA_A     
      {  -1,   -1 },        // kNuisBANFFBeRPA_B     
      {  -1,   -1 },        // kNuisBANFFBeRPA_D     
      {  -1,   -1 },        // kNuisBANFFBeRPA_E     
      {  -1,   -1 },        // kNuisBANFFBeRPA_U     
      {  -1,   -1 },        // kNuisBANFFShape2p2h   
      {  -1,   -1 },        // kNuisBANFFNorm2p2hCtoO
      {  -1,   -1 },        // kNuisBANFFNormCCNu,
      {  -1,   -1 },        // kNuisBANFFNormCCNuBar,
      {  -1,   -1 },        // kNuisBANFFSKFSIAbsorption,
      {  -1,   -1 },        // kNuisBANFFSKFSIQELow,
      {  -1,   -1 },        // kNuisBANFFSKFSICXLow,
      {  -1,   -1 },        // kNuisBANFFSKFSIQEHigh,
      {  -1,   -1 },        // kNuisBANFFSKFSIInelastic
      {  -1,   -1 }         // kNuisBANFFEb

#  endif                        
#endif
    };
    //
        
  } // conv
} // valor

#endif 

