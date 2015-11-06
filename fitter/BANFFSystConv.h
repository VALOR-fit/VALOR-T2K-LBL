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
        
        //
        // Number of BANFF fit parameters
        //
        
        // NEW XSEC NIWG
#if defined(__USE_NIWG2015__)
#if BANFF_BINS == 1 //joint
#ifdef MECBAR_ERROR
        static const int kNumOfParamsInBANFFfit = 40;
        static const int kNumOfFluxParamsInBANFFfit = kNumOfParamsInBANFFfit - 15;
        static const int kNumOfRHCParamsInBANFFfit = 40;
        static const int kNumOfRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - 15;
        // NEW FLUX FHC+RHC
        static const int kNumOfFHCpRHCParamsInBANFFfit = 65;
        static const int kNumOfFHCpRHCFluxParamsInBANFFfit = kNumOfFHCpRHCParamsInBANFFfit - 15;
        
#else
        static const int kNumOfParamsInBANFFfit = 39;
        static const int kNumOfFluxParamsInBANFFfit = kNumOfParamsInBANFFfit - 14;
        static const int kNumOfRHCParamsInBANFFfit = 39;
        static const int kNumOfRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - 14;
        // NEW FLUX FHC+RHC
        static const int kNumOfFHCpRHCParamsInBANFFfit = 64;
        static const int kNumOfFHCpRHCFluxParamsInBANFFfit = kNumOfFHCpRHCParamsInBANFFfit - 14;
        //
#endif
#elif BANFF_BINS == 2 //disap (no scca and sccv)
        static const int kNumOfParamsInBANFFfit = 37;
        static const int kNumOfFluxParamsInBANFFfit = kNumOfParamsInBANFFfit - 12; //25
        static const int kNumOfRHCParamsInBANFFfit = 37;
        static const int kNumOfRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - 12; //25
        // NEW FLUX FHC+RHC
        static const int kNumOfFHCpRHCParamsInBANFFfit = 62;
        static const int kNumOfFHCpRHCFluxParamsInBANFFfit = kNumOfFHCpRHCParamsInBANFFfit - 12; //50;
        //
#elif BANFF_BINS == 3 //app
        static const int kNumOfParamsInBANFFfit = 39;
        static const int kNumOfFluxParamsInBANFFfit = kNumOfParamsInBANFFfit - 14;
        static const int kNumOfRHCParamsInBANFFfit = 39;
        static const int kNumOfRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - 14;
        // NEW FLUX FHC+RHC
        static const int kNumOfFHCpRHCParamsInBANFFfit = 64;
        static const int kNumOfFHCpRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - 14;
        //
#endif
        
#elif defined(__USE_NIWG2013__)
#if BANFF_BINS == 1 //joint
        static const int kNumOfParamsInBANFFfit = 33;
        static const int kNumOfFluxParamsInBANFFfit = kNumOfParamsInBANFFfit - 8;
        static const int kNumOfRHCParamsInBANFFfit = 33;
        static const int kNumOfRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - 8;
        // NEW FLUX FHC+RHC
        static const int kNumOfFHCpRHCParamsInBANFFfit = 58;
        static const int kNumOfFHCpRHCFluxParamsInBANFFfit = kNumOfFHCpRHCParamsInBANFFfit - 8;
        //
#elif BANFF_BINS == 2 //disap
        static const int kNumOfParamsInBANFFfit = 23;
        static const int kNumOfFluxParamsInBANFFfit = kNumOfParamsInBANFFfit - 7;
        static const int kNumOfRHCParamsInBANFFfit = 23;
        static const int kNumOfRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - 7;
        // NEW FLUX FHC+RHC
        static const int kNumOfFHCpRHCParamsInBANFFfit = 39;
        static const int kNumOfFHCpRHCFluxParamsInBANFFfit = kNumOfFHCpRHCParamsInBANFFfit - 7;
        //
#elif BANFF_BINS == 3 //app
        static const int kNumOfParamsInBANFFfit = 30; //27;
        static const int kNumOfFluxParamsInBANFFfit = kNumOfParamsInBANFFfit - 5;
        static const int kNumOfRHCParamsInBANFFfit = 30; //27;
        static const int kNumOfRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - 5;
        // NEW FLUX FHC+RHC
        static const int kNumOfFHCpRHCParamsInBANFFfit = 55; //49;
        static const int kNumOfFHCpRHCFluxParamsInBANFFfit = kNumOfRHCParamsInBANFFfit - 5;
        //
#endif
#endif
        
        
        static const PhysParamId_t kListOfBANFFParams[kNumOfParamsInBANFFfit] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
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
#if (BANFF_BINS != 2)
            kNuisBANFF16 ,
            kNuisBANFF17 ,
            kNuisBANFF18 ,
            kNuisBANFF19 ,
            kNuisBANFF20 ,
            kNuisBANFF21 ,
            // NEW SPLIT ANA
            //#endif
            //#if BANFF_BINS == 1
            kNuisBANFF22 ,
            kNuisBANFF23 ,
            kNuisBANFF24 ,
#endif
            
            kNuisBANFFMAQE   ,
            kNuisBANFFMARES  ,
            kNuisBANFFCCQEE0 ,
#if BANFF_BINS != 3
            kNuisBANFFCCQEE1 ,
            kNuisBANFFCCQEE2 ,
#endif
            kNuisBANFFCC1piE0 ,
#if BANFF_BINS != 3
            kNuisBANFFCC1piE1,
#endif
#if (BANFF_BINS != 2)
            kNuisBANFFNC1pi0 ,
#endif
            
#elif defined(__USE_NIWG2015__)
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
            
            kNuisBANFFMEC          ,
            kNuisBANFFCA5          ,
            kNuisBANFFBgRES        ,
            kNuisBANFFMAQE         ,
            kNuisBANFFMARES        ,
            // NEW SPLIT ANA
#if (BANFF_BINS != 2)
            kNuisBANFFSCCA         ,
            kNuisBANFFSCCV         ,
#endif
            //
            kNuisBANFFFermiMom      ,
            kNuisBANFFShapeCCoth    ,
            kNuisBANFFBindE         ,
            kNuisBANFFNormCCcoh     ,
            kNuisBANFFNormNCcoh     ,
            kNuisBANFFNormNCoth     ,
            kNuisBANFFNormNuEToNuMu ,
#ifdef MECBAR_ERROR
            kNuisBANFFMECBar,
#endif
#endif
        };
        
        static const PhysParamId_t kListOfRHCBANFFParams[kNumOfRHCParamsInBANFFfit] =
        {
#if defined(__USE_NIWG2013__)
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
#if (BANFF_BINS != 2)
            kRHCNuisBANFF16 ,
            kRHCNuisBANFF17 ,
            kRHCNuisBANFF18 ,
            kRHCNuisBANFF19 ,
            kRHCNuisBANFF20 ,
            kRHCNuisBANFF21 ,
            // NEW SPLIT ANA
            //#endif
            //#if BANFF_BINS == 1
            kRHCNuisBANFF22 ,
            kRHCNuisBANFF23 ,
            kRHCNuisBANFF24 ,
#endif
            
            kRHCNuisBANFFMAQE   ,
            kRHCNuisBANFFMARES  ,
            kRHCNuisBANFFCCQEE0 ,
#if BANFF_BINS != 3
            kRHCNuisBANFFCCQEE1 ,
            kRHCNuisBANFFCCQEE2 ,
#endif
            kRHCNuisBANFFCC1piE0 ,
#if BANFF_BINS != 3
            kRHCNuisBANFFCC1piE1,
#endif
#if (BANFF_BINS != 2)
            kRHCNuisBANFFNC1pi0 ,
#endif
            
#elif defined(__USE_NIWG2015__)
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
            
            kRHCNuisBANFFMEC   ,
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
            kRHCNuisBANFFBindE         ,
            kRHCNuisBANFFNormCCcoh     ,
            kRHCNuisBANFFNormNCcoh     ,
            kRHCNuisBANFFNormNCoth     ,
            kRHCNuisBANFFNormNuEToNuMu ,
#ifdef MECBAR_ERROR
            kRHCNuisBANFFMECBar,
#endif
#endif
        };
        
        static const PhysParamId_t kListOfFHCpRHCBANFFParams[kNumOfFHCpRHCParamsInBANFFfit] =
        {
#if defined(__USE_NIWG2013__)
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
#if (BANFF_BINS != 2)
            kNuisBANFF16 ,
            kNuisBANFF17 ,
            kNuisBANFF18 ,
            kNuisBANFF19 ,
            kNuisBANFF20 ,
            kNuisBANFF21 ,
            
            // NEW SPLIT ANA
            //#endif
            //#if BANFF_BINS == 1
            kNuisBANFF22 ,
            kNuisBANFF23 ,
            kNuisBANFF24 ,
#endif
            
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
#if (BANFF_BINS != 2)
            kRHCNuisBANFF16 ,
            kRHCNuisBANFF17 ,
            kRHCNuisBANFF18 ,
            kRHCNuisBANFF19 ,
            kRHCNuisBANFF20 ,
            kRHCNuisBANFF21 ,
            
            // NEW SPLIT ANA
            //#endif
            //#if BANFF_BINS == 1
            kRHCNuisBANFF22 ,
            kRHCNuisBANFF23 ,
            kRHCNuisBANFF24 ,
#endif
            
            //
            kNuisBANFFMAQE   ,
            kNuisBANFFMARES  ,
            kNuisBANFFCCQEE0 ,
#if BANFF_BINS != 3
            kNuisBANFFCCQEE1 ,
            kNuisBANFFCCQEE2 ,
#endif
            kNuisBANFFCC1piE0 ,
#if BANFF_BINS != 3
            kNuisBANFFCC1piE1,
#endif
#if (BANFF_BINS != 2)
            kNuisBANFFNC1pi0 ,
#endif
            
#elif defined(__USE_NIWG2015__)
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
            kNuisBANFFMEC,
            kNuisBANFFCA5,
            kNuisBANFFBgRES,
            kNuisBANFFMAQE   ,
            kNuisBANFFMARES  ,
            // NEW SPLIT ANA
#if (BANFF_BINS != 2)
            kNuisBANFFSCCA   ,
            kNuisBANFFSCCV   ,
#endif
            kNuisBANFFFermiMom      ,
            kNuisBANFFShapeCCoth    ,
            kNuisBANFFBindE         ,
            kNuisBANFFNormCCcoh     ,
            kNuisBANFFNormNCcoh     ,
            kNuisBANFFNormNCoth     ,
            kNuisBANFFNormNuEToNuMu ,
#ifdef MECBAR_ERROR
            kNuisBANFFMECBar,
#endif
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
#if defined(__USE_NIWG2013__)
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
#if (BANFF_BINS != 2)
            true,                 // kNuisBANFF16
            true,                 // kNuisBANFF17
            true,                 // kNuisBANFF18
            true,                 // kNuisBANFF19
            true,                 // kNuisBANFF20
            true,                 // kNuisBANFF21
            // NEW SPLIT ANA
            //#endif
            //#if BANFF_BINS == 1
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
#if (BANFF_BINS != 2)
            false,                 // kNuisBANFFNC1pi0
#endif
            
#elif defined(__USE_NIWG2015__)
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
            false,                // kNuisBANFFMEC
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
            false,                // kNuisBANFFBindE
            false,                // kNuisBANFFNormCCcoh
            false,                // kNuisBANFFNormNCcoh
            false,                // kNuisBANFFNormNCoth
            false,                // kNuisBANFFNormNuEToNuMu
#ifdef MECBAR_ERROR
            false,                  //    kNuisBANFFMECBar
#endif
            //
#endif
        };
        
        static const bool kRHCParamInBANFFfitIsSuperkFluxParam[kNumOfRHCParamsInBANFFfit] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
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
#if (BANFF_BINS != 2)
            true,                 // kRHCNuisBANFF16
            true,                 // kRHCNuisBANFF17
            true,                 // kRHCNuisBANFF18
            true,                 // kRHCNuisBANFF19
            true,                 // kRHCNuisBANFF20
            true,                 // kRHCNuisBANFF21
            // NEW SPLIT ANA
            //#endif
            //#if BANFF_BINS == 1
            true,                 // kRHCNuisBANFF22
            true,                 // kRHCNuisBANFF23
            true,                 // kRHCNuisBANFF24
#endif
            //
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
#if (BANFF_BINS != 2)
            false,                 // kRHCNuisBANFFNC1pi0
#endif
            
#elif defined(__USE_NIWG2015__)
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
            false,                // kRHCNuisBANFFMEC
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
            false,                // kRHCNuisBANFFBindE
            false,                // kRHCNuisBANFFNormCCcoh
            false,                // kRHCNuisBANFFNormNCcoh
            false,                // kRHCNuisBANFFNormNCoth
            false,                // kRHCNuisBANFFNormNuEToNuMu
#ifdef MECBAR_ERROR
            false,                  //     kRHCNuisBANFFMECBar,
            
#endif
            //
#endif
        };
        
        
        
        // NEW FLUX FHC+RHC
        static const bool kFHCpRHCParamInBANFFfitIsSuperkFluxParam[kNumOfFHCpRHCParamsInBANFFfit] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
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
#if (BANFF_BINS != 2)
            true,                 // kNuisBANFF16
            true,                 // kNuisBANFF17
            true,                 // kNuisBANFF18
            true,                 // kNuisBANFF19
            true,                 // kNuisBANFF20
            true,                 // kNuisBANFF21
            // NEW SPLIT ANA
            //#endif
            //#if BANFF_BINS == 1
            true,                 // kNuisBANFF22
            true,                 // kNuisBANFF23
            true,                 // kNuisBANFF24
#endif
            //
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
#if (BANFF_BINS != 2)
            true,                 // kRHCNuisBANFF16
            true,                 // kRHCNuisBANFF17
            true,                 // kRHCNuisBANFF18
            true,                 // kRHCNuisBANFF19
            true,                 // kRHCNuisBANFF20
            true,                 // kRHCNuisBANFF21
            // NEW SPLIT ANA
            //#endif
            //#if BANFF_BINS == 1
            true,                 // kRHCNuisBANFF22
            true,                 // kRHCNuisBANFF23
            true,                 // kRHCNuisBANFF24
#endif
            //
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
#if (BANFF_BINS != 2)
            false,                 // CNuisBANFFNC1pi0
#endif
            
#elif defined(__USE_NIWG2015__)
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
            false,                // kRHCNuisBANFFMEC
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
            false,                // kRHCNuisBANFFBindE
            false,                // kRHCNuisBANFFNormCCcoh
            false,                // kRHCNuisBANFFNormNCcoh
            false,                // kRHCNuisBANFFNormNCoth
            false,                // kRHCNuisBANFFNormNuEToNuMu
#ifdef MECBAR_ERROR
            false,                  //    kRHCNuisBANFFMECBar,
#endif
#endif
        };
        
        
        // flux-related BANFF parameters --> flux neutrino PDG code
        static const int kNuPdgCodeOfFluxParamInBANFFfit[kNumOfParamsInBANFFfit] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
#if (BANFF_BINS == 2)
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
            
            // NEW SPLIT ANA
            //#endif
            /* #if BANFF_BINS == 3 */
            /* kPdgNuMu,             // kNuisBANFF00 */
            /* kPdgNuMu,             // kNuisBANFF01 */
            /* kPdgNuMu,             // kNuisBANFF02 */
            /* kPdgNuMu,             // kNuisBANFF03 */
            /* kPdgNuMu,             // kNuisBANFF04 */
            /* kPdgNuMu,             // kNuisBANFF05 */
            /* kPdgNuMu,             // kNuisBANFF06 */
            /* kPdgNuMu,             // kNuisBANFF07 */
            /* kPdgNuMu,             // kNuisBANFF08 */
            /* kPdgNuMu,             // kNuisBANFF09 */
            /* kPdgNuMu,             // kNuisBANFF10 */
            /* kPdgNuMuBar,          // kNuisBANFF11 */
            /* kPdgNuMuBar,          // kNuisBANFF12 */
            /* kPdgNuE,              // kNuisBANFF13 */
            /* kPdgNuE,              // kNuisBANFF14 */
            /* kPdgNuE,              // kNuisBANFF15 */
            /* kPdgNuE,              // kNuisBANFF16 */
            /* kPdgNuE,              // kNuisBANFF17 */
            /* kPdgNuE,              // kNuisBANFF18 */
            /* kPdgNuE,              // kNuisBANFF19 */
            /* kPdgNuEBar,           // kNuisBANFF20 */
            /* kPdgNuEBar,           // kNuisBANFF21 */
            /* #endif */
            //#if BANFF_BINS == 1
#else
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
#if (BANFF_BINS != 2)
            0,                     // kNuisBANFFNC1pi0, not a flux param
#endif
            
#elif defined(__USE_NIWG2015__)
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
            0,                    // kNuisBANFFMEC
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
            0,                    // kNuisBANFFBindE
            0,                    // kNuisBANFFNormCCcoh
            0,                    // kNuisBANFFNormNCcoh
            0,                    // kNuisBANFFNormNCoth
            0,                    // kNuisBANFFNormNuEToNuMu
#ifdef MECBAR_ERROR
            0,                      //    kNuisBANFFMECBar
#endif
#endif
        };
        
        // NEW FLUX (invert *Nu* w/ *Nu*Bar for RHC mode)
        
        static const int kNuPdgCodeOfRHCFluxParamInBANFFfit[kNumOfRHCParamsInBANFFfit] =
        {
#if defined(__USE_NIWG2013__)
#if (BANFF_BINS == 2)
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
            
            // NEW SPLIT ANA
            //#endif
            /* #if BANFF_BINS == 3 */
            /* kPdgNuMu,             // kRHCNuisBANFF00 */
            /* kPdgNuMu,             // kRHCNuisBANFF01 */
            /* kPdgNuMu,             // kRHCNuisBANFF02 */
            /* kPdgNuMu,             // kRHCNuisBANFF03 */
            /* kPdgNuMu,             // kRHCNuisBANFF04 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF05 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF06 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF07 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF08 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF09 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF10 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF11 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF12 */
            /* kPdgNuE,              // kRHCNuisBANFF13 */
            /* kPdgNuE,              // kRHCNuisBANFF14 */
            /* kPdgNuEBar,           // kRHCNuisBANFF15 */
            /* kPdgNuEBar,           // kRHCNuisBANFF16 */
            /* kPdgNuEBar,           // kRHCNuisBANFF17 */
            /* kPdgNuEBar,           // kRHCNuisBANFF18 */
            /* kPdgNuEBar,           // kRHCNuisBANFF19 */
            /* kPdgNuEBar,           // kRHCNuisBANFF20 */
            /* kPdgNuEBar,           // kRHCNuisBANFF21 */
            /* #endif */
            //#if BANFF_BINS == 1
#else
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
#if (BANFF_BINS != 2)
            0,                     // kRHCNuisBANFFNC1pi0, not a flux param
#ifdef MECBAR_ERROR
            0,                      //     kRHCNuisBANFFMECBar, not a flux param
            
#endif
#endif
            
#elif defined(__USE_NIWG2015__)
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
            0,                    // kRHCNuisBANFFMEC
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
            0,                    // kRHCNuisBANFFBindE
            0,                    // kRHCNuisBANFFNormCCcoh
            0,                    // kRHCNuisBANFFNormNCcoh
            0,                    // kRHCNuisBANFFNormNCoth
            0,                    // kRHCNuisBANFFNormNuEToNuMu
#ifdef MECBAR_ERROR
            0,                      //    kRHCNuisBANFFMECBar,
#endif
#endif
        };
        ///
        
        
        // NEW FLUX FHC+RHC (invert *Nu* w/ *Nu*Bar for RHC mode)
        
        static const int kNuPdgCodeOfFHCpRHCFluxParamInBANFFfit[kNumOfFHCpRHCParamsInBANFFfit] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
#if (BANFF_BINS == 2)
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
            
            // NEW SPLIT ANA
            //#endif
            /* #if BANFF_BINS == 3 */
            /* kPdgNuMu,             // kNuisBANFF00 */
            /* kPdgNuMu,             // kNuisBANFF01 */
            /* kPdgNuMu,             // kNuisBANFF02 */
            /* kPdgNuMu,             // kNuisBANFF03 */
            /* kPdgNuMu,             // kNuisBANFF04 */
            /* kPdgNuMu,             // kNuisBANFF05 */
            /* kPdgNuMu,             // kNuisBANFF06 */
            /* kPdgNuMu,             // kNuisBANFF07 */
            /* kPdgNuMu,             // kNuisBANFF08 */
            /* kPdgNuMu,             // kNuisBANFF09 */
            /* kPdgNuMu,             // kNuisBANFF10 */
            /* kPdgNuMuBar,          // kNuisBANFF11 */
            /* kPdgNuMuBar,          // kNuisBANFF12 */
            /* kPdgNuE,              // kNuisBANFF13 */
            /* kPdgNuE,              // kNuisBANFF14 */
            /* kPdgNuE,              // kNuisBANFF15 */
            /* kPdgNuE,              // kNuisBANFF16 */
            /* kPdgNuE,              // kNuisBANFF17 */
            /* kPdgNuE,              // kNuisBANFF18 */
            /* kPdgNuE,              // kNuisBANFF19 */
            /* kPdgNuEBar,           // kNuisBANFF20 */
            /* kPdgNuEBar,           // kNuisBANFF21 */
            /* #endif */
            //#if BANFF_BINS == 1
            
#else
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
            //
#if (BANFF_BINS == 2)
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
            
            // NEW SPLIT ANA
            //#endif
            /* #if BANFF_BINS == 3 */
            /* kPdgNuMu,             // kRHCNuisBANFF00 */
            /* kPdgNuMu,             // kRHCNuisBANFF01 */
            /* kPdgNuMu,             // kRHCNuisBANFF02 */
            /* kPdgNuMu,             // kRHCNuisBANFF03 */
            /* kPdgNuMu,             // kRHCNuisBANFF04 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF05 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF06 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF07 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF08 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF09 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF10 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF11 */
            /* kPdgNuMuBar,          // kRHCNuisBANFF12 */
            /* kPdgNuE,              // kRHCNuisBANFF13 */
            /* kPdgNuE,              // kRHCNuisBANFF14 */
            /* kPdgNuEBar,           // kRHCNuisBANFF15 */
            /* kPdgNuEBar,           // kRHCNuisBANFF16 */
            /* kPdgNuEBar,           // kRHCNuisBANFF17 */
            /* kPdgNuEBar,           // kRHCNuisBANFF18 */
            /* kPdgNuEBar,           // kRHCNuisBANFF19 */
            /* kPdgNuEBar,           // kRHCNuisBANFF20 */
            /* kPdgNuEBar,           // kRHCNuisBANFF21 */
            /* #endif */
            //#if BANFF_BINS == 1
            
#else
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
#endif
            
#elif defined(__USE_NIWG2015__)
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
            0,                    // kRHCNuisBANFFMEC
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
            0,                    // kRHCNuisBANFFBindE
            0,                    // kRHCNuisBANFFNormCCcoh
            0,                    // kRHCNuisBANFFNormNCcoh
            0,                    // kRHCNuisBANFFNormNCoth
            0,                    // kRHCNuisBANFFNormNuEToNuMu
#ifdef MECBAR_ERROR
            0,                     //    kRHCNuisBANFFMECBar
#endif
#endif
        };
        ///
        
        
        
        // Etrue energy range for flux-related BANFF parameters
        static const double kEtrueRangeOfFluxParamInBANFFfit[kNumOfParamsInBANFFfit][2] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
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
            
            // NEW SPLIT ANA
            //#endif
            /* #if BANFF_BINS == 3 */
            /* { 0.0,  0.4 },        // kNuisBANFF00 */
            /* { 0.4,  0.5 },        // kNuisBANFF01 */
            /* { 0.5,  0.6 },        // kNuisBANFF02 */
            /* { 0.6,  0.7 },        // kNuisBANFF03 */
            /* { 0.7,  1.0 },        // kNuisBANFF04 */
            /* { 1.0,  1.5 },        // kNuisBANFF05 */
            /* { 1.5,  2.5 },        // kNuisBANFF06 */
            /* { 2.5,  3.5 },        // kNuisBANFF07 */
            /* { 3.5,  5.0 },        // kNuisBANFF08 */
            /* { 5.0,  7.0 },        // kNuisBANFF09 */
            /* { 7.0, 30.0 },        // kNuisBANFF10 */
            /* { 0.0,  1.5 },        // kNuisBANFF11 */
            /* { 1.5, 30.0 },        // kNuisBANFF12 */
            /* { 0.0,  0.5 },        // kNuisBANFF13 */
            /* { 0.5,  0.7 },        // kNuisBANFF14 */
            /* { 0.7,  0.8 },        // kNuisBANFF15 */
            /* { 0.8,  1.5 },        // kNuisBANFF16 */
            /* { 1.5,  2.5 },        // kNuisBANFF17 */
            /* { 2.5,  4.0 },        // kNuisBANFF18 */
            /* { 4.0, 30.0 },        // kNuisBANFF19 */
            /* { 0.0,  2.5 },        // kNuisBANFF20 */
            /* { 2.5, 30.0 },        // kNuisBANFF21 */
            /* #endif */
            //#if BANFF_BINS == 1
            
#else
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
#if (BANFF_BINS != 2)
            {  -1,   -1 },         // kNuisBANFFNC1pi0, not a flux param
#endif
            
#elif defined(__USE_NIWG2015__)
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
            {  -1,   -1 },        // kNuisBANFFMEC not a flux param
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
            {  -1,   -1 },        // kNuisBANFFBindE, not a flux param
            {  -1,   -1 },        // kNuisBANFFNormCCcoh, not a flux param
            {  -1,   -1 },        // kNuisBANFFNormNCcoh, not a flux param
            {  -1,   -1 },        // kNuisBANFFNormNCoth, not a flux param
            {  -1,   -1 },        // kNuisBANFFNormNuEToNuMu, not a flux param
#ifdef MECBAR_ERROR
            {  -1,   -1 },        // kNuisBANFFMECBar, not a flux param
#endif
#endif
        };
        
        
        // NEW FLUX
        
        static const double kEtrueRangeOfRHCFluxParamInBANFFfit[kNumOfRHCParamsInBANFFfit][2] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
#if BANFF_BINS == 2
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
            
            // NEW SPLIT ANA
            //#endif
            /* #if BANFF_BINS == 3 */
            /* // numu */
            /* { 0.0,  1.5 },        // kRHCNuisBANFF00 */
            /* { 1.5, 30.0 },        // kRHCNuisBANFF01 */
            /* // anti-numu */
            /* { 0.0,  0.4 },        // kRHCNuisBANFF02 */
            /* { 0.4,  0.5 },        // kRHCNuisBANFF03 */
            /* { 0.5,  0.6 },        // kRHCNuisBANFF04 */
            /* { 0.6,  0.7 },        // kRHCNuisBANFF05 */
            /* { 0.7,  1.0 },        // kRHCNuisBANFF06 */
            /* { 1.0,  1.5 },        // kRHCNuisBANFF07 */
            /* { 1.5,  2.5 },        // kRHCNuisBANFF08 */
            /* { 2.5,  3.5 },        // kRHCNuisBANFF09 */
            /* { 3.5,  5.0 },        // kRHCNuisBANFF10 */
            /* { 5.0,  7.0 },        // kRHCNuisBANFF11 */
            /* { 7.0, 30.0 },        // kRHCNuisBANFF12 */
            /* // nue */
            /* { 0.0,  2.5 },        // kRHCNuisBANFF13 */
            /* { 2.5, 30.0 },        // kRHCNuisBANFF14 */
            /* // anti-nue */
            /* { 0.0,  0.5 },        // kRHCNuisBANFF15 */
            /* { 0.5,  0.7 },        // kRHCNuisBANFF16 */
            /* { 0.7,  0.8 },        // kRHCNuisBANFF17 */
            /* { 0.8,  1.5 },        // kRHCNuisBANFF18 */
            /* { 1.5,  2.5 },        // kRHCNuisBANFF19 */
            /* { 2.5,  4.0 },        // kRHCNuisBANFF20 */
            /* { 4.0, 30.0 },        // kRHCNuisBANFF21 */
            /* #endif */
            //#if BANFF_BINS == 1
            
#else
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
#if (BANFF_BINS != 2)
            {  -1,   -1 },         // kRHCNuisBANFFNC1pi0, not a flux param
#endif
            
#elif defined(__USE_NIWG2015__)
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
            {  -1,   -1 },        // kRHCNuisBANFFMEC not a flux param
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
            {  -1,   -1 },        // kRHCNuisBANFFBindE, not a flux param
            {  -1,   -1 },        // kRHCNuisBANFFNormCCcoh, not a flux param
            {  -1,   -1 },        // kRHCNuisBANFFNormNCcoh, not a flux param
            {  -1,   -1 },        // kRHCNuisBANFFNormNCoth, not a flux param
            {  -1,   -1 },        // kRHCNuisBANFFNormNuEToNuMu, not a flux param
#ifdef MECBAR_ERROR
            {  -1,   -1 },        // kRHCNuisMECBar, not a flux param
#endif
#endif
        };
        
        
        // NEW FLUX FHC+RHC
        // NEW FLUX
        
        static const double kEtrueRangeOfFHCpRHCFluxParamInBANFFfit[kNumOfFHCpRHCParamsInBANFFfit][2] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            // FHC
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
            
            // NEW SPLIT ANA
            //#endif	
            /* #if BANFF_BINS == 3 */
            /* { 0.0,  0.4 },        // kNuisBANFF00 */
            /* { 0.4,  0.5 },        // kNuisBANFF01 */
            /* { 0.5,  0.6 },        // kNuisBANFF02 */
            /* { 0.6,  0.7 },        // kNuisBANFF03 */
            /* { 0.7,  1.0 },        // kNuisBANFF04 */
            /* { 1.0,  1.5 },        // kNuisBANFF05 */
            /* { 1.5,  2.5 },        // kNuisBANFF06 */
            /* { 2.5,  3.5 },        // kNuisBANFF07 */
            /* { 3.5,  5.0 },        // kNuisBANFF08 */
            /* { 5.0,  7.0 },        // kNuisBANFF09 */
            /* { 7.0, 30.0 },        // kNuisBANFF10 */
            /* { 0.0,  1.5 },        // kNuisBANFF11 */
            /* { 1.5, 30.0 },        // kNuisBANFF12 */
            /* { 0.0,  0.5 },        // kNuisBANFF13 */
            /* { 0.5,  0.7 },        // kNuisBANFF14 */
            /* { 0.7,  0.8 },        // kNuisBANFF15 */
            /* { 0.8,  1.5 },        // kNuisBANFF16 */
            /* { 1.5,  2.5 },        // kNuisBANFF17 */
            /* { 2.5,  4.0 },        // kNuisBANFF18 */
            /* { 4.0, 30.0 },        // kNuisBANFF19 */
            /* { 0.0,  2.5 },        // kNuisBANFF20 */
            /* { 2.5, 30.0 },        // kNuisBANFF21 */
            /* #endif */
            //#if BANFF_BINS == 1
            
#else 
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
            
            // RHC
#if BANFF_BINS == 2
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
            //#endif
            
            // NEW SPLIT ANA
            /* #if BANFF_BINS == 3 */
            /* // numu */
            /* { 0.0,  1.5 },        // kRHCNuisBANFF00 */
            /* { 1.5, 30.0 },        // kRHCNuisBANFF01 */
            /* // anti-numu */
            /* { 0.0,  0.4 },        // kRHCNuisBANFF02 */
            /* { 0.4,  0.5 },        // kRHCNuisBANFF03 */
            /* { 0.5,  0.6 },        // kRHCNuisBANFF04 */
            /* { 0.6,  0.7 },        // kRHCNuisBANFF05 */
            /* { 0.7,  1.0 },        // kRHCNuisBANFF06 */
            /* { 1.0,  1.5 },        // kRHCNuisBANFF07 */
            /* { 1.5,  2.5 },        // kRHCNuisBANFF08 */
            /* { 2.5,  3.5 },        // kRHCNuisBANFF09 */
            /* { 3.5,  5.0 },        // kRHCNuisBANFF10 */
            /* { 5.0,  7.0 },        // kRHCNuisBANFF11 */
            /* { 7.0, 30.0 },        // kRHCNuisBANFF12 */
            /* // nue */
            /* { 0.0,  2.5 },        // kRHCNuisBANFF13 */
            /* { 2.5, 30.0 },        // kRHCNuisBANFF14 */
            /* // anti-nue */
            /* { 0.0,  0.5 },        // kRHCNuisBANFF15 */
            /* { 0.5,  0.7 },        // kRHCNuisBANFF16 */
            /* { 0.7,  0.8 },        // kRHCNuisBANFF17 */
            /* { 0.8,  1.5 },        // kRHCNuisBANFF18 */
            /* { 1.5,  2.5 },        // kRHCNuisBANFF19 */
            /* { 2.5,  4.0 },        // kRHCNuisBANFF20 */
            /* { 4.0, 30.0 },        // kRHCNuisBANFF21 */
            /* #endif */
            //#if BANFF_BINS == 1
            
#else 
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
#if (BANFF_BINS != 2)
            {  -1,   -1 },         // kNuisBANFFNC1pi0, not a flux param
#endif
            
#elif defined(__USE_NIWG2015__)
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
            {  -1,   -1 },        // kNuisBANFFMEC not a flux param
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
            {  -1,   -1 },        // kNuisBANFFBindE, not a flux param
            {  -1,   -1 },        // kNuisBANFFNormCCcoh, not a flux param
            {  -1,   -1 },        // kNuisBANFFNormNCcoh, not a flux param
            {  -1,   -1 },        // kNuisBANFFNormNCoth, not a flux param
            {  -1,   -1 },        // kNuisBANFFNormNuEToNuMu, not a flux param
#ifdef MECBAR_ERROR
            {  -1,   -1 },        // kNuisBANFFMECBar, not a flux param
#endif
#endif
        };
        //
        
    } // conv
} // valor

#endif 

