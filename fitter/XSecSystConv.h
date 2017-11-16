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
 
 */
//____________________________________________________________________________

#ifndef xsec_syst_conv_h
#define xsec_syst_conv_h

#include "AnalysisOptions.h"
#include "Enumerations.h"
#include "SamplesAndModes.h"

namespace valor {
  namespace conv   {
        
    //
    // ** Cross-section systematics
    //
        
    // Nominal values of cross-section related params
    static const double kNomMAQE  = 1.21;
    static const double kNomMARES = 1.21;
        
    // Number of cross-section related systematics
        
    // NEW XSEC NIWG
      
#if defined(__USE_NIWG2017__)
    // NEW XSEC
#  if BANFF_BINS == 1
    static const int kNumOfXSecSyst = 31;
    static const int kNumOfRHCXSecSyst = 31;
#  endif
#  ifdef __NUMUONLY__
    static const int kNumOfUncXSecSyst    = 0; // no uncorr xsec
    static const int kNumOfRHCUncXSecSyst = 0; // no uncorr xsec
#  elif __NUEONLY__
    static const int kNumOfUncXSecSyst    = 0; // no uncorr xsec
    static const int kNumOfRHCUncXSecSyst = 0; // no uncorr xsec
#  else
    static const int kNumOfUncXSecSyst    = 0; // no uncorr xsec
    static const int kNumOfRHCUncXSecSyst = 0; // no uncorr xsec
#  endif      

#endif
      

        //
        
        
        // The following is the list of the above kNumOfXSecSyst parameters
    static const PhysParamId_t kListOfXSecSyst[kNumOfXSecSyst] =
    {
            // NEW XSEC NIWG
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )
      // NEW XSEC
      kNuisXSecRPA,               // RPA tuning
        
      kNuisBANFFNorm2p2h,              // BANFF, 2p2h normalization
      kNuisBANFFCA5,              // BANFF, CA5 RES
      kNuisBANFFBgRES,            // BANFF, Isospin 1/2 normalization
      kNuisBANFFMAQE,             // BANFF, M_{A}  QE
      kNuisBANFFMARES,            // BANFF, M_{A}  RES
      // NEW SPLIT ANA
#  if (BANFF_BINS != 2)
      kNuisBANFFSCCA,             // BANFF, SCCA
      kNuisBANFFSCCV,             // BANFF, SCCV
#  endif
      kNuisBANFFFermiMom,         // BANFF, Fermi momentum
      kNuisBANFFShapeCCoth,       // BANFF, CC other shape
      kNuisBANFFNormCCcoh,        // BANFF, CC coherent norm
      kNuisBANFFNormNCcoh,        // BANFF, NC coh norm
      kNuisBANFFNormNCoth,        // BANFF, NC other norm
      kNuisBANFFNormNuEToNuMu,    // BANFF, nue / numu cross section
      kNuisBANFFNormNC1Gamma,		    // BANFF, NC 1 gamma
      kNuisBANFFNormNuEBarToNuMuBar,
      kNuisBANFFNorm2p2hBar,              // BANFF, 2p2hBar normalization
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
      kNuisBANFFSKFSIInelastic

#  endif             
#endif
    };
        
    static const PhysParamId_t kListOfRHCXSecSyst[kNumOfRHCXSecSyst] =
    {
            // NEW XSEC NIWG
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )
      // NEW XSEC
      kRHCNuisXSecRPA,               // RPA tuning
            
      kRHCNuisBANFF2p2h,              // BANFF, 2p2h normalization
      kRHCNuisBANFFCA5,              // BANFF, CA5 RES
      kRHCNuisBANFFBgRES,            // BANFF, Isospin 1/2 normalization
      kRHCNuisBANFFMAQE,             // BANFF, M_{A}  QE
      kRHCNuisBANFFMARES,            // BANFF, M_{A}  RES
      // NEW SPLIT ANA
#  if (BANFF_BINS != 2)
      kRHCNuisBANFFSCCA,             // BANFF, SCCA
      kRHCNuisBANFFSCCV,             // BANFF, SCCV
#  endif
      kRHCNuisBANFFFermiMom,         // BANFF, Fermi momentum
      kRHCNuisBANFFShapeCCoth,       // BANFF, CC other shape
      kRHCNuisBANFFNormCCcoh,        // BANFF, CC coherent norm
      kRHCNuisBANFFNormNCcoh,        // BANFF, NC coh norm
      kRHCNuisBANFFNormNCoth,        // BANFF, NC other norm
      kRHCNuisBANFFNormNuEToNuMu,    // BANFF, nue / numu cross section
      kRHCNuisBANFFNormNC1Gamma,		    // BANFF, NC 1 gamma
      kRHCNuisBANFFNormNuEBarToNuMuBar,
      kRHCNuisBANFF2p2hBar,              // BANFF, 2p2hBar normalization
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
      kRHCNuisBANFFSKFSIInelastic

#  endif             
#endif
        };
        
        
        // The following list specifies which cross-section related systematic has its
        // effect on MC templates pre-computed and stored in cubic splines.
        // The order of systematic parameters should match the one in kListOfXSecSyst.
    static const bool kXSecSystHasWghtSplines[kNumOfXSecSyst] =
    {
      // NEW XSEC NIWG            
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )
      // NEW XSEC
      false,                     // RPA tuning kNuisXSecRPA
            
      false,                     // kNuisBANFFNorm2p2h
      true,                      // kNuisBANFFCA5
      true,                      // kNuisBANFFBgRES
      true,                      // kNuisBANFFMAQE
      true,                      // kNuisBANFFMARES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      true,                      // kNuisBANFFSCCA
      true,                      // kNuisBANFFSCCV
#endif
      true,                      // kNuisBANFFFermiMom
      true,                      // kNuisBANFFShapeCCoth
      false,                     // kNuisBANFFNormCCcoh
      false,                     // kNuisBANFFNormNCcoh
      false,                     // kNuisBANFFNormNCoth
      false,                     // kNuisBANFFNormNuEToNuMu
      false,                     // kNuisBANFFNormNC1Gamma
      false,                     // kNuisBANFFNormNuEBarToNuMuBar
      false,                     // kNuisBANFFNorm2p2hBar
#  if defined(__USE_NIWG2017__)      
      true,                      //kNuisBANFFBeRPA_A     
      true,                      //kNuisBANFFBeRPA_B     
      true,                      //kNuisBANFFBeRPA_D     
      true,                      //kNuisBANFFBeRPA_E     
      false,                     //kNuisBANFFBeRPA_U     
      true,                      //kNuisBANFFShape2p2h   
      false,                     //kNuisBANFFNorm2p2hCtoO
      false,                     // kNuisBANFFNormCCNu,
      false,                     // kNuisBANFFNormCCNuBar,
      true,                      // kNuisBANFFSKFSIAbsorption,
      true,                      // kNuisBANFFSKFSIQELow,
      true,                      // kNuisBANFFSKFSICXLow,
      true,                      // kNuisBANFFSKFSIQEHigh,
      true                       // kNuisBANFFSKFSIInelastic

#  endif             
#endif
          };
        
    static const bool kRHCXSecSystHasWghtSplines[kNumOfRHCXSecSyst] =
    {
            // NEW XSEC NIWG
#if ( defined(__USE_NIWG2015__) || defined (__USE_NIWG2017__) )
      // NEW XSEC
      false,                      // RPA tuning kRHCNuisXSecRPA
            
      false,                     // kRHCNuisBANFF2p2h
      true,                      // kRHCNuisBANFFCA5
      true,                      // kRHCNuisBANFFBgRES
      true,                      // kRHCNuisBANFFMAQE
      true,                      // kRHCNuisBANFFMARES
      // NEW SPLIT ANA
#  if (BANFF_BINS != 2)
      true,                      // kRHCNuisBANFFSCCA
      true,                      // kRHCNuisBANFFSCCV
#  endif
      true,                      // kRHCNuisBANFFFermiMom
      true,                      // kRHCNuisBANFFShapeCCoth
      false,                     // kRHCNuisBANFFNormCCcoh
      false,                     // kRHCNuisBANFFNormNCcoh
      false,                     // kRHCNuisBANFFNormNCoth
      false,                     // kRHCNuisBANFFNormNuEToNuMu
      false,                     // kRHCNuisBANFFNormNC1Gamma
      false,                     // kRHCNuisBANFFNormNuEBarToNuMuBar
      false,                     // kRHCNuisBANFF2p2hBar
#  if defined(__USE_NIWG2017__)      
      true,                      // kRHCNuisBANFFBeRPA_A     
      true,                      // kRHCNuisBANFFBeRPA_B     
      true,                      // kRHCNuisBANFFBeRPA_D     
      true,                      // kRHCNuisBANFFBeRPA_E     
      false,                     // kRHCNuisBANFFBeRPA_U     
      true,                      // kRHCNuisBANFFShape2p2h   
      false,                     // kRHCNuisBANFFNorm2p2hCtoO
      false,                     // kNuisBANFFNormCCNu,
      false,                     // kNuisBANFFNormCCNuBar,
      true,                      // kNuisBANFFSKFSIAbsorption,
      true,                      // kNuisBANFFSKFSIQELow,
      true,                      // kNuisBANFFSKFSICXLow,
      true,                      // kNuisBANFFSKFSIQEHigh,
      true                       // kNuisBANFFSKFSIInelastic

#  endif             

#endif
    };
        
    // The following list specifies which cross-section related systematic
    // parameters were included in the BANFF fit.
    // The order of systematic parameters should match the one in kListOfXSecSyst.
    static const bool kXSecSystIncludedInBANFF[kNumOfXSecSyst] =
    {
            // NEW XSEC
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )
      // NEW XSEC
      false,                     // RPA tuning kNuisXSecRPA
      
      true,                      // kNuisBANFFNorm2p2h
      true,                      // kNuisBANFFCA5
      true,                      // kNuisBANFFBgRES
      true,                      // kNuisBANFFMAQE
      true,                      // kNuisBANFFMARES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      true,                      // kNuisBANFFSCCA
      true,                      // kNuisBANFFSCCV
#endif
      true,                      // kNuisBANFFFermiMom
      true,                      // kNuisBANFFShapeCCoth
      true,                      // kNuisBANFFNormCCcoh
      true,                      // kNuisBANFFNormNCcoh
      true,                      // kNuisBANFFNormNCoth
      true,                      // kNuisBANFFNormNuEToNuMu
      true,                      // kNuisBANFFNormNC1Gamma
      true,                      // kNuisBANFFNormNuEBarToNuMuBar
      true,                      // kNuisBANFFNorm2p2hBar
#  if defined(__USE_NIWG2017__)      
      true,                      // kNuisBANFFBeRPA_A     
      true,                      // kNuisBANFFBeRPA_B     
      true,                      // kNuisBANFFBeRPA_D     
      true,                      // kNuisBANFFBeRPA_E     
      true,                      // kNuisBANFFBeRPA_U     
      true,                      // kNuisBANFFShape2p2h   
      true,                      // kNuisBANFFNorm2p2hCtoO
      true,                      // kNuisBANFFNormCCNu,
      true,                      // kNuisBANFFNormCCNuBar,
      true,                      // kNuisBANFFSKFSIAbsorption,
      true,                      // kNuisBANFFSKFSIQELow,
      true,                      // kNuisBANFFSKFSICXLow,
      true,                      // kNuisBANFFSKFSIQEHigh,
      true                       // kNuisBANFFSKFSIInelastic

#  endif             

#endif
    };
        
        
    static const bool kRHCXSecSystIncludedInBANFF[kNumOfRHCXSecSyst] =
    {
            // NEW XSEC NIWG
#if ( defined(__USE_NIWG2015__)|| defined(__USE_NIWG2017__) )
      // NEW XSEC
      false,                     // RPA tuning kRHCNuisXSecRPA
            
      true,                      // kRHCNuisBANFF2p2h
      true,                      // kRHCNuisBANFFCA5
      true,                      // kRHCNuisBANFFBgRES
      true,                      // kRHCNuisBANFFMAQE
      true,                      // kRHCNuisBANFFMARES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      true,                      // kRHCNuisBANFFSCCA
      true,                      // kRHCNuisBANFFSCCV
#endif
      true,                      // kRHCNuisBANFFFermiMom
      true,                      // kRHCNuisBANFFShapeCCoth
      true,                      // kRHCNuisBANFFNormCCcoh
      true,                      // kRHCNuisBANFFNormNCcoh
      true,                      // kRHCNuisBANFFNormNCoth
      true,                      // kRHCNuisBANFFNormNuEToNuMu
      true,                      // kRHCNuisBANFFNormNC1Gamma
      true,                      // kRHCNuisBANFFNormNuEBarToNuMuBar
      true,                      // kRHCNuisBANFF2p2hBar
#  if defined(__USE_NIWG2017__)      
      true,                      //kNuisBANFFBeRPA_A     
      true,                      //kNuisBANFFBeRPA_B     
      true,                      //kNuisBANFFBeRPA_D     
      true,                      //kNuisBANFFBeRPA_E     
      true,                      //kNuisBANFFBeRPA_U     
      true,                      //kNuisBANFFShape2p2h   
      true,                      //kNuisBANFFNorm2p2hCtoO
      true,                      // kRHCNuisBANFFNormCCNu,
      true,                      // kRHCNuisBANFFNormCCNuBar,
      true,                      // kRHCNuisBANFFSKFSIAbsorption,
      true,                      // kRHCNuisBANFFSKFSIQELow,
      true,                      // kRHCNuisBANFFSKFSICXLow,
      true,                      // kRHCNuisBANFFSKFSIQEHigh,
      true                       // kRHCNuisBANFFSKFSIInelastic

#  endif             

#endif
        };
        
        // Ereco and Etrue cuts for the application of cross-section systematics.
        // The order of systematic parameters should match the one in kListOfXSecSyst.
    static const double kXSecSystErecoCut[kNumOfXSecSyst][2] =
    {
      // NEW XSEC NIWG
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )
      // NEW XSEC
      { 0.0, 100.0 },            // RPA tuning kNuisXSec RPA
            
      { 0.0, 100.0 },            // kNuisBANFFNorm2p2h
      { 0.0, 100.0 },            // kNuisBANFFCA5
      { 0.0, 100.0 },            // kNuisBANFFBgRES
      { 0.0, 100.0 },            // kNuisBANFFMAQE
      { 0.0, 100.0 },            // kNuisBANFFMARES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      { 0.0, 100.0 },            // kNuisBANFFSCCA
      { 0.0, 100.0 },            // kNuisBANFFSCCV
#endif
      { 0.0, 100.0 },            // kNuisBANFFFermiMom
      { 0.0, 100.0 },            // kNuisBANFFShapeCCoth
      { 0.0, 100.0 },            // kNuisBANFFNormCCcoh
      { 0.0, 100.0 },            // kNuisBANFFNormNCcoh
      { 0.0, 100.0 },            // kNuisBANFFNormNCoth
      { 0.0, 100.0 },            // kNuisBANFFNormNuEToNuMu
      { 0.0, 100.0 },            // kNuisBANFFNormNC1Gamma
      { 0.0, 100.0 },            // kNuisBANFFNormNuEBarToNuMuBar
      { 0.0, 100.0 },            // kNuisBANFFNorm2p2hBar
#  if defined(__USE_NIWG2017__)      
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_A     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_B     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_D     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_E     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_U     
      { 0.0, 100.0 },            //kNuisBANFFShape2p2h   
      { 0.0, 100.0 },            //kNuisBANFFNorm2p2hCtoO
      { 0.0, 100.0 },            // kNuisBANFFNormCCNu,
      { 0.0, 100.0 },            // kNuisBANFFNormCCNuBar,
      { 0.0, 100.0 },            // kNuisBANFFSKFSIAbsorption,
      { 0.0, 100.0 },            // kNuisBANFFSKFSIQELow,
      { 0.0, 100.0 },            // kNuisBANFFSKFSICXLow,
      { 0.0, 100.0 },            // kNuisBANFFSKFSIQEHigh,
      { 0.0, 100.0 }             // kNuisBANFFSKFSIInelastic
#  endif             

#endif
    };
        
    static const double kRHCXSecSystErecoCut[kNumOfRHCXSecSyst][2] =
    {
      // NEW XSEC NIWG
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )
      // NEW XSEC
      { 0.0, 100.0 },            // RPA tuning kRHCNuisXSecRPA
            
      { 0.0, 100.0 },            // kRHCNuisBANFF2p2h
      { 0.0, 100.0 },            // kRHCNuisBANFFCA5
      { 0.0, 100.0 },            // kRHCNuisBANFFBgRES
      { 0.0, 100.0 },            // kRHCNuisBANFFMAQE
      { 0.0, 100.0 },            // kRHCNuisBANFFMARES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      { 0.0, 100.0 },            // kRHCNuisBANFFSCCA
      { 0.0, 100.0 },            // kRHCNuisBANFFSCCV
#endif
      { 0.0, 100.0 },            // kRHCNuisBANFFFermiMom
      { 0.0, 100.0 },            // kRHCNuisBANFFShapeCCoth
      { 0.0, 100.0 },            // kRHCNuisBANFFNormCCcoh
      { 0.0, 100.0 },            // kRHCNuisBANFFNormNCcoh
      { 0.0, 100.0 },            // kRHCNuisBANFFNormNCoth
      { 0.0, 100.0 },            // kRHCNuisBANFFNormNuEToNuMu
      { 0.0, 100.0 },            // kRHCNuisBANFFNormNC1Gamma
      { 0.0, 100.0 },            // kRHCNuisBANFFNormNuEBarToNuMuBar
      { 0.0, 100.0 },            // kRHCNuisBANFF2p2hBar
#  if defined(__USE_NIWG2017__)      
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_A     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_B     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_D     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_E     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_U     
      { 0.0, 100.0 },            //kNuisBANFFShape2p2h   
      { 0.0, 100.0 },            //kNuisBANFFNorm2p2hCtoO
      { 0.0, 100.0 },            // kRHCNuisBANFFNormCCNu,
      { 0.0, 100.0 },            // kRHCNuisBANFFNormCCNuBar,
      { 0.0, 100.0 },            // kRHCNuisBANFFSKFSIAbsorption,
      { 0.0, 100.0 },            // kRHCNuisBANFFSKFSIQELow,
      { 0.0, 100.0 },            // kRHCNuisBANFFSKFSICXLow,
      { 0.0, 100.0 },            // kRHCNuisBANFFSKFSIQEHigh,
      { 0.0, 100.0 }             // kRHCNuisBANFFSKFSIInelastic

#  endif             

#endif
    };
        
    static const double kXSecSystEtrueCut[kNumOfXSecSyst][2] =
    {
            // NEW XSEC NIWG
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )
      // NEW XSEC
      { 0.0, 100.0 },            // RPA tuning kNuisXSecRPA
            
      { 0.0, 100.0 },            // kNuisBANFFNorm2p2h
      { 0.0, 100.0 },            // kNuisBANFFCA5
      { 0.0, 100.0 },            // kNuisBANFFBgRES
      { 0.0, 100.0 },            // kNuisBANFFMAQE
      { 0.0, 100.0 },            // kNuisBANFFMARES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      { 0.0, 100.0 },            // kNuisBANFFSCCA
      { 0.0, 100.0 },            // kNuisBANFFSCCV
#endif
      { 0.0, 100.0 },            // kNuisBANFFFermiMom
      { 0.0, 100.0 },            // kNuisBANFFShapeCCoth
      { 0.0, 100.0 },            // kNuisBANFFNormCCcoh
      { 0.0, 100.0 },            // kNuisBANFFNormNCcoh
      { 0.0, 100.0 },            // kNuisBANFFNormNCoth
      { 0.0, 100.0 },            // kNuisBANFFNormNuEToNuMu
      { 0.0, 100.0 },            // kNuisBANFFNormNC1Gamma
      { 0.0, 100.0 },            // kNuisBANFFNormNuEBarToNuMuBar
      { 0.0, 100.0 },            // kNuisBANFFNorm2p2hBar
#  if defined(__USE_NIWG2017__)      
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_A     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_B     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_D     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_E     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_U     
      { 0.0, 100.0 },            //kNuisBANFFShape2p2h   
      { 0.0, 100.0 },            //kNuisBANFFNorm2p2hCtoO
      { 0.0, 100.0 },            // kNuisBANFFNormCCNu,
      { 0.0, 100.0 },            // kNuisBANFFNormCCNuBar,
      { 0.0, 100.0 },            // kNuisBANFFSKFSIAbsorption,
      { 0.0, 100.0 },            // kNuisBANFFSKFSIQELow,
      { 0.0, 100.0 },            // kNuisBANFFSKFSICXLow,
      { 0.0, 100.0 },            // kNuisBANFFSKFSIQEHigh,
      { 0.0, 100.0 }             // kNuisBANFFSKFSIInelastic

#  endif             
#endif
    };
        
    static const double kRHCXSecSystEtrueCut[kNumOfRHCXSecSyst][2] =
    {
      // NEW XSEC NIWG            
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )
      // NEW XSEC
      { 0.0, 100.0 },            // RPA tuning kRHCNuisXSecRPA
            
      { 0.0, 100.0 },            // kRHCNuisBANFF2p2h
      { 0.0, 100.0 },            // kRHCNuisBANFFCA5
      { 0.0, 100.0 },            // kRHCNuisBANFFBgRES
      { 0.0, 100.0 },            // kRHCNuisBANFFMAQE
      { 0.0, 100.0 },            // kRHCNuisBANFFMARES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      { 0.0, 100.0 },            // kRHCNuisBANFFSCCA
      { 0.0, 100.0 },            // kRHCNuisBANFFSCCV
#endif
      { 0.0, 100.0 },            // kRHCNuisBANFFFermiMom
      { 0.0, 100.0 },            // kRHCNuisBANFFShapeCCoth
      { 0.0, 100.0 },            // kRHCNuisBANFFNormCCcoh
      { 0.0, 100.0 },            // kRHCNuisBANFFNormNCcoh
      { 0.0, 100.0 },            // kRHCNuisBANFFNormNCoth
      { 0.0, 100.0 },            // kRHCNuisBANFFNormNuEToNuMu
      { 0.0, 100.0 },            // kRHCNuisBANFFNormNC1Gamma
      { 0.0, 100.0 },            // kRHCNuisBANFFNormNuEBarToNuMuBar
      { 0.0, 100.0 },            // kRHCNuisBANFF2p2hBar
#  if defined(__USE_NIWG2017__)      
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_A     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_B     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_D     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_E     
      { 0.0, 100.0 },            //kNuisBANFFBeRPA_U     
      { 0.0, 100.0 },            //kNuisBANFFShape2p2h   
      { 0.0, 100.0 },            //kNuisBANFFNorm2p2hCtoO
      { 0.0, 100.0 },            // kRHCNuisBANFFNormCCNu,
      { 0.0, 100.0 },            // kRHCNuisBANFFNormCCNuBar,
      { 0.0, 100.0 },            // kRHCNuisBANFFSKFSIAbsorption,
      { 0.0, 100.0 },            // kRHCNuisBANFFSKFSIQELow,
      { 0.0, 100.0 },            // kRHCNuisBANFFSKFSICXLow,
      { 0.0, 100.0 },            // kRHCNuisBANFFSKFSIQEHigh,
      { 0.0, 100.0 }             // kRHCNuisBANFFSKFSIInelastic

#  endif             
#endif
    };
        
        
    // The following list describles in which MC templates we are supposed to apply the
    // weights as a result of tweaking any of the above systematic parameters.
    // The order of systematic parameters should match the one in kListOfXSecSyst.
#if BANFF_BINS == 1
        
    //Todo: make this work if kNumOfRHCXSecSyst is different from kNumOfXSecSyst.
    
    // NEW XSEC
    static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {
          
    // > 1 muon-like ring
    {
      // NEW XSEC NIWG
#if defined(__USE_NIWG2017__)
      /// rpa
      /// banff2p2h banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
      /// banffpF  banffccothshape 
      /// banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
      /// banffnuebarnumubar banff2p2hbar
      /// BeRPA_A  BeRPA_B  BeRPA_D  BeRPA_E  BeRPA_U
      /// 2p2hshapeO  2p2hCtoO
      /// BANFFNormCCNu, BANFFNormCCNuBar
      /// SKFSIAbsorption, SKFSIQELow, SKFSICXLow, SKFSIQEHigh, SKFSIInelastic
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma

      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma

      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCoth
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCoth
#endif
    }//1Rmu
    ,
    // > 1 electron-like ring
    {
#if defined(__USE_NIWG2017__)
      /// rpa
      /// banff2p2h banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
      /// banffpF  banffccothshape 
      /// banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
      /// banffnuebarnumubar banff2p2hbar
      /// BeRPA_A  BeRPA_B  BeRPA_D  BeRPA_E  BeRPA_U
      /// 2p2hshapeO  2p2hCtoO
      /// BANFFNormCCNu, BANFFNormCCNuBar
      /// SKFSIAbsorption, SKFSIQELow, SKFSICXLow, SKFSIQEHigh, SKFSIInelastic
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma

      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma

      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENCoth_e
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCoth
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCoth
#endif
    }//1Re
    ,
    // > 1 muon-like ring RHC
    {
#if defined(__USE_NIWG2017__)
      /// rpa
      /// banff2p2h banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
      /// banffpF  banffccothshape 
      /// banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
      /// banffnuebarnumubar banff2p2hbar
      /// BeRPA_A  BeRPA_B  BeRPA_D  BeRPA_E  BeRPA_U
      /// 2p2hshapeO  2p2hCtoO
      /// BANFFNormCCNu, BANFFNormCCNuBar
      /// SKFSIAbsorption, SKFSIQELow, SKFSICXLow, SKFSIQEHigh, SKFSIInelastic
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma

      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma

      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENC1Gamma
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCoth
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNC1Gamma
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCoth
#endif      
    }//1Rmu RHC
    ,
    // > 1 electron-like ring RHC
    {
#if defined(__USE_NIWG2017__)
      /// rpa
      /// banff2p2h banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
      /// banffpF  banffccothshape 
      /// banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
      /// banffnuebarnumubar banff2p2hbar
      /// BeRPA_A  BeRPA_B  BeRPA_D  BeRPA_E  BeRPA_U
      /// 2p2hshapeO  2p2hCtoO
      /// BANFFNormCCNu, BANFFNormCCNuBar
      /// SKFSIAbsorption, SKFSIQELow, SKFSICXLow, SKFSIQEHigh, SKFSIInelastic
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma

      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma

      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCoth
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCoth
#endif
    }// 1Re RHC
    ,
    // > muon-like multiring
    {
#if defined(__USE_NIWG2017__)
      /// rpa
      /// banff2p2h banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
      /// banffpF  banffccothshape 
      /// banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
      /// banffnuebarnumubar banff2p2hbar
      /// BeRPA_A  BeRPA_B  BeRPA_D  BeRPA_E  BeRPA_U
      /// 2p2hshapeO  2p2hCtoO
      /// BANFFNormCCNu, BANFFNormCCNuBar
      /// SKFSIAbsorption, SKFSIQELow, SKFSICXLow, SKFSIQEHigh, SKFSIInelastic
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma

      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma

      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma
               
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCoth
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCoth
#endif
    }//MultiRmu
    ,
    // > electron-like multiring
    {                
#if defined(__USE_NIWG2017__)
      /// rpa
      /// banff2p2h banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
      /// banffpF  banffccothshape 
      /// banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
      /// banffnuebarnumubar banff2p2hbar
      /// BeRPA_A  BeRPA_B  BeRPA_D  BeRPA_E  BeRPA_U
      /// 2p2hshapeO  2p2hCtoO
      /// BANFFNormCCNu, BANFFNormCCNuBar
      /// SKFSIAbsorption, SKFSIQELow, SKFSICXLow, SKFSIQEHigh, SKFSIInelastic
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma

      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuBarNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma

      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuECCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuENCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCcoh
      { 0,  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuECCoth
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmNuEBarCCoth
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNC1piC
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNC1pi0
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuEBarNCoth
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  0, 0,  1, 1, 1, 1, 1 }, /// kSKmNuMuNC1Gamma
                
      { 0,  0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCQE
      { 0,  0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCC1pi
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCcoh
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCC2p2h
      { 0,  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1, 1,  1, 1, 1, 1, 1 }, /// kSKmOscNuEBarCCoth
#endif
    }//MultiRe
    };

#endif //BANFF_BINS
        
  } // conv
} // valor

#endif 
