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
#if defined(__USE_NIWG2015__)
    // NEW XSEC
#  if BANFF_BINS == 1
    static const int kNumOfXSecSyst = 18;
    static const int kNumOfRHCXSecSyst = 18;
#  elif BANFF_BINS == 2
    static const int kNumOfXSecSyst = 15; //15   // no scca and sccv
    static const int kNumOfRHCXSecSyst = 15; //15   // no scca and sccv
#  elif BANFF_BINS == 3
    static const int kNumOfXSecSyst = 17;
    static const int kNumOfRHCXSecSyst = 17;
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
      
#if defined(__USE_NIWG2017__)
    // NEW XSEC
#  if BANFF_BINS == 1
    static const int kNumOfXSecSyst = 25 + kNumOfExtraXSecParamsInBANFFfit;
    static const int kNumOfRHCXSecSyst = 25 + kNumOfExtraXSecParamsInBANFFfit;
    static const int kNumOfXSecSystOrig = 25;
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
      kNuisBANFFEb
#  endif             

#ifdef EXPANDEDBANFF
      ,
#ifdef NEWBANFF00
      kNuisBANFFNew00,
#endif //#NEWBANFF00
#ifdef NEWBANFF01
      kNuisBANFFNew01,
#endif //#NEWBANFF01
#ifdef NEWBANFF02
      kNuisBANFFNew02,
#endif //#NEWBANFF02
#ifdef NEWBANFF03
      kNuisBANFFNew03,
#endif //#NEWBANFF03
#ifdef NEWBANFF04
      kNuisBANFFNew04,
#endif //#NEWBANFF04
#ifdef NEWBANFF05
      kNuisBANFFNew05,
#endif //#NEWBANFF05
#ifdef NEWBANFF06
      kNuisBANFFNew06,
#endif //#NEWBANFF06
#ifdef NEWBANFF07
      kNuisBANFFNew07,
#endif //#NEWBANFF07
#ifdef NEWBANFF08
      kNuisBANFFNew08,
#endif //#NEWBANFF08
#ifdef NEWBANFF09
      kNuisBANFFNew09
#endif //#NEWBANFF09
#endif //#EXPANDEDBANFF

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
      kRHCNuisBANFFEb
#  endif             

#ifdef EXPANDEDBANFF
      ,
#ifdef NEWBANFF00
      kNuisBANFFNew00,
#endif //#NEWBANFF00
#ifdef NEWBANFF01
      kNuisBANFFNew01,
#endif //#NEWBANFF01
#ifdef NEWBANFF02
      kNuisBANFFNew02,
#endif //#NEWBANFF02
#ifdef NEWBANFF03
      kNuisBANFFNew03,
#endif //#NEWBANFF03
#ifdef NEWBANFF04
      kNuisBANFFNew04,
#endif //#NEWBANFF04
#ifdef NEWBANFF05
      kNuisBANFFNew05,
#endif //#NEWBANFF05
#ifdef NEWBANFF06
      kNuisBANFFNew06,
#endif //#NEWBANFF06
#ifdef NEWBANFF07
      kNuisBANFFNew07,
#endif //#NEWBANFF07
#ifdef NEWBANFF08
      kNuisBANFFNew08,
#endif //#NEWBANFF08
#ifdef NEWBANFF09
      kNuisBANFFNew09
#endif //#NEWBANFF09
#endif //#EXPANDEDBANFF

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
      true,                      //kNuisBANFFEb
#  endif             

#ifdef EXPANDEDBANFF
      ,
#ifdef NEWBANFF00
      false,       //BANFFNew00,
#endif //#NEWBANFF00
#ifdef NEWBANFF01
      false,       //BANFFNew01,
#endif //#NEWBANFF01
#ifdef NEWBANFF02
      false,       //BANFFNew02,
#endif //#NEWBANFF02
#ifdef NEWBANFF03
      false,       //BANFFNew03,
#endif //#NEWBANFF03
#ifdef NEWBANFF04
      false,       //BANFFNew04,
#endif //#NEWBANFF04
#ifdef NEWBANFF05
      false,       //BANFFNew05,
#endif //#NEWBANFF05
#ifdef NEWBANFF06
      false,       //BANFFNew06,
#endif //#NEWBANFF06
#ifdef NEWBANFF07
      false,       //BANFFNew07,
#endif //#NEWBANFF07
#ifdef NEWBANFF08
      false,       //BANFFNew08,
#endif //#NEWBANFF08
#ifdef NEWBANFF09
      false,       //BANFFNew09
#endif //#NEWBANFF09
#endif //#EXPANDEDBANFF
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
      true,                      // kNuisBANFFBeRPA_A     
      true,                      // kNuisBANFFBeRPA_B     
      true,                      // kNuisBANFFBeRPA_D     
      true,                      // kNuisBANFFBeRPA_E     
      false,                     // kNuisBANFFBeRPA_U     
      true,                      // kNuisBANFFShape2p2h   
      false,                      // kNuisBANFFNorm2p2hCtoO
      true,                      // kNuisBANFFEb
#  endif             

#ifdef EXPANDEDBANFF
      ,
#ifdef NEWBANFF00
      false,       //BANFFNew00,
#endif //#NEWBANFF00
#ifdef NEWBANFF01
      false,       //BANFFNew01,
#endif //#NEWBANFF01
#ifdef NEWBANFF02
      false,       //BANFFNew02,
#endif //#NEWBANFF02
#ifdef NEWBANFF03
      false,       //BANFFNew03,
#endif //#NEWBANFF03
#ifdef NEWBANFF04
      false,       //BANFFNew04,
#endif //#NEWBANFF04
#ifdef NEWBANFF05
      false,       //BANFFNew05,
#endif //#NEWBANFF05
#ifdef NEWBANFF06
      false,       //BANFFNew06,
#endif //#NEWBANFF06
#ifdef NEWBANFF07
      false,       //BANFFNew07,
#endif //#NEWBANFF07
#ifdef NEWBANFF08
      false,       //BANFFNew08,
#endif //#NEWBANFF08
#ifdef NEWBANFF09
      false,       //BANFFNew09
#endif //#NEWBANFF09
#endif //#EXPANDEDBANFF

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
      true,                      // kNuisBANFFEb
#  endif             

#ifdef EXPANDEDBANFF
      ,
#ifdef NEWBANFF00
      true,       //BANFFNew00,
#endif //#NEWBANFF00
#ifdef NEWBANFF01
      true,       //BANFFNew01,
#endif //#NEWBANFF01
#ifdef NEWBANFF02
      true,       //BANFFNew02,
#endif //#NEWBANFF02
#ifdef NEWBANFF03
      true,       //BANFFNew03,
#endif //#NEWBANFF03
#ifdef NEWBANFF04
      true,       //BANFFNew04,
#endif //#NEWBANFF04
#ifdef NEWBANFF05
      true,       //BANFFNew05,
#endif //#NEWBANFF05
#ifdef NEWBANFF06
      true,       //BANFFNew06,
#endif //#NEWBANFF06
#ifdef NEWBANFF07
      true,       //BANFFNew07,
#endif //#NEWBANFF07
#ifdef NEWBANFF08
      true,       //BANFFNew08,
#endif //#NEWBANFF08
#ifdef NEWBANFF09
      true,       //BANFFNew09
#endif //#NEWBANFF09
#endif //#EXPANDEDBANFF

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
      true,                      //kNuisBANFFEb
#  endif             

#ifdef EXPANDEDBANFF
      ,
#ifdef NEWBANFF00
      true,       //BANFFNew00,
#endif //#NEWBANFF00
#ifdef NEWBANFF01
      true,       //BANFFNew01,
#endif //#NEWBANFF01
#ifdef NEWBANFF02
      true,       //BANFFNew02,
#endif //#NEWBANFF02
#ifdef NEWBANFF03
      true,       //BANFFNew03,
#endif //#NEWBANFF03
#ifdef NEWBANFF04
      true,       //BANFFNew04,
#endif //#NEWBANFF04
#ifdef NEWBANFF05
      true,       //BANFFNew05,
#endif //#NEWBANFF05
#ifdef NEWBANFF06
      true,       //BANFFNew06,
#endif //#NEWBANFF06
#ifdef NEWBANFF07
      true,       //BANFFNew07,
#endif //#NEWBANFF07
#ifdef NEWBANFF08
      true,       //BANFFNew08,
#endif //#NEWBANFF08
#ifdef NEWBANFF09
      true        //BANFFNew09
#endif //#NEWBANFF09
#endif //#EXPANDEDBANFF

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
      { 0.0, 100.0 },            //kNuisBANFFEb
#  endif             

#ifdef EXPANDEDBANFF
      ,
#ifdef NEWBANFF00
      { 0.0, 100.0 },       //BANFFNew00,
#endif //#NEWBANFF00
#ifdef NEWBANFF01
      { 0.0, 100.0 },       //BANFFNew01,
#endif //#NEWBANFF01
#ifdef NEWBANFF02
      { 0.0, 100.0 },       //BANFFNew02,
#endif //#NEWBANFF02
#ifdef NEWBANFF03
      { 0.0, 100.0 },       //BANFFNew03,
#endif //#NEWBANFF03
#ifdef NEWBANFF04
      { 0.0, 100.0 },       //BANFFNew04,
#endif //#NEWBANFF04
#ifdef NEWBANFF05
      { 0.0, 100.0 },       //BANFFNew05,
#endif //#NEWBANFF05
#ifdef NEWBANFF06
      { 0.0, 100.0 },       //BANFFNew06,
#endif //#NEWBANFF06
#ifdef NEWBANFF07
      { 0.0, 100.0 },       //BANFFNew07,
#endif //#NEWBANFF07
#ifdef NEWBANFF08
      { 0.0, 100.0 },       //BANFFNew08,
#endif //#NEWBANFF08
#ifdef NEWBANFF09
      { 0.0, 100.0 }       //BANFFNew09
#endif //#NEWBANFF09
#endif //#EXPANDEDBANFF

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
      { 0.0, 100.0 },            //kNuisBANFFEb
#  endif             

#ifdef EXPANDEDBANFF
      ,
#ifdef NEWBANFF00
      { 0.0, 100.0 },       //BANFFNew00,
#endif //#NEWBANFF00
#ifdef NEWBANFF01
      { 0.0, 100.0 },       //BANFFNew01,
#endif //#NEWBANFF01
#ifdef NEWBANFF02
      { 0.0, 100.0 },       //BANFFNew02,
#endif //#NEWBANFF02
#ifdef NEWBANFF03
      { 0.0, 100.0 },       //BANFFNew03,
#endif //#NEWBANFF03
#ifdef NEWBANFF04
      { 0.0, 100.0 },       //BANFFNew04,
#endif //#NEWBANFF04
#ifdef NEWBANFF05
      { 0.0, 100.0 },       //BANFFNew05,
#endif //#NEWBANFF05
#ifdef NEWBANFF06
      { 0.0, 100.0 },       //BANFFNew06,
#endif //#NEWBANFF06
#ifdef NEWBANFF07
      { 0.0, 100.0 },       //BANFFNew07,
#endif //#NEWBANFF07
#ifdef NEWBANFF08
      { 0.0, 100.0 },       //BANFFNew08,
#endif //#NEWBANFF08
#ifdef NEWBANFF09
      { 0.0, 100.0 }       //BANFFNew09
#endif //#NEWBANFF09
#endif //#EXPANDEDBANFF

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
      { 0.0, 100.0 },            //kNuisBANFFEb
#  endif             

#ifdef EXPANDEDBANFF
      ,
#ifdef NEWBANFF00
      { 0.0, 100.0 },       //BANFFNew00,
#endif //#NEWBANFF00
#ifdef NEWBANFF01
      { 0.0, 100.0 },       //BANFFNew01,
#endif //#NEWBANFF01
#ifdef NEWBANFF02
      { 0.0, 100.0 },       //BANFFNew02,
#endif //#NEWBANFF02
#ifdef NEWBANFF03
      { 0.0, 100.0 },       //BANFFNew03,
#endif //#NEWBANFF03
#ifdef NEWBANFF04
      { 0.0, 100.0 },       //BANFFNew04,
#endif //#NEWBANFF04
#ifdef NEWBANFF05
      { 0.0, 100.0 },       //BANFFNew05,
#endif //#NEWBANFF05
#ifdef NEWBANFF06
      { 0.0, 100.0 },       //BANFFNew06,
#endif //#NEWBANFF06
#ifdef NEWBANFF07
      { 0.0, 100.0 },       //BANFFNew07,
#endif //#NEWBANFF07
#ifdef NEWBANFF08
      { 0.0, 100.0 },       //BANFFNew08,
#endif //#NEWBANFF08
#ifdef NEWBANFF09
      { 0.0, 100.0 }       //BANFFNew09
#endif //#NEWBANFF09
#endif //#EXPANDEDBANFF
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
      { 0.0, 100.0 },            //kNuisBANFFEb
#  endif             
#ifdef EXPANDEDBANFF
      ,
#ifdef NEWBANFF00
      { 0.0, 100.0 },       //BANFFNew00,
#endif //#NEWBANFF00
#ifdef NEWBANFF01
      { 0.0, 100.0 },       //BANFFNew01,
#endif //#NEWBANFF01
#ifdef NEWBANFF02
      { 0.0, 100.0 },       //BANFFNew02,
#endif //#NEWBANFF02
#ifdef NEWBANFF03
      { 0.0, 100.0 },       //BANFFNew03,
#endif //#NEWBANFF03
#ifdef NEWBANFF04
      { 0.0, 100.0 },       //BANFFNew04,
#endif //#NEWBANFF04
#ifdef NEWBANFF05
      { 0.0, 100.0 },       //BANFFNew05,
#endif //#NEWBANFF05
#ifdef NEWBANFF06
      { 0.0, 100.0 },       //BANFFNew06,
#endif //#NEWBANFF06
#ifdef NEWBANFF07
      { 0.0, 100.0 },       //BANFFNew07,
#endif //#NEWBANFF07
#ifdef NEWBANFF08
      { 0.0, 100.0 },       //BANFFNew08,
#endif //#NEWBANFF08
#ifdef NEWBANFF09
      { 0.0, 100.0 }       //BANFFNew09
#endif //#NEWBANFF09
#endif //#EXPANDEDBANFF

#endif
    };
        
        
    //Fractional error of systematic parameters, as of BANFF 2012v8
    // ONLY USE FOR: macros/misc/test_fitter_inputs/test_valor_kdi_splines.C
    static const double kXSecSystFracErr[kNumOfXSecSyst] =
    {
          // NEW XSEC NIWG
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )
      // NEW XSEC
      0.1,            // RPA tuning
            
      0.1,            // BANFF, 2p2h
      0.1,            // BANFF, CA5
      0.1,            // BANFF, BgRES
      0.0716,         // BANFF, M_{A}QE
      0.0675,         // BANFF, M_{A}RES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      0.1,            // BANFF, SCCA
      0.1,            // BANFF, SCCV
#endif
      0.13333,        // BANFF, Fermi momentum
      0.4,            // BANFF, CC other shape
      //0.33333,        // BANFF, Binding enery
      1.0,            // BANFF, CC coherent norm
      0.3,            // BANFF, NC coh norm
      0.3,            // BANFF, NC other norm
      0.03,           // BANFF, nue / numu cross section
      0.5,	      // BANFF, NC 1 gamma
      0.03            // BANFF, nuebar / numubar cross section
#  if defined(__USE_NIWG2017__)      
     -999,            //kNuisBANFFBeRPA_A     
     -999,            //kNuisBANFFBeRPA_B     
     -999,            //kNuisBANFFBeRPA_D     
     -999,            //kNuisBANFFBeRPA_E     
     -999,            //kNuisBANFFBeRPA_U     
     -999,            //kNuisBANFFShape2p2h   
     -999,            //kNuisBANFFNorm2p2hCtoO
     -999,            //kNuisBANFFEb
#  endif             
#ifdef EXPANDEDBANFF
      ,
#ifdef NEWBANFF00
      -999,       //BANFFNew00,
#endif //#NEWBANFF00
#ifdef NEWBANFF01
      -999,       //BANFFNew01,
#endif //#NEWBANFF01
#ifdef NEWBANFF02
      -999,       //BANFFNew02,
#endif //#NEWBANFF02
#ifdef NEWBANFF03
      -999,       //BANFFNew03,
#endif //#NEWBANFF03
#ifdef NEWBANFF04
      -999,       //BANFFNew04,
#endif //#NEWBANFF04
#ifdef NEWBANFF05
      -999,       //BANFFNew05,
#endif //#NEWBANFF05
#ifdef NEWBANFF06
      -999,       //BANFFNew06,
#endif //#NEWBANFF06
#ifdef NEWBANFF07
      -999,       //BANFFNew07,
#endif //#NEWBANFF07
#ifdef NEWBANFF08
      -999,       //BANFFNew08,
#endif //#NEWBANFF08
#ifdef NEWBANFF09
      -999       //BANFFNew09
#endif //#NEWBANFF09
#endif //#EXPANDEDBANFF

#endif
    };
        
    static const double kRHCXSecSystFracErr[kNumOfRHCXSecSyst] =
    {
      // NEW XSEC NIWG            
#if ( defined(__USE_NIWG2015__) || defined (__USE_NIWG2017__) )
      // NEW XSEC
      0.1,            // RPA tuning
            
      0.1,            // BANFF, 2p2h
      0.1,            // BANFF, CA5
      0.1,            // BANFF, BgRES
      0.0716,         // BANFF, M_{A}QE
      0.0675,         // BANFF, M_{A}RES
      // NEW SPLIT ANA
#if (BANFF_BINS != 2)
      0.1,            // BANFF, SCCA
      0.1,            // BANFF, SCCV
#endif
      0.13333,        // BANFF, Fermi momentum
      0.4,            // BANFF, CC other shape
      0.33333,        // BANFF, Binding enery
      1.0,            // BANFF, CC coherent norm
      0.3,            // BANFF, NC coh norm
      0.3,            // BANFF, NC other norm
      0.03,           // BANFF, nue / numu cross section
      0.5,	      // BANFF, NC 1 gamma
      0.03            // BANFF, nuebar / numubar cross section
#  if defined(__USE_NIWG2017__)      
      -999,            //BANFF, BeRPA_A     
      -999,            //BANFF, BeRPA_B     
      -999,            //BANFF, BeRPA_D     
      -999,            //BANFF, BeRPA_E     
      -999,            //BANFF, BeRPA_U     
      -999,            //BANFF, Shape2p2h   
      -999,            //BANFF, Norm2p2h CtoO
      -999,            //kNuisBANFFEb
#  endif             
#ifdef EXPANDEDBANFF
      ,
#ifdef NEWBANFF00
      -999,       //BANFFNew00,
#endif //#NEWBANFF00
#ifdef NEWBANFF01
      -999,       //BANFFNew01,
#endif //#NEWBANFF01
#ifdef NEWBANFF02
      -999,       //BANFFNew02,
#endif //#NEWBANFF02
#ifdef NEWBANFF03
      -999,       //BANFFNew03,
#endif //#NEWBANFF03
#ifdef NEWBANFF04
      -999,       //BANFFNew04,
#endif //#NEWBANFF04
#ifdef NEWBANFF05
      -999,       //BANFFNew05,
#endif //#NEWBANFF05
#ifdef NEWBANFF06
      -999,       //BANFFNew06,
#endif //#NEWBANFF06
#ifdef NEWBANFF07
      -999,       //BANFFNew07,
#endif //#NEWBANFF07
#ifdef NEWBANFF08
      -999,       //BANFFNew08,
#endif //#NEWBANFF08
#ifdef NEWBANFF09
      -999       //BANFFNew09
#endif //#NEWBANFF09
#endif //#EXPANDEDBANFF

#endif
    };
        
        
        
#ifdef EXPANDEDBANFF
    static const bool kXsApplyXSecSystToSKModeNewBANFF[kns][kNumOfExtraXSecParamsInBANFFfit][knm] = {
      //
      // FHC 1-ring mu-like
      //
      {
	//kSKmNuMuCCQE kSKmNuMuCC1pi kSKmNuMuCCcoh kSKmNuMuCC2p2h kSKmNuMuCCoth kSKmNuMuNC1piC kSKmNuMuNC1pi0 kSKmNuMuNCcoh kSKmNuMuNCoth kSKmNuMuNC1Gamma
	//kSKmNuMuBarCCQE kSKmNuMuBarCC1pi kSKmNuMuBarCCcoh kSKmNuMuBarCC2p2h kSKmNuMuBarCCoth kSKmNuMuBarNC1piC kSKmNuMuBarNC1pi0 kSKmNuMuBarNCcoh kSKmNuMuBarNCoth kSKmNuMuBarNC1Gamma
	//kSKmNuECCQE kSKmNuECC1pi kSKmNuECCcoh kSKmNuECC2p2h kSKmNuECCoth kSKmNuENC1piC kSKmNuENC1pi0 kSKmNuENCcoh kSKmNuENCoth kSKmNuENC1Gamma
	//kSKmSigNuECCQE kSKmSigNuECC1pi kSKmSigNuECCcoh kSKmSigNuECC2p2h kSKmSigNuECCoth
	//kSKmNuEBarCCQE kSKmNuEBarCC1pi kSKmNuEBarCCcoh kSKmNuEBarCC2p2h kSKmNuEBarCCoth kSKmNuEBarNC1piC kSKmNuEBarNC1pi0 kSKmNuEBarNCcoh kSKmNuEBarNCoth kSKmNuEBarNC1Gamma
	//kSKmSigNuEBarCCQE kSKmSigNuEBarCC1pi kSKmSigNuEBarCCcoh kSKmSigNuEBarCC2p2h kSKmSigNuEBarCCoth
#ifdef BANFF00
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF00
#endif //BANFF00
#ifdef BANFF01
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF01
#endif //BANFF01
#ifdef BANFF02
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF02
#endif //BANFF02
#ifdef BANFF03
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF03
#endif //BANFF03
#ifdef BANFF04
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF04
#endif //BANFF04
#ifdef BANFF05
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF05
#endif //BANFF05
#ifdef BANFF06
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF06
#endif //BANFF06
#ifdef BANFF07
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF07
#endif //BANFF07
#ifdef BANFF08
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF08
#endif //BANFF08
#ifdef BANFF09
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }  //BANFF09
#endif //BANFF09
      }, //// FHC 1-ring mu-like

      //
      // FHC 1-ring e-like
      //
      {
	//kSKmNuMuCCQE kSKmNuMuCC1pi kSKmNuMuCCcoh kSKmNuMuCC2p2h kSKmNuMuCCoth kSKmNuMuNC1piC kSKmNuMuNC1pi0 kSKmNuMuNCcoh kSKmNuMuNCoth kSKmNuMuNC1Gamma
	//kSKmNuMuBarCCQE kSKmNuMuBarCC1pi kSKmNuMuBarCCcoh kSKmNuMuBarCC2p2h kSKmNuMuBarCCoth kSKmNuMuBarNC1piC kSKmNuMuBarNC1pi0 kSKmNuMuBarNCcoh kSKmNuMuBarNCoth kSKmNuMuBarNC1Gamma
	//kSKmNuECCQE kSKmNuECC1pi kSKmNuECCcoh kSKmNuECC2p2h kSKmNuECCoth kSKmNuENC1piC kSKmNuENC1pi0 kSKmNuENCcoh kSKmNuENCoth kSKmNuENC1Gamma
	//kSKmSigNuECCQE kSKmSigNuECC1pi kSKmSigNuECCcoh kSKmSigNuECC2p2h kSKmSigNuECCoth
	//kSKmNuEBarCCQE kSKmNuEBarCC1pi kSKmNuEBarCCcoh kSKmNuEBarCC2p2h kSKmNuEBarCCoth kSKmNuEBarNC1piC kSKmNuEBarNC1pi0 kSKmNuEBarNCcoh kSKmNuEBarNCoth kSKmNuEBarNC1Gamma
	//kSKmSigNuEBarCCQE kSKmSigNuEBarCC1pi kSKmSigNuEBarCCcoh kSKmSigNuEBarCC2p2h kSKmSigNuEBarCCoth
#ifdef BANFF00
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF00
#endif //BANFF00
#ifdef BANFF01
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF01
#endif //BANFF01
#ifdef BANFF02
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF02
#endif //BANFF02
#ifdef BANFF03
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF03
#endif //BANFF03
#ifdef BANFF04
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF04
#endif //BANFF04
#ifdef BANFF05
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF05
#endif //BANFF05
#ifdef BANFF06
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF06
#endif //BANFF06
#ifdef BANFF07
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF07
#endif //BANFF07
#ifdef BANFF08
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF08
#endif //BANFF08
#ifdef BANFF09
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }  //BANFF09
#endif //BANFF09
      }, // FHC 1-ring e-like

      //
      // RHC 1-ring mu-like
      //
      {
	//kSKmNuMuCCQE kSKmNuMuCC1pi kSKmNuMuCCcoh kSKmNuMuCC2p2h kSKmNuMuCCoth kSKmNuMuNC1piC kSKmNuMuNC1pi0 kSKmNuMuNCcoh kSKmNuMuNCoth kSKmNuMuNC1Gamma
	//kSKmNuMuBarCCQE kSKmNuMuBarCC1pi kSKmNuMuBarCCcoh kSKmNuMuBarCC2p2h kSKmNuMuBarCCoth kSKmNuMuBarNC1piC kSKmNuMuBarNC1pi0 kSKmNuMuBarNCcoh kSKmNuMuBarNCoth kSKmNuMuBarNC1Gamma
	//kSKmNuECCQE kSKmNuECC1pi kSKmNuECCcoh kSKmNuECC2p2h kSKmNuECCoth kSKmNuENC1piC kSKmNuENC1pi0 kSKmNuENCcoh kSKmNuENCoth kSKmNuENC1Gamma
	//kSKmSigNuECCQE kSKmSigNuECC1pi kSKmSigNuECCcoh kSKmSigNuECC2p2h kSKmSigNuECCoth
	//kSKmNuEBarCCQE kSKmNuEBarCC1pi kSKmNuEBarCCcoh kSKmNuEBarCC2p2h kSKmNuEBarCCoth kSKmNuEBarNC1piC kSKmNuEBarNC1pi0 kSKmNuEBarNCcoh kSKmNuEBarNCoth kSKmNuEBarNC1Gamma
	//kSKmSigNuEBarCCQE kSKmSigNuEBarCC1pi kSKmSigNuEBarCCcoh kSKmSigNuEBarCC2p2h kSKmSigNuEBarCCoth
#ifdef BANFF00
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF00
#endif //BANFF00
#ifdef BANFF01
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF01
#endif //BANFF01
#ifdef BANFF02
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF02
#endif //BANFF02
#ifdef BANFF03
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF03
#endif //BANFF03
#ifdef BANFF04
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF04
#endif //BANFF04
#ifdef BANFF05
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF05
#endif //BANFF05
#ifdef BANFF06
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF06
#endif //BANFF06
#ifdef BANFF07
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF07
#endif //BANFF07
#ifdef BANFF08
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF08
#endif //BANFF08
#ifdef BANFF09
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }  //BANFF09
#endif //BANFF09
      }, // RHC 1-ring mu-like

      //
      // RHC 1-ring e-like
      //
      {
	//kSKmNuMuCCQE kSKmNuMuCC1pi kSKmNuMuCCcoh kSKmNuMuCC2p2h kSKmNuMuCCoth kSKmNuMuNC1piC kSKmNuMuNC1pi0 kSKmNuMuNCcoh kSKmNuMuNCoth kSKmNuMuNC1Gamma
	//kSKmNuMuBarCCQE kSKmNuMuBarCC1pi kSKmNuMuBarCCcoh kSKmNuMuBarCC2p2h kSKmNuMuBarCCoth kSKmNuMuBarNC1piC kSKmNuMuBarNC1pi0 kSKmNuMuBarNCcoh kSKmNuMuBarNCoth kSKmNuMuBarNC1Gamma
	//kSKmNuECCQE kSKmNuECC1pi kSKmNuECCcoh kSKmNuECC2p2h kSKmNuECCoth kSKmNuENC1piC kSKmNuENC1pi0 kSKmNuENCcoh kSKmNuENCoth kSKmNuENC1Gamma
	//kSKmSigNuECCQE kSKmSigNuECC1pi kSKmSigNuECCcoh kSKmSigNuECC2p2h kSKmSigNuECCoth
	//kSKmNuEBarCCQE kSKmNuEBarCC1pi kSKmNuEBarCCcoh kSKmNuEBarCC2p2h kSKmNuEBarCCoth kSKmNuEBarNC1piC kSKmNuEBarNC1pi0 kSKmNuEBarNCcoh kSKmNuEBarNCoth kSKmNuEBarNC1Gamma
	//kSKmSigNuEBarCCQE kSKmSigNuEBarCC1pi kSKmSigNuEBarCCcoh kSKmSigNuEBarCC2p2h kSKmSigNuEBarCCoth
#ifdef BANFF00
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF00
#endif //BANFF00
#ifdef BANFF01
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF01
#endif //BANFF01
#ifdef BANFF02
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF02
#endif //BANFF02
#ifdef BANFF03
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF03
#endif //BANFF03
#ifdef BANFF04
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF04
#endif //BANFF04
#ifdef BANFF05
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF05
#endif //BANFF05
#ifdef BANFF06
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF06
#endif //BANFF06
#ifdef BANFF07
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF07
#endif //BANFF07
#ifdef BANFF08
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF08
#endif //BANFF08
#ifdef BANFF09
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }  //BANFF09
#endif //BANFF09
      }, // RHC 1-ring e-like

      //
      // FHC multi-ring mu-like
      //
      {
	//kSKmNuMuCCQE kSKmNuMuCC1pi kSKmNuMuCCcoh kSKmNuMuCC2p2h kSKmNuMuCCoth kSKmNuMuNC1piC kSKmNuMuNC1pi0 kSKmNuMuNCcoh kSKmNuMuNCoth kSKmNuMuNC1Gamma
	//kSKmNuMuBarCCQE kSKmNuMuBarCC1pi kSKmNuMuBarCCcoh kSKmNuMuBarCC2p2h kSKmNuMuBarCCoth kSKmNuMuBarNC1piC kSKmNuMuBarNC1pi0 kSKmNuMuBarNCcoh kSKmNuMuBarNCoth kSKmNuMuBarNC1Gamma
	//kSKmNuECCQE kSKmNuECC1pi kSKmNuECCcoh kSKmNuECC2p2h kSKmNuECCoth kSKmNuENC1piC kSKmNuENC1pi0 kSKmNuENCcoh kSKmNuENCoth kSKmNuENC1Gamma
	//kSKmSigNuECCQE kSKmSigNuECC1pi kSKmSigNuECCcoh kSKmSigNuECC2p2h kSKmSigNuECCoth
	//kSKmNuEBarCCQE kSKmNuEBarCC1pi kSKmNuEBarCCcoh kSKmNuEBarCC2p2h kSKmNuEBarCCoth kSKmNuEBarNC1piC kSKmNuEBarNC1pi0 kSKmNuEBarNCcoh kSKmNuEBarNCoth kSKmNuEBarNC1Gamma
	//kSKmSigNuEBarCCQE kSKmSigNuEBarCC1pi kSKmSigNuEBarCCcoh kSKmSigNuEBarCC2p2h kSKmSigNuEBarCCoth
#ifdef BANFF00
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF00
#endif //BANFF00
#ifdef BANFF01
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF01
#endif //BANFF01
#ifdef BANFF02
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF02
#endif //BANFF02
#ifdef BANFF03
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF03
#endif //BANFF03
#ifdef BANFF04
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF04
#endif //BANFF04
#ifdef BANFF05
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF05
#endif //BANFF05
#ifdef BANFF06
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF06
#endif //BANFF06
#ifdef BANFF07
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF07
#endif //BANFF07
#ifdef BANFF08
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF08
#endif //BANFF08
#ifdef BANFF09
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }  //BANFF09
#endif //BANFF09
      }, //FHC multi-ring mu-like

      //
      // FHC multi-ring e-like
      //
      {
	//kSKmNuMuCCQE kSKmNuMuCC1pi kSKmNuMuCCcoh kSKmNuMuCC2p2h kSKmNuMuCCoth kSKmNuMuNC1piC kSKmNuMuNC1pi0 kSKmNuMuNCcoh kSKmNuMuNCoth kSKmNuMuNC1Gamma
	//kSKmNuMuBarCCQE kSKmNuMuBarCC1pi kSKmNuMuBarCCcoh kSKmNuMuBarCC2p2h kSKmNuMuBarCCoth kSKmNuMuBarNC1piC kSKmNuMuBarNC1pi0 kSKmNuMuBarNCcoh kSKmNuMuBarNCoth kSKmNuMuBarNC1Gamma
	//kSKmNuECCQE kSKmNuECC1pi kSKmNuECCcoh kSKmNuECC2p2h kSKmNuECCoth kSKmNuENC1piC kSKmNuENC1pi0 kSKmNuENCcoh kSKmNuENCoth kSKmNuENC1Gamma
	//kSKmSigNuECCQE kSKmSigNuECC1pi kSKmSigNuECCcoh kSKmSigNuECC2p2h kSKmSigNuECCoth
	//kSKmNuEBarCCQE kSKmNuEBarCC1pi kSKmNuEBarCCcoh kSKmNuEBarCC2p2h kSKmNuEBarCCoth kSKmNuEBarNC1piC kSKmNuEBarNC1pi0 kSKmNuEBarNCcoh kSKmNuEBarNCoth kSKmNuEBarNC1Gamma
	//kSKmSigNuEBarCCQE kSKmSigNuEBarCC1pi kSKmSigNuEBarCCcoh kSKmSigNuEBarCC2p2h kSKmSigNuEBarCCoth
#ifdef BANFF00
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF00
#endif //BANFF00
#ifdef BANFF01
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF01
#endif //BANFF01
#ifdef BANFF02
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF02
#endif //BANFF02
#ifdef BANFF03
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF03
#endif //BANFF03
#ifdef BANFF04
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF04
#endif //BANFF04
#ifdef BANFF05
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF05
#endif //BANFF05
#ifdef BANFF06
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF06
#endif //BANFF06
#ifdef BANFF07
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF07
#endif //BANFF07
#ifdef BANFF08
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }, //BANFF08
#endif //BANFF08
#ifdef BANFF09
	{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0,     1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     1, 0, 0, 0, 0 }  //BANFF09
#endif //BANFF09
      } //FHC multi-ring e-like

#endif //EXPANDEDBANFF

    // The following list describes in which MC templates we are supposed to apply the
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
      /// Eb
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuMuCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuMuCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
#  ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuMuCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuMuCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCCoth
#  endif
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuMuBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuMuBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
#  ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuMuBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuMuBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCCoth
#  endif
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuECCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuECC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuECCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuECC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECCoth
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuEBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuEBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuEBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuEBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCCoth
#  ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
#  endif
#  ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
#  endif
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
      /// Eb
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuMuCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuMuCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
#  ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuMuCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuMuCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOScNuMuCCoth
#  endif
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuMuBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuMuBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
#  ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuMuBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuMuBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCCoth
#  endif
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuECCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuECC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENCoth_e
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuECCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuECC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECCoth
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuEBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuEBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuEBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuEBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCCoth
#  ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
#  endif
#  ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
#  endif
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
      /// Eb
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuMuCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuMuCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
#  ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuMuCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuMuCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCCoth
#  endif
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuMuBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuMuBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
#  ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuMuBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuMuBar2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCCoth
#  endif
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuECCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuECC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuECCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuECC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECCoth
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuEBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuEBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuEBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuEBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCCoth
#  ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
#  endif
#  ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
#  endif
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
      /// Eb
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuMuCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuMuCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
#  ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuMuCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuMuCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOScNuMuCCoth
#  endif
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuMuBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuMuBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
#  ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuMuBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuMuBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCCoth
#  endif
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuECCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuECC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuECCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuECC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECCoth
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuEBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuEBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuEBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuEBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCCoth
#  ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
#  endif
#  ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
#  endif
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
      /// Eb
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuMuCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuMuCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
#  ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuMuCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuMuCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCCoth
#  endif
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuMuBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuMuBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
#  ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuMuBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuMuBar2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCCoth
#  endif
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuECCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuECC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuECCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuECC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECCoth
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuEBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuEBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuEBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuEBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCCoth
#  ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
#  endif
#  ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
#  endif
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
      /// Eb
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuMuCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuMuCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
#  ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuMuCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuCCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuMuCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOScNuMuCCoth
#  endif
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuMuBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuMuBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuBarNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
#  ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuMuBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuMuBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuMuBarCCoth
#endif
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuECCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuECC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuECCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuENCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuECCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECCcoh
      {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuECC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuECCoth
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmNuEBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmNuEBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarCCoth
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNC1piC
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNC1pi0
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuEBarNCoth
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmNuMuNC1Gamma
                
      {  0,   0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,  1, 1, 1, 1, 1, 0, 0,  1 }, /// kSKmOscNuEBarCCQE
      {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCC1pi
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCCcoh
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,  0, 0, 0, 0, 0, 1, 1,  1 }, /// kSKmOscNuEBarCC2p2h
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmOscNuEBarCCoth
#  ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
#  endif
#  ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
      {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0,  1 }, /// kSKmUndefined
#  endif
#endif
    }//MultiRe
    };

#endif //BANFF_BINS
        
  } // conv
} // valor

#endif 
