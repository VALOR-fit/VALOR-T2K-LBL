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
#if BANFF_BINS == 1
        static const int kNumOfXSecSyst = 18;
        static const int kNumOfRHCXSecSyst = 18;
#elif BANFF_BINS == 2
        static const int kNumOfXSecSyst = 15; //15   // no scca and sccv
        static const int kNumOfRHCXSecSyst = 15; //15   // no scca and sccv
#elif BANFF_BINS == 3
        static const int kNumOfXSecSyst = 17;
        static const int kNumOfRHCXSecSyst = 17;
#endif
#ifdef __NUEONLY__
        static const int kNumOfUncXSecSyst    = 0; // no uncorr xsec
        static const int kNumOfRHCUncXSecSyst = 0; // no uncorr xsec
#else
        static const int kNumOfUncXSecSyst    = 0; // no uncorr xsec
        static const int kNumOfRHCUncXSecSyst = 0; // no uncorr xsec
#endif
      
#elif defined(__USE_NIWG2013__)
#if BANFF_BINS == 1
        static const int kNumOfXSecSyst = 21;
        static const int kNumOfRHCXSecSyst = 21;
#elif BANFF_BINS == 2
        static const int kNumOfXSecSyst = 19;
        static const int kNumOfRHCXSecSyst = 19;
#elif BANFF_BINS == 3
        static const int kNumOfXSecSyst = 15;
        static const int kNumOfRHCXSecSyst = 15;
#endif
#ifdef __NUMUONLY__
        static const int kNumOfUncXSecSyst = 11;
        static const int kNumOfRHCUncXSecSyst = 11;
#else
#ifdef __NUEONLY__
        static const int kNumOfUncXSecSyst = 9;
        static const int kNumOfRHCUncXSecSyst = 9;
#else
        static const int kNumOfUncXSecSyst = 12;
        static const int kNumOfRHCUncXSecSyst = 12;
#endif
#endif
#endif
        //
        
        
        // The following is the list of the above kNumOfXSecSyst parameters
        static const PhysParamId_t kListOfXSecSyst[kNumOfXSecSyst] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            kNuisBANFFMAQE,             // BANFF, M_{A}  QE
            kNuisBANFFMARES,            // BANFF, M_{A}  RES
            kNuisBANFFCCQEE0,           // BANFF, CCQE   normalization, 0 - 1.5 GeV
#if BANFF_BINS != 3
            kNuisBANFFCCQEE1,           // BANFF, CCQE   normalization, 1.5 - 3.5 GeV
            kNuisBANFFCCQEE2,           // BANFF, CCQE   normalization, 3.5 - 30.0 GeV
#endif
            kNuisBANFFCC1piE0,          // BANFF, CC1pi  normalization, 0 - 2.5 GeV
#if BANFF_BINS != 3
            kNuisBANFFCC1piE1,          // BANFF, CC1pi  normalization, 2.5 - 30.0 GeV
#endif
#if (BANFF_BINS != 2)
            kNuisBANFFNC1pi0,           // BANFF, NC1pi0 normalization, 0.0 - 30.0 GeV
#endif
            kNuisXSecFermiMom,          // Fermi momentum
            kNuisXSecSpecFunc,          // spectral function on/off
            kNuisXSecWshape,            // W-shape
            kNuisXSecShapeCCoth,        // CC other shape
            kNuisXSecPiLessDelta,       // Pion-less delta decay
#ifndef __NUEONLY__
            kNuisXSecBindE,             // Binding energy
#endif
            kNuisXSecNormCCcoh,         // CC coherent norm
#ifndef __NUEONLY__
            kNuisXSecNormNC1piC,        // NC 1piC norm
#endif
#ifndef __NUMUONLY__
            kNuisXSecNormNCcoh,         // NC coh norm
#endif
            kNuisXSecNormNCoth,         // NC other norm
#ifndef __NUEONLY__
            kNuisXSecNormAntiNuToNu,    // antinu / nu cross section
#endif
            kNuisXSecNormNuEToNuMu,     // nue / numu cross section
            kNuisRHCToFHC               // RHC / FHC normalisation
            
#elif defined(__USE_NIWG2015__)
            // NEW XSEC
            kNuisXSecRPA,               // RPA tuning
            
            kNuisBANFFMEC,              // BANFF, MEC normalization
            kNuisBANFFCA5,              // BANFF, CA5 RES
            kNuisBANFFBgRES,            // BANFF, Isospin 1/2 normalization
            kNuisBANFFMAQE,             // BANFF, M_{A}  QE
            kNuisBANFFMARES,            // BANFF, M_{A}  RES
            // NEW SPLIT ANA
#if (BANFF_BINS != 2)
            kNuisBANFFSCCA,             // BANFF, SCCA
            kNuisBANFFSCCV,             // BANFF, SCCV
#endif
            kNuisBANFFFermiMom,         // BANFF, Fermi momentum
            kNuisBANFFShapeCCoth,       // BANFF, CC other shape
            kNuisBANFFBindE,            // BANFF, Binding energy
            kNuisBANFFNormCCcoh,        // BANFF, CC coherent norm
            kNuisBANFFNormNCcoh,        // BANFF, NC coh norm
            kNuisBANFFNormNCoth,        // BANFF, NC other norm
            kNuisBANFFNormNuEToNuMu,    // BANFF, nue / numu cross section
            kNuisBANFFNormNC1Gamma,		    // BANFF, NC 1 gamma
            kNuisBANFFNormNuEBarToNuMuBar,
            kNuisBANFFMECBar,              // BANFF, MECBar normalization
#endif
        };
        
        static const PhysParamId_t kListOfRHCXSecSyst[kNumOfRHCXSecSyst] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            kRHCNuisBANFFMAQE,             // BANFF, M_{A}  QE
            kRHCNuisBANFFMARES,            // BANFF, M_{A}  RES
            kRHCNuisBANFFCCQEE0,           // BANFF, CCQE   normalization, 0 - 1.5 GeV
#if BANFF_BINS != 3
            kRHCNuisBANFFCCQEE1,           // BANFF, CCQE   normalization, 1.5 - 3.5 GeV
            kRHCNuisBANFFCCQEE2,           // BANFF, CCQE   normalization, 3.5 - 30.0 GeV
#endif
            kRHCNuisBANFFCC1piE0,          // BANFF, CC1pi  normalization, 0 - 2.5 GeV
#if BANFF_BINS != 3
            kRHCNuisBANFFCC1piE1,          // BANFF, CC1pi  normalization, 2.5 - 30.0 GeV
#endif
#if (BANFF_BINS != 2)
            kRHCNuisBANFFNC1pi0,           // BANFF, NC1pi0 normalization, 0.0 - 30.0 GeV
#endif
            kRHCNuisXSecFermiMom,          // Fermi momentum
            kRHCNuisXSecSpecFunc,          // spectral function on/off
            kRHCNuisXSecWshape,            // W-shape
            kRHCNuisXSecShapeCCoth,        // CC other shape
            kRHCNuisXSecPiLessDelta,       // Pion-less delta decay
#ifndef __NUEONLY__
            kRHCNuisXSecBindE,             // Binding energy
#endif
            kRHCNuisXSecNormCCcoh,         // CC coherent norm
#ifndef __NUEONLY__
            kRHCNuisXSecNormNC1piC,        // NC 1piC norm
#endif
#ifndef __NUMUONLY__
            kRHCNuisXSecNormNCcoh,         // NC coh norm
#endif
            kRHCNuisXSecNormNCoth,         // NC other norm
#ifndef __NUEONLY__
            kRHCNuisXSecNormAntiNuToNu,    // antinu / nu cross section
#endif
            kRHCNuisXSecNormNuEToNuMu,     // nue / numu cross section
            kNuisRHCToFHC,                 // RHC / FHC normalisation
            
#elif defined(__USE_NIWG2015__)
            // NEW XSEC
            kRHCNuisXSecRPA,               // RPA tuning
            
            kRHCNuisBANFFMEC,              // BANFF, MEC normalization
            kRHCNuisBANFFCA5,              // BANFF, CA5 RES
            kRHCNuisBANFFBgRES,            // BANFF, Isospin 1/2 normalization
            kRHCNuisBANFFMAQE,             // BANFF, M_{A}  QE
            kRHCNuisBANFFMARES,            // BANFF, M_{A}  RES
            // NEW SPLIT ANA
#if (BANFF_BINS != 2)
            kRHCNuisBANFFSCCA,             // BANFF, SCCA
            kRHCNuisBANFFSCCV,             // BANFF, SCCV
#endif
            kRHCNuisBANFFFermiMom,         // BANFF, Fermi momentum
            kRHCNuisBANFFShapeCCoth,       // BANFF, CC other shape
            kRHCNuisBANFFBindE,            // BANFF, Binding energy
            kRHCNuisBANFFNormCCcoh,        // BANFF, CC coherent norm
            kRHCNuisBANFFNormNCcoh,        // BANFF, NC coh norm
            kRHCNuisBANFFNormNCoth,        // BANFF, NC other norm
            kRHCNuisBANFFNormNuEToNuMu,    // BANFF, nue / numu cross section
            kRHCNuisBANFFNormNC1Gamma,		    // BANFF, NC 1 gamma
            kRHCNuisBANFFNormNuEBarToNuMuBar,
            kRHCNuisBANFFMECBar,              // BANFF, MECBar normalization

#endif
        };
        
        
        // The following list specifies which cross-section related systematic has its
        // effect on MC templates pre-computed and stored in cubic splines.
        // The order of systematic parameters should match the one in kListOfXSecSyst.
        static const bool kXSecSystHasWghtSplines[kNumOfXSecSyst] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            true,                      // kNuisBANFFMAQE
            true,                      // kNuisBANFFMARES
            false,                     // kNuisBANFFCCQEE0
#if BANFF_BINS != 3
            false,                     // kNuisBANFFCCQEE1
            false,                     // kNuisBANFFCCQEE2
#endif
            false,                     // kNuisBANFFCC1piE0
#if BANFF_BINS != 3
            false,                     // kNuisBANFFCC1piE1
#endif
#if (BANFF_BINS != 2)
            false,                     // kNuisBANFFNC1pi0
#endif
            true,                      // Fermi momentum
            true,                      // spectral function on/off
            true,                      // W-shape
            true,                      // CC other shape
            true,                      // Pion-less delta decay
#ifndef __NUEONLY__
            true,                      // Binding energy
#endif
            false,                     // CC coherent norm
#ifndef __NUEONLY__
            false,                     // NC 1piC norm
#endif
#ifndef __NUMUONLY__
            false,                     // NC coh norm
#endif
            false,                     // NC other norm
#ifndef __NUEONLY__
            false,                     // antinu / nu cross section
#endif
            false,                     // nue / numu cross section
            false                      // RHC / FHC normalisation
            
#elif defined(__USE_NIWG2015__)
            // NEW XSEC
            false,                      // RPA tuning kNuisXSecRPA
            
            false,                     // kNuisBANFFMEC
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
            true,                      // kNuisBANFFBindE
            false,                     // kNuisBANFFNormCCcoh
            false,                     // kNuisBANFFNormNCcoh
            false,                     // kNuisBANFFNormNCoth
            false,                     // kNuisBANFFNormNuEToNuMu
            false,                     // kNuisBANFFNormNC1Gamma
            false,                     // kNuisBANFFNormNuEBarToNuMuBar
            false,                     // kNuisBANFFMECBar

#endif
        };
        
        static const bool kRHCXSecSystHasWghtSplines[kNumOfRHCXSecSyst] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            true,                      // kNuisBANFFMAQE
            true,                      // kNuisBANFFMARES
            false,                     // kNuisBANFFCCQEE0
#if BANFF_BINS != 3
            false,                     // kNuisBANFFCCQEE1
            false,                     // kNuisBANFFCCQEE2
#endif
            false,                     // kNuisBANFFCC1piE0
#if BANFF_BINS != 3
            false,                     // kNuisBANFFCC1piE1
#endif
#if (BANFF_BINS != 2)
            false,                     // kNuisBANFFNC1pi0
#endif
            true,                      // Fermi momentum
            true,                      // spectral function on/off
            true,                      // W-shape
            true,                      // CC other shape
            true,                      // Pion-less delta decay
#ifndef __NUEONLY__
            true,                      // Binding energy
#endif
            false,                     // CC coherent norm
#ifndef __NUEONLY__
            false,                     // NC 1piC norm
#endif
#ifndef __NUMUONLY__
            false,                     // NC coh norm
#endif
            false,                     // NC other norm
#ifndef __NUEONLY__
            false,                     // antinu / nu cross section
#endif
            false,                     // nue / numu cross section
            false                      // RHC / FHC normalisation
            
#elif defined(__USE_NIWG2015__)
            // NEW XSEC
            false,                      // RPA tuning kRHCNuisXSecRPA
            
            false,                     // kRHCNuisBANFFMEC
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
            true,                      // kRHCNuisBANFFBindE
            false,                     // kRHCNuisBANFFNormCCcoh
            false,                     // kRHCNuisBANFFNormNCcoh
            false,                     // kRHCNuisBANFFNormNCoth
            false,                     // kRHCNuisBANFFNormNuEToNuMu
            false,                     // kRHCNuisBANFFNormNC1Gamma
            false,                     // kRHCNuisBANFFNormNuEBarToNuMuBar
            false,                     // kRHCNuisBANFFMECBar

#endif
        };
        
        // The following list specifies which cross-section related systematic
        // parameters were included in the BANFF fit.
        // The order of systematic parameters should match the one in kListOfXSecSyst.
        static const bool kXSecSystIncludedInBANFF[kNumOfXSecSyst] =
        {
            // NEW XSEC
#if defined(__USE_NIWG2013__)
            true,                      // kNuisBANFFMAQE
            true,                      // kNuisBANFFMARES
            true,                      // kNuisBANFFCCQEE0
#if BANFF_BINS != 3
            true,                      // kNuisBANFFCCQEE1
            true,                      // kNuisBANFFCCQEE2
#endif
            true,                      // kNuisBANFFCC1piE0
#if BANFF_BINS != 3
            true,                      // kNuisBANFFCC1piE1
#endif
#if (BANFF_BINS != 2)
            true,                      // kNuisBANFFNC1pi0
#endif
            false,                     // Fermi momentum
            false,                     // spectral function on/off
            false,                     // W-shape
            false,                     // CC other shape
            false,                     // Pion-less delta decay
#ifndef __NUEONLY__
            false,                     // Binding energy
#endif
            false,                     // CC coherent norm
#ifndef __NUEONLY__
            false,                     // NC 1piC norm
#endif
#ifndef __NUMUONLY__
            false,                     // NC coh norm
#endif
            false,                     // NC other norm
#ifndef __NUEONLY__
            false,                     // antinu / nu cross section
#endif
            false,                     // nue / numu cross section
            false                      // RHC / FHC normalisation
            
#elif defined(__USE_NIWG2015__)
            // NEW XSEC
            false,                     // RPA tuning kNuisXSecRPA
            
            true,                      // kNuisBANFFMEC
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
            true,                      // kNuisBANFFBindE
            true,                      // kNuisBANFFNormCCcoh
            true,                      // kNuisBANFFNormNCcoh
            true,                      // kNuisBANFFNormNCoth
            true,                      // kNuisBANFFNormNuEToNuMu
            true,                     // kNuisBANFFNormNC1Gamma
            true,                     // kNuisBANFFNormNuEBarToNuMuBar
            true,                      // kNuisBANFFMECBar

#endif
        };
        
        
        static const bool kRHCXSecSystIncludedInBANFF[kNumOfRHCXSecSyst] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            true,                      // kNuisBANFFMAQE
            true,                      // kNuisBANFFMARES
            true,                      // kNuisBANFFCCQEE0
#if BANFF_BINS != 3
            true,                      // kNuisBANFFCCQEE1
            true,                      // kNuisBANFFCCQEE2
#endif
            true,                      // kNuisBANFFCC1piE0
#if BANFF_BINS != 3
            true,                      // kNuisBANFFCC1piE1
#endif
#if (BANFF_BINS != 2)
            true,                      // kNuisBANFFNC1pi0
#endif
            false,                     // Fermi momentum
            false,                     // spectral function on/off
            false,                     // W-shape
            false,                     // CC other shape
            false,                     // Pion-less delta decay
#ifndef __NUEONLY__
            false,                     // Binding energy
#endif
            false,                     // CC coherent norm
#ifndef __NUEONLY__
            false,                     // NC 1piC norm
#endif
#ifndef __NUMUONLY__
            false,                     // NC coh norm
#endif
            false,                     // NC other norm
#ifndef __NUEONLY__
            false,                     // antinu / nu cross section
#endif
            false,                     // nue / numu cross section
            false                      // RHC / FHC normalisation
            
#elif defined(__USE_NIWG2015__)
            // NEW XSEC
            false,                     // RPA tuning kRHCNuisXSecRPA
            
            true,                      // kRHCNuisBANFFMEC
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
            true,                      // kRHCNuisBANFFBindE
            true,                      // kRHCNuisBANFFNormCCcoh
            true,                      // kRHCNuisBANFFNormNCcoh
            true,                      // kRHCNuisBANFFNormNCoth
            true,                      // kRHCNuisBANFFNormNuEToNuMu
            true,                     // kRHCNuisBANFFNormNC1Gamma
            true,                     // kRHCNuisBANFFNormNuEBarToNuMuBar
            true,                      // kRHCNuisBANFFMECBar

#endif
        };
        
        // Ereco and Etrue cuts for the application of cross-section systematics.
        // The order of systematic parameters should match the one in kListOfXSecSyst.
        static const double kXSecSystErecoCut[kNumOfXSecSyst][2] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            { 0.0, 100.0 },            // kNuisBANFFMAQE
            { 0.0, 100.0 },            // kNuisBANFFMARES
            { 0.0, 100.0 },            // kNuisBANFFCCQEE0
#if BANFF_BINS != 3
            { 0.0, 100.0 },            // kNuisBANFFCCQEE1
            { 0.0, 100.0 },            // kNuisBANFFCCQEE2
#endif
            { 0.0, 100.0 },            // kNuisBANFFCC1piE0
#if BANFF_BINS != 3
            { 0.0, 100.0 },            // kNuisBANFFCC1piE1
#endif
#if (BANFF_BINS != 2)
            { 0.0, 100.0 },            // kNuisBANFFNC1pi0
#endif
            { 0.0, 100.0 },            // Fermi momentum
            { 0.0, 100.0 },            // spectral function on/off
            { 0.0, 100.0 },            // W-shape
            { 0.0, 100.0 },            // CC other shape
            { 0.0, 100.0 },            // Pion-less delta decay
#ifndef __NUEONLY__
            { 0.0, 100.0 },            // Binding energy
#endif
            { 0.0, 100.0 },            // CC coherent norm
#ifndef __NUEONLY__
            { 0.0, 100.0 },            // NC 1piC norm
#endif
#ifndef __NUMUONLY__
            { 0.0, 100.0 },            // NC coh norm
#endif
            { 0.0, 100.0 },            // NC other norm
#ifndef __NUEONLY__
            { 0.0, 100.0 },            // antinu / nu cross section
#endif
            { 0.0, 100.0 },            // nue / numu cross section
            { 0.0, 100.0 }             // RHC / FHC normalisation
            
#elif defined(__USE_NIWG2015__)
            // NEW XSEC
            { 0.0, 100.0 },            // RPA tuning kNuisXSec RPA
            
            { 0.0, 100.0 },            // kNuisBANFFMEC
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
            { 0.0, 100.0 },            // kNuisBANFFBindE
            { 0.0, 100.0 },            // kNuisBANFFNormCCcoh
            { 0.0, 100.0 },            // kNuisBANFFNormNCcoh
            { 0.0, 100.0 },            // kNuisBANFFNormNCoth
            { 0.0, 100.0 },            // kNuisBANFFNormNuEToNuMu
            { 0.0, 100.0 },            // kNuisBANFFNormNC1Gamma
            { 0.0, 100.0 },            // kNuisBANFFNormNuEBarToNuMuBar
            { 0.0, 100.0 },            // kNuisBANFFMECBar

#endif
        };
        
        static const double kRHCXSecSystErecoCut[kNumOfRHCXSecSyst][2] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            { 0.0, 100.0 },            // kNuisBANFFMAQE
            { 0.0, 100.0 },            // kNuisBANFFMARES
            { 0.0, 100.0 },            // kNuisBANFFCCQEE0
#if BANFF_BINS != 3
            { 0.0, 100.0 },            // kNuisBANFFCCQEE1
            { 0.0, 100.0 },            // kNuisBANFFCCQEE2
#endif
            { 0.0, 100.0 },            // kNuisBANFFCC1piE0
#if BANFF_BINS != 3
            { 0.0, 100.0 },            // kNuisBANFFCC1piE1
#endif
#if (BANFF_BINS != 2)
            { 0.0, 100.0 },            // kNuisBANFFNC1pi0
#endif
            { 0.0, 100.0 },            // Fermi momentum
            { 0.0, 100.0 },            // spectral function on/off
            { 0.0, 100.0 },            // W-shape
            { 0.0, 100.0 },            // CC other shape
            { 0.0, 100.0 },            // Pion-less delta decay
#ifndef __NUEONLY__
            { 0.0, 100.0 },            // Binding energy
#endif
            { 0.0, 100.0 },            // CC coherent norm
#ifndef __NUEONLY__
            { 0.0, 100.0 },            // NC 1piC norm
#endif
#ifndef __NUMUONLY__
            { 0.0, 100.0 },            // NC coh norm
#endif
            { 0.0, 100.0 },            // NC other norm
#ifndef __NUEONLY__
            { 0.0, 100.0 },            // antinu / nu cross section
#endif
            { 0.0, 100.0 },            // nue / numu cross section
            { 0.0, 100.0 }             // RHC / FHC normalisation
            
#elif defined(__USE_NIWG2015__)
            // NEW XSEC
            { 0.0, 100.0 },            // RPA tuning kRHCNuisXSecRPA
            
            { 0.0, 100.0 },            // kRHCNuisBANFFMEC
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
            { 0.0, 100.0 },            // kRHCNuisBANFFBindE
            { 0.0, 100.0 },            // kRHCNuisBANFFNormCCcoh
            { 0.0, 100.0 },            // kRHCNuisBANFFNormNCcoh
            { 0.0, 100.0 },            // kRHCNuisBANFFNormNCoth
            { 0.0, 100.0 },            // kRHCNuisBANFFNormNuEToNuMu
            { 0.0, 100.0 },            // kRHCNuisBANFFNormNC1Gamma
            { 0.0, 100.0 },            // kRHCNuisBANFFNormNuEBarToNuMuBar
            { 0.0, 100.0 },            // kRHCNuisBANFFMECBar

#endif
        };
        
        static const double kXSecSystEtrueCut[kNumOfXSecSyst][2] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            { 0.0, 100.0 },            // kNuisBANFFMAQE
            { 0.0, 100.0 },            // kNuisBANFFMARES
            { 0.0,   1.5 },            // kNuisBANFFCCQEE0
#if BANFF_BINS != 3
            { 1.5,   3.5 },            // kNuisBANFFCCQEE1
            { 3.5,  30.0 },            // kNuisBANFFCCQEE2
#endif
            { 0.0,   2.5 },            // kNuisBANFFCC1piE0
#if BANFF_BINS != 3
            { 2.5,  30.0 },            // kNuisBANFFCC1piE1
#endif
#if (BANFF_BINS != 2)
            { 0.0, 100.0 },            // kNuisBANFFNC1pi0
#endif
            { 0.0, 100.0 },            // Fermi momentum
            { 0.0, 100.0 },            // spectral function on/off
            { 0.0, 100.0 },            // W-shape
            { 0.0, 100.0 },            // CC other shape
            { 0.0, 100.0 },            // Pion-less delta decay
#ifndef __NUEONLY__
            { 0.0, 100.0 },            // Binding energy
#endif
            { 0.0, 100.0 },            // CC coherent norm
#ifndef __NUEONLY__
            { 0.0, 100.0 },            // NC 1piC norm
#endif
#ifndef __NUMUONLY__
            { 0.0, 100.0 },            // NC coh norm
#endif
            { 0.0, 100.0 },            // NC other norm
#ifndef __NUEONLY__
            { 0.0, 100.0 },            // antinu / nu cross section
#endif
            { 0.0, 100.0 },            // nue / numu cross section
            { 0.0, 100.0 }             // RHC / FHC normalisation
            
#elif defined(__USE_NIWG2015__)
            // NEW XSEC
            { 0.0, 100.0 },            // RPA tuning kNuisXSecRPA
            
            { 0.0, 100.0 },            // kNuisBANFFMEC
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
            { 0.0, 100.0 },            // kNuisBANFFBindE
            { 0.0, 100.0 },            // kNuisBANFFNormCCcoh
            { 0.0, 100.0 },            // kNuisBANFFNormNCcoh
            { 0.0, 100.0 },            // kNuisBANFFNormNCoth
            { 0.0, 100.0 },            // kNuisBANFFNormNuEToNuMu
            { 0.0, 100.0 },            // kNuisBANFFNormNC1Gamma
            { 0.0, 100.0 },            // kNuisBANFFNormNuEBarToNuMuBar
            { 0.0, 100.0 },            // kNuisBANFFMECBar

#endif
        };
        
        static const double kRHCXSecSystEtrueCut[kNumOfRHCXSecSyst][2] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            { 0.0, 100.0 },            // kNuisBANFFMAQE
            { 0.0, 100.0 },            // kNuisBANFFMARES
            { 0.0,   1.5 },            // kNuisBANFFCCQEE0
#if BANFF_BINS != 3
            { 1.5,   3.5 },            // kNuisBANFFCCQEE1
            { 3.5,  30.0 },            // kNuisBANFFCCQEE2
#endif
            { 0.0,   2.5 },            // kNuisBANFFCC1piE0
#if BANFF_BINS != 3
            { 2.5,  30.0 },            // kNuisBANFFCC1piE1
#endif
#if (BANFF_BINS != 2)
            { 0.0, 100.0 },            // kNuisBANFFNC1pi0
#endif
            { 0.0, 100.0 },            // Fermi momentum
            { 0.0, 100.0 },            // spectral function on/off
            { 0.0, 100.0 },            // W-shape
            { 0.0, 100.0 },            // CC other shape
            { 0.0, 100.0 },            // Pion-less delta decay
#ifndef __NUEONLY__
            { 0.0, 100.0 },            // Binding energy
#endif
            { 0.0, 100.0 },            // CC coherent norm
#ifndef __NUEONLY__
            { 0.0, 100.0 },            // NC 1piC norm
#endif
#ifndef __NUMUONLY__
            { 0.0, 100.0 },            // NC coh norm
#endif
            { 0.0, 100.0 },            // NC other norm
#ifndef __NUEONLY__
            { 0.0, 100.0 },            // antinu / nu cross section
#endif
            { 0.0, 100.0 },            // nue / numu cross section
            { 0.0, 100.0 }             // RHC / FHC normalisation
            
#elif defined(__USE_NIWG2015__)
            // NEW XSEC
            { 0.0, 100.0 },            // RPA tuning kRHCNuisXSecRPA
            
            { 0.0, 100.0 },            // kRHCNuisBANFFMEC
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
            { 0.0, 100.0 },            // kRHCNuisBANFFBindE
            { 0.0, 100.0 },            // kRHCNuisBANFFNormCCcoh
            { 0.0, 100.0 },            // kRHCNuisBANFFNormNCcoh
            { 0.0, 100.0 },            // kRHCNuisBANFFNormNCoth
            { 0.0, 100.0 },            // kRHCNuisBANFFNormNuEToNuMu
            { 0.0, 100.0 },            // kRHCNuisBANFFNormNC1Gamma
            { 0.0, 100.0 },            // kRHCNuisBANFFNormNuEBarToNuMuBar
            { 0.0, 100.0 },            // kRHCNuisBANFFMECBar

#endif
        };
        
        
        //Fractional error of systematic parameters, as of BANFF 2012v8
        // ONLY USE FOR: macros/misc/test_fitter_inputs/test_valor_kdi_splines.C
        static const double kXSecSystFracErr[kNumOfXSecSyst] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            0.0716,          // BANFF, M_{A}QE
            0.0675,          // BANFF, M_{A}RES
            0.09,           // BANFF, CCQE   normalization, 0 - 1.5 GeV
#if BANFF_BINS != 3
            0.3,            // BANFF, CCQE   normalization, 1.5 - 3.5 GeV
            0.3,            // BANFF, CCQE   normalization, 3.5 - 30.0 GeV
#endif
            0.149,          // BANFF, CC1pi  normalization, 0 - 2.5 GeV
#if BANFF_BINS != 3
            0.4,            // BANFF, CC1pi  normalization, 2.5 - 30.0 GeV
#endif
#if (BANFF_BINS != 2)
            0.3,            // BANFF, NC1pi0 normalization, 0.0 - 30.0 GeV
#endif
            0.13333,        // Fermi momentum
            1.0,            // spectral function on/off
            0.516,          // W-shape
            0.4,            // CC other shape
            0.2,            // Pion-less delta decay
#ifndef __NUEONLY__
            0.33333,        // Binding enery
#endif
            1.0,            // CC coherent norm
#ifndef __NUEONLY__
            0.3,            // NC 1piC norm
#endif
#ifndef __NUMUONLY__
            0.3,            // NC coh norm
#endif
            0.3,            // NC other norm
#ifndef __NUEONLY__
            0.4,            // antinu / nu cross section
#endif
            0.03,           // nue / numu cross section
            0.10            // RHC / FHC normalisation
            
#elif defined(__USE_NIWG2015__)
            // NEW XSEC
            0.1,            // RPA tuning
            
            0.1,            // BANFF, MEC
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
            0.5,			// BANFF, NC 1 gamma
            0.03            // BANFF, nuebar / numubar cross section
#endif
        };
        
        static const double kRHCXSecSystFracErr[kNumOfRHCXSecSyst] =
        {
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            0.0716,         // BANFF, M_{A}QE
            0.0675,         // BANFF, M_{A}RES
            0.09,           // BANFF, CCQE   normalization, 0 - 1.5 GeV
#if BANFF_BINS != 3
            0.3,            // BANFF, CCQE   normalization, 1.5 - 3.5 GeV
            0.3,            // BANFF, CCQE   normalization, 3.5 - 30.0 GeV
#endif
            0.149,          // BANFF, CC1pi  normalization, 0 - 2.5 GeV
#if BANFF_BINS != 3
            0.4,            // BANFF, CC1pi  normalization, 2.5 - 30.0 GeV
#endif
#if (BANFF_BINS != 2)
            0.3,            // BANFF, NC1pi0 normalization, 0.0 - 30.0 GeV
#endif
            0.13333,        // Fermi momentum
            1.0,            // spectral function on/off
            0.516,          // W-shape
            0.4,            // CC other shape
            0.2,            // Pion-less delta decay
#ifndef __NUEONLY__
            0.33333,        // Binding enery
#endif
            1.0,            // CC coherent norm
#ifndef __NUEONLY__
            0.3,            // NC 1piC norm
#endif
#ifndef __NUMUONLY__
            0.3,            // NC coh norm
#endif
            0.3,            // NC other norm
#ifndef __NUEONLY__
            0.4,            // antinu / nu cross section
#endif
            0.03,           // nue / numu cross section
            0.10            // RHC / FHC normalisation
            
#elif defined(__USE_NIWG2015__)
            // NEW XSEC
            0.1,            // RPA tuning
            
            0.1,            // BANFF, MEC
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
            0.5,			// BANFF, NC 1 gamma
            0.03            // BANFF, nuebar / numubar cross section

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
#if defined(__USE_NIWG2013__)
                // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
                // banffcc1piE0 banffcc1piE1 banffnc1pi0
                // pF  spectral   W-shape  ccothshape pilessdelta bindE
                // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
                // RHC/FHC
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuBarNC1coh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuBarMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1  }, // kSKmNuECCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1  }, // kSKmNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1  }, // kSKmNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuECCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1  }, // kSKmNuECCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0  }, // kSKmNuENC1piC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENC1pi0
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENCoth
                
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuECCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECCoth
                
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuEBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNCoth
#ifdef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
                
#elif defined(__USE_NIWG2015__)
                /// NEW XSEC
                /// rpa
                /// banffmec banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
                /// banffpF  banffccothshape banffbindE
                /// banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, /// kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, /// kSKmOscNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuENC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCoth
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, /// kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, /// kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
#endif
#endif
            }//1Rmu
            ,
            // > 1 electron-like ring
            {
                // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
                // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
                // banffcc1piE0 banffcc1piE1 banffnc1pi0
                // pF  spectral   W-shape  ccothshape pilessdelta bindE
                // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
                // RHC/FHC
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOScNuMuCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1  }, // kSKmNuECCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1  }, // kSKmNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1  }, // kSKmNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuECCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1  }, // kSKmNuECCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENC1piC
                {  0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuENC1pi0
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0  }, // kSKmNuENCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENCoth
                
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuECCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECCoth
                
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuEBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNCoth
#ifdef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
                
#elif defined(__USE_NIWG2015__)
                /// NEW XSEC
                /// rpa
                /// banffmec banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
                /// banffpF  banffccothshape banffbindE
                /// banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOScNuMuCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, /// kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, /// kSKmOscNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuENC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuENCoth_e
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCoth
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, /// kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, /// kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
#endif
#endif
            }//1Re
            ,
            // > 1 muon-like ring RHC
            {
                // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
                // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
                // banffcc1piE0 banffcc1piE1 banffnc1pi0
                // pF  spectral   W-shape  ccothshape pilessdelta bindE
                // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
                // RHC/FHC
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuBarMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1  }, // kSKmNuECCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1  }, // kSKmNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1  }, // kSKmNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuECCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1  }, // kSKmNuECCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0  }, // kSKmNuENC1piC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENC1pi0
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENCoth
                
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuECCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECCoth
                
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuEBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNCoth
#ifdef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
                
#elif defined(__USE_NIWG2015__)
                /// NEW XSEC
                /// rpa
                /// banffmec banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
                /// banffpF  banffccothshape banffbindE
                /// banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, /// kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, /// kSKmOscNuMuBarMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuENC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCoth
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, /// kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, /// kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
#endif
#endif
            }//1Rmu RHC
            ,
            // > 1 electron-like ring RHC
            {
                // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
                // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
                // banffcc1piE0 banffcc1piE1 banffnc1pi0
                // pF  spectral   W-shape  ccothshape pilessdelta bindE
                // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
                // RHC/FHC
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOScNuMuCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1  }, // kSKmNuECCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1  }, // kSKmNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1  }, // kSKmNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuECCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1  }, // kSKmNuECCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENC1piC
                {  0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuENC1pi0
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0  }, // kSKmNuENCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENCoth
                
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuECCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECCoth
                
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuEBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNCoth
#ifdef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
                
#elif defined(__USE_NIWG2015__)
                /// NEW XSEC
                /// rpa
                /// banffmec banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
                /// banffpF  banffccothshape banffbindE
                /// banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOScNuMuCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, /// kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, /// kSKmOscNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuENC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCoth
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, /// kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, /// kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
#endif
#endif
            }// 1Re RHC
            ,
            // > muon-like multiring
            {
                // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
                // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
                // banffcc1piE0 banffcc1piE1 banffnc1pi0
                // pF  spectral   W-shape  ccothshape pilessdelta bindE
                // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
                // RHC/FHC
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuBarMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1  }, // kSKmNuECCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1  }, // kSKmNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1  }, // kSKmNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuECCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1  }, // kSKmNuECCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0  }, // kSKmNuENC1piC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENC1pi0
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENCoth
                
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuECCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECCoth
                
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuEBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNCoth
#ifdef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
                
#elif defined(__USE_NIWG2015__)
                /// NEW XSEC
                /// rpa
                /// banffmec banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
                /// banffpF  banffccothshape banffbindE
                /// banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, /// kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, /// kSKmOscNuMuBarMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuENC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCoth
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, /// kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, /// kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
#endif
#endif
            }//1Rmu
            ,
            // > electron-like multiring
            {
                // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
                // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
                // banffcc1piE0 banffcc1piE1 banffnc1pi0
                // pF  spectral   W-shape  ccothshape pilessdelta bindE
                // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
                // RHC/FHC
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOScNuMuCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0  }, // kSKmNuMuBarCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1  }, // kSKmNuECCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1  }, // kSKmNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1  }, // kSKmNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuECCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1  }, // kSKmNuECCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENC1piC
                {  0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuENC1pi0
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0  }, // kSKmNuENCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuENCoth
                
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuECCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1  }, // kSKmOscNuECCoth
                
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuEBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1  }, // kSKmNuEBarCCoth
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0  }, // kSKmNuEBarNCoth
#ifdef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCMEC
                {  0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
                
#elif defined(__USE_NIWG2015__)
                /// NEW XSEC
                /// rpa
                /// banffmec banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
                /// banffpF  banffccothshape banffbindE
                /// banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOScNuMuCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, /// kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, /// kSKmOscNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmOscNuMuBarCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmNuECCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuENC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 }, /// kSKmOscNuECCoth
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, /// kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmNuEBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, /// kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, /// kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, /// kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, /// kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0 }, /// kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, /// kSKmUndefined
#endif
#endif
            }// MR e
        };
        
        
#elif BANFF_BINS == 2
        
        //Todo: make this work if kNumOfRHCXSecSyst is different from kNumOfXSecSyst.
        static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {
            
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            
            // > 1 muon-like ring
            {
                // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
                // banffcc1piE0 banffcc1piE1
                // pF  spectral   W-shape  ccothshape pilessdelta bindE
                // cccoh    nc1piC     ncoth    antinu/nu   nue/numu
                // RHC/FHC
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0  }, // kSKmNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0  }, // kSKmNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0  }, // kSKmNuMuBarCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0  }, // kSKmNuMuBarCCoth
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1  }, // kSKmNuECCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1  }, // kSKmNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1  }, // kSKmNuECCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1  }, // kSKmNuECCoth
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0  }, // kSKmNuENC1piC
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuENC1pi0
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuENCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuENCoth
                
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1  }, // kSKmOscNuECCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1  }, // kSKmOscNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1  }, // kSKmOscNuECCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1  }, // kSKmOscNuECCoth
                
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1  }, // kSKmNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1  }, // kSKmNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1  }, // kSKmNuEBarCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1  }, // kSKmNuEBarCCoth
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuEBarNCoth
                
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1  }, // kSKmOscNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1  }, // kSKmOscNuEBarCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCoth
                
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
            }//1Rmu
            ,
            // > 1 electron-like ring
            {
                //if BANFF_BINS == 2 - NUMUONLY!
                // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
                // banffcc1piE0 banffcc1piE1 banffnc1pi0
                // pF  spectral   W-shape  ccothshape pilessdelta bindE
                // cccoh    nc1piC     ncoth    antinu/nu  nue/numu
                // RHC/FHC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCoth
#endif
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarCC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarCCoth
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuECCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuECC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuECCoth
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuENC1piC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuENC1pi0
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuENCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuENCoth
                
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuECCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuECC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuECCoth
                
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarCC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarCCoth
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarNCoth
                
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuEBarCC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
            }//1Re
            ,
            // > 1 muon-like ring RHC
            {
                // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
                // banffcc1piE0 banffcc1piE1
                // pF  spectral   W-shape  ccothshape pilessdelta bindE
                // cccoh    nc1piC     ncoth    antinu/nu   nue/numu
                // RHC/FHC
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0  }, // kSKmNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0  }, // kSKmNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0  }, // kSKmNuMuBarCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0  }, // kSKmNuMuBarCCoth
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1  }, // kSKmNuECCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1  }, // kSKmNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1  }, // kSKmNuECCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1  }, // kSKmNuECCoth
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0  }, // kSKmNuENC1piC
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuENC1pi0
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuENCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuENCoth
                
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1  }, // kSKmOscNuECCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1  }, // kSKmOscNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1  }, // kSKmOscNuECCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1  }, // kSKmOscNuECCoth
                
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1  }, // kSKmNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1  }, // kSKmNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1  }, // kSKmNuEBarCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1  }, // kSKmNuEBarCCoth
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuEBarNCoth
                
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1  }, // kSKmOscNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1  }, // kSKmOscNuEBarCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCoth
                
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
            }//1RmuRHC
            ,
            // > 1 electron-like ring RHC
            {
                //if BANFF_BINS == 2 - NUMUONLY!
                // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
                // banffcc1piE0 banffcc1piE1 banffnc1pi0
                // pF  spectral   W-shape  ccothshape pilessdelta bindE
                // cccoh    nc1piC     ncoth    antinu/nu  nue/numu
                // RHC/FHC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCoth
#endif
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarCC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarCCoth
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuECCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuECC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuECCoth
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuENC1piC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuENC1pi0
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuENCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuENCoth
                
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuECCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuECC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuECCoth
                
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarCC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarCCoth
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarNCoth
                
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuEBarCC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
            }//1ReRHC
            
            /// NEW MULTIRING
            ,
            // > Multi muon-like ring
            {
                // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
                // banffcc1piE0 banffcc1piE1
                // pF  spectral   W-shape  ccothshape pilessdelta bindE
                // cccoh    nc1piC     ncoth    antinu/nu   nue/numu
                // RHC/FHC
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0  }, // kSKmNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0  }, // kSKmNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0  }, // kSKmNuMuBarCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0  }, // kSKmNuMuBarCCoth
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1  }, // kSKmNuECCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1  }, // kSKmNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1  }, // kSKmNuECCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1  }, // kSKmNuECCoth
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0  }, // kSKmNuENC1piC
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuENC1pi0
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuENCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuENCoth
                
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1  }, // kSKmOscNuECCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1  }, // kSKmOscNuECC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1  }, // kSKmOscNuECCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1  }, // kSKmOscNuECCoth
                
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1  }, // kSKmNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1  }, // kSKmNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1  }, // kSKmNuEBarCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1  }, // kSKmNuEBarCCoth
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0  }, // kSKmNuEBarNCoth
                
                {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCQE
                {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1  }, // kSKmOscNuEBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1  }, // kSKmOscNuEBarCCcoh
                {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1  }, // kSKmOscNuEBarCCoth
                
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
            }//MultiRmu
            ,
            // > Multi electron-like ring
            {
                //if BANFF_BINS == 2 - NUMUONLY!
                // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
                // banffcc1piE0 banffcc1piE1 banffnc1pi0
                // pF  spectral   W-shape  ccothshape pilessdelta bindE
                // cccoh    nc1piC     ncoth    antinu/nu  nue/numu
                // RHC/FHC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuCCoth
#endif
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarCC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarCCoth
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuECCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuECC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuECCoth
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuENC1piC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuENC1pi0
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuENCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuENCoth
                
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuECCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuECC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuECCoth
                
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarCC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarCCoth
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmNuEBarNCoth
                
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCQE
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuEBarCC0pi
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCcoh
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0  }, // kSKmUndefined
#endif
            }//MultiRe
            
            
#elif defined(__USE_NIWG2015__)
            
            // > 1 muon-like ring
            {
                // NEW XSEC
                // NB: banffscca and banffsccv are not applied to muon-like ring
                // rpa
                // banffmec banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
                // banffpF  banffccothshape banffbindE
                // banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 }, // kSKmNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 }, // kSKmOscNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCoth
                
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
            }//1Rmu
            ,
            // > 1 electron-like ring
            {
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1Gamma
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1Gamma
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
            }//1Re
            
            ,
            // > 1 muon-like ring RHC
            {
                // NEW XSEC
                // NB: banffscca and banffsccv are not applied to muon-like ring
                // rpa
                // banffmec banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
                // banffpF  banffccothshape banffbindE
                // banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarMEC
                {  0,   0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 }, // kSKmNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 }, // kSKmOscNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCoth
                
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
            }//1RmuRHC
            ,
            // > 1 electron-like ring RHC
            {
                // NEW XSEC
                // NB: banffscca and banffsccv are not applied to muon-like ring
                // rpa
                // banffmec banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
                // banffpF  banffccothshape banffbindE
                // banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOScNuMuCCoth
#endif
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCoth
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
            }//1ReRHC
            
            /// NEW MULTIRING
            ,
            // > Multi muon-like ring
            {
                // NEW XSEC
                // NB: banffscca and banffsccv are not applied to muon-like ring
                // rpa
                // banffmec banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
                // banffpF  banffccothshape banffbindE
                // banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarMEC
                {  0,   0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 }, // kSKmNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 }, // kSKmOscNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCoth
                
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCQE
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
            }//MultiRmu
            ,
            // > Multi electron-like ring
            {
                // NEW XSEC
                // NB: banffscca and banffsccv are not applied to muon-like ring
                // rpa
                // banffmec banffca5 banffbgres banffmaqe banffmares banffscca banffsccv
                // banffpF  banffccothshape banffbindE
                // banffcccoh  banffnccoh  banffncoth  banffnue/numu banffNC1gamma
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOScNuMuCCoth
#endif
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC0pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuECCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmOscNuECCoth
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuMuNC1Gamma
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCQE
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
            }//MultiRe
            
#endif
        };
        
#elif BANFF_BINS == 3 // NUE APP ONLY
        
        // Todo: add new systematics and remove old ones
        // Todo: make this work if kNumOfRHCXSecSyst is different from kNumOfXSecSyst.
        
        static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {
            
            // NEW XSEC NIWG
#if defined(__USE_NIWG2013__)
            // > 1 muon-like ring
            {
                // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc1pi0
                // pF  spectral   W-shape  ccothshape pilessdelta
                // cccoh    nccoh   ncoth    nue/numu
                // RHC/FHC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuCCoth
#endif
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCoth
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuECCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuECC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuECCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuECCoth
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuENC1piC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuENC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuENCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuENCoth
                
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCoth
                
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarCCoth
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarNCoth
                
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCoth
                
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
#endif
                
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
            }//1Rmu
            ,
            // > 1 electron-like ring
            {
                // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc0pi0
                // pF  spectral   W-shape  ccothshape pilessdelta
                // cccoh   nccoh   ncoth    nue/numu
                // RHC/FHC
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuMuNC1piC
                {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0  }, // kSKmNuMuCCoth
#endif
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuBarCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0  }, // kSKmNuMuBarCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCoth
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuMuBarNC1piC
                {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1  }, // kSKmNuECCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1  }, // kSKmNuECC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1  }, // kSKmNuECCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1  }, // kSKmNuECCoth
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuENC1piC
                {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0  }, // kSKmNuENC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0  }, // kSKmNuENCcoh
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuENCoth
                
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1  }, // kSKmOscNuECCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1  }, // kSKmOscNuECC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1  }, // kSKmOscNuECCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1  }, // kSKmOscNuECCoth
                
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1  }, // kSKmNuEBarCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1  }, // kSKmNuEBarCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1  }, // kSKmNuEBarCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1  }, // kSKmNuEBarCCoth
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuEBarNC1piC
                {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuEBarNCoth
                
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1  }, // kSKmOscNuEBarCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1  }, // kSKmOscNuEBarCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1  }, // kSKmOscNuEBarCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1  }, // kSKmOscNuEBarCCoth
                
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
#endif
            }//1Re
            ,
            // > 1 muon-like ring RHC
            {
                // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc1pi0
                // pF  spectral   W-shape  ccothshape pilessdelta
                // cccoh    nccoh   ncoth    nue/numu
                // RHC/FHC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuCCoth
#endif
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCoth
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuECCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuECC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuECCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuECCoth
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuENC1piC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuENC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuENCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuENCoth
                
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCoth
                
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarCCoth
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarNCoth
                
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCoth
                
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
#endif
            }//1RmuRHC
            ,
            // > 1 electron-like ring RHC
            {
                // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc0pi0
                // pF  spectral   W-shape  ccothshape pilessdelta
                // cccoh   nccoh   ncoth    nue/numu
                // RHC/FHC
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuMuNC1piC
                {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0  }, // kSKmNuMuCCoth
#endif
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuBarCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0  }, // kSKmNuMuBarCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCoth
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuMuBarNC1piC
                {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1  }, // kSKmNuECCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1  }, // kSKmNuECC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1  }, // kSKmNuECCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1  }, // kSKmNuECCoth
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuENC1piC
                {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0  }, // kSKmNuENC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0  }, // kSKmNuENCcoh
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuENCoth
                
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1  }, // kSKmOscNuECCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1  }, // kSKmOscNuECC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1  }, // kSKmOscNuECCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1  }, // kSKmOscNuECCoth
                
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1  }, // kSKmNuEBarCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1  }, // kSKmNuEBarCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1  }, // kSKmNuEBarCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1  }, // kSKmNuEBarCCoth
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuEBarNC1piC
                {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuEBarNCoth
                
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1  }, // kSKmOscNuEBarCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1  }, // kSKmOscNuEBarCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1  }, // kSKmOscNuEBarCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1  }, // kSKmOscNuEBarCCoth
                
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
#endif
            }//1ReRHC
            /// NEW MULTIRING
            ,
            // > Multi muon-like ring
            {
                // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc1pi0
                // pF  spectral   W-shape  ccothshape pilessdelta
                // cccoh    nccoh   ncoth    nue/numu
                // RHC/FHC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuNC1piC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuCCoth
#endif
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCoth
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarNC1piC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarNCoth
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuECCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuECC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuECCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuECCoth
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuENC1piC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuENC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuENCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuENCoth
                
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCoth
                
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarCCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarCC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarCCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarCCoth
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarNC1piC
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarNCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuEBarNCoth
                
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCQE
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECC0pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCcoh
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuECCoth
                
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
#endif
            }//MultiRmu
            ,
            // > Multi electron-like ring
            {
                // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc0pi0
                // pF  spectral   W-shape  ccothshape pilessdelta
                // cccoh   nccoh   ncoth    nue/numu
                // RHC/FHC
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0  }, // kSKmNuMuCCoth
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuMuNC1piC
                {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0  }, // kSKmNuMuNCcoh
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuMuNCoth
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0  }, // kSKmNuMuCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0  }, // kSKmNuMuCCoth
#endif
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuBarCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0  }, // kSKmNuMuBarCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0  }, // kSKmNuMuBarCCoth
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuMuBarNC1piC
                {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0  }, // kSKmNuMuBarNC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0  }, // kSKmNuMuBarNCcoh
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0  }, // kSKmOscNuMuBarCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0  }, // kSKmOscNuMuBarCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0  }, // kSKmOscNuMuBarCCoth
#endif
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1  }, // kSKmNuECCQE
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1  }, // kSKmNuECC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1  }, // kSKmNuECCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1  }, // kSKmNuECCoth
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuENC1piC
                {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0  }, // kSKmNuENC1pi0
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0  }, // kSKmNuENCcoh
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuENCoth
                
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1  }, // kSKmOscNuECCQE        
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1  }, // kSKmOscNuECC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1  }, // kSKmOscNuECCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1  }, // kSKmOscNuECCoth
                
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1  }, // kSKmNuEBarCCQE        
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1  }, // kSKmNuEBarCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1  }, // kSKmNuEBarCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1  }, // kSKmNuEBarCCoth
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuEBarNC1piC
                {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0  }, // kSKmNuEBarNC1pi0 
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0  }, // kSKmNuEBarNCcoh 
                {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0  }, // kSKmNuEBarNCoth
                
                {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1  }, // kSKmOscNuEBarCCQE        
                {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1  }, // kSKmOscNuEBarCC1pi
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1  }, // kSKmOscNuEBarCCcoh
                {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1  }, // kSKmOscNuEBarCCoth
                
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
                {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0  }, // kSKmUndefined
#endif
            }//MultiRe
            
            // NEW XSEC
#elif defined(__USE_NIWG2015__)
            // > 1 muon-like ring 
            {
                // rpa 
                // banffmec banffca5 banffbgres
                // banffmaqe banffmares banffscca banffsccv 
                // banffpf banffccothshape banffebo banffcccoh banffnccoh banffncoth banffnue/numu banffNC1gamma
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuENC1Gamma
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1Gamma
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
            }//1Rmu
            ,
            // > 1 electron-like ring
            {
                // rpa 
                // banffmec banffca5 banffbgres
                // banffmaqe banffmares banffscca banffsccv 
                // banffpf banffccothshape banffebo banffcccoh banffnccoh banffncoth banffnue/numu banffNC1gamma
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuNCcoh 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_e
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 }, // kSKmNuECCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 }, // kSKmNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, // kSKmNuECCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuENCcoh 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuENC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 }, // kSKmOscNuECCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 }, // kSKmOscNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCoth
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuEBarNCcoh 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 }, // kSKmOscNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif 
            }//1Re
            ,
            // > 1 muon-like ring RHC
            {
                // rpa 
                // banffmec banffca5 banffbgres
                // banffmaqe banffmares banffscca banffsccv 
                // banffpf banffccothshape banffebo banffcccoh banffnccoh banffncoth banffnue/numu banffNC1gamma
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC0pi   
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuENC1Gamma
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1Gamma
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
            }//1Rmu RHC
            ,
            // > 1 electron-like ring RHC
            {
                // rpa 
                // banffmec banffca5 banffbgres
                // banffmaqe banffmares banffscca banffsccv 
                // banffpf banffccothshape banffebo banffcccoh banffnccoh banffncoth banffnue/numu banffNC1gamma
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuNCcoh 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_e
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_e
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif 
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 }, // kSKmNuECCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 }, // kSKmNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, // kSKmNuECCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuENCcoh 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuENC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 }, // kSKmOscNuECCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 }, // kSKmOscNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCoth
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuEBarNCcoh 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_e
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 }, // kSKmOscNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCoth    
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
            }//1Re RHC
            ,
            // > Multiple muon-like ring 
            {
                // rpa 
                // banffmec banffca5 banffbgres
                // banffmaqe banffmares banffscca banffsccv 
                // banffpf banffccothshape banffebo banffcccoh banffnccoh banffncoth banffnue/numu banffNC1gamma
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC0pi  
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuENC1Gamma
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1Gamma
                
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC0pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
            }//MultiRmu
            ,
            // > Multiple electron-like ring 
            {
                // rpa 
                // banffmec banffca5 banffbgres
                // banffmaqe banffmares banffscca banffsccv 
                // banffpf banffccothshape banffebo banffcccoh banffnccoh banffncoth banffnue/numu banffNC1gamma
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuNCcoh 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuNC1Gamma
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuMuBarNC1Gamma
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 }, // kSKmNuECCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 }, // kSKmNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, // kSKmNuECCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuENCcoh 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuENC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 }, // kSKmOscNuECCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuECC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 }, // kSKmOscNuECCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuECCoth
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 }, // kSKmNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, // kSKmNuEBarCCoth
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1pi0 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, // kSKmNuEBarNCcoh 
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_e
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 }, // kSKmNuEBarNC1Gamma
                
                {  1,   0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCQE        
                {  0,   0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCC1pi
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0 }, // kSKmOscNuEBarCCcoh
                {  0,   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCMEC
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 }, // kSKmOscNuEBarCCoth 
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
                {  0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
            }//MultiRe
#endif  
        };
        
#endif //BANFF_BINS
        
    } // conv
} // valor

#endif 
