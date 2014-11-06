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
#if BANFF_BINS == 1
    static const int kNumOfXSecSyst = 14; //21; // NEW XSEC 
    static const int kNumOfRHCXSecSyst = 14; //21; // NEW XSEC
#elif BANFF_BINS == 2
    static const int kNumOfXSecSyst = 14; //19; // NEW XSEC
    static const int kNumOfRHCXSecSyst = 14; //19; // NEW XSEC
#elif BANFF_BINS == 3
    static const int kNumOfXSecSyst = 14; //15; // NEW XSEC
    static const int kNumOfRHCXSecSyst = 14; //15; // NEW XSEC
#endif

#ifdef __NUMUONLY__
    static const int kNumOfUncXSecSyst = 0; //11; // NEW XSEC (no uncorr xsec)
    static const int kNumOfRHCUncXSecSyst = 0; //11; // NEW XSEC (no uncorr xsec)
#else
#ifdef __NUEONLY__
    static const int kNumOfUncXSecSyst = 0; //9; // NEW XSEC (no uncorr xsec)
    static const int kNumOfRHCUncXSecSyst = 0; //9; // NEW XSEC (no uncorr xsec)
#else
    static const int kNumOfUncXSecSyst = 0; //12; // NEW XSEC (no uncorr xsec)
    static const int kNumOfRHCUncXSecSyst = 0; //12; // NEW XSEC (no uncorr xsec)
#endif
#endif

    // The following is the list of the above kNumOfXSecSyst parameters
    static const PhysParamId_t kListOfXSecSyst[kNumOfXSecSyst] = 
      {
	/* 	// NEW XSEC  */

	/* 	kNuisBANFFMAQE,             // BANFF, M_{A}  QE */
	/* 	kNuisBANFFMARES,            // BANFF, M_{A}  RES */
	/* 	kNuisBANFFCCQEE0,           // BANFF, CCQE   normalization, 0 - 1.5 GeV */
	/* #if BANFF_BINS != 3 */
	/* 	kNuisBANFFCCQEE1,           // BANFF, CCQE   normalization, 1.5 - 3.5 GeV */
	/* 	kNuisBANFFCCQEE2,           // BANFF, CCQE   normalization, 3.5 - 30.0 GeV */
	/* #endif */
	/* 	kNuisBANFFCC1piE0,          // BANFF, CC1pi  normalization, 0 - 2.5 GeV */
	/* #if BANFF_BINS != 3 */
	/* 	kNuisBANFFCC1piE1,          // BANFF, CC1pi  normalization, 2.5 - 30.0 GeV */
	/* #endif */
	/* #if (BANFF_BINS != 2) */
	/* 	kNuisBANFFNC1pi0,           // BANFF, NC1pi0 normalization, 0.0 - 30.0 GeV */
	/* #endif */
	/* 	kNuisXSecFermiMom,          // Fermi momentum */
	/* 	kNuisXSecSpecFunc,          // spectral function on/off */
	/* 	kNuisXSecWshape,            // W-shape */
	/* 	kNuisXSecShapeCCoth,        // CC other shape */
	/* 	kNuisXSecPiLessDelta,       // Pion-less delta decay */
	/* #ifndef __NUEONLY__ */
	/* 	kNuisXSecBindE,             // Binding energy */
	/* #endif */
	/* 	kNuisXSecNormCCcoh,         // CC coherent norm */
	/* #ifndef __NUEONLY__ */
	/* 	kNuisXSecNormNC1piC,        // NC 1piC norm */
	/* #endif */
	/* #ifndef __NUMUONLY__ */
	/* 	kNuisXSecNormNCcoh,         // NC coh norm */
	/* #endif */
	/* 	kNuisXSecNormNCoth,         // NC other norm */
	/* #ifndef __NUEONLY__ */
	/* 	kNuisXSecNormAntiNuToNu,    // antinu / nu cross section */
	/* #endif */
	/* 	kNuisXSecNormNuEToNuMu,     // nue / numu cross section */
	/* 	kNuisRHCToFHC               // RHC / FHC normalisation */
  
	// NEW XSEC
	kNuisBANFFMEC,              // BANFF, MEC normalization
	kNuisBANFFCA5,              // BANFF, CA5 RES
	kNuisBANFFBgRES,            // BANFF, Isospin 1/2 normalization
	kNuisBANFFMAQE,             // BANFF, M_{A}  QE
	kNuisBANFFMARES,            // BANFF, M_{A}  RES
	kNuisBANFFSCCA,             // BANFF, SCCA
	kNuisBANFFSCCV,             // BANFF, SCCV
	kNuisXSecFermiMom,          // Fermi momentum
	kNuisXSecShapeCCoth,        // CC other shape
	kNuisXSecBindE,             // Binding energy
	kNuisXSecNormCCcoh,         // CC coherent norm
	kNuisXSecNormNCcoh,         // NC coh norm
	kNuisXSecNormNCoth,         // NC other norm
	kNuisXSecNormNuEToNuMu,     // nue / numu cross section
      };

    static const PhysParamId_t kListOfRHCXSecSyst[kNumOfRHCXSecSyst] = 
      {
	/* 	// NEW XSEC */

	/* 	kRHCNuisBANFFMAQE,             // BANFF, M_{A}  QE */
	/*         kRHCNuisBANFFMARES,            // BANFF, M_{A}  RES */
	/*         kRHCNuisBANFFCCQEE0,           // BANFF, CCQE   normalization, 0 - 1.5 GeV */
	/* #if BANFF_BINS != 3 */
	/*         kRHCNuisBANFFCCQEE1,           // BANFF, CCQE   normalization, 1.5 - 3.5 GeV */
	/*         kRHCNuisBANFFCCQEE2,           // BANFF, CCQE   normalization, 3.5 - 30.0 GeV */
	/* #endif */
	/*         kRHCNuisBANFFCC1piE0,          // BANFF, CC1pi  normalization, 0 - 2.5 GeV */
	/* #if BANFF_BINS != 3 */
	/*         kRHCNuisBANFFCC1piE1,          // BANFF, CC1pi  normalization, 2.5 - 30.0 GeV */
	/* #endif */
	/* #if (BANFF_BINS != 2) */
	/*         kRHCNuisBANFFNC1pi0,           // BANFF, NC1pi0 normalization, 0.0 - 30.0 GeV */
	/* #endif */
	/*         kRHCNuisXSecFermiMom,          // Fermi momentum */
	/*         kRHCNuisXSecSpecFunc,          // spectral function on/off */
	/*         kRHCNuisXSecWshape,            // W-shape */
	/*         kRHCNuisXSecShapeCCoth,        // CC other shape */
	/*         kRHCNuisXSecPiLessDelta,       // Pion-less delta decay */
	/* #ifndef __NUEONLY__ */
	/*         kRHCNuisXSecBindE,             // Binding energy */
	/* #endif */
	/*         kRHCNuisXSecNormCCcoh,         // CC coherent norm */
	/* #ifndef __NUEONLY__ */
	/*         kRHCNuisXSecNormNC1piC,        // NC 1piC norm */
	/* #endif */
	/* #ifndef __NUMUONLY__ */
	/*         kRHCNuisXSecNormNCcoh,         // NC coh norm */
	/* #endif */
	/*         kRHCNuisXSecNormNCoth,         // NC other norm */
	/* #ifndef __NUEONLY__ */
	/*         kRHCNuisXSecNormAntiNuToNu,    // antinu / nu cross section */
	/* #endif */
	/*         kRHCNuisXSecNormNuEToNuMu,     // nue / numu cross section */
	/*         kNuisRHCToFHC,                 // RHC / FHC normalisation */
    
	// NEW XSEC
	kRHCNuisBANFFMEC,              // BANFF, MEC normalization
	kRHCNuisBANFFCA5,              // BANFF, CA5 RES
	kRHCNuisBANFFBgRES,            // BANFF, Isospin 1/2 normalization
	kRHCNuisBANFFMAQE,             // BANFF, M_{A}  QE
	kRHCNuisBANFFMARES,            // BANFF, M_{A}  RES
	kRHCNuisBANFFSCCA,             // BANFF, SCCA
	kRHCNuisBANFFSCCV,             // BANFF, SCCV
	kRHCNuisXSecFermiMom,          // Fermi momentum
	kRHCNuisXSecShapeCCoth,        // CC other shape
	kRHCNuisXSecBindE,             // Binding energy
	kRHCNuisXSecNormCCcoh,         // CC coherent norm
	kRHCNuisXSecNormNCcoh,         // NC coh norm
	kRHCNuisXSecNormNCoth,         // NC other norm
	kRHCNuisXSecNormNuEToNuMu,     // nue / numu cross section
      };
 
 
    // The following list specifies which cross-section related systematic has its
    // effect on MC templates pre-computed and stored in cubic splines.
    // The order of systematic parameters should match the one in kListOfXSecSyst.
    static const bool kXSecSystHasWghtSplines[kNumOfXSecSyst] = 
      {
	/* 	// NEW XSEC */
	/* 	true,                      // kNuisBANFFMAQE */
	/* 	true,                      // kNuisBANFFMARES */
	/* 	false,                     // kNuisBANFFCCQEE0 */
	/* #if BANFF_BINS != 3 */
	/* 	false,                     // kNuisBANFFCCQEE1 */
	/* 	false,                     // kNuisBANFFCCQEE2 */
	/* #endif */
	/* 	false,                     // kNuisBANFFCC1piE0 */
	/* #if BANFF_BINS != 3 */
	/* 	false,                     // kNuisBANFFCC1piE1 */
	/* #endif */
	/* #if (BANFF_BINS != 2) */
	/* 	false,                     // kNuisBANFFNC1pi0 */
	/* #endif */
	/* 	true,                      // Fermi momentum */
	/* 	true,                      // spectral function on/off */
	/* 	true,                      // W-shape */
	/* 	true,                      // CC other shape */
	/* 	true,                      // Pion-less delta decay */
	/* #ifndef __NUEONLY__ */
	/* 	true,                      // Binding energy */
	/* #endif */
	/* 	false,                     // CC coherent norm */
	/* #ifndef __NUEONLY__ */
	/* 	false,                     // NC 1piC norm */
	/* #endif */
	/* #ifndef __NUMUONLY__ */
	/* 	false,                     // NC coh norm */
	/* #endif */
	/* 	false,                     // NC other norm */
	/* #ifndef __NUEONLY__ */
	/* 	false,                     // antinu / nu cross section */
	/* #endif */
	/* 	false,                     // nue / numu cross section */
	/* 	false                      // RHC / FHC normalisation */
  

	// NEW XSEC
	false,                     // kNuisBANFFMEC
	true,                      // kNuisBANFFCA5
	false,                     // kNuisBANFFBgRES
	true,                      // kNuisBANFFMAQE
	true,                      // kNuisBANFFMARES
	true,                      // kNuisBANFFSCCA
	true,                      // kNuisBANFFSCCV

	true,                      // Fermi momentum
	true,                      // CC other shape
	true,                      // Binding energy
	false,                     // CC coherent norm
	false,                     // NC coh norm
	false,                     // NC other norm
	false,                     // nue / numu cross section
      };

    static const bool kRHCXSecSystHasWghtSplines[kNumOfRHCXSecSyst] = 
      {
	/* 	// NEW XSEC */
	/* 	true,                      // kNuisBANFFMAQE */
	/* 	true,                      // kNuisBANFFMARES */
	/* 	false,                     // kNuisBANFFCCQEE0 */
	/* #if BANFF_BINS != 3 */
	/* 	false,                     // kNuisBANFFCCQEE1 */
	/* 	false,                     // kNuisBANFFCCQEE2 */
	/* #endif */
	/* 	false,                     // kNuisBANFFCC1piE0 */
	/* #if BANFF_BINS != 3 */
	/* 	false,                     // kNuisBANFFCC1piE1 */
	/* #endif */
	/* #if (BANFF_BINS != 2) */
	/* 	false,                     // kNuisBANFFNC1pi0 */
	/* #endif */
	/* 	true,                      // Fermi momentum */
	/* 	true,                      // spectral function on/off */
	/* 	true,                      // W-shape */
	/* 	true,                      // CC other shape */
	/* 	true,                      // Pion-less delta decay */
	/* #ifndef __NUEONLY__ */
	/* 	true,                      // Binding energy */
	/* #endif */
	/* 	false,                     // CC coherent norm */
	/* #ifndef __NUEONLY__ */
	/* 	false,                     // NC 1piC norm */
	/* #endif */
	/* #ifndef __NUMUONLY__ */
	/* 	false,                     // NC coh norm */
	/* #endif */
	/* 	false,                     // NC other norm */
	/* #ifndef __NUEONLY__ */
	/* 	false,                     // antinu / nu cross section */
	/* #endif */
	/* 	false,                     // nue / numu cross section */
	/* 	false                      // RHC / FHC normalisation */


	// NEW XSEC
	false,                     // kRHCNuisBANFFMEC
	true,                      // kRHCNuisBANFFCA5
	false,                     // kRHCNuisBANFFBgRES
	true,                      // kNuisBANFFMAQE
	true,                      // kNuisBANFFMARES
	true,                      // kNuisBANFFSCCA             
	true,                      // kNuisBANFFSCCV             

	true,                      // Fermi momentum
	true,                      // CC other shape
	true,                      // Binding energy
	false,                     // CC coherent norm
	false,                     // NC coh norm
	false,                     // NC other norm
	false,                     // nue / numu cross section
      };

    // The following list specifies which cross-section related systematic 
    // parameters were included in the BANFF fit.
    // The order of systematic parameters should match the one in kListOfXSecSyst.
    static const bool kXSecSystIncludedInBANFF[kNumOfXSecSyst] = 
      {
	/* 	// NEW XSEC */
	/* 	true,                      // kNuisBANFFMAQE */
	/* 	true,                      // kNuisBANFFMARES */
	/* 	true,                      // kNuisBANFFCCQEE0 */
	/* #if BANFF_BINS != 3 */
	/* 	true,                      // kNuisBANFFCCQEE1 */
	/* 	true,                      // kNuisBANFFCCQEE2 */
	/* #endif */
	/* 	true,                      // kNuisBANFFCC1piE0 */
	/* #if BANFF_BINS != 3 */
	/* 	true,                      // kNuisBANFFCC1piE1 */
	/* #endif */
	/* #if (BANFF_BINS != 2) */
	/* 	true,                      // kNuisBANFFNC1pi0 */
	/* #endif */
	/* 	false,                     // Fermi momentum */
	/* 	false,                     // spectral function on/off */
	/* 	false,                     // W-shape */
	/* 	false,                     // CC other shape */
	/* 	false,                     // Pion-less delta decay */
	/* #ifndef __NUEONLY__ */
	/* 	false,                     // Binding energy */
	/* #endif */
	/* 	false,                     // CC coherent norm */
	/* #ifndef __NUEONLY__ */
	/* 	false,                     // NC 1piC norm */
	/* #endif */
	/* #ifndef __NUMUONLY__ */
	/* 	false,                     // NC coh norm */
	/* #endif */
	/* 	false,                     // NC other norm */
	/* #ifndef __NUEONLY__ */
	/* 	false,                     // antinu / nu cross section */
	/* #endif */
	/* 	false,                     // nue / numu cross section */
	/* 	false                      // RHC / FHC normalisation */

	// NEW XSEC
	true,                      // kNuisBANFFMEC
	true,                      // kNuisBANFFCA5
	true,                      // kNuisBANFFBgRES
	true,                      // kNuisBANFFMAQE
	true,                      // kNuisBANFFMARES
	true,                      // kNuisBANFFSCCA             
	true,                      // kNuisBANFFSCCV             

	true,                     // Fermi momentum // constrained from banff in O
	true,                     // CC other shape // constrained from banff in O
	true,                     // Binding energy
	true,                     // CC coherent norm // constrained from banff in O
	true,                     // NC coh norm // constrained from banff in O
	true,                     // NC other norm
	true,                     // nue / numu cross section
      };

    static const bool kRHCXSecSystIncludedInBANFF[kNumOfRHCXSecSyst] = 
      {
	/* 	// NEW XSEC */
	/* 	true,                      // kNuisBANFFMAQE */
	/* 	true,                      // kNuisBANFFMARES */
	/* 	true,                      // kNuisBANFFCCQEE0 */
	/* #if BANFF_BINS != 3 */
	/* 	true,                      // kNuisBANFFCCQEE1 */
	/* 	true,                      // kNuisBANFFCCQEE2 */
	/* #endif */
	/* 	true,                      // kNuisBANFFCC1piE0 */
	/* #if BANFF_BINS != 3 */
	/* 	true,                      // kNuisBANFFCC1piE1 */
	/* #endif */
	/* #if (BANFF_BINS != 2) */
	/* 	true,                      // kNuisBANFFNC1pi0 */
	/* #endif */
	/* 	false,                     // Fermi momentum */
	/* 	false,                     // spectral function on/off */
	/* 	false,                     // W-shape */
	/* 	false,                     // CC other shape */
	/* 	false,                     // Pion-less delta decay */
	/* #ifndef __NUEONLY__ */
	/* 	false,                     // Binding energy */
	/* #endif */
	/* 	false,                     // CC coherent norm */
	/* #ifndef __NUEONLY__ */
	/* 	false,                     // NC 1piC norm */
	/* #endif */
	/* #ifndef __NUMUONLY__ */
	/* 	false,                     // NC coh norm */
	/* #endif */
	/* 	false,                     // NC other norm */
	/* #ifndef __NUEONLY__ */
	/* 	false,                     // antinu / nu cross section */
	/* #endif */
	/* 	false,                     // nue / numu cross section */
	/* 	false                      // RHC / FHC normalisation */

	// NEW XSEC
	true,                      // kNuisBANFFMEC
	true,                      // kNuisBANFFCA5
	true,                      // kNuisBANFFBgRES
	true,                      // kNuisBANFFMAQE
	true,                      // kNuisBANFFMARES
	true,                      // kNuisBANFFSCCA             
	true,                      // kNuisBANFFSCCV             

	true,                     // Fermi momentum // constrained from banff in O
	true,                     // CC other shape // constrained from banff in O
	true,                     // Binding energy
	true,                     // CC coherent norm // constrained from banff in O
	true,                     // NC coh norm // constrained from banff in O
	true,                     // NC other norm
	true,                     // nue / numu cross section
      };

    // Ereco and Etrue cuts for the application of cross-section systematics.
    // The order of systematic parameters should match the one in kListOfXSecSyst.
    static const double kXSecSystErecoCut[kNumOfXSecSyst][2] =
      {
	/* 	// NEW XSEC */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFMAQE */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFMARES */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFCCQEE0 */
	/* #if BANFF_BINS != 3 */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFCCQEE1 */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFCCQEE2 */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFCC1piE0 */
	/* #if BANFF_BINS != 3 */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFCC1piE1 */
	/* #endif */
	/* #if (BANFF_BINS != 2) */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFNC1pi0 */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // Fermi momentum */
	/* 	{ 0.0, 100.0 },            // spectral function on/off */
	/* 	{ 0.0, 100.0 },            // W-shape */
	/* 	{ 0.0, 100.0 },            // CC other shape */
	/* 	{ 0.0, 100.0 },            // Pion-less delta decay */
	/* #ifndef __NUEONLY__ */
	/* 	{ 0.0, 100.0 },            // Binding energy */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // CC coherent norm */
	/* #ifndef __NUEONLY__ */
	/* 	{ 0.0, 100.0 },            // NC 1piC norm */
	/* #endif */
	/* #ifndef __NUMUONLY__ */
	/* 	{ 0.0, 100.0 },            // NC coh norm */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // NC other norm */
	/* #ifndef __NUEONLY__ */
	/* 	{ 0.0, 100.0 },            // antinu / nu cross section */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // nue / numu cross section */
	/* 	{ 0.0, 100.0 }             // RHC / FHC normalisation */

	// NEW XSEC
	{ 0.0, 100.0 },            // kNuisBANFFMEC
	{ 0.0, 100.0 },            // kNuisBANFFCA5
	{ 0.0, 100.0 },            // kNuisBANFFBgRES
	{ 0.0, 100.0 },            // kNuisBANFFMAQE
	{ 0.0, 100.0 },            // kNuisBANFFMARES
	{ 0.0, 100.0 },            // kNuisBANFFSCCA             
	{ 0.0, 100.0 },            // kNuisBANFFSCCV             

	{ 0.0, 100.0 },            // Fermi momentum
	{ 0.0, 100.0 },            // CC other shape
	{ 0.0, 100.0 },            // Binding energy
	{ 0.0, 100.0 },            // CC coherent norm
	{ 0.0, 100.0 },            // NC coh norm
	{ 0.0, 100.0 },            // NC other norm
	{ 0.0, 100.0 },            // nue / numu cross section
      };

    static const double kRHCXSecSystErecoCut[kNumOfRHCXSecSyst][2] =
      {
	/* 	// NEW XSEC */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFMAQE */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFMARES */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFCCQEE0 */
	/* #if BANFF_BINS != 3 */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFCCQEE1 */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFCCQEE2 */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFCC1piE0 */
	/* #if BANFF_BINS != 3 */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFCC1piE1 */
	/* #endif */
	/* #if (BANFF_BINS != 2) */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFNC1pi0 */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // Fermi momentum */
	/* 	{ 0.0, 100.0 },            // spectral function on/off */
	/* 	{ 0.0, 100.0 },            // W-shape */
	/* 	{ 0.0, 100.0 },            // CC other shape */
	/* 	{ 0.0, 100.0 },            // Pion-less delta decay */
	/* #ifndef __NUEONLY__ */
	/* 	{ 0.0, 100.0 },            // Binding energy */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // CC coherent norm */
	/* #ifndef __NUEONLY__ */
	/* 	{ 0.0, 100.0 },            // NC 1piC norm */
	/* #endif */
	/* #ifndef __NUMUONLY__ */
	/* 	{ 0.0, 100.0 },            // NC coh norm */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // NC other norm */
	/* #ifndef __NUEONLY__ */
	/* 	{ 0.0, 100.0 },            // antinu / nu cross section */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // nue / numu cross section */
	/* 	{ 0.0, 100.0 }             // RHC / FHC normalisation */

	// NEW XSEC
	{ 0.0, 100.0 },            // kNuisBANFFMEC
	{ 0.0, 100.0 },            // kNuisBANFFCA5
	{ 0.0, 100.0 },            // kNuisBANFFBgRES
	{ 0.0, 100.0 },            // kNuisBANFFMAQE
	{ 0.0, 100.0 },            // kNuisBANFFMARES
	{ 0.0, 100.0 },            // kNuisBANFFSCCA             
	{ 0.0, 100.0 },            // kNuisBANFFSCCV             

	{ 0.0, 100.0 },            // Fermi momentum
	{ 0.0, 100.0 },            // CC other shape
	{ 0.0, 100.0 },            // Binding energy
	{ 0.0, 100.0 },            // CC coherent norm
	{ 0.0, 100.0 },            // NC coh norm
	{ 0.0, 100.0 },            // NC other norm
	{ 0.0, 100.0 },            // nue / numu cross section
      };

    static const double kXSecSystEtrueCut[kNumOfXSecSyst][2] =
      {
	/* 	// NEW XSEC */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFMAQE */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFMARES */
	/* 	{ 0.0,   1.5 },            // kNuisBANFFCCQEE0 */
	/* #if BANFF_BINS != 3 */
	/* 	{ 1.5,   3.5 },            // kNuisBANFFCCQEE1 */
	/* 	{ 3.5,  30.0 },            // kNuisBANFFCCQEE2 */
	/* #endif */
	/* 	{ 0.0,   2.5 },            // kNuisBANFFCC1piE0 */
	/* #if BANFF_BINS != 3 */
	/* 	{ 2.5,  30.0 },            // kNuisBANFFCC1piE1 */
	/* #endif */
	/* #if (BANFF_BINS != 2) */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFNC1pi0 */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // Fermi momentum */
	/* 	{ 0.0, 100.0 },            // spectral function on/off */
	/* 	{ 0.0, 100.0 },            // W-shape */
	/* 	{ 0.0, 100.0 },            // CC other shape */
	/* 	{ 0.0, 100.0 },            // Pion-less delta decay */
	/* #ifndef __NUEONLY__ */
	/* 	{ 0.0, 100.0 },            // Binding energy */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // CC coherent norm */
	/* #ifndef __NUEONLY__ */
	/* 	{ 0.0, 100.0 },            // NC 1piC norm */
	/* #endif */
	/* #ifndef __NUMUONLY__ */
	/* 	{ 0.0, 100.0 },            // NC coh norm */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // NC other norm */
	/* #ifndef __NUEONLY__ */
	/* 	{ 0.0, 100.0 },            // antinu / nu cross section */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // nue / numu cross section */
	/* 	{ 0.0, 100.0 }             // RHC / FHC normalisation */


	// NEW XSEC
	{ 0.0, 100.0 },            // kNuisBANFFMEC
	{ 0.0, 100.0 },            // kNuisBANFFCA5
	{ 0.0, 100.0 },            // kNuisBANFFBgRES
	{ 0.0, 100.0 },            // kNuisBANFFMAQE
	{ 0.0, 100.0 },            // kNuisBANFFMARES
	{ 0.0, 100.0 },            // kNuisBANFFSCCA             
	{ 0.0, 100.0 },            // kNuisBANFFSCCV             

	{ 0.0, 100.0 },            // Fermi momentum
	{ 0.0, 100.0 },            // CC other shape
	{ 0.0, 100.0 },            // Binding energy
	{ 0.0, 100.0 },            // CC coherent norm
	{ 0.0, 100.0 },            // NC coh norm
	{ 0.0, 100.0 },            // NC other norm
	{ 0.0, 100.0 },            // nue / numu cross section
      };

    static const double kRHCXSecSystEtrueCut[kNumOfRHCXSecSyst][2] =
      {
	/* 	// NEW XSEC */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFMAQE */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFMARES */
	/* 	{ 0.0,   1.5 },            // kNuisBANFFCCQEE0 */
	/* #if BANFF_BINS != 3 */
	/* 	{ 1.5,   3.5 },            // kNuisBANFFCCQEE1 */
	/* 	{ 3.5,  30.0 },            // kNuisBANFFCCQEE2 */
	/* #endif */
	/* 	{ 0.0,   2.5 },            // kNuisBANFFCC1piE0 */
	/* #if BANFF_BINS != 3 */
	/* 	{ 2.5,  30.0 },            // kNuisBANFFCC1piE1 */
	/* #endif */
	/* #if (BANFF_BINS != 2) */
	/* 	{ 0.0, 100.0 },            // kNuisBANFFNC1pi0 */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // Fermi momentum */
	/* 	{ 0.0, 100.0 },            // spectral function on/off */
	/* 	{ 0.0, 100.0 },            // W-shape */
	/* 	{ 0.0, 100.0 },            // CC other shape */
	/* 	{ 0.0, 100.0 },            // Pion-less delta decay */
	/* #ifndef __NUEONLY__ */
	/* 	{ 0.0, 100.0 },            // Binding energy */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // CC coherent norm */
	/* #ifndef __NUEONLY__ */
	/* 	{ 0.0, 100.0 },            // NC 1piC norm */
	/* #endif */
	/* #ifndef __NUMUONLY__ */
	/* 	{ 0.0, 100.0 },            // NC coh norm */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // NC other norm */
	/* #ifndef __NUEONLY__ */
	/* 	{ 0.0, 100.0 },            // antinu / nu cross section */
	/* #endif */
	/* 	{ 0.0, 100.0 },            // nue / numu cross section */
	/* 	{ 0.0, 100.0 }             // RHC / FHC normalisation */



	// NEW XSEC
	{ 0.0, 100.0 },            // kNuisBANFFMEC
	{ 0.0, 100.0 },            // kNuisBANFFCA5
	{ 0.0, 100.0 },            // kNuisBANFFBgRES
	{ 0.0, 100.0 },            // kNuisBANFFMAQE
	{ 0.0, 100.0 },            // kNuisBANFFMARES
	{ 0.0, 100.0 },            // kNuisBANFFSCCA             
	{ 0.0, 100.0 },            // kNuisBANFFSCCV             

	{ 0.0, 100.0 },            // Fermi momentum
	{ 0.0, 100.0 },            // CC other shape
	{ 0.0, 100.0 },            // Binding energy
	{ 0.0, 100.0 },            // CC coherent norm
	{ 0.0, 100.0 },            // NC coh norm
	{ 0.0, 100.0 },            // NC other norm
	{ 0.0, 100.0 },            // nue / numu cross section
      };


    //Fractional error of systematic parameters, as of BANFF 2012v8
    // ONLY USE FOR: macros/misc/test_fitter_inputs/test_valor_kdi_splines.C
    static const double kXSecSystFracErr[kNumOfXSecSyst] =
      {
	/* 	// NEW XSEC */
	/* 	0.0716,          // BANFF, M_{A}QE */
	/* 	0.0675,          // BANFF, M_{A}RES */
	/* 	0.09,           // BANFF, CCQE   normalization, 0 - 1.5 GeV */
	/* #if BANFF_BINS != 3 */
	/* 	0.3,            // BANFF, CCQE   normalization, 1.5 - 3.5 GeV */
	/* 	0.3,            // BANFF, CCQE   normalization, 3.5 - 30.0 GeV */
	/* #endif */
	/* 	0.149,          // BANFF, CC1pi  normalization, 0 - 2.5 GeV */
	/* #if BANFF_BINS != 3 */
	/* 	0.4,            // BANFF, CC1pi  normalization, 2.5 - 30.0 GeV */
	/* #endif */
	/* #if (BANFF_BINS != 2) */
	/* 	0.3,            // BANFF, NC1pi0 normalization, 0.0 - 30.0 GeV */
	/* #endif */
	/* 	0.13333,        // Fermi momentum */
	/* 	1.0,            // spectral function on/off */
	/* 	0.516,          // W-shape */
	/* 	0.4,            // CC other shape */
	/* 	0.2,            // Pion-less delta decay */
	/* #ifndef __NUEONLY__ */
	/* 	0.33333,        // Binding enery */
	/* #endif */
	/* 	1.0,            // CC coherent norm */
	/* #ifndef __NUEONLY__ */
	/* 	0.3,            // NC 1piC norm */
	/* #endif */
	/* #ifndef __NUMUONLY__ */
	/* 	0.3,            // NC coh norm */
	/* #endif */
	/* 	0.3,            // NC other norm */
	/* #ifndef __NUEONLY__ */
	/* 	0.4,            // antinu / nu cross section */
	/* #endif */
	/* 	0.03,           // nue / numu cross section */
	/* 	0.10            // RHC / FHC normalisation */


	// NEW XSEC
	1.0,            // BANFF, MEC
	1.0,            // BANFF, CA5
	1.0,            // BANFF, BgRES
	0.0716,         // BANFF, M_{A}QE
	0.0675,         // BANFF, M_{A}RES
	1.0,            // BANFF, SCCA
	1.0,            // BANFF, SCCV

	0.13333,        // Fermi momentum
	0.4,            // CC other shape
	0.33333,        // Binding enery
	1.0,            // CC coherent norm
	0.3,            // NC coh norm
	0.3,            // NC other norm
	0.03,           // nue / numu cross section
      };

    static const double kRHCXSecSystFracErr[kNumOfRHCXSecSyst] =
      {
	/* 	// NEW XSEC */
	/* 	0.0716,         // BANFF, M_{A}QE */
	/* 	0.0675,         // BANFF, M_{A}RES */
	/* 	0.09,           // BANFF, CCQE   normalization, 0 - 1.5 GeV */
	/* #if BANFF_BINS != 3 */
	/* 	0.3,            // BANFF, CCQE   normalization, 1.5 - 3.5 GeV */
	/* 	0.3,            // BANFF, CCQE   normalization, 3.5 - 30.0 GeV */
	/* #endif */
	/* 	0.149,          // BANFF, CC1pi  normalization, 0 - 2.5 GeV */
	/* #if BANFF_BINS != 3 */
	/* 	0.4,            // BANFF, CC1pi  normalization, 2.5 - 30.0 GeV */
	/* #endif */
	/* #if (BANFF_BINS != 2) */
	/* 	0.3,            // BANFF, NC1pi0 normalization, 0.0 - 30.0 GeV */
	/* #endif */
	/* 	0.13333,        // Fermi momentum */
	/* 	1.0,            // spectral function on/off */
	/* 	0.516,          // W-shape */
	/* 	0.4,            // CC other shape */
	/* 	0.2,            // Pion-less delta decay */
	/* #ifndef __NUEONLY__ */
	/* 	0.33333,        // Binding enery */
	/* #endif */
	/* 	1.0,            // CC coherent norm */
	/* #ifndef __NUEONLY__ */
	/* 	0.3,            // NC 1piC norm */
	/* #endif */
	/* #ifndef __NUMUONLY__ */
	/* 	0.3,            // NC coh norm */
	/* #endif */
	/* 	0.3,            // NC other norm */
	/* #ifndef __NUEONLY__ */
	/* 	0.4,            // antinu / nu cross section */
	/* #endif */
	/* 	0.03,           // nue / numu cross section */
	/* 	0.10            // RHC / FHC normalisation */


	// NEW XSEC
	0.1,            // BANFF, MEC
	0.1,            // BANFF, CA5
	0.1,            // BANFF, BgRES
	0.0716,         // BANFF, M_{A}QE
	0.0675,         // BANFF, M_{A}RES
	0.1,            // BANFF, SCCA
	0.1,            // BANFF, SCCV

	0.13333,        // Fermi momentum
	0.4,            // CC other shape
	0.33333,        // Binding enery
	1.0,            // CC coherent norm
	0.3,            // NC coh norm
	0.3,            // NC other norm
	0.03,           // nue / numu cross section
      };

    // The following list describles in which MC templates we are supposed to apply the
    // weights as a result of tweaking any of the above systematic parameters.
    // The order of systematic parameters should match the one in kListOfXSecSyst.
#if BANFF_BINS == 1
#ifndef HKHACK

    //Todo: make this work if kNumOfRHCXSecSyst is different from kNumOfXSecSyst.

    // NEW XSEC
    static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {

      // > 1 muon-like ring 
      {
	/* 	// MEC CA5 BgRES */
	/* 	// banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2 */
	/* 	// banffcc1piE0 banffcc1piE1 banffnc1pi0 */
	/* 	// pF  spectral   W-shape  ccothshape pilessdelta bindE */
	/* 	// cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu */
	/* 	// RHC/FHC */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1piC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_mu */
	/* #ifdef __USE_OSC_NUMU_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth */
	/* #endif */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1piC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu */
	/* #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuOscMuBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuOscMuBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0,  0 }, // kSKmNuOscMuBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuOscMuBarCCoth */
	/* #endif */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuENC1piC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_mu */

	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCoth */

	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarNC1piC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_mu */
	/* #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCoth */
	/* #endif */
	/* #ifndef __USE_OSC_NUMU_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }  // kSKmUndefined */




	// NEW XSEC

	// MEC CA5 BgRES banffmaqe banffmares banffscca banffsccv
	// pF  ccothshape bindE
	// cccoh  nccoh  ncoth  nue/numu
	{  0, 0, 0, 1, 0, 0, 0,  1, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 0 }, // kSKmNuMuCCoth
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNC1piC
	{  0, 0, 0, 0, 1, 0, 0,  0, 0, 0,   0, 0, 1, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0,  1, 0, 1,   0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 0 }, // kSKmOscNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 1, 0, 0, 0,  1, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 0 }, // kSKmNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCoth
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
	{  0, 0, 0, 0, 1, 0, 0,  0, 0, 0,   0, 0, 1, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0,  1, 0, 1,   0, 0, 0, 0 }, // kSKmNuOscMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmNuOscMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 0 }, // kSKmNuOscMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 0 }, // kSKmNuOscMuBarCCoth
#endif
	{  0, 0, 0, 1, 0, 1, 1,  1, 0, 1,   0, 0, 0, 1 }, // kSKmNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 1 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 1 }, // kSKmNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 1 }, // kSKmNuECCoth
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENC1piC
	{  0, 0, 0, 0, 1, 0, 0,  0, 0, 0,   0, 0, 1, 0 }, // kSKmNuENCoth_mu

	{  0, 0, 0, 1, 0, 1, 1,  1, 0, 1,   0, 0, 0, 1 }, // kSKmOscNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 1 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 1 }, // kSKmOscNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 1 }, // kSKmOscNuECCoth

	{  0, 0, 0, 1, 0, 1, 1,  1, 0, 1,   0, 0, 0, 1 }, // kSKmNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 1 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 1 }, // kSKmNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 1 }, // kSKmNuEBarCCoth
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarNC1piC
	{  0, 0, 0, 0, 1, 0, 0,  0, 0, 0,   0, 0, 1, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 1, 1,  1, 0, 1,   0, 0, 0, 1 }, // kSKmOscNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 1 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 1 }, // kSKmOscNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 1 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }  // kSKmUndefined
      }//1Rmu
      ,
      // > 1 electron-like ring 
      {
	/* 	// MEC CA5 BgRES */
	/* 	// banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2 */
	/* 	// banffcc1piE0 banffcc1piE1 banffnc1pi0 */
	/* 	// pF  spectral   W-shape  ccothshape pilessdelta bindE */
	/* 	// cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu */
	/* 	// RHC/FHC */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0  */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmNuMuNCcoh  */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_e */
	/* #ifdef __USE_OSC_NUMU_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOScNuMuCCoth */
	/* #endif */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1pi0  */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmNuMuBarNCcoh  */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_e */
	/* #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCoth */
	/* #endif */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuENC1pi0  */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmNuENCcoh  */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_e */

	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCoth */

	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuEBarC1pi0  */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmNuEBarNCcoh  */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_e */
	/* #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCoth */
	/* #endif */
	/* #ifndef __USE_OSC_NUMU_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */





      

	// NEW XSEC

	// MEC CA5 BgRES banffmaqe banffmares scca sccv
	// pF  ccothshape  bindE
	// cccoh  nccoh   ncoth  nue/numu
	{  0, 0, 0, 1, 0, 0, 0,  1, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 0 }, // kSKmNuMuCCoth
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 1, 0, 0 }, // kSKmNuMuNCcoh 
	{  0, 0, 0, 0, 1, 0, 0,  0, 0, 0,   0, 0, 1, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0,  1, 0, 1,   0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 0 }, // kSKmOscNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 0 }, // kSKmOScNuMuCCoth
#endif
	{  0, 0, 0, 1, 0, 0, 0,  1, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 0 }, // kSKmNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCoth
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 1, 0, 0 }, // kSKmNuMuBarNCcoh 
	{  0, 0, 0, 0, 1, 0, 0,  0, 0, 0,   0, 0, 1, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0,  1, 0, 1,   0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 1, 0, 1, 1,  1, 0, 1,   0, 0, 0, 1 }, // kSKmNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 1 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 1 }, // kSKmNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 1 }, // kSKmNuECCoth
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENC1pi0 
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 1, 0, 0 }, // kSKmNuENCcoh 
	{  0, 0, 0, 0, 1, 0, 0,  0, 0, 0,   0, 0, 1, 0 }, // kSKmNuENCoth_e

	{  0, 0, 0, 1, 0, 1, 1,  1, 0, 1,   0, 0, 0, 1 }, // kSKmOscNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 1 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 1 }, // kSKmOscNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 1 }, // kSKmOscNuECCoth

	{  0, 0, 0, 1, 0, 1, 1,  1, 0, 1,   0, 0, 0, 1 }, // kSKmNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 1 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 1 }, // kSKmNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 1 }, // kSKmNuEBarCCoth
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarC1pi0 
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 1, 0, 0 }, // kSKmNuEBarNCcoh 
	{  0, 0, 0, 0, 1, 0, 0,  0, 0, 0,   0, 0, 1, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 1, 1,  1, 0, 1,   0, 0, 0, 1 }, // kSKmOscNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,  0, 0, 0,   0, 0, 0, 1 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   1, 0, 0, 1 }, // kSKmOscNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,  0, 1, 0,   0, 0, 0, 1 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
      }//1Re
      ,
      // > 1 muon-like ring RHC
      {
	/* 	// MEC CA5 BgRES */
	/* 	// banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2 */
	/* 	// banffcc1piE0 banffcc1piE1 banffnc1pi0 */
	/* 	// pF  spectral   W-shape  ccothshape pilessdelta bindE */
	/* 	// cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu */
	/* 	// RHC/FHC */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1piC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_mu */
	/* #ifdef __USE_OSC_NUMU_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth */
	/* #endif */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1piC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu */
	/* #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuOscMuBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuOscMuBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0,  0 }, // kSKmNuOscMuBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuOscMuBarCCoth */
	/* #endif */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuENC1piC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_mu */

	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCoth */

	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarNC1piC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_mu */
	/* #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCoth */
	/* #endif */
	/* #ifndef __USE_OSC_NUMU_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }  // kSKmUndefined */
      



	// NEW XSEC

	// MEC CA5 BgRES banffmaqe banffmares scca sccv 
	// pF  ccothshape bindE
	// cccoh  nccoh  ncoth  nue/numu
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmNuMuCCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNC1piC
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmOscNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmNuOscMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuOscMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmNuOscMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmNuOscMuBarCCoth
#endif
	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmNuECCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENC1piC
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuENCoth_mu

	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmOscNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmOscNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmOscNuECCoth

	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmNuEBarCCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarNC1piC
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmOscNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmOscNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }  // kSKmUndefined



      }//1Rmu RHC
      ,
      // > 1 electron-like ring RHC
      {
	/* 	// MEC CA5 BgRES */
	/* 	// banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2 */
	/* 	// banffcc1piE0 banffcc1piE1 banffnc1pi0 */
	/* 	// pF  spectral   W-shape  ccothshape pilessdelta bindE */
	/* 	// cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu */
	/* 	// RHC/FHC */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0  */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmNuMuNCcoh  */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_e */
	/* #ifdef __USE_OSC_NUMU_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOScNuMuCCoth */
	/* #endif */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1pi0  */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmNuMuBarNCcoh  */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_e */
	/* #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCoth */
	/* #endif */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuENC1pi0  */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmNuENCcoh  */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_e */

	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCoth */

	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuEBarC1pi0  */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmNuEBarNCcoh  */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_e */
	/* #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCoth */
	/* #endif */
	/* #ifndef __USE_OSC_NUMU_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */



	// NEW XSEC

	// MEC CA5 BgRES banffmaqe banffmares scca sccv
	// pF  ccothshape  bindE
	// cccoh  nccoh   ncoth  nue/numu
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmNuMuCCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 1, 0, 0 }, // kSKmNuMuNCcoh 
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmOscNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmOScNuMuCCoth
#endif
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 1, 0, 0 }, // kSKmNuMuBarNCcoh 
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmNuECCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENC1pi0 
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 1, 0, 0 }, // kSKmNuENCcoh 
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuENCoth_e

	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmOscNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmOscNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmOscNuECCoth

	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmNuEBarCCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarC1pi0 
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 1, 0, 0 }, // kSKmNuEBarNCcoh 
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmOscNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmOscNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined



      }// 1Re RHC
      ,
      // > muon-like multiring
      {
	/* 	// MEC CA5 BgRES */
	/* 	// banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2 */
	/* 	// banffcc1piE0 banffcc1piE1 banffnc1pi0 */
	/* 	// pF  spectral   W-shape  ccothshape pilessdelta bindE */
	/* 	// cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu */
	/* 	// RHC/FHC */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1piC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_mu */
	/* #ifdef __USE_OSC_NUMU_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth */
	/* #endif */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1piC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu */
	/* #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuOscMuBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuOscMuBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0,  0 }, // kSKmNuOscMuBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuOscMuBarCCoth */
	/* #endif */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuENC1piC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_mu */

	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCoth */

	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarNC1piC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_mu */
	/* #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCoth */
	/* #endif */
	/* #ifndef __USE_OSC_NUMU_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }  // kSKmUndefined */



	// NEW XSEC

	// MEC CA5 BgRES banffmaqe banffmares scca sccv
	// pF  ccothshape bindE
	// cccoh  nccoh  ncoth  nue/numu
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmNuMuCCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNC1piC
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmOscNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmNuOscMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuOscMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmNuOscMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmNuOscMuBarCCoth
#endif
	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmNuECCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENC1piC
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuENCoth_mu

	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmOscNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmOscNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmOscNuECCoth

	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmNuEBarCCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarNC1piC
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmOscNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmOscNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }  // kSKmUndefined




      }//1Rmu
      ,
      // > electron-like multiring 
      {
	/* 	// MEC CA5 BgRES */
	/* 	// banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2 */
	/* 	// banffcc1piE0 banffcc1piE1 banffnc1pi0 */
	/* 	// pF  spectral   W-shape  ccothshape pilessdelta bindE */
	/* 	// cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu */
	/* 	// RHC/FHC */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0  */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmNuMuNCcoh  */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_e */
	/* #ifdef __USE_OSC_NUMU_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmOScNuMuCCoth */
	/* #endif */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1pi0  */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmNuMuBarNCcoh  */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_e */
	/* #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCoth */
	/* #endif */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmNuECCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuENC1pi0  */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmNuENCcoh  */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_e */

	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCoth */

	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCoth */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmNuEBarC1pi0  */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmNuEBarNCcoh  */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_e */
	/* #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 1, 0, 1, 1, 1,  0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCQE         */
	/* 	{  0, 1, 1, 0, 1, 0, 0, 0,  1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCC1pi */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCcoh */
	/* 	{  1, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // MEC */
	/* 	{  0, 0, 0, 0, 1, 0, 0, 0,  0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCoth */
	/* #endif */
	/* #ifndef __USE_OSC_NUMU_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__ */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* #endif */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */
	/* 	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined */




	// NEW XSEC

	// MEC CA5 BgRES banffmaqe banffmares scca sccv 
	// pF  ccothshape  bindE
	// cccoh  nccoh   ncoth  nue/numu
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmNuMuCCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 1, 0, 0 }, // kSKmNuMuNCcoh 
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmOscNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmOScNuMuCCoth
#endif
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 1, 0, 0 }, // kSKmNuMuBarNCcoh 
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0,   1, 0, 1,   0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmNuECCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENC1pi0 
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 1, 0, 0 }, // kSKmNuENCcoh 
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuENCoth_e

	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmOscNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmOscNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmOscNuECCoth

	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmNuEBarCCoth
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarC1pi0 
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 1, 0, 0 }, // kSKmNuEBarNCcoh 
	{  0, 0, 0, 0, 1, 0, 0,   0, 0, 0,   0, 0, 1, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 1, 1,   1, 0, 1,   0, 0, 0, 1 }, // kSKmOscNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0,   0, 0, 0,   0, 0, 0, 1 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   1, 0, 0, 1 }, // kSKmOscNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0,   0, 1, 0,   0, 0, 0, 1 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0,   0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined

      }// MR e

    };
#else //HKHACK 
    ///////////////////////////
    ///////HACK
    ///////////////////////////


    //Todo: make this work if kNumOfRHCXSecSyst is different from kNumOfXSecSyst.


    // NEW XSEC
    static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {

      // NEW XSEC
      
      /* 
      // > 1 muon-like ring 
      {
      // banff & spline sections all 0
      // numusignal numubackground nuesignal nuebackground, aka:
      // cccoh    nc1piC   nccoh   ncoth    0 0
      // RHC/FHC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1piC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_mu
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1piC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuOscMuBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuOscMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuOscMuBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuOscMuBarCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuECCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuECC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuECCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuECCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuENC1piC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuENCoth_mu

      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuECC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarNC1piC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmMEC
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1mu
      ,
      // > 1 electron-like ring 
      {
      // banff & spline sections all 0
      // numusignal numubackground nuesignal nuebackground, aka:
      // cccoh    nc1piC   nccoh   ncoth    0 0
      // RHC/FHC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNC1pi0 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCcoh 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_e
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOScNuMuCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNC1pi0 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_e
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuECCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuECC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuECCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuECCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENC1pi0 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENCcoh 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_e

      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarC1pi0 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCcoh 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_e
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuEBarCCoth
      #endif

      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuECC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuECCoth
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1e
      ,
      // > 1 muon-like ring RHC
      {
      // banff & spline sections all 0
      // numusignal numubackground nuesignal nuebackground, aka:
      // cccoh    nc1piC   nccoh   ncoth    0 0
      // RHC/FHC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1piC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_mu
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1piC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuOscMuBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuOscMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuOscMuBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuOscMuBarCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuECCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuECC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuECCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuECCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuENC1piC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuENCoth_mu

      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuECC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarNC1piC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmMEC
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1muRHC
      ,
      // > 1 electron-like ring RHC
      {
      // banff & spline sections all 0
      // numusignal numubackground nuesignal nuebackground, aka:
      // cccoh    nc1piC   nccoh   ncoth    0 0
      // RHC/FHC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNC1pi0 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCcoh 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_e
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOScNuMuCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNC1pi0 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_e
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuECCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuECC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuECCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuECCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENC1pi0 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENCcoh 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_e

      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarC1pi0 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCcoh 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_e
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuEBarCCoth
      #endif

      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuECC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuECCoth
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1eRHC  


      /// NEW MULTIRING
      ,
      // > Multi muon-like ring 
      {
      // banff & spline sections all 0
      // numusignal numubackground nuesignal nuebackground, aka:
      // cccoh    nc1piC   nccoh   ncoth    0 0
      // RHC/FHC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1piC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_mu
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1piC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0,  0 }, // kSKmNuOscMuBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuOscMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuOscMuBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuOscMuBarCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuECCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuECC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuECCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuECCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuENC1piC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuENCoth_mu

      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuECC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarNC1piC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0,  0 }, // kSKmMEC
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//MultiRmu
      ,
      // > Multi electron-like ring 
      {
      // banff & spline sections all 0
      // numusignal numubackground nuesignal nuebackground, aka:
      // cccoh    nc1piC   nccoh   ncoth    0 0
      // RHC/FHC
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNC1pi0 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCcoh 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_e
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOScNuMuCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNC1pi0 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_e
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuECCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuECC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuECCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuECCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENC1pi0 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENCcoh 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_e

      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarCCoth
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarC1pi0 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCcoh 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_e
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuEBarCCoth
      #endif

      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuECC1pi
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0,  0 }, // kSKmOscNuECCoth
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif 
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//MultiRe  

      */





      // > 1 muon-like ring 
      {
	// banff & spline sections all 0
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
      }//1mu
      ,
      // > 1 electron-like ring 
      {
	// banff & spline sections all 0
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOScNuMuCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi   
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#endif

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
      }//1e
      ,
      // > 1 muon-like ring RHC
      {
	// banff & spline sections all 0
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
      }//1muRHC
      ,
      // > 1 electron-like ring RHC
      {
	// banff & spline sections all 0
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOScNuMuCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#endif

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
      }//1eRHC  


      /// NEW MULTIRING
      ,
      // > Multi muon-like ring 
      {
	// banff & spline sections all 0
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth
	
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
      }//MultiRmu
      ,
      // > Multi electron-like ring 
      {
	// banff & spline sections all 0
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOScNuMuCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuECCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth
#endif

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmOscNuECCoth
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  // kSKmUndefined
      }//MultiRe  

    };
    //////////////////
    ///END HACK
    //////////////////
#endif //HKHACK
#elif BANFF_BINS == 2






    
    //Todo: make this work if kNumOfRHCXSecSyst is different from kNumOfXSecSyst.
 
    // NEW XSEC
    static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {

      /*
      // > 1 muon-like ring 
      {
      // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
      // banffcc1piE0 banffcc1piE1
      // pF  spectral   W-shape  ccothshape pilessdelta bindE
      // cccoh    nc1piC     ncoth    antinu/nu   nue/numu
      // RHC/FHC
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuMuNC1piC
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_mu
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
      #endif
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuMuBarNC1piC
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1,  0 }, // kSKmNuECCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,  0 }, // kSKmNuECC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1,  0 }, // kSKmNuECCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1,  0 }, // kSKmNuECCoth
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuENC1piC
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_mu

      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1,  0 }, // kSKmOscNuECCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCoth

      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCoth
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuEBarNC1piC
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCoth
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmMEC
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }  // kSKmUndefined
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
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNCcoh 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_e
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC0pi0 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_e
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCoth
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENC0pi0 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENCcoh 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENCoth_e

      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarC0pi0 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCcoh 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_e
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCoth
      #endif
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1Re
      ,
      // > 1 muon-like ring RHC
      {
      // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
      // banffcc1piE0 banffcc1piE1
      // pF  spectral   W-shape  ccothshape pilessdelta bindE
      // cccoh    nc1piC     ncoth    antinu/nu   nue/numu
      // RHC/FHC
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuMuNC1piC
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_mu
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
      #endif
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuMuBarNC1piC
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1,  0 }, // kSKmNuECCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,  0 }, // kSKmNuECC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1,  0 }, // kSKmNuECCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1,  0 }, // kSKmNuECCoth
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuENC1piC
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_mu

      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1,  0 }, // kSKmOscNuECCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCoth

      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCoth
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuEBarNC1piC
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCoth
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmMEC
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }  // kSKmUndefined
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
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNCcoh 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_e
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC0pi0 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_e
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCoth
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENC0pi0 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENCcoh 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENCoth_e

      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarC0pi0 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCcoh 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_e
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCoth
      #endif
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }  // kSKmUndefined
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
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuMuNC1piC
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_mu
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
      #endif
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuMuBarNC1piC
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1,  0 }, // kSKmNuECCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,  0 }, // kSKmNuECC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1,  0 }, // kSKmNuECCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1,  0 }, // kSKmNuECCoth
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuENC1piC
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_mu

      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1,  0 }, // kSKmOscNuECCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCoth

      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCoth
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuEBarNC1piC
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  1, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 1, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 1, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCoth
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmMEC
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }  // kSKmUndefined
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
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNCcoh 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_e
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC0pi0 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_e
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCoth
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENC0pi0 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENCcoh 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENCoth_e

      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarC0pi0 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCcoh 
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_e
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCC0pi
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCoth
      #endif
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//MultiRe

      */




      // Todo: remove old systematics


      // > 1 muon-like ring 
      {
	// MEC CA5 BgRES
	// banffmaqe banffmares 
	// banffscca banffsccv
	// banffccqeE0 banffccqeE1 banffccqeE2
	// banffcc1piE0 banffcc1piE1
	// pF  spectral   W-shape  ccothshape pilessdelta bindE
	// cccoh    nc1piC     ncoth    antinu/nu   nue/numu
	// RHC/FHC
	{  0, 0, 0, 1, 0, 0, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuMuNC1piC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 1, 0, 0, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuMuBarNC1piC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 1, 0, 1, 1, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1,  0 }, // kSKmNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,  0 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1,  0 }, // kSKmNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1,  0 }, // kSKmNuECCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuENC1piC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_mu

	{  0, 0, 0, 1, 0, 1, 1, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1,  0 }, // kSKmOscNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 1, 0, 1, 1, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuEBarNC1piC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 1, 1, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1Rmu
      ,
      // > 1 electron-like ring 
      {
	//if BANFF_BINS == 2 - NUMUONLY!
	// MEC CA5 BgRES
	// banffmaqe banffmares 
	// banffscca banffsccv
	// banffccqeE0 banffccqeE1 banffccqeE2
	// banffcc1piE0 banffcc1piE1 banffnc1pi0
	// pF  spectral   W-shape  ccothshape pilessdelta bindE
	// cccoh    nc1piC     ncoth    antinu/nu  nue/numu
	// RHC/FHC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC0pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENC0pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENCoth_e

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarC0pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1Re
      ,
      // > 1 muon-like ring RHC
      {
	// MEC CA5 BgRES
	// banffmaqe banffmares 
	// banffscca banffsccv
	// banffccqeE0 banffccqeE1 banffccqeE2
	// banffcc1piE0 banffcc1piE1
	// pF  spectral   W-shape  ccothshape pilessdelta bindE
	// cccoh    nc1piC     ncoth    antinu/nu   nue/numu
	// RHC/FHC
	{  0, 0, 0, 1, 0, 0, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuMuNC1piC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 1, 0, 0, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuMuBarNC1piC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 1, 0, 1, 1, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1,  0 }, // kSKmNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,  0 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1,  0 }, // kSKmNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1,  0 }, // kSKmNuECCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuENC1piC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_mu

	{  0, 0, 0, 1, 0, 1, 1, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1,  0 }, // kSKmOscNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 1, 0, 1, 1, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuEBarNC1piC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 1, 1, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1Rmu RHC
      ,
      // > 1 electron-like ring RHC
      {
	//if BANFF_BINS == 2 - NUMUONLY!
	// MEC CA5 BgRES
	// banffmaqe banffmares 
	// banffscca banffsccv
	// banffccqeE0 banffccqeE1 banffccqeE2
	// banffcc1piE0 banffcc1piE1 banffnc1pi0
	// pF  spectral   W-shape  ccothshape pilessdelta bindE
	// cccoh    nc1piC     ncoth    antinu/nu  nue/numu
	// RHC/FHC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC0pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENC0pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENCoth_e

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarC0pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1Re RHC
      ,

      // > Multi muon-like ring 
      {
	// MEC CA5 BgRES
	// banffmaqe banffmares 
	// banffscca banffsccv 
	// banffccqeE0 banffccqeE1 banffccqeE2
	// banffcc1piE0 banffcc1piE1
	// pF  spectral   W-shape  ccothshape pilessdelta bindE
	// cccoh    nc1piC     ncoth    antinu/nu   nue/numu
	// RHC/FHC
	{  0, 0, 0, 1, 0, 0, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuMuNC1piC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 1, 0, 0, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0,  0 }, // kSKmNuMuBarCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuMuBarNC1piC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0,  0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 1, 0, 1, 1, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1,  0 }, // kSKmNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,  0 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1,  0 }, // kSKmNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1,  0 }, // kSKmNuECCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuENC1piC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuENCoth_mu

	{  0, 0, 0, 1, 0, 1, 1, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1,  0 }, // kSKmOscNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1,  0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 1, 0, 1, 1, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1,  0 }, // kSKmNuEBarCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0,  0 }, // kSKmNuEBarNC1piC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 1, 1, 1, 1, 1,  0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 0,  1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1,  0 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//MultiRmu
      ,
      // > Multi electron-like ring 
      {
	//if BANFF_BINS == 2 - NUMUONLY!
	// MEC CA5 BgRES
	// banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
	// banffcc1piE0 banffcc1piE1 banffnc1pi0
	// pF  spectral   W-shape  ccothshape pilessdelta bindE
	// cccoh    nc1piC     ncoth    antinu/nu  nue/numu
	// RHC/FHC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNC0pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuECCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENC0pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuENCoth_e

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarC0pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCcoh 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // MEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,  0 }  // kSKmUndefined
      }//MultiRe



    };
#elif BANFF_BINS == 3 // NUE APP ONLY

    // Todo: add new systematics and remove old ones

    //Todo: make this work if kNumOfRHCXSecSyst is different from kNumOfXSecSyst.

    // NEW XSEC
    static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {


      // NEW XSEC

      /*
      // > 1 muon-like ring 
      {
      // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc1pi0
      // pF  spectral   W-shape  ccothshape pilessdelta
      // cccoh    nccoh   ncoth    nue/numu
      // RHC/FHC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuNC0piC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_mu
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNC0piC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCoth
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuENC0piC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuENCoth_mu

      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarNC0piC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCoth
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmMEC
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1Rmu
      ,
      // > 1 electron-like ring 
      {
      // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc0pi0
      // pF  spectral   W-shape  ccothshape pilessdelta 
      // cccoh   nccoh   ncoth    nue/numu 
      // RHC/FHC
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0 
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuMuNCcoh 
      {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuMuNCoth_e
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      #endif
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1pi0 
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
      {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuMuBarNCoth_e
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmNuECCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmNuECC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmNuECCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmNuECCoth
      {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuENC1pi0 
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuENCcoh 
      {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuENCoth_e

      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuECCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuECCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuECCoth

      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmNuEBarCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmNuEBarCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmNuEBarCCoth
      {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuEBarC1pi0 
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuEBarNCcoh 
      {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuEBarNCoth_e
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCCoth
      #endif
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1Re
      ,
      // > 1 muon-like ring RHC
      {
      // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc1pi0
      // pF  spectral   W-shape  ccothshape pilessdelta
      // cccoh    nccoh   ncoth    nue/numu
      // RHC/FHC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuNC0piC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_mu
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNC0piC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCoth
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuENC0piC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuENCoth_mu

      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarNC0piC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCoth
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmMEC
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1RmuRHC
      ,
      // > 1 electron-like ring RHC
      {
      // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc0pi0
      // pF  spectral   W-shape  ccothshape pilessdelta 
      // cccoh   nccoh   ncoth    nue/numu 
      // RHC/FHC
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0 
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuMuNCcoh 
      {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuMuNCoth_e
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      #endif
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1pi0 
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
      {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuMuBarNCoth_e
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmNuECCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmNuECC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmNuECCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmNuECCoth
      {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuENC1pi0 
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuENCcoh 
      {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuENCoth_e

      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuECCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuECCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuECCoth

      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmNuEBarCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmNuEBarCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmNuEBarCCoth
      {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuEBarC1pi0 
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuEBarNCcoh 
      {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuEBarNCoth_e
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCCoth
      #endif
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1ReRHC



      /// NEW MULTIRING
      ,
      // > Multi muon-like ring 
      {
      // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc1pi0
      // pF  spectral   W-shape  ccothshape pilessdelta
      // cccoh    nccoh   ncoth    nue/numu
      // RHC/FHC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuNC0piC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_mu
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNC0piC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCoth
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuENC0piC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuENCoth_mu

      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarNC0piC
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCoth
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmMEC
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }  // kSKmUndefined
      }//MultiRmu
      ,
      // > Multi electron-like ring 
      {
      // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc0pi0
      // pF  spectral   W-shape  ccothshape pilessdelta 
      // cccoh   nccoh   ncoth    nue/numu 
      // RHC/FHC
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0 
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuMuNCcoh 
      {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuMuNCoth_e
      #ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
      #endif
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
      {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1pi0 
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
      {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuMuBarNCoth_e
      #ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
      #endif
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmNuECCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmNuECC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmNuECCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmNuECCoth
      {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuENC1pi0 
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuENCcoh 
      {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuENCoth_e

      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuECCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuECCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuECCoth

      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmNuEBarCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmNuEBarCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmNuEBarCCoth
      {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuEBarC1pi0 
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuEBarNCcoh 
      {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuEBarNCoth_e
      #ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCCQE        
      {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCC1pi
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuEBarCCcoh
      {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCCoth
      #endif
      #ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      #ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      #endif
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }  // kSKmUndefined
      }//MultiRe



      };

      */



      // Todo: remove old systematics


      // > 1 muon-like ring 
      {
	// MEC CA5 BgRES
	// banffmaqe banffmares 
	// banffscca banffsccv
	// banffccqeE0 banffcc1piE0 banffnc1pi0
	// pF  spectral   W-shape  ccothshape pilessdelta
	// cccoh    nccoh   ncoth    nue/numu
	// RHC/FHC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuNC0piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNC0piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuENC0piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuENCoth_mu

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarNC0piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1Rmu
      ,
      // > 1 electron-like ring
      {
	// MEC CA5 BgRES
	// banffmaqe banffmares 
	// banffscca banffsccv
	// banffccqeE0 banffcc1piE0 banffnc0pi0
	// pF  spectral   W-shape  ccothshape pilessdelta 
	// cccoh   nccoh   ncoth    nue/numu 
	// RHC/FHC
	{  0, 0, 0, 1, 0, 0, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuMuNCcoh 
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
#endif
	{  0, 0, 0, 1, 0, 0, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmNuECCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuENC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuENCcoh 
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuENCoth_e

	{  0, 0, 0, 1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmNuEBarCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuEBarC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuEBarNCcoh 
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuEBarCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1Re
      ,
      // > 1 muon-like ring RHC
      {
	// banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc1pi0
	// pF  spectral   W-shape  ccothshape pilessdelta
	// cccoh    nccoh   ncoth    nue/numu
	// RHC/FHC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuNC0piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCC0pi   
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNC0piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuENC0piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuENCoth_mu

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarNC0piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1Rmu RHC
      ,
      // > 1 electron-like ring RHC
      {
	// MEC CA5 BgRES
	// banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc0pi0
	// pF  spectral   W-shape  ccothshape pilessdelta 
	// cccoh   nccoh   ncoth    nue/numu 
	// RHC/FHC
	{  0, 0, 0, 1, 0, 0, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuMuNCcoh 
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
#endif
	{  0, 0, 0, 1, 0, 0, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmNuECCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuENC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuENCcoh 
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuENCoth_e

	{  0, 0, 0, 1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmNuEBarCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuEBarC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuEBarNCcoh 
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuEBarCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }  // kSKmUndefined
      }//1Re RHC
      ,
      // > Multiple muon-like ring 
      {
	// MEC CA5 BgRES
	// banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc1pi0
	// pF  spectral   W-shape  ccothshape pilessdelta
	// cccoh    nccoh   ncoth    nue/numu
	// RHC/FHC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuNC0piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCC0pi  
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNC0piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuENC0piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuENCoth_mu

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCoth
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarNC0piC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCQE        
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECC0pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCcoh
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCMEC
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }  // kSKmUndefined
      }//MultiRmu
      ,
      // > Multiple electron-like ring 
      {
	// MEC CA5 BgRES
	// banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc0pi0
	// pF  spectral   W-shape  ccothshape pilessdelta 
	// cccoh   nccoh   ncoth    nue/numu 
	// RHC/FHC
	{  0, 0, 0, 1, 0, 0, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuMuNCcoh 
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuCCoth
#endif
	{  0, 0, 0, 1, 0, 0, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmNuMuBarCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuMuBarNC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuMuBarNCcoh 
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 0, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuMuBarCCoth
#endif
	{  0, 0, 0, 1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuECCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmNuECCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuENC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuENCcoh 
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuENCoth_e

	{  0, 0, 0, 1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuECCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmOscNuECC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuECCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmOscNuECCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuECCoth

	{  0, 0, 0, 1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmNuEBarCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmNuEBarCCoth
	{  0, 1, 1, 0, 1, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0,  0 }, // kSKmNuEBarC1pi0 
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0,  0 }, // kSKmNuEBarNCcoh 
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0,  0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCCQE        
	{  0, 1, 1, 0, 1, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCC1pi
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuEBarCCcoh
	{  1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1,  0 }, // kSKmOscNuEBarCCMEC
	{  0, 0, 0, 0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1,  0 }, // kSKmOscNuEBarCCoth
#endif
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
#endif
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }, // kSKmUndefined
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0,  0 }  // kSKmUndefined
      }//MultiRe

    };

#endif //BANFF_BINS

  } // conv
} // valor
  
#endif 

