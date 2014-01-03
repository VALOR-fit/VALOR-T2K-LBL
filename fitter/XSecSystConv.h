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
Added extra arrays for reverse horn current modes

*/
//____________________________________________________________________________

#ifndef xsec_syst_conv_h
#define xsec_syst_conv_h

#include "AnalysisOptions.h"
#include "Enumerations.h"
#include "SamplesAndModes.h"

namespace t2kana {
namespace conv   {

//
// ** Cross-section systematics
//

// Nominal values of cross-section related params 
static const double kNomMAQE  = 1.21;
static const double kNomMARES = 1.21;

// Number of cross-section related systematics 
#if BANFF_BINS == 1
static const int kNumOfXSecSyst = 20;
static const int kNumOfRHCXSecSyst = 20;
#elif BANFF_BINS == 2
static const int kNumOfXSecSyst = 18;
static const int kNumOfRHCXSecSyst = 18;
#elif BANFF_BINS == 3
static const int kNumOfXSecSyst = 14;
static const int kNumOfRHCXSecSyst = 14;
#elif BANFF_BINS == 4
static const int kNumOfXSecSyst = 20;
static const int kNumOfRHCXSecSyst = 20;
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

// The following is the list of the above kNumOfXSecSyst parameters
static const PhysParamId_t kListOfXSecSyst[kNumOfXSecSyst] = 
{
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
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
   kNuisBANFFNC1pi0,           // BANFF, NC1pi0 normalization, 0.0 - 30.0 GeV
#endif
#if BANFF_BINS == 4
   kNuisBANFFMECNorm,          // BANFF, MEC template normalisation
   kNuisBANFFPDDDoubleCount,   // BANFF, PDD double counting systematic
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
   kNuisXSecNormAntiNuToNu,     // antinu / nu cross section
#endif
   kNuisXSecNormNuEToNuMu     // nue / numu cross section

};

static const PhysParamId_t kListOfRHCXSecSyst[kNumOfRHCXSecSyst] = 
  {
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
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    kRHCNuisBANFFNC1pi0,           // BANFF, NC1pi0 normalization, 0.0 - 30.0 GeV
#endif
#if BANFF_BINS == 4
    kRHCNuisBANFFMECNorm,          // BANFF, MEC template normalisation
    kRHCNuisBANFFPDDDoubleCount,   // BANFF, PDD double counting systematic
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
    kRHCNuisXSecNormAntiNuToNu,     // antinu / nu cross section
#endif
    kRHCNuisXSecNormNuEToNuMu     // nue / numu cross section
  };


// The following list specifies which cross-section related systematic has its
// effect on MC templates pre-computed and stored in cubic splines.
// The order of systematic parameters should match the one in kListOfXSecSyst.
static const bool kXSecSystHasWghtSplines[kNumOfXSecSyst] = 
{
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
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
   false,                     // kNuisBANFFNC1pi0
#endif
#if BANFF_BINS == 4
   false,                     // kNuisBANFFMECNorm
   true,                      // kNuisBANFFPDDDoubleCount
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
   false                      // nue / numu cross section
};

static const bool kRHCXSecSystHasWghtSplines[kNumOfRHCXSecSyst] = 
  {
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
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    false,                     // kNuisBANFFNC1pi0
#endif
#if BANFF_BINS == 4
    false,                     // kNuisBANFFMECNorm
    true,                      // kNuisBANFFPDDDoubleCount
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
    false                      // nue / numu cross section
  };

// The following list specifies which cross-section related systematic 
// parameters were included in the BANFF fit.
// The order of systematic parameters should match the one in kListOfXSecSyst.
static const bool kXSecSystIncludedInBANFF[kNumOfXSecSyst] = 
{
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
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
   true,                      // kNuisBANFFNC1pi0
#endif
#if BANFF_BINS == 4
   true,                      // kNuisBANFFMECNorm
   true,                      // kNuisBANFFPDDDoubleCount
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
   false                      // nue / numu cross section
};

static const bool kRHCXSecSystIncludedInBANFF[kNumOfRHCXSecSyst] = 
  {
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
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    true,                      // kNuisBANFFNC1pi0
#endif
#if BANFF_BINS == 4
    true,                      // kNuisBANFFMECNorm
    true,                      // kNuisBANFFPDDDoubleCount
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
    false                      // nue / numu cross section
  };

// Ereco and Etrue cuts for the application of cross-section systematics.
// The order of systematic parameters should match the one in kListOfXSecSyst.
static const double kXSecSystErecoCut[kNumOfXSecSyst][2] =
{
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
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
   { 0.0, 100.0 },            // kNuisBANFFNC1pi0
#endif
#if BANFF_BINS == 4
   { 0.0, 100.0 },            // kNuisBANFFMECNorm
   { 0.0, 100.0 },            // kNuisBANFFPDDDoubleCount
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
   { 0.0, 100.0 }             // nue / numu cross section
};

static const double kRHCXSecSystErecoCut[kNumOfRHCXSecSyst][2] =
  {
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
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    { 0.0, 100.0 },            // kNuisBANFFNC1pi0
#endif
#if BANFF_BINS == 4
    { 0.0, 100.0 },            // kNuisBANFFMECNorm
    { 0.0, 100.0 },            // kNuisBANFFPDDDoubleCount
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
    { 0.0, 100.0 }             // nue / numu cross section
  };

static const double kXSecSystEtrueCut[kNumOfXSecSyst][2] =
{
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
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
   { 0.0, 100.0 },            // kNuisBANFFNC1pi0
#endif
#if BANFF_BINS == 4
   { 0.0, 100.0 },            // kNuisBANFFMECNorm
   { 0.0, 100.0 },            // kNuisBANFFPDDDoubleCount
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
   { 0.0, 100.0 }             // nue / numu cross section
};

static const double kRHCXSecSystEtrueCut[kNumOfRHCXSecSyst][2] =
  {
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
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    { 0.0, 100.0 },            // kNuisBANFFNC1pi0
#endif
#if BANFF_BINS == 4
    { 0.0, 100.0 },            // kNuisBANFFMECNorm
    { 0.0, 100.0 },            // kNuisBANFFPDDDoubleCount
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
    { 0.0, 100.0 }             // nue / numu cross section
  };


//Fractional error of systematic parameters, as of BANFF 2012v8
// ONLY USE FOR: macros/misc/test_fitter_inputs/test_valor_kdi_splines.C
static const double kXSecSystFracErr[kNumOfXSecSyst] =
{
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
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
  0.3,            // BANFF, NC1pi0 normalization, 0.0 - 30.0 GeV 
#endif
#if BANFF_BINS == 4
  -9999,          // BANFF, MEC template normalisation
  0.2,          // BANFF, PDD double counting systematic
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
  0.4,             // antinu / nu cross section
#endif
  0.03           // nue / numu cross section
};

static const double kRHCXSecSystFracErr[kNumOfRHCXSecSyst] =
  {
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
#if (BANFF_BINS != 2) && (BANFF_BINS != 4)
    0.3,            // BANFF, NC1pi0 normalization, 0.0 - 30.0 GeV 
#endif
#if BANFF_BINS == 4
    -9999,          // BANFF, MEC template normalisation
    0.2,          // BANFF, PDD double counting systematic
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
    0.4,             // antinu / nu cross section
#endif
    0.03           // nue / numu cross section
  };

// The following list describles in which MC templates we are supposed to apply the
// weights as a result of tweaking any of the above systematic parameters.
// The order of systematic parameters should match the one in kListOfXSecSyst.
#if BANFF_BINS == 1
#ifndef HKHACK

//Todo: make this work if kNumOfRHCXSecSyst is different from kNumOfXSecSyst.

static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {

#ifdef __FHC_RUNNING__

  // > 1 muon-like ring 
  {
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_mu

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_mu

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCCoth_mu
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth_mu

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_mu

#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCoth_mu
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmSigNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmSigNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmSigNuECCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmMEC

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring 
  {
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOScNuMuCCoth_e
#endif
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCCoth_e
#endif
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth_e

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuEBarC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCoth_e
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmSigNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmSigNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmSigNuECCoth_e
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined

  }

#else

  // > 1 muon-like ring 
  {
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_mu

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_mu

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCCoth_mu
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth_mu

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_mu

#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCoth_mu
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmSigNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmSigNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmSigNuECCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmMEC

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring 
  {

    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOScNuMuCCoth_e
#endif
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCCoth_e
#endif
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth_e

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuEBarC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCoth_e
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmSigNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmSigNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmSigNuECCoth_e
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined

  }

#endif
  ,

#ifdef __RHC_RUNNING__

  // > 1 muon-like ring with reverse horn current 
  {
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_mu

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_mu

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCCoth_mu
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth_mu

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_mu

#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCoth_mu
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmMEC

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring with reverse horn current 
  {

    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOScNuMuCCoth_e
#endif
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCCoth_e
#endif
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth_e

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuEBarC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCoth_e
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCoth_e
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined

  }

#else

  // > 1 muon-like ring with reverse horn current 
  {
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_mu

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_mu

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuOscMuBarCCoth_mu
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth_mu

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_mu

#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCoth_mu
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmMEC

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring with reverse horn current 
  {

    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC   nccoh   ncoth    antinu/nu   nue/numu
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmOScNuMuCCoth_e
#endif
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 0 }, // kSKmOscNuMuBarCCoth_e
#endif
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmNuECCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth_e

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_e
    {  0, 1, 0, 0, 0, 0, 0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 0, 0, 0, 0 }, // kSKmNuEBarC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCoth_e
#endif

    {  1, 0, 1, 1, 1, 0, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1, 0,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCoth_e
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined

  }

#endif
  
};
#else //HKHACK
///////////////////////////
///////HACK
///////////////////////////

//Todo: make this work if kNumOfRHCXSecSyst is different from kNumOfXSecSyst.

static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {

#ifdef __FHC_RUNNING__

  // > 1 muon-like ring 
  {
    // first sections all 0
    // numusignal numubackground nuesignal nuebackground, aka:
    // cccoh    nc1piC   nccoh   ncoth    0 0
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuOscMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCoth_mu
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu

#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth_mu
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmSigNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmSigNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmSigNuECCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmMEC

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring 
  {
    // first sections all 0
    // numusignal numubackground nuesignal nuebackground, aka:
    // cccoh    nc1piC   nccoh   ncoth    0 0
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOScNuMuCCoth_e
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCoth_e
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmOscNuEBarCCoth_e
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmSigNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmSigNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmSigNuECCoth_e
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined

  }

#else

  // > 1 muon-like ring 
  {
    // first sections all 0
    // numusignal numubackground nuesignal nuebackground, aka:
    // cccoh    nc1piC   nccoh   ncoth    0 0
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuOscMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCoth_mu
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu

#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuEBarCCoth_mu
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmSigNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmSigNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmSigNuECCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmMEC

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring 
  {
    // first sections all 0
    // numusignal numubackground nuesignal nuebackground, aka:
    // cccoh    nc1piC   nccoh   ncoth    0 0
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOScNuMuCCoth_e
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCoth_e
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmOscNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmOscNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmOscNuEBarCCoth_e
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmSigNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmSigNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmSigNuECCoth_e
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined

  }

#endif
  ,

#ifdef __RHC_RUNNING__

  // > 1 muon-like ring with reverse horn current 
  {
    // first sections all 0
    // numusignal numubackground nuesignal nuebackground, aka:
    // cccoh    nc1piC   nccoh   ncoth    0 0
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuOscMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCoth_mu
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu

#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCOscNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCOscNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCOscNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCOscNuECCoth_mu
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmMEC

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring with reverse horn current 
  {
    // first sections all 0
    // numusignal numubackground nuesignal nuebackground, aka:
    // cccoh    nc1piC   nccoh   ncoth    0 0
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOScNuMuCCoth_e
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCoth_e
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCOscNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCOscNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCOscNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCOscNuECCoth_e
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCSigNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCSigNuEBarCCoth_e
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined

  }

#else

  // > 1 muon-like ring with reverse horn current 
  {
    // first sections all 0
    // numusignal numubackground nuesignal nuebackground, aka:
    // cccoh    nc1piC   nccoh   ncoth    0 0
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuNCoth_mu

#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu

#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuOscMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuOscMuBarCCoth_mu
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuECCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarCCoth_mu
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu

#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCOscNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCOscNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCOscNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCOscNuECCoth_mu
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCoth_mu

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0 }, // kSKmMEC

#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring with reverse horn current 
  {
    // first sections all 0
    // numusignal numubackground nuesignal nuebackground, aka:
    // cccoh    nc1piC   nccoh   ncoth    0 0
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOScNuMuCCoth_e
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmOscNuMuBarCCoth_e
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarC1pi0 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCOscNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCOscNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCOscNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCOscNuECCoth_e
#endif

    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCSigNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 1, 0, 0, 0 }, // kSKmRHCSigNuEBarCCoth_e
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0 }  // kSKmUndefined

  }

#endif
  
};
//////////////////
///END HACK
//////////////////
#endif //HKHACK
#elif BANFF_BINS == 2

//Todo: make this work if kNumOfRHCXSecSyst is different from kNumOfXSecSyst.

static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {

#ifdef __FHC_RUNNING__

  // > 1 muon-like ring 
  {
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC     ncoth    antinu/nu   nue/numu
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmNuECCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuENCoth_mu

    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmSigNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmSigNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmSigNuECCoth_mu

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring 
  {
    //if BANFF_BINS == 2 - NUMUONLY!
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC     ncoth    antinu/nu  nue/numu
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 } // kSKmUndefined
  }

#else

  // > 1 muon-like ring 
  {
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC     ncoth    antinu/nu   nue/numu
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmNuECCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuENCoth_mu

    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmSigNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmSigNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmSigNuECCoth_mu

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring 
  {
    //if BANFF_BINS == 2 - NUMUONLY!
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC     ncoth    antinu/nu  nue/numu
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 } // kSKmUndefined
  }  

#endif
  ,
#ifdef __RHC_RUNNING__

  // > 1 muon-like ring with reverse horn current 
  {
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC     ncoth    antinu/nu   nue/numu
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmNuECCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuENCoth_mu

    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCoth_mu

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring with reverse horn current
  {
    //if BANFF_BINS == 2 - NUMUONLY!
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC     ncoth    antinu/nu  nue/numu
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 } // kSKmUndefined
  }

#else

    // > 1 muon-like ring with reverse horn current 
    {
      // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
      // banffcc1piE0 banffcc1piE1
      // pF  spectral   W-shape  ccothshape pilessdelta bindE
      // cccoh    nc1piC     ncoth    antinu/nu   nue/numu
      {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
      {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuNC1piC
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
      {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
      {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
#endif
      {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
      {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_mu
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuBarNC1piC
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCQE        
      {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmNuMuBarCCcoh
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmNuMuBarCCoth_mu
#endif
      {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmNuECCQE        
      {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmNuECC1pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmNuECCcoh
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmNuECCoth_mu
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuENC1piC
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuENCoth_mu

      {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
      {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuEBarNC1piC
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
      {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
#endif
      {  1, 0, 1, 1, 1, 0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCQE        
      {  0, 1, 0, 0, 0, 1, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCC1pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCcoh
      {  0, 1, 0, 0, 0, 0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmRHCSigNuEBarCCoth_mu

      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }  // kSKmUndefined
    }
    ,
    // > 1 electron-like ring with reverse horn current 
    {
      //if BANFF_BINS == 2 - NUMUONLY!
      // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
      // banffcc1piE0 banffcc1piE1 banffnc1pi0
      // pF  spectral   W-shape  ccothshape pilessdelta bindE
      // cccoh    nc1piC     ncoth    antinu/nu  nue/numu
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNC1pi0 
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e

      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e

      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_e
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNC0pi0 
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e

      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCQE        
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECC0pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCcoh
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCoth_e
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENC0pi0 
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCC0pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth_e
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarC0pi0 
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e

      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCQE        
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCC0pi
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCcoh
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCoth_e

      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
      {  0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 } // kSKmUndefined
    }  

#endif
  
};
#elif BANFF_BINS == 3

//Todo: make this work if kNumOfRHCXSecSyst is different from kNumOfXSecSyst.

static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {

#ifdef __FHC_RUNNING__

  // > 1 muon-like ring 
  {
    // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta
    // cccoh    nccoh   ncoth    nue/numu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarNC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmSigNuECCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmSigNuECC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmSigNuECCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmSigNuECCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmMEC

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 0 electron-like ring 
  {
    // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc0pi0
    // pF  spectral   W-shape  ccothshape pilessdelta 
    // cccoh   nccoh   ncoth    nue/numu 
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmNuMuCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmNuMuCCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuMuNCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmNuMuCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmNuMuCCoth_e
#endif
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmNuMuBarCCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuMuBarNCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmOscNuMuBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmOscNuMuBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmOscNuMuBarCCoth_e
#endif
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmNuECCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmNuECCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuENCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuENCoth_e

    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmNuEBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmNuEBarCCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuEBarC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuEBarNCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmOscNuEBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmOscNuEBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmOscNuEBarCCoth_e
#endif
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmSigNuECCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmSigNuECCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmSigNuECCoth_e
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined

  }

#else

  // > 1 muon-like ring 
  {
    // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta
    // cccoh    nccoh   ncoth    nue/numu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarNC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmSigNuECCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmSigNuECC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmSigNuECCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmSigNuECCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmMEC

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 0 electron-like ring 
  {

    // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc0pi0
    // pF  spectral   W-shape  ccothshape pilessdelta 
    // cccoh   nccoh   ncoth    nue/numu 
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmNuMuCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmNuMuCCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuMuNCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmNuMuCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmNuMuCCoth_e
#endif
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmNuMuBarCCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuMuBarNCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmOscNuMuBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmOscNuMuBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmOscNuMuBarCCoth_e
#endif
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmNuECCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmNuECCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuENCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuENCoth_e

    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmNuEBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmNuEBarCCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuEBarC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuEBarNCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmOscNuEBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmOscNuEBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmOscNuEBarCCoth_e
#endif
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmSigNuECCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmSigNuECCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmSigNuECCoth_e
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined

  }

#endif
  ,
#ifdef __RHC_RUNNING__

  // > 1 muon-like ring with reverse horn current 
  {
    // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta
    // cccoh    nccoh   ncoth    nue/numu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarNC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmRHCSigNuEBarCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmMEC

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 0 electron-like ring with reverse horn current 
  {
    // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc0pi0
    // pF  spectral   W-shape  ccothshape pilessdelta 
    // cccoh   nccoh   ncoth    nue/numu 
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmNuMuCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmNuMuCCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuMuNCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmNuMuCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmNuMuCCoth_e
#endif
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmNuMuBarCCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuMuBarNCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmOscNuMuBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmOscNuMuBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmOscNuMuBarCCoth_e
#endif
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmNuECCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmNuECCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuENCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuENCoth_e

    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmNuEBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmNuEBarCCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuEBarC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuEBarNCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmRHCOscNuECCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmRHCOscNuECC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmRHCOscNuECCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmRHCOscNuECCoth_e
#endif
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmRHCSigNuEBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmRHCSigNuEBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmRHCSigNuEBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmRHCSigNuEBarCCoth_e
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined

  }

#else

  // > 1 muon-like ring with reverse horn current
  {
    // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta
    // cccoh    nccoh   ncoth    nue/numu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuNCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarCCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuMuBarNCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuECCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuENCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarCCoth_mu
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarNC0piC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmNuEBarNCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmRHCSigNuEBarCC0pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCoth_mu

    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmMEC
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 0 electron-like ring with reverse horn current 
  {

    // banffmaqe banffmares banffccqeE0 banffcc1piE0 banffnc0pi0
    // pF  spectral   W-shape  ccothshape pilessdelta 
    // cccoh   nccoh   ncoth    nue/numu 
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmNuMuCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmNuMuCCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuMuNC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuMuNCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuMuNCoth_e
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmNuMuCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmNuMuCCoth_e
#endif
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmNuMuBarCCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuMuBarNC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuMuBarNCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuMuBarNCoth_e
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 0}, // kSKmOscNuMuBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 0}, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 0}, // kSKmOscNuMuBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 0}, // kSKmOscNuMuBarCCoth_e
#endif
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmNuECCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmNuECCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuENC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuENCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuENCoth_e

    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmNuEBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmNuEBarCCoth_e
    {  0, 1, 0, 0, 1,   0, 0, 1, 0, 1,   0, 0, 0, 0}, // kSKmNuEBarC1pi0 
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 1, 0, 0}, // kSKmNuEBarNCcoh 
    {  0, 1, 0, 0, 0,   0, 0, 1, 0, 1,   0, 0, 1, 0}, // kSKmNuEBarNCoth_e
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmRHCOscNuECCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmRHCOscNuECC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmRHCOscNuECCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmRHCOscNuECCoth_e
#endif
    {  1, 0, 1, 0, 0,   1, 1, 0, 0, 0,   0, 0, 0, 1}, // kSKmRHCSigNuEBarCCQE        
    {  0, 1, 0, 1, 0,   0, 0, 0, 0, 1,   0, 0, 0, 1}, // kSKmRHCSigNuEBarCC1pi
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   1, 0, 0, 1}, // kSKmRHCSigNuEBarCCcoh
    {  0, 1, 0, 0, 0,   0, 0, 0, 1, 0,   0, 0, 0, 1}, // kSKmRHCSigNuEBarCCoth_e
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined
    {  0, 0, 0, 0, 0,   0, 0, 0, 0, 0,   0, 0, 0, 0}, // kSKmUndefined

  }

#endif
  
};
#elif BANFF_BINS == 4

//Todo: make this work if kNumOfRHCXSecSyst is different from kNumOfXSecSyst.

static const bool kXsApplyXSecSystToSKMode[kns][knm][kNumOfXSecSyst] = {

#ifdef __FHC_RUNNING__

  // > 1 muon-like ring 
  {
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1
    // banff MEC, banff PDDDC
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC    ncoth    nue/numu antinu/nu
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
 #endif
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmNuMuBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmOscNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmOscNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmOscNuMuBarCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmNuECCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuENCoth_mu

    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmOscNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmSigNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmSigNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmSigNuECCoth_mu

    //note that BANFFMECNorm is not handled by NuXSec (instead by SuperkMC) so it is set to false here
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring 
  {

    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc0piC   nccoh   ncoth    nue/numu antinu/nu
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }  // kSKmUndefined

  }

#else

  // > 1 muon-like ring 
  {
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1
    // banff MEC, banff PDDDC
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC    ncoth    nue/numu antinu/nu
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
 #endif
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmNuMuBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmOscNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmOscNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmOscNuMuBarCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmNuECCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuENCoth_mu

    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmOscNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmOscNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmOscNuEBarCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmSigNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmSigNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmSigNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmSigNuECCoth_mu

    //note that BANFFMECNorm is not handled by NuXSec (instead by SuperkMC) so it is set to false here
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring 
  {
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc0piC   nccoh   ncoth    nue/numu antinu/nu
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmSigNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }  // kSKmUndefined

  }

#endif
  ,
#ifdef __RHC_RUNNING__

  // > 1 muon-like ring with reverse horn current
  {
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1
    // banff MEC, banff PDDDC
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC    ncoth    nue/numu antinu/nu
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
 #endif
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmNuMuBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmOscNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmOscNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmOscNuMuBarCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmNuECCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuENCoth_mu

    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmRHCOscNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmRHCOscNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmRHCSigNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmRHCSigNuEBarCCoth_mu

    //note that BANFFMECNorm is not handled by NuXSec (instead by SuperkMC) so it is set to false here
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring with reverse horn current
  {

    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc0piC   nccoh   ncoth    nue/numu antinu/nu
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }  // kSKmUndefined

  }

#else

  // > 1 muon-like ring with reverse horn current 
  {
    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1
    // banff MEC, banff PDDDC
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc1piC    ncoth    nue/numu antinu/nu
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuNCoth_mu
#ifdef __USE_OSC_NUMU_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_mu
 #endif
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmNuMuBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuMuBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuMuBarNCoth_mu
#ifdef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 0, 1 }, // kSKmOscNuMuBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 0, 1 }, // kSKmOscNuMuBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 0, 1 }, // kSKmOscNuMuBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 0, 1 }, // kSKmOscNuMuBarCCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmNuECCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuENC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuENCoth_mu

    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmNuEBarCCoth_mu
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 1, 0, 0, 0 }, // kSKmNuEBarNC1piC
    {  0, 1, 0, 0, 0, 0, 0,   0, 1,   0, 0, 1, 0, 1, 0,   0, 0, 1, 0, 0 }, // kSKmNuEBarNCoth_mu
#ifdef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 1 }, // kSKmRHCOscNuECC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 1 }, // kSKmRHCOscNuECCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 1 }, // kSKmRHCOscNuECCoth_mu
#endif
    {  1, 0, 1, 1, 1, 0, 0,   0, 0,   1, 1, 0, 0, 0, 1,   0, 0, 0, 1, 0 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 1, 0, 0, 0, 1, 1,   0, 1,   0, 0, 0, 0, 1, 0,   0, 0, 0, 1, 0 }, // kSKmRHCSigNuEBarCC1pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   1, 0, 0, 1, 0 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 1, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 1, 0, 0,   0, 0, 0, 1, 0 }, // kSKmRHCSigNuEBarCCoth_mu

    //note that BANFFMECNorm is not handled by NuXSec (instead by SuperkMC) so it is set to false here
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmMEC
#ifndef __USE_OSC_NUMU_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUMUBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
#ifndef __USE_OSC_NUEBAR_MC_TEMPLATES__
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
#endif
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }  // kSKmUndefined
  }
  ,
  // > 1 electron-like ring with reverse horn current
  {

    // banffmaqe banffmares banffccqeE0 banffccqeE1 banffccqeE2
    // banffcc1piE0 banffcc1piE1 banffnc1pi0
    // pF  spectral   W-shape  ccothshape pilessdelta bindE
    // cccoh    nc0piC   nccoh   ncoth    nue/numu antinu/nu
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmOscNuMuCCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuMuBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuECCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuENCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarCCoth_e
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarC0pi0 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCcoh 
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmNuEBarNCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCQE        
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCC0pi
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCcoh
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmRHCSigNuEBarCCoth_e

    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }, // kSKmUndefined
    {  0, 0, 0, 0, 0, 0, 0,   0, 0,   0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0 }  // kSKmUndefined

  }

#endif
  
};
#endif

} // conv
} // t2kana
  
#endif 

