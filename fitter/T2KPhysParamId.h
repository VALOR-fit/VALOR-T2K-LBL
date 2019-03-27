//____________________________________________________________________________
/*!

\class     valor::conv::PhysParamId_t;

\brief     Conventions

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory
          
\created   Oct 20, 2010

\copyright (c) 2010-2014, The VALOR group

*/
//____________________________________________________________________________

#ifndef _t2k_phys_param_id_h_
#define _t2k_phys_param_id_h_

#include <string>

#include "AnalysisOptions.h"

#include "valor/Conventions/PhysParamId.h"

using std::string;

namespace valor {

namespace conv {
   
   
   /*
    * List of 2017 analysis nuisance parameters 
    * -----------------------------------------
   */

   //********* FORWARD HORN CURRENT ************************ 

   /*
    * SuperK detector errors (energy scale and efficiencies)
   */   
static const PhysParamId_t kNuisSuperkEnergyScale     = kNuis000;     // SK energy scale
#if FSI_BINS != 3
static const PhysParamId_t kNuisSuperkEff1RmuCCQEoth  = kNuis001; // SK 1Rmu CCQE efficiency (other)
static const PhysParamId_t kNuisSuperkEff1RmuCCQErc0u = kNuis002; // SK 1Rmu CCQE efficiency (ring counting; Ereco bin 1)
static const PhysParamId_t kNuisSuperkEff1RmuCCQErc1u = kNuis003; // SK 1Rmu CCQE efficiency (ring counting; Ereco bin 2)
static const PhysParamId_t kNuisSuperkEff1RmuCCQErc2u = kNuis004; // SK 1Rmu CCQE efficiency (ring counting; Ereco bin 3)
static const PhysParamId_t kNuisSuperkEff1RmuCCnonQE  = kNuis005; // SK 1Rmu CCnonQE efficiency
static const PhysParamId_t kNuisSuperkEff1RmuNC       = kNuis006; // SK 1Rmu NC efficiency
static const PhysParamId_t kNuisSuperkEff1RmuNueCC    = kNuis007; // SK 1Rmu nue CC efficiency
#endif

#if FSI_BINS != 2
static const PhysParamId_t kNuisSuperkEff1ReNueSigCC0 = kNuis001;    // SK 1Re CC eff for osc nue Ereco bin 0
static const PhysParamId_t kNuisSuperkEff1ReNueSigCC1 = kNuis002;    // SK 1Re CC eff for osc nue Ereco bin 1
static const PhysParamId_t kNuisSuperkEff1ReNueSigCC2 = kNuis003;    // SK 1Re CC eff for osc nue Ereco bin 2
static const PhysParamId_t kNuisSuperkEff1ReNumuCC0   = kNuis004;    // SK 1Re CC eff for numu Ereco bin 0
static const PhysParamId_t kNuisSuperkEff1ReNumuCC1   = kNuis005;    // SK 1Re CC eff for numu Ereco bin 1
static const PhysParamId_t kNuisSuperkEff1ReNumuCC2   = kNuis006;    // SK 1Re CC eff for numu Ereco bin 2
static const PhysParamId_t kNuisSuperkEff1ReNueCC0    = kNuis007;    // SK 1Re CC eff for nue Ereco bin 0
static const PhysParamId_t kNuisSuperkEff1ReNueCC1    = kNuis008;    // SK 1Re CC eff for nue Ereco bin 1
static const PhysParamId_t kNuisSuperkEff1ReNueCC2    = kNuis009;    // SK 1Re CC eff for nue Ereco bin 2
static const PhysParamId_t kNuisSuperkEff1ReNC0       = kNuis010;    // SK 1Re NC eff Ereco bin 0
static const PhysParamId_t kNuisSuperkEff1ReNC1       = kNuis011;    // SK 1Re NC eff Ereco bin 1
static const PhysParamId_t kNuisSuperkEff1ReNC2       = kNuis012;    // SK 1Re NC eff Ereco bin 2
#endif

   /*
    * (flux and cross-section) correlated parameters used in the BANFF fit
   */   
static const PhysParamId_t kNuisBANFF00 = kNuis020; // # 00; SK     numu    flux;  E =  0.0 -  0.4 GeV
static const PhysParamId_t kNuisBANFF01 = kNuis021; // # 01; SK     numu    flux;  E =  0.4 -  0.5 GeV
static const PhysParamId_t kNuisBANFF02 = kNuis022; // # 02; SK     numu    flux;  E =  0.5 -  0.6 GeV
static const PhysParamId_t kNuisBANFF03 = kNuis023; // # 03; SK     numu    flux;  E =  0.6 -  0.7 GeV
static const PhysParamId_t kNuisBANFF04 = kNuis024; // # 04; SK     numu    flux;  E =  0.7 -  1.0 GeV
static const PhysParamId_t kNuisBANFF05 = kNuis025; // # 05; SK     numu    flux;  E =  1.0 -  1.5 GeV
static const PhysParamId_t kNuisBANFF06 = kNuis026; // # 06; SK     numu    flux;  E =  1.5 -  2.5 GeV
static const PhysParamId_t kNuisBANFF07 = kNuis027; // # 07; SK     numu    flux;  E =  2.5 -  3.5 GeV
static const PhysParamId_t kNuisBANFF08 = kNuis028; // # 08; SK     numu    flux;  E =  3.5 -  5.0 GeV
static const PhysParamId_t kNuisBANFF09 = kNuis029; // # 09; SK     numu    flux;  E =  5.0 -  7.0 GeV
static const PhysParamId_t kNuisBANFF10 = kNuis030; // # 10; SK     numu    flux;  E =  7.0 - 30.0 GeV
static const PhysParamId_t kNuisBANFF11 = kNuis031; // # 11; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF12 = kNuis032; // # 12; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF13 = kNuis033; // # 13; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF14 = kNuis034; // # 14; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF15 = kNuis035; // # 15; different meanings depending on samples defined
#if !(BANFF_BINS == 2)
static const PhysParamId_t kNuisBANFF16 = kNuis036; // # 16; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF17 = kNuis037; // # 17; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF18 = kNuis038; // # 18; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF19 = kNuis039; // # 19; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF20 = kNuis040; // # 20; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF21 = kNuis041; // # 21; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF22 = kNuis042; // # 22; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF23 = kNuis043; // # 23; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF24 = kNuis044; // # 24; different meanings depending on samples defined
#endif

static const PhysParamId_t kNuisBANFFMAQE    = kNuis045;  // # M_{A}  QE
static const PhysParamId_t kNuisBANFFMARES   = kNuis046;  // # M_{A}  RES
static const PhysParamId_t kNuisBANFFCCQEE0  = kNuis047;  // # CCQE   normalization; 0 - 1.5 GeV
#if BANFF_BINS != 3
static const PhysParamId_t kNuisBANFFCCQEE1  = kNuis048;  // # CCQE   normalization; 1.5 -  3.5 GeV
static const PhysParamId_t kNuisBANFFCCQEE2  = kNuis049;  // # CCQE   normalization; 3.5 - 30.0 GeV
#endif
static const PhysParamId_t kNuisBANFFCC1piE0 = kNuis050; // # CC1pi  normalization; 0 - 2.5 GeV
#if BANFF_BINS != 3
static const PhysParamId_t kNuisBANFFCC1piE1 = kNuis051; // # CC1pi  normalization; 2.5 - 30.0 GeV
#endif
#if (BANFF_BINS != 2)
static const PhysParamId_t kNuisBANFFNC1pi0  = kNuis052; // # NC1pi0 normalization; 0 - 30.0 GeV
#endif

static const PhysParamId_t kNuisBANFFNorm2p2h = kNuis053;  // 2p2h nu normalization
static const PhysParamId_t kNuisBANFFCA5      = kNuis054;  // CA5 RES 
static const PhysParamId_t kNuisBANFFBgRES    = kNuis055;  // Isospin 1/2 normalization
// NEW SPLIT ANA
#if (BANFF_BINS != 2)
static const PhysParamId_t kNuisBANFFSCCA     = kNuis056;  // Second current axial
static const PhysParamId_t kNuisBANFFSCCV     = kNuis057;  // Second current vector
#endif
static const PhysParamId_t kNuisBANFFFermiMom            = kNuis058;   // Fermi momentum
static const PhysParamId_t kNuisBANFFShapeCCoth          = kNuis059;   // CC other shape
  //static const PhysParamId_t kNuisBANFFBindE               = kNuis060;   // Binding energy, not in 2017 summer analysis
static const PhysParamId_t kNuisBANFFNormCCcoh           = kNuis060;   // CC coherent norm
static const PhysParamId_t kNuisBANFFNormNCcoh           = kNuis061;   // NC coh norm
static const PhysParamId_t kNuisBANFFNormNCoth           = kNuis062;   // NC other norm
static const PhysParamId_t kNuisBANFFNormNuEToNuMu       = kNuis063;   // nue / numu cross section
static const PhysParamId_t kNuisBANFFNormNC1Gamma        = kNuis064;   // NC 1 gamma
static const PhysParamId_t kNuisBANFFNormNuEBarToNuMuBar = kNuis065;   // nuebar / numubar 
static const PhysParamId_t kNuisBANFFNorm2p2hBar         = kNuis066;   // 2p2h nubar normalization
static const PhysParamId_t kNuisBANFFBeRPA_A             = kNuis067;   // BeRPA A
static const PhysParamId_t kNuisBANFFBeRPA_B             = kNuis068;   // BeRPA B
static const PhysParamId_t kNuisBANFFBeRPA_D             = kNuis069;   // BeRPA D
static const PhysParamId_t kNuisBANFFBeRPA_E             = kNuis070;   // BeRPA E
static const PhysParamId_t kNuisBANFFBeRPA_U             = kNuis071;   // BeRPA U
static const PhysParamId_t kNuisBANFFShape2p2h           = kNuis072;   // 2p2h nu/nubar shape
static const PhysParamId_t kNuisBANFFNorm2p2hCtoO        = kNuis073;   // 2p2h CtoO norm
static const PhysParamId_t kNuisBANFFEb                  = kNuis074;   // Eb
static const PhysParamId_t kNuisBANFFMultiReNorm         = kNuis075;   // Normalisation of all events in the FHC MultiRe sample
//

static const PhysParamId_t kLastBANFFParam  = kNuis076; // Utility, should always come after last FHC BANFF

   /*
    * Cross-section parameters (excluding cross-section params which were included in BANFF fit and are correlated with the flux params)
    * These other paramers are taken to be uncorrelated and are applied to SK only.
   */   
   //systematics which use splines

// NEW XSEC
static const PhysParamId_t kNuisXSecRPA            = kNuis080;   // RPA tuning

static const PhysParamId_t kNuisXSecFermiMom       = kNuis081;   // Fermi momentum
static const PhysParamId_t kNuisXSecSpecFunc       = kNuis082;   // spectral function on/off
static const PhysParamId_t kNuisXSecWshape         = kNuis083;   // W-shape
static const PhysParamId_t kNuisXSecShapeCCoth     = kNuis084;   // CC other shape
static const PhysParamId_t kNuisXSecPiLessDelta    = kNuis085;   // Pion-less delta decay
static const PhysParamId_t kNuisXSecBindE          = kNuis086;   // Binding energy
   //systematics that are normalisations
static const PhysParamId_t kNuisXSecNormCCcoh      = kNuis087;   // CC coherent norm
static const PhysParamId_t kNuisXSecNormNC1piC     = kNuis088;   // NC 1piC norm
static const PhysParamId_t kNuisXSecNormNCcoh      = kNuis089;   // NC coh norm
static const PhysParamId_t kNuisXSecNormNCoth      = kNuis090;   // NC other norm
static const PhysParamId_t kNuisXSecNormNuEToNuMu  = kNuis091;   // nue / numu cross section
static const PhysParamId_t kNuisXSecNormAntiNuToNu = kNuis092;   // antinu / nu cross section

static const PhysParamId_t kLastXSecParam = kNuis093;   // Utility, should always come after last FHC Xsec param

   /*
    * FSI parameters.
    * Assuming here a 10x10 covariance matrix will be made available.
    * Each parameter is a scaling factor to apply at some Ereco range (to be defined) for some
    * SK reaction mode (or the total? to be defined)
    */

// NEW SKDETFSI

#if defined(__USE_SKDETFSI2013__)
 // SKDET+FSI matrix 18x18
 #if FSI_BINS != 3
 static const PhysParamId_t kNuisSKDetFSI000 = kNuis101;
 static const PhysParamId_t kNuisSKDetFSI001 = kNuis102;
 static const PhysParamId_t kNuisSKDetFSI002 = kNuis103;
 static const PhysParamId_t kNuisSKDetFSI003 = kNuis104;
 static const PhysParamId_t kNuisSKDetFSI004 = kNuis105;
 static const PhysParamId_t kNuisSKDetFSI005 = kNuis106;
 #endif
 #if FSI_BINS != 2
 static const PhysParamId_t kNuisSKDetFSI006 = kNuis107;
 static const PhysParamId_t kNuisSKDetFSI007 = kNuis108;
 static const PhysParamId_t kNuisSKDetFSI008 = kNuis109;
 static const PhysParamId_t kNuisSKDetFSI009 = kNuis110;
 static const PhysParamId_t kNuisSKDetFSI010 = kNuis111;
 static const PhysParamId_t kNuisSKDetFSI011 = kNuis112;
 static const PhysParamId_t kNuisSKDetFSI012 = kNuis113;
 static const PhysParamId_t kNuisSKDetFSI013 = kNuis114;
 static const PhysParamId_t kNuisSKDetFSI014 = kNuis115;
 static const PhysParamId_t kNuisSKDetFSI015 = kNuis116;
 static const PhysParamId_t kNuisSKDetFSI016 = kNuis117;
 static const PhysParamId_t kNuisSKDetFSI017 = kNuis118;
 #endif

 // NEW SKDETFSI
 #if FSI_BINS != 3
 static const PhysParamId_t kNuisSKDetFSI_first = kNuisSKDetFSI000; 
 #else
 static const PhysParamId_t kNuisSKDetFSI_first = kNuisSKDetFSI006; 
 #endif

#endif



   //********* REVERSE HORN CURRENT **************************

   /*
    * Special RHC to FHC errors
    */

static const PhysParamId_t kNuisRHCToFHC = kNuis125;

   /*
    * SuperK detector errors (energy scale and efficiencies)
    */
   
   /*
    * SuperK detector errors (energy scale and efficiencies)
   */   
static const PhysParamId_t kRHCNuisSuperkEnergyScale     = kNuis150;     // SK energy scale
#if FSI_BINS != 3
static const PhysParamId_t kRHCNuisSuperkEff1RmuCCQEoth  = kNuis151; // SK 1Rmu CCQE efficiency (other)
static const PhysParamId_t kRHCNuisSuperkEff1RmuCCQErc0u = kNuis152; // SK 1Rmu CCQE efficiency (ring counting; Ereco bin 1)
static const PhysParamId_t kRHCNuisSuperkEff1RmuCCQErc1u = kNuis153; // SK 1Rmu CCQE efficiency (ring counting; Ereco bin 2)
static const PhysParamId_t kRHCNuisSuperkEff1RmuCCQErc2u = kNuis154; // SK 1Rmu CCQE efficiency (ring counting; Ereco bin 3)
static const PhysParamId_t kRHCNuisSuperkEff1RmuCCnonQE  = kNuis155; // SK 1Rmu CCnonQE efficiency
static const PhysParamId_t kRHCNuisSuperkEff1RmuNC       = kNuis156; // SK 1Rmu NC efficiency
static const PhysParamId_t kRHCNuisSuperkEff1RmuNueCC    = kNuis157; // SK 1Rmu nue CC efficiency
#endif

#if FSI_BINS != 2
static const PhysParamId_t kRHCNuisSuperkEff1ReNueSigCC0 = kNuis151;    // SK 1Re CC eff for osc nue Ereco bin 0
static const PhysParamId_t kRHCNuisSuperkEff1ReNueSigCC1 = kNuis152;    // SK 1Re CC eff for osc nue Ereco bin 1
static const PhysParamId_t kRHCNuisSuperkEff1ReNueSigCC2 = kNuis153;    // SK 1Re CC eff for osc nue Ereco bin 2
static const PhysParamId_t kRHCNuisSuperkEff1ReNumuCC0   = kNuis154;    // SK 1Re CC eff for numu Ereco bin 0
static const PhysParamId_t kRHCNuisSuperkEff1ReNumuCC1   = kNuis155;    // SK 1Re CC eff for numu Ereco bin 1
static const PhysParamId_t kRHCNuisSuperkEff1ReNumuCC2   = kNuis156;    // SK 1Re CC eff for numu Ereco bin 2
static const PhysParamId_t kRHCNuisSuperkEff1ReNueCC0    = kNuis157;    // SK 1Re CC eff for nue Ereco bin 0
static const PhysParamId_t kRHCNuisSuperkEff1ReNueCC1    = kNuis158;    // SK 1Re CC eff for nue Ereco bin 1
static const PhysParamId_t kRHCNuisSuperkEff1ReNueCC2    = kNuis159;    // SK 1Re CC eff for nue Ereco bin 2
static const PhysParamId_t kRHCNuisSuperkEff1ReNC0       = kNuis160;    // SK 1Re NC eff Ereco bin 0
static const PhysParamId_t kRHCNuisSuperkEff1ReNC1       = kNuis161;    // SK 1Re NC eff Ereco bin 1
static const PhysParamId_t kRHCNuisSuperkEff1ReNC2       = kNuis162;    // SK 1Re NC eff Ereco bin 2
#endif

   /*
    * (flux and cross-section) correlated parameters used in the BANFF fit
   */   
static const PhysParamId_t kRHCNuisBANFF00 = kNuis170; // # 00; SK numu
static const PhysParamId_t kRHCNuisBANFF01 = kNuis171; // # 01; SK numu
static const PhysParamId_t kRHCNuisBANFF02 = kNuis172; // # 02; SK numu
static const PhysParamId_t kRHCNuisBANFF03 = kNuis173; // # 03; SK numu 
static const PhysParamId_t kRHCNuisBANFF04 = kNuis174; // # 04; SK numu
static const PhysParamId_t kRHCNuisBANFF05 = kNuis175; // # 05; SK numu 
static const PhysParamId_t kRHCNuisBANFF06 = kNuis176; // # 06; SK anti-numu
static const PhysParamId_t kRHCNuisBANFF07 = kNuis177; // # 07; SK anti-numu
static const PhysParamId_t kRHCNuisBANFF08 = kNuis178; // # 08; SK anti-numu 
static const PhysParamId_t kRHCNuisBANFF09 = kNuis179; // # 09; SK anti-numu 
static const PhysParamId_t kRHCNuisBANFF10 = kNuis180; // # 10; SK anti-numu 
static const PhysParamId_t kRHCNuisBANFF11 = kNuis181; // # 11; SK anti-numu (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF12 = kNuis182; // # 12; SK anti-numu (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF13 = kNuis183; // # 13; SK anti-numu (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF14 = kNuis184; // # 14; SK anti-numu (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF15 = kNuis185; // # 15; SK anti-numu (different meanings depending on samples defined)
#if !(BANFF_BINS == 2)
static const PhysParamId_t kRHCNuisBANFF16 = kNuis186; // # 16; SK nue (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF17 = kNuis187; // # 17; SK nue (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF18 = kNuis188; // # 18; SK anti-nue (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF19 = kNuis189; // # 19; SK anti-nu (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF20 = kNuis190; // # 20; SK anti-nue (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF21 = kNuis191; // # 21; SK anti-nue (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF22 = kNuis192; // # 22; SK anti-nue (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF23 = kNuis193; // # 23; SK anti-nue (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF24 = kNuis194; // # 24; SK anti-nue (different meanings depending on samples defined)
#endif

static const PhysParamId_t kRHCNuisBANFFMAQE    = kNuis195;  // # M_{A}  QE
static const PhysParamId_t kRHCNuisBANFFMARES   = kNuis196;  // # M_{A}  RES
static const PhysParamId_t kRHCNuisBANFFCCQEE0  = kNuis197;  // # CCQE   normalization; 0 - 1.5 GeV
#if BANFF_BINS != 3
static const PhysParamId_t kRHCNuisBANFFCCQEE1  = kNuis198;  // # CCQE   normalization; 1.5 -  3.5 GeV
static const PhysParamId_t kRHCNuisBANFFCCQEE2  = kNuis199;  // # CCQE   normalization; 3.5 - 30.0 GeV
#endif
static const PhysParamId_t kRHCNuisBANFFCC1piE0 = kNuis200; // # CC1pi  normalization; 0 - 2.5 GeV
#if BANFF_BINS != 3
static const PhysParamId_t kRHCNuisBANFFCC1piE1 = kNuis201; // # CC1pi  normalization; 2.5 - 30.0 GeV
#endif
#if (BANFF_BINS != 2)
static const PhysParamId_t kRHCNuisBANFFNC1pi0  = kNuis202; // # NC1pi0 normalization; 0 - 30.0 GeV
#endif

// NEW XSEC
static const PhysParamId_t kRHCNuisBANFF2p2h     = kNuis203;  // 2p2h normalization
static const PhysParamId_t kRHCNuisBANFFCA5     = kNuis204;  // CA5 RES 
static const PhysParamId_t kRHCNuisBANFFBgRES   = kNuis205;  // Isospin 1/2 normalization
// NEW SPLIT ANA
#if (BANFF_BINS != 2)
static const PhysParamId_t kRHCNuisBANFFSCCA    = kNuis206;  // Second current axial
static const PhysParamId_t kRHCNuisBANFFSCCV    = kNuis207;  // Second current vector
#endif
static const PhysParamId_t kRHCNuisBANFFFermiMom            = kNuis208;   // Fermi momentum
static const PhysParamId_t kRHCNuisBANFFShapeCCoth          = kNuis209;   // CC other shape
  //static const PhysParamId_t kRHCNuisBANFFBindE               = kNuis210;   // Binding energy
static const PhysParamId_t kRHCNuisBANFFNormCCcoh           = kNuis210;   // CC coherent norm
static const PhysParamId_t kRHCNuisBANFFNormNCcoh           = kNuis211;   // NC coh norm
static const PhysParamId_t kRHCNuisBANFFNormNCoth           = kNuis212;   // NC other norm
static const PhysParamId_t kRHCNuisBANFFNormNuEToNuMu       = kNuis213;   // nue / numu cross section
static const PhysParamId_t kRHCNuisBANFFNormNC1Gamma        = kNuis214;   // nue / numu cross section
static const PhysParamId_t kRHCNuisBANFFNormNuEBarToNuMuBar = kNuis215;   // nue / numu cross section
static const PhysParamId_t kRHCNuisBANFF2p2hBar             = kNuis216;  // 2p2h normalization
static const PhysParamId_t kRHCNuisBANFFBeRPA_A             = kNuis217;   // BeRPA A
static const PhysParamId_t kRHCNuisBANFFBeRPA_B             = kNuis218;   // BeRPA B
static const PhysParamId_t kRHCNuisBANFFBeRPA_D             = kNuis219;   // BeRPA D
static const PhysParamId_t kRHCNuisBANFFBeRPA_E             = kNuis220;   // BeRPA E
static const PhysParamId_t kRHCNuisBANFFBeRPA_U             = kNuis221;   // BeRPA U
static const PhysParamId_t kRHCNuisBANFFShape2p2h           = kNuis222;   // 2p2h nu shape
static const PhysParamId_t kRHCNuisBANFFNorm2p2hCtoO        = kNuis223;   // 2p2h CtoO norm
static const PhysParamId_t kRHCNuisBANFFEb                  = kNuis224;   // 2p2h CtoO norm
static const PhysParamId_t kRHCNuisBANFFMultiReNorm         = kNuis225;   // RHC MultiRe normalisation. NOT USED

//

static const PhysParamId_t kRHCLastBANFFParam  = kNuis226; // Utility, should always come after last RHC BANFF

   /*
    * Cross-section parameters (excluding cross-section params which were included in BANFF fit and are correlated with the flux params)
    * These other paramers are taken to be uncorrelated and are applied to SK only.
   */   
   //systematics which use splines

// NEW XSEC
static const PhysParamId_t kRHCNuisXSecRPA            = kNuis257;   // RPA tuning

static const PhysParamId_t kRHCNuisXSecFermiMom       = kNuis258;   // Fermi momentum
static const PhysParamId_t kRHCNuisXSecSpecFunc       = kNuis259;   // spectral function on/off
static const PhysParamId_t kRHCNuisXSecWshape         = kNuis260;   // W-shape
static const PhysParamId_t kRHCNuisXSecShapeCCoth     = kNuis261;   // CC other shape
static const PhysParamId_t kRHCNuisXSecPiLessDelta    = kNuis262;   // Pion-less delta decay
static const PhysParamId_t kRHCNuisXSecBindE          = kNuis263;   // Binding energy
   //systematics that are normalisations
static const PhysParamId_t kRHCNuisXSecNormCCcoh      = kNuis264;   // CC coherent norm
static const PhysParamId_t kRHCNuisXSecNormNC1piC     = kNuis265;   // NC 1piC norm
static const PhysParamId_t kRHCNuisXSecNormNCcoh      = kNuis266;   // NC coh norm
static const PhysParamId_t kRHCNuisXSecNormNCoth      = kNuis267;   // NC other norm
static const PhysParamId_t kRHCNuisXSecNormNuEToNuMu  = kNuis268;   // nue / numu cross section
static const PhysParamId_t kRHCNuisXSecNormAntiNuToNu = kNuis269;   // antinu / nu cross section

static const PhysParamId_t kRHCLastXSecParam = kNuis270;   // Utility, should always come after last RHC Xsec param

   /*
    * FSI parameters.
    * Assuming here a 10x10 covariance matrix will be made available.
    * Each parameter is a scaling factor to apply at some Ereco range (to be defined) for some
    * SK reaction mode (or the total? to be defined)
    */

// NEW SKDETFSI
 

#if defined(__USE_SKDETFSI2013__)
 // SKDET+FSI matrix 18x18
 #if FSI_BINS != 3
 static const PhysParamId_t kRHCNuisSKDetFSI000 = kNuis271;
 static const PhysParamId_t kRHCNuisSKDetFSI001 = kNuis272;
 static const PhysParamId_t kRHCNuisSKDetFSI002 = kNuis273;
 static const PhysParamId_t kRHCNuisSKDetFSI003 = kNuis274;
 static const PhysParamId_t kRHCNuisSKDetFSI004 = kNuis275;
 static const PhysParamId_t kRHCNuisSKDetFSI005 = kNuis276;
 #endif
 #if FSI_BINS != 2
 static const PhysParamId_t kRHCNuisSKDetFSI006 = kNuis277;
 static const PhysParamId_t kRHCNuisSKDetFSI007 = kNuis278;
 static const PhysParamId_t kRHCNuisSKDetFSI008 = kNuis279;
 static const PhysParamId_t kRHCNuisSKDetFSI009 = kNuis280;
 static const PhysParamId_t kRHCNuisSKDetFSI010 = kNuis281;
 static const PhysParamId_t kRHCNuisSKDetFSI011 = kNuis282;
 static const PhysParamId_t kRHCNuisSKDetFSI012 = kNuis283;
 static const PhysParamId_t kRHCNuisSKDetFSI013 = kNuis284;
 static const PhysParamId_t kRHCNuisSKDetFSI014 = kNuis285;
 static const PhysParamId_t kRHCNuisSKDetFSI015 = kNuis286;
 static const PhysParamId_t kRHCNuisSKDetFSI016 = kNuis287;
 static const PhysParamId_t kRHCNuisSKDetFSI017 = kNuis288;
#endif

 // NEW SKDETFSI
 #if FSI_BINS != 3
 static const PhysParamId_t kRHCNuisSKDetFSI_first = kRHCNuisSKDetFSI000; 
 #else
 static const PhysParamId_t kRHCNuisSKDetFSI_first = kRHCNuisSKDetFSI006; 
 #endif

#endif 

static const PhysParamId_t kRHCLastParam = kNuis295;






//*********  MULTIRING FORWARD HORN CURRENT **************************


/*
 * SuperK detector errors (energy scale and efficiencies)
 */   

/// NEW MULTIRING 3
static const PhysParamId_t kMultiRingNuisSuperkEnergyScale = kNuis306;   // SK energy scale
///

/// NEW MULTIRING 2

// NEW SKDETFSI

#if defined(__USE_SKDETFSI2013__)
 // SKDET+FSI matrix 18x18 
 #if MR_FSI_BINS != 3
 static const PhysParamId_t kMultiRingNuisSKDetFSI000 = kNuis307;
 static const PhysParamId_t kMultiRingNuisSKDetFSI001 = kNuis308;
 static const PhysParamId_t kMultiRingNuisSKDetFSI002 = kNuis309;
 static const PhysParamId_t kMultiRingNuisSKDetFSI003 = kNuis310;
 static const PhysParamId_t kMultiRingNuisSKDetFSI004 = kNuis311;
 static const PhysParamId_t kMultiRingNuisSKDetFSI005 = kNuis312;
#endif
 #if MR_FSI_BINS != 2
 static const PhysParamId_t kMultiRingNuisSKDetFSI006 = kNuis313;
 static const PhysParamId_t kMultiRingNuisSKDetFSI007 = kNuis314;
 static const PhysParamId_t kMultiRingNuisSKDetFSI008 = kNuis315;
 static const PhysParamId_t kMultiRingNuisSKDetFSI009 = kNuis316;
 static const PhysParamId_t kMultiRingNuisSKDetFSI010 = kNuis317;
 static const PhysParamId_t kMultiRingNuisSKDetFSI011 = kNuis318;
 static const PhysParamId_t kMultiRingNuisSKDetFSI012 = kNuis319;
 static const PhysParamId_t kMultiRingNuisSKDetFSI013 = kNuis320;
 static const PhysParamId_t kMultiRingNuisSKDetFSI014 = kNuis321;
 static const PhysParamId_t kMultiRingNuisSKDetFSI015 = kNuis322;
 static const PhysParamId_t kMultiRingNuisSKDetFSI016 = kNuis323;
 static const PhysParamId_t kMultiRingNuisSKDetFSI017 = kNuis324;
 #endif

 // NEW SKDETFSI
 #if MR_FSI_BINS != 3
 static const PhysParamId_t kMultiRingNuisSKDetFSI_first = kMultiRingNuisSKDetFSI000; 
 #else
  // Would normally be kMultiRingNuisSKDetFSI006 but don't have lowest energy bin
  // Set to kMultiRingNuisSKDetFSI007 to return correct description
 static const PhysParamId_t kMultiRingNuisSKDetFSI_first = kMultiRingNuisSKDetFSI007; 
 #endif
  
#endif

 
static const PhysParamId_t kMultiRingLastParam = kNuis331;
////

} // conv




namespace util {

  string AsString                   (conv::PhysParamId_t p);
  string AsStringTex                (conv::PhysParamId_t p);
  string AsStringPaperTexNoDollar   (conv::PhysParamId_t p);
  string AsStringPaperTex           (conv::PhysParamId_t p);
  string AsStringPaperRooTex        (conv::PhysParamId_t p);
  string AsStringPaperTable         (conv::PhysParamId_t p);
  string AsStringSimple             (conv::PhysParamId_t p);
  string AsStringCompact            (conv::PhysParamId_t p);
  string AsStringCompact_splinehack (conv::PhysParamId_t p);
  string AsStringCompact_KDI        (conv::PhysParamId_t p);
  string AsStringSpline             (conv::PhysParamId_t p);
  string AsStringNIWG               (conv::PhysParamId_t p);
  string AsStringOfficialAxis       (conv::PhysParamId_t p);

  bool IsSystematic  (conv::PhysParamId_t p);
  bool IsMixingAngle (conv::PhysParamId_t p);
  
  bool IsSystematic       (conv::PhysParamId_t p);
  bool IsFHCSystematic    (conv::PhysParamId_t p);
  bool IsRHCSystematic    (conv::PhysParamId_t p);
  /// NEW MULTIRING 2
  bool IsMultiRingSystematic    (conv::PhysParamId_t p);
  /// 
  bool IsMixingAngle          (conv::PhysParamId_t p);
  bool IsBANFF                (conv::PhysParamId_t p);
  bool IsRHCBANFF             (conv::PhysParamId_t p); 
  bool IsXSec                 (conv::PhysParamId_t p);
  bool IsUsingSplines         (conv::PhysParamId_t p);
  bool IsEff                  (conv::PhysParamId_t p);
  bool IsEff1Rmu              (conv::PhysParamId_t p);
  bool IsEff1Re               (conv::PhysParamId_t p);
  bool IsFSI                  (conv::PhysParamId_t p);
  bool IsSKDetFSI             (conv::PhysParamId_t p);
  bool IsFlux                 (conv::PhysParamId_t p);
  bool IsNormalisation        (conv::PhysParamId_t p);
  bool IsEnergyScale          (conv::PhysParamId_t p);
  bool IsBANFFUnconstrained   (conv::PhysParamId_t p);
  bool IsBANFFNotUnconstrained(conv::PhysParamId_t p);  
  bool IsRPA                  (conv::PhysParamId_t p);
  bool IsBeRPA                (conv::PhysParamId_t p);  
  bool IsSCC                  (conv::PhysParamId_t p);
  bool IsNormXSecParam        (conv::PhysParamId_t p);
  bool IsXSecNumuNue          (conv::PhysParamId_t p);
  bool IsEScale               (conv::PhysParamId_t p); 
  
  conv::PhysParamId_t      PhysParamIdFromString      (string h);
  conv::PhysParamId_t      PhysParamIdFromLongString  (string h);
  conv::PhysParamId_t      PhysParamIdFromSplineString(string h);  
  
  // Aliases to avoid overloading with enums
  string ParamAsString                   (conv::PhysParamId_t p);
  string ParamAsStringTex                (conv::PhysParamId_t p);
  string ParamAsStringPaperTexNoDollar   (conv::PhysParamId_t p);
  string ParamAsStringPaperTex           (conv::PhysParamId_t p);
  string ParamAsStringPaperRooTex        (conv::PhysParamId_t p);
  string ParamAsStringPaperTable         (conv::PhysParamId_t p);
  string ParamAsStringSimple             (conv::PhysParamId_t p);
  string ParamAsStringCompact            (conv::PhysParamId_t p);
  string ParamAsStringCompact_splinehack (conv::PhysParamId_t p);
  string ParamAsStringCompact_KDI        (conv::PhysParamId_t p);
  string ParamAsStringSpline             (conv::PhysParamId_t p);
  string ParamAsStringNIWG               (conv::PhysParamId_t p);

  conv::PhysParamId_t GetBarredVersionOfPhysParam(conv::PhysParamId_t p);
  
  // Tell the framework classes what our parameter names are
  void PopulateNameRegistry(void);
  
} // util

} // valor
  
#endif // _phys_param_id_h_

