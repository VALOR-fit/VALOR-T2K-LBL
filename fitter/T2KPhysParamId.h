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
#include <AnalysisOptions.h>

#include "Conventions/PhysParamId.h"

using std::string;

namespace valor {

namespace conv {
   
   
   /*
    * List of 2011 analysis nuisance parameters 
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
#if (BANFF_BINS != 2)
static const PhysParamId_t kNuisBANFF16 = kNuis036; // # 16; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF17 = kNuis037; // # 17; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF18 = kNuis038; // # 18; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF19 = kNuis039; // # 19; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF20 = kNuis040; // # 20; different meanings depending on samples defined
static const PhysParamId_t kNuisBANFF21 = kNuis041; // # 21; different meanings depending on samples defined
#endif
#if BANFF_BINS == 1
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

// NEW XSEC
static const PhysParamId_t kNuisBANFFMEC      = kNuis053;  // MEC normalization
static const PhysParamId_t kNuisBANFFCA5      = kNuis054;  // CA5 RES 
static const PhysParamId_t kNuisBANFFBgRES    = kNuis055;  // Isospin 1/2 normalization
//

static const PhysParamId_t kLastBANFFParam  = kNuis056; // Utility, should always come after last FHC BANFF


 


   /*
    * Cross-section parameters (excluding cross-section params which were included in BANFF fit and are correlated with the flux params)
    * These other paramers are taken to be uncorrelated and are applied to SK only.
   */   
   //systematics which use splines
static const PhysParamId_t kNuisXSecFermiMom       = kNuis060;   // Fermi momentum
static const PhysParamId_t kNuisXSecSpecFunc       = kNuis061;   // spectral function on/off
static const PhysParamId_t kNuisXSecWshape         = kNuis062;   // W-shape
static const PhysParamId_t kNuisXSecShapeCCoth     = kNuis063;   // CC other shape
static const PhysParamId_t kNuisXSecPiLessDelta    = kNuis064;   // Pion-less delta decay
static const PhysParamId_t kNuisXSecBindE          = kNuis065;   // Binding energy
   //systematics that are normalisations
static const PhysParamId_t kNuisXSecNormCCcoh      = kNuis066;   // CC coherent norm
static const PhysParamId_t kNuisXSecNormNC1piC     = kNuis067;   // NC 1piC norm
static const PhysParamId_t kNuisXSecNormNCcoh      = kNuis068;   // NC coh norm
static const PhysParamId_t kNuisXSecNormNCoth      = kNuis069;   // NC other norm
static const PhysParamId_t kNuisXSecNormNuEToNuMu  = kNuis070;   // nue / numu cross section
static const PhysParamId_t kNuisXSecNormAntiNuToNu = kNuis071;   // antinu / nu cross section


static const PhysParamId_t kLastXSecParam = kNuis072;   // Utility, should always come after last FHC Xsec param

   /*
    * FSI parameters.
    * Assuming here a 10x10 covariance matrix will be made available.
    * Each parameter is a scaling factor to apply at some Ereco range (to be defined) for some
    * SK reaction mode (or the total? to be defined)
    */

   // SKDET+FSI matrix 18x18
#if FSI_BINS != 3
static const PhysParamId_t kNuisSKDetFSI000 = kNuis080;
static const PhysParamId_t kNuisSKDetFSI001 = kNuis081;
static const PhysParamId_t kNuisSKDetFSI002 = kNuis082;
static const PhysParamId_t kNuisSKDetFSI003 = kNuis083;
static const PhysParamId_t kNuisSKDetFSI004 = kNuis084;
static const PhysParamId_t kNuisSKDetFSI005 = kNuis085;
#endif
#if FSI_BINS != 2
static const PhysParamId_t kNuisSKDetFSI006 = kNuis086;
static const PhysParamId_t kNuisSKDetFSI007 = kNuis087;
static const PhysParamId_t kNuisSKDetFSI008 = kNuis088;
static const PhysParamId_t kNuisSKDetFSI009 = kNuis089;
static const PhysParamId_t kNuisSKDetFSI010 = kNuis090;
static const PhysParamId_t kNuisSKDetFSI011 = kNuis091;
static const PhysParamId_t kNuisSKDetFSI012 = kNuis092;
static const PhysParamId_t kNuisSKDetFSI013 = kNuis093;
static const PhysParamId_t kNuisSKDetFSI014 = kNuis094;
static const PhysParamId_t kNuisSKDetFSI015 = kNuis095;
static const PhysParamId_t kNuisSKDetFSI016 = kNuis096;
static const PhysParamId_t kNuisSKDetFSI017 = kNuis097;
#endif

   //********* REVERSE HORN CURRENT **************************

   /*
    * Special RHC to FHC errors
    */

static const PhysParamId_t kNuisRHCToFHC = kNuis099;

   /*
    * SuperK detector errors (energy scale and efficiencies)
    */
   
   /*
    * SuperK detector errors (energy scale and efficiencies)
   */   
static const PhysParamId_t kRHCNuisSuperkEnergyScale     = kNuis100;     // SK energy scale
#if FSI_BINS != 3
static const PhysParamId_t kRHCNuisSuperkEff1RmuCCQEoth  = kNuis101; // SK 1Rmu CCQE efficiency (other)
static const PhysParamId_t kRHCNuisSuperkEff1RmuCCQErc0u = kNuis102; // SK 1Rmu CCQE efficiency (ring counting; Ereco bin 1)
static const PhysParamId_t kRHCNuisSuperkEff1RmuCCQErc1u = kNuis103; // SK 1Rmu CCQE efficiency (ring counting; Ereco bin 2)
static const PhysParamId_t kRHCNuisSuperkEff1RmuCCQErc2u = kNuis104; // SK 1Rmu CCQE efficiency (ring counting; Ereco bin 3)
static const PhysParamId_t kRHCNuisSuperkEff1RmuCCnonQE  = kNuis105; // SK 1Rmu CCnonQE efficiency
static const PhysParamId_t kRHCNuisSuperkEff1RmuNC       = kNuis106; // SK 1Rmu NC efficiency
static const PhysParamId_t kRHCNuisSuperkEff1RmuNueCC    = kNuis107; // SK 1Rmu nue CC efficiency
#endif

#if FSI_BINS != 2
static const PhysParamId_t kRHCNuisSuperkEff1ReNueSigCC0 = kNuis101;    // SK 1Re CC eff for osc nue Ereco bin 0
static const PhysParamId_t kRHCNuisSuperkEff1ReNueSigCC1 = kNuis102;    // SK 1Re CC eff for osc nue Ereco bin 1
static const PhysParamId_t kRHCNuisSuperkEff1ReNueSigCC2 = kNuis103;    // SK 1Re CC eff for osc nue Ereco bin 2
static const PhysParamId_t kRHCNuisSuperkEff1ReNumuCC0   = kNuis104;    // SK 1Re CC eff for numu Ereco bin 0
static const PhysParamId_t kRHCNuisSuperkEff1ReNumuCC1   = kNuis105;    // SK 1Re CC eff for numu Ereco bin 1
static const PhysParamId_t kRHCNuisSuperkEff1ReNumuCC2   = kNuis106;    // SK 1Re CC eff for numu Ereco bin 2
static const PhysParamId_t kRHCNuisSuperkEff1ReNueCC0    = kNuis107;    // SK 1Re CC eff for nue Ereco bin 0
static const PhysParamId_t kRHCNuisSuperkEff1ReNueCC1    = kNuis108;    // SK 1Re CC eff for nue Ereco bin 1
static const PhysParamId_t kRHCNuisSuperkEff1ReNueCC2    = kNuis109;    // SK 1Re CC eff for nue Ereco bin 2
static const PhysParamId_t kRHCNuisSuperkEff1ReNC0       = kNuis110;    // SK 1Re NC eff Ereco bin 0
static const PhysParamId_t kRHCNuisSuperkEff1ReNC1       = kNuis111;    // SK 1Re NC eff Ereco bin 1
static const PhysParamId_t kRHCNuisSuperkEff1ReNC2       = kNuis112;    // SK 1Re NC eff Ereco bin 2
#endif

   /*
    * (flux and cross-section) correlated parameters used in the BANFF fit
   */   
static const PhysParamId_t kRHCNuisBANFF00 = kNuis120; // # 00; SK numu
static const PhysParamId_t kRHCNuisBANFF01 = kNuis121; // # 01; SK numu
static const PhysParamId_t kRHCNuisBANFF02 = kNuis122; // # 02; SK numu
static const PhysParamId_t kRHCNuisBANFF03 = kNuis123; // # 03; SK numu 
static const PhysParamId_t kRHCNuisBANFF04 = kNuis124; // # 04; SK numu
static const PhysParamId_t kRHCNuisBANFF05 = kNuis125; // # 05; SK numu 
static const PhysParamId_t kRHCNuisBANFF06 = kNuis126; // # 06; SK anti-numu
static const PhysParamId_t kRHCNuisBANFF07 = kNuis127; // # 07; SK anti-numu
static const PhysParamId_t kRHCNuisBANFF08 = kNuis128; // # 08; SK anti-numu 
static const PhysParamId_t kRHCNuisBANFF09 = kNuis129; // # 09; SK anti-numu 
static const PhysParamId_t kRHCNuisBANFF10 = kNuis130; // # 10; SK anti-numu 
static const PhysParamId_t kRHCNuisBANFF11 = kNuis131; // # 11; SK anti-numu (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF12 = kNuis132; // # 12; SK anti-numu (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF13 = kNuis133; // # 13; SK anti-numu (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF14 = kNuis134; // # 14; SK anti-numu (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF15 = kNuis135; // # 15; SK anti-numu (different meanings depending on samples defined)
#if (BANFF_BINS != 2)
static const PhysParamId_t kRHCNuisBANFF16 = kNuis136; // # 16; SK nue (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF17 = kNuis137; // # 17; SK nue (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF18 = kNuis138; // # 18; SK anti-nue (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF19 = kNuis139; // # 19; SK anti-nu (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF20 = kNuis140; // # 20; SK anti-nue (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF21 = kNuis141; // # 21; SK anti-nue (different meanings depending on samples defined)
#endif
#if BANFF_BINS == 1
static const PhysParamId_t kRHCNuisBANFF22 = kNuis142; // # 22; SK anti-nue (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF23 = kNuis143; // # 23; SK anti-nue (different meanings depending on samples defined)
static const PhysParamId_t kRHCNuisBANFF24 = kNuis144; // # 24; SK anti-nue (different meanings depending on samples defined)
#endif
static const PhysParamId_t kRHCNuisBANFFMAQE    = kNuis145;  // # M_{A}  QE
static const PhysParamId_t kRHCNuisBANFFMARES   = kNuis146;  // # M_{A}  RES
static const PhysParamId_t kRHCNuisBANFFCCQEE0  = kNuis147;  // # CCQE   normalization; 0 - 1.5 GeV
#if BANFF_BINS != 3
static const PhysParamId_t kRHCNuisBANFFCCQEE1  = kNuis148;  // # CCQE   normalization; 1.5 -  3.5 GeV
static const PhysParamId_t kRHCNuisBANFFCCQEE2  = kNuis149;  // # CCQE   normalization; 3.5 - 30.0 GeV
#endif
static const PhysParamId_t kRHCNuisBANFFCC1piE0 = kNuis150; // # CC1pi  normalization; 0 - 2.5 GeV
#if BANFF_BINS != 3
static const PhysParamId_t kRHCNuisBANFFCC1piE1 = kNuis151; // # CC1pi  normalization; 2.5 - 30.0 GeV
#endif
#if (BANFF_BINS != 2)
static const PhysParamId_t kRHCNuisBANFFNC1pi0  = kNuis152; // # NC1pi0 normalization; 0 - 30.0 GeV
#endif

// NEW XSEC
static const PhysParamId_t kRHCNuisBANFFMEC     = kNuis153;  // MEC normalization
static const PhysParamId_t kRHCNuisBANFFCA5     = kNuis154;  // CA5 RES 
static const PhysParamId_t kRHCNuisBANFFBgRES   = kNuis155;  // Isospin 1/2 normalization
//

static const PhysParamId_t kRHCLastBANFFParam  = kNuis156; // Utility, should always come after last RHC BANFF


   /*
    * Cross-section parameters (excluding cross-section params which were included in BANFF fit and are correlated with the flux params)
    * These other paramers are taken to be uncorrelated and are applied to SK only.
   */   
   //systematics which use splines
static const PhysParamId_t kRHCNuisXSecFermiMom       = kNuis160;   // Fermi momentum
static const PhysParamId_t kRHCNuisXSecSpecFunc       = kNuis161;   // spectral function on/off
static const PhysParamId_t kRHCNuisXSecWshape         = kNuis162;   // W-shape
static const PhysParamId_t kRHCNuisXSecShapeCCoth     = kNuis163;   // CC other shape
static const PhysParamId_t kRHCNuisXSecPiLessDelta    = kNuis164;   // Pion-less delta decay
static const PhysParamId_t kRHCNuisXSecBindE          = kNuis165;   // Binding energy
   //systematics that are normalisations
static const PhysParamId_t kRHCNuisXSecNormCCcoh      = kNuis166;   // CC coherent norm
static const PhysParamId_t kRHCNuisXSecNormNC1piC     = kNuis167;   // NC 1piC norm
static const PhysParamId_t kRHCNuisXSecNormNCcoh      = kNuis168;   // NC coh norm
static const PhysParamId_t kRHCNuisXSecNormNCoth      = kNuis169;   // NC other norm
static const PhysParamId_t kRHCNuisXSecNormNuEToNuMu  = kNuis170;   // nue / numu cross section
static const PhysParamId_t kRHCNuisXSecNormAntiNuToNu = kNuis171;   // antinu / nu cross section

static const PhysParamId_t kRHCLastXSecParam = kNuis172;   // Utility, should always come after last RHC Xsec param

   /*
    * FSI parameters.
    * Assuming here a 10x10 covariance matrix will be made available.
    * Each parameter is a scaling factor to apply at some Ereco range (to be defined) for some
    * SK reaction mode (or the total? to be defined)
    */

   // SKDET+FSI matrix 18x18
#if FSI_BINS != 3
static const PhysParamId_t kRHCNuisSKDetFSI000 = kNuis180;
static const PhysParamId_t kRHCNuisSKDetFSI001 = kNuis181;
static const PhysParamId_t kRHCNuisSKDetFSI002 = kNuis182;
static const PhysParamId_t kRHCNuisSKDetFSI003 = kNuis183;
static const PhysParamId_t kRHCNuisSKDetFSI004 = kNuis184;
static const PhysParamId_t kRHCNuisSKDetFSI005 = kNuis185;
#endif
#if FSI_BINS != 2
static const PhysParamId_t kRHCNuisSKDetFSI006 = kNuis186;
static const PhysParamId_t kRHCNuisSKDetFSI007 = kNuis187;
static const PhysParamId_t kRHCNuisSKDetFSI008 = kNuis188;
static const PhysParamId_t kRHCNuisSKDetFSI009 = kNuis189;
static const PhysParamId_t kRHCNuisSKDetFSI010 = kNuis190;
static const PhysParamId_t kRHCNuisSKDetFSI011 = kNuis191;
static const PhysParamId_t kRHCNuisSKDetFSI012 = kNuis192;
static const PhysParamId_t kRHCNuisSKDetFSI013 = kNuis193;
static const PhysParamId_t kRHCNuisSKDetFSI014 = kNuis194;
static const PhysParamId_t kRHCNuisSKDetFSI015 = kNuis195;
static const PhysParamId_t kRHCNuisSKDetFSI016 = kNuis196;
static const PhysParamId_t kRHCNuisSKDetFSI017 = kNuis197;
#endif

static const PhysParamId_t kRHCLastParam = kNuis198;






//*********  MULTIRING FORWARD HORN CURRENT **************************


/*
 * SuperK detector errors (energy scale and efficiencies)
 */   

/// NEW MULTIRING 3
static const PhysParamId_t kMultiRingNuisSuperkEnergyScale = kNuis199;   // SK energy scale
///

/// NEW MULTIRING 2

// SKDET+FSI matrix 18x18 
#if FSI_BINS != 3
static const PhysParamId_t kMultiRingNuisSKDetFSI000 = kNuis200;
static const PhysParamId_t kMultiRingNuisSKDetFSI001 = kNuis201;
static const PhysParamId_t kMultiRingNuisSKDetFSI002 = kNuis202;
static const PhysParamId_t kMultiRingNuisSKDetFSI003 = kNuis203;
static const PhysParamId_t kMultiRingNuisSKDetFSI004 = kNuis204;
static const PhysParamId_t kMultiRingNuisSKDetFSI005 = kNuis205;
#endif
#if FSI_BINS != 2
static const PhysParamId_t kMultiRingNuisSKDetFSI006 = kNuis206;
static const PhysParamId_t kMultiRingNuisSKDetFSI007 = kNuis207;
static const PhysParamId_t kMultiRingNuisSKDetFSI008 = kNuis208;
static const PhysParamId_t kMultiRingNuisSKDetFSI009 = kNuis209;
static const PhysParamId_t kMultiRingNuisSKDetFSI010 = kNuis210;
static const PhysParamId_t kMultiRingNuisSKDetFSI011 = kNuis211;
static const PhysParamId_t kMultiRingNuisSKDetFSI012 = kNuis212;
static const PhysParamId_t kMultiRingNuisSKDetFSI013 = kNuis213;
static const PhysParamId_t kMultiRingNuisSKDetFSI014 = kNuis214;
static const PhysParamId_t kMultiRingNuisSKDetFSI015 = kNuis215;
static const PhysParamId_t kMultiRingNuisSKDetFSI016 = kNuis216;
static const PhysParamId_t kMultiRingNuisSKDetFSI017 = kNuis217;
#endif

static const PhysParamId_t kMultiRingLastParam = kNuis218;
////

} // conv




namespace util {

  string AsString                   (conv::PhysParamId_t p);
  string AsStringTex                (conv::PhysParamId_t p);
  string AsStringPaperTexNoDollar   (conv::PhysParamId_t p);
  string AsStringPaperTex           (conv::PhysParamId_t p);
  string AsStringPaperRooTex        (conv::PhysParamId_t p);
  string AsStringCompact            (conv::PhysParamId_t p);
  string AsStringCompact_splinehack (conv::PhysParamId_t p);
  string AsStringCompact_KDI        (conv::PhysParamId_t p);

  bool IsSystematic  (conv::PhysParamId_t p);
  bool IsMixingAngle (conv::PhysParamId_t p);
  
  bool IsSystematic       (conv::PhysParamId_t p);
  bool IsFHCSystematic    (conv::PhysParamId_t p);
  bool IsRHCSystematic    (conv::PhysParamId_t p);
  /// NEW MULTIRING 2
  bool IsMultiRingSystematic    (conv::PhysParamId_t p);
  /// 
  bool IsMixingAngle      (conv::PhysParamId_t p);
  bool IsBANFF            (conv::PhysParamId_t p);
  bool IsXSec             (conv::PhysParamId_t p);
  bool IsUsingSplines     (conv::PhysParamId_t p);
  bool IsEff              (conv::PhysParamId_t p);
  bool IsEff1Rmu          (conv::PhysParamId_t p);
  bool IsEff1Re           (conv::PhysParamId_t p);
  bool IsFSI              (conv::PhysParamId_t p);
  bool IsSKDetFSI         (conv::PhysParamId_t p);
  
  conv::PhysParamId_t      PhysParamIdFromString     (string h);
  conv::PhysParamId_t      PhysParamIdFromLongString (string h);
  
  // Aliases to avoid overloading with enums
  string ParamAsString                   (conv::PhysParamId_t p);
  string ParamAsStringTex                (conv::PhysParamId_t p);
  string ParamAsStringPaperTexNoDollar   (conv::PhysParamId_t p);
  string ParamAsStringPaperTex           (conv::PhysParamId_t p);
  string ParamAsStringPaperRooTex        (conv::PhysParamId_t p);
  string ParamAsStringCompact            (conv::PhysParamId_t p);
  string ParamAsStringCompact_splinehack (conv::PhysParamId_t p);
  string ParamAsStringCompact_KDI        (conv::PhysParamId_t p);
  
} // util

} // valor
  
#endif // _phys_param_id_h_

