//____________________________________________________________________________
/*!

\class     t2kana::conv::PhysParamId_t;

\brief     Conventions

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory
          
\created   Oct 20, 2010

\copyright (c) 2010-2014, The VALOR group

*/
//____________________________________________________________________________

#ifndef _phys_param_id_h_
#define _phys_param_id_h_

#include <string>
#include <AnalysisOptions.h>

using std::string;

namespace t2kana {

namespace conv {
   
// Max number of fit  parameters listed below (for booking large enough param arrays)
// Fit parameters
static const int knppmax = 200; // maximum number of physics params listed below (for booking large enough param arrays)

typedef enum EPhysParamId 
{   
   kPhypUndefined = -1, 
   kPhypMixing23,              // 23-mixing parameter (either sin^2(2*\theta_{23}) or sin^2(\theta_{23}) depending on config choices)
   kPhypMixing13,              // 13-mixing parameter (either sin^2(2*\theta_{13}) or sin^2(\theta_{13}) depending on config choices)
   kPhypMixing12,              // 12-mixing parameter (either sin^2(2*\theta_{12}) or sin^2(\theta_{12}) depending on config choices)
   kPhypCPV,                   // CP-invariance violating param (either deltacp, or sin(deltacp) depending on config choices)
   kPhypDm23sq,                // modulus of large squared-mass splitting (Dm^2_32, Dm^2_31 or Fogli-Lisi Dm^2 depending on config choices)
   kPhypDm12sq,                // small squared-mass splitting (Dm^2_21 or Fogli-Lisi dm^2 depending on config choices)
   kPhypAlphaDecay,            //
   kPhypSinsq2ThetaDecay,      // 
   kPhypMusqDecoh,             //
   kPhypSinsq2ThetaDecoh,      //
   kPhypMixing14,              // Sterile Mixings
   kPhypMixing15,              // Sterile Mixings
   kPhypMixing24,              // Sterile Mixings
   kPhypMixing25,              // Sterile Mixings
   kPhypMixing34,              // Sterile Mixings
   kPhypMixing35,              // Sterile Mixings
   kPhypDm41sq,                // Sterile mass splitting
   kPhypDm51sq,                // Sterile mass splitting
   kPhypDelta12,               // Sterile phase
   kPhypDelta13,               // Sterile phase - same as kPhypCPV, but separate for constistency.
   kPhypDelta15,               // Sterile phase
   kPhypDelta34,               // Sterile phase
   kPhypDelta35,               // Sterile phase
   kPhypNSIepsMuTau,           // Nonstandard matter interactions

   /*
    * List of 2011 analysis nuisance parameters 
    * -----------------------------------------
   */

   //********* FORWARD HORN CURRENT ************************ 

   /*
    * SuperK detector errors (energy scale and efficiencies)
   */   
   kNuisSuperkEnergyScale,     // SK energy scale
#if FSI_BINS != 3
   kNuisSuperkEff1RmuCCQEoth,  // SK 1Rmu CCQE efficiency (other)
   kNuisSuperkEff1RmuCCQErc0u, // SK 1Rmu CCQE efficiency (ring counting, Ereco bin 1)
   kNuisSuperkEff1RmuCCQErc1u, // SK 1Rmu CCQE efficiency (ring counting, Ereco bin 2)
   kNuisSuperkEff1RmuCCQErc2u, // SK 1Rmu CCQE efficiency (ring counting, Ereco bin 3)
   kNuisSuperkEff1RmuCCnonQE,  // SK 1Rmu CCnonQE efficiency
   kNuisSuperkEff1RmuNC,       // SK 1Rmu NC efficiency
   kNuisSuperkEff1RmuNueCC,    // SK 1Rmu nue CC efficiency
#endif

#if FSI_BINS != 2
   kNuisSuperkEff1ReNueSigCC0,    // SK 1Re CC eff for osc nue Ereco bin 0
   kNuisSuperkEff1ReNueSigCC1,    // SK 1Re CC eff for osc nue Ereco bin 1
   kNuisSuperkEff1ReNueSigCC2,    // SK 1Re CC eff for osc nue Ereco bin 2
   kNuisSuperkEff1ReNumuCC0,      // SK 1Re CC eff for numu Ereco bin 0
   kNuisSuperkEff1ReNumuCC1,      // SK 1Re CC eff for numu Ereco bin 1
   kNuisSuperkEff1ReNumuCC2,      // SK 1Re CC eff for numu Ereco bin 2
   kNuisSuperkEff1ReNueCC0,       // SK 1Re CC eff for nue Ereco bin 0
   kNuisSuperkEff1ReNueCC1,       // SK 1Re CC eff for nue Ereco bin 1
   kNuisSuperkEff1ReNueCC2,       // SK 1Re CC eff for nue Ereco bin 2
   kNuisSuperkEff1ReNC0,          // SK 1Re NC eff Ereco bin 0
   kNuisSuperkEff1ReNC1,          // SK 1Re NC eff Ereco bin 1
   kNuisSuperkEff1ReNC2,          // SK 1Re NC eff Ereco bin 2
#endif

   /*
    * (flux and cross-section) correlated parameters used in the BANFF fit
   */   
   kNuisBANFF00,               // # 00, SK     numu    flux,  E =  0.0 -  0.4 GeV
   kNuisBANFF01,               // # 01, SK     numu    flux,  E =  0.4 -  0.5 GeV
   kNuisBANFF02,               // # 02, SK     numu    flux,  E =  0.5 -  0.6 GeV
   kNuisBANFF03,               // # 03, SK     numu    flux,  E =  0.6 -  0.7 GeV
   kNuisBANFF04,               // # 04, SK     numu    flux,  E =  0.7 -  1.0 GeV
   kNuisBANFF05,               // # 05, SK     numu    flux,  E =  1.0 -  1.5 GeV
   kNuisBANFF06,               // # 06, SK     numu    flux,  E =  1.5 -  2.5 GeV
   kNuisBANFF07,               // # 07, SK     numu    flux,  E =  2.5 -  3.5 GeV
   kNuisBANFF08,               // # 08, SK     numu    flux,  E =  3.5 -  5.0 GeV
   kNuisBANFF09,               // # 09, SK     numu    flux,  E =  5.0 -  7.0 GeV
   kNuisBANFF10,               // # 10, SK     numu    flux,  E =  7.0 - 30.0 GeV
   kNuisBANFF11,               // # 11, different meanings depending on samples defined
   kNuisBANFF12,               // # 12, different meanings depending on samples defined
   kNuisBANFF13,               // # 13, different meanings depending on samples defined
   kNuisBANFF14,               // # 14, different meanings depending on samples defined
   kNuisBANFF15,               // # 15, different meanings depending on samples defined
#if (BANFF_BINS != 2)
   kNuisBANFF16,               // # 16, different meanings depending on samples defined
   kNuisBANFF17,               // # 17, different meanings depending on samples defined
   kNuisBANFF18,               // # 18, different meanings depending on samples defined
   kNuisBANFF19,               // # 19, different meanings depending on samples defined
   kNuisBANFF20,               // # 20, different meanings depending on samples defined
   kNuisBANFF21,               // # 21, different meanings depending on samples defined
#endif
#if BANFF_BINS == 1
   kNuisBANFF22,               // # 22, different meanings depending on samples defined
   kNuisBANFF23,               // # 23, different meanings depending on samples defined
   kNuisBANFF24,               // # 24, different meanings depending on samples defined
#endif
   kNuisBANFFMAQE,             // # M_{A}  QE
   kNuisBANFFMARES,            // # M_{A}  RES
   kNuisBANFFCCQEE0,           // # CCQE   normalization, 0 - 1.5 GeV
#if BANFF_BINS != 3
   kNuisBANFFCCQEE1,           // # CCQE   normalization, 1.5 -  3.5 GeV
   kNuisBANFFCCQEE2,           // # CCQE   normalization, 3.5 - 30.0 GeV
#endif
   kNuisBANFFCC1piE0,          // # CC1pi  normalization, 0 - 2.5 GeV
#if BANFF_BINS != 3
   kNuisBANFFCC1piE1,          // # CC1pi  normalization, 2.5 - 30.0 GeV
#endif
#if (BANFF_BINS != 2)
   kNuisBANFFNC1pi0,           // # NC1pi0 normalization, 0 - 30.0 GeV
#endif

   /*
    * Cross-section parameters (excluding cross-section params which were included in BANFF fit and are correlated with the flux params)
    * These other paramers are taken to be uncorrelated and are applied to SK only.
   */   
   //systematics which use splines
   kNuisXSecFermiMom,          // Fermi momentum
   kNuisXSecSpecFunc,          // spectral function on/off
   kNuisXSecWshape,            // W-shape
   kNuisXSecShapeCCoth,        // CC other shape
   kNuisXSecPiLessDelta,       // Pion-less delta decay
   kNuisXSecBindE,             // Binding energy
   //systematics that are normalisations
   kNuisXSecNormCCcoh,         // CC coherent norm
   kNuisXSecNormNC1piC,        // NC 1piC norm
   kNuisXSecNormNCcoh,         // NC coh norm
   kNuisXSecNormNCoth,         // NC other norm
   kNuisXSecNormNuEToNuMu,     // nue / numu cross section
   kNuisXSecNormAntiNuToNu,    // antinu / nu cross section

   /*
    * FSI parameters.
    * Assuming here a 10x10 covariance matrix will be made available.
    * Each parameter is a scaling factor to apply at some Ereco range (to be defined) for some
    * SK reaction mode (or the total? to be defined)
    */

   // FSI matrix 6x6 for 1Rmu
#if FSI_BINS != 3
   kNuisFSI000,
   kNuisFSI001,
   kNuisFSI002,
   kNuisFSI003,
   kNuisFSI004,
   kNuisFSI005,
#endif
   // FSI matrix 12x12 for 1Re
#if FSI_BINS != 2
   kNuisFSI006,
   kNuisFSI007,
   kNuisFSI008,
   kNuisFSI009,
   kNuisFSI010,
   kNuisFSI011,
   kNuisFSI012,
   kNuisFSI013,
   kNuisFSI014,
   kNuisFSI015,
   kNuisFSI016,
   kNuisFSI017,
#endif

   // SKDET+FSI matrix 18x18
#if FSI_BINS != 3
   kNuisSKDetFSI000,
   kNuisSKDetFSI001,
   kNuisSKDetFSI002,
   kNuisSKDetFSI003,
   kNuisSKDetFSI004,
   kNuisSKDetFSI005,
#endif
#if FSI_BINS != 2
   kNuisSKDetFSI006,
   kNuisSKDetFSI007,
   kNuisSKDetFSI008,
   kNuisSKDetFSI009,
   kNuisSKDetFSI010,
   kNuisSKDetFSI011,
   kNuisSKDetFSI012,
   kNuisSKDetFSI013,
   kNuisSKDetFSI014,
   kNuisSKDetFSI015,
   kNuisSKDetFSI016,
   kNuisSKDetFSI017,
#endif

   //********* REVERSE HORN CURRENT **************************

   /*
    * Special RHC to FHC errors
    */

   kNuisRHCToFHC,

   /*
    * SuperK detector errors (energy scale and efficiencies)
    */
   
   kRHCNuisSuperkEnergyScale,     // SK energy scale
   
   /*                                               
    * (flux and cross-section) correlated parameters used in the BANFF fit               
    */
   
   kRHCNuisBANFF00,               // # 00, SK     numu    flux,  E =  0.0 -  0.4 GeV
   kRHCNuisBANFF01,               // # 01, SK     numu    flux,  E =  0.4 -  0.5 GeV
   kRHCNuisBANFF02,               // # 02, SK     numu    flux,  E =  0.5 -  0.6 GeV   
   kRHCNuisBANFF03,               // # 03, SK     numu    flux,  E =  0.6 -  0.7 GeV
   kRHCNuisBANFF04,               // # 04, SK     numu    flux,  E =  0.7 -  1.0 GeV  
   kRHCNuisBANFF05,               // # 05, SK     numu    flux,  E =  1.0 -  1.5 GeV
   kRHCNuisBANFF06,               // # 06, SK     numu    flux,  E =  1.5 -  2.5 GeV 
   kRHCNuisBANFF07,               // # 07, SK     numu    flux,  E =  2.5 -  3.5 GeV
   kRHCNuisBANFF08,               // # 08, SK     numu    flux,  E =  3.5 -  5.0 GeV
   kRHCNuisBANFF09,               // # 09, SK     numu    flux,  E =  5.0 -  7.0 GeV
   kRHCNuisBANFF10,               // # 10, SK     numu    flux,  E =  7.0 - 30.0 GeV 
   kRHCNuisBANFF11,               // # 11, different meanings depending on samples defined
   kRHCNuisBANFF12,               // # 12, different meanings depending on samples defined
   kRHCNuisBANFF13,               // # 13, different meanings depending on samples defined       
   kRHCNuisBANFF14,               // # 14, different meanings depending on samples defined      
   kRHCNuisBANFF15,               // # 15, different meanings depending on samples defined       
#if (BANFF_BINS != 2)
   kRHCNuisBANFF16,               // # 16, different meanings depending on samples defined
   kRHCNuisBANFF17,               // # 17, different meanings depending on samples defined
   kRHCNuisBANFF18,               // # 18, different meanings depending on samples defined
   kRHCNuisBANFF19,               // # 19, different meanings depending on samples defined       
   kRHCNuisBANFF20,               // # 20, different meanings depending on samples defined
   kRHCNuisBANFF21,               // # 21, different meanings depending on samples defined      
#endif
#if BANFF_BINS == 1
   kRHCNuisBANFF22,               // # 22, different meanings depending on samples defined
   kRHCNuisBANFF23,               // # 23, different meanings depending on samples defined       
   kRHCNuisBANFF24,               // # 24, different meanings depending on samples defined
#endif
   kRHCNuisBANFFMAQE,             // # M_{A}  QE                                     
   kRHCNuisBANFFMARES,            // # M_{A}  RES                                   
   kRHCNuisBANFFCCQEE0,           // # CCQE   normalization, 0 - 1.5 GeV   
#if BANFF_BINS != 3
   kRHCNuisBANFFCCQEE1,           // # CCQE   normalization, 1.5 -  3.5 GeV        
   kRHCNuisBANFFCCQEE2,           // # CCQE   normalization, 3.5 - 30.0 GeV        
#endif
   kRHCNuisBANFFCC1piE0,          // # CC1pi  normalization, 0 - 2.5 GeV           
#if BANFF_BINS != 3
   kRHCNuisBANFFCC1piE1,          // # CC1pi  normalization, 2.5 - 30.0 GeV       
#endif
#if (BANFF_BINS != 2)
   kRHCNuisBANFFNC1pi0,           // # NC1pi0 normalization, 0 - 30.0 GeV         
#endif
   
   /*                                                                             
    * Cross-section parameters (excluding cross-section params which were included in BANFF fit and are correlated with the flux params)
    * These other paramers are taken to be uncorrelated and are applied to SK only.     
    */
   
   //systematics which use splines                                               
   kRHCNuisXSecFermiMom,          // Fermi momentum                                 
   kRHCNuisXSecSpecFunc,          // spectral function on/off         
   kRHCNuisXSecWshape,            // W-shape                          
   kRHCNuisXSecShapeCCoth,        // CC other shape               
   kRHCNuisXSecPiLessDelta,       // Pion-less delta decay                
   kRHCNuisXSecBindE,             // Binding energy               

   //systematics that are normalisations                                         
   kRHCNuisXSecNormCCcoh,         // CC coherent norm                              
   kRHCNuisXSecNormNC1piC,        // NC 1piC norm                                   
   kRHCNuisXSecNormNCcoh,         // NC coh norm                                   
   kRHCNuisXSecNormNCoth,         // NC other norm                       
   kRHCNuisXSecNormNuEToNuMu,     // nue / numu cross section                  
   kRHCNuisXSecNormAntiNuToNu,    // antinu / nu cross section        
   
   /*                                                                          
    * SK and FSI parameters.                                 
    * Assuming here a 10x10 covariance matrix will be made available.            
    * Each parameter is a scaling factor to apply at some Ereco range (to be defined) for some
    * SK reaction mode (or the total? to be defined)                              
    */
   
   // SKDET+FSI matrix 18x18                                                  
#if FSI_BINS != 3
   kRHCNuisSKDetFSI000,
   kRHCNuisSKDetFSI001,
   kRHCNuisSKDetFSI002,
   kRHCNuisSKDetFSI003,
   kRHCNuisSKDetFSI004,
   kRHCNuisSKDetFSI005,
#endif
#if FSI_BINS != 2
   kRHCNuisSKDetFSI006,
   kRHCNuisSKDetFSI007,
   kRHCNuisSKDetFSI008,
   kRHCNuisSKDetFSI009,
   kRHCNuisSKDetFSI010,
   kRHCNuisSKDetFSI011,
   kRHCNuisSKDetFSI012,
   kRHCNuisSKDetFSI013,
   kRHCNuisSKDetFSI014,
   kRHCNuisSKDetFSI015,
   kRHCNuisSKDetFSI016,
   kRHCNuisSKDetFSI017,
#endif
   kRHCLastParam
   
} 
PhysParamId_t;

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
  
  bool GlobalFit(conv::PhysParamId_t p, double & bestfit_out, double & error_out);


} // util

} // t2kana
  
#endif // _phys_param_id_h_

