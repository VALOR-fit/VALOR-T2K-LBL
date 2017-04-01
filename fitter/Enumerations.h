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
Add samples, modes and systematic parameters for reverse horn current running

21-7-2014 Davide
Added samples for multiring forward horn current


*/
//____________________________________________________________________________

#ifndef enumerations_h
#define enumerations_h

#include "AnalysisOptions.h"

#include "T2KPhysParamId.h"

#include "valor/Conventions/ConfLevel.h"
#include "valor/Conventions/NuMassHierarchy.h"
#include "valor/Conventions/NuOscDeltaCP.h"
#include "valor/Conventions/NuOscDmsqScheme.h"
#include "valor/Conventions/NuOscMixParam.h"
#include "valor/Conventions/NuOscOctant.h"
#include "valor/Conventions/OscChannelId.h"
#include "valor/Conventions/PDG.h"
#include "valor/Conventions/PhysHypothesisId.h"
//#include "valor/Conventions/PhysParamId.h"
//#include "valor/Conventions/ReacMode.h"

#include "SuperkModes.h"  
#include "SuperkSamples.h"

namespace valor {
namespace conv   {

//
// ------------------------------------------------------------------------
// Enumerations
// ------------------------------------------------------------------------
//

// Combined or single analysis 
typedef enum EOscAnaType 
{   
   kAnaUndefined = -1, 

   kAnaCombined,
   kAnaNuMuDisappOnly,
   kAnaNuEAppOnly,

   kAnaRHCCombined,
   kAnaRHCNuMuDisappOnly,
   kAnaRHCNuEAppOnly,

   kAnaFHCAndRHCNuMuDisappOnly,
   kAnaFHCAndRHCNuEAppOnly,
   kAnaFHCAndRHCCombined,

   // MultiRing-FHC
   kAnaMultiRingCombined,
   kAnaMultiRingNuMuDisappOnly,
   kAnaMultiRingNuEAppOnly,

   // MultiRing-FHC + 1Ring-FHC
   kAnaMultiRingAnd1RingCombined,
   kAnaMultiRingAnd1RingNuMuDisappOnly,
   kAnaMultiRingAnd1RingNuEAppOnly,

   // MultiRing-FHC + 1Ring-FHC&RHC
   kAnaMultiRingFHCAnd1RingFHCAndRHCCombined,
   kAnaMultiRingFHCAnd1RingFHCAndRHCNuMuDisappOnly,
   kAnaMultiRingFHCAnd1RingFHCAndRHCNuEAppOnly,

} 
OscAnaType_t;

// Type of SK MC template
typedef enum ESKMCTemplateType
{   
   kSKMCtUndefined = -1, 
   kSKMCtOri, // nominal MC
   kSKMCtPOT, // POT and nominal flux weights applied. No NIWG tune.
   kSKMCtNom, // nominal MC + tuning (ext) weights 
   kSKMCtRPA, // same as above, but with (relativistic rather than nonrelativistic RPA) || (Bernstein polynomial RPA)
   kSKMCtTun, // nominal MC + tuning (ext) weights + tuning (T2K) weights
   kSKMCtSys, // nominal MC + tuning (ext) weights + tuning (T2K) weights + systematic variations
   kSKMCtFin  // nominal MC + tuning (ext) weights + tuning (T2K) weights + systematic variations + physics hypothesis 
} 
SKMCTemplateType_t;

// Enumeration of flux histograms types held by NuFlux
typedef enum EFluxType
{   
   kFlxUndefined      = -1, 
   kFlxNominal        =  0,  // flux used for MC production
   kFlxExtTuning      =  1,  // flux used for MC production + tuning wghts based on external data (NA61,...)
   kFlxFinalTuning    =  2,  // flux used for MC production + tuning wghts based on external data (NA61,...) + tuning wghts based on T2K data (ND280 data fit by BANFF,...)
   kFlxCurrent        =  3,  // current flux from systematic param variation of kFlxFinalTuning flux
   kRHCFlxNominal     =  4,  // RHC flux used for MC production
   kRHCFlxExtTuning   =  5,  // RHC flux used for MC production + tuning wghts based on external data (NA61,...)
   kRHCFlxFinalTuning =  6,  // RHC flux used for MC production + tuning wghts based on external data (NA61,...) + tuning wghts based on T2K data (ND280 data fit by BANFF,...)
   kRHCFlxCurrent     =  7,   // RHC current flux from systematic param variation of kFlxFinalTuning flux
   
   /* /// NEW MULTIRING */
   /* kMultiRingFlxNominal     =  8,  // RHC flux used for MC production */
   /* kMultiRingFlxExtTuning   =  9,  // RHC flux used for MC production + tuning wghts based on external data (NA61,...) */
   /* kMultiRingFlxFinalTuning =  10,  // RHC flux used for MC production + tuning wghts based on external data (NA61,...) + tuning wghts based on T2K data (ND280 data fit by BANFF,...) */
   /* kMultiRingFlxCurrent     =  11   // RHC current flux from systematic param variation of kFlxFinalTuning flux */
   /* /// */
} 
FluxType_t;

} // conv
} // valor
  
#endif 

