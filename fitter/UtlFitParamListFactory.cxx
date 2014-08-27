//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

21-7-2014 Davide
Added samples for multiring forward horn current

*/
//________________________________________________________________________________________

#include <cassert>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>

#include <TSystem.h>

#include "Enumerations.h"
#include "Conventions.h"
#include "Utilities.h"

#include "UtlFitParamListFactory.h"
#include "FitCore/FitParamList.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;

using namespace valor;
using namespace valor::conv;

FitParamList * util::NamedParameterList(string name)
{

  //
  // 3-flavour oscillation fit with/without systematics included 
  //

  if(name == "nuosc3f") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  if(name == "nuosc3f+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "nuosc3f_12constrain+escale") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq,   true);
    int np = fp->GetNumberParameters();
    for(int ip = 0; ip < np; ip++) {
      conv::PhysParamId_t p = fp->WhichParam(ip);
      if(!fp->IsNuisance(p)) continue; 
      if(p != kNuisSuperkEnergyScale &&
	 p != kPhypMixing12 &&
	 p != kPhypDm12sq ) 
	fp -> FixParameter(p);
    }
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  
  //
  // 3-flavour oscillation fit with various oscillation param subsets either fixed or
  // treated as nuisance paramers.
  //

  if(name == "nuosc3f_12constrain") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  if(name == "nuosc3f_12constrain+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "nuosc3f_12fix") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }

  if(name == "nuosc3f_12fix+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
   }

//   /// NEW MULTIRING 2 (to add rhc sample)
//   if(name == "nuosc3f_12fix+nuisance") {
//     FitParamList * fp;  
//     //#ifdef __FHC_RUNNING__
// #if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)
//     #ifdef __MULTIRING_FHC_RUNNING__
//     fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParamsFHC_MultiR);   // fhc all + skdetfsi multiring (no escale)
//     #else
//     fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);  // fhc only 
//     #endif
// #else
//     #ifdef __MULTIRING_FHC_RUNNING__
//     fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParamsMultiR);  // multiring fhc only 
//     #else
//     cerr << "No samples have been chosen for the analysis!!!" << endl;
//     exit(1);    
//     #endif
// #endif // sample analysis option
//     util::AppendAllSystematicsTo(fp, 6);
//     fp -> FixParameter (kPhypMixing12);
//     fp -> FixParameter (kPhypDm12sq  );
//     fp -> SetParListName(name);
//     fp -> SetSystMode(1);
//     return fp;
//    }
//   ///








  if(name == "nuosc3f_numudisapp_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  
  if(name == "nuosc3f_numudisapp_fix4+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypCPV);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  
  if(name == "nuosc3f_numudisapp_nuis3fixdcp+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance (kPhypMixing13, true);
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq, true);
    fp -> FixParameter (kPhypCPV);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }

  if(name == "nuosc3f_nueapp_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm23sq);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  if(name == "nuosc3f_nueapp_fix4+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm23sq);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity_fix5") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDm23sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity_fix3") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(4);//use future osc errors
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(3);//use current osc errors
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity_fix3+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(5);//use future osc errors
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
    util::AppendAllSystematicsTo(fp, 6);
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> SetNuisance  (kPhypMixing13, true);
    fp -> SetNuisance  (kPhypMixing12, true);
    fp -> SetNuisance  (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(6);//use current osc errors
    return fp;
  }
    
    if(name == "nuosc3f_dcpsensitivity_allsyst") {
        FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams);
        util::AppendAllSystematicsTo(fp, 6);
	//        fp -> SetNuisance  (kPhypDm23sq,   true);
	//        fp -> SetNuisance  (kPhypMixing13, true);
        fp -> FixParameter (kPhypMixing12      );
        fp -> FixParameter (kPhypDm12sq        );
        //fp -> FixParameter (kPhypCPV);
        fp -> SetParListName(name);
        fp -> SetSystMode(1);//future oscillation
        return fp;
    }
    
    if(name == "nuosc3f_dcpsensitivity_HK") {
        FitParamList * fp = FitParamList::NuOsc3FlavourFit(kNumOfNuisanceParams + kNumOfParamsInBANFFfit);
        AppendAllSuperkDetectorFSISystematics      (fp, 6);
        AppendAllUncorrelatedXSecSystematics       (fp, 6 + kNumOfSKDetFSIParams+1);
        AppendAllCorrelatedBANFFSystematics        (fp,6 + kNumOfSKDetFSIParams+1 + kNumOfUncXSecSyst);
        AppendAllCorrelatedBANFFSystematicsRHC     (fp,6 + kNumOfSKDetFSIParams+1 + kNumOfUncXSecSyst + kNumOfParamsInBANFFfit);
        fp -> SetNuisance  (kPhypMixing13, true);
        fp -> FixParameter (kPhypMixing12      );
        fp -> FixParameter (kPhypDm12sq        );
        fp -> SetParListName(name);
        fp -> SetSystMode(3);//RHC banff
        return fp;
    }
  
  cerr<<"FitParamListFactory/NamedParamList >> Unable to make: "<<name<<endl;
  exit(-1); 
}

//________________________________________________________________________________________
void util::AppendAllSystematicsTo(FitParamList * fparm, int start)
{
// Add nuisance params
//

  /// NEW MULTIRING 3
#if defined(__NUMUONLY__)

  #if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)  
    #ifdef __MULTIRING_FHC_RUNNING__
    // fhc all + skdetfsi multiring
    util::AppendCorrelatedBANFFSystematicsDisap  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsDisap(fparm, start+kNumOfParamsInBANFFfit);
    util::AppendUncorrelatedXSecSystematicsDisap (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRing+1);
    #else
    // fhc only 
    util::AppendCorrelatedBANFFSystematicsDisap  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsDisap(fparm, start+kNumOfParamsInBANFFfit);
    util::AppendUncorrelatedXSecSystematicsDisap (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams);
    #endif
  #else
    #ifdef __MULTIRING_FHC_RUNNING__
    // multiring fhc only 
    util::AppendCorrelatedBANFFSystematicsDisap  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsDisap(fparm, start+kNumOfParamsInBANFFfit);
    util::AppendUncorrelatedXSecSystematicsDisap (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParamsMultiRing);
    #endif
  #endif // sample analysis option

#elif defined(__NUEONLY__)
  #if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)  
    #ifdef __MULTIRING_FHC_RUNNING__
    // fhc all + skdetfsi multiring
    util::AppendCorrelatedBANFFSystematicsApp  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsApp(fparm, start+kNumOfParamsInBANFFfit);
    util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRing+1);
    #else
    // fhc only 
    util::AppendCorrelatedBANFFSystematicsApp  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsApp(fparm, start+kNumOfParamsInBANFFfit);
    util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams);
    #endif
  #else
    #ifdef __MULTIRING_FHC_RUNNING__
    // multiring fhc only 
    util::AppendCorrelatedBANFFSystematicsApp  (fparm, start);
    util::AppendSuperkDetectorFSISystematicsApp(fparm, start+kNumOfParamsInBANFFfit);
    util::AppendUncorrelatedXSecSystematicsApp (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParamsMultiRing);
    #endif
  #endif // sample analysis option
#else
    ////


    /// NEW MULTIRING 2

  //#ifdef __FHC_RUNNING__  
#if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)  
  #ifdef __MULTIRING_FHC_RUNNING__
    // fhc all + skdetfsi multiring 
  util::AppendAllCorrelatedBANFFSystematics    (fparm, start);
  util::AppendAllSuperkDetectorFSISystematics  (fparm, start+kNumOfParamsInBANFFfit);          
  /// NEW MULTIRING 3
  util::AppendAllUncorrelatedXSecSystematics   (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams+kNumOfSKDetFSIParamsMultiRing+1);  
  ///
  #else
  // fhc only 
  util::AppendAllCorrelatedBANFFSystematics    (fparm, start);
  util::AppendAllSuperkDetectorFSISystematics  (fparm, start+kNumOfParamsInBANFFfit); 
  util::AppendAllUncorrelatedXSecSystematics   (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams);
  #endif
#else
  #ifdef __MULTIRING_FHC_RUNNING__
  // multiring fhc only 
  util::AppendAllCorrelatedBANFFSystematics    (fparm, start);
  util::AppendAllSuperkDetectorFSISystematics  (fparm, start+kNumOfParamsInBANFFfit);  
  util::AppendAllUncorrelatedXSecSystematics   (fparm, start+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParamsMultiRing); 
#endif
#endif // sample analysis option

  ///

  
#endif // oscillation analysis option

}
//________________________________________________________________________________________

//________________________________________________________________________________________
void util::AppendAllUncorrelatedXSecSystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecFermiMom));
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, 0.,               1.,              SYST_STEP,  false, true, AsString(kNuisXSecSpecFunc));
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecWshape));
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecShapeCCoth));
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecPiLessDelta));
#ifndef __NUEONLY__
  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecBindE));
#endif
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormCCcoh));
#ifndef __NUEONLY__
  fparm->AddParam( ip++, kNuisXSecNormNC1piC,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNC1piC));
#endif
#ifndef __NUMUONLY__
  fparm->AddParam( ip++, kNuisXSecNormNCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP*0.3,  false, true, AsString(kNuisXSecNormNCcoh));
#endif
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNCoth));
#ifndef __NUEONLY__
  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormAntiNuToNu));
#endif
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNuEToNuMu));


}
#ifdef __NUMUONLY__
//________________________________________________________________________________________
void util::AppendUncorrelatedXSecSystematicsDisap(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list
  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecFermiMom) );
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, 0.,               1.,              SYST_STEP,  false, true, AsString(kNuisXSecSpecFunc) );
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecWshape) );
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecShapeCCoth) );
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecPiLessDelta) );
  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecBindE) );
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormCCcoh) );
  fparm->AddParam( ip++, kNuisXSecNormNC1piC,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNC1piC) );
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNCoth) );
  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormAntiNuToNu) );
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNuEToNuMu) );
}
#endif
#ifdef __NUEONLY__
//________________________________________________________________________________________
void util::AppendUncorrelatedXSecSystematicsApp(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list
  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecFermiMom) );
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, 0.,               1.,              SYST_STEP,  false, true, AsString(kNuisXSecSpecFunc) );
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecWshape) );
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecShapeCCoth) );
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecPiLessDelta) );
  //  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecBindE) );
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormCCcoh) );
  fparm->AddParam( ip++, kNuisXSecNormNCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNCcoh) );
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNCoth) );
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormNuEToNuMu) );
  //  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisXSecNormAntiNuToNu) );
}
#endif
//________________________________________________________________________________________

void util::AppendAllSuperkDetectorFSISystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  // OLD (UNTIL NEW MULTIRING 2)
  // // SK detector energy scale
  // fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale) );

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get FSI covariance matrix
#if FSI_BINS == 1
  TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,"skdetfsi/covariance_matrix");
#elif FSI_BINS == 2
  TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,"skdetfsi/disap/covariance_matrix");
#elif FSI_BINS == 3
  TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,"skdetfsi/app/covariance_matrix");
#endif
  if(!SKDetFSICovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No SKDET+FSI covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & SKDetFSICov = *SKDetFSICovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;







//   // FSI 
//   for(int i = 0; i < kNumOfSKDetFSIParams; i++) {
// #if FSI_BINS != 3
//     PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + i);
// #else
//     PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + i);
// #endif
//     double err = sqrt(SKDetFSICov[i][i]);
//     double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
//     fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
//   }




  
  ///////////////////////////////////////////////////////////////
  //                                                           //
  // N.B. SAME COVARIANCE MATRIX FOR FHC AND MULTIRING SAMPLES //
  //                                                           //
  ///////////////////////////////////////////////////////////////

  ////////////////////////////////
  //                            //
  // N.B. No RHC sample for now //
  //                            //
  ////////////////////////////////




  /// NEW MULTIRING 2
  
  // FSI 
  //#ifdef __FHC_RUNNING__  
#if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)  
  #ifdef __MULTIRING_FHC_RUNNING__
  // fhc all + skdetfsi multiring
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale) );   // FHC SK escale
  for(int i = 0; i < kNumOfSKDetFSIParams; i++) {
    #if FSI_BINS != 3
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + i);
    #else
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + i);
    #endif
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  //
  fparm->AddParam( ip++, kMultiRingNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kMultiRingNuisSuperkEnergyScale) );   // MultiRing SK escale
  for(int i = 0; i < kNumOfSKDetFSIParamsMultiRing; i++) {
    #if FSI_BINS != 3
    PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI000 + i);
    #else
    PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI006 + i);
    #endif
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  #else
  // fhc only 
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale) );   // FHC SK escale
  for(int i = 0; i < kNumOfSKDetFSIParams; i++) {
    #if FSI_BINS != 3
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + i);
    #else
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + i);
    #endif
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
#endif
#else
  #ifdef __MULTIRING_FHC_RUNNING__
  // multiring fhc only 
  fparm->AddParam( ip++, kMultiRingNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kMultiRingNuisSuperkEnergyScale) );   // MultiRing SK escale
  for(int i = 0; i < kNumOfSKDetFSIParamsMultiRing; i++) {
    #if FSI_BINS != 3
    PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI000 + i);
    #else
    PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI006 + i);
    #endif
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  #endif
#endif // sample analysis option
///

}

#ifdef __NUMUONLY__
//________________________________________________________________________________________
void util::AppendSuperkDetectorFSISystematicsDisap(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list

  // OLD (UNTIL NEW MULTIRING 2)
  // SK detector energy scale
  //fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get FSI covariance matrix
  TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,"skdetfsi/disap/covariance_matrix");
  if(!SKDetFSICovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No SKDET+FSI covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & SKDetFSICov = *SKDetFSICovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;





  // // FSI
  // for(int i = 0; i < kNumOfSKDetFSIParams1Rmu; i++) {
  //   PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + i);
  //   double err = sqrt(SKDetFSICov[i][i]);
  //   double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
  //   fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  // }
  
  ///////////////////////////////////////////////////////////////
  //                                                           //
  // N.B. SAME COVARIANCE MATRIX FOR FHC AND MULTIRING SAMPLES //
  //                                                           //
  ///////////////////////////////////////////////////////////////

  ////////////////////////////////
  //                            //
  // N.B. No RHC sample for now //
  //                            //
  ////////////////////////////////


  /// NEW MULTIRING 2
  
  // FSI
  //#ifdef __FHC_RUNNING__  
#if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)  
  #ifdef __MULTIRING_FHC_RUNNING__
  // fhc all + skdetfsi multiring
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));  // FHC SK escale
  for(int i = 0; i < kNumOfSKDetFSIParams1Rmu; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  //
  fparm->AddParam( ip++, kMultiRingNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kMultiRingNuisSuperkEnergyScale));  // MultiRing SK escale
  for(int i = 0; i < kNumOfSKDetFSIParamsMultiRmu; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI000 + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  #else
  // fhc only
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));  // FHC SK escale
  for(int i = 0; i < kNumOfSKDetFSIParams1Rmu; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  } 
#endif
#else
  #ifdef __MULTIRING_FHC_RUNNING__
  // multiring fhc only 
  fparm->AddParam( ip++, kMultiRingNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kMultiRingNuisSuperkEnergyScale));  // MultiRing SK escale
  for(int i = 0; i < kNumOfSKDetFSIParamsMultiRmu; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI000 + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  #endif
#endif // sample analysis option
  ///



}
#endif // #ifdef __NUMUONLY__

#ifdef __NUEONLY__
//________________________________________________________________________________________
void util::AppendSuperkDetectorFSISystematicsApp(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list
  
  // OLD (UNTIL NEW MULTIRING 2)
  // SK detector energy scale
  //fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get FSI covariance matrix
  TMatrixD * SKDetFSICovMtrx = util::XML2TMatrixD(xml,root_node,"skdetfsi/app/covariance_matrix");
  if(!SKDetFSICovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No SKDET+FSI covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & SKDetFSICov = *SKDetFSICovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;




  // // FSI
  // for(int i = 0; i < kNumOfSKDetFSIParams1Re; i++) {
  //   PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + i);
  //   double err = sqrt(SKDetFSICov[i][i]);
  //   double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
  //   fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  // }
  
  ///////////////////////////////////////////////////////////////
  //                                                           //
  // N.B. SAME COVARIANCE MATRIX FOR FHC AND MULTIRING SAMPLES //
  //                                                           //
  ///////////////////////////////////////////////////////////////

  ////////////////////////////////
  //                            //
  // N.B. No RHC sample for now //
  //                            //
  ////////////////////////////////


  /// NEW MULTIRING 2 
  
  // FSI
  //#ifdef __FHC_RUNNING__  
#if defined(__FHC_RUNNING__) || defined(__RHC_RUNNING__)  
  #ifdef __MULTIRING_FHC_RUNNING__
  // fhc all + skdetfsi multiring
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));  // FHC SK escale
  for(int i = 0; i < kNumOfSKDetFSIParams1Re; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  //
  fparm->AddParam( ip++, kMultiRingNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kMultiRingNuisSuperkEnergyScale));  // MultiRing SK escale
  for(int i = 0; i < kNumOfSKDetFSIParamsMultiRe; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI006 + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  #else
  // fhc only
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kNuisSuperkEnergyScale));  // FHC SK escale
  for(int i = 0; i < kNumOfSKDetFSIParams1Re; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
#endif
#else
  #ifdef __MULTIRING_FHC_RUNNING__
  // multiring fhc only 
  fparm->AddParam( ip++, kMultiRingNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true, AsString(kMultiRingNuisSuperkEnergyScale));  // MultiRing SK escale
  for(int i = 0; i < kNumOfSKDetFSIParamsMultiRe; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kMultiRingNuisSKDetFSI006 + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
  #endif
#endif // sample analysis option
 ///


}
#endif // #ifdef __NUEONLY__

//________________________________________________________________________________________
void util::AppendAllCorrelatedBANFFSystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get BANFF fit central values
#if BANFF_BINS == 1
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/best_fit_values");
#elif BANFF_BINS == 2
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/disap/best_fit_values");
#elif BANFF_BINS == 3
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/app/best_fit_values");
#endif
  if(!BANFFBestFitValues) {
    cerr << "FitParamList/AppendSystematics >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
  const TVectorD & BANFFBestFit = *BANFFBestFitValues;
  // Get BANFF fit covariance matrix
#if BANFF_BINS == 1
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/covariance_matrix");
#elif BANFF_BINS == 2
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/disap/covariance_matrix");
#elif BANFF_BINS == 3
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/app/covariance_matrix");
#endif
  if(!BANFFCovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No BANFF covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & BANFFCov = *BANFFCovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // BANFF 
  for(int i = 0; i < kNumOfParamsInBANFFfit; i++) {
    PhysParamId_t p = kListOfBANFFParams[i];
    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }

}
//_______________________________________________________added by raj_________________________________
void util::AppendAllCorrelatedBANFFSystematicsRHC(FitParamList * fparm, int start)
{
    // Add nuisance params
    //
    
    int ip = start; // first unallocated element in the param list
    
    //    Get the covariance matrices & best fits of FSI & BANFF from xml,
    //    to have something sensible to initialise to
    // Parse the input XML configuration file
    const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
    TXMLEngine* xml = new TXMLEngine;
    XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
    if (xmldoc==0) {
        cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
        delete xml;
        exit(1);
    }
    XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
    // Get BANFF fit central values
#if BANFF_BINS == 1
    TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/best_fit_values");
#elif BANFF_BINS == 2
    TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/disap/best_fit_values");
#elif BANFF_BINS == 3
    TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/app/best_fit_values");
#endif
    if(!BANFFBestFitValues) {
        cerr << "FitParamList/AppendSystematics >> No BANFF best-fit parameter vector was found." << endl;
        exit(1);
    }
    const TVectorD & BANFFBestFit = *BANFFBestFitValues;
    // Get BANFF fit covariance matrix
#if BANFF_BINS == 1
    TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/covariance_matrix");
#elif BANFF_BINS == 2
    TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/disap/covariance_matrix");
#elif BANFF_BINS == 3
    TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/app/covariance_matrix");
#endif
    if(!BANFFCovMtrx) {
        cerr << "FitParamList/AppendSystematics >> No BANFF covariance matrix was found." << endl;
        exit(1);
    }
    const TMatrixD & BANFFCov = *BANFFCovMtrx;
    // Release memory
    //  xml->FreeDoc(xmldoc);
    xml->FreeNode(root_node);
    delete xml;
    
    // BANFF
    for(int i = 0; i < kNumOfParamsInBANFFfit; i++) {
        PhysParamId_t p = PhysParamId_t ((int) kRHCNuisBANFF00 + i);
        double bf = BANFFBestFit[i];
        double err = sqrt(BANFFCov[i][i]);
        double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
        fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true,AsString(p));
    }
    
}
#ifdef __NUMUONLY__
//________________________________________________________________________________________
void util::AppendCorrelatedBANFFSystematicsDisap(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get BANFF fit central values
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/disap/best_fit_values");
  if(!BANFFBestFitValues) {
    cerr << "FitParamList/AppendSystematics >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
  const TVectorD & BANFFBestFit = *BANFFBestFitValues;
  // Get BANFF fit covariance matrix
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/disap/covariance_matrix");
  if(!BANFFCovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No BANFF covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & BANFFCov = *BANFFCovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // BANFF 
  for(int i = 0; i < kNumOfParamsInBANFFfit; i++) {
    PhysParamId_t p = kListOfBANFFParams[i];
    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
}
#endif
#ifdef __NUEONLY__
//________________________________________________________________________________________
void util::AppendCorrelatedBANFFSystematicsApp(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list

  //    Get the covariance matrices & best fits of FSI & BANFF from xml, 
  //    to have something sensible to initialise to
  // Parse the input XML configuration file
  const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");

  TXMLEngine* xml = new TXMLEngine;   
  XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
  if (xmldoc==0) {
    cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
    delete xml;
    exit(1);  
  }
  XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
  // Get BANFF fit central values
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/app/best_fit_values");
  if(!BANFFBestFitValues) {
    cerr << "FitParamList/AppendSystematics >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
  const TVectorD & BANFFBestFit = *BANFFBestFitValues;
  // Get BANFF fit covariance matrix
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/app/covariance_matrix");
  if(!BANFFCovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No BANFF covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & BANFFCov = *BANFFCovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // BANFF 
  for(int i = 0; i < kNumOfParamsInBANFFfit; i++) {
    PhysParamId_t p = kListOfBANFFParams[i];
    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true, AsString(p));
  }
}
#endif
