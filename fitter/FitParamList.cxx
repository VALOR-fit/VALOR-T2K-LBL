//________________________________________________________________________________________
/*
 Changes since v1.0.0:

 @ 05/07/2011 - CA
   Was renamed from PhysParamList to FitParamList so that distinction between the two
   SuperkExpectation interfaces is clear. Added named ctor to return the numu-disappearance
   fit param list including all 34 nuisance params of the 2010 analysis. Added named ctor
   which returns a fit param list by specifying a name. Added boolean array (fIsNuisance) 
   to mark which of the fit params are nuisance params. Added ChisquarePenaltyTerm() 
   method to return the penalty chisq term when nuisance params are included in the fit.
   All nuisance params are de-correlated (see Cholesky)
 @ 15/07/2011 - CA
   Added 3-flavour fit parameter lists and added NamedParameterList(string name) method.
 @ 21/07/2011 - CA
   Added Reset() and WhichParam(int id). Implemented ChisquarePenaltyTerm() for 
   marginalised oscillation parameters.
 @ 29/07/2011 - CA
   Added GetNumberFreeParameters()
 @ 15/09/2011 - CA
   Renamed Reset() to ResetCurrToInitValues() to avoid confusion.
   Added mutable flag fShowCurrValuesOnly and a setter (ShowOnlyCurrValuesInNextPrint())
   method. The flag affects the 'next' print-out only and then it is reset to false.
   Method SetInitValue(PhysParamId_t p, double value) also gets the current value in sync
   with the initial value. Calling TVirtualFitter::Clear() before setting the fit params
   to the fitter in SetTo(TVirtualFitter *). Modified the default, named fit param lists.

*/
//________________________________________________________________________________________

//#define _debug_

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "FitParamList.h"
#include "Utilities.h"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;

using namespace t2kana;
using namespace t2kana::conv;

//________________________________________________________________________________________
namespace t2kana {
  ostream & operator << (ostream & stream, const FitParamList & plst)
  {
    plst.Print(stream);
    return stream;
  }
}
//________________________________________________________________________________________
FitParamList::FitParamList(int nparams)
{
  this->Init(nparams);
}
//________________________________________________________________________________________
FitParamList::~FitParamList()
{
  this->CleanUp();
}
//________________________________________________________________________________________
bool FitParamList::AddParam (int ip, PhysParamId_t p, 
  double init, double min, double max, double step, bool fix, bool nuis, string name)
{
  if(ip < 0 || ip >= fNParams) return false;

  fIdMap.insert(map<PhysParamId_t, int>::value_type(p,ip));

  if(name.size()==0) {
    name = util::AsString(p);
  }

  fCurrValue [ip] = init;
  fInitValue [ip] = init;
  fMinValue  [ip] = min; 
  fMaxValue  [ip] = max; 
  fStepSize  [ip] = step; 
  fIsFixed   [ip] = fix;  
  fIsNuisance[ip] = nuis;  
  fName      [ip] = ( (name == "") ? util::AsString(p) : name );

  return true;
}
//________________________________________________________________________________________
bool FitParamList::Exists(PhysParamId_t p) const
{
  map<PhysParamId_t, int>::const_iterator iter = fIdMap.find(p);
  if(iter != fIdMap.end()) return true;
  return false;
}
//________________________________________________________________________________________
void FitParamList::SetTo(TVirtualFitter * fitter)
{
  if(!fitter) return;
  fitter->Clear();
  map<t2kana::conv::PhysParamId_t, int>::const_iterator iter = fIdMap.begin();
  for ( ; iter != fIdMap.end(); ++iter) {
    PhysParamId_t p = iter->first;
    int iparam = iter->second;
    fitter->SetParameter(
      iparam,
      this->Name(p).c_str(),
      this->InitValue(p),
      this->StepSize(p),
      this->MinValue(p),
      this->MaxValue(p) );

    if(this->IsFixed(p)) {
     fitter->FixParameter(iparam);
    }
  }
}
//________________________________________________________________________________________
void FitParamList::Sync(double * par)
{
  if(!par) return;
  for(int ip = 0; ip < fNParams; ip++) {
     fCurrValue[ip] = par[ip];
  }
}
//________________________________________________________________________________________
void FitParamList::ResetCurrToInitValues(void)
{
  for(int ip = 0; ip < fNParams; ip++) {
     fCurrValue[ip] = fInitValue[ip];
  }
}
//________________________________________________________________________________________
double FitParamList::ChisquarePenaltyTerm(void)
{
  double chisq_penalty = 0;
  for(int ip=0; ip<fNParams;  ip++) {
    // is the current param marked as a nuisance one?
    if( fIsNuisance[ip] ) {
        // what is the actual param at array slot ip?
        PhysParamId_t p = this->WhichParam(ip);
        assert(p != kPhypUndefined);
        // if it is one of the systematics, include the (tweaking dial)^2 in the
        // penalty term; tweaking dials are given in units of the corresponding standard
        // deviation, 0 corresponds to nominal and they are the un-correllated inputs to
        // Cholesky method (wherever correlations are taken into account).
        if( util::IsSystematic(p) ) {
           double twk = fCurrValue[ip];
           chisq_penalty += (twk*twk);
        }
        // if the nuiance param is one of the oscillation params that get marginalised
        // in this particular fit, do one of the following:
        else if(p == kPhypSinsq2Theta23) 
        { 
           // penalize deviations from global analysis best-fit value
           double sinsq2theta23 = fCurrValue[ip];
           double twk = (sinsq2theta23 - kGlobFitSinsq2Theta23) / 
                        (kGlobFitSinsq2Theta23 * kGlobFitErrSinsq2Theta23);
           chisq_penalty += (twk*twk);
        }
        else if(p == kPhypSinsq2Theta13) 
        { 
          // no contribution to \chi^2_{penalty}
        }
        else if(p == kPhypSinsq2Theta12) 
        { 
           // penalize deviations from global analysis best-fit value
           double sinsq2theta12 = fCurrValue[ip];
           double twk = (sinsq2theta12 - kGlobFitSinsq2Theta12) / 
                        (kGlobFitSinsq2Theta12 * kGlobFitErrSinsq2Theta12);
           chisq_penalty += (twk*twk);
        }
        else if(p == kPhypDeltaCP) 
        { 
          // no contribution to \chi^2_{penalty}
        }
        else if(p == kPhypDm23sq) 
        { 
           // penalize deviations from global analysis best-fit value
           double dmsq = fCurrValue[ip];     
           double twk = (dmsq - kGlobFitDm23sq) / (kGlobFitDm23sq * kGlobFitErrDm23sq);
           chisq_penalty += (twk*twk);
        }
        else if(p == kPhypDm12sq) 
        { 
           // penalize deviations from global analysis best-fit value
           double dmsq = fCurrValue[ip];     
           double twk = (dmsq - kGlobFitDm12sq) / (kGlobFitDm12sq * kGlobFitErrDm12sq);
           chisq_penalty += (twk*twk);
        }

    }//nuisance
  }
  return chisq_penalty;
}
//________________________________________________________________________________________
int FitParamList::GetNumberFreeParameters(void)
{
  int nparams = 0;

  for(int ip=0; ip<fNParams;  ip++) {
    bool is_fixed = fIsFixed[ip];
    if(!is_fixed) {
      nparams++;
    }
  }
  return nparams;
}
//________________________________________________________________________________________
void FitParamList::Print(ostream & stream) const
{
  if(fShowCurrValuesOnly) {
     stream << endl
            << setw(6)  << "id"
            << setw(20) << "name"
            << setw(15) << "curr"
            << endl;
     for(int ip=0; ip<fNParams;  ip++) {
       stream << setw(6)  << ip 
              << setw(20) << fName[ip]
              << setprecision(6)
              << setw(15) << fCurrValue[ip]
              << endl;
     }
  } else {
     stream << endl
            << setw(6)  << "id"
            << setw(20) << "name"
            << setw(15) << "curr"
            << setw(15) << "init"
            << setw(15) << "min"
            << setw(15) << "max"
            << setw(8)  << "step"
            << setw(8)  << "fixed?"
            << setw(10) << "nuisance?"
            << endl;
     for(int ip=0; ip<fNParams;  ip++) {
       stream << setw(6)  << ip 
              << setw(20) << fName[ip]
              << setprecision(6)
              << setw(15) << fCurrValue[ip]
              << setw(15) << fInitValue[ip]
              << setw(15) << fMinValue[ip]
              << setw(15) << fMaxValue[ip]
              << setw(8) << fStepSize[ip]
              << setw(8) << ((fIsFixed[ip])    ? " **yes** " : "no")
              << setw(8) << ((fIsNuisance[ip]) ? " **yes** " : "no")
              << endl;
     }
  }

  fShowCurrValuesOnly = false; // has one-off effect
}
//________________________________________________________________________________________
FitParamList * FitParamList::NamedParameterList(string name)
{
  if(name == "nuosc3f") {
    return FitParamList::NuOsc3FlavourFitNoSyst();
  }
  if(name == "nuosc3f_12constrain") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitNoSyst();
    fp -> SetNuisance (kPhypSinsq2Theta12, true);
    fp -> SetNuisance (kPhypDm12sq, true);
  }
  if(name == "nuosc3f_numudisapp_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitNoSyst();
    fp -> FixParameter (kPhypSinsq2Theta13);
    fp -> FixParameter (kPhypSinsq2Theta12);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> FixParameter (kPhypDm12sq);
    return fp;
  }
  if(name == "nuosc3f_nueapp_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitNoSyst();
    fp -> FixParameter (kPhypSinsq2Theta23);
    fp -> FixParameter (kPhypSinsq2Theta12);
    fp -> FixParameter (kPhypDm23sq);
    fp -> FixParameter (kPhypDm12sq);
    return fp;
  }
  if(name == "numudisapp2f") {
    return FitParamList::NumuDisapFitNoSyst();
  }
  if(name == "numudisapp2f+nuisance") {
    return FitParamList::NumuDisapFitAllSyst();
  }
  if(name == "decay") {
    return FitParamList::NumuDecayFitNoSyst();
  }
  if(name == "decoh") {
    return FitParamList::NumuDecohFitNoSyst();
  }

  return 0;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitNoSyst(void)
{
  FitParamList * fparm = new FitParamList(6);

  //                                     init            min     max     step     fix?   nuis?
  fparm->AddParam(0, kPhypSinsq2Theta23, 0.95,           0.50,   1.00,   1.0E-2,  false, false);
  fparm->AddParam(1, kPhypSinsq2Theta13, 0.10,           0.00,   0.20,   1.0E-3,  false, false);
  fparm->AddParam(2, kPhypSinsq2Theta12, kGlobFitDm12sq, 0.70,   1.00,   1.0E-2,  false, false);
  fparm->AddParam(3, kPhypDeltaCP,       kpi,            0.00,   2*kpi,  1.0E-2,  false, false);
  fparm->AddParam(4, kPhypDm23sq,        kGlobFitDm23sq, 1.0E-3, 4.0E-3, 1.0E-4,  false, false);
  fparm->AddParam(5, kPhypDm12sq,        kGlobFitDm12sq, 7.0E-5, 9.0E-5, 2.0E-6,  false, false);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDisapFitNoSyst(void)
{
  FitParamList * fparm = new FitParamList(2);

  //                                     init     min   max    step  fix?   nuis?
  fparm->AddParam(0, kPhypSinsq2Theta23, 0.99,    0.0,  1.0,  0.01,  false, false);
  fparm->AddParam(1, kPhypDm23sq,        2.32E-3, 1E-3, 6E-3, 1E-4,  false, false);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDisapFitAllSyst(void)
{
  FitParamList * fparm = new FitParamList(36);

  //                                          init     min   max    step  fix?   nuis?
  fparm->AddParam( 0, kPhypSinsq2Theta23,     0.99,    0.0,  1.0,  0.01,  false, false);
  fparm->AddParam( 1, kPhypDm23sq,            2.32E-3, 1E-3, 6E-3, 1E-4,  false, false);
  fparm->AddParam( 2, kNuisSuperkEffCCQEoth,  0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam( 3, kNuisSuperkEffCCQErc0u, 0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam( 4, kNuisSuperkEffCCQErc1u, 0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam( 5, kNuisSuperkEffCCQErc2u, 0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam( 6, kNuisSuperkEffCCnonQE,  0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam( 7, kNuisSuperkEffNC,       0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam( 8, kNuisSuperkEffNueCC,    0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam( 9, kNuisSuperkEnergyScale, 0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(10, kNuisNd280Eff,          0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(11, kNuisFluxNorm,          0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(12, kNuisFluxShapeNuMu0u,   0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(13, kNuisFluxShapeNuMu1u,   0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(14, kNuisFluxShapeNuMu2u,   0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(15, kNuisFluxShapeNuMu3u,   0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(16, kNuisFluxShapeNuMu4u,   0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(17, kNuisFluxShapeNuMu5u,   0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(18, kNuisFluxShapeNuMu6u,   0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(19, kNuisFluxShapeNuMu7u,   0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(20, kNuisFluxShapeNuMu8u,   0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(21, kNuisFluxShapeNuMu9u,   0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(22, kNuisFluxShapeNuMu10u,  0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(23, kNuisFluxShapeNuMu11u,  0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(24, kNuisFluxShapeNuMu12u,  0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(25, kNuisFluxShapeNuMu13u,  0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(26, kNuisFluxShapeNuMu14u,  0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(27, kNuisFluxShapeNuMu15u,  0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(28, kNuisFluxShapeNuMu16u,  0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(29, kNuisFluxShapeNuMu17u,  0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(30, kNuisFluxShapeNuMu18u,  0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(31, kNuisXSecCCQE,          0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(32, kNuisXSecCC1piToCCQE,   0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(33, kNuisXSecCCothToCCQE,   0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(34, kNuisXSecNCToCCQE,      0.0,    -1,    1,    1E-1,  false, true);
  fparm->AddParam(35, kNuisXSecNuclEffect,    0.0,    -1,    1,    1E-1,  false, true);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDecayFitNoSyst(void)
{
  FitParamList * fparm = new FitParamList(2);

  //                                     init    min      max    step   fix?   nuis?
  fparm->AddParam(0, kPhypSinsq2Theta23, 1.0,    0.0,     1.0,   0.02,  false, false);
  fparm->AddParam(1, kPhypAlphaDecay,    0.01,   1E-4,   1E-3,   1E-5,  false, false);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDecohFitNoSyst(void)
{
  FitParamList * fparm = new FitParamList(2);

  //                                       init    min     max     step    fix?   nuis?  
  fparm->AddParam(0, kPhypSinsq2Theta23,   1.0,    0.0,    1.0,    0.02,   false, false);
  fparm->AddParam(1, kPhypMusqDecoherence, 0.01,   0.001,  0.04,   0.001,  false, false);

  return fparm;
}
//________________________________________________________________________________________
PhysHypothesis_t FitParamList::PhysicsHypothesis(void) const
{
  // neutrino oscillations (3-flavour framework)
  if( this->Exists(kPhypDm12sq) &&
      this->Exists(kPhypDm23sq) &&
      this->Exists(kPhypSinsq2Theta13) &&
      this->Exists(kPhypSinsq2Theta12) &&
      this->Exists(kPhypSinsq2Theta23) &&
      this->Exists(kPhypDeltaCP) )
  {
    return kHypOsc3F;
  }

  // numu disappearance, 2-flavour neutrino oscillation framework
  else
  if( this->Exists(kPhypDm23sq) &&
      this->Exists(kPhypSinsq2Theta23) )
  {
    return kHypNuMuDisapp2F;
  }
  
  // neutrino decay   
  else
  if( this->Exists(kPhypAlphaDecay) &&
      this->Exists(kPhypSinsq2Theta23) )
  {
    return kHypDecay;
  }
   
  // neutrino decoherence
  else
  if( this->Exists(kPhypMusqDecoherence) &&
      this->Exists(kPhypSinsq2Theta23) )
  {
    return kHypDecoh;
  }

  // undefined
  return kHypUndefined;
}
//________________________________________________________________________________________
void FitParamList::Init(int nparams)
{
   assert(nparams>0);

   fNParams = nparams; 

   fCurrValue  = new double[fNParams];
   fInitValue  = new double[fNParams];
   fMinValue   = new double[fNParams];
   fMaxValue   = new double[fNParams];
   fStepSize   = new double[fNParams];
   fName       = new string[fNParams];
   fIsFixed    = new bool  [fNParams];
   fIsNuisance = new bool  [fNParams];

   fShowCurrValuesOnly = false;
}
//________________________________________________________________________________________
void FitParamList::CleanUp(void)
{
   fIdMap.clear();
   delete [] fCurrValue;
   delete [] fInitValue;
   delete [] fMinValue;
   delete [] fMaxValue;
   delete [] fStepSize;
   delete [] fName;
   delete [] fIsFixed;
   delete [] fIsNuisance;
}
//________________________________________________________________________________________
int FitParamList::Id(PhysParamId_t p) const
{
   map<PhysParamId_t, int>::const_iterator it = fIdMap.find(p);
   if(it == fIdMap.end()) {
     cerr << "Parameter : " << util::AsString(p) 
          << " was not included in current FitParamList" 
         << endl;
     exit(1);
   }
   int id = it->second;
   return id;
}
//________________________________________________________________________________________
t2kana::conv::PhysParamId_t FitParamList::WhichParam(int id) const
{
   map<PhysParamId_t, int>::const_iterator it = fIdMap.begin();
   for( ; it != fIdMap.end(); ++it) {
      int curr_id = it->second;
      if(curr_id == id) return it->first;
   }
   return kPhypUndefined;
}
//________________________________________________________________________________________

