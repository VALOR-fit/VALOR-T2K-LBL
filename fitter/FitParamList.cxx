//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/
//________________________________________________________________________________________

//#define __show_contributions_to_penalty_term__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "FitParamList.h"
#include "Utilities.h"
#include "Opt.h"

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

void FitParamList::SetTheSample (Option_t * sample = ""){
      fSample = sample;
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

  fBestFitValue [ip] = -99999999;
  fCurrValue    [ip] = init;
  fInitValue    [ip] = init;
  fMinValue     [ip] = min; 
  fMaxValue     [ip] = max; 
  fStepSize     [ip] = step; 
  fIsFixed      [ip] = fix;  
  fIsNuisance   [ip] = nuis;  
  fName         [ip] = ( (name == "") ? util::AsString(p) : name );

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
void FitParamList::Update(PhysParamId_t p, double v)
{
  fCurrValue [this->Id(p)] = v;
}
//________________________________________________________________________________________
void FitParamList::SetInitValue(PhysParamId_t p, double v)
{
  fInitValue [this->Id(p)] = v;
  fCurrValue [this->Id(p)] = v;
}                  
//________________________________________________________________________________________
void FitParamList::ResetCurrToInitValues(void)
{
  for(int ip = 0; ip < fNParams; ip++) {
     fCurrValue[ip] = fInitValue[ip];
  }
}
//________________________________________________________________________________________
void FitParamList::CopyCurrToBestFitValues(void)
{
  for(int ip = 0; ip < fNParams; ip++) {
     fBestFitValue[ip] = fCurrValue[ip];
  }
}
//________________________________________________________________________________________
void FitParamList::CopyBestFitToCurrValues(void)
{
  for(int ip = 0; ip < fNParams; ip++) {
     fCurrValue[ip] = fBestFitValue[ip];
  }
}
//________________________________________________________________________________________
void FitParamList::FixParameter(PhysParamId_t p)
{
  fIsFixed [this->Id(p)] = true;
}
//________________________________________________________________________________________
void FitParamList::FreeParameter(PhysParamId_t p)
{
  fIsFixed [this->Id(p)] = false;
}
//________________________________________________________________________________________
void FitParamList::SetNuisance(PhysParamId_t p, bool tf)
{
  fIsNuisance [this->Id(p)] = tf;
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
           double dchisq_penalty = twk*twk;
           chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
           cout << "FitParamList::ChisquarePenaltyTerm >> " 
                << util::AsString(p) << " = " << twk 
                << " ==> deltachi^2_{penalty} = " << dchisq_penalty << endl;
#endif
        }
        // if the nuisance param is one of the oscillation params that get marginalised
        // in this particular fit, do one of the following:
        else if(p == kPhypMixing23) 
        { 
           // penalize deviations from global analysis best-fit value
           double mix = fCurrValue[ip];
           double nom = 0;
           double err = 0;
           NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
           switch(omp) {
              case(kNuMixSinsq2Theta):
                nom = kGlobFitSinsq2Theta23;
                err = kGlobFitSinsq2Theta23 * kGlobFitErrSinsq2Theta23;
                break;
              case(kNuMixSinsqTheta):
                nom = kGlobFitSinsqTheta23;
                err = kGlobFitSinsqTheta23 * kGlobFitErrSinsqTheta23;
                break;
              default:
                exit(1);
           }
           double twk = (mix - nom) / err;
           double dchisq_penalty = twk*twk;
           chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
           cout << "FitParamList::ChisquarePenaltyTerm >> mixing23 = " << mix
                << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
                << " ==> deltachi^2_{penalty} = " << dchisq_penalty << endl;
#endif
        }
        else if(p == kPhypMixing13) 
        { 
          // added contribution to \chi^2_{penalty}
	  // penalize deviations from global analysis best-fit value                                                                                                                  
	  double mix = fCurrValue[ip];
	  double nom = 0;
	  double err = 0;
	  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
	  switch(omp) {
	  case(kNuMixSinsq2Theta):
	    nom = kGlobFitSinsq2Theta13;
	    err = kGlobFitSinsq2Theta13 * kGlobFitErrSinsq2Theta13;
	    break;
	  case(kNuMixSinsqTheta):
	    nom = kGlobFitSinsqTheta13;
	    err = kGlobFitSinsqTheta13 * kGlobFitErrSinsqTheta13;
	    break;
	  default:
	    exit(1);
	  }
	  double twk = (mix - nom) / err;
	  double dchisq_penalty = twk*twk;
	  chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
	  cout << "FitParamList::ChisquarePenaltyTerm >> mixing13 = " << mix
	       << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
	       << " ==> deltachi^2_{penalty} = " << dchisq_penalty << endl;
#endif

        }
        else if(p == kPhypMixing12) 
        { 
           // penalize deviations from global analysis best-fit value
           double mix = fCurrValue[ip];
           double nom = 0;
           double err = 0;
           NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
           switch(omp) {
              case(kNuMixSinsq2Theta):
                  nom = kGlobFitSinsq2Theta12;
                  err = kGlobFitSinsq2Theta12 * kGlobFitErrSinsq2Theta12;
                  break;
              case(kNuMixSinsqTheta):
                  nom = kGlobFitSinsqTheta12;
                  err = kGlobFitSinsqTheta12 * kGlobFitErrSinsqTheta12;
                  break;
              default:
                exit(1);
           }
           double twk = (mix - nom) / err;
           double dchisq_penalty = twk*twk;
           chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
           cout << "FitParamList::ChisquarePenaltyTerm >> mixing12 = " << mix
                << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
                << " ==> deltachi^2_{penalty} = " << dchisq_penalty << endl;
#endif
        }
        else if(p == kPhypDeltaCP) 
        { 
          // added contribution to \chi^2_{penalty}
	  // penalize deviations from global analysis best-fit value                                                                                                                  
	  double dcp  = fCurrValue[ip];
	  double nom  = kGlobFitDeltaCP;
	  double err  = kGlobFitErrDeltaCP;
	  double twk  = (dcp - nom) / err;
	  double dchisq_penalty = twk*twk;
	  chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
	  cout << "FitParamList::ChisquarePenaltyTerm >> dcp = " << dcp
	       << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
	       << " ==> deltachi^2_{penalty} = " << dchisq_penalty << endl;
#endif
        }
        else if(p == kPhypDm12sq)
	  {
	    // penalize deviations from global analysis best-fit value                                                                                                                  
	    double dmsq = fCurrValue[ip];
	    double nom  = kGlobFitDm12sq;
	    double err  = kGlobFitDm12sq * kGlobFitErrDm12sq;
	    double twk  = (dmsq - nom) / err;
	    double dchisq_penalty = twk*twk;
	    chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
	    cout << "FitParamList::ChisquarePenaltyTerm >> dm12sq = " << dmsq
		 << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
		 << " ==> deltachi^2_{penalty} = " << dchisq_penalty << endl;
#endif
      
  }
        else if(p == kPhypDm23sq) 
        { 
           // penalize deviations from global analysis best-fit value
           double dmsq = fCurrValue[ip];    
           double nom  = kGlobFitDm23sq; 
           double err  = kGlobFitDm23sq * kGlobFitErrDm23sq;
           double twk  = (dmsq - nom) / err;
           double dchisq_penalty = twk*twk;
           chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
           cout << "FitParamList::ChisquarePenaltyTerm >> |dm23sq| = " << dmsq
                << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
                << " ==> deltachi^2_{penalty} = " << dchisq_penalty << endl;
#endif
        }
        else if(p == kPhypDm12sq) 
        { 
           // penalize deviations from global analysis best-fit value
           double dmsq = fCurrValue[ip];     
           double nom  = kGlobFitDm12sq;
           double err  = kGlobFitDm12sq * kGlobFitErrDm12sq;
           double twk  = (dmsq - nom) / err;
           double dchisq_penalty = twk*twk;
           chisq_penalty += dchisq_penalty;
#ifdef __show_contributions_to_penalty_term__
           cout << "FitParamList::ChisquarePenaltyTerm >> dm12sq = " << dmsq
                << " (nominal = " << nom << ", error = " << err << ", twk = " << twk << ")"
                << " ==> deltachi^2_{penalty} = " << dchisq_penalty << endl;
#endif
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
  if(fShowBestFitValuesOnly) {
     fShowBestFitValuesOnly = false; // has one-off effect
     stream << endl
            << setw(6)  << "id"
            << setw(20) << "name"
            << setw(15) << "best-fit"
            << endl;
     for(int ip=0; ip<fNParams;  ip++) {
       stream << setw(6)  << ip 
              << setw(20) << fName[ip]
              << setprecision(6)
              << setw(15) << fBestFitValue[ip]
              << endl;
     }
  } 
  else
  if(fShowCurrValuesOnly) {
     fShowCurrValuesOnly = false; // has one-off effect
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
  } 
  else 
  {
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

}
//________________________________________________________________________________________
FitParamList * FitParamList::NamedParameterList(string name)
{
  if(name == "nuosc3f") {
    return FitParamList::NuOsc3FlavourFitNoSyst();
  }
  if(name == "nuosc3f_12constrain") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitNoSyst();
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq,   true);
    return fp;
  }
  if(name == "nuosc3f_12fix") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitNoSyst();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    return fp;
  }
  if(name == "nuosc3f_12fix+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitAllSyst();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    return fp;
  }
  if(name == "nuosc3f_numudisapp_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitNoSyst();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> FixParameter (kPhypDm12sq);
    return fp;
  }
  if(name == "nuosc3f_nueapp_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitNoSyst();
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm23sq);
    fp -> FixParameter (kPhypDm12sq);
    return fp;
  }
  if(name == "nuosc3f_nueapp_fix4+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitNuESyst();
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
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

  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
  case(kNuMixSinsq2Theta):
    fparm->AddParam(0, kPhypMixing23,  kInitFitSinsq2Theta23,  kMinFitSinsq2Theta23,  kMaxFitSinsq2Theta23,  1.0E-2,  false, false);
    fparm->AddParam(1, kPhypMixing13,  kInitFitSinsq2Theta13,  kMinFitSinsq2Theta13,  kMaxFitSinsq2Theta13,  1.0E-3,  false, false);
    fparm->AddParam(2, kPhypMixing12,  kInitFitSinsq2Theta12,  kMinFitSinsq2Theta12,  kMaxFitSinsq2Theta12,  1.0E-2,  false, false);
    break;
  case(kNuMixSinsqTheta):
    fparm->AddParam(0, kPhypMixing23,  kInitFitSinsqTheta23,   kMinFitSinsqTheta23,   kMaxFitSinsqTheta23,   1.0E-2,  false, false);
    fparm->AddParam(1, kPhypMixing13,  kInitFitSinsqTheta13,   kMinFitSinsqTheta13,   kMaxFitSinsqTheta13,   1.0E-3,  false, false);
    fparm->AddParam(2, kPhypMixing12,  kInitFitSinsqTheta12,   kMinFitSinsqTheta12,   kMaxFitSinsqTheta12,   1.0E-2,  false, false);
    break;
  default:
    delete fparm;
    return 0;
  }

  fparm->AddParam(3, kPhypDeltaCP,   0,                -kpi,          kpi,           1.0E-2,  false, false);
  fparm->AddParam(4, kPhypDm23sq,    kInitFitDm23sq,   kMinFitDm23sq, kMaxFitDm23sq, 1.0E-4,  false, false);
  fparm->AddParam(5, kPhypDm12sq,    kInitFitDm12sq,   kMinFitDm12sq, kMaxFitDm12sq, 2.0E-6,  false, false);

  return fparm;
}
//________________________________________________________________________________________

FitParamList * FitParamList::NuOsc3FlavourFitNuESyst(void)
{
      FitParamList * fparm = new FitParamList(27);

      NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
      switch(omp) {
            case(kNuMixSinsq2Theta):
                  fparm->AddParam(0, kPhypMixing23,  kInitFitSinsq2Theta23,  kMinFitSinsq2Theta23,  kMaxFitSinsq2Theta23,  1.0E-2,  false, false);
                  fparm->AddParam(1, kPhypMixing13,  kInitFitSinsq2Theta13,  kMinFitSinsq2Theta13,  kMaxFitSinsq2Theta13,  1.0E-3,  false, false);
                  fparm->AddParam(2, kPhypMixing12,  kInitFitSinsq2Theta12,  kMinFitSinsq2Theta12,  kMaxFitSinsq2Theta12,  1.0E-2,  false, false);
                  break;
            case(kNuMixSinsqTheta):
                  fparm->AddParam(0, kPhypMixing23,  kInitFitSinsqTheta23,   kMinFitSinsqTheta23,   kMaxFitSinsqTheta23,   1.0E-2,  false, false);
                  fparm->AddParam(1, kPhypMixing13,  kInitFitSinsqTheta13,   kMinFitSinsqTheta13,   kMaxFitSinsqTheta13,   1.0E-3,  false, false);
                  fparm->AddParam(2, kPhypMixing12,  kInitFitSinsqTheta12,   kMinFitSinsqTheta12,   kMaxFitSinsqTheta12,   1.0E-2,  false, false);
                  break;
            default:
                  delete fparm;
                  return 0;
      }

      fparm->AddParam(3, kPhypDeltaCP,   0,                -kpi,          kpi,           1.0E-2,  false, false);
      fparm->AddParam(4, kPhypDm23sq,    kInitFitDm23sq,   kMinFitDm23sq, kMaxFitDm23sq, 1.0E-4,  false, false);
      fparm->AddParam(5, kPhypDm12sq,    kInitFitDm12sq,   kMinFitDm12sq, kMaxFitDm12sq, 2.0E-6,  false, false);

      //
      // add nuisance params
      //
/*
      fparm->AddParam( 6, kNuisSuperkEff1RmuCCQEoth,  0.0,           -1,    1,    1E-1,  false, true); 
      fparm->AddParam( 7, kNuisSuperkEff1RmuCCQErc0u, 0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam( 8, kNuisSuperkEff1RmuCCQErc1u, 0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam( 9, kNuisSuperkEff1RmuCCQErc2u, 0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(10, kNuisSuperkEff1RmuCCnonQE,  0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(11, kNuisSuperkEff1RmuNC,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(12, kNuisSuperkEff1RmuNueCC,    0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(13, kNuisSuperkEnergyScale,     0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(14, kNuisNd280Eff,              0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(15, kNuisFluxNorm,              0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(16, kNuisXSecCCQE,              0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(17, kNuisXSecCC1piToCCQE,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(18, kNuisXSecCCothToCCQE,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(19, kNuisXSecNCToCCQE,          0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(20, kNuisXSecNuclEffect,        0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(21,kNuisSuperkEff1ReSig,        0.0,           -1,    1,    1E-1,  false, true);    
      fparm->AddParam(22,kNuisSuperkEff1ReBkg,        0.0,           -1,    1,    1E-1,  false, true);    
      fparm->AddParam(22,kNuisXSecCCcohToCCQE,        0.0,           -1,    1,    1E-1,  false, true);    
      fparm->AddParam(23,kNuisXSecNC1pi0ToCCQE,       0.0,           -1,    1,    1E-1,  false, true);    
      fparm->AddParam(24,kNuisXSecNCcohToCCQE,        0.0,           -1,    1,    1E-1,  false, true);    
      fparm->AddParam(25,kNuisXSecNCothToCCQE,        0.0,           -1,    1,    1E-1,  false, true);    
      fparm->AddParam(26,kNuisXSecNuEToNuMu,          0.0,           -1,    1,    1E-1,  false, true);    
*/

      return fparm;
}
//________________________________________________________________________________________

FitParamList * FitParamList::NuOsc3FlavourFitAllSyst(void)
{
      FitParamList * fparm = new FitParamList(47);

      NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
      switch(omp) {
            case(kNuMixSinsq2Theta):
                  fparm->AddParam(0, kPhypMixing23,  kInitFitSinsq2Theta23,  kMinFitSinsq2Theta23,  kMaxFitSinsq2Theta23,  1.0E-2,  false, false);
                  fparm->AddParam(1, kPhypMixing13,  kInitFitSinsq2Theta13,  kMinFitSinsq2Theta13,  kMaxFitSinsq2Theta13,  1.0E-3,  false, false);
                  fparm->AddParam(2, kPhypMixing12,  kInitFitSinsq2Theta12,  kMinFitSinsq2Theta12,  kMaxFitSinsq2Theta12,  1.0E-2,  false, false);
                  break;
            case(kNuMixSinsqTheta):
                  fparm->AddParam(0, kPhypMixing23,  kInitFitSinsqTheta23,   kMinFitSinsqTheta23,   kMaxFitSinsqTheta23,   1.0E-2,  false, false);
                  fparm->AddParam(1, kPhypMixing13,  kInitFitSinsqTheta13,   kMinFitSinsqTheta13,   kMaxFitSinsqTheta13,   1.0E-3,  false, false);
                  fparm->AddParam(2, kPhypMixing12,  kInitFitSinsqTheta12,   kMinFitSinsqTheta12,   kMaxFitSinsqTheta12,   1.0E-2,  false, false);
                  break;
            default:
                  delete fparm;
                  return 0;
      }

      fparm->AddParam(3, kPhypDeltaCP,   0,                -kpi,          kpi,           1.0E-2,  false, false);
      fparm->AddParam(4, kPhypDm23sq,    kInitFitDm23sq,   kMinFitDm23sq, kMaxFitDm23sq, 1.0E-4,  false, false);
      fparm->AddParam(5, kPhypDm12sq,    kInitFitDm12sq,   kMinFitDm12sq, kMaxFitDm12sq, 2.0E-6,  false, false);

      //
      // add nuisance params
      //
/*
      //numu
      fparm->AddParam( 6, kNuisSuperkEff1RmuCCQEoth,  0.0,           -1,    1,    1E-1,  false, true); 
      fparm->AddParam( 7, kNuisSuperkEff1RmuCCQErc0u, 0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam( 8, kNuisSuperkEff1RmuCCQErc1u, 0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam( 9, kNuisSuperkEff1RmuCCQErc2u, 0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(10, kNuisSuperkEff1RmuCCnonQE,  0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(11, kNuisSuperkEff1RmuNC,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(12, kNuisSuperkEff1RmuNueCC,    0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(13, kNuisSuperkEnergyScale,     0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(14, kNuisNd280Eff,              0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(15, kNuisFluxNorm,              0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(16, kNuisFluxShapeNuMu0u,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(17, kNuisFluxShapeNuMu1u,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(18, kNuisFluxShapeNuMu2u,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(19, kNuisFluxShapeNuMu3u,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(20, kNuisFluxShapeNuMu4u,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(21, kNuisFluxShapeNuMu5u,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(22, kNuisFluxShapeNuMu6u,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(23, kNuisFluxShapeNuMu7u,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(24, kNuisFluxShapeNuMu8u,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(25, kNuisFluxShapeNuMu9u,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(26, kNuisFluxShapeNuMu10u,      0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(27, kNuisFluxShapeNuMu11u,      0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(28, kNuisFluxShapeNuMu12u,      0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(29, kNuisFluxShapeNuMu13u,      0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(30, kNuisFluxShapeNuMu14u,      0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(31, kNuisFluxShapeNuMu15u,      0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(32, kNuisFluxShapeNuMu16u,      0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(33, kNuisFluxShapeNuMu17u,      0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(34, kNuisFluxShapeNuMu18u,      0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(35, kNuisXSecCCQE,              0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(36, kNuisXSecCC1piToCCQE,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(37, kNuisXSecCCothToCCQE,       0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(38, kNuisXSecNCToCCQE,          0.0,           -1,    1,    1E-1,  false, true);
      fparm->AddParam(39, kNuisXSecNuclEffect,        0.0,           -1,    1,    1E-1,  false, true);
      //additional, only important for nue
      fparm->AddParam(40,kNuisSuperkEff1ReSig,        0.0,           -1,    1,    1E-1,  false, true);    
      fparm->AddParam(41,kNuisSuperkEff1ReBkg,        0.0,           -1,    1,    1E-1,  false, true);    
      fparm->AddParam(42,kNuisXSecCCcohToCCQE,        0.0,           -1,    1,    1E-1,  false, true);    
      fparm->AddParam(43,kNuisXSecNC1pi0ToCCQE,       0.0,           -1,    1,    1E-1,  false, true);    
      fparm->AddParam(44,kNuisXSecNCcohToCCQE,        0.0,           -1,    1,    1E-1,  false, true);    
      fparm->AddParam(45,kNuisXSecNCothToCCQE,        0.0,           -1,    1,    1E-1,  false, true);    
      fparm->AddParam(46,kNuisXSecNuEToNuMu,          0.0,           -1,    1,    1E-1,  false, true);    
*/
      return fparm;
}
//________________________________________________________________________________________

FitParamList * FitParamList::NumuDisapFitNoSyst(void)
{
  FitParamList * fparm = new FitParamList(2);

  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
  case(kNuMixSinsq2Theta):
    fparm->AddParam(0, kPhypMixing23, kInitFitSinsq2Theta23, kMinFitSinsq2Theta23,  kMaxFitSinsq2Theta23,  0.01,  false, false);
    break;
  case(kNuMixSinsqTheta):
    fparm->AddParam(0, kPhypMixing23, kInitFitSinsqTheta23,  kMinFitSinsqTheta23,   kMaxFitSinsqTheta23,   0.01,  false, false);
    break;
  default:
    delete fparm;
    return 0;
  }

  fparm->AddParam(1, kPhypDm23sq, kInitFitDm23sq, kMinFitDm23sq, kMaxFitDm23sq, 1E-4, false, false);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDisapFitAllSyst(void)
{
  FitParamList * fparm = new FitParamList(36);

  //
  // add osc params
  //
  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
  case(kNuMixSinsq2Theta):
    fparm->AddParam(0, kPhypMixing23, kInitFitSinsq2Theta23, kMinFitSinsq2Theta23,  kMaxFitSinsq2Theta23,  0.01,  false, false);
    break;
  case(kNuMixSinsqTheta):
    fparm->AddParam(0, kPhypMixing23, kInitFitSinsqTheta23,  kMinFitSinsqTheta23,   kMaxFitSinsqTheta23,   0.01,  false, false);
    break;
  default:
    delete fparm;
    return 0;
  }

  fparm->AddParam(1, kPhypDm23sq, kInitFitDm23sq, kMinFitDm23sq, kMaxFitDm23sq, 1E-4, false, false);

  //
  // add nuisance params
  //
/*
  fparm->AddParam( 2, kNuisSuperkEff1RmuCCQEoth,  0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam( 3, kNuisSuperkEff1RmuCCQErc0u, 0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam( 4, kNuisSuperkEff1RmuCCQErc1u, 0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam( 5, kNuisSuperkEff1RmuCCQErc2u, 0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam( 6, kNuisSuperkEff1RmuCCnonQE,  0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam( 7, kNuisSuperkEff1RmuNC,       0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam( 8, kNuisSuperkEff1RmuNueCC,    0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam( 9, kNuisSuperkEnergyScale,     0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(10, kNuisNd280Eff,              0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(11, kNuisFluxNorm,              0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(12, kNuisFluxShapeNuMu0u,       0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(13, kNuisFluxShapeNuMu1u,       0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(14, kNuisFluxShapeNuMu2u,       0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(15, kNuisFluxShapeNuMu3u,       0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(16, kNuisFluxShapeNuMu4u,       0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(17, kNuisFluxShapeNuMu5u,       0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(18, kNuisFluxShapeNuMu6u,       0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(19, kNuisFluxShapeNuMu7u,       0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(20, kNuisFluxShapeNuMu8u,       0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(21, kNuisFluxShapeNuMu9u,       0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(22, kNuisFluxShapeNuMu10u,      0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(23, kNuisFluxShapeNuMu11u,      0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(24, kNuisFluxShapeNuMu12u,      0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(25, kNuisFluxShapeNuMu13u,      0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(26, kNuisFluxShapeNuMu14u,      0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(27, kNuisFluxShapeNuMu15u,      0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(28, kNuisFluxShapeNuMu16u,      0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(29, kNuisFluxShapeNuMu17u,      0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(30, kNuisFluxShapeNuMu18u,      0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(31, kNuisXSecCCQE,              0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(32, kNuisXSecCC1piToCCQE,       0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(33, kNuisXSecCCothToCCQE,       0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(34, kNuisXSecNCToCCQE,          0.0,           -1,    1,    1E-1,  false, true);
  fparm->AddParam(35, kNuisXSecNuclEffect,        0.0,           -1,    1,    1E-1,  false, true);
*/
  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDecayFitNoSyst(void)
{
  FitParamList * fparm = new FitParamList(2);

  //                                         init    min      max    step   fix?   nuis?
  fparm->AddParam(0, kPhypSinsq2ThetaDecay,  1.0,    0.0,     1.0,   0.02,  false, false);
  fparm->AddParam(1, kPhypAlphaDecay,        0.01,   1E-4,   1E-3,   1E-5,  false, false);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDecohFitNoSyst(void)
{
  FitParamList * fparm = new FitParamList(2);

  //                                          init    min     max     step    fix?   nuis?  
  fparm->AddParam(0, kPhypSinsq2ThetaDecoh,   1.0,    0.0,    1.0,    0.02,   false, false);
  fparm->AddParam(1, kPhypMusqDecoh,          0.01,   0.001,  0.04,   0.001,  false, false);

  return fparm;
}
//________________________________________________________________________________________
PhysHypothesis_t FitParamList::PhysicsHypothesis(void) const
{
  // neutrino oscillations (3-flavour framework)
  if( this->Exists(kPhypDm12sq)   &&
      this->Exists(kPhypDm23sq)   &&
      this->Exists(kPhypMixing13) &&
      this->Exists(kPhypMixing12) &&
      this->Exists(kPhypMixing23) &&
      this->Exists(kPhypDeltaCP) )
  {
    return kHypOsc3F;
  }

  // numu disappearance, 2-flavour neutrino oscillation framework
  else
  if( this->Exists(kPhypDm23sq) &&
      this->Exists(kPhypMixing23) )
  {
    return kHypNuMuDisapp2F;
  }
  
  // neutrino decay   
  else
  if( this->Exists(kPhypAlphaDecay) &&
      this->Exists(kPhypSinsq2ThetaDecay) )
  {
    return kHypDecay;
  }
   
  // neutrino decoherence
  else
  if( this->Exists(kPhypMusqDecoh) &&
      this->Exists(kPhypSinsq2ThetaDecoh) )
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

   fBestFitValue = new double[fNParams];
   fCurrValue    = new double[fNParams];
   fInitValue    = new double[fNParams];
   fMinValue     = new double[fNParams];
   fMaxValue     = new double[fNParams];
   fStepSize     = new double[fNParams];
   fName         = new string[fNParams];
   fIsFixed      = new bool  [fNParams];
   fIsNuisance   = new bool  [fNParams];

   fShowCurrValuesOnly    = false;
   fShowBestFitValuesOnly = false;
}
//________________________________________________________________________________________
void FitParamList::CleanUp(void)
{
   fIdMap.clear();
   delete [] fBestFitValue;
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
     cerr << "FitParamList::Id() >> Parameter: " << util::AsString(p) 
          << " was not included in current fit parameter list" << endl;
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

