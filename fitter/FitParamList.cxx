//________________________________________________________________________________________
/*

 \copyright (c) 2010-2013, VaLOR group.

*/
//________________________________________________________________________________________

#define SYST_NSIGMA 5   // number of sigma away from nominal nuisance parameters are allowed to be
#define SYST_STEP   0.1 // step size (in units of sigma)
#define NSELECTED_NONBANFF 10
#define NSELECTED_NONBANFF_DISAP 6
#define NSELECTED_BANFF    3
#define NSELECTED_BANFF_DISAP    2

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>

#include "PhysParamList.h"
#include "FitParamList.h"
#include "Utilities.h"
#include "Opt.h"

#include <TSystem.h>

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;

using namespace t2kana;
using namespace t2kana::conv;

ClassImp(FitParamList)

//________________________________________________________________________________________
namespace t2kana {
  ostream & operator << (ostream & stream, const FitParamList & plst)
  {
    plst.Print(stream);
    return stream;
  }
}
//________________________________________________________________________________________
FitParamList::FitParamList() :
TObject()
{
  this->Init(0);
}
//________________________________________________________________________________________
FitParamList::FitParamList(int nparams) :
TObject()
{
  this->Init(nparams);
}
//________________________________________________________________________________________
FitParamList::FitParamList(const FitParamList & plst) :
TObject()
{
  this->Init(plst.fNParams);
  this->Copy(plst);
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

  fName         [ip].SetString(name.c_str());
  fBestFitValue [ip] = -99999999;
  fCurrValue    [ip] = init;
  fInitValue    [ip] = init;
  fMinValue     [ip] = min; 
  fMaxValue     [ip] = max; 
  fStepSize     [ip] = step; 
  fIsFixed      [ip] = fix;  
  fIsNuisance   [ip] = nuis;  

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
void FitParamList::CopyBestFitToInitValues(void)
{
  for(int ip = 0; ip < fNParams; ip++) {
     fInitValue[ip] = fBestFitValue[ip];
  }
}
//________________________________________________________________________________________
void FitParamList::CopyBestFitToInitValues(const FitParamList & plst)
{
  if(plst.GetParListName() != fParListName) {
    cerr << "Cannot copy BestFit values in FitParamList of type " << plst.GetParListName()
	 << " to InitValues in FitParamList of type " << fParListName << endl;
    exit(1);
  }
  for(int ip = 0; ip < fNParams; ip++) {
    fInitValue[ip] = plst.BestFitValue(ip);
  }
}
//________________________________________________________________________________________
void FitParamList::SetErrors(double * errplus, double * errminus, double * errmigrad)
{
  for(int ip = 0; ip < fNParams; ip++) {
     fPlusError  [ip] = errplus  [ip];
     fMinusError [ip] = errminus [ip];
     fMigradError[ip] = errmigrad[ip];
     if(errplus[ip] > kASmallNum || errminus[ip] > kASmallNum)
       fMinosErrorsSet = true;
  }
  fMigradErrorsSet = true;
}
//________________________________________________________________________________________
void FitParamList::SetMigradErrors(vector<double> err)
{
  if(err.size() != (unsigned int)this->Size()) {
    cerr << "FitterKernelMinuit2/SetMigradErrors >> Error vector size: " << err.size()
	 << " is not the same size as the number of free parameters: "
	 << this->GetNumberFreeParameters() << endl;
    exit(1);
  }
  // vector only contains free parameters, so need to make sure we handle the
  // fixed parameters in the parameter list correctly
  /*
  int ifree = 0;
  for(int ip = 0; ip < fNParams; ip++) {
    if(fIsFixed[ip]) {
      fMigradError [ip] = 0;
    }
    else {
      fMigradError[ip] = err[ifree];
      ifree++;
    }
  }
  */
  for(int ip = 0; ip < fNParams; ip++) {
    fMigradError [ip] = err[ip];
  }
  fMigradErrorsSet = true;
}
//________________________________________________________________________________________
void FitParamList::SetHesseErrors(vector<double> err)
{
  if(err.size() != (unsigned int)this->Size()) {
    cerr << "FitterKernelMinuit2/SetHesseErrors >> Error vector size: " << err.size()
	 << " is not the same size as the number of free parameters: "
	 << this->GetNumberFreeParameters() << endl;
    exit(1);
  }
  /*
  // vector only contains free parameters, so need to make sure we handle the
  // fixed parameters in the parameter list correctly
  int ifree = 0;
  for(int ip = 0; ip < fNParams; ip++) {
    if(fIsFixed[ip]) {
      fHesseError [ip] = 0;
    }
    else {
      fHesseError[ip] = err[ifree];
      ifree++;
    }
  }
  */
  for(int ip = 0; ip < fNParams; ip++) {
    fHesseError [ip] = err[ip];
  }
  fHesseErrorsSet = true;
}
//________________________________________________________________________________________
void FitParamList::SetMinosErrors(vector<double> plus, vector<double> minus)
{
  if((int)plus.size() != this->GetNumberFreeParameters()) {
    cerr << "FitterKernelMinuit2/SetMinosErrors >> Error vector size: " << plus.size()
	 << " is not the same size as the number of free parameters: "
	 << this->GetNumberFreeParameters() << endl;
    exit(1);
  }
  if(plus.size() != minus.size()) {
    cerr << "FitterKernelMinuit2/SetMinosErrors >> Error vector size (plus): " << plus.size()
	 << " is not the same size as error vector size (minus): " << minus.size() << endl;
    exit(1);
  }

  // vector only contains free parameters, so need to make sure we handle the
  // fixed parameters in the parameter list correctly
  int ifree = 0;
  for(int ip = 0; ip < fNParams; ip++) {
    if(fIsFixed[ip]) {
      fPlusError   [ip] = 0;
      fMinusError  [ip] = 0;
    }
    else {
      fPlusError [ip] = plus [ifree];
      fMinusError[ip] = minus[ifree];
      ifree++;
    }
  }
  /*
  for(int ip = 0; ip < fNParams; ip++) {
    fPlusError [ip] = plus [ip];
    fMinusError[ip] = minus[ip];
  }
  */
  fMinosErrorsSet = true;
}
//________________________________________________________________________________________
void FitParamList::FixParameter(PhysParamId_t p)
{
  fIsFixed [this->Id(p)] = true;
}
//________________________________________________________________________________________
void FitParamList::FixParameter(int ip)
{
  fIsFixed [ip] = true;
}
//________________________________________________________________________________________
void FitParamList::FreeParameter(PhysParamId_t p)
{
  fIsFixed [this->Id(p)] = false;
}
//________________________________________________________________________________________
void FitParamList::FreeParameter(int ip)
{
  fIsFixed [ip] = false;
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
void FitParamList::SetTo(ROOT::Minuit2::MnUserParameters & params)
{
  //if(!params) return;
  //params.Clear();
  for(int ip = 0; ip < fNParams; ip++) {
    t2kana::conv::PhysParamId_t p = WhichParam(ip);
    //keep mixing angles in physical region
    if(util::IsMixingAngle(p)) {
      params.Add(this->Name(p).c_str(), this->InitValue(p), this->StepSize(p), 0, 1);
    }
    /*
    //stop splines going too far from nominal
    //(normalisations will eventually make poisson term grow too, not necessarily true for splines)
    else if(util::ArrIdx::Instance()->XSecSystIdx(p) >= 0) {
      if(kXSecSystHasWghtSplines[util::ArrIdx::Instance()->XSecSystIdx(p)])
	params.Add(this->Name(p).c_str(), this->InitValue(p), this->StepSize(p), this->MinValue(p), this->MaxValue(p));
      else
	params.Add(this->Name(p).c_str(), this->InitValue(p), this->StepSize(p));
    }
    //sk energy scale is also special
    if(p == kNuisSuperkEnergyScale)
	params.Add(this->Name(p).c_str(), this->InitValue(p), this->StepSize(p), this->MinValue(p), this->MaxValue(p));
    */
    // confine everything else
    else {
      params.Add(this->Name(p).c_str(), this->InitValue(p), this->StepSize(p), this->MinValue(p), this->MaxValue(p));
    }
    /*
    //everything else is free
    else {
      params.Add(this->Name(p).c_str(), this->InitValue(p), this->StepSize(p));
    }
    */
    if(this->IsFixed(p)) {
     params.Fix(ip);
    }
  }//loop over params
}
//________________________________________________________________________________________
void FitParamList::SetTo(ROOT::Math::Minimizer * fitter)
{
  if(!fitter) return;
  fitter->Clear();
  map<t2kana::conv::PhysParamId_t, int>::const_iterator iter = fIdMap.begin();
  // for ( ; iter != fIdMap.end(); ++iter) {
  for(int ip = 0; ip < fNParams; ip++) {  
    //    PhysParamId_t p = iter->first;
    PhysParamId_t p = WhichParam(ip);
    //    int iparam = iter->second;
    int iparam = ip;

    if(this->IsFixed(p)) {
      fitter->SetFixedVariable(
       iparam,
       this->Name(p).c_str(), 
       this->InitValue(p));
    }
    else{
      fitter->SetLimitedVariable(
	  iparam, 
	  this->Name(p).c_str(), 
	  this->InitValue(p), 
	  this->StepSize(p),
	  this->MinValue(p),
	  this->MaxValue(p) );
    }
  }
}
//________________________________________________________________________________________
void FitParamList::SetTo(ROOT::Fit::Fitter * fitter)
{
  if(!fitter) return;
  map<t2kana::conv::PhysParamId_t, int>::const_iterator iter = fIdMap.begin();
  // for ( ; iter != fIdMap.end(); ++iter) {
  double par0[fNParams];
  for(int ip = 0; ip < fNParams; ip++) {
    PhysParamId_t p = WhichParam(ip);
    par0[ip] = this->InitValue(p);
  }
  (fitter->Config()).SetParamsSettings(fNParams,par0);

  for(int ip = 0; ip < fNParams; ip++) {  
    //    PhysParamId_t p = iter->first;
    PhysParamId_t p = WhichParam(ip);
    //    int iparam = iter->second;
    int iparam = ip;

    if(this->IsFixed(p)) {
      (fitter->Config()).ParSettings(iparam).Set(this->Name(p).c_str(),this->InitValue(p));
    }
    else{
      (fitter->Config()).ParSettings(iparam).Set(this->Name(p).c_str(),this->InitValue(p),this->StepSize(p),this->MinValue(p),this->MaxValue(p) );
    }
  }
}
//________________________________________________________________________________________
void FitParamList::Sync(std::vector<double> par)
{
  if(par.size() <= 0) return;
  for(int ip = 0; ip < fNParams; ip++) {
     fCurrValue[ip] = par[ip];
  }
}
//________________________________________________________________________________________
void FitParamList::Sync(const double * par)
{
  if(!par) return;
  for(int ip = 0; ip < fNParams; ip++) {
     fCurrValue[ip] = par[ip];
  }
}
//________________________________________________________________________________________
PhysParamList FitParamList::AsPhysParamList(Option_t * option) const
{
  PhysParamList pl;
  TString opt(option);
  map<PhysParamId_t, int>::const_iterator iter = fIdMap.begin();
  for( ; iter != fIdMap.end(); ++iter)
  {
    PhysParamId_t p = iter->first;
    int ip = iter->second;
    double v = 0;
    if      ( opt.Contains("initial")  ) v = fInitValue[ip];
    else if ( opt.Contains("current")  ) v = fCurrValue[ip];
    else if ( opt.Contains("best-fit") ) v = fBestFitValue[ip];
    pl.SetValue(p,v);
  }
  return pl;
}
//________________________________________________________________________________________
void FitParamList::FromPhysParamList(t2kana::PhysParamList pl, Option_t * option, bool skip_osc_params)
{
  TString opt(option);
  map<PhysParamId_t, int>::const_iterator iter = fIdMap.begin();
  for( ; iter != fIdMap.end(); ++iter)
  {
    PhysParamId_t p = iter->first;
    if(skip_osc_params && !util::IsSystematic(p))
      continue;
    int ip = iter->second;
    double input = pl.Value(p);
    if(input > -99998) {
      if      ( opt.Contains("initial")  ) fInitValue[ip]    = input;
      else if ( opt.Contains("current")  ) fCurrValue[ip]    = input;
      else if ( opt.Contains("best-fit") ) fBestFitValue[ip] = input;
    }
  }
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
int FitParamList::GetNumberParameters(void)
{
  int nparams = 0;

  for(int ip=0; ip<fNParams;  ip++) {
    nparams++;
  }
  return nparams;
}
//________________________________________________________________________________________
void FitParamList::ClearStatusCodes(void)
{
  fStatusCodes.clear();
}
//________________________________________________________________________________________
void FitParamList::SetStatusCode(string st, int code)
{
  fStatusCodes[st] = code;
}
//________________________________________________________________________________________
int FitParamList::GetStatusCode(string st) const
{
  map<string,int>::const_iterator sciter = fStatusCodes.find(st);
  if(sciter != fStatusCodes.end()) {
    return sciter->second;
  }
  return -9876;
}
//________________________________________________________________________________________
void FitParamList::Print(ostream & stream) const
{
  if(fShowBestFitValuesOnly) {
     fShowBestFitValuesOnly = false; // has one-off effect
     stream << endl
            << setw(6)  << "id"
            << setw(80) << "name"
            << setw(15) << "best-fit"
            << endl;
     for(int ip=0; ip<fNParams;  ip++) {
       stream << setw(6)  << ip 
              << setw(80) << fName[ip].GetString().Data()
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
            << setw(80) << "name"
            << setw(15) << "curr"
            << endl;
     for(int ip=0; ip<fNParams;  ip++) {
       stream << setw(6)  << ip 
              << setw(80) << fName[ip].GetString().Data()
              << setprecision(6)
              << setw(15) << fCurrValue[ip]
              << endl;
     }
  } 
  else 
  {
     stream << endl
            << setw(6)  << "id"
            << setw(80) << "name"
            << setw(15) << "curr"
            << setw(15) << "init"
            << setw(15) << "min"
            << setw(15) << "max"
            << setw(15) << "step"
            << setw(8)  << "fixed?"
            << setw(10) << "nuisance?"
            << endl;
     for(int ip=0; ip<fNParams;  ip++) {
       stream << setw(6)  << ip 
              << setw(80) << fName[ip].GetString().Data()
              << setprecision(6)
              << setw(15) << fCurrValue[ip]
              << setw(15) << fInitValue[ip]
              << setw(15) << fMinValue[ip]
              << setw(15) << fMaxValue[ip]
              << setw(15) << fStepSize[ip]
              << setw(8) << ((fIsFixed[ip])    ? "yes" : "no")
              << setw(8) << ((fIsNuisance[ip]) ? "yes" : "no")
              << endl;
     }
  }

  stream << "Status codes: " << endl;  
  map<string,int>::const_iterator sciter = fStatusCodes.begin();
  for( ; sciter != fStatusCodes.end(); ++sciter) {
    string name = sciter->first;
    int code = sciter->second;
    stream << setw(20) << name << setw(10) << code << endl;  
  }
  stream << "Fit time:  " << fFitTime  << " sec" << endl;
  stream << "Fit tries: " << fFitTries << " tries" << endl;
  stream << "Best-fit chisq: " << fBestFitChisq << endl;  
}
//________________________________________________________________________________________
void FitParamList::Print() const
{
  this->Print(cout);
}
//________________________________________________________________________________________
FitParamList * FitParamList::NamedParameterList(string name)
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
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSyst();
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "nuosc3f+xsec") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithXSec();
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "nuosc3f+nuisanceselection") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithSelectedSyst();
    fp -> SetParListName(name);
    fp -> SetSystMode(2);
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
  if(name == "nuosc3f_12fix") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  if(name == "nuosc3f_12fix+escale") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithEScale();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "nuosc3f_12fix+xsec") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithXSec();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "nuosc3f_12fix+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSyst();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "nuosc3f_12constrain+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSyst();
    fp -> SetNuisance (kPhypMixing12, true);
    fp -> SetNuisance (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "nuosc3f_12fix+nuisance_woskdetfsi") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSystButSKDetFSI();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "nuosc3f_12fix+nuisanceselection") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithSelectedSyst();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "nuosc3f_numudisapp_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  if(name == "nuosc3f_numudisapp_fix4+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSystDisap();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "nuosc3f_numudisapp_fix4+nuisanceselection") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithSelectedSyst();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(2);
    return fp;
  }
  if(name == "nuosc3f_numudisapp_nuis4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetNuisance  (kPhypMixing13, true);
    fp -> SetNuisance  (kPhypMixing12, true);
    fp -> SetNuisance  (kPhypDeltaCP,  true);
    fp -> SetNuisance  (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(3);//special value (get's chisq for osc params, but not for systematics)
    return fp;
  }
  if(name == "nuosc3f_numudisapp_nuis4+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSystDisap();
    fp -> SetNuisance  (kPhypMixing13, true);
    fp -> SetNuisance  (kPhypMixing12, true);
    fp -> SetNuisance  (kPhypDeltaCP,  true);
    fp -> SetNuisance  (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "nuosc3f_numudisapp_nuis4+nuisanceselection") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithSelectedSyst();
    fp -> SetNuisance  (kPhypMixing13, true);
    fp -> SetNuisance  (kPhypMixing12, true);
    fp -> SetNuisance  (kPhypDeltaCP,  true);
    fp -> SetNuisance  (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(2);
    return fp;
  }
  if(name == "nuosc3f_numudisapp_nuis3fixdcp+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSystDisap();
    fp -> SetNuisance  (kPhypMixing13, true);
    fp -> SetNuisance  (kPhypMixing12, true);
    fp -> SetNuisance  (kPhypDm12sq,   true);
    fp -> FixParameter (kPhypDeltaCP);
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
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSystApp();
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm23sq);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }

  /*
  if(name == "nuosc3f_numudisapp_fix4+nuisanceselectionnobanff") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithSelectedSystNoBANFF();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode();
    return fp;
  }
  if(name == "nuosc3f_numudisapp_fix4+nuisancefullbanffselection") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithSelectedSystFullBANFF();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode();
    return fp;
  }
  */

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
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> SetNuisance  (kPhypMixing13, true);
    fp -> SetNuisance  (kPhypMixing12, true);
    fp -> SetNuisance  (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(3);//use current osc errors
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity_fix3+nuisance") {
    //FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSyst();
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithHackSyst();
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    /*
    fp -> FixParameter (kNuisXSecFermiMom);
    fp -> FixParameter (kNuisXSecSpecFunc);
    fp -> FixParameter (kNuisXSecBindE);
    */
    fp -> SetParListName(name);
    fp -> SetSystMode(5);//use future osc errors
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity+nuisance") {
    //FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSyst();
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithHackSyst();
    fp -> SetNuisance  (kPhypMixing23, true);
    fp -> SetNuisance  (kPhypDm23sq,   true);
    fp -> SetNuisance  (kPhypMixing13, true);
    fp -> SetNuisance  (kPhypMixing12, true);
    fp -> SetNuisance  (kPhypDm12sq,   true);
    fp -> SetParListName(name);
    fp -> SetSystMode(6);//use current osc errors
    return fp;
  }

  if(name == "nuosc3f_dcpsensitivity_13error")
    {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetNuisance  (kPhypMixing13, true);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> SetParListName(name);
    fp -> SetSystMode(4);
    return fp;
    }
  //FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSystApp();

  //
  // 2-flavour numu-disappearance oscillation fit with/without systematics included 
  //
  if(name == "numudisapp2f") {
    FitParamList * fp = FitParamList::NumuDisapFit();
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }
  if(name == "numudisapp2f+nuisance") {
    FitParamList * fp = FitParamList::NumuDisapFitWithAllSyst();
    fp -> SetParListName(name);
    fp -> SetSystMode(1);
    return fp;
  }
  if(name == "numudisapp2f+nuisanceselection") {
    FitParamList * fp = FitParamList::NumuDisapFitWithSelectedSyst();
    fp -> SetParListName(name);
    fp -> SetSystMode(2);
    return fp;
  }

  //
  // numu decay fit
  //
  if(name == "decay") {
    FitParamList * fp = FitParamList::NumuDecayFit();
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }

  //
  // numu decoherence fit
  //
  if(name == "decoh") {
    FitParamList * fp = FitParamList::NumuDecohFit();
    fp -> SetParListName(name);
    fp -> SetSystMode(0);
    return fp;
  }

  return 0;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFit(void)
{
  FitParamList * fparm = new FitParamList(6);

  FitParamList::Append3FlavourOscillation(fparm, 0);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithAllSyst(void)
{
 FitParamList * fparm = new FitParamList(6+kNumOfNuisanceParams); // NOTE: This will have to change for the 2012a 1Re and joint analysis

  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendAllCorrelatedBANFFSystematics (fparm, 6);
#ifndef SKFSITOGETHER
  FitParamList::AppendAllSuperkDetectorSystematics  (fparm, 6+kNumOfParamsInBANFFfit);
  FitParamList::AppendAllFSISystematics             (fparm, 6+kNumOfParamsInBANFFfit+kNumOfSKDetParams);
  FitParamList::AppendAllUncorrelatedXSecSystematics(fparm, 6+kNumOfParamsInBANFFfit+kNumOfSKDetParams+kNumOfFSIParams);
#endif
#ifdef SKFSITOGETHER
  FitParamList::AppendAllSuperkDetectorFSISystematics  (fparm, 6+kNumOfParamsInBANFFfit);
  FitParamList::AppendAllUncorrelatedXSecSystematics(fparm, 6+kNumOfParamsInBANFFfit+1+kNumOfSKDetFSIParams);
#endif
  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithHackSyst(void)
{
  FitParamList * fparm = new FitParamList(6+4);

  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  int ip = 6;
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormNC1piC,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormNCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithEScale(void)
{
 FitParamList * fparm = new FitParamList(6+1); // NOTE: This will have to change for the 2012a 1Re and joint analysis

  FitParamList::Append3FlavourOscillation(fparm, 0);
  fparm->AddParam( 6, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithXSec(void)
{
 FitParamList * fparm = new FitParamList(6+12); // NOTE: This will have to change for the 2012a 1Re and joint analysis

  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendAllUncorrelatedXSecSystematics(fparm, 6);
  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithAllSystButSKDetFSI(void)
{
 FitParamList * fparm = new FitParamList(6+kNumOfNuisanceParams-kNumOfSKDetFSIParams-1); // NOTE: This will have to change for the 2012a 1Re and joint analysis

  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendAllCorrelatedBANFFSystematics (fparm, 6);
  FitParamList::AppendAllUncorrelatedXSecSystematics(fparm, 6+kNumOfParamsInBANFFfit);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithAllSystDisap(void)
{
  //only way to work correctly with stand-alone modes is to have defined
  // __NUMUONLY__ or __NUEONLY__, to have correct matrices, bins, etc.
  // But if that is defined, number of params should already be correct for mode

  FitParamList * fparm = new FitParamList(6+kNumOfNuisanceParams); // NOTE: This will have to change for the 2012a 1Re and joint analysis
  FitParamList::Append3FlavourOscillation(fparm, 0);
#ifdef __NUMUONLY__  
  //
  // add nuisance params
  //
  FitParamList::AppendAllCorrelatedBANFFSystematics (fparm, 6);
#ifndef SKFSITOGETHER
  FitParamList::AppendAllSuperkDetectorSystematics  (fparm, 6 +kNumOfParamsInBANFFfit);
  FitParamList::AppendAllFSISystematics             (fparm, 6 +kNumOfParamsInBANFFfit+kNumOfSKDetParams1Rmu);
  FitParamList::AppendAllUncorrelatedXSecSystematics(fparm, 6 +kNumOfParamsInBANFFfit+kNumOfSKDetParams1Rmu+kNumOfFSIParams1Rmu);
#endif
#ifdef SKFSITOGETHER
  FitParamList::AppendAllSuperkDetectorFSISystematics(fparm, 6 +kNumOfParamsInBANFFfit);
  FitParamList::AppendAllUncorrelatedXSecSystematics (fparm, 6 +kNumOfParamsInBANFFfit+kNumOfSKDetFSIParams + 1 /*energy scale*/);
#endif
  //old method: not necessary, now previous methods work with stand-alone as well,
  //number of parameters in each category is already well defined depending on analysis mode
  //  FitParamList::AppendCorrelatedBANFFSystematicsDisap (fparm, 6);
  //#ifndef SKFSITOGETHER
  //FitParamList::AppendSuperkDetectorSystematicsDisap  (fparm, 6 +kNumOfParamsInBANFFfit);
  //FitParamList::AppendFSISystematicsDisap             (fparm, 6 +kNumOfParamsInBANFFfit+kNumOfSKDetParams1Rmu);
  //FitParamList::AppendUncorrelatedXSecSystematicsDisap(fparm, 6 +kNumOfParamsInBANFFfit+kNumOfSKDetParams1Rmu+kNumOfFSIParams1Rmu);
  //#endif
  //#ifdef SKFSITOGETHER
  //FitParamList::AppendSuperkDetectorFSISystematicsDisap(fparm, 6 +kNumOfParamsInBANFFfit);
  //FitParamList::AppendUncorrelatedXSecSystematicsDisap (fparm, 6 +kNumOfParamsInBANFFfit+kNumOfSKDetFSIParams1Rmu + 1 /*energy scale*/);
  //#endif
#endif
  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithAllSystApp(void)
{
  //only way to work correctly with stand-alone modes is to have defined                    
  // __NUMUONLY__ or __NUEONLY__, to have correct matrices, bins, etc.                                                                                          
  // But if that is defined, number of params should already be correct for mode                                                                             

  FitParamList * fparm = new FitParamList(6+kNumOfNuisanceParams); // NOTE: This will have to change for the 2012a 1Re and joint analysis
#if FSI_BINS == 3
  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendAllCorrelatedBANFFSystematics   (fparm, 6);
#ifndef SKFSITOGETHER
  FitParamList::AppendAllSuperkDetectorSystematics    (fparm, 6+kNumOfParamsInBANFFfit);
  FitParamList::AppendAllFSISystematics               (fparm, 6+kNumOfParamsInBANFFfit+1+kNumOfSKDetParams1Re);
  FitParamList::AppendAllUncorrelatedXSecSystematics  (fparm, 6+kNumOfParamsInBANFFfit+1+kNumOfSKDetParams1Re+kNumOfFSIParams1Re);
#endif
#ifdef SKFSITOGETHER
  FitParamList::AppendAllSuperkDetectorFSISystematics    (fparm, 6+kNumOfParamsInBANFFfit);
  FitParamList::AppendAllUncorrelatedXSecSystematics  (fparm, 6+kNumOfParamsInBANFFfit+kNumOfSKDetFSIParams1Re+1);
#endif
  //old method: not necessary, now previous methods work with stand-alone as well,           
  //number of parameters in each category is already well defined depending on analysis mode 
  //#ifndef SKFSITOGETHER
  //FitParamList::AppendSuperkDetectorSystematicsApp    (fparm, 6+kNumOfParamsInBANFFfit);
  //FitParamList::AppendFSISystematicsApp               (fparm, 6+kNumOfParamsInBANFFfit+1+kNumOfSKDetParams1Re);
  //FitParamList::AppendUncorrelatedXSecSystematicsApp  (fparm, 6+kNumOfParamsInBANFFfit+1+kNumOfSKDetParams1Re+kNumOfFSIParams1Re);
  //#endif
  //#ifdef SKFSITOGETHER
  //FitParamList::AppendSuperkDetectorFSISystematicsApp    (fparm, 6+kNumOfParamsInBANFFfit);
  //FitParamList::AppendUncorrelatedXSecSystematicsApp  (fparm, 6+kNumOfParamsInBANFFfit+kNumOfSKDetFSIParams1Re+1);
  //#endif
#endif
  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithSelectedSyst(void)
{
  FitParamList * fparm = 0;
  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaCombined)
    fparm = new FitParamList(6+NSELECTED_NONBANFF+NSELECTED_BANFF); // NOTE: This will have to change for the 2012a 1Re and joint analysis
  else 
    fparm = new FitParamList(6+NSELECTED_NONBANFF_DISAP+NSELECTED_BANFF_DISAP);

  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendSelectedNonBANFFSystematics       (fparm, 6);
  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaCombined){
    FitParamList::AppendSelectedCorrelatedBANFFSystematics(fparm, 6+NSELECTED_NONBANFF);
  }
  else FitParamList::AppendSelectedCorrelatedBANFFSystematics(fparm, 6+NSELECTED_NONBANFF_DISAP);
  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithSelectedSystFullBANFF(void)
{
  FitParamList * fparm = 0 ;
  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaCombined)
    fparm = new FitParamList(6+NSELECTED_NONBANFF+kNumOfParamsInBANFFfit); // NOTE: This will have to change for the 2012a 1Re and joint analysis
  else
    fparm = new FitParamList(6+NSELECTED_NONBANFF_DISAP+kNumOfParamsInBANFFfit); // NOTE: This will have to change for the 2012a 1Re and joint analysis

  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendSelectedNonBANFFSystematics  (fparm, 6);
  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaCombined){
    FitParamList::AppendAllCorrelatedBANFFSystematics(fparm, 6+NSELECTED_NONBANFF);
  }
  else   FitParamList::AppendAllCorrelatedBANFFSystematics(fparm, 6+NSELECTED_NONBANFF_DISAP);
  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithSelectedSystNoBANFF(void)
{
  FitParamList * fparm = 0 ;
  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaCombined)
   fparm = new FitParamList(6+NSELECTED_NONBANFF); // NOTE: This will have to change for the 2012a 1Re and joint analysis
  else
   fparm = new FitParamList(6+NSELECTED_NONBANFF_DISAP); // NOTE: This will have to change for the 2012a 1Re and joint analysis

  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendSelectedNonBANFFSystematics         (fparm, 6);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDisapFit(void)
{
  FitParamList * fparm = new FitParamList(2);

  FitParamList::Append2FlavourOscillation(fparm, 0);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDisapFitWithAllSyst(void)
{

  FitParamList * fparm = new FitParamList(2+kNumOfNuisanceParams);
#ifdef __NUMUONLY__
  FitParamList::Append2FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendAllCorrelatedBANFFSystematics (fparm, 2);
#ifndef SKFSITOGETHER
  FitParamList::AppendAllSuperkDetectorSystematics  (fparm, 2+kNumOfParamsInBANFFfit);
  FitParamList::AppendAllFSISystematics             (fparm, 2+kNumOfParamsInBANFFfit+kNumOfSKDetParams1Rmu);
  FitParamList::AppendAllUncorrelatedXSecSystematics(fparm, 2+kNumOfParamsInBANFFfit+kNumOfSKDetParams1Rmu+kNumOfFSIParams1Rmu);
#endif
#ifdef SKFSITOGETHER
  FitParamList::AppendAllSuperkDetectorFSISystematics(fparm, 2+kNumOfParamsInBANFFfit);
  FitParamList::AppendAllUncorrelatedXSecSystematics (fparm, 2+kNumOfParamsInBANFFfit+kNumOfSKDetFSIParams1Rmu+ 1/*energy scale*/);
#endif
  //old method: not necessary, now previous methods work with stand-alone as well,           
  //number of parameters in each category is already well defined depending on analysis mode 
  //  FitParamList::AppendCorrelatedBANFFSystematicsDisap (fparm, 2);
  //#ifndef SKFSITOGETHER
  //FitParamList::AppendSuperkDetectorSystematicsDisap  (fparm, 2+kNumOfParamsInBANFFfit);
  //FitParamList::AppendFSISystematicsDisap             (fparm, 2+kNumOfParamsInBANFFfit+kNumOfSKDetParams1Rmu);
  //FitParamList::AppendUncorrelatedXSecSystematicsDisap(fparm, 2+kNumOfParamsInBANFFfit+kNumOfSKDetParams1Rmu+kNumOfFSIParams1Rmu);
  //#endif
  //#ifdef SKFSITOGETHER
  //FitParamList::AppendSuperkDetectorFSISystematicsDisap(fparm, 2+kNumOfParamsInBANFFfit);
  //FitParamList::AppendUncorrelatedXSecSystematicsDisap (fparm, 2+kNumOfParamsInBANFFfit+kNumOfSKDetFSIParams1Rmu+ 1/*energy scale*/);
  //#endif

#endif
  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDisapFitWithSelectedSyst(void)
{
  FitParamList * fparm = new FitParamList(2+NSELECTED_NONBANFF_DISAP+NSELECTED_BANFF_DISAP);

  FitParamList::Append2FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendSelectedNonBANFFSystematics       (fparm, 2);
  FitParamList::AppendSelectedCorrelatedBANFFSystematics(fparm, 2+NSELECTED_NONBANFF_DISAP);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDecayFit(void)
{
  FitParamList * fparm = new FitParamList(2);

  //                                         init    min      max    step   fix?   nuis?
  fparm->AddParam(0, kPhypSinsq2ThetaDecay,  1.0,    0.0,     1.0,   0.02,  false, false);
  fparm->AddParam(1, kPhypAlphaDecay,        0.01,   1E-4,   1E-3,   1E-5,  false, false);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDecohFit(void)
{
  FitParamList * fparm = new FitParamList(2);

  //                                          init    min     max     step    fix?   nuis?  
  fparm->AddParam(0, kPhypSinsq2ThetaDecoh,   1.0,    0.0,    1.0,    0.02,   false, false);
  fparm->AddParam(1, kPhypMusqDecoh,          0.01,   0.001,  0.04,   0.001,  false, false);

  return fparm;
}
//________________________________________________________________________________________
void FitParamList::Append3FlavourOscillation(FitParamList * fparm, int start)
{

  int ip = start; // first unallocated element in the param list

  //
  // add osc params
  //
  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
  case(kNuMixSinsq2Theta):
     fparm->AddParam( ip++, kPhypMixing23,  kInitFitSinsq2Theta23,  kMinFitSinsq2Theta23,  kMaxFitSinsq2Theta23,  kStepSinsq2Theta23,  false, false);
     fparm->AddParam( ip++, kPhypMixing13,  kInitFitSinsq2Theta13,  kMinFitSinsq2Theta13,  kMaxFitSinsq2Theta13,  kStepSinsq2Theta13,  false, false);
     fparm->AddParam( ip++, kPhypMixing12,  kInitFitSinsq2Theta12,  kMinFitSinsq2Theta12,  kMaxFitSinsq2Theta12,  kStepSinsq2Theta12,  false, false);
     break;
  case(kNuMixSinsqTheta):
     fparm->AddParam( ip++, kPhypMixing23,  kInitFitSinsqTheta23,   kMinFitSinsqTheta23,   kMaxFitSinsqTheta23,   kStepSinsqTheta23,  false, false);
     fparm->AddParam( ip++, kPhypMixing13,  kInitFitSinsqTheta13,   kMinFitSinsqTheta13,   kMaxFitSinsqTheta13,   kStepSinsqTheta13,  false, false);
     fparm->AddParam( ip++, kPhypMixing12,  kInitFitSinsqTheta12,   kMinFitSinsqTheta12,   kMaxFitSinsqTheta12,   kStepSinsqTheta12,  false, false);
     break;
  default:
     delete fparm;
     return;
  }
  
  NuOscDeltaCP_t dcp = Opt::Instance()->gOscDeltaCP;
  switch(dcp){
  case(kNuDeltaCPPhase):
  fparm->AddParam( ip++, kPhypDeltaCP,   0,                -kpi,          kpi,           1.0E-2,  false, false);
  break;
  case(kNuDeltaCPSin):
  fparm->AddParam( ip++, kPhypDeltaCP,   0,                  -1,            1,           1.0E-2,  false, false);
  break;
  default:
    delete fparm;
    return;
  }

  //  fparm->AddParam( ip++, kPhypDm23sq,    kInitFitDm23sq,   kMinFitDm23sq, kMaxFitDm23sq, 1.0E-4,  false, false);
  NuOscDmsqScheme_t odmsq = Opt::Instance()->gOscDmsqScheme;
  switch (odmsq) {
  case(kNuDmsqFogliLisi):
    fparm->AddParam( ip++, kPhypDm23sq,    kGlobFitDmsqFL,   kMinFitDm23sq, kMaxFitDm23sq, kStepDm23sq,  false, false);
    break;
  case(kNuDmsq12andAbs23):
    fparm->AddParam( ip++, kPhypDm23sq,    kInitFitDm23sq,   kMinFitDm23sq, kMaxFitDm23sq, kStepDm23sq,  false, false);
    break;
  case(kNuDmsq12andAbs31):
    if(Opt::Instance()->gMassHierarchyAssumed == kNmhNormal)
      fparm->AddParam( ip++, kPhypDm23sq,    kInitFitDm23sq + kInitFitDm12sq,   kMinFitDm23sq, kMaxFitDm23sq, kStepDm23sq,  false, false);
    else if(Opt::Instance()->gMassHierarchyAssumed == kNmhInverted)
      fparm->AddParam( ip++, kPhypDm23sq,    kInitFitDm23sq - kInitFitDm12sq,   kMinFitDm23sq, kMaxFitDm23sq, kStepDm23sq,  false, false);
    break;
  default:
    delete fparm;
    return;
  }
  
  fparm->AddParam( ip++, kPhypDm12sq,    kInitFitDm12sq,   kMinFitDm12sq, kMaxFitDm12sq, kStepDm12sq,  false, false);

}
//________________________________________________________________________________________
void FitParamList::Append2FlavourOscillation(FitParamList * fparm, int start)
{

  int ip = start; // first unallocated element in the param list

  //
  // add osc params
  //
  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
  case(kNuMixSinsq2Theta):
    fparm->AddParam( ip++, kPhypMixing23, kInitFitSinsq2Theta23, kMinFitSinsq2Theta23,  kMaxFitSinsq2Theta23,  kStepSinsq2Theta23,  false, false);
    break;
  case(kNuMixSinsqTheta):
    fparm->AddParam( ip++, kPhypMixing23, kInitFitSinsqTheta23,  kMinFitSinsqTheta23,   kMaxFitSinsqTheta23,   kStepSinsq2Theta23,  false, false);
    break;
  default:
    delete fparm;
    return;
  }
  fparm->AddParam( ip++, kPhypDm23sq, kInitFitDm23sq, kMinFitDm23sq, kMaxFitDm23sq, kStepDm23sq, false, false);

}
//________________________________________________________________________________________
void FitParamList::AppendAllUncorrelatedXSecSystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, 0.,               1.,              SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
#ifndef __NUEONLY__
  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
#endif
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
#ifndef __NUEONLY__
  fparm->AddParam( ip++, kNuisXSecNormNC1piC,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
#endif
#ifndef __NUMUONLY__
  fparm->AddParam( ip++, kNuisXSecNormNCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP*0.3,  false, true);
#endif
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
#ifndef __NUEONLY__
  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
#endif
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);


}
#ifdef __NUMUONLY__
//________________________________________________________________________________________
void FitParamList::AppendUncorrelatedXSecSystematicsDisap(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list
  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, 0.,               1.,              SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormNC1piC,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
}
#endif
#ifdef __NUEONLY__
//________________________________________________________________________________________
void FitParamList::AppendUncorrelatedXSecSystematicsApp(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list
  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, 0.,               1.,              SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  //  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormNCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  //  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
}
#endif
//________________________________________________________________________________________
void FitParamList::AppendAllSuperkDetectorSystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);

#if FSI_BINS != 3
  // SK detector
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCQEoth,  0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCQErc0u, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCQErc1u, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCQErc2u, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCnonQE,  0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuNC,       0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuNueCC,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
#endif

#if FSI_BINS != 2
  // Adding the part corresponding to 1Re efficiencies (with covariance matrix) 

  //    Get the covariance matrices & best fits of SKDet
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
  TMatrixD * SKDet1ReCovMtrx = util::XML2TMatrixD(xml,root_node,"detector_systematics/superk/efficiency_relative_error/one_e_ring");
  if(!SKDet1ReCovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No SKDet 1Re eff covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & SKDet1ReCov = *SKDet1ReCovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // SKDet 1Re eff
  for(int i = 0; i < kNumOfSKDetParams1Re; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSuperkEff1ReNueSigCC0 + i);
    double err = sqrt(SKDet1ReCov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }
#endif

}
#ifdef __NUMUONLY__
//________________________________________________________________________________________
void FitParamList::AppendSuperkDetectorSystematicsDisap(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list

  // SK detector
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCQEoth,  0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCQErc0u, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCQErc1u, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCQErc2u, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCnonQE,  0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuNC,       0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuNueCC,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
}
#endif
#ifdef __NUEONLY__
//________________________________________________________________________________________
void FitParamList::AppendSuperkDetectorSystematicsApp(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list

  // SK detector
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);

  // Adding the part corresponding to 1Re efficiencies (with covariance matrix) 

  //    Get the covariance matrices & best fits of SKDet
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
  TMatrixD * SKDet1ReCovMtrx = util::XML2TMatrixD(xml,root_node,"detector_systematics/superk/efficiency_relative_error/one_e_ring");
  if(!SKDet1ReCovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No SKDet 1Re eff covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & SKDet1ReCov = *SKDet1ReCovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // SKDet 1Re eff
  for(int i = 0; i < kNumOfSKDetParams1Re; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSuperkEff1ReNueSigCC0 + i);
    double err = sqrt(SKDet1ReCov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }
}
#endif 
//________________________________________________________________________________________
void FitParamList::AppendAllSuperkDetectorFSISystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list
  // SK detector energy scale
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);

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

  // FSI
  for(int i = 0; i < kNumOfSKDetFSIParams; i++) {
#if FSI_BINS != 3
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + i);
#else
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + i);
#endif
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }
}
#ifdef __NUMUONLY__
//________________________________________________________________________________________
void FitParamList::AppendSuperkDetectorFSISystematicsDisap(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list
  // SK detector energy scale
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);

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

  // FSI
  for(int i = 0; i < kNumOfSKDetFSIParams1Rmu; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI000 + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }
}
#endif
#ifdef __NUEONLY__
//________________________________________________________________________________________
void FitParamList::AppendSuperkDetectorFSISystematicsApp(FitParamList * fparm, int start)
{
// Add nuisance params
//
  int ip = start; // first unallocated element in the param list
  // SK detector energy scale
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);

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

  // FSI
  for(int i = 0; i < kNumOfSKDetFSIParams1Re; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisSKDetFSI006 + i);
    double err = sqrt(SKDetFSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }
}
#endif
//________________________________________________________________________________________
void FitParamList::AppendAllFSISystematics(FitParamList * fparm, int start)
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
  // Get FSI covariance matrix
#if FSI_BINS == 1
  TMatrixD * FSICovMtrx = util::XML2TMatrixD(xml,root_node,"fsi/covariance_matrix");
#elif FSI_BINS == 2
  TMatrixD * FSICovMtrx = util::XML2TMatrixD(xml,root_node,"fsi/disap/covariance_matrix");
#elif FSI_BINS == 3
  TMatrixD * FSICovMtrx = util::XML2TMatrixD(xml,root_node,"fsi/app/covariance_matrix");
#endif
  if(!FSICovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No FSI covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & FSICov = *FSICovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // FSI
  for(int i = 0; i < kNumOfFSIParams; i++) {
#if FSI_BINS != 3
    PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + i);
#else
    PhysParamId_t p = PhysParamId_t ((int) kNuisFSI006 + i);
#endif
    double err = sqrt(FSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }
}
#ifdef __NUMUONLY__
//________________________________________________________________________________________
void FitParamList::AppendFSISystematicsDisap(FitParamList * fparm, int start)
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
  // Get FSI covariance matrix
  TMatrixD * FSICovMtrx = util::XML2TMatrixD(xml,root_node,"fsi/disap/covariance_matrix");
  if(!FSICovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No FSI covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & FSICov = *FSICovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // FSI
  for(int i = 0; i < kNumOfFSIParams1Rmu; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + i);
    double err = sqrt(FSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }
}
#endif
#ifdef __NUEONLY__
//________________________________________________________________________________________
void FitParamList::AppendFSISystematicsApp(FitParamList * fparm, int start)
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
  // Get FSI covariance matrix
  TMatrixD * FSICovMtrx = util::XML2TMatrixD(xml,root_node,"fsi/app/covariance_matrix");
  if(!FSICovMtrx) {
    cerr << "FitParamList/AppendSystematics >> No FSI covariance matrix was found." << endl;
    exit(1);
  }
  const TMatrixD & FSICov = *FSICovMtrx;
  // Release memory
  //  xml->FreeDoc(xmldoc);
  xml->FreeNode(root_node);
  delete xml;

  // FSI
  for(int i = 0; i < kNumOfFSIParams1Re; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisFSI006 + i);
    double err = sqrt(FSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }

}
#endif
//________________________________________________________________________________________
void FitParamList::AppendAllCorrelatedBANFFSystematics(FitParamList * fparm, int start)
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
#elif BANFF_BINS == 4
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/disapmec/best_fit_values");
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
#elif BANFF_BINS == 4
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/disapmec/covariance_matrix");
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
    PhysParamId_t p = PhysParamId_t ((int) kNuisBANFF00 + i);
    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }

}
#ifdef __NUMUONLY__
//________________________________________________________________________________________
void FitParamList::AppendCorrelatedBANFFSystematicsDisap(FitParamList * fparm, int start)
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
    PhysParamId_t p = PhysParamId_t ((int) kNuisBANFF00 + i);
    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }
}
#endif
#ifdef __NUEONLY__
//________________________________________________________________________________________
void FitParamList::AppendCorrelatedBANFFSystematicsApp(FitParamList * fparm, int start)
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
    PhysParamId_t p = PhysParamId_t ((int) kNuisBANFF00 + i);
    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }
}
#endif
//________________________________________________________________________________________
void FitParamList::AppendSelectedNonBANFFSystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  // XSec
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, 0.,               1.,              SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaCombined){
    fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true); //new - joint
  }
#if FSI_BINS != 3
  // SK detector
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCnonQE,  0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuNC,       0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
#endif
  //new - joint
#if FSI_BINS == 1
  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaCombined){
    const char * config_file = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), "analysis_configuration.xml");
    TXMLEngine* xml = new TXMLEngine;
    XMLDocPointer_t xmldoc = xml->ParseFile(config_file);
    if (xmldoc==0) {
      cerr << "Exp >> " << "Can not parse input: " << config_file << endl;
      delete xml;
      exit(1);
    }
    XMLNodePointer_t root_node = xml->DocGetRootElement(xmldoc);
    TMatrixD * SKDet1ReCovMtrx = util::XML2TMatrixD(xml,root_node,"detector_systematics/superk/efficiency_relative_error/one_e_ring");
    if(!SKDet1ReCovMtrx) {
      cerr << "FitParamList/AppendSystematics >> No SKDet 1Re eff covariance matrix was found." << endl;
      exit(1);
    }
    const TMatrixD & SKDet1ReCov = *SKDet1ReCovMtrx;
    // Release memory 
    xml->FreeNode(root_node);
    delete xml;
    for(int i = 0; i < kNumOfSKDetParams1Re; i++) {
      PhysParamId_t p = PhysParamId_t ((int) kNuisSuperkEff1ReNueSigCC0 + i);
      double err = sqrt(SKDet1ReCov[i][i]);
      double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
      if ((p == kNuisSuperkEff1ReNueSigCC2) || (p == kNuisSuperkEff1ReNueCC2) || (p == kNuisSuperkEff1ReNC2) ) 
	fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
    }
  }
#endif

}
//________________________________________________________________________________________
void FitParamList::AppendSelectedCorrelatedBANFFSystematics(FitParamList * fparm, int start)
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
  TVectorD * BANFFBestFitValues = util::XML2TVectorD(xml,root_node,"banff/best_fit_values");
  if(!BANFFBestFitValues) {
    cerr << "FitParamList/AppendSystematics >> No BANFF best-fit parameter vector was found." << endl;
    exit(1);
  }
  const TVectorD & BANFFBestFit = *BANFFBestFitValues;
  // Get BANFF fit covariance matrix
  TMatrixD * BANFFCovMtrx = util::XML2TMatrixD(xml,root_node,"banff/covariance_matrix");
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
  int bidx = util::ArrIdx::Instance()->BANFFSystIdx(kNuisBANFFMAQE);
  double bf  = BANFFBestFit[bidx];
  double err = sqrt(BANFFCov[bidx][bidx]);
  double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
  fparm->AddParam( ip++, kNuisBANFFMAQE,    bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  bidx = util::ArrIdx::Instance()->BANFFSystIdx(kNuisBANFFCCQEE0);
  bf  = BANFFBestFit[bidx];
  err = sqrt(BANFFCov[bidx][bidx]);
  lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
  fparm->AddParam( ip++, kNuisBANFFCCQEE0,  bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  //new - joint
  if (Opt::Instance()->gOscAnaTypeChoice == conv::kAnaCombined){
    bidx = util::ArrIdx::Instance()->BANFFSystIdx(kNuisBANFF01);
    bf  = BANFFBestFit[bidx];
    err = sqrt(BANFFCov[bidx][bidx]);
    lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    fparm->AddParam( ip++, kNuisBANFF01,  bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }
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
void FitParamList::Copy(const FitParamList & plst)
{
  if(this->fNParams != plst.fNParams)
  {
     this->CleanUp();
     this->Init(plst.fNParams);
  }

  this->fParListName = plst.fParListName;
  this->fSystMode    = plst.fSystMode;

  this->fIdMap.insert(plst.fIdMap.begin(),  plst.fIdMap.end());

  for(int i = 0; i < plst.fNParams; i++) {
     this->fName[i].SetString(plst.fName[i].GetString().Data());
     this->fBestFitValue[i] = plst.fBestFitValue[i];
     this->fCurrValue[i]    = plst.fCurrValue[i];
     this->fInitValue[i]    = plst.fInitValue[i];
     this->fMinValue[i]     = plst.fMinValue[i];
     this->fMaxValue[i]     = plst.fMaxValue[i];
     this->fStepSize[i]     = plst.fStepSize[i];
     this->fIsFixed[i]      = plst.fIsFixed[i];
     this->fIsNuisance[i]   = plst.fIsNuisance[i];
     this->fMigradError[i]  = plst.fMigradError[i];
     this->fHesseError[i]   = plst.fHesseError[i];
     this->fPlusError[i]    = plst.fPlusError[i];
     this->fMinusError[i]   = plst.fMinusError[i];
   } 

   this->fBestFitChisq = plst.fBestFitChisq;
   this->fFitTime      = plst.fFitTime;
   this->fFitTries     = plst.fFitTries;

   this->fStatusCodes.insert(plst.fStatusCodes.begin(),  plst.fStatusCodes.end());

   this->fMigradErrorsSet = plst.fMigradErrorsSet;
   this->fHesseErrorsSet  = plst.fHesseErrorsSet;
   this->fMinosErrorsSet  = plst.fMinosErrorsSet;

   this->fShowCurrValuesOnly    = plst.fShowCurrValuesOnly;
   this->fShowBestFitValuesOnly = plst.fShowBestFitValuesOnly;
}
//________________________________________________________________________________________
void FitParamList::Init(int nparams)
{
   fIdMap.clear();

   fParListName = "";
   fSystMode = -1;

   fNParams = nparams; 

   if(nparams>0) {
     fName         = new TObjString[fNParams];
     fBestFitValue = new double[fNParams];
     fCurrValue    = new double[fNParams];
     fInitValue    = new double[fNParams];
     fMinValue     = new double[fNParams];
     fMaxValue     = new double[fNParams];
     fStepSize     = new double[fNParams];
     fIsFixed      = new bool  [fNParams];
     fIsNuisance   = new bool  [fNParams];
     fMigradError  = new double[fNParams];
     fHesseError   = new double[fNParams];
     fPlusError    = new double[fNParams];
     fMinusError   = new double[fNParams];
   } else {
     fName         = 0;
     fBestFitValue = 0;
     fCurrValue    = 0;
     fInitValue    = 0;
     fMinValue     = 0;
     fMaxValue     = 0;
     fStepSize     = 0;
     fIsFixed      = 0;
     fIsNuisance   = 0;
     fMigradError  = 0;
     fHesseError   = 0;
     fPlusError    = 0;
     fMinusError   = 0;
   }

   fShowCurrValuesOnly    = false;
   fShowBestFitValuesOnly = false;

   fBestFitChisq = 0.;
   fFitTime      = -1;
   fFitTries     = 0;

   fStatusCodes.clear();

   fMigradErrorsSet = false;
   fHesseErrorsSet  = false;
   fMinosErrorsSet  = false;
}
//________________________________________________________________________________________
void FitParamList::CleanUp(void)
{
   fIdMap.clear();
   fStatusCodes.clear();

   delete [] fName;
   delete [] fBestFitValue;
   delete [] fCurrValue;
   delete [] fInitValue;
   delete [] fMinValue;
   delete [] fMaxValue;
   delete [] fStepSize;
   delete [] fIsFixed;
   delete [] fIsNuisance;
   delete [] fMigradError;
   delete [] fHesseError;
   delete [] fPlusError;
   delete [] fMinusError;

   fName         = 0;
   fBestFitValue = 0;
   fCurrValue    = 0;
   fInitValue    = 0;
   fMinValue     = 0;
   fMaxValue     = 0;
   fStepSize     = 0;
   fIsFixed      = 0;
   fIsNuisance   = 0;
   fMigradError  = 0;
   fHesseError   = 0;
   fPlusError    = 0;
   fMinusError   = 0;

   fBestFitChisq = 0.;
   fFitTime      = 0.;
   fFitTries     = -1;

   fParListName = "";
   fSystMode = -1;
}
//________________________________________________________________________________________
int FitParamList::Id(PhysParamId_t p) const
{
   map<PhysParamId_t, int>::const_iterator it = fIdMap.find(p);
   if(it == fIdMap.end()) {
     cerr << "FitParamList::Id() >> Parameter: " << util::AsString(p) << " PhysParamId_t p = " << p 
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

