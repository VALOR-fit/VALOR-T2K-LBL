//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/
//________________________________________________________________________________________

#define SYST_NSIGMA 5   // number of sigma away from nominal nuisance parameters are allowed to be
#define SYST_STEP   0.1 // step size (in units of sigma)
#define NSELECTED_NONBANFF 5
#define NSELECTED_BANFF    2
#define NSELECTED_UNBANFF  0

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <math.h>

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
void FitParamList::Sync(std::vector<double> par)
{
  if(par.size() <= 0) return;
  for(int ip = 0; ip < fNParams; ip++) {
     fCurrValue[ip] = par[ip];
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
            << setw(50) << "name"
            << setw(15) << "best-fit"
            << endl;
     for(int ip=0; ip<fNParams;  ip++) {
       stream << setw(6)  << ip 
              << setw(50) << fName[ip]
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
            << setw(50) << "name"
            << setw(15) << "curr"
            << endl;
     for(int ip=0; ip<fNParams;  ip++) {
       stream << setw(6)  << ip 
              << setw(50) << fName[ip]
              << setprecision(6)
              << setw(15) << fCurrValue[ip]
              << endl;
     }
  } 
  else 
  {
     stream << endl
            << setw(6)  << "id"
            << setw(50) << "name"
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
              << setw(50) << fName[ip]
              << setprecision(6)
              << setw(15) << fCurrValue[ip]
              << setw(15) << fInitValue[ip]
              << setw(15) << fMinValue[ip]
              << setw(15) << fMaxValue[ip]
              << setw(15) << fStepSize[ip]
              << setw(8) << ((fIsFixed[ip])    ? " **yes** " : "no")
              << setw(8) << ((fIsNuisance[ip]) ? " **yes** " : "no")
              << endl;
     }
  }

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
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSyst();
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
  if(name == "nuosc3f_numudisapp_fix4+nuisanceunbanff") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSystUnBANFF();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode();
    return fp;
  }
  if(name == "nuosc3f_numudisapp_fix4+nuisanceunbanffselection") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithSelectedSystUnBANFF();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode();
    return fp;
  }
  if(name == "nuosc3f_numudisapp_fix4+nuisanceunbanfffullbanffselection") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithSelectedSystFullBANFFUnBANFF();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    fp -> SetSystMode();
    return fp;
  }
  */
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
  FitParamList * fparm = new FitParamList(6+kNumOfNuisanceParams1Rmu); // NOTE: This will have to change for the 2012a 1Re and joint analysis

  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendAllCorrelatedBANFFSystematics (fparm, 6);
  FitParamList::AppendAllSuperkDetectorSystematics  (fparm, 6+kNumOfParamsInBANFFfit);
  FitParamList::AppendAllFSISystematics             (fparm, 6+kNumOfParamsInBANFFfit+kNumOfSKDetParams);
  FitParamList::AppendAllUncorrelatedXSecSystematics(fparm, 6+kNumOfParamsInBANFFfit+kNumOfSKDetParams+kNumOfFSIParams);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithSelectedSyst(void)
{
  FitParamList * fparm = new FitParamList(6+NSELECTED_NONBANFF+NSELECTED_BANFF); // NOTE: This will have to change for the 2012a 1Re and joint analysis

  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendSelectedNonBANFFSystematics       (fparm, 6);
  FitParamList::AppendSelectedCorrelatedBANFFSystematics(fparm, 6+NSELECTED_NONBANFF);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithSelectedSystFullBANFF(void)
{
  FitParamList * fparm = new FitParamList(6+NSELECTED_NONBANFF+kNumOfParamsInBANFFfit); // NOTE: This will have to change for the 2012a 1Re and joint analysis

  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendSelectedNonBANFFSystematics  (fparm, 6);
  FitParamList::AppendAllCorrelatedBANFFSystematics(fparm, 6+NSELECTED_NONBANFF);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithAllSystUnBANFF(void)
{
  FitParamList * fparm = new FitParamList(6+kNumOfNuisanceParams1Rmu); // NOTE: This will have to change for the 2012a 1Re and joint analysis

  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendAllUncorrelatedBANFFSystematics(fparm, 6);
  FitParamList::AppendAllSuperkDetectorSystematics   (fparm, 6+kNumOfParamsInBANFFfit);
  FitParamList::AppendAllFSISystematics              (fparm, 6+kNumOfParamsInBANFFfit+kNumOfSKDetParams);
  FitParamList::AppendAllUncorrelatedXSecSystematics (fparm, 6+kNumOfParamsInBANFFfit+kNumOfSKDetParams+kNumOfFSIParams);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithSelectedSystUnBANFF(void)
{
  FitParamList * fparm = new FitParamList(6+NSELECTED_NONBANFF+NSELECTED_UNBANFF); // NOTE: This will have to change for the 2012a 1Re and joint analysis

  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendSelectedNonBANFFSystematics         (fparm, 6);
  FitParamList::AppendSelectedUncorrelatedBANFFSystematics(fparm, 6+NSELECTED_NONBANFF);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithSelectedSystFullBANFFUnBANFF(void)
{
  FitParamList * fparm = new FitParamList(6+NSELECTED_NONBANFF+kNumOfParamsInBANFFfit); // NOTE: This will have to change for the 2012a 1Re and joint analysis

  FitParamList::Append3FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendSelectedNonBANFFSystematics    (fparm, 6);
  FitParamList::AppendAllUncorrelatedBANFFSystematics(fparm, 6+NSELECTED_NONBANFF);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFitWithSelectedSystNoBANFF(void)
{
  FitParamList * fparm = new FitParamList(6+NSELECTED_NONBANFF); // NOTE: This will have to change for the 2012a 1Re and joint analysis

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
  FitParamList * fparm = new FitParamList(2+kNumOfNuisanceParams1Rmu);

  FitParamList::Append2FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendAllCorrelatedBANFFSystematics (fparm, 2);
  FitParamList::AppendAllSuperkDetectorSystematics  (fparm, 2+kNumOfParamsInBANFFfit);
  FitParamList::AppendAllFSISystematics             (fparm, 2+kNumOfParamsInBANFFfit+kNumOfSKDetParams);
  FitParamList::AppendAllUncorrelatedXSecSystematics(fparm, 2+kNumOfParamsInBANFFfit+kNumOfSKDetParams+kNumOfFSIParams);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDisapFitWithSelectedSyst(void)
{
  FitParamList * fparm = new FitParamList(2+NSELECTED_NONBANFF+NSELECTED_BANFF);

  FitParamList::Append2FlavourOscillation(fparm, 0);

  //
  // add nuisance params
  //
  FitParamList::AppendSelectedNonBANFFSystematics       (fparm, 2);
  FitParamList::AppendSelectedCorrelatedBANFFSystematics(fparm, 2+NSELECTED_NONBANFF);

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
     fparm->AddParam( ip++, kPhypMixing23,  kInitFitSinsq2Theta23,  kMinFitSinsq2Theta23,  kMaxFitSinsq2Theta23,  1.0E-2,  false, false);
     fparm->AddParam( ip++, kPhypMixing13,  kInitFitSinsq2Theta13,  kMinFitSinsq2Theta13,  kMaxFitSinsq2Theta13,  1.0E-3,  false, false);
     fparm->AddParam( ip++, kPhypMixing12,  kInitFitSinsq2Theta12,  kMinFitSinsq2Theta12,  kMaxFitSinsq2Theta12,  1.0E-2,  false, false);
     break;
  case(kNuMixSinsqTheta):
     fparm->AddParam( ip++, kPhypMixing23,  kInitFitSinsqTheta23,   kMinFitSinsqTheta23,   kMaxFitSinsqTheta23,   1.0E-2,  false, false);
     fparm->AddParam( ip++, kPhypMixing13,  kInitFitSinsqTheta13,   kMinFitSinsqTheta13,   kMaxFitSinsqTheta13,   1.0E-3,  false, false);
     fparm->AddParam( ip++, kPhypMixing12,  kInitFitSinsqTheta12,   kMinFitSinsqTheta12,   kMaxFitSinsqTheta12,   1.0E-2,  false, false);
     break;
  default:
     delete fparm;
     return;
  }
  fparm->AddParam( ip++, kPhypDeltaCP,   0,                -kpi,          kpi,           1.0E-2,  false, false);
  fparm->AddParam( ip++, kPhypDm23sq,    kInitFitDm23sq,   kMinFitDm23sq, kMaxFitDm23sq, 1.0E-4,  false, false);
  fparm->AddParam( ip++, kPhypDm12sq,    kInitFitDm12sq,   kMinFitDm12sq, kMaxFitDm12sq, 2.0E-6,  false, false);

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
    fparm->AddParam( ip++, kPhypMixing23, kInitFitSinsq2Theta23, kMinFitSinsq2Theta23,  kMaxFitSinsq2Theta23,  0.01,  false, false);
    break;
  case(kNuMixSinsqTheta):
    fparm->AddParam( ip++, kPhypMixing23, kInitFitSinsqTheta23,  kMinFitSinsqTheta23,   kMaxFitSinsqTheta23,   0.01,  false, false);
    break;
  default:
    delete fparm;
    return;
  }
  fparm->AddParam( ip++, kPhypDm23sq, kInitFitDm23sq, kMinFitDm23sq, kMaxFitDm23sq, 1E-4, false, false);

}
//________________________________________________________________________________________
void FitParamList::AppendAllUncorrelatedXSecSystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.5, -0.,              1.,              SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  //fparm->AddParam( ip++, kNuisXSecENuShape,      0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true); // not used in numu analysis!
  fparm->AddParam( ip++, kNuisXSecBindE,         0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
#if BANFF_BINS <= 1
  fparm->AddParam( ip++, kNuisXSecNormCCQEE1,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormCCQEE2,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormCC1piE1,   0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
#endif
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormNC1piC,    0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  //fparm->AddParam( ip++, kNuisXSecNormNCcoh,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true); // not used in numu analysis!
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);

}
//________________________________________________________________________________________
void FitParamList::AppendAllSuperkDetectorSystematics(FitParamList * fparm, int start)
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

} //________________________________________________________________________________________
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
  TMatrixD * FSICovMtrx = util::XML2TMatrixD(xml,root_node,"fsi/covariance_matrix");
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
    PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + i);
    double err = sqrt(FSICov[i][i]);
    double lowerr = TMath::Max(0., 1.0-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p,  1.0, lowerr, 1.0+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }

}
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
  for(int i = 0; i < kNumOfParamsInBANFFfit; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisBANFF00 + i);
    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    double lowerr = TMath::Max(0., bf-(err*SYST_NSIGMA));
    fparm->AddParam(ip++, p, bf, lowerr, bf+(err*SYST_NSIGMA),  err*SYST_STEP,  false, true);
  }
#if BANFF_BINS <= 2
  fparm -> SetNuisance(kNuisBANFFNC1pi0 ,false); // not used in numu analysis!
#endif

}
//________________________________________________________________________________________
void FitParamList::AppendAllUncorrelatedBANFFSystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  // Uncorrelated BANFF 
  for(int i = 0; i < kNumOfParamsInBANFFfit; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisUnBANFF00 + i);
    fparm->AddParam( ip++, p, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  }

}
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

  // SK detector
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCnonQE,  0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuNC,       0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);
  //fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);

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

}
//________________________________________________________________________________________
void FitParamList::AppendSelectedUncorrelatedBANFFSystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//
  return;

  //int ip = start; // first unallocated element in the param list

  // Uncorrelated BANFF 
  //don't have any selected, but if we did, it'd be lines like:
  //fparm->AddParam( ip++, kNuisUnBANFF00, 0.0, -1.0*SYST_NSIGMA, 1.0*SYST_NSIGMA, SYST_STEP,  false, true);

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

   fParListName = "";
   fSystMode = -1;

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
   fMigradError  = new double[fNParams];
   fHesseError   = new double[fNParams];
   fPlusError    = new double[fNParams];
   fMinusError   = new double[fNParams];

   fShowCurrValuesOnly    = false;
   fShowBestFitValuesOnly = false;

   fMigradErrorsSet = false;
   fHesseErrorsSet  = false;
   fMinosErrorsSet  = false;
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
   delete [] fMigradError;
   delete [] fHesseError;
   delete [] fPlusError;
   delete [] fMinusError;

   fParListName = "";
   fSystMode = -1;
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

