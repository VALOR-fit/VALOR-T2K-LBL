//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/
//________________________________________________________________________________________

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
void FitParamList::SetErrors(double * errplus, double * errminus, double * errpara)
{
  for(int ip = 0; ip < fNParams; ip++) {
     fPlusError [ip] = errplus [ip];
     fMinusError[ip] = errminus[ip];
     fParaError [ip] = errpara [ip];
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
  //
  // 3-flavour oscillation fit with/without systematics included 
  //
  if(name == "nuosc3f") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> SetParListName(name);
    return fp;
  }
  if(name == "nuosc3f+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSyst();
    fp -> SetParListName(name);
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
    return fp;
  }
  if(name == "nuosc3f_12fix") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm12sq  );
    fp -> SetParListName(name);
    return fp;
  }
  if(name == "nuosc3f_numudisapp_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    return fp;
  }
  if(name == "nuosc3f_numudisapp_fix4+nuisance") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSyst();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    return fp;
  }
  if(name == "nuosc3f_numudisapp_fix4+nuisanceselection") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFitWithAllSyst();
    fp -> FixParameter (kPhypMixing13);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDeltaCP);
    fp -> FixParameter (kPhypDm12sq);
    for(int i = 0; i < kNumOfFSIParams; i++) {
      PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + i);
      fp -> FixParameter (p);
      fp -> SetNuisance  (p, false);
    }
    fp -> FixParameter (kNuisSuperkEff1RmuCCQEoth);
    fp -> SetNuisance  (kNuisSuperkEff1RmuCCQEoth ,false);
    fp -> FixParameter (kNuisSuperkEff1RmuCCQErc0u);
    fp -> SetNuisance  (kNuisSuperkEff1RmuCCQErc0u ,false);
    fp -> FixParameter (kNuisSuperkEff1RmuCCQErc1u);
    fp -> SetNuisance  (kNuisSuperkEff1RmuCCQErc1u ,false);
    fp -> FixParameter (kNuisSuperkEff1RmuCCQErc2u);
    fp -> SetNuisance  (kNuisSuperkEff1RmuCCQErc2u ,false);
    fp -> FixParameter (kNuisSuperkEff1RmuNueCC);
    fp -> SetNuisance  (kNuisSuperkEff1RmuNueCC ,false);
    fp -> FixParameter (kNuisSuperkEnergyScale);
    fp -> SetNuisance  (kNuisSuperkEnergyScale ,false);
    fp -> FixParameter (kNuisXSecFermiMom);
    fp -> SetNuisance  (kNuisXSecFermiMom ,false);
    fp -> FixParameter (kNuisXSecSpecFunc);
    fp -> SetNuisance  (kNuisXSecSpecFunc ,false);
    fp -> FixParameter (kNuisXSecWshape);
    fp -> SetNuisance  (kNuisXSecWshape ,false);
    fp -> FixParameter (kNuisXSecShapeCCoth);
    fp -> SetNuisance  (kNuisXSecShapeCCoth ,false);
#if BANFF_BINS != 2
    fp -> FixParameter (kNuisXSecNormCCQEE1);
    fp -> SetNuisance  (kNuisXSecNormCCQEE1 ,false);
    fp -> FixParameter (kNuisXSecNormCCQEE2);
    fp -> SetNuisance  (kNuisXSecNormCCQEE2 ,false);
    fp -> FixParameter (kNuisXSecNormCC1piE1);
    fp -> SetNuisance  (kNuisXSecNormCC1piE1 ,false);
#endif
    fp -> FixParameter (kNuisXSecNormCCcoh);
    fp -> SetNuisance  (kNuisXSecNormCCcoh ,false);
    fp -> FixParameter (kNuisXSecNormNuEToNuMu);
    fp -> SetNuisance  (kNuisXSecNormNuEToNuMu ,false);
    fp -> FixParameter (kNuisXSecNormAntiNuToNu);
    fp -> SetNuisance  (kNuisXSecNormAntiNuToNu ,false);
    fp -> SetParListName(name);
    return fp;
  }
  if(name == "nuosc3f_nueapp_fix4") {
    FitParamList * fp = FitParamList::NuOsc3FlavourFit();
    fp -> FixParameter (kPhypMixing23);
    fp -> FixParameter (kPhypMixing12);
    fp -> FixParameter (kPhypDm23sq);
    fp -> FixParameter (kPhypDm12sq);
    fp -> SetParListName(name);
    return fp;
  }

  //
  // 2-flavour numu-disappearance oscillation fit with/without systematics included 
  //
  if(name == "numudisapp2f") {
    FitParamList * fp = FitParamList::NumuDisapFit();
    fp -> SetParListName(name);
    return fp;
  }
  if(name == "numudisapp2f+nuisance") {
    FitParamList * fp = FitParamList::NumuDisapFitWithAllSyst();
    fp -> SetParListName(name);
    return fp;
  }

  //
  // numu decay fit
  //
  if(name == "decay") {
    FitParamList * fp = FitParamList::NumuDecayFit();
    fp -> SetParListName(name);
    return fp;
  }

  //
  // numu decoherence fit
  //
  if(name == "decoh") {
    FitParamList * fp = FitParamList::NumuDecohFit();
    fp -> SetParListName(name);
    return fp;
  }

  return 0;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NuOsc3FlavourFit(void)
{
  FitParamList * fparm = new FitParamList(6);

  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  switch(omp) {
  case(kNuMixSinsq2Theta):
    fparm->AddParam(0, kPhypMixing23,  kInitFitSinsq2Theta23,  kMinFitSinsq2Theta23,  10,  1.0E-2,  false, false);
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
FitParamList * FitParamList::NuOsc3FlavourFitWithAllSyst(void)
{
  FitParamList * fparm = new FitParamList(6+kNumOfNuisanceParams1Rmu); // NOTE: This will have to change for the 2012a 1Re and joint analysis

  //
  // add osc params
  //
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
  FitParamList::AppendSystematics(fparm, 6);

  return fparm;
}
//________________________________________________________________________________________
FitParamList * FitParamList::NumuDisapFit(void)
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
FitParamList * FitParamList::NumuDisapFitWithAllSyst(void)
{
  FitParamList * fparm = new FitParamList(2+kNumOfNuisanceParams1Rmu);

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
  FitParamList::AppendSystematics(fparm, 2);

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
void FitParamList::AppendSystematics(FitParamList * fparm, int start)
{
// Add nuisance params
//

  int ip = start; // first unallocated element in the param list

  int nsigma = 3; // number of sigma away from nominal nuisance parameters are allowed to be

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

  // SK detector
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCQEoth,  0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCQErc0u, 0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCQErc1u, 0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCQErc2u, 0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuCCnonQE,  0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuNC,       0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEff1RmuNueCC,    0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisSuperkEnergyScale,     0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);

  // BANFF 
  for(int i = 0; i < kNumOfParamsInBANFFfit; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisBANFF00 + i);
    double bf = BANFFBestFit[i];
    double err = sqrt(BANFFCov[i][i]);
    fparm->AddParam(ip++, p, bf, bf-(err*nsigma), bf+(err*nsigma),  err/2,  false, true);
  }
  fparm -> SetNuisance(kNuisBANFFNC1pi0 ,false); // not used in numu analysis!

  // XSec
  fparm->AddParam( ip++, kNuisXSecFermiMom,      0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisXSecSpecFunc,      0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisXSecWshape,        0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisXSecShapeCCoth,    0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisXSecPiLessDelta,   0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  //fparm->AddParam( ip++, kNuisXSecENuShape,      0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true); // not used in numu analysis!
#if BANFF_BINS != 2
  fparm->AddParam( ip++, kNuisXSecNormCCQEE1,    0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormCCQEE2,    0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormCC1piE1,   0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
#endif
  fparm->AddParam( ip++, kNuisXSecNormCCcoh,     0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormNC1piC,    0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  //fparm->AddParam( ip++, kNuisXSecNormNCcoh,     0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true); // not used in numu analysis!
  fparm->AddParam( ip++, kNuisXSecNormNCoth,     0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormNuEToNuMu, 0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);
  fparm->AddParam( ip++, kNuisXSecNormAntiNuToNu,0.0, -1.0*nsigma, 1.0*nsigma, 5E-1,  false, true);

  // FSI
  for(int i = 0; i < kNumOfFSIParams; i++) {
    PhysParamId_t p = PhysParamId_t ((int) kNuisFSI000 + i);
    double err = sqrt(FSICov[i][i]);
    fparm->AddParam(ip++, p,  1.0, 1.0-(err*nsigma), 1.0+(err*nsigma),  err/2,  false, true);
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
void FitParamList::Init(int nparams)
{
   assert(nparams>0);

   fParListName = "";

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
   fPlusError    = new double[fNParams];
   fMinusError   = new double[fNParams];
   fParaError    = new double[fNParams];

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
   delete [] fPlusError;
   delete [] fMinusError;
   delete [] fParaError;

   fParListName = "";
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

