//________________________________________________________________________________________
/*

Copyright (c) 2010-2014, The VALOR group. All rights reserved.

ChangeLog:

*/
//________________________________________________________________________________________

#include <string>
#include <cstdlib>

#include "Opt.h"
#include "Rnd.h"
#include "ParameterStorage.h"
#include "Utilities.h"

using std::endl;
using std::string;
using std::cerr;
using std::cout;

using namespace t2kana;
using namespace t2kana::conv;

ClassImp(ParameterStorage);

//________________________________________________________________________________________
namespace t2kana {
  ostream & operator<<(ostream & stream, const ParameterStorage & obj)
  {
    obj.Print(stream);
    return stream;
  }
}
//________________________________________________________________________________________
ParameterStorage * ParameterStorage::fInstance = 0;
//________________________________________________________________________________________
ParameterStorage::ParameterStorage()
{
  this->SetDefaults();
  fInstance = 0;
}
//________________________________________________________________________________________
ParameterStorage::ParameterStorage(const ParameterStorage &) 
{

}
//________________________________________________________________________________________
ParameterStorage::~ParameterStorage()
{

}
//________________________________________________________________________________________
ParameterStorage * ParameterStorage::Instance()
{
  if(fInstance == 0) {
    static ParameterStorage::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new ParameterStorage;
  }
  return fInstance;
}
//________________________________________________________________________________________
double ParameterStorage::GetGlobalBestFit(PhysParamId_t p) const
{
  assert(p!=kPhypUndefined);
  return fGlobalBestFit[p];
}
//________________________________________________________________________________________
double ParameterStorage::GetGlobalFitError(PhysParamId_t p) const
{
  assert(p!=kPhypUndefined);
  return fGlobalFitError[p];
}
//________________________________________________________________________________________
double ParameterStorage::GetInitialValue(PhysParamId_t p) const
{
  assert(p!=kPhypUndefined);
  return fInitialValue[p];
}
//________________________________________________________________________________________
double ParameterStorage::GetRandomGaus(PhysParamId_t p) const
{
  assert(p!=kPhypUndefined);
  
  double mean  = fGlobalBestFit[p]; 
  double sigma = fGlobalFitError[p]*mean;
  
  assert(mean!=kBadResult);
  double randval = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
  cerr<<util::AsString(p)<<"\t"<<mean<<"\t"<<sigma<<"\t"<<randval<<endl;
  return randval;
}
//________________________________________________________________________________________
void ParameterStorage::SetDefaults(void) 
{
  // set empty values
  for (int ip = 0 ; ip < knppmax ; ip++ ) {
    fGlobalBestFit [ip] = kBadResult;
    fGlobalFitError[ip] = kBadResult;
    fInitialValue  [ip] = 0.;
    //~ fInitialError  [ip] = 1.;
  }
  
  assert(not Opt::Instance()->Undefined());
  
  NuOscDmsqScheme_t odmsq = Opt::Instance()->gOscDmsqScheme;
  NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
  
  double mean;
  
  //dm12
  fGlobalBestFit  [kPhypDm12sq] = kGlobFitDm12sq;
  fGlobalFitError [kPhypDm12sq] = kGlobFitErrDm12sq;
  fInitialValue   [kPhypDm12sq] = kInitFitDm12sq;
 
  //dm23
  switch (odmsq) {
    case(kNuDmsqFogliLisi):
      mean  =  kGlobFitDmsqFL;
      break;
    case(kNuDmsq12andAbs23):
      mean  = kGlobFitDm23sq; 
      break;
    case(kNuDmsq12andAbs31):
      if(Opt::Instance()->gMassHierarchyAssumed == kNmhNormal)
        mean  = kGlobFitDm23sq + kGlobFitDm12sq; 
      else if(Opt::Instance()->gMassHierarchyAssumed == kNmhInverted)
        mean  = kGlobFitDm23sq - kGlobFitDm12sq; 
      break;
    default:
      mean  = kGlobFitDm23sq;  
  }
  fGlobalBestFit  [kPhypDm23sq] = mean;
  fGlobalFitError [kPhypDm23sq] = kGlobFitErrDm23sq;
  fInitialValue   [kPhypDm23sq] = fGlobalBestFit[kPhypDm23sq];
  
  // mixing angles
  switch(omp) {
    case(kNuMixSinsq2Theta):
    {
      fGlobalBestFit  [kPhypMixing12] = kGlobFitSinsq2Theta12;
      fGlobalFitError [kPhypMixing12] = kGlobFitErrSinsq2Theta12;
      fInitialValue   [kPhypMixing12] = fGlobalBestFit[kPhypMixing12];
      
      fGlobalBestFit  [kPhypMixing23] = kGlobFitSinsq2Theta23;
      fGlobalFitError [kPhypMixing23] = kGlobFitErrSinsq2Theta23;
      fInitialValue   [kPhypMixing23] = fGlobalBestFit[kPhypMixing23];
      
      fGlobalBestFit  [kPhypMixing13] = kGlobFitSinsq2Theta13;
      fGlobalFitError [kPhypMixing13] = kGlobFitErrSinsq2Theta13;
      fInitialValue   [kPhypMixing13] = fGlobalBestFit[kPhypMixing13];
      break;
    }
    case(kNuMixSinsqTheta):
    {
      fGlobalBestFit  [kPhypMixing12] = kGlobFitSinsqTheta12;
      fGlobalFitError [kPhypMixing12] = kGlobFitErrSinsqTheta12;
      fInitialValue   [kPhypMixing12] = fGlobalBestFit[kPhypMixing12];
      
      fGlobalBestFit  [kPhypMixing23] = kGlobFitSinsqTheta23;
      fGlobalFitError [kPhypMixing23] = kGlobFitErrSinsqTheta23;
      fInitialValue   [kPhypMixing23] = fGlobalBestFit[kPhypMixing23];
      
      fGlobalBestFit  [kPhypMixing13] = kGlobFitSinsqTheta13;
      fGlobalFitError [kPhypMixing13] = kGlobFitErrSinsqTheta13;
      fInitialValue   [kPhypMixing13] = fGlobalBestFit[kPhypMixing13];
      break;
    }
    default: exit(-1);
  }
}
//________________________________________________________________________________________
void ParameterStorage::Print(ostream & stream) const
{
  stream << "ParameterStorage/Print >> Unimplemented"
         << endl;
}
//________________________________________________________________________________________
void ParameterStorage::SetGlobalBestFit  (conv::PhysParamId_t p, double value)
{
  assert(p>=0);
  assert(p<knppmax);
  fGlobalBestFit[p] = value;
}
//________________________________________________________________________________________
void ParameterStorage::SetGlobalFitError (conv::PhysParamId_t p, double value)
{
  assert(p>=0);
  assert(p<knppmax);
  fGlobalFitError[p] = value;
}
//________________________________________________________________________________________
void ParameterStorage::SetInitialValue   (conv::PhysParamId_t p, double value)
{
  assert(p>=0);
  assert(p<knppmax);
  fInitialValue[p] = value;
}
//________________________________________________________________________________________
