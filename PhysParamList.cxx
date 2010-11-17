//________________________________________________________________________________________
/*

*/
//________________________________________________________________________________________

//#define _debug_

#include <cassert>
#include <cstdlib>
#include <iomanip>

#include "year1_numu_osc_analysis/PhysParamList.h"
#include "year1_numu_osc_analysis/Utilities.h"

using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;

using namespace t2kana;
using namespace t2kana::yr1numu;
using namespace t2kana::yr1numu::conv;

//________________________________________________________________________________________
namespace t2kana {
namespace yr1numu {
  ostream & operator << (ostream & stream, const PhysParamList & plst)
  {
    plst.Print(stream);
    return stream;
  }
}
}
//________________________________________________________________________________________
PhysParamList::PhysParamList(int nparams)
{
  this->Init(nparams);
}
//________________________________________________________________________________________
PhysParamList::~PhysParamList()
{
  this->CleanUp();
}
//________________________________________________________________________________________
bool PhysParamList::AddParam (int ip, PhysParamId_t p, 
        double init, double min, double max, double step, bool fix, string name)
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
  fName      [ip] = name;     

  return true;
}
//________________________________________________________________________________________
bool PhysParamList::Exists(PhysParamId_t p) const
{
  map<PhysParamId_t, int>::const_iterator iter = fIdMap.find(p);
  if(iter != fIdMap.end()) return true;
  return false;
}
//________________________________________________________________________________________
void PhysParamList::SetTo(TVirtualFitter * fitter)
{
  if(!fitter) return;

  for(int iparam = 0; iparam < this->Size(); iparam++) {
    PhysParamId_t p = (PhysParamId_t) iparam;
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
void PhysParamList::Sync(double * par)
{
  if(!par) return;
  for(int ip = 0; ip < fNParams; ip++) {
     fCurrValue[ip] = par[ip];
  }
}
//________________________________________________________________________________________
void PhysParamList::Print(ostream & stream) const
{
  stream << endl
         << setw(6)  << "id"
         << setw(20) << "name"
         << setw(10) << "curr"
         << setw(10) << "init"
         << setw(10) << "min"
         << setw(10) << "max"
         << setw(8)  << "step"
         << setw(8)  << "fixed?"
         << endl;
  for(int ip=0; ip<fNParams;  ip++) {
    stream << setw(6)  << ip 
           << setw(20) << fName[ip]
           << setprecision(6)
           << setw(10) << fCurrValue[ip]
           << setw(10) << fInitValue[ip]
           << setw(10) << fMinValue[ip]
           << setw(10) << fMaxValue[ip]
           << setw(8) << fStepSize[ip]
           << setw(8) << ((fIsFixed[ip]) ? " **yes** " : "no")
           << endl;
  }
}
//________________________________________________________________________________________
PhysParamList * PhysParamList::NumuDisapFitNoSyst(void)
{
  PhysParamList * fparm = new PhysParamList(2);

  //                                     init   min   max  step  fix?   name
  fparm->AddParam(0, kPhypSinsq2Theta23, 0.9,   0.7,  1.0, 0.02, false, "sinsq2theta23");
  fparm->AddParam(1, kPhypDm23sq,        2E-3, 1E-4, 8E-3, 1E-5, false, "dmsq23");

  return fparm;
}
//________________________________________________________________________________________
PhysParamList * PhysParamList::NumuDisapFitShapeOnlyNoSyst(void)
{
  PhysParamList * fparm = new PhysParamList(3);

  //                                     init   min   max  step  fix?   name
  fparm->AddParam(0, kPhypSinsq2Theta23, 0.9,   0.7,  1.0, 0.02, false, "sinsq2theta23");
  fparm->AddParam(1, kPhypDm23sq,        2E-3, 1E-4, 8E-3, 1E-5, false, "dmsq23");
  fparm->AddParam(2, kPhypNorm,          1, 1E-2, 1E+2, 0.1 , false, "norm");

  return fparm;
}
//________________________________________________________________________________________
void PhysParamList::Init(int nparams)
{
   assert(nparams>0);

   fNParams = nparams; 

   fCurrValue = new double[fNParams];
   fInitValue = new double[fNParams];
   fMinValue  = new double[fNParams];
   fMaxValue  = new double[fNParams];
   fStepSize  = new double[fNParams];
   fName      = new string[fNParams];
   fIsFixed   = new bool  [fNParams];
}
//________________________________________________________________________________________
void PhysParamList::CleanUp(void)
{
   fIdMap.clear();
   delete [] fCurrValue;
   delete [] fInitValue;
   delete [] fMinValue;
   delete [] fMaxValue;
   delete [] fStepSize;
   delete [] fName;
   delete [] fIsFixed;
}
//________________________________________________________________________________________
int PhysParamList::Id(PhysParamId_t p) const
{
   map<PhysParamId_t, int>::const_iterator it = fIdMap.find(p);
   if(it == fIdMap.end()) {
     cerr << "Parameter : " << util::AsString(p) 
          << " was not included in current PhysParamList" 
         << endl;
     exit(1);
   }
   int id = it->second;
   return id;
}
//________________________________________________________________________________________
