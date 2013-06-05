//________________________________________________________________________________________
/*

 \copyright (c) 2010-2013, VaLOR group.

*/
//________________________________________________________________________________________

#include <iomanip>
#include <cstdlib>

#include "PhysParamList.h"
#include "Utilities.h"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;

using namespace t2kana;
using namespace t2kana::conv;

ClassImp(PhysParamList)

//________________________________________________________________________________________
namespace t2kana {
  ostream & operator << (ostream & stream, const PhysParamList & pl)
  {
    pl.Print(stream);
    return stream;
  }
}
//________________________________________________________________________________________
PhysParamList::PhysParamList() :
TObject()
{
  this->fParams.clear();
}
//________________________________________________________________________________________
PhysParamList::PhysParamList(const PhysParamList & pl) :
TObject()
{
  this->fParams.clear();

  map<int,double>::const_iterator it = pl.fParams.begin();
  for( ; it != pl.fParams.end(); ++it) {
     PhysParamId_t p = (PhysParamId_t) it->first;
     double        v = it->second;
     this->SetValue(p,v);
  }
}
//________________________________________________________________________________________
PhysParamList::~PhysParamList()
{
  fParams.clear();
}
//________________________________________________________________________________________
int PhysParamList::NParams(void) const
{
  return (int) fParams.size();
}
//________________________________________________________________________________________
double PhysParamList::Value(PhysParamId_t p) const
{
  map<int,double>::const_iterator it = fParams.find(p);
  if(it == fParams.end()) {
     cerr << "PhysParamList >> Parameter [" 
          << util::AsString(p) << "] was not set!" << endl;
     return -99999;
  } 
  return it->second;
}
//________________________________________________________________________________________
void PhysParamList::SetValue(PhysParamId_t p, double v)
{
  fParams[p] = v;
}
//________________________________________________________________________________________
void PhysParamList::Print(ostream & stream) const
{
  int i=0;
  stream << endl;
  map<int,double>::const_iterator it = fParams.begin();
  for( ; it != fParams.end(); ++it) {
     PhysParamId_t p = (PhysParamId_t) it->first;
     double        v = it->second;
     stream << setw(5) << i++
            << setw(75) << util::AsString(p) << " = "
            << setprecision(8) << setw(12) << v
            << endl;      
  }
}
//________________________________________________________________________________________
double PhysParamList::operator () (conv::PhysParamId_t p) const
{
  return this->Value(p);
}
//________________________________________________________________________________________

