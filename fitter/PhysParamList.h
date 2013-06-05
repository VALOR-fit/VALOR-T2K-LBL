//____________________________________________________________________________
/*!
   
\class    t2kana::PhysParamList
      
\brief    A set of physics parameter values.
   
\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  June 05, 2013

\copyright (c) 2010-2013, VaLOR group.

*/
//____________________________________________________________________________
   
#ifndef _t2kana_phys_param_list_h_
#define _t2kana_phys_param_list_h_

#include <iostream>
#include <string>
#include <map>

#include <TObject.h>

#include "Conventions.h"
#include "FitParamList.h"

using std::ostream;
using std::string;
using std::map;

namespace t2kana {

class PhysParamList: public TObject
{
public:
   PhysParamList();
   PhysParamList(const PhysParamList & pl);
  ~PhysParamList();

   int    NParams  (void) const;
   double Value    (conv::PhysParamId_t p) const;
   void   SetValue (conv::PhysParamId_t p, double v);
   void   Clear    (void);
   void   Append   (const PhysParamList & pl);
   void   Print    (ostream & stream) const;

   double           operator () (conv::PhysParamId_t p) const;
   friend ostream & operator << (ostream & stream, const PhysParamList & pl);

   const map<int,double> & Map(void) const { return fParams; }

private:

   map<int,double> fParams;

ClassDef(PhysParamList,1)
};

} // t2kana

#endif
