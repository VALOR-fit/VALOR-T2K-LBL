//____________________________________________________________________________
/*!

\class    t2kana::FitParamList

\brief    List of parameters included in the fit

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_fit_param_list_h_
#define _t2kana_fit_param_list_h_

#include <iostream>
#include <string>
#include <map>

#include <TVirtualFitter.h>

#include "Conventions.h"

using std::string;
using std::map;

namespace t2kana {

class FitParamList {

public:
   FitParamList(int nparams);
  ~FitParamList();

   // Get information for the fit parameters
   int    Size       (void)                          const { return fNParams;                   }
   string Name       (t2kana::conv::PhysParamId_t p) const { return fName      [this->Id(p)];   }
   double CurrValue  (t2kana::conv::PhysParamId_t p) const { return fCurrValue [this->Id(p)];   }
   double InitValue  (t2kana::conv::PhysParamId_t p) const { return fInitValue [this->Id(p)];   }
   double MinValue   (t2kana::conv::PhysParamId_t p) const { return fMinValue  [this->Id(p)];   }
   double MaxValue   (t2kana::conv::PhysParamId_t p) const { return fMaxValue  [this->Id(p)];   }
   double StepSize   (t2kana::conv::PhysParamId_t p) const { return fStepSize  [this->Id(p)];   }
   bool   IsFixed    (t2kana::conv::PhysParamId_t p) const { return fIsFixed   [this->Id(p)];   }
   bool   IsNuisance (t2kana::conv::PhysParamId_t p) const { return fIsNuisance[this->Id(p)];   }

   // Add a parameter
   bool AddParam (
     int ip, t2kana::conv::PhysParamId_t p, double init, 
     double min, double max, double step, bool fix=false, bool nuisance=false, string name="");

   // Check whether the input parameter was included in the physics param list
   bool Exists (t2kana::conv::PhysParamId_t p) const;

   // Update current/running value for the given parameter
   void Update(t2kana::conv::PhysParamId_t p, double v) 
   { 
     fCurrValue [this->Id(p)] = v; 
   }

   // allow modifying the initial value (minimization may be sensitive to this choice)
   void SetInitValue(t2kana::conv::PhysParamId_t p, double v) 
   { 
     fInitValue [this->Id(p)] = v; 
   }

   // Fix a parameter
   void FixParameter (t2kana::conv::PhysParamId_t p) 
   { 
     fIsFixed [this->Id(p)] = true;  
   }

   // Set parameter list to input TVirtualFitter 
   void SetTo (TVirtualFitter * fitter);
   
   //
   void Sync (double * par);

   // Chisq penalty term
   double ChisquarePenaltyTerm(void);

   // Print
   void Print (ostream & stream) const;
   friend ostream & operator << (ostream & stream, const FitParamList & plst);

   // Using the "Named Constructor" C++ idiom for fast creation of typical parameter lists.
   static FitParamList * NamedParameterList (string name);
   static FitParamList * NumuDisapFitNoSyst (void);
   static FitParamList * NumuDisapFitAllSyst(void);
   static FitParamList * NumuDecayFitNoSyst (void);
   static FitParamList * NumuDecohFitNoSyst (void);

private:

   void CleanUp (void);
   void Init    (int nparams);
   int  Id      (t2kana::conv::PhysParamId_t p) const; // PhysParamId_t --> array slot

   map<t2kana::conv::PhysParamId_t, int> fIdMap;

   int      fNParams;   
   double * fCurrValue;  //[fNParams]
   double * fInitValue;  //[fNParams]
   double * fMinValue;   //[fNParams]
   double * fMaxValue;   //[fNParams]
   double * fStepSize;   //[fNParams]
   string * fName;       //[fNParams]
   bool *   fIsFixed;    //[fNParams]
   bool *   fIsNuisance; //[fNParams]
};

} // t2kana
  
#endif


