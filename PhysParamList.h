//____________________________________________________________________________
/*!

\class    t2kana::yr1numu::PhysParamList

\brief    List of physics parameters included in the fit

\author   STFC, Rutherford Appleton Laboratory / Oxford T2K groups

          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          Jelena Ilic         <Jelena.Ilic \at stfc.ac.uk>
          Nick Grant          <nick.grant \at stfc.ac.uk>
          Geoff Pearce        <geoff.pearce \at stfc.ac.uk>
          Alfons Weber        <alfons.weber \at stfc.ac.uk>

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_phys_param_list_h_
#define _t2kana_phys_param_list_h_

#include <iostream>
#include <string>
#include <map>

#include <TVirtualFitter.h>

#include "year1_numu_osc_analysis/Conventions.h"

using std::string;
using std::map;

namespace t2kana {
namespace yr1numu {

class PhysParamList {

public:
   PhysParamList(int nparams);
  ~PhysParamList();

   // Get information for the physics parameter list
   int    Size      (void)                                                const { return fNParams;                  }
   string Name      (t2kana::yr1numu::conv::PhysParamId_t p) const { return fName     [this->Id(p)];   }
   double CurrValue (t2kana::yr1numu::conv::PhysParamId_t p) const { return fCurrValue[this->Id(p)];   }
   double InitValue (t2kana::yr1numu::conv::PhysParamId_t p) const { return fInitValue[this->Id(p)];   }
   double MinValue  (t2kana::yr1numu::conv::PhysParamId_t p) const { return fMinValue [this->Id(p)];   }
   double MaxValue  (t2kana::yr1numu::conv::PhysParamId_t p) const { return fMaxValue [this->Id(p)];   }
   double StepSize  (t2kana::yr1numu::conv::PhysParamId_t p) const { return fStepSize [this->Id(p)];   }
   bool   IsFixed   (t2kana::yr1numu::conv::PhysParamId_t p) const { return fIsFixed  [this->Id(p)];   }

   // Add a parameter
   bool AddParam (
     int ip, t2kana::yr1numu::conv::PhysParamId_t p, double init, 
     double min, double max, double step, bool fix=false, string name="");

   // Update current/running value for the given parameter
   void Update(t2kana::yr1numu::conv::PhysParamId_t p, double v) 
   { 
     fCurrValue [this->Id(p)] = v; 
   }

   // allow modifying the initial value (minimization may be sensitive to this choice)
   void SetInitValue(t2kana::yr1numu::conv::PhysParamId_t p, double v) 
   { 
     fInitValue [this->Id(p)] = v; 
   }

   // Fix a parameter
   void FixParameter (t2kana::yr1numu::conv::PhysParamId_t p) 
   { 
     fIsFixed [this->Id(p)] = true;  
   }

   // Set parameter list to input TVirtualFitter 
   void SetTo (TVirtualFitter * fitter);
   
   //
   void Sync (double * par);

   // Print
   void Print (ostream & stream) const;
   friend ostream & operator << (ostream & stream, const PhysParamList & plst);

private:

   void CleanUp (void);
   void Init    (int nparams);
   int  Id      (t2kana::yr1numu::conv::PhysParamId_t p) const; // PhysParamId_t --> array slot

   map<t2kana::yr1numu::conv::PhysParamId_t, int> fIdMap;

   int      fNParams;   
   double * fCurrValue; //[fNParams]
   double * fInitValue; //[fNParams]
   double * fMinValue;  //[fNParams]
   double * fMaxValue;  //[fNParams]
   double * fStepSize;  //[fNParams]
   string * fName;      //[fNParams]
   bool *   fIsFixed;   //[fNParams]
};

} // yr1numu
} // t2kana
  
#endif


