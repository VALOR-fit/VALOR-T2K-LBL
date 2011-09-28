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
   int    Size         (void)                          const { return fNParams;                      }
   string Name         (t2kana::conv::PhysParamId_t p) const { return fName         [this->Id(p)];   }
   double BestFitValue (t2kana::conv::PhysParamId_t p) const { return fBestFitValue [this->Id(p)];   }
   double CurrValue    (t2kana::conv::PhysParamId_t p) const { return fCurrValue    [this->Id(p)];   }
   double InitValue    (t2kana::conv::PhysParamId_t p) const { return fInitValue    [this->Id(p)];   }
   double MinValue     (t2kana::conv::PhysParamId_t p) const { return fMinValue     [this->Id(p)];   }
   double MaxValue     (t2kana::conv::PhysParamId_t p) const { return fMaxValue     [this->Id(p)];   }
   double StepSize     (t2kana::conv::PhysParamId_t p) const { return fStepSize     [this->Id(p)];   }
   bool   IsFixed      (t2kana::conv::PhysParamId_t p) const { return fIsFixed      [this->Id(p)];   }
   bool   IsNuisance   (t2kana::conv::PhysParamId_t p) const { return fIsNuisance   [this->Id(p)];   }

   // Add a parameter
   bool AddParam (
     int ip, t2kana::conv::PhysParamId_t p, double init, 
     double min, double max, double step, bool fix=false, bool nuisance=false, string name="");

   // Check whether the input parameter was included in the physics param list
   bool Exists (t2kana::conv::PhysParamId_t p) const;

   // Update current/running value for the given parameter
   void Update(t2kana::conv::PhysParamId_t p, double v);

   // Allow modifying the initial value (minimization may be sensitive to this choice)
   void SetInitValue(t2kana::conv::PhysParamId_t p, double v);

   // Reset all current param values to their specified initial values
   void ResetCurrToInitValues (void);

   // Copy current param values to best-fit values & vice-versa
   void CopyCurrToBestFitValues (void);
   void CopyBestFitToCurrValues (void);

   // Fix a parameter
   void FixParameter (t2kana::conv::PhysParamId_t p); 

   // Free a parameter
   void FreeParameter (t2kana::conv::PhysParamId_t p); 

   // Set as nuisance parameter (contributes to \chi^{2}_{penalty})
   void SetNuisance (t2kana::conv::PhysParamId_t p, bool tf);

   // Set parameter list to input TVirtualFitter 
   void SetTo (TVirtualFitter * fitter);
   
   // Sync fit param values stored this object with the values stored in the input array.
   void Sync (double * par);

   // Chisq penalty term
   double ChisquarePenaltyTerm(void);

   // Get number of free (not fixed) parameters
   int GetNumberFreeParameters(void);

   // Print
   void Print (ostream & stream) const;
   void ShowOnlyCurrValuesInNextPrint(void) const { fShowCurrValuesOnly = true; }
   friend ostream & operator << (ostream & stream, const FitParamList & plst);

   // Using the "Named Constructor" C++ idiom for fast creation of typical parameter lists.
   // The string argument of NamedParameterList() can take the following values:
   // - "nuosc3f"                        : Fit for sin^2(2theta12), sin^2(2theta13), sin^2(2theta23), delta_CP, |Dm23^2| and Dm12^2 (3-flavour formalism).
   // - "nuosc3f_12constrain"            : As "nuosc3f", but penalizing sin^2(2theta12), Dm12^2 devations from values measured in solar/reactor expts.
   // - "nuosc3f_numudisapp_fix4"        : As "nuosc3f" but fixing  sin^2(2theta12), sin^2(2theta13), delta_CP and Dm12^2.
   // - "nuosc3f_nueapp_fix4"            : As "nuosc3f" but fixing  sin^2(2theta12), sin^2(2theta23), |Dm23^2| and Dm12^2.
   // - "numudisapp2f"                   : Fit for sin^2(2theta23) and |Dm23^2| (2-flavour formalism)
   // - "numudisapp2f+nuisance"          : As "numudisapp2f" + the 34 nuisance parameters of 2010 analysis
   // - "decay"                          : Fit for sin^2(2theta) and \alpha_{decay}
   // - "decoh"                          : Fit for sin^2(2theta) and \mu_{decoherence}
   static FitParamList * NamedParameterList     (string name);
   static FitParamList * NuOsc3FlavourFitNoSyst (void);
   static FitParamList * NumuDisapFitNoSyst     (void);
   static FitParamList * NumuDisapFitAllSyst    (void);
   static FitParamList * NumuDecayFitNoSyst     (void);
   static FitParamList * NumuDecohFitNoSyst     (void);

   // Work-out the physics hypothesis
   conv::PhysHypothesis_t PhysicsHypothesis(void) const;

private:

   void CleanUp (void);
   void Init    (int nparams);
   int  Id      (t2kana::conv::PhysParamId_t p) const; // PhysParamId_t --> array slot

   t2kana::conv::PhysParamId_t WhichParam(int id) const; // array slot --> PhysParamId_t

   map<t2kana::conv::PhysParamId_t, int> fIdMap;

   int      fNParams;   
   double * fBestFitValue; //[fNParams]
   double * fCurrValue;    //[fNParams]
   double * fInitValue;    //[fNParams]
   double * fMinValue;     //[fNParams]
   double * fMaxValue;     //[fNParams]
   double * fStepSize;     //[fNParams]
   string * fName;         //[fNParams]
   bool *   fIsFixed;      //[fNParams]
   bool *   fIsNuisance;   //[fNParams]

   mutable bool fShowCurrValuesOnly;
};

} // t2kana
  
#endif


