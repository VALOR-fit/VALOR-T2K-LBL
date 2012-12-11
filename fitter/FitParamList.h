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
#include <vector>

#include <TVirtualFitter.h>
#include <Minuit2/MnUserParameters.h>

#include "Conventions.h"

using std::string;
using std::map;

namespace t2kana {

class FitParamList {

public:
   FitParamList(int nparams);
  ~FitParamList();

   // Get information for the fit parameters
   int    Size           (void)                          const { return fNParams;                      }
   string Name           (t2kana::conv::PhysParamId_t p) const { return fName         [this->Id(p)];   }
   double BestFitValue   (t2kana::conv::PhysParamId_t p) const { return fBestFitValue [this->Id(p)];   }
   double CurrValue      (t2kana::conv::PhysParamId_t p) const { return fCurrValue    [this->Id(p)];   }
   double InitValue      (t2kana::conv::PhysParamId_t p) const { return fInitValue    [this->Id(p)];   }
   double MinValue       (t2kana::conv::PhysParamId_t p) const { return fMinValue     [this->Id(p)];   }
   double MaxValue       (t2kana::conv::PhysParamId_t p) const { return fMaxValue     [this->Id(p)];   }
   double StepSize       (t2kana::conv::PhysParamId_t p) const { return fStepSize     [this->Id(p)];   }
   bool   IsFixed        (t2kana::conv::PhysParamId_t p) const { return fIsFixed      [this->Id(p)];   }
   bool   IsNuisance     (t2kana::conv::PhysParamId_t p) const { return fIsNuisance   [this->Id(p)];   }
   double MigradError    (t2kana::conv::PhysParamId_t p) const { return fMigradError  [this->Id(p)];   }
   double HesseError     (t2kana::conv::PhysParamId_t p) const { return fHesseError   [this->Id(p)];   }
   double PlusError      (t2kana::conv::PhysParamId_t p) const { return fPlusError    [this->Id(p)];   }
   double MinusError     (t2kana::conv::PhysParamId_t p) const { return fMinusError   [this->Id(p)];   }

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

   // Set the errors
   // used with FitterKernel
   void SetErrors(double * errplus, double * errminus, double * errmigrad);
   // used with FitterKernelMinuit2
   void SetMigradErrors(std::vector<double> err);
   void SetHesseErrors (std::vector<double> err);
   void SetMinosErrors (std::vector<double> plus, std::vector<double> minus);

   // Fix a parameter
   void FixParameter (t2kana::conv::PhysParamId_t p); 

   // Free a parameter
   void FreeParameter (t2kana::conv::PhysParamId_t p); 

   // Set as nuisance parameter (contributes to \chi^{2}_{penalty})
   void SetNuisance (t2kana::conv::PhysParamId_t p, bool tf);

   // Set parameter list to input TVirtualFitter 
   void SetTo (TVirtualFitter * fitter);
   void SetTo (ROOT::Minuit2::MnUserParameters & params);
   
   // Sync fit param values stored this object with the values stored in the input array.
   void Sync (double * par);
   void Sync (std::vector<double> par);

   // Get number of free (not fixed) parameters
   int GetNumberFreeParameters(void);

   // Print
   void Print (ostream & stream) const;
   friend ostream & operator << (ostream & stream, const FitParamList & plst);

   // Methods with an one-off effect in Print()
   void ShowOnlyCurrValuesInNextPrint   (void) const { fShowCurrValuesOnly    = true; }
   void ShowOnlyBestFitValuesInNextPrint(void) const { fShowBestFitValuesOnly = true; }


   // Using the "Named Constructor" C++ idiom for fast creation of typical parameter lists.
   // The string argument of NamedParameterList() can take the following values:
   // - "nuosc3f"                        : Fit for sin^2(2theta12), sin^2(2theta13), sin^2(2theta23), delta_CP, |Dm23^2| and Dm12^2 (3-flavour formalism).
   // - "nuosc3f+nuisance"               : As "nuosc3f" + all nuisance parameters
   // - "nuosc3f+nuisanceselection"      : As "nuosc3f" + a selection of nuisance parameters
   // - "nuosc3f_12constrain"            : As "nuosc3f", but penalizing sin^2(2theta12), Dm12^2 devations from values measured in solar/reactor expts.
   // - "nuosc3f_12fix"                  : As "nuosc3f", but fixing sin^2(2theta12), Dm12^2 devations to values measured in solar/reactor expts.
   // - "nuosc3f_numudisapp_fix4"        : As "nuosc3f" but fixing  sin^2(2theta12), sin^2(2theta13), delta_CP and Dm12^2.
   // - "nuosc3f_numudisapp_fix4+nuisanceselection": As "nuosc3f" but fixing  sin^2(2theta12), sin^2(2theta13), delta_CP and Dm12^2 + a selection of nuisance parameters.
   // - "nuosc3f_numudisapp_fix4+nuisanceselectionnobanff": As "nuosc3f" but fixing  sin^2(2theta12), sin^2(2theta13), delta_CP and Dm12^2 + a selection of nuisance parameters (no BANFF"s included).
   // - "nuosc3f_numudisapp_fix4+nuisancefullbanffselection": As "nuosc3f" but fixing  sin^2(2theta12), sin^2(2theta13), delta_CP and Dm12^2 + a selection of nuisance parameters (including the full BANFF matrix).
   // - "nuosc3f_numudisapp_fix4+nuisance": As "nuosc3f" but fixing  sin^2(2theta12), sin^2(2theta13), delta_CP and Dm12^2 + all nuisance parameters.
   // - "nuosc3f_nueapp_fix4"            : As "nuosc3f" but fixing  sin^2(2theta12), sin^2(2theta23), |Dm23^2| and Dm12^2.
   // - "numudisapp2f"                   : Fit for sin^2(2theta23) and |Dm23^2| (2-flavour formalism)
   // - "numudisapp2f+nuisanceselection" : As "numudisapp2f" + selected nuisance parameters 
   // - "numudisapp2f+nuisance"          : As "numudisapp2f" + all nuisance parameters 
   // - "decay"                          : Fit for sin^2(2theta) and \alpha_{decay}
   // - "decoh"                          : Fit for sin^2(2theta) and \mu_{decoherence}
   static FitParamList * NamedParameterList (string name);

   // Get the name used in NamedParameterList
   string GetParListName (void) const { return fParListName; }
   int    GetSystMode     (void) const { return fSystMode;     }

   // Work-out the physics hypothesis
   conv::PhysHypothesis_t PhysicsHypothesis(void) const;

   //Are the errors set?
   bool MigradErrorsSet (void) const { return fMigradErrorsSet; }
   bool HesseErrorsSet  (void) const { return fHesseErrorsSet;  }
   bool MinosErrorsSet  (void) const { return fMinosErrorsSet;  }

private:

   void CleanUp (void);
   void Init    (int nparams);
   int  Id      (t2kana::conv::PhysParamId_t p) const; // PhysParamId_t --> array slot
   void SetParListName (string s)   { fParListName = s; }
   void SetSystMode     (int m)      { fSystMode     = m; }

   // Various functions called by NamedParameterList(name)
   static FitParamList * NuOsc3FlavourFit                 (void);
   static FitParamList * NuOsc3FlavourFitWithAllSyst      (void);
   static FitParamList * NuOsc3FlavourFitWithSelectedSyst (void);
   static FitParamList * NuOsc3FlavourFitWithSelectedSystFullBANFF (void);
   static FitParamList * NuOsc3FlavourFitWithSelectedSystNoBANFF (void);
   static FitParamList * NumuDisapFit                     (void);
   static FitParamList * NumuDisapFitWithSelectedSyst     (void);
   static FitParamList * NumuDisapFitWithAllSyst          (void);
   static FitParamList * NumuDecayFit                     (void);
   static FitParamList * NumuDecohFit                     (void);
   static void           Append3FlavourOscillation                  (FitParamList * fp, int start);
   static void           Append2FlavourOscillation                  (FitParamList * fp, int start);
   static void           AppendAllUncorrelatedXSecSystematics       (FitParamList * fp, int start);
   static void           AppendAllSuperkDetectorSystematics         (FitParamList * fp, int start);
   static void           AppendAllFSISystematics                    (FitParamList * fp, int start);
   static void           AppendAllCorrelatedBANFFSystematics        (FitParamList * fp, int start);
   static void           AppendAllUncorrelatedBANFFSystematics      (FitParamList * fp, int start);
   static void           AppendSelectedNonBANFFSystematics          (FitParamList * fp, int start);
   static void           AppendSelectedCorrelatedBANFFSystematics   (FitParamList * fp, int start);

   //
   // Private data members
   //

   t2kana::conv::PhysParamId_t WhichParam(int id) const; // array slot --> PhysParamId_t

   map<t2kana::conv::PhysParamId_t, int> fIdMap;

   string   fParListName;
   int      fSystMode; // used in FitInputs to make deciding on what to include in the chisq penalty faster. 0 no syst, 1 all syst, 2 selected syst, ...

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
   double * fMigradError;  //[fNParams]
   double * fHesseError;   //[fNParams]
   double * fPlusError;    //[fNParams]
   double * fMinusError;   //[fNParams]

   bool fMigradErrorsSet;
   bool fHesseErrorsSet;
   bool fMinosErrorsSet;

   mutable bool fShowCurrValuesOnly;
   mutable bool fShowBestFitValuesOnly;
};

} // t2kana
  
#endif


