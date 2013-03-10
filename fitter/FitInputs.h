//____________________________________________________________________________
/*!

\class    t2kana::FitInputs

\brief    Bundles fitter inputs into a single TObject-derived class so that 
          they can be passed to the minuit fitting function via 
          TMinuit::Set/GetObjectFit();

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 20, 2010

*/
//____________________________________________________________________________

#ifndef _t2kana_fit_inputs_h_
#define _t2kana_fit_inputs_h_

#include <Minuit2/FCNBase.h>
#include <Math/IFunction.h>
#include "Math/FitMethodFunction.h"
#include "HFitInterface.h"

#include "FitParamList.h"
#include "Conventions.h"

namespace t2kana {
  
class SuperkObservation;
class SuperkExpectation;
 
class FitInputs;
  
 //class FitFunctor: public ROOT::Math::IBaseFunctionMultiDim 
class FitFunctor: public ROOT::Math::IMultiGenFunction
{
public:
  FitFunctor(FitInputs * inputs);
  ~FitFunctor();
  double operator() (const double * x) const;
  double DoEval(const double * x) const;
  //ROOT::Math::IBaseFunctionMultiDim * Clone(void) const;
  ROOT::Math::IMultiGenFunction * Clone(void) const;
  unsigned int	NDim(void) const;
 private:
  FitInputs * fFitInputs;
};
  
class FitMethodFunc : public ROOT::Math::FitMethodFunction 
{
 public:
  //(ROOT::Math::FitMethodFunction & f, FitInputs * inputs);//(int dim, int npoints, FitInputs * inputs);
 FitMethodFunc( int dim, int npoints,  FitInputs * inputs) :
  ROOT::Math::FitMethodFunction(dim,npoints),
    fFitInputs2(inputs) {}
  ~FitMethodFunc();
  double operator() (const double * x) const;
  ROOT::Math::IMultiGenFunction * Clone(void) const{
    return new FitMethodFunc(*this);
  }
  virtual double DataElement(const double* x, unsigned int i, double* g = 0) const;
  unsigned int   NDim(void) const;
  unsigned int   NCalls(void) const;
  unsigned int   NPoints(void) const;
  virtual void   ResetNCalls(void) const;
  virtual void   UpdateNCalls(void) const;
  virtual  ROOT::Math::BasicFitMethodFunction<ROOT::Math::IBaseFunctionMultiDim>::Type_t Type(void) const ;
  
 private:
  FitInputs * fFitInputs2;
  double DoEval(const double * x) const;
  mutable unsigned int fNCalls;
  
  //int fNDim;
  //int fNPoints;
  //ROOT::Math::FitMethodFunction * fChi2;
};
  
class FitInputs: public ROOT::Minuit2::FCNBase, public TObject
{
public :
  FitInputs(SuperkObservation * o, SuperkExpectation * e, FitParamList * p, string opt = "poisson", double errdef = 1, bool goodnessoffit = false);
  ~FitInputs();
  
  SuperkObservation & Observation(void);
  SuperkExpectation & Expectation(void);
  FitParamList &      FitParams  (void);
  
  // set the chisq definition 5 bins vs 73 bins
  void SetChisqDefinition (bool goodnessoffit) { fGoodnessOfFit = goodnessoffit; }
  
  // judge the agreement between the input data and the current prediction
  double Chisquare (Option_t * opt = "poisson", conv::SuperkSample_t s = conv::kSKsUndefined) const;
  
  // get chisquare for the best fit prediction
  double BestFitChisquare (Option_t * opt = "poisson");
  
  // Set the chisquare type (e.g. poisson) to be used in the fit
  void SetChisquareFitOption(Option_t * opt) { fOpt = opt; }
  
  // Minuit2: Minimum chi2 fit (delta chi2 = 1 for 1 sigma error)
  double Up          (void) const    { return fErrDef;   }
  double ErrorDef    (void) const    { return Up();   }
  void   SetErrorDef (double errdef) { fErrDef = errdef; }
  
  // Minuit2: Returns chi^2
  double operator() (const std::vector<double>& x) const;
  
  // GSL
  // ROOT::Math::IBaseFunctionMultiDim * AsFunctor(void) { return new FitFunctor(this); }
  ROOT::Math::IMultiGenFunction * AsFunctor(void) { return new FitFunctor(this); }
  ROOT::Math::FitMethodFunction * AsFitMethodFunction(void) { return new FitMethodFunc(this->GetNDim(),this->GetNPoints(),this); }
  
  int GetNDim(void) const {return fNDim;}
  int GetNPoints (void) const {return fNPoints;}
  
private:
  
  double ChisquarePenaltyTerm(void) const;
  
  SuperkObservation *         fObs;
  mutable SuperkExpectation * fExp;
  mutable FitParamList *      fFitParams;
  string              fOpt;
  double              fErrDef;
  bool                fAllBANFFFit;
  bool                fGoodnessOfFit;
  int                 fNPoints;
  int                 fNDim;
  
  ClassDef(FitInputs, 4)
    
};
  
//
// MINUIT, GSL fit functions
//

void FitFunc (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);

} // t2kana

#endif

