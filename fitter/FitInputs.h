//____________________________________________________________________________
/*!

\class     t2kana::FitInputs

\brief     Bundles fitter inputs into a single TObject-derived class so that 
           they can be passed to the minuit fitting function via 
           TMinuit::Set/GetObjectFit();

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Laboratory

\created   October 20, 2010

\copyright (c) 2010-2014, The VALOR group

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
  
class FitDataABC;
class FitModelABC; 
class FitInputs;
//____________________________________________________________________________  
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
//____________________________________________________________________________
class FitMethodFunc : public ROOT::Math::FitMethodFunction 
{
 public:
 FitMethodFunc( int dim, int npoints,  FitInputs * inputs) :
  ROOT::Math::FitMethodFunction(dim,npoints),
    fFitInputs(inputs) {}
  ~FitMethodFunc();
  double operator() (const double * x) const;
  ROOT::Math::IMultiGenFunction * Clone(void) const{
    return new FitMethodFunc(*this);
  }
  
  double DataElement(const double* x, unsigned int i, double* g = 0) const;
  
  virtual  ROOT::Math::BasicFitMethodFunction<ROOT::Math::IBaseFunctionMultiDim>::Type_t Type(void) const ;
  
 private:
  mutable FitInputs * fFitInputs;
  virtual double DoEval(const double * x) const;
  
};
//____________________________________________________________________________
class FitInputs: public ROOT::Minuit2::FCNBase, public TObject
{
public :
  FitInputs(FitDataABC * o, FitModelABC * e, FitParamList * p, string opt = "poisson", double errdef = 1, bool goodnessoffit = false);
  ~FitInputs();
  
  FitDataABC &    Observation (void);
  FitModelABC &   Expectation (void);
  FitParamList &  FitParams   (void);
  void UpdateParamList(FitParamList * fpar);
  
  double Statistic(Option_t * stat_opt = "likelihood-ratio") const;

  // set the chisq definition 5 bins vs 73 bins
  void SetChisqDefinition (bool goodnessoffit) { fGoodnessOfFit = goodnessoffit; }
  
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
  
  
  // Old interface
  double BestFitChisquare(void) const { return fFitParams->GetBestFitStatistic(); }
  // Set the chisquare type (e.g. poisson) to be used in the fit
  void SetChisquareFitOption(Option_t * opt) { fOpt = opt; }

  double ChisquarePenaltyTerm(void) const;
  // judge the agreement between the input data and the current prediction
  double Chisquare (Option_t * opt,      conv::SuperkSample_t s = conv::kSKsUndefined, int point = -1) const;
  double Chisquare (int        mode = 1, conv::SuperkSample_t s = conv::kSKsUndefined, int point = -1) const; //mode 1 is poisson, mode 2 is gauss

private:
  
  FitDataABC *            fObs;
  mutable FitModelABC *   fExp;
  mutable FitParamList *  fFitParams;
  string                  fOpt;
  int                     fMode;
  double                  fErrDef;
  bool                    fAllBANFFFit;
  bool                    fAll1ReEff;
  bool                    fGoodnessOfFit;
  int                     fNPoints;
  int                     fNDim;
  
  ClassDef(FitInputs, 5)    
};
  
//
// MINUIT, GSL fit functions
//

void FitFunc (Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t iflag);

} // t2kana

#endif

