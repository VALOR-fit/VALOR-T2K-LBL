//____________________________________________________________________________
/*!

\class    valor::FitMarginalization

\brief    Code to perform the fit marginalizing the nuisance parameters

\author   Davide Sgalaberna <davide.sgalaberna \at cern.ch>
          ETH Zurich
          
\created  Feb 16, 2015

\copyright (c) 2010-2015, VaLOR group.

*/
//____________________________________________________________________________

#ifndef _valor_fitmarginalization__h_
#define _valor_fitmarginalization__h_


#include <string>

#include <TH1F.h>
#include <TH2F.h>

#include "Enumerations.h"

#include "Utilities/VALORArrayI.h"
#include "Utilities/VALORArrayF.h"

// NuBar mode

#define mixing13Bar_min 0.          
#define mixing13Bar_max 1.           

#define mixing23Bar_min 0.           // both sinsq2theta and sinsqtheta
#define mixing23Bar_max 1.           

#define mixing12Bar_min 0.           // both sinsq2theta and sinsqtheta
#define mixing12Bar_max 1.           

#define dmsq32Bar_min valor::conv::kMinFitDm23sqBar  
#define dmsq32Bar_max valor::conv::kMaxFitDm23sqBar 

#define dmsq12Bar_min valor::conv::kMinFitDm12sqBar  
#define dmsq12Bar_max valor::conv::kMaxFitDm12sqBar 

#define deltacpBar_min -3.14159265358979312
#define deltacpBar_max +3.14159265358979312

// Nu mode

#define mixing13_min 0.          
#define mixing13_max 1.           

#define mixing23_min 0.           // both sinsq2theta and sinsqtheta
#define mixing23_max 1.           

#define dmsq32_min valor::conv::kMinFitDm23sq  
#define dmsq32_max valor::conv::kMaxFitDm23sq 

#define mixing12_min 0.           // both sinsq2theta and sinsqtheta
#define mixing12_max 1.           

#define dmsq12_min valor::conv::kMinFitDm12sq  
#define dmsq12_max valor::conv::kMaxFitDm12sq 

#define deltacp_min -3.14159265358979312
#define deltacp_max +3.14159265358979312

namespace valor {
  
  class FitParamList;
  class PhysParamList;
  class SuperkObservation;
  class SuperkExpectation;
  
class FitMarginalization
{
 public:

  FitMarginalization();
  ~FitMarginalization(void);

  void Init();

  void CleanUp();
  
  //void SetConfigFile();

  void SetPhysParFree(conv::PhysParamId_t *parfree,int ntoys_marg,int *nbins_param,double *minpar,double *maxpar);
  void SetFreeParam(int ParID,conv::PhysParamId_t *parfree,double *minpar,double *maxpar);
  void SetSubRange(int *submin,int *submax);

   
  void SetFitParamList(const char * fit_param_list);
  FitParamList * GetFitParamList(void) {return fparm;};
  
  void PrintPhysHypothesis(void);
  //void PrintFitParamList(void);

  TH1F * GetChisqMarg   (void);
  TH2F * GetChisqMarg2D (void);

  bool IsFit2D(){ return fIsFit2D;}

  void MakeToysMarginalization(int tweak_systematics,SuperkExpectation *skexp);

  void WriteToys2Binary(int idfile, const char *tag);
  void WriteToys2BinaryMerge(const char *tag);

  void ReadToysFromBinary(int idfile, const char *tag);
  
  void ResetChisqMarg(int *min_bin_sub_param,int *max_bin_sub_param);

  void FillChisqMarg(SuperkObservation *skobs_data); // Loop over the free params and call CalcChisqMargPt()
  float CalcChisqMargPt(SuperkObservation *skobs_data,int iterpar0,int iterpar1);
  
  int GetNCallsCalcChisqMargPt(void) { return fNCallsCalcChisqMargPt;};

  void Minimize(TH1F *hInChisqMarg);
  void Minimize(TH2F *hInChisqMarg);
  void Minimize(void);

  float GetChisqBestFit(void){ return fchisq_bf;};
  float GetParamBestFit(int ParId){ return fparam_bf[ParId];};
  
  void PrintBestFit(void);

  VALORArrayF * GetToyMargErecoVArray() {return fToyMargErecoVArray;};
  TH1F * GetToyMargHisto(conv::SuperkSample_t s,int itoy,int ifreepar0,int ifreepar1,int iMassHierarchy);

  void SetMargOscParam();

  void SetPrior(conv::PhysParamId_t physpar, std::string priorname,double par1,double par2);
  void SetPrior(conv::PhysParamId_t physpar, std::string priorname,TH1F *hprior);
  std::string GetPrior(conv::PhysParamId_t physpar);

  double Randomize(conv::PhysParamId_t physpar);
  bool RethrowOscUnphysical(conv::PhysParamId_t p,double oscval); // NEW 
  
  void SetMargMassHierarchy(bool domargMH);  
  void SetRethrowSystUnphysical(bool dorethrow);
  
  void MergeInputFiles();

  void SetNSKsamples();
  int GetNSKsamples(){return fNSKsamples;};
  conv::SuperkSample_t GetSKSamples(int isksample); 
  
  void SetSaveParamLists(bool save_lists = true);
  PhysParamList * const * GetParamListArray() const;
  TObjArray * GetToyParamThrowHistograms(int nbins = -1) const; // Array needs deleting after call, histograms don't.





  // Histogram minimization

  void MinimizeChisqHisto(SuperkObservation *skobs_data);
  //void MinimizeChisqHisto(); // Minimize the filled histogram
  void SetMinimInterval(double XvalFirst,double XvalLast,double YvalFirst,double YvalLast);
  void SetMinimPrecision(double eps);
  void SetMinimUncorrXY(bool flaguncorr);
  void SetMaxIterationSingleFit(int maxiter); // by default is the tot # of histo bins
  void SetMaxFitNumber(int nfit);
  void SetBorderDoFit(bool dofit);
  void SetBorderMaxFitNumber(int nfit);
  TH1F *GetChisqStep_1D(void){ return hchisq_marg_1D_Step;};  
  TH2F *GetChisqStep_2D(void){ return hchisq_marg_2D_Step;};  
  //





  
 private:
  
  TH1F *hchisq_marg_1D; 
  TH2F *hchisq_marg_2D;  
  
  std::string fnamefreeparam[2]; // Use param 1 only for 2D fits

  conv::PhysParamId_t fPhysParamFree[2]; 

  bool fIsFit2D;

  bool doMargMassHierarchy; // marginalize MH. Not by default     
  bool fRethrowSystUnphysical; // rethrow toys with unphysical systematics. True by default

  char fnamefilebinary[200];

  double fchisq_bf;
  double fparam_bf[2]; 
  
  UInt_t fSeed;

  int fnbins_param[2]; 
  int fnbins_sub_param[2]; 
  
  int fntoys_marg;

  int fmin_bin_sub_param[2]; // first bin of sub-range to scan (for partial toys and chisq marg)
  int fmax_bin_sub_param[2]; // last bin of sub-range to scan (for partial toys and chisq marg)

  double fmin_sub_param[2];  // first value of sub-range to scan (for partial toys and chisq marg)
  double fmax_sub_param[2];  // last value of sub-range to scan (for partial toys and chisq marg)

  double fmin_param[2];  
  double fmax_param[2];  
  double fstep_param[2]; 

  int fNCallsCalcChisqMargPt; // # of calls of CalcChisqMargPt
  
  VALORArrayF *flikeVArray;
  VALORArrayF *fChisqMargVArray;
  
  int fNSKsamples;
  VALORArrayI *fnerecobins;
  VALORArrayI *fsksamples; // pointer of SuperkSample_t doesn't work. Then use kSKSamples[]
  VALORArrayF *fToyMargErecoVArray;
  VALORArrayF *fToyMargErecoVArray_Clone;
  TH1F *hToyMargCurr[conv::kns]; // a histogram for each sk sample. Filled using fsksamples and fNSKsamples, not as usal by looping over all the samples and skipping the undefined ones 
  int fMaxSampleErecoBins;

  int fnerecobins_RHCNumu;
  VALORArrayF *fToyMargErecoVArray_RHCNumu_Clone; 
  VALORArrayF *fToyMargErecoVArray_RHCNue_Clone; 
  VALORArrayF *fToyMargErecoVArray_Numu_Clone; 
  VALORArrayF *fToyMargErecoVArray_Nue_Clone; 
  bool MergeInputFiles_Called;
  int fNmergedFiles;

  TH1F *hToyMargCurr_RHCNumu;

  VALORArrayF *fToyMargErecoVArray_RHCNumu; 
  VALORArrayF *fToyMargErecoVArray_RHCNue;  
  VALORArrayF *fToyMargErecoVArray_Numu;    
  VALORArrayF *fToyMargErecoVArray_Nue;     

  FitParamList * fparm;
  conv::PhysHypothesis_t fhypothesis;
  
  bool             fShouldSaveToyParamLists;
  PhysParamList ** fSavedToyParamLists; //[fntoys_marg]







  // Histogram minimization
  
  TH1F *hchisq_marg_1D_Step; // Get track of the minimization steps
  TH2F *hchisq_marg_2D_Step; // Get track of the minimization steps
  
  int fMaxIterationSingleFit; // Partial # of iterations before restarting the fit
  int fMaxFitNumber;        // Maximum # of fits

  int fBorderDoFit; // if the best-fit is on the border, the fit is done again (default = true)
  int fBorderMaxFitNumber; // if the best-fit is on the border, the fit is redone fBorderMaxFitNumber times (default = # of borders = 2 x NDim)

  bool fUncorrXY; // if minimum found in Y, keep it fixed and vary only X (viceversa)
  double fEPS; // precision of the fit
  
  // Set the interval where to find the minimum (bin numbers)
  double fFirstVal_X;
  double fLastVal_X;
  double fFirstVal_Y;
  double fLastVal_Y;

  int fFirstBin_X;
  int fLastBin_X;
  int fFirstBin_Y;
  int fLastBin_Y;
  
  //









};
 
} // valor
#endif
