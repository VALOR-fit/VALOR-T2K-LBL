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
#include <TH1D.h> // NEW PRECISION
#include <TH2D.h> // NEW PRECISION
#include <TTree.h>
#include <TFile.h>
#include <TStopwatch.h>

#include "Enumerations.h"

#include "valor/VArray/VALORArrayI.h"
#include "valor/VArray/VALORArrayF.h"
#include "valor/VArray/VALORArrayD.h" // NEW PRECISION

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

// NEW BETABAR
#define betaBar_min 0.
#define betaBar_max 20.
//

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
  virtual ~FitMarginalization(void);

  void Init(SuperkExpectation *skexp);

  void CleanUp();
  
  //void SetConfigFile();

  void SetPhysParFree(conv::PhysParamId_t *parfree,int ntoys_marg,int *nbins_param,double *minpar,double *maxpar);
  void SetFreeParam(int ParID,conv::PhysParamId_t *parfree,double *minpar,double *maxpar);
  void SetFullRange();
  void SetSubRange(int *submin,int *submax);
  
  int GetFirstBinInSubRange(int par=0) const {return fmin_bin_sub_param[par];}
  int GetLastBinInSubRange (int par=0) const {return fmax_bin_sub_param[par];}

  int GetFirstBinInFullRange(int par=0) const {return fmin_bin_param[par];}
  int GetLastBinInFullRange (int par=0) const {return fmax_bin_param[par];}
  
  double GetBinParameterValue(int par, int ibin) const {return fmin_param[par] + ibin*fstep_param[par];} // TO CHECK
   
  void SetFitParamList(const char * fit_param_list);

  FitParamList * GetFitParamList(void) {return fparm;};
  // NEW MARG MH
  FitParamList * GetFitParamList_NH(void) {return fparm_nh;};
  FitParamList * GetFitParamList_IH(void) {return fparm_ih;};
  //
  
  //Store parameter variations in a root file
  void SaveTreeOfToyParamLists(TFile* file);
  
  const FitParamList * GetFitParamList(void) const {return fparm;};
  // NEW MARG MH
  const FitParamList * GetFitParamList_NH(void) const {return fparm_nh;};
  const FitParamList * GetFitParamList_IH(void) const {return fparm_ih;};
  //

  // NEW FHC+RHC JOINT 
  double GetMinOscParam(conv::PhysParamId_t physpar) const;
  double GetMaxOscParam(conv::PhysParamId_t physpar) const;
  //


  // NEW 2D TEMPLATES (VARRAY)nuosc3fbar_numudissap_allsyst
  double GetNormPOT(conv::SuperkSample_t s) const;
  double GetPOT(conv::SuperkSample_t s) const;
  //

  void PrintPhysHypothesis(void);

  bool IsFit2D(){ return fIsFit2D;}

  void MakeToysMarginalization(int tweak_systematics,SuperkExpectation *skexp);

  void WriteToys2Binary(int idfile, const char *tag);
  void WriteToys2BinaryMerge(const char *tag);
  void ReadToysFromBinary(int idfile, const char *tag);

  // NEW SPECTRA MARG
  void WriteLikeToys2Binary(const char *tag); 
  void ReadLikeToysFromBinary(const char *tag);
  void SetStoreLikeWeights(bool doStoreLikeWeights);
  //
  // NEW OSC TREE
  void WriteOscToys2Tree(int idfile, const char *tag); 
  void InitTreeOscToys();
  //

  void ResetChisqMarg(int *min_bin_sub_param,int *max_bin_sub_param);

  void BuildDeltaChisqMarg(void);

  void FillChisqMarg(SuperkObservation *skobs_data,bool readbinary = false); // Loop over the free params and call CalcChisqMargPt()
    

  // NEW PRECISION
  //double CalcChisqMargPt(SuperkObservation *skobs_data,int iterpar0,int iterpar1, bool allow_cached = true);
  double CalcChisqMargPtRead(SuperkObservation *skobs_data,int iterpar0,int iterpar1,bool readbinary = false,bool allow_cached = true); // Use the map of the binary files
  //

  // NEW PRECISION
  //double GetChisqMargTrue(double valpar0,double valpar1,SuperkObservation *skobs);
  double GetChisqMargTrueRead(double valpar0,double valpar1,SuperkObservation *skobs,bool readbinary = false); // Use the map of the binary files
  //
  
  int GetNCallsCalcChisqMargPt(void) const { return fNCallsCalcChisqMargPt;};

  void CalcBestFitErrors(SuperkObservation *skobs,string method="MINOS",bool readbinary=false); // calculate the best-fit errors
  void CalcBestFitErrorsMINOS(SuperkObservation *skobs,bool readbinary = false); // calculate the best-fit errors
  void CalcBestFitErrorsHESSE(SuperkObservation *skobs,bool readbinary = false,int Nh = 1); // calculate the best-fit errors

  void MinimizeHisto(TH1D *hInChisqMarg); // NEW PRECISION
  void MinimizeHisto(TH2D *hInChisqMarg); // NEW PRECISION
  void MinimizeHisto(void);
  
  double GetChisqBestFit(void) const      { return fchisq_bf;};
  double GetParamBestFit(int ParId) const { return fparam_bf[ParId];};
  double GetParamBestFitErrNeg(int ParId) const { return fparam_bf_errneg[ParId];};
  double GetParamBestFitErrPos(int ParId) const { return fparam_bf_errpos[ParId];};
  double GetParamBestFitErrSym(int ParId) const { return fparam_bf_errsym[ParId];};
  
  void PrintBestFit(void);
  
  void SetChisqBestFit(double chisq) { fchisq_bf = chisq;};
  void SetParamBestFit(int ParId,double val) { fparam_bf[ParId] = val;};
  void SetParamBestFitErrNeg(int ParId,double val) { fparam_bf_errneg[ParId] = val;};
  void SetParamBestFitErrPos(int ParId,double val) { fparam_bf_errpos[ParId] = val;};
  void SetParamBestFitErrSym(int ParId,double val) { fparam_bf_errsym[ParId] = val;};
  

  VALORArrayF * GetToyMargVArray()             {return fToyMargVArray;};
  const VALORArrayF * GetToyMargVArray() const {return fToyMargVArray;};

  // NEW SPECTRA MARG
  double GetLikeToyRead(double valpar0,double valpar1,int iexperiment,bool readbinary = false); // Use the map of the binary files
  //

  // NEW 2D TEMPLATES (VARRAY)
  VALORArrayF * GetToyMargVArraySample(conv::SuperkSample_t s);
  VALORArrayF * GetToyMargVArraySample_Clone(conv::SuperkSample_t s);
  //

  TH1F * GetToyMargHisto(conv::SuperkSample_t s,int itoy,int ifreepar0,int ifreepar1); 
  
  // NEW 2D TEMPLATES
  TH2F * GetToyMargHisto2D(conv::SuperkSample_t s,int itoy,int ifreepar0,int ifreepar1); 
  //
  
  // THIS ISN"T IMPLEMENTED ANYMORE! void SetMargOscParam();
  void SetMargOscParam(FitParamList * infparm);
  
  // THIS ISN"T IMPLEMENTED ANYMORE! void SetPrior(conv::PhysParamId_t physpar, std::string priorname,double par1,double par2);
  // THIS ISN"T IMPLEMENTED ANYMORE! void SetPrior(conv::PhysParamId_t physpar, std::string priorname,TH1D *hprior);
  // THIS ISN"T IMPLEMENTED ANYMORE! void SetPrior(conv::PhysParamId_t physpar, std::string priorname,TH2D *hprior);
  
  void SetPrior_NH(conv::PhysParamId_t physpar, std::string priorname,double par1,double par2);
  void SetPrior_NH(conv::PhysParamId_t physpar, std::string priorname,TH2D *hprior);
  void SetPrior_NH(conv::PhysParamId_t physpar, std::string priorname,TH1D *hprior);
  void SetPrior_IH(conv::PhysParamId_t physpar, std::string priorname,double par1,double par2);
  void SetPrior_IH(conv::PhysParamId_t physpar, std::string priorname,TH1D *hprior);
  void SetPrior_IH(conv::PhysParamId_t physpar, std::string priorname,TH2D *hprior);
  //
  
  // THIS ISN"T IMPLEMENTED ANYMORE! std::string GetPrior(conv::PhysParamId_t physpar) const;
  // NEW MARG MH
  std::string GetPrior_NH(conv::PhysParamId_t physpar) const;
  std::string GetPrior_IH(conv::PhysParamId_t physpar) const;
  //
  
  void BuildDefaultPriors(void);

  // THIS ISN"T IMPLEMENTED ANYMORE! double Randomize(conv::PhysParamId_t physpar) const;
  // NEW MARG MH
  double Randomize(conv::PhysParamId_t physpar, FitParamList * infparm) const;
  //
  
  // NEW MARG MH 2
  int RandomizeMH(const double probNH) const; // input is the prob for NH. P(IH) is 1-P(NH) 
  void SetProbNormalHierarchy(double probnh);
  double GetProbNormalHierarchy(void) const { return fProbNormalHierarchy;};

  TH1D *GetHistoPriorMH(void) const {return fhPriorMH;};
  TH1D *GetHistoToyMH(void) const {return fhToyMH;};
  //


  bool RethrowOscUnphysical(conv::PhysParamId_t p,double oscval) const; // NEW 
  
  void SetMargMassHierarchy(bool domargMH);  
  void SetRethrowSystUnphysical(bool dorethrow);
  
  double OscPhysicalRegion(conv::PhysParamId_t p,double oscval);

  void MergeInputFiles();

  void SetNSKsamples();
  int GetNSKsamples ()  const {return fNSKsamples;};
  conv::SuperkSample_t GetSKSamples (int isksample) const ; 


  // NEW FHC+RHC JOINT
  bool UseSKSample(int isksample) const;
  bool UseSKSample(conv::SuperkSample_t s) const;
  //



  void SetSaveParamLists(bool save_lists = true);
  PhysParamList * const * GetParamListArray() const;
  TObjArray * GetToyParamThrowHistograms(int nbins = -1) const; // Array needs deleting after call, histograms don't.

  // Histogram minimization
  void MinimizeChisqIntervalRead(SuperkObservation *skobs_data,bool readbinary = false);
  void MinimizeChisqStepRead(SuperkObservation *skobs_data,bool readbinary = false);  // MAP  
  void SetMinimIntervalRead(double XvalFirst,double XvalLast,double YvalFirst,double YvalLast,bool readbinary = false); // MAP  
  void SetMaxIterationSingleFit(int maxiter); // by default is the tot # of histo bins
  void SetMaxFitNumber(int nfit);

  // Step (descent gradient)
  void GridStartBinRead(SuperkObservation *skobs_data,bool readbinary = false); // MAP
  void SetGridStartRange(double minX,double maxX,double minY,double maxY);
  //

  void SetStartBin(double start_X,double start_Y); // if <0 --> throw random
  void SetNumStartBin(int Nthrows_X,int Nthrows_Y);
  void SetExtractStartPt(string doextract);
  
  // Interval (descent gradient)
  void SetMinimUncorrXY(bool flaguncorr);
  void SetBorderDoFit(bool dofit);
  void SetBorderMaxFitNumber(int nfit);

  TH1D *GetChisqMarg1D (void){ return fhchisq_marg_1D;}; // NEW PRECISION
  TH2D *GetChisqMarg2D (void){ return fhchisq_marg_2D;}; // NEW PRECISION
  
  // NEW PENALTY TERM
  double CalcPenaltyTerm (const double par0, const double par1=0.);
  void SetPenaltyTerm (const double mean0, const double sigma0, const double mean1=0., const double sigma1=0.);
  void UsePenaltyTerm (const bool usepenalty0, const bool usepenalty1);
  const TH1D *GetPenaltyTerm1D (void){ return fhpenaltyterm_1D;}; 
  const TH2D *GetPenaltyTerm2D (void){ return fhpenaltyterm_2D;}; 
  //

  TH1D *GetDeltaChisqMarg1D (void){ return fhdeltachisq_marg_1D;}; // NEW PRECISION 
  TH2D *GetDeltaChisqMarg2D (void){ return fhdeltachisq_marg_2D;}; // NEW PRECISION

  TH1D *GetChisqStep_1D(void){ return fhchisq_marg_1D_Step;}; // NEW PRECISION 
  TH2D *GetChisqStep_2D(void){ return fhchisq_marg_2D_Step;}; // NEW PRECISION

  int GetNumFit(void){ return fNFits;};
  int GetNumIterations(void){ return fNiterations;};
  int GetFitStatus(void){ return fFitStatus;};

  double GetInterpolationX(double X1,double Y1,double X2, double Y2,double Yin);

  // MAP
  TH1D *GetMapChisqRangeFileID_1D(void){ return fhmap_chisqrange_fileid_1D;};  
  TH2D *GetMapChisqRangeFileID_2D(void){ return fhmap_chisqrange_fileid_2D;};  
  void SetTagFileName(const char *tag);
  void SetFileID(int nfiles);
  int GetFileIDRead(int binpar0,int binpar1); // need histogram bin number as input
  int GetCurrentFileID(void){ return fCurrentFileIDminim;};
  //
  
  void SetTestStatistic(string teststatistic); // "rateonly" or "rate+shape". Always use the log-likelihood ratio
  void SetSeed(UInt_t seed) {fSeed = seed;};
  
  //NuEBar app
  double FillExpectedNuEBarAppHisto(SuperkObservation *skobs_data,int iterpar0,bool readbinary);
  TH1D* GetNuEBarHisto(int i){return fNuEAppHisto[i];};
  // Histo below is required for p-value statistical uncertainty
  TH1D* GetNuEBarHistoUnscaled(int i){return fNuEAppHistoUnscaled[i];};

private:
  void SetNoOscillation(FitParamList* fparm);
  bool TweakSystematics(SuperkExpectation* skexp, int tweak_systematics);
  bool ProduceToys(SuperkObservation* skobs, SuperkExpectation* skexp,
        FitParamList* fparm, int iexperiment);

private:

  //NuEBar app histogram one for each value of Beta (is only initialised and 
  //filled if beta is the only parameter of intereset 
  TH1D* fNuEAppHisto[2];
  // histo below is required for p-value statistical uncertainty
  TH1D* fNuEAppHistoUnscaled[2];

  // NEW PRECISION
  TH1D *fhchisq_marg_1D; 
  TH2D *fhchisq_marg_2D;  
  TH1D *fhdeltachisq_marg_1D; 
  TH2D *fhdeltachisq_marg_2D;  
  //  

  // NEW PENALTY TERM
  TH1D *fhpenaltyterm_1D; 
  TH2D *fhpenaltyterm_2D; 
  double fpenalty_mean0;
  double fpenalty_mean1;
  double fpenalty_sigma0;
  double fpenalty_sigma1;
  bool fusepenalty0;
  bool fusepenalty1;
  //

  // MAP
  TH1D *fhmap_chisqrange_fileid_1D;
  TH2D *fhmap_chisqrange_fileid_2D;

  int fCurrentFileIDminim;

  string ftagfile;       // tag of the binary file name
  string *fnamefile; // all the binary file names (done in SetFileID), read during minimization

  double fGridStartMinBin_X;
  double fGridStartMaxBin_X;
  double fGridStartMinBin_Y;
  double fGridStartMaxBin_Y;

  double fGridStartMinVal_X;
  double fGridStartMaxVal_X;
  double fGridStartMinVal_Y;
  double fGridStartMaxVal_Y;
 
  //

  std::string fnamefreeparam[2]; // Use param 1 only for 2D fits

  conv::PhysParamId_t fPhysParamFree[2]; 

  bool fIsFit2D;

  bool fdoMargMassHierarchy; // marginalize MH. Not by default     
  bool fRethrowSystUnphysical; // rethrow toys with unphysical systematics. True by default

  char fnamefilebinary[200];

  // NEW SPECTRA MARG
  char fnamelikefilebinary[200]; // to store toy likelihood from the fit
  //

  double fchisq_bf;
  double fparam_bf[2]; 
  double fparam_bf_errpos[2]; 
  double fparam_bf_errneg[2]; 
  double fparam_bf_errsym[2]; 
  
  UInt_t fSeed;

  int fnbins_param[2];
  int fnbins_sub_param[2]; 
  
  int fntoys_marg;
  
  int fmin_bin_param[2]; // first bin of sub-range to scan (for partial toys and chisq marg)
  int fmax_bin_param[2]; // last bin of sub-range to scan (for partial toys and chisq marg)

  int fmin_bin_sub_param[2]; // first bin of sub-range to scan (for partial toys and chisq marg)
  int fmax_bin_sub_param[2]; // last bin of sub-range to scan (for partial toys and chisq marg)

  double fmin_sub_param[2];  // first value of sub-range to scan (for partial toys and chisq marg)
  double fmax_sub_param[2];  // last value of sub-range to scan (for partial toys and chisq marg)

  double fmin_param[2];  
  double fmax_param[2];  
  double fstep_param[2]; 

  int fNCallsCalcChisqMargPt; // # of calls of CalcChisqMargPt
  
  // NEW PRECISION
  //VALORArrayD *flikeVArray; //UNUSED
  VALORArrayD *fChisqMargVArray;
  //

  VALORArrayD *flikeToyVArray; // Fit likelihood of each toy  // NEW SPECTRA MARG

  
  int fNSKsamples;
  VALORArrayI *fnerecobins;
  VALORArrayI *fsksamples; // pointer of SuperkSample_t doesn't work. Then use kSKSamples[]
  
  // NEW 2D TEMPLATES (VARRAY)
  VALORArrayF *fToyMargVArray; // TO REMOVE
  //VALORArrayF *fToyMargVArray_Clone; // TO REMOVE
  VALORArrayF *fToyMargVArray_1Rmu;
  VALORArrayF *fToyMargVArray_1Re;
  VALORArrayF *fToyMargVArray_RHC1Rmu;
  VALORArrayF *fToyMargVArray_RHC1Re;

  VALORArrayF *fToyMargVArray_MultiRe;
  
  VALORArrayF *fToyMargVArray_Clone_1Rmu;
  VALORArrayF *fToyMargVArray_Clone_1Re;
  VALORArrayF *fToyMargVArray_Clone_RHC1Rmu;
  VALORArrayF *fToyMargVArray_Clone_RHC1Re;

  VALORArrayF *fToyMargVArray_Clone_MultiRe;
  //
  
  TH1F *fhToyMargCurr[conv::kns]; // a histogram for each sk sample. Filled using fsksamples and fNSKsamples, not as usual by looping over all the samples and skipping the undefined ones 
  int fMaxSampleErecoBins;

  // NEW 2D TEMPLATES
  VALORArrayI *fnq2bins; 
  TH2F *fhToyMargCurr2D[conv::kns]; // a histogram for each sk sample. Filled using fsksamples and fNSKsamples, not as usual by looping over all the samples and skipping the undefined ones 
  int fMaxSampleQ2Bins;
  //
  
  bool fMergeInputFiles_Called;
  int fNmergedFiles;

  FitParamList * fparm;
  // NEW MARG MH
  FitParamList * fparm_nh;
  FitParamList * fparm_ih;
  //

  // NEW MARG MH 2
  double fProbNormalHierarchy; // Probability for Normal Hierarchy (for not flat prior)
  //
  
  conv::PhysHypothesis_t fhypothesis;
  
  bool             fShouldSaveToyParamLists;
  PhysParamList ** fSavedToyParamLists; //[fntoys_marg]

  SuperkExpectation *fskexp;

  // Histogram minimization
  
  // NEW PRECISION
  TH1D *fhchisq_marg_1D_Step; // Get track of the minimization steps
  TH2D *fhchisq_marg_2D_Step; // Get track of the minimization steps
  //
  
  int fMaxIterationSingleFit; // Partial # of iterations before restarting the fit
  int fMaxFitNumber;        // Maximum # of fits

  int fBorderDoFit; // if the best-fit is on the border, the fit is done again (default = true)
  int fBorderMaxFitNumber; // if the best-fit is on the border, the fit is redone fBorderMaxFitNumber times (default = # of borders = 2 x NDim)

  bool fUncorrXY; // if minimum found in Y, keep it fixed and vary only X (viceversa)
  
  // Set the interval where to find the minimum (bin numbers)
  double fFirstMinimVal_X;
  double fLastMinimVal_X;
  double fFirstMinimVal_Y;
  double fLastMinimVal_Y;

  int fFirstMinimBin_X;
  int fLastMinimBin_X;
  int fFirstMinimBin_Y;
  int fLastMinimBin_Y;
  
  // Histogram minimization with steps (descent gradient)
  int fStep;
  int fStartBin_X;
  int fStartBin_Y;
  int fNumStartBin_X; 
  int fNumStartBin_Y;

  string fExtractStartPt;
  
  int fNFits;
  int fNiterations;
  int fFitStatus; // 1 --> OK, 0 --> it didn't converge

  double fPOT[conv::kns]; // POT of each sk sample to print in the output text file with marg toys
  double fNormPOT[conv::kns]; // POT of each sk sample to print in the output text file with marg toys. Normalized to 1E+20
 
  string fteststatistic; // "rateonly" or "rate+shape"

  // NEW MARG MH 2
  TH1D *fhPriorMH;
  TH1D *fhToyMH;
  //

  

  // NEW OSC TREE
  bool fStoreLikeWeights;
  char fnameoscfiletree[200]; // to store toy likelihood from the fit
  TTree *treetoys;
  TFile *filetreetoys;
  //
  
  ClassDef(FitMarginalization,9)
};
 
} // valor
#endif
