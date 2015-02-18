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





#define mixing23Bar_min 0.           // first value of param1 
#define mixing23Bar_max 1.           // last value of param1 
#define nbins_tot_mixing23Bar 21     // tot # grid values for param1

#define dmsq32Bar_min 2.0E-4           // first value of param1 
#define dmsq32Bar_max 2.0E-2           // last value of param1 
#define nbins_tot_dmsq32Bar 41     // tot # grid values for param1



#define NFILES_MARG 1     // number of files where marginalization toys are stored
#define NTOYS_MARG 10     // total number of marginalization toys
//#define NBINS_TOT_PARAM1 nbins_tot_mixing23Bar




#include <string>

#include <TObject.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSpline.h>
#include <TVirtualPad.h>
#include <TMatrixTSym.h>
#include <TVectorD.h>
#include <TXMLEngine.h>

#include "AnalysisOptions.h"
#include "Enumerations.h"
#include "Conventions.h"
#include "Utilities.h"
#include "SuperkExpectation.h"
#include "SuperkObservation.h"

#include "FitCore/FitModelABC.h"
#include "UtlFitParamListFactory.h"
#include "FitCore/FitParamList.h"
#include "FitCore/PhysParamList.h"
#include "Utilities/VALORArray.h"
#include "Utilities/VALORArrayI.h"
#include "Utilities/VALORArrayF.h"
#include "FitCore/PhysHypothesis.h"

using std::string;
using std::vector;

using namespace valor;
using namespace valor::conv;

namespace valor {
  
  class FitParamList;
  class SuperkObservation;
  
class FitMarginalization: public TObject
{
 public:
  
  FitMarginalization(PhysParamId_t parfree, int ntoys_marg, int nbins_param);
  ~FitMarginalization(void);
  
  void Init(PhysParamId_t parfree, int ntoys_marg, int nbins_param);
  void CleanUp();
  
  void SetConfigFile();

  void SetFreeParam(PhysParamId_t parfree);
  
  void SetFitParamList(char * fit_param_list);
  FitParamList * GetFitParamList(void) {return fparm;};
  
  void PrintPhysHypothesis(void);
  void PrintFitParamList(void);

  TH1F * GetChisqMarg (void) { return hchisq_marg; };

  void MakeToysMarginalization(int tweak_systematics,SuperkExpectation *skexp,int min_bin_sub_param,int max_bin_sub_param);

  void WriteToys2Binary(int idfile,char *tag);
  void ReadToysFromBinary(int idfile, char *tag);
  
  void ResetChisqMarg(int min_bin_sub_param,int max_bin_sub_param);
  void FillChisqMarg(SuperkObservation *skobs_data);
  
  void Minimize(TH1F *hInChisqMarg);
  void Minimize(void);

  double GetChisqBestFit(void){ return fchisq_bf;};
  double GetParamBestFit(void){ return fparam_bf;};
  
  void SetSubRange();

  void PrintBestFit(void);

  VALORArrayF * GetToyMargErecoVArray(SuperkSample_t s);



  void SetMargOscParam();

  void SetPrior(PhysParamId_t physpar,string priorname,double par1,double par2);
  void SetPrior(PhysParamId_t physpar,string priorname,TH1F *hprior);
  string GetPrior(PhysParamId_t physpar);

  double Randomize(PhysParamId_t physpar);
  

 private:
  
  TH1F *hchisq_marg; 
  
  string fnamefreeparam;
  PhysParamId_t fPhysParamFree;

  double fchisq_bf;
  double fparam_bf;

  int fnbins_param;
  int fntoys_marg;
  int fnerecobins_RHCNumu;
  
  int fmin_bin_sub_param; // first bin of sub-range to scan (for partial toys and chisq marg)
  int fmax_bin_sub_param; // last bin of sub-range to scan (for partial toys and chisq marg)

  double fmin_sub_param; // first value of sub-range to scan (for partial toys and chisq marg)
  double fmax_sub_param; // last value of sub-range to scan (for partial toys and chisq marg)

  double fmin_param;
  double fmax_param;
  double fstep_param;
  
  VALORArrayF *flikeVArray;
  VALORArrayF *fChisqMargVArray;
  VALORArrayF *fToyMargErecoVArray_RHCNumu;

  FitParamList * fparm;
  PhysHypothesis_t fhypothesis;

  ClassDef(FitMarginalization, 1)
    };
 
} // valor
#endif
