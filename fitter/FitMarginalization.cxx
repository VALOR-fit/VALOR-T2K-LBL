//________________________________________________________________________________________
/*
 
 \copyright (c) 2010-2015, VaLOR group.
 
 */
//________________________________________________________________________________________

//#define __print_oscparam_all__ // print extracted osc param for both NH and IH
//#define __print_oscparam_all_verbose__ // print extracted osc param used in osc prob
//#define __print_masshierarchy__
//#define __print_oscparam_marginalized__
//#define __print_build_map_fileid__
//#define __print_array_toy_maketoys__
//#define __print_array_toy_maketoys_verbose__
// #define __print_like_toy_fillchisq__
// #define __print_array_toy_fillchisq__
// #define __print_array_toy_fillchisq_verbose__
//#define __print_array_toy_merge__
//#define __print_systematic_variations__
//#define __print_chisqmarg__
//#define __print_chisqmarg_ptval__
//#define __print_chisqmarg_tstopwatch__
//#define __print_chisqmarg_changebinary__
// #define __readtoys_debug__
//#define __print_minimization_steps__ //
//#define __print_hesseerror__
//#define __print_hesseerror_verbose__
//#define __show_each_experiment__

//#define FILL_TREE_OSCPARAM_TOYS

#define NUMOFBINSAROUNDMINTOCHECK 3  // Check bins both on the right and left Emac

#define SCALE_PRECISION 80. // Rescale the likelihood by exp(scale) when sum over the toys, then scale back at the end (CalcChisqMargPtRead())
                            // Avoid precision problems w/ floats or double when chisq~200 and like=exp(-0.5*chisq) becomes too small
                            // floats support numbers until exp(+88). The likelihood becomes like=exp(-0.5*chisq+88) >= exp(+88)
                            // When chisq gets larger, the scaling factor can be increased, paying attention that -200 < -0.5*chisq+scale < 89
                            // Now double are used for chisq and likelihood, but the scale factor could be needed anyway

//#define __MARGINALIZE_23PARAMS_CORR__     // Extract mixing23 and dmsq23 from as correlated the same 2D prior
//#define __MARGINALIZE_TH13_DELTACP_CORR__ // Extract mixing13 and deltaCP from as correlated the same 2D prior

//#define __TOYS_SYNCHRONIZE_IH_TO_NH__ // Make toys using the same extracted oscillated parameters for NH and IH

#define __FILL_OSCPARAM_TOY_HISTO__


// NEW OSC PRIOR
#define __USE_PRIOR_SINSQ2TH13__    // PDG gives sinsq2th
#define __USE_PRIOR_SINSQ2TH12__    // PDG gives sinsq2th
                                    //#define __USE_PRIOR_SINSQ2TH23__  // We don't use PDG
                                    //



#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <TSystem.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TString.h>
#include <TLegend.h>
#include <TSpline.h>
#include <TRandom3.h>
#include <TMath.h>

#include "valor/Conventions/Opt.h"
#include "valor/Utilities/Rnd.h"
#include "valor/Utilities/UtlHist.h"

#include "AnalysisOptions.h"
#include "Build.h"
#include "T2KPhysicalParameterLimits.h"
#include "SuperkExpectation.h"
#include "SuperkObservation.h"
#include "SuperkMC.h"
#include "SuperkDetectorFSI.h"
#include "NuFlux.h"
#include "NuXSec.h"
#include "BANFF.h"
#include "Utilities.h"
#include "UtlFitParamListFactory.h"
#include "FitMarginalization.h"

#include <TSystemDirectory.h>

// NEW FHC+RHC JOINT
#if defined(__MULTIRING_FHC_RUNNING__) && defined(__FHC_RUNNING__) && defined(__RHC_RUNNING__) && defined(__MULTIRING_FHC_NUE__)
//CCPi fit
#define __USE_1RMU_IN_CHISQ__
#define __USE_RHC1RMU_IN_CHISQ__
#define __USE_1RE_IN_CHISQ__
#define __USE_RHC1RE_IN_CHISQ__
#define __USE_MULTIRE_IN_CHISQ__
#else //JointFit
#define __USE_1RMU_IN_CHISQ__
#define __USE_RHC1RMU_IN_CHISQ__
#define __USE_1RE_IN_CHISQ__
#define __USE_RHC1RE_IN_CHISQ__
#endif

// 3 FHC sample testing
//#undef __RHC_RUNNING__
//#undef __USE_RHC1RMU_IN_CHISQ__
//#undef __USE_RHC1RE_IN_CHISQ__

// No CC1pi or FHC 1Rmu
//#undef __MULTIRING_FHC_RUNNING__
//#undef __MULTIRING_FHC_NUE__
//#undef __USE_MULTIRE_IN_CHISQ__
//#undef __USE_1RMU_IN_CHISQ__
//#undef __USE_RHC1RE_IN_CHISQ__

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::fstream;
using std::ofstream;

using namespace valor;
using namespace valor::conv;

ClassImp(FitMarginalization)
//________________________________________________________________________________________
FitMarginalization::FitMarginalization()
{
  // NEW SPECTRA MARG
  //flikeVArray                        = 0;
  flikeToyVArray                     = 0;
  //
  fChisqMargVArray                   = 0;
  
  //
  fhchisq_marg_1D                    = 0;
  fhchisq_marg_2D                    = 0;
  fhchisq_marg_1D_Step               = 0;
  fhchisq_marg_2D_Step               = 0;
  fhdeltachisq_marg_1D               = 0;
  fhdeltachisq_marg_2D               = 0;
  fhmap_chisqrange_fileid_1D         = 0;
  fhmap_chisqrange_fileid_2D         = 0;
  fparm                              = 0;
  // NEW MARG MH
  fparm_nh                           = 0;
  fparm_ih                           = 0;
  //
  fSavedToyParamLists                = 0;
  
  // to be removed once the sk sample issue is fixed
  //for(int ikns=0;ikns<conv::kns;ikns++) fhToyMargCurr[ikns] = 0;
  
  fShouldSaveToyParamLists           = false;
  fMaxSampleErecoBins                = 1;
  
  // NEW MARG MH 2
  fCurrentFileIDminim = conv::kBadResult;
  ftagfile = "";
  fnamefile = 0;
  fGridStartMinBin_X = conv::kBadResult;
  fGridStartMaxBin_X = conv::kBadResult;
  fGridStartMinBin_Y = conv::kBadResult;
  fGridStartMaxBin_Y = conv::kBadResult;
  fGridStartMinVal_X = conv::kBadResult;
  fGridStartMaxVal_X = conv::kBadResult;
  fGridStartMinVal_Y = conv::kBadResult;
  fGridStartMaxVal_Y = conv::kBadResult;
  fnamefreeparam[0] = "";
  fnamefreeparam[1] = "";
  fPhysParamFree[0] = conv::kPhypUndefined;
  fPhysParamFree[1] = conv::kPhypUndefined;
  fIsFit2D = false;
  fdoMargMassHierarchy = false;
  fRethrowSystUnphysical = false;
  //fnamefilebinary = ;
  fchisq_bf = conv::kALargeNum;
  fparam_bf[0] = conv::kBadResult;
  fparam_bf[1] = conv::kBadResult;
  fparam_bf_errpos[0] = conv::kBadResult;
  fparam_bf_errpos[1] = conv::kBadResult;
  fparam_bf_errneg[0] = conv::kBadResult;
  fparam_bf_errneg[1] = conv::kBadResult;
  fparam_bf_errsym[0] = conv::kBadResult;
  fparam_bf_errsym[1] = conv::kBadResult;
  fSeed = 999;
  fnbins_param[0] = 0;
  fnbins_param[1] = 0;
  fnbins_sub_param[0] = 0;
  fnbins_sub_param[1] = 0;
  fntoys_marg = 0;
  fmin_bin_param[0] = 0;
  fmin_bin_param[1] = 0;
  fmin_bin_sub_param[0] = 0;
  fmin_bin_sub_param[1] = 0;
  fmin_sub_param[0] = 0;
  fmin_sub_param[1] = 0;
  fmin_param[0] = 0;
  fmin_param[1] = 0;
  fstep_param[0] = 0;
  fstep_param[1] = 0;
  fNCallsCalcChisqMargPt = 0;
  // NEW SPECTRA MARG
  //flikeVArray = 0;
  flikeToyVArray = 0;
  //
  fChisqMargVArray = 0;
  fNSKsamples = 0;
  fsksamples = 0;
  
  fnerecobins = 0;
  
  
  
  // NEW 2D TEMPLATES (VARRAY)
  
  //fToyMargVArray = 0;
  //fToyMargVArray_Clone = 0;
  fToyMargVArray_1Re = 0;
  fToyMargVArray_1Rmu = 0;
  fToyMargVArray_RHC1Rmu = 0;
  fToyMargVArray_RHC1Re = 0;
  fToyMargVArray_MultiRe = 0;
  
  fToyMargVArray_Clone_1Re = 0;
  fToyMargVArray_Clone_1Rmu = 0;
  fToyMargVArray_Clone_RHC1Rmu = 0;
  fToyMargVArray_Clone_RHC1Re = 0;
  fToyMargVArray_Clone_MultiRe = 0;
  //
  
  for(int i=0;i<conv::kns;i++) fhToyMargCurr[i] = 0;
  fMaxSampleErecoBins = 0;
  
  // NEW 2D TEMPLATES
  fnq2bins = 0;
  for(int i=0;i<conv::kns;i++) fhToyMargCurr2D[i] = 0;
  fMaxSampleQ2Bins = 0;
  //
  
  fMergeInputFiles_Called = true;
  fNmergedFiles = 0;
  
  fparm = 0;
  fparm_nh = 0;
  fparm_ih = 0;
  fProbNormalHierarchy = 0;
  //fhypothesis
  fskexp = 0;
  fhchisq_marg_1D_Step = 0;
  fhchisq_marg_2D_Step = 0;
  fMaxIterationSingleFit = 0;
  fMaxFitNumber = 0;
  fBorderDoFit = 0;
  fBorderMaxFitNumber = 0;
  fUncorrXY = false;
  fFirstMinimVal_X = 0;
  fLastMinimVal_X = 0;
  fFirstMinimVal_Y = 0;
  fLastMinimVal_Y = 0;
  fFirstMinimBin_X = 0;
  fLastMinimBin_X = 0;
  fFirstMinimBin_Y = 0;
  fLastMinimBin_Y = 0;
  fStep = 0;
  fStartBin_X = 0;
  fStartBin_Y = 0;
  fNumStartBin_X = 0;
  fNumStartBin_Y = 0;
  fExtractStartPt = "";
  fNFits = 0;
  fNiterations = 0;
  fFitStatus = 0;
  for(int i=0;i<conv::kns;i++) fPOT[i] = 0;
  for(int i=0;i<conv::kns;i++) fNormPOT[i] = 0;
  fteststatistic = "";
  fhPriorMH = 0;
  fhToyMH = 0;
  
  //NuEBarAppHisto
  for(int i = 0;i<2;i++){
    fNuEAppHisto[i] = 0;
    fNuEAppHistoUnscaled[i] = 0;
  }
  
  // NEW PENALTY TERM
  fpenalty_mean0 = 0.;
  fpenalty_mean1 = 0.;
  fpenalty_sigma0 = 0.;
  fpenalty_sigma1 = 0.;
  fhpenaltyterm_1D = 0;
  fhpenaltyterm_2D = 0;
  fusepenalty0 = false;
  fusepenalty1 = false;
  //
  
  //
}
//________________________________________________________________________________________
FitMarginalization::~FitMarginalization(void)
{
  this->CleanUp();
}

//________________________________________________________________________________________
void FitMarginalization::Init(SuperkExpectation *skexp){
  
  // Test statistic to be used in the fit
  if( fteststatistic!="rate+shape" && fteststatistic!="rate-only" && fteststatistic!=""){
    cerr << endl;
    cerr << "FitMarginalization::Init" << endl;
    cerr << "The test statistic is not available: " << fteststatistic << endl;
    cerr << endl;
    util::Die(1);
  }
  else if( fteststatistic=="" ){
    fteststatistic = "rate+shape";
  }
  cout << endl;
  cout << "The used test statistic is: " << fteststatistic << endl;
  cout << endl;
  
  
  // Before Init() we need SetPhysParFree()
  
  fdoMargMassHierarchy = false;
  fRethrowSystUnphysical = true; // rethrow toys if systematics are unphysical
  
  fIsFit2D = false; // 1D fit by default
  
  fSeed = 999;
  
  fchisq_bf = conv::kALargeNum;
  fparam_bf_errpos[0] = conv::kBadResult;
  fparam_bf_errpos[1] = conv::kBadResult;
  fparam_bf_errneg[0] = conv::kBadResult;
  fparam_bf_errneg[1] = conv::kBadResult;
  fparam_bf_errsym[0] = conv::kBadResult;
  fparam_bf_errsym[1] = conv::kBadResult;
  fparam_bf[0] = conv::kBadResult;
  fparam_bf[1] = conv::kBadResult;
  
  sprintf(fnamefilebinary,"%s","%s_file%d_toy%d_x_oscpar%d_oscpar%d_MH%s%s%s");
  
  // NEW SPECTRA MARG
  sprintf(fnamelikefilebinary,"%s","liketoy_%s_toy%d_x_oscpar%d_oscpar%d_MH%s%s");
  sprintf(fnameoscfiletree,"%s","%s_osctoys_file%d_toy%d_x_oscpar%d_oscpar%d_MH%s%s");
  //
  
  fMergeInputFiles_Called = true; // false the first time MergeInputFiles() is called
  fNmergedFiles = 0; // fNmergedFiles++ everytime MergeInputFiles() is called (after reading each file)
  
  fNCallsCalcChisqMargPt = 0; // Get track of the # of calls to CalcChisqMargPt
  
  const int nmh=2; // mass hierarchy
  
  // Get number of sk samples used in the fit
  // Needed to initialize the Ereco arrays
  
  this->SetNSKsamples(); // also defines fnerecobins
  
  
  
  
  
  // NEW 2D TEMPLATES (VARRAY)
  
  // #ifndef __USE_2D_TEMPLATES__
  //   fToyMargVArray = new VALORArrayF(1,1,1,1,1);
  //   fToyMargVArray->SetAll(-1.);   // needed for MergeInputFiles()
  //   fToyMargVArray->SetDimensionName(0,"# toys marginalization");
  //   fToyMargVArray->SetDimensionName(1,"# bins free parameter 0");
  //   fToyMargVArray->SetDimensionName(2,"# bins free parameter 1");
  //   fToyMargVArray->SetDimensionName(3,"# sk samples");
  //   fToyMargVArray->SetDimensionName(4,"# bins ereco");
  // #else
  //   // NEW 2D TEMPLATES
  //   fToyMargVArray = new VALORArrayF(1,1,1,1,1,1);
  //   fToyMargVArray->SetAll(-1.);   // needed for MergeInputFiles()
  //   fToyMargVArray->SetDimensionName(0,"# toys marginalization");
  //   fToyMargVArray->SetDimensionName(1,"# bins free parameter 0");
  //   fToyMargVArray->SetDimensionName(2,"# bins free parameter 1");
  //   fToyMargVArray->SetDimensionName(3,"# sk samples");
  //   fToyMargVArray->SetDimensionName(4,"# bins ereco");
  //   fToyMargVArray->SetDimensionName(5,"# bins q2"); // Follow same convention as SuperkMC. Q2 is the last element, opposite to histograms
  //   //
  // #endif
  
#ifndef __USE_2D_TEMPLATES__
  fToyMargVArray_1Re     = new VALORArrayF(1,1,1,1);
  fToyMargVArray_1Rmu    = new VALORArrayF(1,1,1,1);
  fToyMargVArray_RHC1Re  = new VALORArrayF(1,1,1,1);
  fToyMargVArray_RHC1Rmu = new VALORArrayF(1,1,1,1);
  fToyMargVArray_MultiRe = new VALORArrayF(1,1,1,1);
#else
  fToyMargVArray_1Re     = new VALORArrayF(1,1,1,1,1);
  fToyMargVArray_1Rmu    = new VALORArrayF(1,1,1,1,1);
  fToyMargVArray_RHC1Re  = new VALORArrayF(1,1,1,1,1);
  fToyMargVArray_RHC1Rmu = new VALORArrayF(1,1,1,1,1);
  fToyMargVArray_MultiRe = new VALORArrayF(1,1,1,1,1);
#endif
  
  fToyMargVArray_1Re->SetAll(-1.);   // needed for MergeInputFiles()
  fToyMargVArray_1Re->SetDimensionName(0,"# toys marginalization");
  fToyMargVArray_1Re->SetDimensionName(1,"# bins free parameter 0");
  fToyMargVArray_1Re->SetDimensionName(2,"# bins free parameter 1");
  fToyMargVArray_1Re->SetDimensionName(3,"# bins ereco");
  
  fToyMargVArray_1Rmu->SetAll(-1.);   // needed for MergeInputFiles()
  fToyMargVArray_1Rmu->SetDimensionName(0,"# toys marginalization");
  fToyMargVArray_1Rmu->SetDimensionName(1,"# bins free parameter 0");
  fToyMargVArray_1Rmu->SetDimensionName(2,"# bins free parameter 1");
  fToyMargVArray_1Rmu->SetDimensionName(3,"# bins ereco");
  
  fToyMargVArray_RHC1Re->SetAll(-1.);   // needed for MergeInputFiles()
  fToyMargVArray_RHC1Re->SetDimensionName(0,"# toys marginalization");
  fToyMargVArray_RHC1Re->SetDimensionName(1,"# bins free parameter 0");
  fToyMargVArray_RHC1Re->SetDimensionName(2,"# bins free parameter 1");
  fToyMargVArray_RHC1Re->SetDimensionName(3,"# bins ereco");
  
  fToyMargVArray_RHC1Rmu->SetAll(-1.);   // needed for MergeInputFiles()
  fToyMargVArray_RHC1Rmu->SetDimensionName(0,"# toys marginalization");
  fToyMargVArray_RHC1Rmu->SetDimensionName(1,"# bins free parameter 0");
  fToyMargVArray_RHC1Rmu->SetDimensionName(2,"# bins free parameter 1");
  fToyMargVArray_RHC1Rmu->SetDimensionName(3,"# bins ereco");
  
  fToyMargVArray_MultiRe->SetAll(-1.);   // needed for MergeInputFiles()
  fToyMargVArray_MultiRe->SetDimensionName(0,"# toys marginalization");
  fToyMargVArray_MultiRe->SetDimensionName(1,"# bins free parameter 0");
  fToyMargVArray_MultiRe->SetDimensionName(2,"# bins free parameter 1");
  fToyMargVArray_MultiRe->SetDimensionName(3,"# bins ereco");
  
#ifdef __USE_2D_TEMPLATES__
  // Follow same convention as SuperkMC. Q2 is the last element, opposite to histograms
  fToyMargVArray_1Re    ->SetDimensionName(4,"# bins q2");
  fToyMargVArray_1Rmu   ->SetDimensionName(4,"# bins q2");
  fToyMargVArray_RHC1Re ->SetDimensionName(4,"# bins q2");
  fToyMargVArray_RHC1Rmu->SetDimensionName(4,"# bins q2");
  fToyMargVArray_MultiRe->SetDimensionName(4,"# bins q2");
#endif
  
  //
  
  // need to declare array of histograms (then set to 0) because
  // problems with ROOT when declaring double pointers as global in header file
  // # of elements = conv::kns, but the used ones are is<fNSKsamples
  // Set to 0 if is>=fNSKsamples
  
#ifndef __USE_2D_TEMPLATES__
  for(int is=0;is<conv::kns;is++){
    if(is<fNSKsamples){
      char name[256];
      SuperkSample_t s = this->GetSKSamples(is);
      sprintf(name,"hToyMargCurr_%s",util::AsStringCompact(s).c_str());
      fhToyMargCurr[is] = new TH1F(name,name,util::NumErecoBins(s),util::ErecoBinEdges(s));
      cout << name << endl;
    }
    else{
      fhToyMargCurr[is] = 0;
    }
  }
#else
  // NEW 2D TEMPLATES
  for(int is=0;is<conv::kns;is++){
    if(is<fNSKsamples){
      char name[256];
      SuperkSample_t s = this->GetSKSamples(is);
      sprintf(name,"hToyMargCurr2D_%s",util::AsStringCompact(s).c_str());
      fhToyMargCurr2D[is] = new TH2F(name,name,util::NumQ2Bins(s),util::Q2BinEdges(s),util::NumErecoBins(s),util::ErecoBinEdges(s));
      cout << name << endl;
    }
    else{
      fhToyMargCurr2D[is] = 0;
    }
  }
  //
#endif
  
  //
  
  fChisqMargVArray            = new VALORArrayD(fnbins_param[0], fnbins_param[1]); // NEW PRECISION
  fChisqMargVArray->SetAll(0.);
  // NEW SPECTRA MARG
  //flikeVArray                 = new VALORArrayD(fnbins_param[0], fnbins_param[1]); // NEW PRECISION
  //flikeVArray->SetAll(0.);
  
  //
  
  char Xtitle[256];
  char Ytitle[256];
  sprintf(Xtitle,"%s",util::AsStringPaperRooTex(fPhysParamFree[0]).c_str());
  sprintf(Ytitle,"%s",util::AsStringPaperRooTex(fPhysParamFree[1]).c_str());
  
  // cout << "fmin_param[0] = " << fmin_param[0] << endl;
  // cout << "fmax_param[0] = " << fmax_param[0] << endl;
  // cout << "fstep_param[0] = " << fstep_param[0] << endl;
  // cout << "fmin_param[1] = " << fmin_param[1] << endl;
  // cout << "fmax_param[1] = " << fmax_param[1] << endl;
  // cout << "fstep_param[1] = " << fstep_param[1] << endl;
  
  fhchisq_marg_1D = new TH1D("hchisq_marg_1D","hchisq_marg_1D",fnbins_param[0],fmin_param[0]-fstep_param[0]/2.,fmax_param[0]+fstep_param[0]/2.);  // NEW PRECISION
  fhchisq_marg_1D->GetXaxis()->SetTitle(Xtitle);
  fhchisq_marg_1D->GetYaxis()->SetTitle("#Delta #chi^{2}");
  
  cout << endl;
  cout << "Initialize histogram fchisq_marg_1D:" << endl;
  cout << "BinCenter(1) = " << fhchisq_marg_1D->GetXaxis()->GetBinCenter(1) << endl;
  cout << "BinCenter(Nbins) = " << fhchisq_marg_1D->GetXaxis()->GetBinCenter(fnbins_param[0]) << endl;
  
  // NEW PRECISION
  fhchisq_marg_2D = new TH2D("hchisq_marg_2D","hchisq_marg_2D",
                             fnbins_param[0],fmin_param[0]-fstep_param[0]/2.,fmax_param[0]+fstep_param[0]/2.,
                             fnbins_param[1],fmin_param[1]-fstep_param[1]/2.,fmax_param[1]+fstep_param[1]/2. );
  //
  fhchisq_marg_2D->GetXaxis()->SetTitle(Xtitle);
  fhchisq_marg_2D->GetYaxis()->SetTitle(Ytitle);
  
  cout << endl;
  cout << "Initialize histogram fchisq_marg_2D:" << endl;
  cout << "X: BinCenter(1) = " << fhchisq_marg_2D->GetXaxis()->GetBinCenter(1) << endl;
  cout << "X: BinCenter(Nbins) = " << fhchisq_marg_2D->GetXaxis()->GetBinCenter(fnbins_param[1]) << endl;
  cout << "Y: BinCenter(1) = " << fhchisq_marg_2D->GetYaxis()->GetBinCenter(1) << endl;
  cout << "Y: BinCenter(Nbins) = " << fhchisq_marg_2D->GetYaxis()->GetBinCenter(fnbins_param[1]) << endl;
  
  
  //fhdeltachisq_marg_1D = 0;
  //fhdeltachisq_marg_2D = 0;
  fhdeltachisq_marg_1D = (TH1D*) fhchisq_marg_1D->Clone("fhdeltachisq_marg_1D");
  fhdeltachisq_marg_1D->Reset();
  fhdeltachisq_marg_1D->SetName("fhdeltachisq_marg_1D");
  fhdeltachisq_marg_1D->SetTitle("fhdeltachisq_marg_1D");
  fhdeltachisq_marg_2D = (TH2D*) fhchisq_marg_2D->Clone("fhdeltachisq_marg_2D");
  fhdeltachisq_marg_2D->Reset();
  fhdeltachisq_marg_2D->SetName("fhdeltachisq_marg_2D");
  fhdeltachisq_marg_2D->SetTitle("fhdeltachisq_marg_2D");
  
  fhchisq_marg_1D_Step = (TH1D*)fhchisq_marg_1D->Clone("hchisq_marg_1D_Step");
  fhchisq_marg_1D_Step->Reset();
  fhchisq_marg_1D_Step->SetName("hchisq_marg_1D_Step");
  fhchisq_marg_1D_Step->SetTitle("hchisq_marg_1D_Step");
  fhchisq_marg_2D_Step = (TH2D*)fhchisq_marg_2D->Clone("hchisq_marg_2D_Step");
  fhchisq_marg_2D_Step->Reset();
  fhchisq_marg_2D_Step->SetName("hchisq_marg_2D_Step");
  fhchisq_marg_2D_Step->SetTitle("hchisq_marg_2D_Step");
  
  // NEW PENALTY TERM
  fhpenaltyterm_1D = (TH1D*) fhchisq_marg_1D->Clone("hpenaltyterm_1D");
  fhpenaltyterm_1D->Reset();
  fhpenaltyterm_1D->SetName("hpenaltyterm_1D");
  fhpenaltyterm_1D->SetTitle("hpenaltyterm_1D");
  fhpenaltyterm_2D = (TH2D*) fhchisq_marg_2D->Clone("hpenaltyterm_2D");
  fhpenaltyterm_2D->Reset();
  fhpenaltyterm_2D->SetName("hpenaltyterm_2D");
  fhpenaltyterm_2D->SetTitle("hpenaltyterm_2D");
  
  //
  
  
  cout << endl;
  cout << "Initialize VALORArray fToyMargVArray with:"<< endl;
  cout << "fntoys_marg = " << fntoys_marg << endl;
  cout << endl;
  cout << "fNSKsamples = " << fNSKsamples << endl;
  cout << "Used sub-samples:" << endl;
  
  //for(int is = 0; is < kns; is++) {
  //SuperkSample_t s = kSKSamples[is];
  //if(s == kSKsUndefined) continue;
  
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s = this->GetSKSamples(is);
    cout << "   - " << util::AsStringCompact(s)
    << " - # bins ereco = " << this->fnerecobins->GetEntry(is)
#ifdef __USE_2D_TEMPLATES__
    << " - # bins q2 = " << this->fnq2bins->GetEntry(is) // NEW 2D TEMPLATES
#endif
    << endl;
  }
  cout << "Max # of Ereco bins = " << fMaxSampleErecoBins << endl;
#ifdef __USE_2D_TEMPLATES__
  cout << "Max # of Q2 bins = " << fMaxSampleQ2Bins << endl; // NEW 2D TEMPLATES
#endif
  cout << endl;
  cout << "fnbins_sub_param[0] = " << fnbins_sub_param[0] << endl;
  cout << "fnbins_sub_param[1] = " << fnbins_sub_param[1] << endl;
  cout << "nmh = " << nmh << endl;
  cout << endl;
  
  
  // NEW FHC+RHC JOINT
  cout << "Build chisq using the following SK samples:" << endl;
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s = this->GetSKSamples(is);
    bool usesample = this->UseSKSample(is);
    if(!usesample) continue;
    cout << util::AsString(s) << endl;
  }
  cout << endl;
  //
  
  
  fmin_sub_param[0] = fhchisq_marg_2D->GetXaxis()->GetBinCenter(fmin_bin_sub_param[0]+1);
  fmax_sub_param[0] = fhchisq_marg_2D->GetXaxis()->GetBinCenter(fmax_bin_sub_param[0]+1);
  fmin_sub_param[1] = fhchisq_marg_2D->GetYaxis()->GetBinCenter(fmin_bin_sub_param[1]+1);
  fmax_sub_param[1] = fhchisq_marg_2D->GetYaxis()->GetBinCenter(fmax_bin_sub_param[1]+1);
  
  cout << endl;
  cout << "Free parameter 0 sub-range: ";
  cout << "[" << fmin_sub_param[0] << "," << fmax_sub_param[0] << "]" << endl;
  cout << "# of bins parameter 0 = " << fnbins_sub_param[0] << endl;
  cout << endl;
  cout << "Free parameter 1 sub-range: ";
  cout << "[" << fmin_sub_param[1] << "," << fmax_sub_param[1] << "]" << endl;
  cout << "# of bins parameter 1 = " << fnbins_sub_param[1] << endl;
  cout << endl;
  
  cout << "fmin_bin_sub_param[0] = " << fmin_bin_sub_param[0] << endl;
  cout << "fmax_bin_sub_param[0] = " << fmax_bin_sub_param[0] << endl;
  cout << "fmin_bin_sub_param[1] = " << fmin_bin_sub_param[1] << endl;
  cout << "fmax_bin_sub_param[1] = " << fmax_bin_sub_param[1] << endl;
  
  
  
  // Variables for minimization of the marginalized Chisq
  fUncorrXY = true; // by default if minimum found in Y keep it fixed and vary only X (viceversa)
  fMaxIterationSingleFit = 1000;
  fMaxFitNumber = 10;
  
  fBorderDoFit = true;
  if(fIsFit2D) fBorderMaxFitNumber = 4; // 2dim x 2 borders
  else         fBorderMaxFitNumber = 2; // 1dim x 2 borders
  
  fStartBin_X = 1;
  fStartBin_Y = 1;
  fNumStartBin_X = 1;
  fNumStartBin_Y = 1;
  
  fExtractStartPt = "";
  fNFits = 0;
  fNiterations = 0;
  fFitStatus = 0;
  
  // Set the interval where to find the minimum (bin numbers), when minimizing the chisq marg
  fFirstMinimBin_X = 1;
  fLastMinimBin_X = fhchisq_marg_2D->GetNbinsX();
  fFirstMinimBin_Y = 1;
  fLastMinimBin_Y = fhchisq_marg_2D->GetNbinsY();
  
  // Set the default map
  fhmap_chisqrange_fileid_1D = (TH1D*) fhchisq_marg_1D->Clone("hmap_chisqrange_fileid_1D");
  fhmap_chisqrange_fileid_1D->Reset();
  fhmap_chisqrange_fileid_1D->SetName("hmap_chisqrange_fileid_1D");
  fhmap_chisqrange_fileid_1D->SetTitle("hmap_chisqrange_fileid_1D");
  fhmap_chisqrange_fileid_2D = (TH2D*) fhchisq_marg_2D->Clone("hmap_chisqrange_fileid_2D");
  fhmap_chisqrange_fileid_2D->Reset();
  fhmap_chisqrange_fileid_2D->SetName("hmap_chisqrange_fileid_2D");
  fhmap_chisqrange_fileid_2D->SetTitle("hmap_chisqrange_fileid_2D");
  
  fCurrentFileIDminim = conv::kBadResult; // current file ID during the chisq minimization process
  
  
  // Store # of POT for each used sample
  
  //for(int is = 0; is < kns; is++) {
  //SuperkSample_t s = kSKSamples[is];
  //if(s == kSKsUndefined) continue;
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s = this->GetSKSamples(is);
    fPOT[is] = skexp->GetSuperkMC()->AnalyzedSamplePOT(s);
    fNormPOT[is] = fPOT[is] / (1E+20);
    cout << "POT ( " << util::AsString(s) << " ) = " << fPOT[is] << endl;
    cout << "Norm POT ( " << util::AsString(s) << " ) = " << fNormPOT[is] << endl;
  }
  
  fProbNormalHierarchy = 0.5; // As default initialized to flat prior
  fhPriorMH = new TH1D("hPriorMH","Prior MH: bin 1 --> P(NH), bin 2 --> P(IH)",2,0,2);
  fhPriorMH->SetBinContent(1,fProbNormalHierarchy);
  fhPriorMH->SetBinContent(2,1.-fProbNormalHierarchy);
  
  fhToyMH = new TH1D("hToyMH","Toy MH: bin 1 --> NH, bin 2 --> IH",2,0,2);
  
  fStoreLikeWeights = false; // don't store the full likelihood by default
  
  //NuEApp Raj
  for(int i = 0;i<2;i++){
    if(fPhysParamFree[0] == kPhypBetaBar && fPhysParamFree[1] == kPhypUndefined){
      fNuEAppHisto[i] = new TH1D(Form("NuEBarAppHisto_%d",i),Form("Expected NuEBar Events; BetaBar bin %d",i),25000,0,25);
      fNuEAppHistoUnscaled[i] = new TH1D(Form("NuEBarAppHistoUnscaled_%d",i),Form("Expected NuEBar Events Unscaled; BetaBar bin %d",i),25000,0,25);
    }
  }
  
}
//________________________________________________________________________________________
void FitMarginalization::SetStoreLikeWeights(bool doStoreLikeWeights){
  fStoreLikeWeights = doStoreLikeWeights;
  if(this->fStoreLikeWeights){
    flikeToyVArray              = new VALORArrayD(fnbins_param[0], fnbins_param[1], fntoys_marg);
    flikeToyVArray->SetAll(0.);
  }
}
//________________________________________________________________________________________
void FitMarginalization::CleanUp(void){
  if (fSavedToyParamLists) {
    for (int i = 0 ; i < fntoys_marg ; i++) {
      if (fSavedToyParamLists[i]) { delete fSavedToyParamLists[i];}
    }
    delete fSavedToyParamLists;
    fSavedToyParamLists = 0;
  }
  
  // NEW SPECTRA MARG
  //if (flikeVArray                       ) { delete flikeVArray                       ; flikeVArray                   = 0;}
  if (flikeToyVArray                    ) { delete flikeToyVArray                    ; flikeToyVArray                = 0;}
  //
  if (fChisqMargVArray                  ) { delete fChisqMargVArray                  ; fChisqMargVArray              = 0;}
  if (fnerecobins                       ) { delete fnerecobins                       ; fnerecobins                   = 0;}
  // NEW 2D TEMPLATES
  if (fnq2bins                          ) { delete fnq2bins                          ; fnq2bins                      = 0;}
  //
  if (fsksamples                        ) { delete fsksamples                        ; fsksamples                    = 0;}
  
  
  // NEW 2D TEMPLATES (VARRAY)
  //if (fToyMargVArray                    ) { delete fToyMargVArray                    ; fToyMargVArray                = 0;}
  //if (fToyMargVArray_Clone              ) { delete fToyMargVArray_Clone              ; fToyMargVArray_Clone          = 0;}
  if (fToyMargVArray_1Re                ) { delete fToyMargVArray_1Re                ; fToyMargVArray_1Re            = 0;}
  if (fToyMargVArray_1Rmu               ) { delete fToyMargVArray_1Rmu               ; fToyMargVArray_1Rmu           = 0;}
  if (fToyMargVArray_RHC1Re             ) { delete fToyMargVArray_RHC1Re             ; fToyMargVArray_RHC1Re         = 0;}
  if (fToyMargVArray_RHC1Rmu            ) { delete fToyMargVArray_RHC1Rmu            ; fToyMargVArray_RHC1Rmu        = 0;}
  if (fToyMargVArray_Clone_1Re          ) { delete fToyMargVArray_Clone_1Re          ; fToyMargVArray_Clone_1Re      = 0;}
  if (fToyMargVArray_Clone_1Rmu         ) { delete fToyMargVArray_Clone_1Rmu         ; fToyMargVArray_Clone_1Rmu     = 0;}
  if (fToyMargVArray_Clone_RHC1Re       ) { delete fToyMargVArray_Clone_RHC1Re       ; fToyMargVArray_Clone_RHC1Re   = 0;}
  if (fToyMargVArray_Clone_RHC1Rmu      ) { delete fToyMargVArray_Clone_RHC1Rmu      ; fToyMargVArray_Clone_RHC1Rmu  = 0;}
  //
  
  // NEW PENALTY TERM
  if (fhpenaltyterm_1D                  ) { delete fhpenaltyterm_1D                  ; fhpenaltyterm_1D                   = 0;}
  if (fhpenaltyterm_2D                  ) { delete fhpenaltyterm_2D                  ; fhpenaltyterm_2D                   = 0;}
  //
  
  if (fhchisq_marg_1D                   ) { delete fhchisq_marg_1D                   ; fhchisq_marg_1D                    = 0;}
  if (fhchisq_marg_2D                   ) { delete fhchisq_marg_2D                   ; fhchisq_marg_2D                    = 0;}
  if (fhchisq_marg_1D_Step              ) { delete fhchisq_marg_1D_Step              ; fhchisq_marg_1D_Step               = 0;}
  if (fhchisq_marg_2D_Step              ) { delete fhchisq_marg_2D_Step              ; fhchisq_marg_2D_Step               = 0;}
  if (fhdeltachisq_marg_1D              ) { delete fhdeltachisq_marg_1D              ; fhdeltachisq_marg_1D               = 0;}
  if (fhdeltachisq_marg_2D              ) { delete fhdeltachisq_marg_2D              ; fhdeltachisq_marg_2D               = 0;}
  if (fhmap_chisqrange_fileid_1D        ) { delete fhmap_chisqrange_fileid_1D        ; fhmap_chisqrange_fileid_1D         = 0;}
  if (fhmap_chisqrange_fileid_2D        ) { delete fhmap_chisqrange_fileid_2D        ; fhmap_chisqrange_fileid_2D         = 0;}
  if (fparm                             ) { delete fparm                             ; fparm                              = 0;}
  if (fparm_nh                          ) { delete fparm_nh                          ; fparm_nh                           = 0;}
  if (fparm_ih                          ) { delete fparm_ih                          ; fparm_ih                           = 0;}
  if (fSavedToyParamLists               ) { delete fSavedToyParamLists               ; fSavedToyParamLists                = 0;}
  if (fhPriorMH                         ) { delete fhPriorMH                         ; fhPriorMH                          = 0;}
  if (fhToyMH                           ) { delete fhToyMH                           ; fhToyMH                            = 0;}
  
  // to be removed once the sk sample issue is fixed
#ifndef __USE_2D_TEMPLATES__
  for(int ikns=0;ikns<conv::kns;ikns++){
    if (fhToyMargCurr[ikns]              ) { delete fhToyMargCurr[ikns]               ; fhToyMargCurr[ikns]                = 0;}
  }
  fMaxSampleErecoBins = 0;
#else
  // NEW 2D TEMPLATES
  for(int ikns=0;ikns<conv::kns;ikns++){
    if (fhToyMargCurr2D[ikns]            ) { delete fhToyMargCurr2D[ikns]           ; fhToyMargCurr2D[ikns]                = 0;}
  }
  fMaxSampleQ2Bins = 0;
  //
#endif
  
#ifdef FILL_TREE_OSCPARAM_TOYS
  if (treetoys                         ) { delete treetoys                         ; treetoys                = 0;}  // NEW OSC TREE
  if (filetreetoys                     ) { delete filetreetoys                     ; filetreetoys            = 0;}  // NEW OSC TREE
#endif
  
  //NuEBarApp
  for(int i=0;i<2;i++){
    if(fNuEAppHisto[i]){
      delete fNuEAppHisto[i];
      fNuEAppHisto[i] = 0;
    }
    if(fNuEAppHistoUnscaled[i]){
      delete fNuEAppHistoUnscaled[i];
      fNuEAppHistoUnscaled[i] = 0;
    }
  }
}
//________________________________________________________________________________________
void FitMarginalization::SetFullRange(){
  
  // Default values
  
  fmin_bin_param[0] = 0;
  fmax_bin_param[0] = fnbins_param[0] - 1;
  fmin_bin_param[1] = 0;
  fmax_bin_param[1] = fnbins_param[1] - 1;
  
  if(fPhysParamFree[0] == conv::kPhypUndefined){
    fmin_bin_param[0] = 0.;
    fmax_bin_param[0] = 0.;
    fnbins_param[0] = 1; // dimension of parameter ParId
    cerr << endl;
    cerr << "FitMarginalization::SetFullRange" << endl;
    cerr << "Param 0 is undefined" << endl;
    cerr << endl;
    util::Die(1);
  }
  if(fPhysParamFree[1] == conv::kPhypUndefined){
    fmin_bin_param[1] = 0.;
    fmax_bin_param[1] = 0.;
    fnbins_param[1] = 1; // dimension of parameter ParId
  }
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetSubRange(int *submin,int *submax){
  
  // sub-range: scan all the range by default
  fmin_bin_sub_param[0] = submin[0];
  fmax_bin_sub_param[0] = submax[0];
  fmin_bin_sub_param[1] = submin[1];
  fmax_bin_sub_param[1] = submax[1];
  
  // default values
  if(fmin_bin_sub_param[0] == conv::kBadResult) fmin_bin_sub_param[0] = 0;
  if(fmax_bin_sub_param[0] == conv::kBadResult) fmax_bin_sub_param[0] = fnbins_param[0] - 1;
  if(fmin_bin_sub_param[1] == conv::kBadResult) fmin_bin_sub_param[1] = 0;
  if(fmax_bin_sub_param[1] == conv::kBadResult) fmax_bin_sub_param[1] = fnbins_param[1] - 1;
  
  // # of bins
  fnbins_sub_param[0] = fmax_bin_sub_param[0] - fmin_bin_sub_param[0] + 1;
  fnbins_sub_param[1] = fmax_bin_sub_param[1] - fmin_bin_sub_param[1] + 1;
  
  if(fPhysParamFree[0] == conv::kPhypUndefined){
    fmin_bin_sub_param[0] = 0.;
    fmax_bin_sub_param[0] = 0.;
    fnbins_sub_param[0] = 1; // dimension of parameter ParId
  }
  if(fPhysParamFree[1] == conv::kPhypUndefined){
    fmin_bin_sub_param[1] = 0.;
    fmax_bin_sub_param[1] = 0.;
    fnbins_sub_param[1] = 1; // dimension of parameter ParId
  }
  
  // checks
  if( (fPhysParamFree[0] != conv::kPhypUndefined) &&
     ((fnbins_sub_param[0] > fnbins_param[0]) ||
      (fmin_bin_sub_param[0] < 0) ||
      (fmax_bin_sub_param[0] > fnbins_param[0] - 1) ))
  {
    cerr << endl;
    cerr << "FitMarginalization::SetSubRange" << endl;
    cerr << "Param 0 sub-range outside the entire free parameter range" << endl;
    cerr << "fnbins_param[0] = " << fnbins_param[0] << endl;
    cerr << "fnbins_sub_param[0] = " << fnbins_sub_param[0] << endl;
    cerr << "fmin_bin_sub_param[0] = " << fmin_bin_sub_param[0] << endl;
    cerr << "fmax_bin_sub_param[0] = " << fmax_bin_sub_param[0] << endl;
    cerr << endl;
    util::Die(1);
  }
  
  if( (fPhysParamFree[1] != conv::kPhypUndefined) &&
     ((fnbins_sub_param[1] > fnbins_param[1]) ||
      (fmin_bin_sub_param[1] < 0) ||
      (fmax_bin_sub_param[1] > fnbins_param[1] - 1) ))
  {
    cerr << endl;
    cerr << "FitMarginalization::SetSubRange" << endl;
    cerr << "Param 1 sub-range outside the entire free parameter range" << endl;
    cerr << "fnbins_param[1] = " << fnbins_param[1] << endl;
    cerr << "fnbins_sub_param[1] = " << fnbins_sub_param[1] << endl;
    cerr << "fmin_bin_sub_param[1] = " << fmin_bin_sub_param[1] << endl;
    cerr << "fmax_bin_sub_param[1] = " << fmax_bin_sub_param[1] << endl;
    cerr << endl;
    util::Die(1);
  }
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetPhysParFree(PhysParamId_t *parfree,int ntoys_marg,int *nbins_param,double *minpar,double *maxpar){
  
  // Set parameters
  fntoys_marg = ntoys_marg;
  fnbins_param[0] = nbins_param[0];
  fnbins_param[1] = nbins_param[1];
  
  // Set free parameters
  this->SetFreeParam(0,parfree,minpar,maxpar);
  this->SetFreeParam(1,parfree,minpar,maxpar);
  
  this->SetFullRange();
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetFreeParam(int ParId,conv::PhysParamId_t *parfree,double *minpar,double *maxpar){
  
  fPhysParamFree[ParId] = parfree[ParId];
  fnamefreeparam[ParId] = util::AsString(fPhysParamFree[ParId]);
  
  fmin_param[ParId] = minpar[ParId];
  fmax_param[ParId] = maxpar[ParId];
  
  // NuBar Osc parameters
  if(fPhysParamFree[ParId]==kPhypMixing13Bar){ // mixing13Bar
    if(minpar[ParId]==conv::kBadResult && maxpar[ParId]==conv::kBadResult){ // if limits are not given use the default ones
      fmin_param[ParId] = mixing13Bar_min;
      fmax_param[ParId] = mixing13Bar_max;
    }
  }
  else if(fPhysParamFree[ParId]==kPhypMixing23Bar){ // mixing23Bar
    if(minpar[ParId]==conv::kBadResult && maxpar[ParId]==conv::kBadResult){ // if limits are not given use the default ones
      fmin_param[ParId] = mixing23Bar_min;
      fmax_param[ParId] = mixing23Bar_max;
    }
  }
  else if(fPhysParamFree[ParId]==kPhypDm23sqBar){ // Dm23sqBar
    if(minpar[ParId]==conv::kBadResult && maxpar[ParId]==conv::kBadResult){ // if limits are not given use the default ones
      fmin_param[ParId] = dmsq32Bar_min;
      fmax_param[ParId] = dmsq32Bar_max;
    }
  }
  else if(fPhysParamFree[ParId]==kPhypCPVBar){ // deltaCPBar
    if(minpar[ParId]==conv::kBadResult && maxpar[ParId]==conv::kBadResult){ // if limits are not given use the default ones
      fmin_param[ParId] = deltacpBar_min;
      fmax_param[ParId] = deltacpBar_max;
    }
  }
  else if(fPhysParamFree[ParId]==kPhypMixing12Bar){ // mixing12Bar
    if(minpar[ParId]==conv::kBadResult && maxpar[ParId]==conv::kBadResult){ // if limits are not given use the default ones
      fmin_param[ParId] = mixing12Bar_min;
      fmax_param[ParId] = mixing12Bar_max;
    }
  }
  else if(fPhysParamFree[ParId]==kPhypDm12sqBar){ // Dm12sqBar
    if(minpar[ParId]==conv::kBadResult && maxpar[ParId]==conv::kBadResult){ // if limits are not given use the default ones
      fmin_param[ParId] = dmsq12Bar_min;
      fmax_param[ParId] = dmsq12Bar_max;
    }
  }
  
  // Nu Osc parameters
  else if(fPhysParamFree[ParId]==kPhypMixing13){ // mixing13
    if(minpar[ParId]==conv::kBadResult && maxpar[ParId]==conv::kBadResult){ // if limits are not given use the default ones
      fmin_param[ParId] = mixing13_min;
      fmax_param[ParId] = mixing13_max;
    }
  }
  else if(fPhysParamFree[ParId]==kPhypMixing23){ // mixing23
    if(minpar[ParId]==conv::kBadResult && maxpar[ParId]==conv::kBadResult){ // if limits are not given use the default ones
      fmin_param[ParId] = mixing23_min;
      fmax_param[ParId] = mixing23_max;
    }
  }
  else if(fPhysParamFree[ParId]==kPhypDm23sq){ // Dm23sq
    if(minpar[ParId]==conv::kBadResult && maxpar[ParId]==conv::kBadResult){ // if limits are not given use the default ones
      fmin_param[ParId] = dmsq32_min;
      fmax_param[ParId] = dmsq32_max;
    }
  }
  else if(fPhysParamFree[ParId]==kPhypCPV){ // deltaCP
    if(minpar[ParId]==conv::kBadResult && maxpar[ParId]==conv::kBadResult){ // if limits are not given use the default ones
      fmin_param[ParId] = deltacp_min;
      fmax_param[ParId] = deltacp_max;
    }
  }
  else if(fPhysParamFree[ParId]==kPhypMixing12){ // mixing12
    if(minpar[ParId]==conv::kBadResult && maxpar[ParId]==conv::kBadResult){ // if limits are not given use the default ones
      fmin_param[ParId] = mixing12_min;
      fmax_param[ParId] = mixing12_max;
    }
  }
  else if(fPhysParamFree[ParId]==kPhypDm12sq){ // Dm12sq
    if(minpar[ParId]==conv::kBadResult && maxpar[ParId]==conv::kBadResult){ // if limits are not given use the default ones
      fmin_param[ParId] = dmsq12_min;
      fmax_param[ParId] = dmsq12_max;
    }
  }
  
  // NEW BETABAR
  else if(fPhysParamFree[ParId]==kPhypBetaBar){ // BetaBar (numubar->nuebar normalization)
    if(minpar[ParId]==conv::kBadResult && maxpar[ParId]==conv::kBadResult){ // if limits are not given use the default ones
      fmin_param[ParId] = betaBar_min;
      fmax_param[ParId] = betaBar_max;
    }
  }
  //
  
  // The parameter is not selected
  else{
    if(ParId==0){
      cerr << "FitMarginalization::SetFreeParam" << endl;
      cerr << endl;
      cerr << "The parameter 0 has not been chosen. Always do 1D fit with parameter 0" << endl;
      cerr << endl;
      util::Die(1);
    }
    else{
      fmin_param[ParId] = 0.;
      fmax_param[ParId] = 0.;
      fnbins_param[ParId] = 1; // dimension of parameter ParId
      fPhysParamFree[ParId] = conv::kPhypUndefined;
    }
  }
  
  // Wrong assignement of the free parameter limit
  if(minpar[ParId] > maxpar[ParId]){
    cerr << endl;
    cerr << "FitMarginalization::SetFreeParam" << endl;
    cerr << "ParID = " << ParId << endl
    << "minpar >= maxpar && minpar!=-999 && maxpar!=-999" << endl
    << " - minpar = " << minpar[ParId] << endl
    << " - maxpar = " << maxpar[ParId] << endl
    << endl;
    util::Die(1);
  }
  
  // Define the grid step size
  if(fPhysParamFree[ParId] == conv::kPhypUndefined) fstep_param[ParId] = 0.;
  else fstep_param[ParId] = (fmax_param[ParId]-fmin_param[ParId])/(double)(fnbins_param[ParId] - 1);
  
  cout << endl;
  cout << "Minimize " << fnamefreeparam[ParId] << endl;
  cout << "# toys for marginalization: " << fntoys_marg << endl;
  cout << "# bins of free parameter to scan: " << fnbins_param[ParId] << endl;
  cout << "# free parameter step size: " << fstep_param[ParId] << endl;
  cout << "Interval: [" << fmin_param[ParId] << "," << fmax_param[ParId] << "]" << endl;
  cout << endl;
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetFitParamList(const char * fit_param_list){
  
  fparm  = util::NamedParameterList(fit_param_list);
  fparm_nh = new FitParamList(*fparm);
  fparm_ih = new FitParamList(*fparm);
  
  fhypothesis = fparm->PhysicsHypothesis();
  
  if( fPhysParamFree[0] == conv::kPhypUndefined &&
     fPhysParamFree[1] == conv::kPhypUndefined ){
    cerr << endl;
    cerr << "FitMarginalization::SetFitParamList" << endl;
    cerr << "Both free parameters are not defined" << endl;
    cerr << endl;
    util::Die(1);
  }
  
  if(fPhysParamFree[0]!=conv::kPhypUndefined && fPhysParamFree[1]!=conv::kPhypUndefined){
    fIsFit2D = true;
    cout << endl;
    cout << "Perform a fit 2D" << endl;
    cout << endl;
  }
  else if(fPhysParamFree[0]!=conv::kPhypUndefined){
    fIsFit2D = false;
    cout << endl;
    cout << "Perform a fit 1D" << endl;
    cout << endl;
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::SetFitParamList" << endl;
    cerr << "Define correctly parameter 0 for 1D fit and parameter 1 for 2D fit" << endl;
    cerr << endl;
    util::Die(1);
  }
  
  
  // Parameter 0
  if( fparm->IsFixed(fPhysParamFree[0]) ){
    cerr << "FitMarginalization::SetFitParamList" << endl;
    cerr << "The free parameter " << fnamefreeparam[0]
    << " is fixed in FitParamList " << fit_param_list
    << endl << *fparm << endl;
    
    util::Die(1);
  }
  else if( fparm->IsNuisance(fPhysParamFree[0]) ){
    cerr << "FitMarginalization::SetFitParamList" << endl;
    cerr << "The free parameter " << fnamefreeparam[0]
    << " is nuisance in FitParamList " << fit_param_list
    << endl << *fparm << endl;
    util::Die(1);
  }
  
  // Parameter 1
  if( fPhysParamFree[1] == conv::kPhypUndefined ){
    cout << endl;
    cout << "Parameter 1 is undefined" << endl;
    cout << endl;
  }
  else if( fparm->IsFixed(fPhysParamFree[1]) ){
    cerr << "FitMarginalization::SetFitParamList" << endl;
    cerr << "The free parameter " << fnamefreeparam[1]
    << " is fixed in FitParamList " << fit_param_list
    << endl << *fparm << endl;
    util::Die(1);
  }
  else if( fparm->IsNuisance(fPhysParamFree[1]) ){
    cerr << "FitMarginalization::SetFitParamList" << endl;
    cerr << "The free parameter " << fnamefreeparam[1]
    << " is nuisance in FitParamList " << fit_param_list
    << endl << *fparm << endl;
    util::Die(1);
  }
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::PrintPhysHypothesis(){
  cout << endl;
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(fhypothesis) << endl;
  cout << endl;
}
//________________________________________________________________________________________
void FitMarginalization::MakeToysMarginalization(int tweak_systematics,SuperkExpectation *skexp){
 
   TFile* toys_file = new TFile("toys.root", "recreate");
  
  
  // NEW 2D TEMPLATES (VARRAY)
  
  // #ifndef __USE_2D_TEMPLATES__
  //   fToyMargVArray->Resize(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], fNSKsamples, fMaxSampleErecoBins);
  // #else
  //   // NEW 2D TEMPLATES
  //   fToyMargVArray->Resize(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], fNSKsamples, fMaxSampleErecoBins, fMaxSampleQ2Bins);
  // #endif
  
  int nbins_ereco_1Re     = util::NumErecoBins(kSKsFCFV1RingELike);
  int nbins_ereco_1Rmu    = util::NumErecoBins(kSKsFCFV1RingMuLike);
  int nbins_ereco_RHC1Re  = util::NumErecoBins(kSKsRHCFCFV1RingELike);
  int nbins_ereco_RHC1Rmu = util::NumErecoBins(kSKsRHCFCFV1RingMuLike);
  int nbins_ereco_MultiRe = util::NumErecoBins(kSKsFCFVMultiRingELike);
  
#ifndef __USE_2D_TEMPLATES__
  fToyMargVArray_1Re    ->Resize(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_1Re);
  fToyMargVArray_1Rmu   ->Resize(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_1Rmu);
  fToyMargVArray_RHC1Re ->Resize(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_RHC1Re);
  fToyMargVArray_RHC1Rmu->Resize(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_RHC1Rmu);
  fToyMargVArray_MultiRe->Resize(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_MultiRe);
#else
  
  int nbins_q2_1Re     = util::NumQ2Bins(kSKsFCFV1RingELike);
  int nbins_q2_1Rmu    = util::NumQ2Bins(kSKsFCFV1RingMuLike);
  int nbins_q2_RHC1Re  = util::NumQ2Bins(kSKsRHCFCFV1RingELike);
  int nbins_q2_RHC1Rmu = util::NumQ2Bins(kSKsRHCFCFV1RingMuLike);
  int nbins_q2_MultiRe = util::NumQ2Bins(kSKsFCFVMultiRingELike);
  
  fToyMargVArray_1Re    ->Resize(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_1Re,     nbins_q2_1Re);
  fToyMargVArray_1Rmu   ->Resize(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_1Rmu,    nbins_q2_1Rmu);
  fToyMargVArray_RHC1Re ->Resize(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_RHC1Re,  nbins_q2_RHC1Re);
  fToyMargVArray_RHC1Rmu->Resize(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_RHC1Rmu, nbins_q2_RHC1Rmu);
  fToyMargVArray_MultiRe->Resize(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_MultiRe, nbins_q2_MultiRe);
#endif
  
  // const int *dimensions_1Re     = fToyMargVArray_1Re    ->GetDimensions();
  // const int *dimensions_1Rmu    = fToyMargVArray_1Rmu   ->GetDimensions();
  // const int *dimensions_RHC1Re  = fToyMargVArray_RHC1Re ->GetDimensions();
  // const int *dimensions_RHC1Rmu = fToyMargVArray_RHC1Rmu->GetDimensions();
  // cout << fToyMargVArray_1Re    ->GetNDim() << " --> " << dimensions_1Re[3]     << " " << dimensions_1Re[4]     << endl;
  // cout << fToyMargVArray_1Rmu   ->GetNDim() << " --> " << dimensions_1Rmu[3]    << " " << dimensions_1Rmu[4]    << endl;
  // cout << fToyMargVArray_RHC1Re ->GetNDim() << " --> " << dimensions_RHC1Re[3]  << " " << dimensions_RHC1Re[4]  << endl;
  // cout << fToyMargVArray_RHC1Rmu->GetNDim() << " --> " << dimensions_RHC1Rmu[3] << " " << dimensions_RHC1Rmu[4] << endl;
  
  //
  
  cout << "////////////////////////////////////////////////////////" << endl;
  cout << endl;
  cout << "Make toys for marginalization and store in a binary file" << endl;
  cout << endl;
  
  cout << endl;
  cout << "Scan the parameter 0 range:" << endl
  << "bin: ["
  << fmin_bin_sub_param[0]
  << ","
  << fmax_bin_sub_param[0]
  << "]"
  << endl
  << fnamefreeparam[0] << ": "
  << "[";
  cout << fhchisq_marg_2D->GetXaxis()->GetBinCenter(fmin_bin_sub_param[0] + 1)
  << ","
  << fhchisq_marg_2D->GetXaxis()->GetBinCenter(fmax_bin_sub_param[0] + 1);
  cout << "]"
  << endl;
  
  cout << endl;
  cout << "Scan the parameter 1 range:" << endl
  << "bin: ["
  << fmin_bin_sub_param[1]
  << ","
  << fmax_bin_sub_param[1]
  << "]"
  << endl
  << fnamefreeparam[1] << ": "
  << "[";
  cout << fhchisq_marg_2D->GetYaxis()->GetBinCenter(fmin_bin_sub_param[1] + 1)
  << ","
  << fhchisq_marg_2D->GetYaxis()->GetBinCenter(fmax_bin_sub_param[1] + 1);
  cout << "]"
  << endl;
  
  // if(!fIsFit2D){
  //   fmin_bin_sub_param[1] = 0.;
  //   fmax_bin_sub_param[1] = 0.;
  // }
  
  cout << endl;
  cout << "fmin_bin_sub_param 0 = " << fmin_bin_sub_param[0]
  << " - fmax_bin_sub_param 0 = " << fmax_bin_sub_param[0] << endl;
  cout << "fmin_bin_sub_param 1 = " << fmin_bin_sub_param[1]
  << " - fmax_bin_sub_param 1 = " << fmax_bin_sub_param[1] << endl;
  cout << endl;
  
  
  // It can extract toys for sub-ranges of the free parameter
  
  for(int iexperiment = 0; iexperiment < fntoys_marg; iexperiment++) {
    
    // Initialize the seed to the experiment number. fSeed is 999
    
    Rnd::Instance()->SetSeeds(fSeed+iexperiment+10,fSeed+iexperiment+20,fSeed+iexperiment+30,fSeed+iexperiment+40,fSeed+iexperiment+50); // stat, syst, osc, fitter, other
    
#ifdef __show_each_experiment__
    cout << endl;
    cout << "iexperiment " << iexperiment << endl;
    cout << "# params 0: " << fnbins_param[0] << endl;
    cout << "# params 1: " << fnbins_param[1] << endl;
#else
    if (iexperiment%100 == 0) { cout << "iexperiment " << iexperiment << endl; }
#endif
    
    fparm->ResetCurrToInitValues(); // Doesn't make sense if SetInitValue is used
                                    // NEW MARG MH
    fparm_nh->ResetCurrToInitValues(); // Doesn't make sense if SetInitValue is used
    fparm_ih->ResetCurrToInitValues(); // Doesn't make sense if SetInitValue is used
                                       //
    
/*    // Extract oscillation parameters to be marginalized
    
    this->SetMargOscParam(fparm_nh);
    this->SetMargOscParam(fparm_ih);*/
    this->SetNoOscillation(fparm_nh);
    this->SetNoOscillation(fparm_ih);
    
#ifdef  __print_oscparam_all__
    
    cout << "toy=" << iexperiment << endl;
    
    cout << "Normal hierarchy parameters:" << endl;
    int np_print=0;
    np_print = fparm_nh->GetNumberParameters();
    for(int ip = 0; ip < np_print; ip++) {
      conv::PhysParamId_t p = fparm_nh->WhichParam(ip);
      if(util::IsSystematic(p)) continue; // not oscillation parameters
      cout << "Osc param " << p << " --> " << util::AsString(p) << " --> " << fparm_nh->CurrValue(p) << endl;
    }
    cout << "Inverted hierarchy parameters:" << endl;
    np_print = fparm_ih->GetNumberParameters();
    for(int ip = 0; ip < np_print; ip++) {
      conv::PhysParamId_t p = fparm_ih->WhichParam(ip);
      if(util::IsSystematic(p)) continue; // not oscillation parameters
      cout << "Osc param " << p << " --> " << util::AsString(p) << " --> " << fparm_ih->CurrValue(p) << endl;
    }
    cout << endl;
#endif
    
    // Set seed after oscillation parameter extraction
    Rnd::Instance()->SetSeeds(
          Rnd::Instance()->RndToyMCstat().GetSeed() + 1,
          Rnd::Instance()->RndToyMCsyst().GetSeed() + 1,
          Rnd::Instance()->RndToyMCosc().GetSeed() + 1,
          Rnd::Instance()->RndFitter().GetSeed() + 1,
          Rnd::Instance()->RndOther().GetSeed()+ 1);
    
    SuperkObservation * skobs = new SuperkObservation;
    
    bool applied_negevents_condition = true;
    bool unphysical_tweak = true;
    while(applied_negevents_condition || unphysical_tweak)
    {
      if(tweak_systematics)
      {
        unphysical_tweak = this->TweakSystematics(skexp, tweak_systematics);
        if(unphysical_tweak)
        { // Repeat loop to get physical systematics
          continue;
        }
      }
      else
      { // no toys
        skexp->ResetSystematicVariations();
        unphysical_tweak = false; // otherwise it doesn't exit the while
      }
      
      // Increment seeds after systematic parameter extraction
      Rnd::Instance()->SetSeeds(
            Rnd::Instance()->RndToyMCstat().GetSeed() + 1,
            Rnd::Instance()->RndToyMCsyst().GetSeed() + 1,
            Rnd::Instance()->RndToyMCosc().GetSeed() + 1,
            Rnd::Instance()->RndFitter().GetSeed() + 1,
            Rnd::Instance()->RndOther().GetSeed() + 1);

      // Marginalize the Mass Hierarchy
      
      conv::NuMassHierarchy_t hier;
      conv::NuMassHierarchy_t hierBar;
      
      int iMH = conv::kBadResult;  // Index of Mass Hierarchy: 0 --> NH, 1 --> IH
      
      // Don't marginalize Mass Hierarchy. Only the assumed one
      if(this->fdoMargMassHierarchy == false){
        
        hier    = Opt::Instance()->gMassHierarchyAssumed;
        hierBar = Opt::Instance()->gMassHierarchyAssumed;
        
        // Used convention in MakeToysMarginalization():
        //  - iMH=0 --> kNmhNormal
        //  - iMH=1 --> kNmhInverted
        if     ( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal   ) iMH = 0;
        else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted ) iMH = 1;
        
        // Mass Hierarchy has to be already defined
        else {
          cerr << endl;
          cerr << "FitMarginalization::MakeToysMarginalization" << endl;
          cerr << "Mass Hierarchy is not defined" << endl;
          cerr << endl;
          util::Die(1);
        }
      }
      
      // Extract randomly the Mass Hierarchy
      else
      iMH = this->RandomizeMH(fProbNormalHierarchy);   // iMH = 0 --> NH, iMH = 1 --> IH
      
      if( iMH==0 ){
        Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("normal");
        hier    = kNmhNormal;
        hierBar = kNmhNormal;
        skexp->GetSuperkMC()->SetMassHierarchy(hier);
        skexp->GetSuperkMC()->SetMassHierarchyBar(hierBar);
        fparm->Copy(*fparm_nh); // NEW MARG MH
      }
      else if( iMH==1 ){
        Opt::Instance()->gMassHierarchyAssumed = util::HierarchyFromString("inverted");
        hier    = kNmhInverted;
        hierBar = kNmhInverted;
        skexp->GetSuperkMC()->SetMassHierarchy(hier);
        skexp->GetSuperkMC()->SetMassHierarchyBar(hierBar);
        fparm->Copy(*fparm_ih); // NEW MARG MH
      }
      
      // If we need to and haven't already save the toy param list
      
      
      if (fShouldSaveToyParamLists && (!fSavedToyParamLists[iexperiment])) {
        // Get Syst tweaks
        PhysParamList true_nuisance_params = skexp->Tweaks();
        PhysParamList * pparm = new PhysParamList (true_nuisance_params);
        
        // Loop over all parameters to check we have enough info
        for (int ip = 0 ; ip < conv::knppmax ; ip++) {
          PhysParamId_t p = (PhysParamId_t)ip;
          if (fparm->Exists(p)) {
            if (not pparm->Contains(p)) {
              // Add parameters that are in fitparamlist but not physparamlist (eg Osc)
              pparm->SetValue(p,fparm->CurrValue(p));
            }
            //~ if (fparm->IsFixed(p)) { // remove fixed parameters
            //~ pparm->Remove(p);
            //~ }
          }
          else { //And remove things not in the fitparamlist
            pparm->Remove(p);
          }
        }
        fSavedToyParamLists[iexperiment] = pparm;
        //~ cerr<<*fparm<<std::endl;
        //~ cerr<<"\n\n===\n\n"<<*pparm<<std::endl;
        //~ exit(-1);
      }

      
#ifdef __FILL_OSCPARAM_TOY_HISTO__
      
      // Fill toy osc parameter histograms
      
      int np = fparm->GetNumberParameters();
      //int npfree = 0;
      for(int ip = 0; ip < np; ip++) {
        conv::PhysParamId_t p = fparm->WhichParam(ip);
        //cout << p << " " << util::AsString(p) << endl;
        
        // NEW FHC+RHC JOINT
        // if(p!=fPhysParamFree[0])
        //   fparm->SetDoMarg(fPhysParamFree[0],true);
        // else{
        //   fparm->SetDoMarg(fPhysParamFree[0],false);  // Do not marginalize if free parameter
        //   continue;
        // }
        // if(fIsFit2D){
        //   if(p!=fPhysParamFree[1])
        //     fparm->SetDoMarg(fPhysParamFree[1],true);
        //   else{
        //     fparm->SetDoMarg(fPhysParamFree[1],false);  // Do not marginalize if free parameter
        //     continue;
        //   }
        // }
        //
        
        if(util::IsSystematic(p)) continue; // not oscillation parameters
        if(fparm->IsFixed(p))     continue; // do not marginalize if fixed
        if(p==conv::kPhypBetaBar) continue; // NEVER marginalize betaBar
        
        
        //cout << "Osc param " << p << " --> " << util::AsString(p) << " --> " << fparm->CurrValue(p) << endl;
        
        
        if(fparm->IsPrior2D(p)){ // Fill it later as correlated (only mix23 and dm23 and/or mix13 and dcp)
          
#ifdef __MARGINALIZE_23PARAMS_CORR__
          if(p==conv::kPhypDm23sq || p==conv::kPhypMixing23) continue;
#endif
#ifdef __MARGINALIZE_TH13_DELTACP_CORR__
          if(p==conv::kPhypCPV || p==conv::kPhypMixing13) continue;
#endif
          cerr << endl;
          cerr << "FitMarginalization::MakeToysMarginalization" << endl;
          cerr << "There are other parameters which are not 2D correlated" << endl;
          cerr << "Now the only 2D correlated parameters are " << endl;
#ifdef __MARGINALIZE_23PARAMS_CORR__
          cerr << util::AsString(conv::kPhypDm23sq)
          << " and "
          << util::AsString(conv::kPhypMixing23)
          << endl;
#endif
#ifdef __MARGINALIZE_TH13_DELTACP_CORR__
          cerr << util::AsString(conv::kPhypCPV)
          << " and "
          << util::AsString(conv::kPhypMixing13)
          << endl;
#endif
          util::Die(1);
        }
        
        
        // NEW FHC+RHC JOINT
        //double usedvalue = fparm->InitValue(p);
        double usedvalue = fparm->CurrValue(p);
        //
        
        fparm->GetToyHisto(p)->Fill(usedvalue);
        
        if(Opt::Instance()->gMassHierarchyAssumed == kNmhNormal){
          
          // NEW FHC+RHC JOINT
          //usedvalue = fparm_nh->InitValue(p);
          usedvalue = fparm_nh->CurrValue(p);
          //
          fparm_nh->GetToyHisto(p)->Fill(usedvalue);
        }
        else if(Opt::Instance()->gMassHierarchyAssumed == kNmhInverted){
          
          // NEW FHC+RHC JOINT
          //usedvalue = fparm_ih->InitValue(p);
          usedvalue = fparm_ih->CurrValue(p);
          //
          fparm_ih->GetToyHisto(p)->Fill(usedvalue);
        }
        
      } // loop over p
      
#if defined(__MARGINALIZE_23PARAMS_CORR__) || defined(__MARGINALIZE_TH13_DELTACP_CORR__)
      double usedvalue_x = conv::kBadResult;
      double usedvalue_y = conv::kBadResult;
#endif
      
#ifdef __MARGINALIZE_23PARAMS_CORR__
      
      // NEW FHC+RHC JOINT
      //usedvalue_x = fparm->InitValue(conv::kPhypMixing23);
      usedvalue_x = fparm->CurrValue(conv::kPhypMixing23);
      //usedvalue_y = fparm->InitValue(conv::kPhypDm23sq);
      usedvalue_y = fparm->CurrValue(conv::kPhypDm23sq);
      //
      
      fparm->GetToyHisto2D(conv::kPhypMixing23)->Fill(usedvalue_x,usedvalue_y);
      fparm->GetToyHisto2D(conv::kPhypDm23sq)->Fill(usedvalue_x,usedvalue_y);
      
      if(Opt::Instance()->gMassHierarchyAssumed == kNmhNormal){
        fparm_nh->GetToyHisto2D(conv::kPhypMixing23)->Fill(usedvalue_x,usedvalue_y);
        fparm_nh->GetToyHisto2D(conv::kPhypDm23sq)->Fill(usedvalue_x,usedvalue_y);
      }
      else if(Opt::Instance()->gMassHierarchyAssumed == kNmhInverted){
        fparm_ih->GetToyHisto2D(conv::kPhypMixing23)->Fill(usedvalue_x,usedvalue_y);
        fparm_ih->GetToyHisto2D(conv::kPhypDm23sq)->Fill(usedvalue_x,usedvalue_y);
      }
#endif
      
#ifdef __MARGINALIZE_TH13_DELTACP_CORR__
      
      // NEW FHC+RHC JOINT
      //usedvalue_x = fparm->InitValue(conv::kPhypMixing13);
      //usedvalue_y = fparm->InitValue(conv::kPhypCPV);
      usedvalue_x = fparm->CurrValue(conv::kPhypMixing13);
      usedvalue_y = fparm->CurrValue(conv::kPhypCPV);
      //
      
      fparm->GetToyHisto2D(conv::kPhypMixing13)->Fill(usedvalue_x,usedvalue_y);
      fparm->GetToyHisto2D(conv::kPhypCPV)->Fill(usedvalue_x,usedvalue_y);
      
      if(Opt::Instance()->gMassHierarchyAssumed == kNmhNormal){
        fparm_nh->GetToyHisto2D(conv::kPhypMixing13)->Fill(usedvalue_x,usedvalue_y);
        fparm_nh->GetToyHisto2D(conv::kPhypCPV)->Fill(usedvalue_x,usedvalue_y);
      }
      else if(Opt::Instance()->gMassHierarchyAssumed == kNmhInverted){
        fparm_ih->GetToyHisto2D(conv::kPhypMixing13)->Fill(usedvalue_x,usedvalue_y);
        fparm_ih->GetToyHisto2D(conv::kPhypCPV)->Fill(usedvalue_x,usedvalue_y);
      }
#endif
      
      if(Opt::Instance()->gMassHierarchyAssumed == kNmhNormal)
      fhToyMH->Fill(0.+conv::kASmallNum);
      else if(Opt::Instance()->gMassHierarchyAssumed == kNmhInverted)
      fhToyMH->Fill(1.+conv::kASmallNum);
      
#endif // __FILL_OSCPARAM_TOY_HISTO__
      
#ifdef __print_masshierarchy__
      cout << endl;
      cout << "iMH = " << iMH << " - MH: " << util::HierarchyAsString(hier)
      << " - MHbar: " << util::HierarchyAsString(hierBar) << endl;
      cout << endl;
#endif
      
#ifdef  __print_oscparam_all_verbose__
      np = fparm->GetNumberParameters();
      for(int ip = 0; ip < np; ip++) {
        conv::PhysParamId_t p = fparm->WhichParam(ip);
        if(util::IsSystematic(p)) continue; // not oscillation parameters
        
        // NEW FHC+RHC JOINT
        //cout << "Osc param " << p << " --> " << util::AsString(p) << " --> " << fparm->InitValue(p) << endl;
        cout << "Osc param " << p << " --> " << util::AsString(p) << " --> " << fparm->CurrValue(p) << endl;
        //
      }
      cout << endl;
#endif
   
      applied_negevents_condition = this->ProduceToys(skobs, skexp, fparm, iexperiment);
    } // applied_negevents_condition

    for (int is = 0 ; is < conv::kns ; is++)
    {
      SuperkSample_t s = conv::kSKSamples[is];
      if (s == conv::kSKsUndefined) continue;
      #ifdef __USE_2D_TEMPLATES__ //1D templates
      TH1D* spectrum_x = skobs->XSpectrum(s);
      spectrum_x->SetDirectory(0);
      spectrum_x->Write(Form("spectrum_x_%s_%d", util::AsStringCompact(s).c_str(), iexperiment));
      delete spectrum_x;
      TH1D* spectrum_y = skobs->YSpectrum(s);
      spectrum_y->SetDirectory(0);
      spectrum_y->Write(Form("spectrum_y_%s_%d", util::AsStringCompact(s).c_str(), iexperiment));
      delete spectrum_y;
      TH2D* spectrum_2D = skobs->TwoDSpectrum(s);
      spectrum_2D->SetDirectory(0);
      spectrum_2D->Write(Form("spectrum_2D_%s_%d", util::AsStringCompact(s).c_str(), iexperiment));
      delete spectrum_2D;
      #endif
    }


#ifdef FILL_TREE_OSCPARAM_TOYS
    treetoys->Fill(); // NEW OSC TREE
#endif
    
    delete skobs;
  } // iexperiment

  toys_file->Close();
  delete toys_file;
  
  cout << endl;
  cout << "Toys for marginalization are extracted" << endl;
  cout << endl;
  
  return;
}
//________________________________________________________________________________________
bool FitMarginalization::TweakSystematics(SuperkExpectation* skexp, int tweak_systematics)
{
   if(tweak_systematics == 1)
   {
      skexp->TweakAllSystematics();
   }
   else if(tweak_systematics == 2)
   {
      skexp->TweakAllBANFFSystematics();
   }
   else if(tweak_systematics == 3)
   {
      skexp->TweakAllUncorrelatedXSecSystematics();
   }
   else if(tweak_systematics == 4)
   {
      skexp->TweakAllSKDetectorFSISystematics();
   }
   else if(tweak_systematics == 5)
   {  // tweak all sources of systematic error but fix 2p2h to its nominal value (only FHC+RHC joint fit)
      #if defined(__RHC_RUNNING__) && defined(__FHC_RUNNING__) && !defined(__FULLY_CORRELATED__)
      skexp->TweakAllSystematics();
      TVectorD twkv(skexp->GetBANFF()->CurrentValues());
      double mecval = skexp->GetBANFF()->InitialValue(conv::kNuisBANFFNorm2p2h);
      int ip = util::ArrIdx::Instance()->FHCpRHCBANFFSystIdx(conv::kNuisBANFFNorm2p2h);
      twkv[ip] = mecval;
      skexp->GetBANFF()->Tweak(true,&twkv);
      #else
      cerr << endl;
      cerr << "FitMarginalization::MakeToysMarginalization" << endl;
      cerr << "Tweaking all systematics but fixing 2p2h is available only for FHC+RHC not fully correlated" << endl;
      cerr << endl;
      util::Die(1);
      #endif
   }
   #ifdef __print_systematic_variations__
   skexp->PrintSystematicVariations();
   #endif

   PhysParamList true_nuisance_params = skexp->Tweaks();
   if(this->fRethrowSystUnphysical)
   {
      return !(valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(true_nuisance_params));
   }
   else
   {
      return false;
   }
}

bool FitMarginalization::ProduceToys(SuperkObservation* skobs,
      SuperkExpectation* skexp, FitParamList* fparm, int iexperiment)
{
   // Update the expectation and produce the mock data
   // Build marginalization distribution w/ systematic fluctations, but no stat fluctuations
   if(skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), false))
   {
      return true;
   }

   // Loop over the used SK samples (This loop can be improved. Loop only over the used SK samples, pre-stored)
   int iscount = 0;
   for(int is = 0; is < fNSKsamples; ++is)
   {
      SuperkSample_t s = this->GetSKSamples(is);
      VALORArrayF* array = this->GetToyMargVArraySample(s);

      #ifndef __USE_2D_TEMPLATES__
      double* obsErec = skobs->ErecoSpectrumArray(s);
      // Build toy Ereco distributions
      for(int ib = 0; ib < util::NumErecoBins(s); ib++)
      {
         float bincont = (float)(obsErec[ib] / fNormPOT[is]);  // NEW PRECISION
         array->SetEntry(bincont, iexperiment, 0, 0, ib);
      }
      #else
      VALORArray* obs = skobs->TwoDSpectrumVArray(s);
      // Build toy 2D distributions (Ereco Vs Q2)
      for(int iereco = 0; iereco < util::NumErecoBins(s); iereco++)
      {
         for(int iq2 = 0; iq2 < util::NumQ2Bins(s); iq2++)
         {
            float bincont = (float) obs->GetEntry(iereco,iq2) / fNormPOT[is];
            array->SetEntry(bincont, iexperiment, 0, 0, iereco, iq2);
         }
      }
      #endif
      iscount++; // index of the toy array
   }

   return false;
}

// NEW 2D TEMPLATES (VARRAY)
VALORArrayF * FitMarginalization::GetToyMargVArraySample(conv::SuperkSample_t s)
{
  switch (s) {
    case ( kSKsFCFV1RingELike     ) : return fToyMargVArray_1Re;
    case ( kSKsFCFV1RingMuLike    ) : return fToyMargVArray_1Rmu;
    case ( kSKsRHCFCFV1RingELike  ) : return fToyMargVArray_RHC1Re;
    case ( kSKsRHCFCFV1RingMuLike ) : return fToyMargVArray_RHC1Rmu;
    case ( kSKsFCFVMultiRingELike ) : return fToyMargVArray_MultiRe;
      
    default                         :
      cerr << endl;
      cerr << "FitMarginalization::GetToyMargVArraySample" << endl;
      cerr << util::AsString(s) << " sample is not implemented" << endl;
      cerr << endl;
      util::Die(1);
  }
  return 0;
}
//________________________________________________________________________________________
VALORArrayF * FitMarginalization::GetToyMargVArraySample_Clone(conv::SuperkSample_t s)
{
  switch (s) {
    case ( kSKsFCFV1RingELike     ) : return fToyMargVArray_Clone_1Re;
    case ( kSKsFCFV1RingMuLike    ) : return fToyMargVArray_Clone_1Rmu;
    case ( kSKsRHCFCFV1RingELike  ) : return fToyMargVArray_Clone_RHC1Re;
    case ( kSKsRHCFCFV1RingMuLike ) : return fToyMargVArray_Clone_RHC1Rmu;
    case ( kSKsFCFVMultiRingELike ) : return fToyMargVArray_Clone_MultiRe;
      
    default                         :
      cerr << endl;
      cerr << "FitMarginalization::GetToyMargVArraySample_Clone" << endl;
      cerr << util::AsString(s) << " sample is not implemented" << endl;
      cerr << endl;
      util::Die(1);
  }
  return 0;
}
//



//________________________________________________________________________________________
void FitMarginalization::FillChisqMarg(SuperkObservation *skobs_data,bool readbinary){
  
  int NtoysTot=0;
  int NparTot0=0;
  int NparTot1=0;
  
  cout << endl;
  cout << "Fill marginalized Chisq in the parameter range:" << endl
  << "bin: ["
  << fmin_bin_sub_param[0]
  << ","
  << fmax_bin_sub_param[0]
  << "]"
  << endl
  << fnamefreeparam[0] << ": "
  << "["
  << fhchisq_marg_2D->GetXaxis()->GetBinCenter(fmin_bin_sub_param[0] + 1)
  << ","
  << fhchisq_marg_2D->GetXaxis()->GetBinCenter(fmax_bin_sub_param[0] + 1)
  << "]"
  << endl;
  cout << "bin: ["
  << fmin_bin_sub_param[1]
  << ","
  << fmax_bin_sub_param[1]
  << "]"
  << endl
  << fnamefreeparam[1] << ": "
  << "["
  << fhchisq_marg_2D->GetYaxis()->GetBinCenter(fmin_bin_sub_param[1] + 1)
  << ","
  << fhchisq_marg_2D->GetYaxis()->GetBinCenter(fmax_bin_sub_param[1] + 1)
  << "]"
  << endl;
  
  
  // It can marginalize sub-ranges of the free parameter
  // fmin_bin_sub_param_1 and fmax_bin_sub_param_1 read from the binary file
  
  
  // Partial Chisq
  
  if(!readbinary){ // old approach: fill range of the read binary file
    NparTot0=0;
    for(int ipar0=fmin_bin_sub_param[0];ipar0<=fmax_bin_sub_param[0];ipar0++){
      NparTot1=0;
      for(int ipar1=fmin_bin_sub_param[1];ipar1<=fmax_bin_sub_param[1];ipar1++){
        
        CalcChisqMargPtRead(skobs_data,ipar0,ipar1,readbinary);
        
        NparTot1++;
      } // ipar1 (free parameter)
      NparTot0++;
    } // ipar0 (free parameter)
  }
  
  
  // Full Chisq loop over files
  
  else{ // new approach: fill full range reading new binary file at each point
    
    NparTot0=0;
    for(int ipar0=fmin_bin_param[0];ipar0<=fmax_bin_param[0];ipar0++){
      NparTot1=0;
      for(int ipar1=fmin_bin_param[1];ipar1<=fmax_bin_param[1];ipar1++){
        
        CalcChisqMargPtRead(skobs_data,ipar0,ipar1,readbinary);
        
        NparTot1++;
      } // ipar1 (free parameter)
      NparTot0++;
    } // ipar0 (free parameter)
  }
  
  NtoysTot = fntoys_marg;
  //NtoysTot = fntoys_marg * 2; // Used when MH was hard-coded set with prior of 0.5
  
  cout << endl;
  cout << "# of scanned free parameters 0: " << NparTot0 << endl;
  cout << "# of scanned free parameters 1: " << NparTot1 << endl;
  cout << "# of toys per point: " << NtoysTot << endl;
  cout << endl;
  
#ifdef __print_chisqmarg__
  cout << endl;
  cout << "Marginalized Chisq distribution:" << endl;
  cout << "fmin_bin_sub_param 0 = " << fmin_bin_sub_param[0] << endl;
  cout << "fmax_bin_sub_param 0 = " << fmax_bin_sub_param[0] << endl;
  cout << "fmin_bin_sub_param 1 = " << fmin_bin_sub_param[1] << endl;
  cout << "fmax_bin_sub_param 1 = " << fmax_bin_sub_param[1] << endl;
  for(int ipar0=fmin_bin_param[0];ipar0<=fmax_bin_param[0];ipar0++){
    for(int ipar1=fmin_bin_param[1];ipar1<=fmax_bin_param[1];ipar1++){
      cout << fChisqMargVArray->GetEntry(ipar0,ipar1) << ", ";
    }
  }
  cout << endl;
#endif
  
  return;
}

//________________________________________________________________________________________
double FitMarginalization::CalcChisqMargPtRead(SuperkObservation *skobs_data,int iterpar0,int iterpar1,bool readbinary,bool allow_cached){ // NEW PRECISION
  
  // If chisq already colculated return the stored value and break
  double chisq = fChisqMargVArray->GetEntry(iterpar0,iterpar1);
  if(allow_cached && chisq > 0.) {return chisq;}
  
  
#ifdef __print_chisqmarg_tstopwatch__
  TStopwatch sw;
#endif
  
  // NEW SPECTRA MARG
  // Set to 0 in case it has been already calculated before, for example when fabs(chisq)<1.E-8.
  //flikeVArray->SetEntry(0.,iterpar0,iterpar1); // not used
  //
  
  fNCallsCalcChisqMargPt++;
  
  
  if(fPhysParamFree[0] == kPhypBetaBar && fPhysParamFree[1] == kPhypUndefined && fteststatistic == "rate-only")
  {
    double bincont =  this->FillExpectedNuEBarAppHisto(skobs_data, iterpar0, readbinary);
    if(bincont<kAVerySmallNum)bincont = kAVerySmallNum;
    fChisqMargVArray->SetEntry(bincont,iterpar0,iterpar1);
    if(!fIsFit2D) fhchisq_marg_1D->SetBinContent(iterpar0+1,bincont);
    return bincont;
  }
  
  // Read binary file
  if(readbinary){
    
    // Get file ID corresponding to the input free parameters
    
    int fileid = GetFileIDRead(iterpar0,iterpar1);
    const char *tag = ftagfile.c_str();
    
    if( (fCurrentFileIDminim==conv::kBadResult) || (fCurrentFileIDminim!=fileid) ){
      //if( (fileid==conv::kBadResult) || (fCurrentFileIDminim!=fileid) ){
      
#ifdef __print_chisqmarg_changebinary__
      cout << endl;
      cout << "Change binary file from " << fCurrentFileIDminim << " to " << fileid << endl;
#endif
      
      fCurrentFileIDminim = fileid;
      this->ReadToysFromBinary(fileid,tag);
    }
  }
  
  //cout << "iterpar0 = " << iterpar0 << " - iterpar1 = " << iterpar1 << "   ";
  //double currval0 = fmin_param[0] + iterpar0*fstep_param[0];
  //double currval1 = fmin_param[1] + iterpar1*fstep_param[1];
  //cout << "currval0 = " << currval0 << endl;
  //cout << "currval1 = " << currval1 << endl;
  
  
  
  // NEW PRECISION
  vector<double> vec_chisq_toy;
  //
  
  
  for(int iexperiment = 0; iexperiment < fntoys_marg; iexperiment++) {
    //cout << "toy=" << iexperiment << " - iterpar0 = " << iterpar0 << " - iterpar1 = " << iterpar1 << endl;
    
#ifdef __print_array_toy_fillchisq__
    cout << endl;
    cout << "Toy " << iexperiment << endl;
    cout << endl;
#endif
    
    // compute the likelihood for a single marginalization toy
    chisq = 0.;
    
    // Loop over the used SK samples
    int iscount=0;
    //for(int is = 0; is < kns; is++) {
    //SuperkSample_t s = kSKSamples[is];
    //if(s == kSKsUndefined) continue;
    for(int is=0;is<fNSKsamples;is++){
      SuperkSample_t s = this->GetSKSamples(is);
      
      
      // NEW FHC+RHC JOINT
      // Skip the sk samples which shouldn't be used in the chisq computation
      bool usesample = this->UseSKSample(is);
      if(!usesample){
        iscount++; // Increase the sk-sample counter, otherwise at the next iteration it takes the wrong sample
        continue;
      }
      //
      
      
      
      if( fteststatistic=="rate+shape"){
        
#ifdef __print_array_toy_fillchisq__
        cout << util::AsString(s) << " --> # Ereco bins = " << util::NumErecoBins(s) << endl;
#ifdef __USE_2D_TEMPLATES__
        // NEW 2D TEMPLATES
        cout << util::AsString(s) << " --> # Q2 bins = " << util::NumQ2Bins(s) << endl;
        //
#endif
        double nexp_sample = 0.;
        double nobs_sample = 0.;
        double chisq_sample = 0.;
#endif
        
        
#ifndef __USE_2D_TEMPLATES__
        
        for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){
          
          // NEW PRECISION
          
          // NEW 2D TEMPLATES (VARRAY)
          //double nexp = (double)fToyMargVArray->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],iscount,ibin);
          double nexp = (double)this->GetToyMargVArraySample(s)->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],ibin);
          //
          
          nexp *= (double)fNormPOT[is];
          double nobs = (double)skobs_data->ErecoSpectrumArray(s)[ibin];
          chisq += valor::util::ChisqPoisson(nobs,nexp); // penalty term is a simple constant --> not needed
                                                         //
          
#ifdef __print_array_toy_fillchisq__
          nexp_sample += nexp;
          nobs_sample += nobs;
          chisq_sample += valor::util::ChisqPoisson(nobs,nexp);
#endif
        } // loop over ereco bins
        
        
#else // __USE_2D_TEMPLATES__
        
        // NEW 2D TEMPLATES
        
        for(int ibinEreco=0;ibinEreco<util::NumErecoBins(s);ibinEreco++){
/*           if(s == kSKsFCFV1RingELike)
           {
              if(kErecoBinEdges1Re[ibinEreco] == 0.100 && kErecoBinEdges1Re[ibinEreco + 1] == 0.150)
              {
                 // Skip bin 0.1-0.15
                 continue;
              }
           }*/

          for(int ibinQ2=0;ibinQ2<util::NumQ2Bins(s);ibinQ2++){
            
            // NEW 2D TEMPLATES (VARRAY)
            //double nexp = (double)fToyMargVArray->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],iscount,ibinEreco,ibinQ2);
            double nexp = (double)this->GetToyMargVArraySample(s)->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],ibinEreco,ibinQ2);
            //
            
            nexp *= (double)fNormPOT[is];
            double nobs = (double)skobs_data->TwoDSpectrumVArray(s)->GetEntry(ibinEreco,ibinQ2); // Maybe it's slower than the simple array --> TO CHECK
            chisq += valor::util::ChisqPoisson(nobs,nexp); // penalty term is a simple constant --> not needed
            
#ifdef __print_array_toy_fillchisq__
            nexp_sample += nexp;
            nobs_sample += nobs;
            chisq_sample += valor::util::ChisqPoisson(nobs,nexp);
#endif
          } // loop over q2 bins
        } // loop over ereco bins
          //
#endif // __USE_2D_TEMPLATES__
        
#ifdef __print_array_toy_fillchisq__
        cout << "# exp sample = " << nexp_sample << " - # obs sample = " << nobs_sample << " --> chisq sample = " << chisq_sample << " --> chisq tot = " << chisq << endl;
#endif
        
      } // rate+shape
      
      else if( fteststatistic=="rate-only" ){
        
        // Take the integral of the ereco distributions (rate-only)
        double nexp = 0.;
        double nobs = 0.;
        
#ifndef __USE_2D_TEMPLATES__
        for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){
          
          // NEW 2D TEMPLATES (VARRAY)
          //nexp += (double)fToyMargVArray->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],iscount,ibin);
          nexp += (double)this->GetToyMargVArraySample(s)->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],ibin);
          //
          
          nobs += (double)skobs_data->ErecoSpectrumArray(s)[ibin];
        }
#else
        // NEW 2D TEMPLATES
        for(int ibinEreco=0;ibinEreco<util::NumErecoBins(s);ibinEreco++){
          for(int ibinQ2=0;ibinQ2<util::NumQ2Bins(s);ibinQ2++){
            
            // NEW 2D TEMPLATES (VARRAY)
            //nexp += (double)fToyMargVArray->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],iscount,ibinEreco,ibinQ2);
            nexp += (double)this->GetToyMargVArraySample(s)->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],ibinEreco,ibinQ2);
            //
            
            nobs += (double)skobs_data->TwoDSpectrumVArray(s)->GetEntry(ibinEreco,ibinQ2);
          }
        }
        //
#endif
        
        nexp *= (double)fNormPOT[is];
        
        //cout << "toy " << iexperiment << " --> " << nexp << endl;
        //cout << endl;
        
        chisq += valor::util::ChisqPoisson(nobs,nexp); // penalty term is a simple constant --> not needed
      } // if(stat-only)
      
      else{
        cerr << endl;
        cerr << "FitMarginalization::CalcChisqMargPtRead" << endl;
        cerr << "The test statistic is not available: " << fteststatistic << endl;
        cerr << endl;
        util::Die(1);
      }
      
      iscount++;
    } // loop SK samples
    
#ifdef __print_array_toy_fillchisq_verbose__
    cout << endl;
    cout << "toy=" << iexperiment << " - iterpar0=" << iterpar0 << " - iterpar1=" << iterpar1 << " : " << endl;
    
    iscount=0;
    //for(int is = 0; is < kns; is++) {
    //SuperkSample_t s = kSKSamples[is];
    //if(s == kSKsUndefined) continue;
    for(int is=0;is<fNSKsamples;is++){
      SuperkSample_t s = this->GetSKSamples(is);
      
      // NEW FHC+RHC JOINT
      // Skip the sk samples which shouldn't be used in the chisq computation
      bool usesample = this->UseSKSample(is);
      if(!usesample){
        iscount++; // Increase the sk-sample counter, otherwise at the next iteration it takes the wrong sample
        continue;
      }
      //
      
      cout << util::AsStringCompact(s) << endl;
      
      double nexp = 0.;
      
#ifndef __USE_2D_TEMPLATES__
      for(int ib = 0; ib < util::NumErecoBins(s); ib++){
        
        // NEW 2D TEMPLATES (VARRAY)
        //double bincont = (double)fToyMargVArray->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],iscount,ib);
        double bincont = (double)this->GetToyMargVArraySample(s)->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],ib);
        //
        
        cout << bincont * (double)fNormPOT[is] << " ";
        nexp += bincont * (double)fNormPOT[is];
      }
#else
      // NEW 2D TEMPLATES
      for(int ibEreco = 0; ibEreco < util::NumErecoBins(s); ibEreco++){
        for(int ibQ2 = 0; ibQ2 < util::NumQ2Bins(s); ibQ2++){
          
          // NEW 2D TEMPLATES (VARRAY)
          //double bincont = (double)fToyMargVArray->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],iscount,ibEreco,ibQ2);
          double bincont = (double)this->GetToyMargVArraySample(s)->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],ibEreco,ibQ2);
          //
          
          cout << bincont * (double)fNormPOT[is] << " ";
          nexp += bincont * (double)fNormPOT[is];
        }
        cout << endl;
      }
      //
#endif
      
      cout << endl;
      cout << "Tot # exp = " << nexp << endl;
      cout << endl;
      
      iscount++;
    } // loop SK samples
    
#endif
    
    // NEW SPECTRA MARG
    
    // // NEW PRECISION
    // double like = exp(-chisq/2.+SCALE_PRECISION) / fntoys_marg;
    // //double like = exp(-chisq/2.) / fntoys_marg;
    // flikeVArray->AddToEntry(like,iterpar0,iterpar1); // this is correct if binary files differ only for free param sub-ranges
    // //
    //
    
    // Get total likelihood of each toy experiment
    double liketoy = exp(-chisq/2.);
    
    if(this->fStoreLikeWeights){
      flikeToyVArray->AddToEntry(liketoy,iterpar0,iterpar1,iexperiment); // this is correct if binary files differ only for free param sub-ranges
    }
    
    
    
    
#ifdef __print_like_toy_fillchisq__
    double val0 = conv::kBadResult;
    double val1 = conv::kBadResult;
    if(!fIsFit2D){
      val0 = fhchisq_marg_1D->GetXaxis()->GetBinCenter(iterpar0+1);
    }
    else{
      val0 = fhchisq_marg_2D->GetXaxis()->GetBinCenter(iterpar0+1);
      val1 = fhchisq_marg_2D->GetYaxis()->GetBinCenter(iterpar1+1);
    }
    
    //cout << "iterpar0 = " << iterpar0 << " - iterpar1 = " << iterpar1 << endl;
    cout << "val0 = " << val0 << " - val1 = " <<  val1 << " --> ";
    cout << "toy " << iexperiment;
    
    cout << endl;
    
    for(int is=0;is<fNSKsamples;is++){
      SuperkSample_t s = this->GetSKSamples(is);
      
      double nexp = 0.;
      for(int ibEreco = 0; ibEreco < util::NumErecoBins(s); ibEreco++){
        for(int ibQ2 = 0; ibQ2 < util::NumQ2Bins(s); ibQ2++){
          double bincont = (double)this->GetToyMargVArraySample(s)->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],ibEreco,ibQ2);
          //cout << bincont * (double)fNormPOT[is] << " ";
          nexp += bincont;
        }
      }
      nexp *= (double)fNormPOT[is];
      cout << util::AsStringCompact(s) << " --> " << nexp << endl;
    }
    
    cout << " -  chisq = " << chisq << " - like = " << liketoy << endl;
    cout << endl;
#endif
    //
    
    
    
#ifdef __print_array_toy_fillchisq__
    cout << "chisq = " << chisq << " - like = " << liketoy << endl;
    cout << endl;
#endif
    
    vec_chisq_toy.push_back(chisq); // NEW PRECISION
    
  } // iexperiment
  
  
  // Now the marginalized likelihood is computed using toys from all the binary files
  // Compute the chisq for each free osc parameter value and fill the array and histo
  
  
  // NEW PRECISION
  
  double bincont = 0.;
  // double bincont = -2*log( flikeVArray->GetEntry(iterpar0,iterpar1) ) + SCALE_PRECISION*2.;
  // //double bincont = -2*log( flikeVArray->GetEntry(iterpar0,iterpar1) );
  // //cout << "Classical chisq = " << bincont << endl;
  
  std::sort(vec_chisq_toy.rbegin(),vec_chisq_toy.rend());
  
  // Take the average value, in order not to exceed both the max and min allowed values
  
  double diff_chisq = (vec_chisq_toy[0] - vec_chisq_toy[vec_chisq_toy.size()-1]) / 2.;
  diff_chisq /= 2.;
  double scale = diff_chisq;
  
  if(diff_chisq > SCALE_PRECISION) scale = (double)SCALE_PRECISION;
  scale = TMath::Min((double)SCALE_PRECISION,scale);
  
  double like_new = 0.;
  for (std::vector<double>::iterator it=vec_chisq_toy.begin(); it!=vec_chisq_toy.end(); ++it){
    double arg = -*it/2. + scale;
    like_new += exp(arg) / fntoys_marg;
  }
  bincont = 2*scale - 2*log(like_new);
  //cout << "New chisq = " << bincont << endl;
  
  //
  
  
  // //
  // // Use the following approximation: log(sum) = log(a0) + log( 1 + sum(a0/ai) ), a0 > a1 > ... > aN
  // // See https://en.wikipedia.org/wiki/List_of_logarithmic_identities
  // //
  // std::sort(vec_chisq_toy.rbegin(),vec_chisq_toy.rend());
  // double sum_chisq = 0.;
  // for (std::vector<double>::iterator it=vec_chisq_toy.begin()+1; it!=vec_chisq_toy.end(); ++it){
  //   //double arg = -0.5 * (*it - vec_chisq_toy[0]);
  //   double arg = -0.5 * (*it - vec_chisq_toy[0]) - SCALE_PRECISION;
  //   //sum_chisq += exp( arg ) / fntoys_marg;
  //   sum_chisq += exp( arg ) * (exp(SCALE_PRECISION)/fntoys_marg);
  //   //cout << "SCALE_PRECISION = " << SCALE_PRECISION << " - vec_chisq_toy[0] = " << vec_chisq_toy[0]
  //   //<< " - *it = " << *it << " - *it - vec_chisq_toy[0] = " << *it - vec_chisq_toy[0] << " - arg = " << arg
  //   //<< " - sum_chisq = " << sum_chisq << endl;
  // }
  // bincont = vec_chisq_toy[0] - 2*log( 1./fntoys_marg + sum_chisq ); // + 2*log(fntoys_marg);
  // // cout << "sum_chisq = " << sum_chisq << endl;
  // // cout << "2*log( 1 + sum_chisq ) = " << 2*log( 1 + sum_chisq ) << endl;
  // // cout << "vec_chisq_toy[0] = " << vec_chisq_toy[0] << endl;
  // // cout << " 2*log(fntoys_marg) = " <<  2*log(fntoys_marg) << endl;
  // cout << "bincont = " << bincont << endl;
  
  
  
  
  
  // NEW PENALTY TERM
  
  // Calculate the penalty term of the free parameters
  
  //cout << endl;
  //cout << "bincont = " << bincont << endl;
  
  double penaltyterm = 0.;
  double freeparval0 = conv::kBadResult;
  double freeparval1 = conv::kBadResult;
  
  if(fusepenalty0 || fusepenalty1){
    if(!fIsFit2D){
      freeparval0 = fhchisq_marg_1D->GetXaxis()->GetBinCenter(iterpar0+1);
      penaltyterm = CalcPenaltyTerm(freeparval0);
    }
    else{
      freeparval0 = fhchisq_marg_2D->GetXaxis()->GetBinCenter(iterpar0+1);
      freeparval1 = fhchisq_marg_2D->GetYaxis()->GetBinCenter(iterpar1+1);
      penaltyterm = CalcPenaltyTerm(freeparval0,freeparval1);
    }
    
    bincont += penaltyterm;
  }
  
  if(!fIsFit2D) fhpenaltyterm_1D->SetBinContent(iterpar0+1,penaltyterm);
  else          fhpenaltyterm_2D->SetBinContent(iterpar0+1,iterpar1+1,penaltyterm);
  
  //cout << "bincont = " << bincont << endl;
  
  //
  
  fChisqMargVArray->SetEntry(bincont,iterpar0,iterpar1);
  if(!fIsFit2D) fhchisq_marg_1D->SetBinContent(iterpar0+1,bincont);
  else          fhchisq_marg_2D->SetBinContent(iterpar0+1,iterpar1+1,bincont);
  
  
  
#ifdef __print_chisqmarg_ptval__
  freeparval0 = fhchisq_marg_2D->GetXaxis()->GetBinCenter(iterpar0+1);
  freeparval1 = fhchisq_marg_2D->GetYaxis()->GetBinCenter(iterpar1+1);
  cout << "ParVal0 = " << freeparval0
  << " - ParVal1 = " << freeparval1
  << " --> chisq = " << bincont
#ifdef __print_chisqmarg_tstopwatch__
  << " --> cpu time = " << sw.CpuTime() << " sec"
  << " --> real time = " << sw.RealTime() << " sec"
#endif
  << endl;
#endif
  
  return bincont;
}
//________________________________________________________________________________________


// NEW PENALTY TERM

void FitMarginalization::UsePenaltyTerm (const bool usepenalty0, const bool usepenalty1){
  fusepenalty0 = usepenalty0;
  fusepenalty1 = usepenalty1;
}
//________________________________________________________________________________________
void FitMarginalization::SetPenaltyTerm(const double mean0, const double sigma0, const double mean1, const double sigma1){
  
  fpenalty_mean0 = mean0;
  fpenalty_mean1 = mean1;
  fpenalty_sigma0 = sigma0;
  fpenalty_sigma1 = sigma1;
  
  // Check if the parameters of the penalty term are set
  if(fusepenalty0){
    if(fpenalty_sigma0<=0.){
      cerr << endl;
      cerr << "FitMarginalization::SetPenaltyTerm" << endl;
      cerr << "Using penalty term for free parameter 0" << endl;
      cerr << "sigma = " << fpenalty_sigma0 << endl;
      cerr << "The error must be >0" << endl;
      cerr << endl;
      util::Die(1);
    }
  }
  if(fusepenalty1){
    if(fpenalty_sigma1<=0.){
      cerr << endl;
      cerr << "FitMarginalization::SetPenaltyTerm" << endl;
      cerr << "Using penalty term for free parameter 1" << endl;
      cerr << "sigma = " << fpenalty_sigma1 << endl;
      cerr << "The error must be >0" << endl;
      cerr << endl;
      util::Die(1);
    }
  }
  
  // Check if the parameters of the penalty term are in the physical region
  conv::PhysParamId_t parfree0 = fPhysParamFree[0];
  conv::PhysParamId_t parfree1 = fPhysParamFree[1];
  bool IsPhysical0 = this->RethrowOscUnphysical(parfree0, fpenalty_mean0);
  bool IsPhysical1 = this->RethrowOscUnphysical(parfree1, fpenalty_mean1);
  if(IsPhysical0){
    cerr << endl;
    cerr << "FitMarginalization::SetPenaltyTerm" << endl;
    cerr << "Penalty term of parameter 0: " << util::AsString(parfree0) << endl;
    cerr << "The central value is in the unphysical region: " << fpenalty_mean0 << endl;
    cerr << "Not allowed --> Set a new value" << endl;
    cerr << endl;
    util::Die(1);
  }
  
  if (parfree1 != conv::kPhypUndefined && fusepenalty1) {
    bool IsPhysical1 = this->RethrowOscUnphysical(parfree1, fpenalty_mean1);
    if(IsPhysical1){
      cerr << endl;
      cerr << "FitMarginalization::SetPenaltyTerm" << endl;
      cerr << "Penalty term of parameter 1: " << util::AsString(parfree1) << endl;
      cerr << "The central value is in the unphysical region: " << fpenalty_mean1 << endl;
      cerr << "Not allowed --> Set a new value" << endl;
      cerr << endl;
      util::Die(1);
    }
  }
  
  // Print parameters of penalty term
  cout << endl;
  if(fusepenalty0){
    cout << "Use Penalty term on parameter 0:" << endl;
    cout << "  - Mean = " << fpenalty_mean0 << endl;
    cout << "  - Sigma = " << fpenalty_sigma0 << endl;
  }
  if(fusepenalty1){
    cout << "Use Penalty term on parameter 1:" << endl;
    cout << "  - Mean = " << fpenalty_mean1 << endl;
    cout << "  - Sigma = " << fpenalty_sigma1 << endl;
  }
  cout << endl;
  
  return;
}
//________________________________________________________________________________________
double FitMarginalization::CalcPenaltyTerm(const double par0, const double par1){
  
  double penaltyterm = 0.;
  
  conv::PhysParamId_t parfree0 = fPhysParamFree[0];
  conv::PhysParamId_t parfree1 = fPhysParamFree[1];
  
  double parval0 = par0;
  double parval1 = par1;
  
  //
  // mixparamchoice = sinsqtheta
  //
  
  if( Opt::Instance()->gOscMixingParamChoice == kNuMixSinsqTheta){
    
    //
    // sinsqtheta running mode is used
    // If sinsq2theta penalty term convert free parameter from sinsqtheta to sinsq2theta
    // Otherwise don't convert anything
    //
    
#ifdef __USE_PRIOR_SINSQ2TH13__
    if(fusepenalty0){
      if(parfree0 == conv::kPhypMixing13 || parfree0 == conv::kPhypMixing13Bar){
        parval0 = util::FromSinsqthToSinsq2th(parval0);
      }
    }
    if(fusepenalty1){
      if(parfree1 == conv::kPhypMixing13 || parfree1 == conv::kPhypMixing13Bar){
        parval1 = util::FromSinsqthToSinsq2th(parval1);
      }
    }
#endif
#ifdef __USE_PRIOR_SINSQ2TH12__
    if(fusepenalty0){
      if(parfree0 == conv::kPhypMixing12 || parfree0 == conv::kPhypMixing12Bar){
        parval0 = util::FromSinsqthToSinsq2th(parval0);
      }
    }
    if(fusepenalty1){
      if(parfree1 == conv::kPhypMixing12 || parfree1 == conv::kPhypMixing12Bar){
        parval1 = util::FromSinsqthToSinsq2th(parval1);
      }
    }
#endif
#ifdef __USE_PRIOR_SINSQ2TH23__
    if(fusepenalty0){
      if(parfree0 == conv::kPhypMixing23 || parfree0 == conv::kPhypMixing23Bar){
        parval0 = util::FromSinsqthToSinsq2th(parval0);
      }
    }
    if(fusepenalty1){
      if(parfree1 == conv::kPhypMixing23 || parfree1 == conv::kPhypMixing23Bar){
        parval1 = util::FromSinsqthToSinsq2th(parval1);
      }
    }
#endif
  }
  
  else if( Opt::Instance()->gOscMixingParamChoice == kNuMixSinsq2Theta){
    
    //
    // sinsq2theta running mode is used
    // If sinsqtheta penalty term convert free parameter from sinsq2theta to sinsqtheta
    // Otherwise don't convert anything
    //
    
#ifndef __USE_PRIOR_SINSQ2TH13__
    if(fusepenalty0){
      if(parfree0 == conv::kPhypMixing13 || parfree0 == conv::kPhypMixing13Bar){
        parval0 = util::FromSinsq2thToSinsqth(parval0);
      }
    }
    if(fusepenalty1){
      if(parfree1 == conv::kPhypMixing13 || parfree1 == conv::kPhypMixing13Bar){
        parval1 = util::FromSinsq2thToSinsqth(parval1);
      }
    }
#endif
#ifndef __USE_PRIOR_SINSQ2TH12__
    if(fusepenalty0){
      if(parfree0 == conv::kPhypMixing12 || parfree0 == conv::kPhypMixing12Bar){
        parval0 = util::FromSinsq2thToSinsqth(parval0);
      }
    }
    if(fusepenalty1){
      if(parfree1 == conv::kPhypMixing12 || parfree1 == conv::kPhypMixing12Bar){
        parval1 = util::FromSinsq2thToSinsqth(parval1);
      }
    }
#endif
#ifndef __USE_PRIOR_SINSQ2TH23__
    if(fusepenalty0){
      if(parfree0 == conv::kPhypMixing23 || parfree0 == conv::kPhypMixing23Bar){
        parval0 = util::FromSinsq2thToSinsqth(parval0);
      }
    }
    if(fusepenalty1){
      if(parfree1 == conv::kPhypMixing23 || parfree1 == conv::kPhypMixing23Bar){
        parval1 = util::FromSinsq2thToSinsqth(parval1);
      }
    }
#endif
  }
  
  double penalty0 = 0.;
  double penalty1 = 0.;
  
  if(fusepenalty0) penalty0 = (parval0 - fpenalty_mean0)*(parval0 - fpenalty_mean0) /  ( fpenalty_sigma0 * fpenalty_sigma0);
  if(fusepenalty1) penalty1 = (parval1 - fpenalty_mean1)*(parval1 - fpenalty_mean1) /  ( fpenalty_sigma1 * fpenalty_sigma1);
  
  penaltyterm = penalty0 + penalty1;
  
  // cout << endl;
  // cout << "Parameter 0: " << util::AsString(parfree0) << endl;
  // cout << "parval0 = " << par0 << endl;
  // cout << "parval0 (converted) = " << parval0 << endl;
  // cout << "Parameter 1: " << util::AsString(parfree1) << endl;
  // cout << "parval1 = " << par1 << endl;
  // cout << "parval1 (converted) = " << parval1 << endl;
  // cout << endl;
  // cout << "penalty0 = " << penalty0 << endl;
  // cout << "penalty1 = " << penalty1 << endl;
  // cout << "penaltyterm = " << penaltyterm << endl;
  // cout << endl;
  
  return penaltyterm;
}

//


//________________________________________________________________________________________
void FitMarginalization::SetTestStatistic(string teststatistic){
  
  // "rate-only" or "rate+shape". Always use the log-likelihood ratio
  
  if( teststatistic!="rate+shape" && teststatistic!="rate-only" ){
    cerr << endl;
    cerr << "FitMarginalization::SetTestStatistic" << endl;
    cerr << "The test statistic is not available: " << teststatistic << endl;
    cerr << endl;
    util::Die(1);
  }
  
  fteststatistic = teststatistic;
}
//________________________________________________________________________________________
void FitMarginalization::WriteToys2Binary(int idfile, const char *tag){
  
  cout << endl;
  cout << "Write marginalization toys in a binary file:" << endl;
  
  char namemh[256];
  if(fdoMargMassHierarchy) sprintf(namemh,"%s","marg");
  else{
    if( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal )
    sprintf(namemh,"%s","NH");
    else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted )
    sprintf(namemh,"%s","IH");
    else{
      cerr << "FitMarginalization::WriteToys2Binary" << endl;
      cerr << "Mass Hierarchy not defined" << endl;
      util::Die(1);
    }
  }
  
  char outname[256];
  sprintf(outname,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"",".bin");
  std::string soutname(outname);
  
  
  // NEW 2D TEMPLATES (VARRAY)
  
  //fToyMargVArray->WriteFile(soutname);
  
  ofstream os_toys(soutname.c_str(),ios::binary|ios::out);
  assert(os_toys);
  
  int Size = conv::kBadResult;
  int NDim = conv::kBadResult;
  
  // FHC e-like
  NDim = this->GetToyMargVArraySample(kSKsFCFV1RingELike)->GetNDim();
  Size = this->GetToyMargVArraySample(kSKsFCFV1RingELike)->GetSize(-1);
  os_toys.write( (char*)&NDim       ,sizeof(int));
  os_toys.write( (char*) this->GetToyMargVArraySample(kSKsFCFV1RingELike)->GetDimensions() ,sizeof(int)*6); // Max # of dimensions is 6 for VALORArray
  os_toys.write( (char*) this->GetToyMargVArraySample(kSKsFCFV1RingELike)->GetArray(), Size*sizeof(float) ); // It's a VALORArrayF
  
  // FHC mu-like
  NDim = this->GetToyMargVArraySample(kSKsFCFV1RingMuLike)->GetNDim();
  Size = this->GetToyMargVArraySample(kSKsFCFV1RingMuLike)->GetSize(-1);
  os_toys.write( (char*)&NDim       ,sizeof(int));
  os_toys.write( (char*) this->GetToyMargVArraySample(kSKsFCFV1RingMuLike)->GetDimensions() ,sizeof(int)*6); // Max # of dimensions is 6 for VALORArray
  os_toys.write( (char*) this->GetToyMargVArraySample(kSKsFCFV1RingMuLike)->GetArray(), Size*sizeof(float) ); // It's a VALORArrayF
  
  // RHC e-like
  NDim = this->GetToyMargVArraySample(kSKsRHCFCFV1RingELike)->GetNDim();
  Size = this->GetToyMargVArraySample(kSKsRHCFCFV1RingELike)->GetSize(-1);
  os_toys.write( (char*)&NDim       ,sizeof(int));
  os_toys.write( (char*) this->GetToyMargVArraySample(kSKsRHCFCFV1RingELike)->GetDimensions() ,sizeof(int)*6); // Max # of dimensions is 6 for VALORArray
  os_toys.write( (char*) this->GetToyMargVArraySample(kSKsRHCFCFV1RingELike)->GetArray(), Size*sizeof(float) ); // It's a VALORArrayF
  
  // RHC mu-like
  NDim = this->GetToyMargVArraySample(kSKsRHCFCFV1RingMuLike)->GetNDim();
  Size = this->GetToyMargVArraySample(kSKsRHCFCFV1RingMuLike)->GetSize(-1);
  os_toys.write( (char*)&NDim       ,sizeof(int));
  os_toys.write( (char*) this->GetToyMargVArraySample(kSKsRHCFCFV1RingMuLike)->GetDimensions() ,sizeof(int)*6); // Max # of dimensions is 6 for VALORArray
  os_toys.write( (char*) this->GetToyMargVArraySample(kSKsRHCFCFV1RingMuLike)->GetArray(), Size*sizeof(float) ); // It's a VALORArrayF
  
  // FHC MuliRing e-like
  NDim = this->GetToyMargVArraySample(kSKsFCFVMultiRingELike)->GetNDim();
  Size = this->GetToyMargVArraySample(kSKsFCFVMultiRingELike)->GetSize(-1);
  os_toys.write( (char*)&NDim       ,sizeof(int));
  os_toys.write( (char*) this->GetToyMargVArraySample(kSKsFCFVMultiRingELike)->GetDimensions() ,sizeof(int)*6); // Max # of dimensions is 6 for VALORArray
  os_toys.write( (char*) this->GetToyMargVArraySample(kSKsFCFVMultiRingELike)->GetArray(), Size*sizeof(float) ); // It's a VALORArrayF
  
  //
  os_toys.close();
  
  
  cout << " - Toy for marginalization are stored in " << soutname << endl;
  
  // Write store sub-range of the free parameter
  
  //int min_bin_sub_param0 = fmin_bin_sub_param[0];
  //int max_bin_sub_param0 = fmax_bin_sub_param[0];
  //int min_bin_sub_param1 = fmin_bin_sub_param[1];
  //int max_bin_sub_param1 = fmax_bin_sub_param[1];
  
  char outname_range[256];
  sprintf(outname_range,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"_subrange",".bin");
  
  std::string soutname_range(outname_range);
  
  ofstream os(soutname_range.c_str(),ios::binary|ios::out);
  assert(os);
  os.write( (char*)&fmin_bin_sub_param[0],sizeof(int));
  os.write( (char*)&fmax_bin_sub_param[0],sizeof(int));
  os.write( (char*)&fmin_bin_sub_param[1],sizeof(int));
  os.write( (char*)&fmax_bin_sub_param[1],sizeof(int));
  os.close();
  
  cout << " - Free param sub-range for marginalization is stored in " << soutname_range << endl;
  cout << endl;
  
  //
  
  sprintf(outname_range,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"_subrange",".txt");
  
  ofstream ostext(outname_range,ios::out);
  ostext << endl;
  
  const char * config_file = T2KOpt::Instance()->fConfigFileName.c_str();
  ostext << "Configuration file: " << config_file << endl;
  
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s = this->GetSKSamples(is);
    ostext << "  - POT ( " << util::AsString(s) << " ) = " << fPOT[is] << endl;
    ostext << "  - Norm POT ( " << util::AsString(s) << " ) = " << fNormPOT[is] << endl;
  }
  
  ostext << endl;
  ostext << "Used sub-range of the free parameter 0: " << fnamefreeparam[0] << endl;
  ostext << "First bin of free param: " << fmin_bin_sub_param[0] << endl;
  ostext << "Last bin of free param: " << fmax_bin_sub_param[0] << endl;
  ostext << "fnamefreeparam0: " << endl
  << " - " << soutname << endl
  << " - " << soutname_range << endl;
  ostext << endl;
  ostext << "Used sub-range of the free parameter 1: " << fnamefreeparam[1] << endl;
  ostext << "First bin of free param: " << fmin_bin_sub_param[1] << endl;
  ostext << "Last bin of free param: " << fmax_bin_sub_param[1] << endl;
  ostext << "fnamefreeparam1: " << endl
  << " - " << soutname << endl
  << " - " << soutname_range << endl;
  ostext << endl;
  ostext << "Mass hierarchy: ";
  if( this->fdoMargMassHierarchy ){
    ostext << " Marginalize" << endl;
    ostext << " Prob Normal Hierarchy = " << this->GetProbNormalHierarchy() << endl; // NEW MARG MH 2
  }
  else ostext << util::AsString(Opt::Instance()->gMassHierarchyAssumed) << endl;
  ostext << endl;
  if(fRethrowSystUnphysical) ostext << "Rethrow unphysical systematic toys" << endl;
  else ostext << "Don't rethrow unphysical systematic toys" << endl;
  ostext << "fntoys_marg = " << fntoys_marg << endl;
  ostext << endl;
  ostext << "Initialize VALORArray fToyMargVArray with:"<< endl;
  ostext << "fntoys_marg = " << fntoys_marg << endl;
  
  ostext << "fNSKsamples = " << fNSKsamples << endl;
  ostext << "Used sub-samples:" << endl;
  int iscount=0;
  //for(int is = 0; is < kns; is++) {
  //SuperkSample_t s = kSKSamples[is];
  //if(s == kSKsUndefined) continue;
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s = this->GetSKSamples(is);
    ostext << "   - " << util::AsStringCompact(s) << endl;
    ostext << "     fnerecobins = " << fnerecobins->GetEntry(iscount) << endl;
#ifdef __USE_2D_TEMPLATES__
    // NEW 2D TEMPLATES
    ostext << "     fnq2bins = " << fnq2bins->GetEntry(iscount) << endl;
    //
#endif
    iscount++;
  }
  ostext << endl;
  ostext << "fnbins_param0 = " << fnbins_param[0] << endl;
  ostext << "fmin_bin_sub_param0 = " << fmin_bin_sub_param[0] << endl;
  ostext << "fmax_bin_sub_param0 = " << fmax_bin_sub_param[0] << endl;
  ostext << "fmin_sub_param0 = " << fmin_sub_param[0] << endl;
  ostext << "fmax_sub_param0 = " << fmax_sub_param[0] << endl;
  ostext << "fmin_param0 = " << fmin_param[0] << endl;
  ostext << "fmax_param0 = " << fmax_param[0] << endl;
  ostext << "fstep_param0 = " << fstep_param[0] << endl;
  ostext << endl;
  ostext << "fnbins_param1 = " << fnbins_param[1] << endl;
  ostext << "fmin_bin_sub_param1 = " << fmin_bin_sub_param[1] << endl;
  ostext << "fmax_bin_sub_param1 = " << fmax_bin_sub_param[1] << endl;
  ostext << "fmin_sub_param1 = " << fmin_sub_param[1] << endl;
  ostext << "fmax_sub_param1 = " << fmax_sub_param[1] << endl;
  ostext << "fmin_param1 = " << fmin_param[1] << endl;
  ostext << "fmax_param1 = " << fmax_param[1] << endl;
  ostext << "fstep_param1 = " << fstep_param[1] << endl;
  ostext << endl;
  ostext << "fhypothesis = " << util::AsString(fhypothesis) << endl;
  ostext << endl;
  ostext << "FitParamList:" << endl;
  ostext << *fparm << endl;
  ostext << endl;
  ostext << "Priors of oscillation parameters are:" << endl;
  
  for (int ip = 0 ; ip < conv::kNuis000 ; ip++) {
    PhysParamId_t p = (PhysParamId_t)ip;
    if (not fparm->Exists(p)) {continue;}
    if (util::IsMixingAngle(p) || util::IsMassSplitting(p)) {
      ostext << " - " << util::AsString(p) << " (NH) : " << endl;
      ostext << "     Prior = " << this->GetPrior_NH(p) << endl;
      ostext << "     Min Value = " << fparm_nh->GetMinValueFlat(p) << endl;
      ostext << "     Max Value = " << fparm_nh->GetMaxValueFlat(p) << endl;
      ostext << "     Mean Value = " << fparm_nh->GetMeanGaus(p) << endl;
      ostext << "     Sigma Value = " << fparm_nh->GetSigmaGaus(p) << endl;
      ostext << " - " << util::AsString(p) << " (IH) : " << endl;
      ostext << "     Prior = " << this->GetPrior_NH(p) << endl;
      ostext << "     Min Value = " << fparm_nh->GetMinValueFlat(p) << endl;
      ostext << "     Max Value = " << fparm_nh->GetMaxValueFlat(p) << endl;
      ostext << "     Mean Value = " << fparm_nh->GetMeanGaus(p) << endl;
      ostext << "     Sigma Value = " << fparm_nh->GetSigmaGaus(p) << endl;
    }
  }
  ostext << endl;
  
  
  ostext.close();
}
//________________________________________________________________________________________

void FitMarginalization::WriteToys2BinaryMerge(const char *tag){
  
  cout << endl;
  cout << "Write marginalization toys in a binary file:" << endl;
  
  char namemh[256];
  if(fdoMargMassHierarchy) sprintf(namemh,"%s","marg");
  else{
    if( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal )
    sprintf(namemh,"%s","NH");
    else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted )
    sprintf(namemh,"%s","IH");
    else{
      cerr << "FitMarginalization::WriteToys2BinaryMerge" << endl;
      cerr << "Mass Hierarchy not defined" << endl;
      util::Die(1);
    }
  }
  
  char outtag[256] = "%s_%s";
  sprintf(outtag,"%s_%s",tag,"Merged");
  
  char outname[256];
  sprintf(outname,fnamefilebinary,outtag,0,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"",".bin");
  
  std::string soutname(outname);
  
  
  // NEW 2D TEMPLATES (VARRAY)
  
  //fToyMargVArray_Clone->WriteFile(soutname);
  
  ofstream os_toys(soutname.c_str(),ios::binary|ios::out);
  assert(os_toys);
  
  int Size = conv::kBadResult;
  int NDim = conv::kBadResult;
  
  // FHC e-like
  NDim = this->GetToyMargVArraySample_Clone(kSKsFCFV1RingELike)->GetNDim();
  Size = this->GetToyMargVArraySample_Clone(kSKsFCFV1RingELike)->GetSize(-1);
  os_toys.write( (char*)&NDim       ,sizeof(int));
  os_toys.write( (char*) this->GetToyMargVArraySample_Clone(kSKsFCFV1RingELike)->GetDimensions() ,sizeof(int)*6); // Max # of dimensions is 6 for VALORArray
  os_toys.write( (char*) this->GetToyMargVArraySample_Clone(kSKsFCFV1RingELike)->GetArray(), Size*sizeof(float) ); // It's a VALORArrayF
  
  // FHC mu-like
  NDim = this->GetToyMargVArraySample_Clone(kSKsFCFV1RingMuLike)->GetNDim();
  Size = this->GetToyMargVArraySample_Clone(kSKsFCFV1RingMuLike)->GetSize(-1);
  os_toys.write( (char*)&NDim       ,sizeof(int));
  os_toys.write( (char*) this->GetToyMargVArraySample_Clone(kSKsFCFV1RingMuLike)->GetDimensions() ,sizeof(int)*6); // Max # of dimensions is 6 for VALORArray
  os_toys.write( (char*) this->GetToyMargVArraySample_Clone(kSKsFCFV1RingMuLike)->GetArray(), Size*sizeof(float) ); // It's a VALORArrayF
  
  // RHC e-like
  NDim = this->GetToyMargVArraySample_Clone(kSKsRHCFCFV1RingELike)->GetNDim();
  Size = this->GetToyMargVArraySample_Clone(kSKsRHCFCFV1RingELike)->GetSize(-1);
  os_toys.write( (char*)&NDim       ,sizeof(int));
  os_toys.write( (char*) this->GetToyMargVArraySample_Clone(kSKsRHCFCFV1RingELike)->GetDimensions() ,sizeof(int)*6); // Max # of dimensions is 6 for VALORArray
  os_toys.write( (char*) this->GetToyMargVArraySample_Clone(kSKsRHCFCFV1RingELike)->GetArray(), Size*sizeof(float) ); // It's a VALORArrayF
  
  // RHC mu-like
  NDim = this->GetToyMargVArraySample_Clone(kSKsRHCFCFV1RingMuLike)->GetNDim();
  Size = this->GetToyMargVArraySample_Clone(kSKsRHCFCFV1RingMuLike)->GetSize(-1);
  os_toys.write( (char*)&NDim       ,sizeof(int));
  os_toys.write( (char*) this->GetToyMargVArraySample_Clone(kSKsRHCFCFV1RingMuLike)->GetDimensions() ,sizeof(int)*6); // Max # of dimensions is 6 for VALORArray
  os_toys.write( (char*) this->GetToyMargVArraySample_Clone(kSKsRHCFCFV1RingMuLike)->GetArray(), Size*sizeof(float) ); // It's a VALORArrayF
  
  // FHC Multi Ring e-like
  NDim = this->GetToyMargVArraySample_Clone(kSKsFCFVMultiRingELike)->GetNDim();
  Size = this->GetToyMargVArraySample_Clone(kSKsFCFVMultiRingELike)->GetSize(-1);
  os_toys.write( (char*)&NDim       ,sizeof(int));
  os_toys.write( (char*) this->GetToyMargVArraySample_Clone(kSKsFCFVMultiRingELike)->GetDimensions() ,sizeof(int)*6); // Max # of dimensions is 6 for VALORArray
  os_toys.write( (char*) this->GetToyMargVArraySample_Clone(kSKsFCFVMultiRingELike)->GetArray(), Size*sizeof(float) ); // It's a VALORArrayF
  
  
  os_toys.close();
  
  //
  
  
  
  // Write store sub-range of the free parameter
  
  //int min_bin_sub_param0 = fmin_bin_sub_param[0];
  //int max_bin_sub_param0 = fmax_bin_sub_param[0];
  //int min_bin_sub_param1 = fmin_bin_sub_param[1];
  //int max_bin_sub_param1 = fmax_bin_sub_param[1];
  
  char outname_range[256];
  sprintf(outname_range,fnamefilebinary,outtag,0,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"_subrange",".bin");
  
  std::string soutname_range(outname_range);
  
  int min_bin0 = 0;
  int max_bin0 = fnbins_param[0] - 1;
  int min_bin1 = 0;
  int max_bin1 = fnbins_param[1] - 1;
  
  ofstream os(soutname_range.c_str(),ios::binary|ios::out);
  assert(os);
  os.write( (char*)&min_bin0,sizeof(int));
  os.write( (char*)&max_bin0,sizeof(int));
  os.write( (char*)&min_bin1,sizeof(int));
  os.write( (char*)&max_bin1,sizeof(int));
  os.close();
  
  cout << " - Free param sub-range for marginalization is stored in " << soutname_range << endl;
  cout << endl;
  
  //
  
  char outname_text[256];
  std::string soutname_text(outname_text);
  sprintf(outname_text,fnamefilebinary,outtag,0,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"",".txt");
  
  ofstream ostext(outname_text,ios::out);
  ostext << "fnamefreeparam0: " << fnamefreeparam[0] << endl;
  ostext << "fnamefreeparam1: " << fnamefreeparam[1] << endl;
  ostext << "Mass hierarchy: ";
  if( this->fdoMargMassHierarchy ){
    ostext << " Marginalize" << endl;
    ostext << " Prob Normal Hierarchy = " << this->GetProbNormalHierarchy() << endl; // NEW MARG MH 2
  }
  else ostext << util::AsString(Opt::Instance()->gMassHierarchyAssumed) << endl;
  ostext << endl;
  if(fRethrowSystUnphysical) ostext << "Rethrow unphysical systematic toys" << endl;
  else ostext << "Don't rethrow unphysical systematic toys" << endl;
  ostext << "# of merged files: " << fNmergedFiles << endl;
  ostext << "fntoys_marg = " << fntoys_marg << endl;
  
  ostext << "fNSKsamples = " << fNSKsamples << endl;
  ostext << "Used sub-samples:" << endl;
  int iscount=0;
  //for(int is = 0; is < kns; is++) {
  //SuperkSample_t s = kSKSamples[is];
  //if(s == kSKsUndefined) continue;
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s = this->GetSKSamples(is);
    ostext << "   - " << util::AsStringCompact(s) << endl;
    ostext << "     fnerecobins = " << fnerecobins->GetEntry(iscount) << endl;
#ifdef __USE_2D_TEMPLATES__
    // NEW 2D TEMPLATES
    ostext << "     fnq2bins = " << fnq2bins->GetEntry(iscount) << endl;
#endif
    iscount++;
  }
  
  ostext << endl;
  ostext << "fnbins_param0 = " << fnbins_param[0] << endl;
  ostext << "min_bin0 = " << min_bin0 << endl;
  ostext << "max_bin0 = " << max_bin0 << endl;
  ostext << "fmin_param[0] = " << fmin_param[0] << endl;
  ostext << "fmax_param[0] = " << fmax_param[0] << endl;
  ostext << "fstep_param[0] = " << fstep_param[0] << endl;
  ostext << endl;
  ostext << "fnbins_param1 = " << fnbins_param[1] << endl;
  ostext << "min_bin1 = " << min_bin1 << endl;
  ostext << "max_bin1 = " << max_bin1 << endl;
  ostext << "fmin_param[1] = " << fmin_param[1] << endl;
  ostext << "fmax_param[1] = " << fmax_param[1] << endl;
  ostext << "fstep_param[1] = " << fstep_param[1] << endl;
  ostext << endl;
  ostext << "fhypothesis = " << util::AsString(fhypothesis) << endl;
  ostext << endl;
  ostext << "FitParamList:" << endl;
  ostext << *fparm << endl;
  
  ostext.close();
}
//________________________________________________________________________________________
void FitMarginalization::ReadToysFromBinary(int idfile, const char *tag){
  
  char namemh[256];
  if(fdoMargMassHierarchy) sprintf(namemh,"%s","marg");
  else{
    if( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal )
    sprintf(namemh,"%s","NH");
    else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted )
    sprintf(namemh,"%s","IH");
    else{
      cerr << "FitMarginalization::ReadToysFromBinary" << endl;
      cerr << "Mass Hierarchy not defined" << endl;
      util::Die(1);
    }
  }
  
  // Write store sub-range of the free parameter
  
  int min_bin_sub_param0=0;
  int max_bin_sub_param0=0;
  int min_bin_sub_param1=0;
  int max_bin_sub_param1=0;
  //int nbins_param0=0;
  //int nbins_param1=0;
  
  char inname_range[256];
  sprintf(inname_range,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"_subrange",".bin");
  
  std::string sinname_range(inname_range);
  
#ifdef __readtoys_debug__
  cout << "Reading free param sub-range from binary file:" << sinname_range << endl;
#endif
  
  ifstream is(sinname_range.c_str(),ios::binary|ios::in);
  assert(is);
  is.read( (char*)&min_bin_sub_param0,sizeof(int));
  is.read( (char*)&max_bin_sub_param0,sizeof(int));
  is.read( (char*)&min_bin_sub_param1,sizeof(int));
  is.read( (char*)&max_bin_sub_param1,sizeof(int));
  is.close();
  
#ifdef __readtoys_debug__
  cout<<"Done!"<<endl;
#endif
  
  fmin_bin_sub_param[0] = min_bin_sub_param0;
  fmax_bin_sub_param[0] = max_bin_sub_param0;
  fmin_bin_sub_param[1] = min_bin_sub_param1;
  fmax_bin_sub_param[1] = max_bin_sub_param1;
  
  fmin_sub_param[0] = fhchisq_marg_2D->GetXaxis()->GetBinCenter(fmin_bin_sub_param[0] + 1);
  fmax_sub_param[0] = fhchisq_marg_2D->GetXaxis()->GetBinCenter(fmax_bin_sub_param[0] + 1);
  fmin_sub_param[1] = fhchisq_marg_2D->GetYaxis()->GetBinCenter(fmin_bin_sub_param[1] + 1);
  fmax_sub_param[1] = fhchisq_marg_2D->GetYaxis()->GetBinCenter(fmax_bin_sub_param[1] + 1);
  
  // Check inputs
  
  if( !is.good() ){
    cerr << endl;
    cerr << "The input file " << sinname_range << " doesn't exist" << endl;
    cerr << "Check the input parameters: " << endl;
    cerr << "  - tag = " << tag << endl;
    cerr << "  - fntoys_marg = " << fntoys_marg << endl;
    cerr << "  - # of pt param 0 = " << fnbins_param[0] << endl;
    cerr << "  - # of pt param 1 = " << fnbins_param[1] << endl;
    //cerr << "  - MH = " << namemh << endl; // NEW MARG MH 2
    //cerr << "  - # of binary files = " << << endl;
    cerr << endl;
    util::Die(1);
  }
  
  // cout << fToyMargVArray->GetSize() << endl;
  // cout << fToyMargVArray->GetNDim() << endl;
  // cout << fToyMargVArray->GetSize(0) << endl; // # toys marginalization
  // cout << fToyMargVArray->GetSize(1) << endl; // # bins free parameter 0
  // cout << fToyMargVArray->GetSize(2) << endl; // # bins free parameter 1
  // cout << fToyMargVArray->GetSize(3) << endl; // # mass hierarchies
  // cout << fToyMargVArray->GetSize(4) << endl; // # sk samples
  // cout << fToyMargVArray->GetSize(5) << endl; // # bins ereco
  
  // int NtoyVArray = fToyMargVArray->GetSize(0); // # toys marginalization
  // int Nbins0VArray = fToyMargVArray->GetSize(1);  // # bins free parameter 0
  // int Nbins1VArray = fToyMargVArray->GetSize(2);  // # bins free parameter 1
  // int NmhVArray = fToyMargVArray->GetSize(3); // # mass hierarchies
  // int NsksampleVArray = fToyMargVArray->GetSize(4); // # sk samples
  // int NbinserecoVArray = fToyMargVArray->GetSize(5); // # bins ereco
  
  // if(NtoyVArray != fntoys_marg){
  //   cerr << endl;
  //   cerr << "# of toys from file = " << NtoyVArray << endl;
  //   cerr << "# of toys from input = " << fntoys_marg << endl;
  //   cerr << "If the chosen binary file is correct, set # of toys to " << NtoyVArray << endl;
  //   cerr << endl;
  //   util::Die(1);
  // }
  
  char inname[256];
  sprintf(inname,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"",".bin");
  
  std::string sinname(inname);
  //cout << sinname << endl;
  
#ifdef __readtoys_debug__
  cout << "Reading toys for marginalization from binary file: " << sinname << endl;
#endif
  
  
  // NEW 2D TEMPLATES (VARRAY)
  
  //fToyMargVArray->ReadFile(sinname);
  
  ifstream is_toys(sinname.c_str(),ios::binary|ios::in);
  assert(is_toys);
  
  int Size = conv::kBadResult;
  int NDim = conv::kBadResult;
  
  // 1Re
  //this->GetToyMargVArraySample(kSKsFCFV1RingELike)->Cleanup();
#ifndef __USE_2D_TEMPLATES__
  this->GetToyMargVArraySample(kSKsFCFV1RingELike)->Resize(1,1,1,1);
#else
  this->GetToyMargVArraySample(kSKsFCFV1RingELike)->Resize(1,1,1,1,1);
#endif
  int ndim_1Re = 0;
  int dims_1Re[6] = {};
  is_toys.read( (char*)&ndim_1Re  ,sizeof(int));
  is_toys.read( (char*)dims_1Re   ,sizeof(int)*6);
  fToyMargVArray_1Re = new VALORArrayF(ndim_1Re,dims_1Re);
  NDim = this->GetToyMargVArraySample(kSKsFCFV1RingELike)->GetNDim();
  Size = this->GetToyMargVArraySample(kSKsFCFV1RingELike)->GetSize(-1);
  is_toys.read( (char*)GetToyMargVArraySample(kSKsFCFV1RingELike)->GetArray(),Size*sizeof(float));
  
  
  // 1Rmu
#ifndef __USE_2D_TEMPLATES__
  this->GetToyMargVArraySample(kSKsFCFV1RingMuLike)->Resize(1,1,1,1);
#else
  this->GetToyMargVArraySample(kSKsFCFV1RingMuLike)->Resize(1,1,1,1,1);
#endif
  int ndim_1Rmu = 0;
  int dims_1Rmu[6] = {};
  is_toys.read( (char*)&ndim_1Rmu  ,sizeof(int));
  is_toys.read( (char*)dims_1Rmu   ,sizeof(int)*6);
  fToyMargVArray_1Rmu = new VALORArrayF(ndim_1Rmu,dims_1Rmu);
  NDim = this->GetToyMargVArraySample(kSKsFCFV1RingMuLike)->GetNDim();
  Size = this->GetToyMargVArraySample(kSKsFCFV1RingMuLike)->GetSize(-1);
  is_toys.read( (char*)GetToyMargVArraySample(kSKsFCFV1RingMuLike)->GetArray(),Size*sizeof(float));
  
  // RHC1Re
#ifndef __USE_2D_TEMPLATES__
  this->GetToyMargVArraySample(kSKsRHCFCFV1RingELike)->Resize(1,1,1,1);
#else
  this->GetToyMargVArraySample(kSKsRHCFCFV1RingELike)->Resize(1,1,1,1,1);
#endif
  int ndim_RHC1Re = 0;
  int dims_RHC1Re[6] = {};
  is_toys.read( (char*)&ndim_RHC1Re  ,sizeof(int));
  is_toys.read( (char*)dims_RHC1Re   ,sizeof(int)*6);
  fToyMargVArray_RHC1Re = new VALORArrayF(ndim_RHC1Re,dims_RHC1Re);
  NDim = this->GetToyMargVArraySample(kSKsRHCFCFV1RingELike)->GetNDim();
  Size = this->GetToyMargVArraySample(kSKsRHCFCFV1RingELike)->GetSize(-1);
  is_toys.read( (char*)GetToyMargVArraySample(kSKsRHCFCFV1RingELike)->GetArray(),Size*sizeof(float));
  
  // RHC1Rmu
#ifndef __USE_2D_TEMPLATES__
  this->GetToyMargVArraySample(kSKsRHCFCFV1RingMuLike)->Resize(1,1,1,1);
#else
  this->GetToyMargVArraySample(kSKsRHCFCFV1RingMuLike)->Resize(1,1,1,1,1);
#endif
  int ndim_RHC1Rmu = 0;
  int dims_RHC1Rmu[6] = {};
  is_toys.read( (char*)&ndim_RHC1Rmu  ,sizeof(int));
  is_toys.read( (char*)dims_RHC1Rmu   ,sizeof(int)*6);
  fToyMargVArray_RHC1Rmu = new VALORArrayF(ndim_RHC1Rmu,dims_RHC1Rmu);
  NDim = this->GetToyMargVArraySample(kSKsRHCFCFV1RingMuLike)->GetNDim();
  Size = this->GetToyMargVArraySample(kSKsRHCFCFV1RingMuLike)->GetSize(-1);
  is_toys.read( (char*)GetToyMargVArraySample(kSKsRHCFCFV1RingMuLike)->GetArray(),Size*sizeof(float));
  
  //FHC MultiRing E
  //  cout << "MultiRe" << endl;
#ifndef __USE_2D_TEMPLATES__
  this->GetToyMargVArraySample(kSKsFCFVMultiRingELike)->Resize(1,1,1,1);
#else
  this->GetToyMargVArraySample(kSKsFCFVMultiRingELike)->Resize(1,1,1,1,1);
#endif
  int ndim_MultiRe = 0;
  int dims_MultiRe[6] = {};
  is_toys.read( (char*)&ndim_MultiRe  ,sizeof(int));
  is_toys.read( (char*)dims_MultiRe   ,sizeof(int)*6);
  fToyMargVArray_MultiRe = new VALORArrayF(ndim_MultiRe,dims_MultiRe);
  NDim = this->GetToyMargVArraySample(kSKsFCFVMultiRingELike)->GetNDim();
  Size = this->GetToyMargVArraySample(kSKsFCFVMultiRingELike)->GetSize(-1);
  is_toys.read( (char*)GetToyMargVArraySample(kSKsFCFVMultiRingELike)->GetArray(),Size*sizeof(float));
  
  
  is_toys.close();
  
  
  //
  
  
  
  
  
#ifdef __readtoys_debug__
  cout << "Done!"<<endl;
#endif
  
  if(fmin_sub_param[0]<(fmin_param[0]-conv::kAVerySmallNum) || fmax_sub_param[0]>(fmax_param[0]+conv::kAVerySmallNum)){
    cerr << endl;
    cerr << "FitMarginalization::ReadToysFromBinary" << endl;
    cerr << "The free parameter 0 sub-range ["
    << fmin_sub_param[0] << "," << fmax_sub_param[0] << "] "
    << "provided by the files "
    << sinname << endl
    << " and "
    << sinname_range << endl
    << " exceeds the total range set for"
    << fnamefreeparam[0]
    << " --> [" << fmin_param[0] << "," << fmax_param[0] << "]"
    << endl;
    util::Die(1);
  }
  
  if(fmin_sub_param[1]<(fmin_param[1]-conv::kAVerySmallNum) || fmax_sub_param[1]>(fmax_param[1]+conv::kAVerySmallNum)){  // w/o conv::kAVerySmallNum it enters the loop if bounds are identical. To check
    cerr << endl;
    cerr << "FitMarginalization::ReadToysFromBinary" << endl;
    cerr << "The free parameter 1 sub-range ["
    << fmin_sub_param[1] << "," << fmax_sub_param[1] << "] "
    << "provided by the files "
    << sinname << endl
    << " and "
    << sinname_range << endl
    << " exceeds the total range set for"
    << fnamefreeparam[1]
    << " --> [" << fmin_param[1] << "," << fmax_param[1] << "]"
    << endl;
    util::Die(1);
  }
  
#ifdef __readtoys_debug__
  cout << "File sub-range: ["
  << fmin_sub_param[0] << "," << fmax_sub_param[0]
  << "]" << endl;
  cout << "File bin sub-range: ["
  << fmin_bin_sub_param[0]+1 << "," << fmax_bin_sub_param[0]+1
  << "]" << endl << endl;
  
  cout << "File sub-range: ["
  << fmin_sub_param[1] << "," << fmax_sub_param[1]
  << "]" << endl;
  cout << "File bin sub-range: ["
  << fmin_bin_sub_param[1]+1 << "," << fmax_bin_sub_param[1]+1
  << "]" << endl << endl;
#endif
  
  return;
}
//________________________________________________________________________________________



// NEW SPECTRA MARG

void FitMarginalization::WriteLikeToys2Binary(const char *tag){
  
  if(this->fStoreLikeWeights){
    
    cout << endl;
    cout << "Write toy likelihoods in a binary file:" << endl;
    
    char namemh[256];
    if(fdoMargMassHierarchy) sprintf(namemh,"%s","marg");
    else{
      if( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal )
      sprintf(namemh,"%s","NH");
      else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted )
      sprintf(namemh,"%s","IH");
      else{
        cerr << "FitMarginalization::WriteLikeToys2Binary" << endl;
        cerr << "Mass Hierarchy not defined" << endl;
        util::Die(1);
      }
    }
    char outname[256];
    sprintf(outname,fnamelikefilebinary,tag,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,".bin");
    std::string soutname(outname);
    if(this->fStoreLikeWeights) flikeToyVArray->WriteFile(soutname);
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::WriteLikeToys2Binary" << endl;
    cerr << "Storing/Reading of toy likelihood weights is not enabled" << endl;
    cerr << endl;
    util::Die(1);
  }
  
}
//________________________________________________________________________________________
void FitMarginalization::ReadLikeToysFromBinary(const char *tag){
  
  if(this->fStoreLikeWeights){
    
    cout << endl;
    cout << "Read fit toy likelihoods from binary file:" << endl;
    
    std::string sinname(tag);
    cout << sinname << endl;
    
    flikeToyVArray->ReadFile(sinname);
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::WriteLikeToys2Binary" << endl;
    cerr << "Storing/Reading of toy likelihood weights is not enabled" << endl;
    cerr << endl;
  }
  
}











// NEW OSC TREE

//________________________________________________________________________________________
void FitMarginalization::WriteOscToys2Tree(int idfile, const char *tag){
  
#ifdef FILL_TREE_OSCPARAM_TOYS
  
  char namemh[256];
  if(fdoMargMassHierarchy) sprintf(namemh,"%s","marg");
  else{
    if( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal )
    sprintf(namemh,"%s","NH");
    else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted )
    sprintf(namemh,"%s","IH");
    else{
      cerr << "FitMarginalization::WriteOscToys2Tree" << endl;
      cerr << "Mass Hierarchy not defined" << endl;
      util::Die(1);
    }
  }
  char outname[256];
  sprintf(outname,fnameoscfiletree,tag,idfile,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,".root");
  
  filetreetoys = new TFile(outname,"RECREATE");
  
  cout << endl;
  cout << "Write toy oscillation parameters to tree in file" << endl;
  filetreetoys->cd();
  treetoys->Write();
  filetreetoys->Close();
  delete filetreetoys;
  filetreetoys = 0;
  
#else
  cout << "FitMarginalization/WriteOscToys2Tree > The tree with toy oscillation parameters is not stored" << endl;
#endif
}
//________________________________________________________________________________________
void FitMarginalization::InitTreeOscToys(){
  treetoys = new TTree("treemargtoys","FitParamList and osc param values of marg toys");
  treetoys->Branch("best_fit_params",      "valor::FitParamList",  &fparm,      1280000,0);
  treetoys->SetDirectory(0);
}











//________________________________________________________________________________________
void FitMarginalization::MinimizeHisto(){
  
  int minbin = conv::kBadResult;
  int minbinX = conv::kBadResult;
  int minbinY = conv::kBadResult;
  int minbinZ = conv::kBadResult;
  
  cout << endl;
  cout << "Find minimum of marginalized chisq" << endl;
  cout << endl;
  
  if(!fIsFit2D){
    fchisq_bf = fhchisq_marg_1D->GetBinContent(fhchisq_marg_1D->GetMinimumBin());
    minbin = fhchisq_marg_1D->GetMinimumBin();
    fparam_bf[0] = fhchisq_marg_1D->GetBinCenter(minbin);
    fparam_bf[1] = conv::kBadResult;
  }
  else{
    fchisq_bf = fhchisq_marg_2D->GetBinContent(fhchisq_marg_2D->GetMinimumBin());
    fhchisq_marg_2D->GetMinimumBin(minbinX,minbinY,minbinZ);
    fparam_bf[0] = fhchisq_marg_2D->GetXaxis()->GetBinCenter(minbinX);
    fparam_bf[1] = fhchisq_marg_2D->GetYaxis()->GetBinCenter(minbinY);
  }
  return;
}
//________________________________________________________________________________________
void FitMarginalization::MinimizeHisto(TH1D *hInChisqMarg){ // NEW PRECISION
  
  int minbin = conv::kBadResult;
  
  fchisq_bf = hInChisqMarg->GetBinContent(hInChisqMarg->GetMinimumBin());
  minbin = hInChisqMarg->GetMinimumBin();
  fparam_bf[0] = hInChisqMarg->GetBinCenter(minbin);
  fparam_bf[1] = conv::kBadResult;
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::MinimizeHisto(TH2D *hInChisqMarg){ // NEW PRECISION
  
  int minbin = conv::kBadResult;
  
  fchisq_bf = hInChisqMarg->GetBinContent(hInChisqMarg->GetMinimumBin());
  minbin = hInChisqMarg->GetMinimumBin();
  fparam_bf[0] = hInChisqMarg->GetXaxis()->GetBinCenter(minbin);
  fparam_bf[1] = hInChisqMarg->GetYaxis()->GetBinCenter(minbin);
  return;
}
//________________________________________________________________________________________
void FitMarginalization::PrintBestFit(){
  
  cout << endl;
  cout << "Best fit ChisqMarg = " << this->GetChisqBestFit() << endl;
  cout << "Best fit param 0: " << fnamefreeparam[0] << " = " << this->GetParamBestFit(0) << endl;
  if(fIsFit2D)
  cout << "Best fit param 1: " << fnamefreeparam[1] << " = " << this->GetParamBestFit(1) << endl;
  cout << endl;
  
  return;
}
//________________________________________________________________________________________
TH1F * FitMarginalization::GetToyMargHisto(SuperkSample_t sksample,int itoy,int ifreepar0,int ifreepar1){  // NEW MARG MH 2
  
  // Get file ID corresponding to the input free parameters
  
  int fileid = GetFileIDRead(ifreepar0,ifreepar1);
  const char *tag = ftagfile.c_str();
  
  if( (fCurrentFileIDminim==conv::kBadResult) || (fCurrentFileIDminim!=fileid) ){
    fCurrentFileIDminim = fileid;
    this->ReadToysFromBinary(fileid,tag);
  }
  
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s = this->GetSKSamples(is);
    if(sksample == s){
      
      // cout << "sksample = " << sksample << endl;
      // cout << "itoy = " << itoy << endl;
      // cout << "ifreepar0 = " << ifreepar0 << endl;
      // cout << "ifreepar1 = " << ifreepar1 << endl;
      // cout << "bincont:" << endl;
      
      for(int ibin=0;ibin<fnerecobins->GetEntry(is);ibin++){
        
        // NEW 2D TEMPLATES (VARRAY)
        //float bincont = fToyMargVArray->GetEntry(itoy,ifreepar0-fmin_bin_sub_param[0],ifreepar1-fmin_bin_sub_param[1],is,ibin);
        float bincont = this->GetToyMargVArraySample(s)->GetEntry(itoy,ifreepar0-fmin_bin_sub_param[0],ifreepar1-fmin_bin_sub_param[1],ibin);
        //
        
        bincont *= (float)fNormPOT[is];
        
        fhToyMargCurr[is]->SetBinContent(ibin+1,bincont);
        //cout << bincont << " ";
      }
      //cout << endl;
      return fhToyMargCurr[is];
    }
  }
  
  cerr << endl;
  cerr << "FitMarginalization::GetToyMargHisto" << endl;
  cerr << "You have chosen a wrong Superk Sample: "
  << util::AsStringCompact(sksample)
  << endl;
  cerr << "The samples used in the analysis are:" << endl;
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s = this->GetSKSamples(is);
    cerr << "  - " << util::AsStringCompact(s) << endl;
  }
  cerr << endl;
  util::Die(1);
}
//________________________________________________________________________________________


// NEW 2D TEMPLATES

TH2F * FitMarginalization::GetToyMargHisto2D(SuperkSample_t sksample,int itoy,int ifreepar0,int ifreepar1){  // NEW MARG MH 2
  
  // Get file ID corresponding to the input free parameters
  
  int fileid = GetFileIDRead(ifreepar0,ifreepar1);
  const char *tag = ftagfile.c_str();
  
  if( (fCurrentFileIDminim==conv::kBadResult) || (fCurrentFileIDminim!=fileid) ){
    fCurrentFileIDminim = fileid;
    this->ReadToysFromBinary(fileid,tag);
  }
  
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s = this->GetSKSamples(is);
    if(sksample == s){
      
      // cout << "sksample = " << sksample << endl;
      // cout << "itoy = " << itoy << endl;
      // cout << "ifreepar0 = " << ifreepar0 << endl;
      // cout << "ifreepar1 = " << ifreepar1 << endl;
      // cout << "bincont:" << endl;
      
      for(int ibinEreco=0;ibinEreco<fnerecobins->GetEntry(is);ibinEreco++){
        for(int ibinQ2=0;ibinQ2<fnq2bins->GetEntry(is);ibinQ2++){
          
          // NEW 2D TEMPLATES (VARRAY)
          //float bincont = fToyMargVArray->GetEntry(itoy,ifreepar0-fmin_bin_sub_param[0],ifreepar1-fmin_bin_sub_param[1],is,ibinEreco,ibinQ2);
          float bincont = this->GetToyMargVArraySample(s)->GetEntry(itoy,ifreepar0-fmin_bin_sub_param[0],ifreepar1-fmin_bin_sub_param[1],ibinEreco,ibinQ2);
          //
          
          bincont *= (float)fNormPOT[is];
          fhToyMargCurr2D[is]->SetBinContent(ibinQ2+1,ibinEreco+1,bincont);
          //cout << bincont << " ";
        }
      }
      //cout << endl;
      return fhToyMargCurr2D[is];
    }
  }
  
  cerr << endl;
  cerr << "FitMarginalization::GetToyMargHisto" << endl;
  cerr << "You have chosen a wrong Superk Sample: "
  << util::AsStringCompact(sksample)
  << endl;
  cerr << "The samples used in the analysis are:" << endl;
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s = this->GetSKSamples(is);
    cerr << "  - " << util::AsStringCompact(s) << endl;
  }
  cerr << endl;
  util::Die(1);
}
//



//________________________________________________________________________________________




// NEW 2D TEMPLATES (VARRAY)
double FitMarginalization::GetNormPOT(conv::SuperkSample_t s) const
{
  
  double nnormpot = conv::kBadResult;
  
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s_loop = this->GetSKSamples(is);
    //fPOT[is];
    if(s_loop == s) nnormpot = fNormPOT[is];
  }
  
  if(nnormpot < 0.){
    cerr << endl;
    cerr << "The normalization of POT sample " << util::AsString(s) << " is negative:" << endl;
    cerr << "fNormPOT = " << nnormpot << endl;
    cerr << endl;
    util::Die(1);
  }
  
  return nnormpot;
}
//________________________________________________________________________________________
double FitMarginalization::GetPOT(conv::SuperkSample_t s) const
{
  
  double npot = conv::kBadResult;
  
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s_loop = this->GetSKSamples(is);
    //fPOT[is];
    if(s_loop == s) npot = fPOT[is];
  }
  
  if(npot < 0.){
    cerr << endl;
    cerr << "The # of POT sample " << util::AsString(s) << " is negative:" << endl;
    cerr << "fNormPOT = " << npot << endl;
    cerr << endl;
    util::Die(1);
  }
  
  return npot;
}
//


//________________________________________________________________________________________
void FitMarginalization::ResetChisqMarg(int * /*min_bin_sub_param*/ ,int * /*max_bin_sub_param*/)
{
  for(int ipar0=fmin_bin_sub_param[0];ipar0<=fmax_bin_sub_param[0];ipar0++){
    for(int ipar1=fmin_bin_sub_param[1];ipar1<=fmax_bin_sub_param[1];ipar1++){
      double bincont = 0;
      fChisqMargVArray->SetEntry(bincont,ipar0,ipar1);
      fhchisq_marg_1D->SetBinContent(ipar0+1,bincont);
      fhchisq_marg_2D->SetBinContent(ipar0+1,ipar1+1,bincont);
      fhchisq_marg_1D_Step->SetBinContent(ipar0+1,bincont);
      fhchisq_marg_2D_Step->SetBinContent(ipar0+1,ipar1+1,bincont);
      fhdeltachisq_marg_1D->SetBinContent(ipar0+1,bincont);
      fhdeltachisq_marg_2D->SetBinContent(ipar0+1,ipar1+1,bincont);
      
      
      // NEW PENALTY TERM
      fhpenaltyterm_1D->SetBinContent(ipar0+1,bincont);
      fhpenaltyterm_2D->SetBinContent(ipar0+1,ipar1+1,bincont);
      //
      
    }
  }
  return;
}
//________________________________________________________________________________________
void FitMarginalization::BuildDeltaChisqMarg()
{
  if(fchisq_bf==conv::kALargeNum){
    cerr << endl;
    cerr << "FitMarginalization::BuildDeltaChisqMarg" << endl;
    cerr << "The best-fit has not been found" << endl;
    cerr << "Can't get the marginalized DeltaChisq distribution" << endl;
    cerr << endl;
    util::Die(1);
  }
  
  // 1D
  if(!fIsFit2D){
    // fhdeltachisq_marg_1D = (TH1D*) fhchisq_marg_1D->Clone("fhdeltachisq_marg_1D"); // NEW PRECISION
    // //fhdeltachisq_marg_1D->Reset(); // don't reset because must contain all the values of hchisq
    // fhdeltachisq_marg_1D->SetName("fhdeltachisq_marg_1D");
    // fhdeltachisq_marg_1D->SetTitle("fhdeltachisq_marg_1D");
    
    for(int ibin=0;ibin<fhdeltachisq_marg_1D->GetNbinsX();ibin++){
      
      //double bincont = fhdeltachisq_marg_1D->GetBinContent(ibin+1);
      double bincont = fhchisq_marg_1D->GetBinContent(ibin+1);
      
      if(bincont>0.){
        bincont -= this->GetChisqBestFit();
      }
      
      if(bincont<0.){
        cerr << endl;
        cerr << "FitMarginalization::BuildDeltaChisqMarg" << endl;
        cerr << "Delta Chisq content is negative" << endl;
        cerr << "   - X = " << fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(ibin+1);
        cerr << "   - chisq = " << bincont << endl;
        cerr << endl;
        util::Die(1);
      }
      fhdeltachisq_marg_1D->SetBinContent(ibin+1,bincont);
    }
  }
  
  // 2D
  else{
    fhdeltachisq_marg_2D = (TH2D*) fhchisq_marg_2D->Clone("fhdeltachisq_marg_2D"); // NEW PRECISION
    
    for(int ibinX=0;ibinX<fhdeltachisq_marg_2D->GetNbinsX();ibinX++){
      for(int ibinY=0;ibinY<fhdeltachisq_marg_2D->GetNbinsY();ibinY++){
        
        //double bincont = fhdeltachisq_marg_2D->GetBinContent(ibinX+1,ibinY+1);
        double bincont = fhchisq_marg_2D->GetBinContent(ibinX+1,ibinY+1);
        if(bincont>0.){
          bincont -= this->GetChisqBestFit();
        }
        
        if(bincont<0.){
          cerr << endl;
          cerr << "FitMarginalization::BuildDeltaChisqMarg" << endl;
          cerr << "Delta Chisq content is negative" << endl;
          cerr << "   - X = " << fhdeltachisq_marg_2D->GetXaxis()->GetBinCenter(ibinX+1);
          cerr << "   - Y = " << fhdeltachisq_marg_2D->GetYaxis()->GetBinCenter(ibinY+1);
          cerr << "   - chisq = " << bincont << endl;
          cerr << endl;
          util::Die(1);
        }
        fhdeltachisq_marg_2D->SetBinContent(ibinX+1,ibinY+1,bincont);
      }
    }
  }
  
  return;
}


// //________________________________________________________________________________________
// void FitMarginalization::SetPrior(PhysParamId_t physpar,string priorname,double par1,double par2){

//   fparm->SetPriorMode(physpar,priorname);
//   // NEW MARG MH
//   fparm_nh->SetPriorMode(physpar,priorname);
//   fparm_ih->SetPriorMode(physpar,priorname);
//   //

//   if(priorname=="flat"){
//     fparm->SetMinValueFlat(physpar,par1);
//     fparm->SetMaxValueFlat(physpar,par2);
//     // NEW MARG MH
//     fparm_nh->SetMinValueFlat(physpar,par1);
//     fparm_nh->SetMaxValueFlat(physpar,par2);
//     fparm_ih->SetMinValueFlat(physpar,par1);
//     fparm_ih->SetMaxValueFlat(physpar,par2);
//     //
//   }
//   else if(priorname=="gaus"){
//     fparm->SetMeanGaus(physpar,par1);
//     fparm->SetSigmaGaus(physpar,par2);
//     // NEW MARG MH
//     fparm_nh->SetMeanGaus(physpar,par1);
//     fparm_nh->SetSigmaGaus(physpar,par2);
//     fparm_ih->SetMeanGaus(physpar,par1);
//     fparm_ih->SetSigmaGaus(physpar,par2);
//     //
//   }
//   else if(priorname=="user"){
//     cerr << endl;
//     cerr << "FitMarginalization::SetPrior" << endl;
//     cerr << "You chose a user prior but you didn't give any histogram as input" << endl;
//     util::Die(1);
//   }
//   else{
//     cerr << endl;
//     cerr << "FitMarginalization::SetPrior" << endl;
//     cerr << "You chose the wrong prior option" << endl;
//     util::Die(1);
//   }

//   return;
// }
//________________________________________________________________________________________
void FitMarginalization::SetPrior_NH(PhysParamId_t physpar,string priorname,double par1,double par2){
  
  // NEW OSC PRIOR
  cout << "Set NH prior of " << util::AsString(physpar) << " --> " << priorname << endl;
  
#ifdef __USE_PRIOR_SINSQ2TH13__
  if(physpar == conv::kPhypMixing13 || physpar == conv::kPhypMixing13Bar){
    cout << "The prior is defined on sinsq2th13 --> convert to sinsqth13" << endl;
  }
#endif
#ifdef __USE_PRIOR_SINSQ2TH23__
  if(physpar == conv::kPhypMixing23 || physpar == conv::kPhypMixing23Bar){
    cout << "The prior is defined on sinsq2th23 --> convert to sinsqth23" << endl;
  }
#endif
#ifdef __USE_PRIOR_SINSQ2TH12__
  if(physpar == conv::kPhypMixing12 || physpar == conv::kPhypMixing12Bar){
    cout << "The prior is defined on sinsq2th12 --> convert to sinsqth12" << endl;
  }
#endif
  if(priorname == "flat"){
    cout << " - low = " << par1 << endl;
    cout << " - up = " << par2 << endl;
  }
  else if(priorname == "gaus"){
    cout << " - mean = " << par1 << endl;
    cout << " - sigma = " << par2 << endl;
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_NH" << endl;
    cerr << "Used prior is: " << priorname << endl;
    cerr << "Prior not well defined for " << util::AsString(physpar) << endl;
    util::Die(1);
  }
  cout << endl;
  //
  
  fparm->SetPriorMode(physpar,priorname);
  // NEW MARG MH
  fparm_nh->SetPriorMode(physpar,priorname);
  //
  
  
  // NEW FHC+RHC JOINT
  
  double nbins = 100;
  double min = this->GetMinOscParam(physpar);
  double max = this->GetMaxOscParam(physpar);
  
  // NEW OSC PRIOR
  if(physpar == conv::kPhypMixing13 || physpar == conv::kPhypMixing13Bar){ // Define a smaller range for mixing13 histogram
    min = 0.;
    max = 0.4;
    nbins = 500;
  }
  //
  
  char name[256];
  sprintf(name,"h_%s",util::AsString(physpar).c_str());
  TH1D *h = new TH1D(name,name,nbins,min,max);
  fparm->SetToyHisto(physpar,h); // empty histogram
  fparm_nh->SetToyHisto(physpar,h); // empty histogram
                                    //cout << util::AsString(physpar) << " - " << name << " - min = " << min << " - max = " << max << endl;
  
  //
  
  
  if(priorname=="flat"){
    fparm->SetMinValueFlat(physpar,par1);
    fparm->SetMaxValueFlat(physpar,par2);
    // NEW MARG MH
    fparm_nh->SetMinValueFlat(physpar,par1);
    fparm_nh->SetMaxValueFlat(physpar,par2);
    //
  }
  else if(priorname=="gaus"){
    fparm->SetMeanGaus(physpar,par1);
    fparm->SetSigmaGaus(physpar,par2);
    // NEW MARG MH
    fparm_nh->SetMeanGaus(physpar,par1);
    fparm_nh->SetSigmaGaus(physpar,par2);
    //
  }
  else if(priorname=="user"){
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_NH" << endl;
    cerr << "You chose a user prior but you didn't give any histogram as input" << endl;
    util::Die(1);
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_NH" << endl;
    cerr << "You chose the wrong prior option" << endl;
    util::Die(1);
  }
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetPrior_IH(PhysParamId_t physpar,string priorname,double par1,double par2){
  
  // NEW OSC PRIOR
  cout << "Set IH prior of " << util::AsString(physpar) << " --> " << priorname << endl;
  
  // TODO: Running in sinsqtheta mode. The sinsq2theta running mode is not available yet.
  
#ifdef __USE_PRIOR_SINSQ2TH13__
  if(physpar == conv::kPhypMixing13 || physpar == conv::kPhypMixing13Bar){
    cout << "The prior is defined on sinsq2th13 --> convert to sinsqth13" << endl;
  }
#endif
#ifdef __USE_PRIOR_SINSQ2TH23__
  if(physpar == conv::kPhypMixing23 || physpar == conv::kPhypMixing23Bar){
    cout << "The prior is defined on sinsq2th23 --> convert to sinsqth23" << endl;
  }
#endif
#ifdef __USE_PRIOR_SINSQ2TH12__
  if(physpar == conv::kPhypMixing12 || physpar == conv::kPhypMixing12Bar){
    cout << "The prior is defined on sinsq2th12 --> convert to sinsqth12" << endl;
  }
#endif
  if(priorname == "flat"){
    cout << " - low = " << par1 << endl;
    cout << " - up = " << par2 << endl;
  }
  else if(priorname == "gaus"){
    cout << " - mean = " << par1 << endl;
    cout << " - sigma = " << par2 << endl;
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_NH" << endl;
    cerr << "Used prior is: " << priorname << endl;
    cerr << "Prior not well defined for " << util::AsString(physpar) << endl;
    util::Die(1);
  }
  cout << endl;
  //
  
  // NEW MARG MH
  fparm_ih->SetPriorMode(physpar,priorname);
  //
  
  
  // NEW FHC+RHC JOINT
  
  double nbins = 100;
  double min = this->GetMinOscParam(physpar);
  double max = this->GetMaxOscParam(physpar);
  
  // NEW OSC PRIOR
  if(physpar == conv::kPhypMixing13 || physpar == conv::kPhypMixing13Bar){ // Define a smaller range for mixing13 histogram
    min = 0.;
    max = 0.4;
    nbins = 500;
  }
  //
  
  char name[256];
  sprintf(name,"h_%s_ih",util::AsString(physpar).c_str());
  TH1D *h = new TH1D(name,name,nbins,min,max);
  fparm_ih->SetToyHisto(physpar,h); // empty histogram
  
  //cout << util::AsString(physpar) << " - " << name << " - min = " << min << " - max = " << max << endl;
  
  //
  
  
  if(priorname=="flat"){
    // NEW MARG MH
    fparm_ih->SetMinValueFlat(physpar,par1);
    fparm_ih->SetMaxValueFlat(physpar,par2);
    //
  }
  else if(priorname=="gaus"){
    // NEW MARG MH
    fparm_ih->SetMeanGaus(physpar,par1);
    fparm_ih->SetSigmaGaus(physpar,par2);
    //
  }
  else if(priorname=="user"){
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_IH" << endl;
    cerr << "You chose a user prior but you didn't give any histogram as input" << endl;
    util::Die(1);
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_IH" << endl;
    cerr << "You chose the wrong prior option" << endl;
    util::Die(1);
  }
  
  return;
}





// //________________________________________________________________________________________
// void FitMarginalization::SetPrior(PhysParamId_t physpar,string priorname,TH1D *hprior){
//   fparm->SetPriorMode(physpar,priorname);
//   // NEW MARG MH
//   fparm_nh->SetPriorMode(physpar,priorname);
//   fparm_ih->SetPriorMode(physpar,priorname);
//   //

//   if(priorname=="user"){
//     fparm->SetPriorHisto(physpar,hprior);
//     fparm->SetToyHisto(physpar,hprior);
//     // NEW MARG MH
//     fparm_nh->SetPriorHisto(physpar,hprior);
//     fparm_nh->SetToyHisto(physpar,hprior);
//     fparm_ih->SetPriorHisto(physpar,hprior);
//     fparm_ih->SetToyHisto(physpar,hprior);
//     //
//   }
//   else if(priorname=="flat" || priorname=="gaus"){
//     cerr << endl;
//     cerr << "FitMarginalization::SetPrior" << endl;
//     cerr << "You chose a " << priorname << " prior, but you gave as input an histogram" << endl;
//     util::Die(1);
//   }
//   else{
//     cerr << endl;
//     cerr << "FitMarginalization::SetPrior" << endl;
//     cerr << "You chose the wrong prior option" << endl;
//     util::Die(1);
//   }
//   return;
// }
//________________________________________________________________________________________
void FitMarginalization::SetPrior_NH(PhysParamId_t physpar,string priorname,TH1D *hprior){
  fparm->SetPriorMode(physpar,priorname);
  // NEW MARG MH
  fparm_nh->SetPriorMode(physpar,priorname);
  //
  
  if(priorname=="user"){
    fparm->SetPriorHisto(physpar,hprior);
    fparm->SetToyHisto(physpar,hprior);
    // NEW MARG MH
    fparm_nh->SetPriorHisto(physpar,hprior);
    fparm_nh->SetToyHisto(physpar,hprior);
    //
  }
  else if(priorname=="flat" || priorname=="gaus"){
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_NH" << endl;
    cerr << "You chose a " << priorname << " prior, but you gave as input an histogram" << endl;
    util::Die(1);
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_NH" << endl;
    cerr << "You chose the wrong prior option" << endl;
    util::Die(1);
  }
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetPrior_IH(PhysParamId_t physpar,string priorname,TH1D *hprior){
  // NEW MARG MH
  fparm_ih->SetPriorMode(physpar,priorname);
  //
  
  if(priorname=="user"){
    // NEW MARG MH
    fparm_ih->SetPriorHisto(physpar,hprior);
    fparm_ih->SetToyHisto(physpar,hprior);
    //
  }
  else if(priorname=="flat" || priorname=="gaus"){
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_IH" << endl;
    cerr << "You chose a " << priorname << " prior, but you gave as input an histogram" << endl;
    util::Die(1);
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_IH" << endl;
    cerr << "You chose the wrong prior option" << endl;
    util::Die(1);
  }
  return;
}






// //________________________________________________________________________________________
// void FitMarginalization::SetPrior(PhysParamId_t physpar,string priorname,TH2D *hprior){
//   fparm->SetPriorMode(physpar,priorname);
//   // NEW MARG MH
//   fparm_nh->SetPriorMode(physpar,priorname);
//   fparm_ih->SetPriorMode(physpar,priorname);
//   //

//   if(priorname=="user"){
//     fparm->SetPriorHisto2D(physpar,hprior);
//     fparm->SetToyHisto2D(physpar,hprior);
//     // NEW MARG MH
//     fparm_nh->SetPriorHisto2D(physpar,hprior);
//     fparm_nh->SetToyHisto2D(physpar,hprior);
//     fparm_ih->SetPriorHisto2D(physpar,hprior);
//     fparm_ih->SetToyHisto2D(physpar,hprior);
//     //
//   }
//   else if(priorname=="flat" || priorname=="gaus"){
//     cerr << endl;
//     cerr << "FitMarginalization::SetPrior" << endl;
//     cerr << "You chose a " << priorname << " prior, but you gave as input an histogram" << endl;
//     util::Die(1);
//   }
//   else{
//     cerr << endl;
//     cerr << "FitMarginalization::SetPrior" << endl;
//     cerr << "You chose the wrong prior option" << endl;
//     util::Die(1);
//   }
//   return;
// }
//________________________________________________________________________________________
void FitMarginalization::SetPrior_NH(PhysParamId_t physpar,string priorname,TH2D *hprior){
  fparm->SetPriorMode(physpar,priorname);
  // NEW MARG MH
  fparm_nh->SetPriorMode(physpar,priorname);
  //
  
  if(priorname=="user"){
    fparm->SetPriorHisto2D(physpar,hprior);
    fparm->SetToyHisto2D(physpar,hprior);
    // NEW MARG MH
    fparm_nh->SetPriorHisto2D(physpar,hprior);
    fparm_nh->SetToyHisto2D(physpar,hprior);
    //
  }
  else if(priorname=="flat" || priorname=="gaus"){
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_NH" << endl;
    cerr << "You chose a " << priorname << " prior, but you gave as input an histogram" << endl;
    util::Die(1);
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_NH" << endl;
    cerr << "You chose the wrong prior option" << endl;
    util::Die(1);
  }
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetPrior_IH(PhysParamId_t physpar,string priorname,TH2D *hprior){
  // NEW MARG MH
  fparm_ih->SetPriorMode(physpar,priorname);
  //
  
  if(priorname=="user"){
    // NEW MARG MH
    fparm_ih->SetPriorHisto2D(physpar,hprior);
    fparm_ih->SetToyHisto2D(physpar,hprior);
    //
  }
  else if(priorname=="flat" || priorname=="gaus"){
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_IH" << endl;
    cerr << "You chose a " << priorname << " prior, but you gave as input an histogram" << endl;
    util::Die(1);
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::SetPrior_IH" << endl;
    cerr << "You chose the wrong prior option" << endl;
    util::Die(1);
  }
  return;
}

//________________________________________________________________________________________
string FitMarginalization::GetPrior_NH(PhysParamId_t physpar) const
{
  return fparm_nh->GetPriorModeName(physpar);
}
//________________________________________________________________________________________
string FitMarginalization::GetPrior_IH(PhysParamId_t physpar) const
{
  return fparm_ih->GetPriorModeName(physpar);
}


// NEW MARG MH

//________________________________________________________________________________________
void FitMarginalization::SetMargOscParam(FitParamList * infparm)
{
  int Nparam_with_2Dprior = 0;
  
  // Loop over the fit parameter list
  int np = infparm->GetNumberParameters();
  //int npfree = 0;
  for(int ip = 0; ip < np; ip++) {
    conv::PhysParamId_t p = infparm->WhichParam(ip);
    //cout << p << " " << util::AsString(p) << endl;
    
    //bool fdoMarginalize = true;
    
    if(p!=fPhysParamFree[0])
    infparm->SetDoMarg(p,true);
    else{
      infparm->SetDoMarg(p,false);  // Do not marginalize if free parameter
      continue;
    }
    
    if(fIsFit2D){
      if(p!=fPhysParamFree[1])
      infparm->SetDoMarg(p,true);
      else{
        infparm->SetDoMarg(p,false);  // Do not marginalize if free parameter
        continue;
      }
    }
    
    if(util::IsSystematic(p)) continue; // not oscillation parameters
    if(infparm->IsFixed(p))   continue; // do not marginalize if fixed
    
    // Now extract the oscillation parameter to marginalize from 1D prior
    
    if(infparm->IsPrior2D(p)){
      //cout << p << " " << util::AsString(p) << " --> 2D prior --> skip it!" << endl;
      Nparam_with_2Dprior++;
      continue;
    }
    
    double currval = conv::kBadResult;
    
    bool rethrow=true;
    while(rethrow){
      //currval = Randomize(p);
      currval = Randomize(p,infparm); // NEW MARG MH
      rethrow = RethrowOscUnphysical(p,currval);
      if(rethrow) cout << "Rethrow " << util::AsString(p) << " " << currval << endl;
    }
    
#ifdef  __print_oscparam_marginalized__
    cout << "Marginalize " << p << " " << util::AsString(p) << " --> " << currval << endl;
#endif
    
    // Update current values of parameters
    // NEW FHC+RHC JOINT
    infparm->Update(p,currval);
    //infparm->SetInitValue(p,currval);
    //
    
  } // loop over parameters
  
  
  
  // Extract correlated parameters
  
#if defined(__MARGINALIZE_TH13_DELTACP_CORR__) && defined(__MARGINALIZE_23PARAMS_CORR__)
  
  if(Nparam_with_2Dprior!=4){
    cerr << endl;
    cerr << "FitMarginalization::SetMargOscParam" << endl;
    cerr << "# of marginalization parameters with 2D prior != 4" << endl;
    cerr << "Nparam_with_2Dprior = " << Nparam_with_2Dprior << endl;
    cerr << "We don't expect parameters with 2D priors for both 23 parameters and theta13 Vs dcp" << endl;
    cerr << endl;
    util::Die(1);
  }
  
  if( !( infparm->IsPrior2D(conv::kPhypDm23sq) ) ){
    cerr << endl;
    cerr << util::AsString(conv::kPhypDm23sq)
    << " don't have a 2D prior" << endl;
    cerr << "There is something wrong!" << endl;
    cerr << endl;
    util::Die(1);
  }
  if( !( infparm->IsPrior2D(conv::kPhypMixing23) ) ){
    cerr << endl;
    cerr << util::AsString(conv::kPhypMixing23)
    << " don't have a 2D prior" << endl;
    cerr << "There is something wrong!" << endl;
    cerr << endl;
    util::Die(1);
  }
  if( !( infparm->IsPrior2D(conv::kPhypCPV) ) ){
    cerr << endl;
    cerr << util::AsString(conv::kPhypCPV)
    << " don't have a 2D prior" << endl;
    cerr << "There is something wrong!" << endl;
    cerr << endl;
    util::Die(1);
  }
  if( !( infparm->IsPrior2D(conv::kPhypMixing13) ) ){
    cerr << endl;
    cerr << util::AsString(conv::kPhypMixing13)
    << " don't have a 2D prior" << endl;
    cerr << "There is something wrong!" << endl;
    cerr << endl;
    util::Die(1);
  }
  
  
#elif defined(__MARGINALIZE_TH13_DELTACP_CORR__) && !defined(__MARGINALIZE_23PARAMS_CORR__)
  
  if(Nparam_with_2Dprior!=2){
    cerr << endl;
    cerr << "FitMarginalization::SetMargOscParam" << endl;
    cerr << "# of marginalization parameters with 2D prior != 2" << endl;
    cerr << "Nparam_with_2Dprior = " << Nparam_with_2Dprior << endl;
    cerr << "We expect 2 parameters with 2D priors for theta13 Vs dcp" << endl;
    cerr << endl;
    util::Die(1);
  }
  
  if( !( infparm->IsPrior2D(conv::kPhypCPV) ) ){
    cerr << endl;
    cerr << util::AsString(conv::kPhypCPV)
    << " don't have a 2D prior" << endl;
    cerr << "There is something wrong!" << endl;
    cerr << endl;
    util::Die(1);
  }
  if( !( infparm->IsPrior2D(conv::kPhypMixing13) ) ){
    cerr << endl;
    cerr << util::AsString(conv::kPhypMixing13)
    << " don't have a 2D prior" << endl;
    cerr << "There is something wrong!" << endl;
    cerr << endl;
    util::Die(1);
  }
  
#elif !defined(__MARGINALIZE_TH13_DELTACP_CORR__) && defined(__MARGINALIZE_23PARAMS_CORR__)
  if(Nparam_with_2Dprior!=2){
    cerr << endl;
    cerr << "FitMarginalization::SetMargOscParam" << endl;
    cerr << "# of marginalization parameters with 2D prior != 2" << endl;
    cerr << "Nparam_with_2Dprior = " << Nparam_with_2Dprior << endl;
    cerr << "We expect 2 parameters with 2D priors for 23 parameters" << endl;
    cerr << endl;
    util::Die(1);
  }
  
  if( !( infparm->IsPrior2D(conv::kPhypDm23sq) ) ){
    cerr << endl;
    cerr << util::AsString(conv::kPhypDm23sq)
    << " don't have a 2D prior" << endl;
    cerr << "There is something wrong!" << endl;
    cerr << endl;
    util::Die(1);
  }
  if( !( infparm->IsPrior2D(conv::kPhypMixing23) ) ){
    cerr << endl;
    cerr << util::AsString(conv::kPhypMixing23)
    << " don't have a 2D prior" << endl;
    cerr << "There is something wrong!" << endl;
    cerr << endl;
    util::Die(1);
  }
  
#endif
  
  
  
  // Marginalze 23 parameters with 2D priors if correlated
  
#ifdef __MARGINALIZE_23PARAMS_CORR__
  
  if(Nparam_with_2Dprior==0){ // There are not 2D prior parameters
    cerr << endl;
    cerr << "FitMarginalization::SetMargOscParam" << endl;
    cerr << "# of marginalization parameters with 2D prior = 0" << endl;
    cerr << "Nparam_with_2Dprior = " << Nparam_with_2Dprior << endl;
    cerr << "We expect 2 parameters with 2D priors for 23 parameters" << endl;
    cerr << endl;
    util::Die(1);
  }
  else{
    
    double dmsq=conv::kBadResult;
    double mix23=conv::kBadResult;
    infparm->GetPriorHisto2D(conv::kPhypDm23sq)->GetRandom2(mix23,dmsq);
    
#ifdef  __print_oscparam_marginalized__
    cout << "Marginalize " << conv::kPhypMixing23 << " " << util::AsString(conv::kPhypMixing23) << " --> " << mix23 << endl;
    cout << "Marginalize " << conv::kPhypDm23sq << " " << util::AsString(conv::kPhypDm23sq) << " --> " << dmsq << endl;
    //cout << "Nparam_with_2Dprior = " << Nparam_with_2Dprior << endl;
#endif
    
    // NEW FHC+RHC JOINT
    //infparm->SetInitValue(conv::kPhypDm23sq,dmsq);
    //infparm->SetInitValue(conv::kPhypMixing23,mix23);
    infparm->Update(conv::kPhypDm23sq,dmsq);
    infparm->Update(conv::kPhypMixing23,mix23);
    //
  }
  
#endif // __MARGINALIZE_23PARAMS_CORR__
  
  
  
  // Marginalze theta13 and deltaCP parameters with 2D priors if correlated
  
#ifdef __MARGINALIZE_TH13_DELTACP_CORR__
  
  if(Nparam_with_2Dprior==0){ // There are not 2D prior parameters
    cerr << endl;
    cerr << "FitMarginalization::SetMargOscParam" << endl;
    cerr << "# of marginalization parameters with 2D prior = 0" << endl;
    cerr << "Nparam_with_2Dprior = " << Nparam_with_2Dprior << endl;
    cerr << "We expect 2 parameters with 2D priors for 23 parameters" << endl;
    cerr << endl;
    util::Die(1);
  }
  
  else{
    
    double dcp=conv::kBadResult;
    double mix13=conv::kBadResult;
    infparm->GetPriorHisto2D(conv::kPhypMixing13)->GetRandom2(mix13,dcp);
    
#ifdef  __print_oscparam_marginalized__
    cout << "Marginalize " << conv::kPhypMixing13 << " " << util::AsString(conv::kPhypMixing13) << " --> " << mix13 << endl;
    cout << "Marginalize " << conv::kPhypCPV << " " << util::AsString(conv::kPhypCPV) << " --> " << dcp << endl;
    //cout << "Nparam_with_2Dprior = " << Nparam_with_2Dprior << endl;
#endif
    
    // NEW FHC+RHC JOINT
    //infparm->SetInitValue(conv::kPhypMixing13,mix13);
    //infparm->SetInitValue(conv::kPhypCPV,dcp);
    infparm->Update(conv::kPhypMixing13,mix13);
    infparm->Update(conv::kPhypCPV,dcp);
    //
  }
  
#endif // __MARGINALIZE_TH13_DELTACP_CORR__
  
}

void FitMarginalization::SetNoOscillation(FitParamList* fparm)
{
   int N = fparm->GetNumberParameters();
   for(int ip = 0; ip < N; ++ip)
   {
      conv::PhysParamId_t p = fparm->WhichParam(ip);
      if(!util::IsSystematic(p))
      {
         if(util::IsMixingAngle(p) || util::IsMixingAngleBar(p) ||
               p == conv::kPhypCPV || p == conv::kPhypCPVBar)
         {
            cout << "Setting oscillation parameter to zero: " << util::AsString(p) << endl;
            fparm->SetDoMarg(p, false);
            fparm->SetInitValue(p, 0.0);
         }
      }
   }
}

//________________________________________________________________________________________
void FitMarginalization::BuildDefaultPriors(void)
{
  double badmin = -conv::kBadResult;
  double badmax =  conv::kBadResult;
  
  double dcp_min = -conv::kpi;
  double dcp_max =  conv::kpi;
  
  if (Opt::Instance()->gOscDeltaCP==kNuDeltaCPSin) {
    dcp_min = -1;
    dcp_max =  1;
  }
  
  for (int ip = 0 ; ip < conv::kNuis000 ; ip++) {
    PhysParamId_t p = (PhysParamId_t)ip;
    
    if (not fparm->Exists(p)) {continue;}
    
    if (util::IsMixingAngle(p) || util::IsMassSplitting(p)) {
      double mean    = Opt::Instance()->GetGlobalBestFit(p);
      double fracerr = Opt::Instance()->GetGlobalFitError(p);
      
      if ( ( fabs(mean -conv::kBadResult) < conv::kAVerySmallNum ) ||
          ( fabs(fracerr-conv::kBadResult) < conv::kAVerySmallNum ) ) {
        this->SetPrior_NH(p,"flat",badmin,badmax);
        this->SetPrior_IH(p,"flat",badmin,badmax);
      }
      else {
        double sigma = fracerr*mean;
        //~ std::cout<<"Setting default prior: "<<util::AsString(p)<<"\tGaus("<<mean<<","<<sigma<<")"<<std::endl;
        this->SetPrior_NH(p,"gaus",mean,sigma);
        this->SetPrior_IH(p,"gaus",mean,sigma);
      }
    }
    else if (util::IsPhase(p)) {
      this->SetPrior_NH(p,"flat",dcp_min,dcp_max);
      this->SetPrior_IH(p,"flat",dcp_min,dcp_max);
    }
    
    else {
      this->SetPrior_NH(p,"flat",badmin,badmax);
      this->SetPrior_IH(p,"flat",badmin,badmax);
    }
  }
  return;
}
//________________________________________________________________________________________
bool FitMarginalization::RethrowOscUnphysical(conv::PhysParamId_t p,double oscval) const
{
  
  // NuBar
  if(p == conv::kPhypMixing13Bar){
    if(oscval < mixing13Bar_min || oscval > mixing13Bar_max) return true;
  }
  else if(p == conv::kPhypMixing23Bar){
    if(oscval < mixing23Bar_min || oscval > mixing23Bar_max) return true;
  }
  else if(p == conv::kPhypMixing12Bar){
    if(oscval < mixing12Bar_min || oscval > mixing12Bar_max) return true;
  }
  else if(p == conv::kPhypDm23sqBar){
    if(oscval < dmsq32Bar_min || oscval > dmsq32Bar_max) return true;
  }
  else if(p == conv::kPhypDm12sqBar){
    if(oscval < dmsq12Bar_min || oscval > dmsq12Bar_max) return true;
  }
  else if(p == conv::kPhypCPVBar){
    return false;
  }
  
  // NEW BETABAR
  else if(p == conv::kPhypBetaBar){
    if(oscval < betaBar_min || oscval > betaBar_max) return true;
  }
  //
  
  // Nu
  else if(p == conv::kPhypMixing13){
    if(oscval < mixing13_min || oscval > mixing13_max) return true;
  }
  else if(p == conv::kPhypMixing23){
    if(oscval < mixing23_min || oscval > mixing23_max) return true;
  }
  else if(p == conv::kPhypMixing12){
    if(oscval < mixing12_min || oscval > mixing12_max) return true;
  }
  else if(p == conv::kPhypDm23sq){
    if(oscval < dmsq32_min || oscval > dmsq32_max) return true;
  }
  else if(p == conv::kPhypDm12sq){
    if(oscval < dmsq12_min || oscval > dmsq12_max) return true;
  }
  else if(p == conv::kPhypCPV){
    return false;
  }
  
  else{
    cerr << endl;
    cerr << "Parameter " << util::AsString(p) << " is not in the physics list" << endl;
    cerr << endl;
    util::Die(1);
  }
  
  return false;
}





// NEW FHC+RHC JOINT
double FitMarginalization::GetMinOscParam(PhysParamId_t physpar) const
{
  if     (physpar == conv::kPhypMixing12Bar) return mixing12Bar_min; // NEW OSC PRIOR
  else if(physpar == conv::kPhypMixing13Bar) return mixing13Bar_min;
  else if(physpar == conv::kPhypMixing23Bar) return mixing23Bar_min;
  else if(physpar == conv::kPhypDm23sqBar)   return dmsq32Bar_min;
  else if(physpar == conv::kPhypDm12sqBar)   return dmsq12Bar_min;
  else if(physpar == conv::kPhypCPVBar)      return deltacpBar_min;
  else if(physpar == conv::kPhypBetaBar)     return betaBar_min;
  
  else if(physpar == conv::kPhypMixing12)    return mixing12_min; // NEW OSC PRIOR
  else if(physpar == conv::kPhypMixing13)    return mixing13_min;
  else if(physpar == conv::kPhypMixing23)    return mixing23_min;
  else if(physpar == conv::kPhypDm23sq)      return dmsq32_min;
  else if(physpar == conv::kPhypDm12sq)      return dmsq12_min;
  else if(physpar == conv::kPhypCPV)         return deltacp_min;
  
  else return conv::kBadResult;
}

double FitMarginalization::GetMaxOscParam(PhysParamId_t physpar) const
{
  if     (physpar == conv::kPhypMixing12Bar) return mixing12Bar_max; // NEW OSC PRIOR
  else if(physpar == conv::kPhypMixing13Bar) return mixing13Bar_max;
  else if(physpar == conv::kPhypMixing23Bar) return mixing23Bar_max;
  else if(physpar == conv::kPhypDm23sqBar)   return dmsq32Bar_max;
  else if(physpar == conv::kPhypDm12sqBar)   return dmsq12Bar_max;
  else if(physpar == conv::kPhypCPVBar)      return deltacpBar_max;
  else if(physpar == conv::kPhypBetaBar)     return betaBar_max;
  
  else if(physpar == conv::kPhypMixing12)    return mixing12_max; // NEW OSC PRIOR
  else if(physpar == conv::kPhypMixing13)    return mixing13_max;
  else if(physpar == conv::kPhypMixing23)    return mixing23_max;
  else if(physpar == conv::kPhypDm23sq)      return dmsq32_max;
  else if(physpar == conv::kPhypDm12sq)      return dmsq12_max;
  else if(physpar == conv::kPhypCPV)         return deltacp_max;
  
  else return conv::kBadResult;
}
//





// NEW MARG MH

//________________________________________________________________________________________
double FitMarginalization::Randomize(PhysParamId_t physpar, FitParamList * infparm) const
{
  double valrndm = conv::kBadResult;
  
  if(infparm->GetPriorModeName(physpar)=="flat"){
    double min = infparm->GetMinValueFlat(physpar);
    double max = infparm->GetMaxValueFlat(physpar);
    valrndm = Rnd::Instance()->RndToyMCosc().Uniform(min,max);
  }
  
  else if(infparm->GetPriorModeName(physpar)=="gaus"){
    double mean = infparm->GetMeanGaus(physpar);
    double sigma = infparm->GetSigmaGaus(physpar);
    valrndm = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
  }
  
  else if(infparm->GetPriorModeName(physpar)=="user"){
    valrndm =  infparm->GetPriorHisto(physpar)->GetRandom();
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::Randomize" << endl;
    cerr << "Used prior is: " << infparm->GetPriorModeName(physpar) << endl;
    cerr << "Prior not well defined for " << util::AsString(physpar) << endl;
    util::Die(1);
  }
  
  
  
  // NEW OSC PRIOR
  
  // Convert from sinsq2th to sinsqth or viceversa is mixing angle
  
  if(util::IsMixingAngle(physpar)){
    //cout << util::AsString(physpar) << endl;
    
    // mixparamchoice = sinsqtheta
    
    if( Opt::Instance()->gOscMixingParamChoice == kNuMixSinsqTheta){
#ifdef __USE_PRIOR_SINSQ2TH23__
      if(physpar == conv::kPhypMixing23 || physpar == conv::kPhypMixing23Bar)
      valrndm = util::FromSinsq2thToSinsqth(valrndm);
#endif
#ifdef __USE_PRIOR_SINSQ2TH12__
      if(physpar == conv::kPhypMixing12 || physpar == conv::kPhypMixing12Bar)
      valrndm = util::FromSinsq2thToSinsqth(valrndm);
#endif
#ifdef __USE_PRIOR_SINSQ2TH13__
      if(physpar == conv::kPhypMixing13 || physpar == conv::kPhypMixing13Bar)
      valrndm = util::FromSinsq2thToSinsqth(valrndm);
#endif
      return valrndm;
    }
    
    // mixparamchoice = sinsq2theta
    
    else if( Opt::Instance()->gOscMixingParamChoice == kNuMixSinsq2Theta){
#ifndef __USE_PRIOR_SINSQ2TH23__
      if(physpar == conv::kPhypMixing23 || physpar == conv::kPhypMixing23Bar)
      valrndm = util::FromSinsqthToSinsq2th(valrndm);
#endif
#ifndef __USE_PRIOR_SINSQ2TH12__
      if(physpar == conv::kPhypMixing12 || physpar == conv::kPhypMixing12Bar)
      valrndm = util::FromSinsqthToSinsq2th(valrndm);
#endif
#ifndef __USE_PRIOR_SINSQ2TH13__
      if(physpar == conv::kPhypMixing13 || physpar == conv::kPhypMixing13Bar)
      valrndm = util::FromSinsqthToSinsq2th(valrndm);
#endif
      return valrndm;
    }
    
    // mixparamchoice not defined
    
    else {
      cerr << endl;
      cerr << "FitMarginalization::Randomize" << endl;
      cerr << "Mixing param choice is not defined" << endl;
      cerr << endl;
      util::Die(1);
    }
    
  } // IsMixingAngle
  
  //
  
  return valrndm;
}
//________________________________________________________________________________________

// NEW MARG MH 2

int FitMarginalization::RandomizeMH(const double probNH) const
{
  
  const double min = 0.;
  const double max = 1.;
  double Pnh = probNH;
  double Pih = 1 - probNH;
  
  if(Pnh < 0. || Pnh > 1.){
    cerr << endl;
    cerr << "FitMarginalization::RandomizeMH" << endl;
    cerr << "The input probability for NH is not possible" << endl;
    cerr << "P(NH) = " << Pnh << endl;
    cerr << endl;
    util::Die(1);
  }
  if(Pih < 0. || Pih > 1.){
    cerr << endl;
    cerr << "FitMarginalization::RandomizeMH" << endl;
    cerr << "The input probability for IH is not possible" << endl;
    cerr << "P(IH) = " << Pih << endl;
    cerr << endl;
    util::Die(1);
  }
  
  double rndm = Rnd::Instance()->RndToyMCosc().Uniform(min,max);
  
  int MH = conv::kBadResult;
  if(rndm <= Pnh) MH = 0; // Normal Hierarchy
  else            MH = 1; // Inverted Hierarchy
  
  return MH;
}

//
//________________________________________________________________________________________

// NEW MARG MH 2
void FitMarginalization::SetProbNormalHierarchy(double probnh){
  fProbNormalHierarchy = probnh;
  cout << endl;
  cout << "Set prior probability for Normal Hierarchy to " << fProbNormalHierarchy << endl;
  cout << endl;
  
  fhPriorMH->SetBinContent(1,fProbNormalHierarchy);
  fhPriorMH->SetBinContent(2,1.-fProbNormalHierarchy);
}
//

//________________________________________________________________________________________
void FitMarginalization::SetMargMassHierarchy(bool domargMH){
  fdoMargMassHierarchy = domargMH;
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetRethrowSystUnphysical(bool dorethrow){
  fRethrowSystUnphysical = dorethrow;
  return;
}
//________________________________________________________________________________________





void FitMarginalization::MergeInputFiles(){
  
  // Initialize the clone VALORArray only the first time MergeInputFiles() is called
  
  if(fMergeInputFiles_Called){ // Initialize the VALORArray only the first time MergeInputFiles() is called
    
    
    
    
    
    // NEW 2D TEMPLATES (VARRAY)
    
    //     // NEW FHC+RHC JOINT
    //     // Initialize the Clone array to the maximum # of ereco bins among all the sk samples
    // #ifndef __USE_2D_TEMPLATES__
    //     fToyMargVArray_Clone = new VALORArrayF(fntoys_marg, fnbins_param[0],fnbins_param[1], fNSKsamples,fMaxSampleErecoBins); // total range of free parameter
    // #else
    //     // NEW 2D TEMPLATES
    //     fToyMargVArray_Clone = new VALORArrayF(fntoys_marg, fnbins_param[0],fnbins_param[1], fNSKsamples,fMaxSampleErecoBins,fMaxSampleQ2Bins); // total range of free parameter
    //     //
    // #endif
    //     //
    //
    // fToyMargVArray_Clone->SetAll(-1.);
    //
    
    int nbins_ereco_1Re     = util::NumErecoBins(kSKsFCFV1RingELike);
    int nbins_ereco_1Rmu    = util::NumErecoBins(kSKsFCFV1RingMuLike);
    int nbins_ereco_RHC1Re  = util::NumErecoBins(kSKsRHCFCFV1RingELike);
    int nbins_ereco_RHC1Rmu = util::NumErecoBins(kSKsRHCFCFV1RingMuLike);
    int nbins_ereco_MultiRe = util::NumErecoBins(kSKsFCFVMultiRingELike);
    
#ifndef __USE_2D_TEMPLATES__
    
    fToyMargVArray_Clone_1Re     = new VALORArrayF(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_1Re);
    fToyMargVArray_Clone_1Rmu    = new VALORArrayF(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_1Rmu);
    fToyMargVArray_Clone_RHC1Re  = new VALORArrayF(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_RHC1Re);
    fToyMargVArray_Clone_RHC1Rmu = new VALORArrayF(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_RHC1Rmu);
    fToyMargVArray_Clone_MultiRe = new VALORArrayF(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_MultiRe);
    
    fToyMargVArray_Clone_1Re    ->SetAll(-1.);
    fToyMargVArray_Clone_1Rmu   ->SetAll(-1.);
    fToyMargVArray_Clone_RHC1Re ->SetAll(-1.);
    fToyMargVArray_Clone_RHC1Rmu->SetAll(-1.);
    fToyMargVArray_Clone_MultiRe->SetAll(-1.);
#else
    
    int nbins_q2_1Re     = util::NumQ2Bins(kSKsFCFV1RingELike);
    int nbins_q2_1Rmu    = util::NumQ2Bins(kSKsFCFV1RingMuLike);
    int nbins_q2_RHC1Re  = util::NumQ2Bins(kSKsRHCFCFV1RingELike);
    int nbins_q2_RHC1Rmu = util::NumQ2Bins(kSKsRHCFCFV1RingMuLike);
    int nbins_q2_MultiRe = util::NumQ2Bins(kSKsFCFVMultiRingELike);
    
    fToyMargVArray_Clone_1Re     = new VALORArrayF(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_1Re,     nbins_q2_1Re);
    fToyMargVArray_Clone_1Rmu    = new VALORArrayF(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_1Rmu,    nbins_q2_1Rmu);
    fToyMargVArray_Clone_RHC1Re  = new VALORArrayF(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_RHC1Re,  nbins_q2_RHC1Re);
    fToyMargVArray_Clone_RHC1Rmu = new VALORArrayF(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_RHC1Rmu, nbins_q2_RHC1Rmu);
    fToyMargVArray_Clone_MultiRe = new VALORArrayF(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nbins_ereco_MultiRe, nbins_q2_MultiRe);
    
    fToyMargVArray_Clone_1Re     ->SetAll(-1.);
    fToyMargVArray_Clone_1Rmu    ->SetAll(-1.);
    fToyMargVArray_Clone_RHC1Re  ->SetAll(-1.);
    fToyMargVArray_Clone_RHC1Rmu ->SetAll(-1.);
    fToyMargVArray_Clone_MultiRe ->SetAll(-1.);
#endif
    
    //
    
    
    
    
    fMergeInputFiles_Called = false;
  }
  
  cout << endl;
  cout << "Merge toys of input file:" << endl
  << "bin: ["
  << fmin_bin_sub_param[0]
  << ","
  << fmax_bin_sub_param[0]
  << "]"
  << endl
  << fnamefreeparam[0] << ": "
  << "["
  << fhchisq_marg_2D->GetXaxis()->GetBinCenter(fmin_bin_sub_param[0] + 1)
  << ","
  << fhchisq_marg_2D->GetXaxis()->GetBinCenter(fmax_bin_sub_param[0] + 1)
  << "]"
  << endl;
  cout << "bin: ["
  << fmin_bin_sub_param[1]
  << ","
  << fmax_bin_sub_param[1]
  << "]"
  << endl
  << fnamefreeparam[1] << ": "
  << "["
  << fhchisq_marg_2D->GetYaxis()->GetBinCenter(fmin_bin_sub_param[1] + 1)
  << ","
  << fhchisq_marg_2D->GetYaxis()->GetBinCenter(fmax_bin_sub_param[1] + 1)
  << "]"
  << endl;
  
  // Copy read array elements to the clone one. Same loop as in FillChisqMarg()
  for(int iexperiment = 0; iexperiment < fntoys_marg; iexperiment++) {
    
    for(int ipar0=fmin_bin_sub_param[0];ipar0<=fmax_bin_sub_param[0];ipar0++){
      for(int ipar1=fmin_bin_sub_param[1];ipar1<=fmax_bin_sub_param[1];ipar1++){
        
        // Loop over the used SK samples
        int iscount=0;
        //for(int is = 0; is < kns; is++) {
        //SuperkSample_t s = kSKSamples[is];
        //if(s == kSKsUndefined) continue;
        for(int is=0;is<fNSKsamples;is++){
          SuperkSample_t s = this->GetSKSamples(is);
          
          
#ifndef __USE_2D_TEMPLATES__
          
          for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){
            
            // NEW 2D TEMPLATES (VARRAY)
            //float bincont = fToyMargVArray->GetEntry(iexperiment,ipar0-fmin_bin_sub_param[0],ipar1-fmin_bin_sub_param[1],iscount,ibin); // shift to sub-range // bincont = -1 by default
            float bincont = this->GetToyMargVArraySample(s)->GetEntry(iexperiment,ipar0-fmin_bin_sub_param[0],ipar1-fmin_bin_sub_param[1],ibin); // shift to sub-range // bincont = -1 by default
                                                                                                                                                 //
            
            if(bincont >= 0.){
              
              // NEW 2D TEMPLATES (VARRAY)
              //float bincont_clone = fToyMargVArray_Clone->GetEntry(iexperiment,ipar0,ipar1,iscount,ibin);  // total range
              float bincont_clone = this->GetToyMargVArraySample_Clone(s)->GetEntry(iexperiment,ipar0,ipar1,ibin);  // total range
                                                                                                                    //
              
              if(bincont_clone >= 0.){
                cerr << endl;
                cerr << "FitMarginalization::MergeInputFiles" << endl;
                cerr << "Following bin has been already filled:" << endl;
                cerr << " - iexperiment = " << iexperiment << endl;
                cerr << " - ipar0 = " << ipar0 << endl;
                cerr << " - ipar1 = " << ipar1 << endl;
                cerr << " - sample = " << util::AsStringCompact(s) << endl;
                cerr << " - iereco = " << ibin << endl;
                util::Die(1);
              }
              else{
                
                // NEW 2D TEMPLATES (VARRAY)
                //fToyMargVArray_Clone->SetEntry(bincont,iexperiment,ipar0,ipar1,iscount,ibin); // Set to 0 if negative
                this->GetToyMargVArraySample_Clone(s)->SetEntry(bincont,iexperiment,ipar0,ipar1,ibin); // Set to 0 if negative
                                                                                                       //
                
              }
            }
          } // loop ibin
#else
          // NEW 2D TEMPLATES
          for(int ibinEreco=0;ibinEreco<util::NumErecoBins(s);ibinEreco++){
            for(int ibinQ2=0;ibinQ2<util::NumQ2Bins(s);ibinQ2++){
              
              
              // NEW 2D TEMPLATES (VARRAY)
              //float bincont = fToyMargVArray->GetEntry(iexperiment,ipar0-fmin_bin_sub_param[0],ipar1-fmin_bin_sub_param[1],iscount,ibinEreco,ibinQ2); // shift to sub-range // bincont = -1 by default
              float bincont = this->GetToyMargVArraySample(s)->GetEntry(iexperiment,ipar0-fmin_bin_sub_param[0],ipar1-fmin_bin_sub_param[1],ibinEreco,ibinQ2); // shift to sub-range // bincont = -1 by default
                                                                                                                                                               //
              
              
              if(bincont >= 0.){
                
                // NEW 2D TEMPLATES (VARRAY)
                //float bincont_clone = fToyMargVArray_Clone->GetEntry(iexperiment,ipar0,ipar1,iscount,ibinEreco,ibinQ2);  // total range
                float bincont_clone = this->GetToyMargVArraySample_Clone(s)->GetEntry(iexperiment,ipar0,ipar1,ibinEreco,ibinQ2);  // total range
                                                                                                                                  //
                
                
                if(bincont_clone >= 0.){
                  cerr << endl;
                  cerr << "FitMarginalization::MergeInputFiles" << endl;
                  cerr << "Following bin has been already filled:" << endl;
                  cerr << " - iexperiment = " << iexperiment << endl;
                  cerr << " - ipar0 = " << ipar0 << endl;
                  cerr << " - ipar1 = " << ipar1 << endl;
                  cerr << " - sample = " << util::AsStringCompact(s) << endl;
                  cerr << " - iereco = " << ibinEreco << endl;
                  cerr << " - iq2 = " << ibinQ2 << endl;
                  cerr << " - bincont = " << bincont << endl;
                  cerr << " - bincont clone = " << bincont_clone << endl;
                  util::Die(1);
                }
                else{
                  
                  // NEW 2D TEMPLATES (VARRAY)
                  //fToyMargVArray_Clone->SetEntry(bincont,iexperiment,ipar0,ipar1,iscount,ibinEreco,ibinQ2); // Set to 0 if negative
                  this->GetToyMargVArraySample_Clone(s)->SetEntry(bincont,iexperiment,ipar0,ipar1,ibinEreco,ibinQ2); // Set to 0 if negative
                                                                                                                     //
                  
                }
              }
            } // loop ibinQ2
          } // loop ibinEreco
            //
          
#endif // __USE_2D_TEMPLATES__
          
          iscount++;
        } // loop SK samples
        
        
#ifdef __print_array_toy_merge__
        cout << endl;
        cout << "toy=" << iexperiment << " - ipar0 = " << ipar0 << " - ipar1 = " << ipar1 << " : ";
        
        // Loop over the used SK samples
        iscount=0;
        //for(int is = 0; is < kns; is++) {
        //SuperkSample_t s = kSKSamples[is];
        //if(s == kSKsUndefined) continue;
        for(int is=0;is<fNSKsamples;is++){
          SuperkSample_t s = this->GetSKSamples(is);
          
          cout << util::AsStringCompact(s) << endl;
          
#ifndef __USE_2D_TEMPLATES__
          for(int ib = 0; ib < util::NumErecoBins(s); ib++){
            
            // NEW 2D TEMPLATES (VARRAY)
            //double bincont = fToyMargVArray_Clone->GetEntry(iexperiment,ipar0,ipar1,iscount,ib);
            float bincont = this->GetToyMargVArraySample_Clone(s)->GetEntry(iexperiment,ipar0,ipar1,ib);
            //
            
            cout << bincont << " ";
          }
#else
          // NEW 2D TEMPLATES
          for(int ibEreco = 0; ibEreco < util::NumErecoBins(s); ibEreco++){
            for(int ibQ2 = 0; ibQ2 < util::NumQ2Bins(s); ibQ2++){
              
              // NEW 2D TEMPLATES (VARRAY)
              //double bincont = fToyMargVArray_Clone->GetEntry(iexperiment,ipar0,ipar1,iscount,ibEreco,ibQ2);
              float bincont = this->GetToyMargVArraySample_Clone(s)->GetEntry(iexperiment,ipar0,ipar1,ibEreco,ibQ2);
              //
              
              cout << bincont << " ";
            }
          }
#endif
          cout << endl;
          cout << endl;
          
          iscount++;
        } // loop SK samples
#endif
        
      } // loop ipar1
    } // loop ipar0
  } // loop iexperiment
  
  fNmergedFiles++;
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetNSKsamples(){
  
  fNSKsamples = 0; // # of sk samples used in the fit
  fMaxSampleErecoBins = 1;
  // NEW 2D TEMPLATES
  fMaxSampleQ2Bins = 1;
  //
  cout << endl;
  cout << "Use the following SK samples:" << endl;
  
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    cout << "    " << util::AsStringCompact(s) << endl;
    fNSKsamples++;
  }
  cout << "Total # of samples = " << fNSKsamples << endl;
  cout << endl;
  
  // Take the Ereco binning for each SK samples
  this->fnerecobins = new VALORArrayI(fNSKsamples);
  // NEW 2D TEMPLATES
  this->fnq2bins = new VALORArrayI(fNSKsamples);
  //
  this->fsksamples = new VALORArrayI(fNSKsamples);
  
  int iscount=0;
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    int nbins_ereco = util::NumErecoBins(s);
    this->fnerecobins->SetEntry(nbins_ereco,iscount);
    this->fsksamples->SetEntry(is,iscount);
    if (nbins_ereco > fMaxSampleErecoBins) { fMaxSampleErecoBins = nbins_ereco; }
    
#ifdef __USE_2D_TEMPLATES__
    // NEW 2D TEMPLATES
    int nbins_q2    = util::NumQ2Bins(s);
    this->fnq2bins->SetEntry(nbins_q2,iscount);
    if (nbins_q2 > fMaxSampleQ2Bins) { fMaxSampleQ2Bins = nbins_q2; }
    //
#endif
    
    
    //cout << "is " << is << " - iscount " << iscount << endl;
    //cout << util::AsString(s) << " --> # bins ereco = " << fnerecobins->GetEntry(iscount) << endl;
#ifdef __USE_2D_TEMPLATES__
    cout << util::AsString(s) << " --> # bins q2 = "    << fnq2bins->GetEntry(iscount) << endl;
#endif
    
    iscount++;
  }
  
  return;
}
//________________________________________________________________________________________
conv::SuperkSample_t FitMarginalization::GetSKSamples(int isksample) const
{
  int is = this->fsksamples->GetEntry(isksample);
  SuperkSample_t s = kSKSamples[is];
  return s;
}
//________________________________________________________________________________________



// NEW FHC+RHC JOINT
bool FitMarginalization::UseSKSample(int isksample) const
{
  SuperkSample_t s = this->GetSKSamples(isksample);
  
#ifdef __USE_1RMU_IN_CHISQ__
  if(s==conv::kSKsFCFV1RingMuLike) return true;
#endif
  
#ifdef __USE_1RE_IN_CHISQ__
  if(s==conv::kSKsFCFV1RingELike) return true;
#endif
  
#ifdef __USE_RHC1RMU_IN_CHISQ__
  if(s==conv::kSKsRHCFCFV1RingMuLike) return true;
#endif
  
#ifdef __USE_RHC1RE_IN_CHISQ__
  if(s==conv::kSKsRHCFCFV1RingELike) return true;
#endif
  
#ifdef __USE_MULTIRE_IN_CHISQ__
  if(s==conv::kSKsFCFVMultiRingELike) return true;
#endif
  
  
  return false;
}
//________________________________________________________________________________________
bool FitMarginalization::UseSKSample(SuperkSample_t s ) const
{
#ifdef __USE_1RMU_IN_CHISQ__
  if(s==conv::kSKsFCFV1RingMuLike) return true;
#endif
  
#ifdef __USE_1RE_IN_CHISQ__
  if(s==conv::kSKsFCFV1RingELike) return true;
#endif
  
#ifdef __USE_RHC1RMU_IN_CHISQ__
  if(s==conv::kSKsRHCFCFV1RingMuLike) return true;
#endif
  
#ifdef __USE_RHC1RE_IN_CHISQ__
  if(s==conv::kSKsRHCFCFV1RingELike) return true;
#endif
  
#ifdef __USE_MULTIRE_IN_CHISQ__
  if(s==conv::kSKsFCFVMultiRingELike) return true;
#endif
  
  return false;
}
//


//________________________________________________________________________________________
void FitMarginalization::SetSaveParamLists(bool should_save_lists)
{
  if (fSavedToyParamLists) {
    for (int i = 0 ; i < fntoys_marg ; i++) {
      if (fSavedToyParamLists[i]) { delete fSavedToyParamLists[i];}
    }
    delete fSavedToyParamLists;
    fSavedToyParamLists = 0;
  }
  fShouldSaveToyParamLists = should_save_lists;
  if (should_save_lists) {
    fSavedToyParamLists = (PhysParamList **) calloc(fntoys_marg,sizeof(PhysParamList*));
  }
}
//________________________________________________________________________________________
PhysParamList *const* FitMarginalization::GetParamListArray(void) const
{
  return fSavedToyParamLists;
}
//________________________________________________________________________________________
TObjArray * FitMarginalization::GetToyParamThrowHistograms(int nbins) const
{
  if (fShouldSaveToyParamLists) {
    TObjArray * o = util::MakeHistogramsOfParameterValues(fntoys_marg,fSavedToyParamLists,nbins,fparm);
    for (int i = 0 ; i < o->GetEntries();i++) {
      TH1D * h = dynamic_cast<TH1D*>(o->At(i));
      if (not h) continue;
      const char * name = h->GetName();
      string pname(util::ReplaceAll(name,"param_values_",""));
      PhysParamId_t p = util::PhysParamIdFromRawString(pname);
      if (p==conv::kPhypUndefined) continue;
      h->SetName (Form("param_values_%s",util::ParamAsStringCompact(p).c_str()));
      h->SetTitle(Form("%s Values",util::ReplaceAll(util::ParamAsString(p),";","").c_str()));
    }
    return o;
  }
  else {
    return 0;
  }
}
//________________________________________________________________________________________
void FitMarginalization::SetMinimIntervalRead(double XvalFirst,double XvalLast,double YvalFirst,double YvalLast,bool readbinary){
  
  fFirstMinimVal_X = XvalFirst;
  fLastMinimVal_X = XvalLast;
  fFirstMinimVal_Y = YvalFirst;
  fLastMinimVal_Y = YvalLast;
  
  if(fIsFit2D){
    // fFirstMinimBin_X = fhchisq_marg_2D->GetXaxis()->FindBin(fFirstMinimVal_X+conv::kVeryASmallNum);
    // fLastMinimBin_X  = fhchisq_marg_2D->GetXaxis()->FindBin(fLastMinimVal_X+conv::kAVerySmallNum);
    // fFirstMinimBin_Y = fhchisq_marg_2D->GetYaxis()->FindBin(fFirstMinimVal_Y+conv::kAVerySmallNum);
    // fLastMinimBin_Y  = fhchisq_marg_2D->GetYaxis()->FindBin(fLastMinimVal_Y+conv::kAVerySmallNum);
    fFirstMinimBin_X = fhchisq_marg_2D->GetXaxis()->FindBin(fFirstMinimVal_X);
    fLastMinimBin_X  = fhchisq_marg_2D->GetXaxis()->FindBin(fLastMinimVal_X);
    fFirstMinimBin_Y = fhchisq_marg_2D->GetYaxis()->FindBin(fFirstMinimVal_Y);
    fLastMinimBin_Y  = fhchisq_marg_2D->GetYaxis()->FindBin(fLastMinimVal_Y);
  }
  else{
    // fFirstMinimBin_X = fhchisq_marg_1D->FindBin(fFirstMinimVal_X+conv::kAVerySmallNum);
    // fLastMinimBin_X  = fhchisq_marg_1D->FindBin(fLastMinimVal_X+conv::kAVerySmallNum);
    fFirstMinimBin_X = fhchisq_marg_1D->FindBin(fFirstMinimVal_X);
    fLastMinimBin_X  = fhchisq_marg_1D->FindBin(fLastMinimVal_X);
    fFirstMinimBin_Y = conv::kBadResult;
    fLastMinimBin_Y  = conv::kBadResult;
  }
  
  int nbins_X = conv::kBadResult;
  int nbins_Y = conv::kBadResult;
  nbins_X= fnbins_param[0];
  if(fIsFit2D)
  nbins_Y= fnbins_param[1];
  
  // If read binary file at each grid point step, the range must be contained in full range
  // while if binary files are first stored and then the chisq is computed, since a binary
  // file is read and then the corresponding part of chisq is filled, the range must
  // be consistent with the read binary file. The same is when minimization is done without
  // without reading binary at each step: the minimization range must be contained in the
  // binary file one.
  
  if(!readbinary){
    
    // Don't exceed the allowed bins (N.B. Use histogram bins from 1 to NBins!!!)
    if(fFirstMinimBin_X < fmin_bin_sub_param[0]+1){
      cerr << "The first bin of the X interval is outside the sub-range limit" << endl;
      cerr << "Set it to " << fmin_bin_sub_param[0]+1 << endl;
      cerr << endl;
      fFirstMinimBin_X = fmin_bin_sub_param[0]+1;
      if(fIsFit2D) fFirstMinimVal_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(fFirstMinimBin_X);
      else         fFirstMinimVal_X = fhchisq_marg_1D->GetBinCenter(fFirstMinimBin_X);
    }
    if(fLastMinimBin_X > fmax_bin_sub_param[0]+1){
      cerr << endl;
      cerr << "The last bin of the X interval is outside the sub-range limit" << endl;
      cerr << "Set it to " << fmax_bin_sub_param[0]+1 << endl;
      cerr << endl;
      fLastMinimBin_X = fmax_bin_sub_param[0]+1;
      if(fIsFit2D) fLastMinimVal_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(fLastMinimBin_X);
      else         fLastMinimVal_X = fhchisq_marg_1D->GetBinCenter(fLastMinimBin_X);
    }
    if(fIsFit2D){
      if(fFirstMinimBin_Y < fmin_bin_sub_param[1]+1){
        cerr << endl;
        cerr << "The first bin of the Y interval is outside the sub-range limit" << endl;
        cerr << "Set it to " << fmin_bin_sub_param[1]+1 << endl;
        cerr << endl;
        fFirstMinimBin_Y = fmin_bin_sub_param[1]+1;
        if(fIsFit2D) fFirstMinimVal_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(fFirstMinimBin_Y);
        else         fFirstMinimVal_Y = fhchisq_marg_1D->GetBinCenter(fFirstMinimBin_Y);
      }
      if(fLastMinimBin_Y > fmax_bin_sub_param[1]+1){
        cerr << endl;
        cerr << "The last bin of the Y interval is outside the sub-range limit" << endl;
        cerr << "Set it to " << fmax_bin_sub_param[1]+1 << endl;
        cerr << endl;
        fLastMinimBin_Y = fmax_bin_sub_param[1]+1;
        if(fIsFit2D) fLastMinimVal_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(fLastMinimBin_Y);
        else         fLastMinimVal_Y = fhchisq_marg_1D->GetBinCenter(fLastMinimBin_Y);
      }
    }
  }
  
  else{
    // Don't exceed the allowed bins
    if(fFirstMinimBin_X < 1){
      cerr << endl;
      cerr << "The first bin of the X interval is outside the histogram limit" << endl;
      cerr << "Set it to 1" << endl;
      cerr << endl;
      fFirstMinimBin_X = 1;
      fFirstMinimVal_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(fFirstMinimBin_X);
    }
    if(fLastMinimBin_X > nbins_X){
      cerr << endl;
      cerr << "The last bin of the X interval is outside the histogram limit" << endl;
      cerr << "Set it to the number of bins in the histogram:" << nbins_X << endl;
      cerr << endl;
      fLastMinimBin_X = nbins_X;
      fLastMinimVal_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(fLastMinimBin_X);
    }
    if(fIsFit2D){
      if(fFirstMinimBin_Y < 1){
        cerr << endl;
        cerr << "The first bin of the Y interval is outside the histogram limit" << endl;
        cerr << "Set it to 1" << endl;
        cerr << endl;
        fFirstMinimBin_Y = 1;
        fFirstMinimVal_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(fFirstMinimBin_Y);
      }
      if(fLastMinimBin_Y > nbins_Y){
        cerr << endl;
        cerr << "The last bin of the Y interval is outside the histogram limit" << endl;
        cerr << "Set it to the number of bins in the histogram:" << nbins_Y << endl;
        cerr << endl;
        fLastMinimBin_Y = nbins_Y;
        fLastMinimVal_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(fLastMinimBin_Y);
      }
    }
  }
  
  cout << endl;
  cout << "Minimization interval set to:" << endl;
  cout << "  - X: [" << fFirstMinimVal_X << "," << fLastMinimVal_X << "]" << endl;
  cout << "  - X Bin: [" << fFirstMinimBin_X << "," << fLastMinimBin_X << "]" << endl;
  cout << "  - Y: [" << fFirstMinimVal_Y << "," << fLastMinimVal_Y << "]" << endl;
  cout << "  - Y Bin: [" << fFirstMinimBin_Y << "," << fLastMinimBin_Y << "]" << endl;
  cout << endl;
  
  return;
}

//________________________________________________________________________________________
void FitMarginalization::SetMinimUncorrXY(bool flaguncorr){
  fUncorrXY = flaguncorr;
}
//________________________________________________________________________________________
void FitMarginalization::MinimizeChisqIntervalRead(SuperkObservation *skobs_data,bool readbinary){
  
  double bigvalueUOflow = 99999999.;  // NEW PRECISION
  
  // 1D
  // fhchisq_marg_1D_Step = (TH1D*)fhchisq_marg_1D->Clone("hchisq_marg_1D_Step"); // NEW PRECISION
  // fhchisq_marg_1D_Step->Reset();
  // fhchisq_marg_1D_Step->SetName("hchisq_marg_1D_Step");
  // fhchisq_marg_1D_Step->SetTitle("hchisq_marg_1D_Step");
  fhchisq_marg_1D_Step->Reset();
  for(int ibinX=0;ibinX<fhchisq_marg_1D_Step->GetNbinsX();ibinX++){
    fhchisq_marg_1D_Step->SetBinContent(ibinX+1,0.);
  }
  // Set under and over flow of X to very high value
  fhchisq_marg_1D->SetBinContent(0,bigvalueUOflow);
  fhchisq_marg_1D->SetBinContent(fhchisq_marg_1D->GetNbinsX()+1,bigvalueUOflow); // FIX
  
  
  // 2D
  // fhchisq_marg_2D_Step = (TH2D*)fhchisq_marg_2D->Clone("hchisq_marg_2D_Step"); // NEW PRECISION
  // fhchisq_marg_2D_Step->Reset();
  // fhchisq_marg_2D_Step->SetName("hchisq_marg_2D_Step");
  // fhchisq_marg_2D_Step->SetTitle("hchisq_marg_2D_Step");
  fhchisq_marg_2D_Step->Reset();
  for(int ibinX=0;ibinX<fhchisq_marg_2D_Step->GetNbinsX();ibinX++){
    for(int ibinY=0;ibinY<fhchisq_marg_2D_Step->GetNbinsY();ibinY++){
      fhchisq_marg_2D_Step->SetBinContent(ibinX+1,ibinY+1,0.);
    }
  }
  // Set under and over flow of X and Y to very high value
  for(int ibinX=0;ibinX<fhchisq_marg_2D->GetNbinsX();ibinX++){
    fhchisq_marg_2D->SetBinContent(ibinX+1,0,bigvalueUOflow);
    fhchisq_marg_2D->SetBinContent(ibinX+1,fhchisq_marg_2D->GetNbinsY()+1,bigvalueUOflow);
  }
  for(int ibinY=0;ibinY<fhchisq_marg_2D->GetNbinsY();ibinY++){
    fhchisq_marg_2D->SetBinContent(0,ibinY+1,bigvalueUOflow);
    fhchisq_marg_2D->SetBinContent(fhchisq_marg_2D->GetNbinsX()+1,ibinY+1,bigvalueUOflow);
  }
  
  
  TRandom3 *r = new TRandom3(0); // needed to randomly scan the distribution
  
  // Get the interval
  
  double UniFirstBin_X = fFirstMinimBin_X;
  double UniLastBin_X  = fLastMinimBin_X;
  double UniFirstBin_Y = fFirstMinimBin_Y;
  double UniLastBin_Y  = fLastMinimBin_Y;
  
  double binF_X = UniFirstBin_X; // first value of the free parameter
  double binL_X = UniLastBin_X; // last value of the free parameter
  double binF_Y = UniFirstBin_Y; // first value of the free parameter
  double binL_Y = UniLastBin_Y; // last value of the free parameter
  
  double bincont = 0.;
  
  double chisq_X[3] = {0., 0., 0.};  // Chisq(N-1), Chisq(N) and Chisq(N+1)
  double chisq_Y[3] = {0., 0., 0.};  // Chisq(N-1), Chisq(N) and Chisq(N+1)
  
  double chisq[3][3] = {{0.}};  // Chisq(N-1), Chisq(N) and Chisq(N+1)
  
  const int ilike_N = 1;
  const int ilike_Nm1 = 0;
  const int ilike_Np1 = 2;
  
  int Niterations=0;
  int Niterations_Part=0;
  int Nminima=0;
  
#if defined __print_minimization_steps__
  cout << "Extract Uniform [ " << binF_X << " , " << binL_X << " ] " << endl;
  if(fIsFit2D)
  cout << "Extract Uniform [ " << binF_Y << " , " << binL_Y << " ] " << endl;
#endif
  
  // Throw the bin number
  double r_X = r->Uniform(binF_X,binL_X+1); // +1 because it's then truncated into a int. Extract 21th bin as well
  double r_Y = r->Uniform(binF_Y,binL_Y+1); // +1 because it's then truncated into a int. Extract 21th bin as well
  int bin_X = (int)r_X;
  int bin_Y = (int)r_Y;
  
  bool doNewFit = true;
  
  int Nfits=1; // number of fits performed
  
  
  while(1){
    
    if(Nfits > fMaxFitNumber){
      cout << endl;
      cout << "Starting fit # " << Nfits << endl;
      cout << "Maximum allowed # of fits: " << fMaxFitNumber << endl;
      cout << "Convergence has not found" << endl;
      cout << "# of found minima: " << Nminima << endl;
      cout << endl;
      return;
    }
    else if(doNewFit){
      cout << endl;
      cout << "/////////////// Fit " << Nfits << " ///////////////" << endl;
      cout << endl;
      doNewFit = false;
    }
    
    
    bin_X = r->Uniform(binF_X,binL_X+1); // add 1 because it will be truncated to get the bin number
    if(fIsFit2D)
    bin_Y = r->Uniform(binF_Y,binL_Y+1); // add 1 because it will be truncated to get the bin number
    
    int binN_X = conv::kBadResult;
    int binN_Y = conv::kBadResult;
    
    // same point
    binN_X = bin_X;
    if(fIsFit2D) binN_Y = bin_Y;
    
    double binlowedge_X = conv::kBadResult;
    double binlowedge_Y = conv::kBadResult;
    
    // Find direction along X axis
    
    // L(N,N)
    if(fIsFit2D)
    bincont = fChisqMargVArray->GetEntry(binN_X-1,binN_Y-1);
    else
    bincont = fChisqMargVArray->GetEntry(binN_X-1,fmin_bin_sub_param[1]);
    if( bincont<=0.){
      if(fIsFit2D) bincont = CalcChisqMargPtRead(skobs_data,binN_X-1,binN_Y-1,readbinary); // CalcChisqMargPt takes the point not bin
      else bincont = CalcChisqMargPtRead(skobs_data,binN_X-1,fmin_bin_sub_param[1],readbinary);
    }
    chisq_X[ilike_N] = bincont;
    chisq[ilike_N][ilike_N] = bincont;
    
    // L(N-1,N)
    int binNm1_X = binN_X - 1;
    if(binNm1_X < fFirstMinimBin_X)
    binNm1_X = binN_X;  // should happens only when binN = binF
    if(fIsFit2D)
    bincont = fChisqMargVArray->GetEntry(binNm1_X-1,binN_Y-1);
    else
    bincont = fChisqMargVArray->GetEntry(binNm1_X-1,fmin_bin_sub_param[1]);
    if( bincont<=0.){
      if(fIsFit2D) bincont = CalcChisqMargPtRead(skobs_data,binNm1_X-1,binN_Y-1,readbinary);
      else bincont = CalcChisqMargPtRead(skobs_data,binNm1_X-1,fmin_bin_sub_param[1],readbinary);
    }
    chisq_X[ilike_Nm1] = bincont;
    chisq[ilike_Nm1][ilike_N] = bincont;
    
    // L(N+1,N)
    int binNp1_X = binN_X + 1;
    if(binNp1_X > fLastMinimBin_X)
    binNp1_X = binN_X;  // should happens only when binN = binL
    if(fIsFit2D)
    bincont = fChisqMargVArray->GetEntry(binNp1_X-1,binN_Y-1);
    else
    bincont = fChisqMargVArray->GetEntry(binNp1_X-1,fmin_bin_sub_param[1]);
    if( bincont<=0.){
      if(fIsFit2D) bincont = CalcChisqMargPtRead(skobs_data,binNp1_X-1,binN_Y-1,readbinary);
      else bincont = CalcChisqMargPtRead(skobs_data,binNp1_X-1,fmin_bin_sub_param[1],readbinary); // bin_y = 0
    }
    chisq_X[ilike_Np1] = bincont;
    chisq[ilike_Np1][ilike_N] = bincont;
    
    binlowedge_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(bin_X);
    
    
    // Find direction along Y axis
    
    int binNm1_Y = conv::kBadResult;
    int binNp1_Y = conv::kBadResult;
    
    if(fIsFit2D){
      // L(N,N) already computed for 2D chisq
      bincont = fChisqMargVArray->GetEntry(binN_X-1,binN_Y-1);
      if( bincont<=0.)
      bincont = CalcChisqMargPtRead(skobs_data,binN_X-1,binN_Y-1,readbinary);
      chisq_Y[ilike_N] = bincont;
      
      // L(N,N-1)
      binNm1_Y = binN_Y - 1;
      if(binNm1_Y < fFirstMinimBin_Y) binNm1_Y = binN_Y;  // should happens only when binN = binF
      bincont = fChisqMargVArray->GetEntry(binN_X-1,binNm1_Y-1);
      if( bincont<=0.)
      bincont = CalcChisqMargPtRead(skobs_data,binN_X-1,binNm1_Y-1,readbinary);
      chisq_Y[ilike_Nm1] = bincont;
      chisq[ilike_N][ilike_Nm1] = bincont;
      
      // L(N,N+1)
      binNp1_Y = binN_Y + 1;
      if(binNp1_Y > fLastMinimBin_Y) binNp1_Y = binN_Y;  // should happens only when binL = binL
      bincont = fChisqMargVArray->GetEntry(binN_X-1,binNp1_Y-1);
      if( bincont<=0.)
      bincont = CalcChisqMargPtRead(skobs_data,binN_X-1,binNp1_Y-1,readbinary);
      chisq_Y[ilike_Np1] = bincont;
      chisq[ilike_N][ilike_Np1] = bincont;
      
      binlowedge_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(bin_Y);
    }
    
    // Find off-diagonal directions
    
    if(fIsFit2D){
      
      // L(N-1,N-1)
      bincont = fChisqMargVArray->GetEntry(binNm1_X-1,binNm1_Y-1);
      if( bincont<=0.)
      bincont = CalcChisqMargPtRead(skobs_data,binNm1_X-1,binNm1_Y-1,readbinary);
      chisq[ilike_Nm1][ilike_Nm1] = bincont;
      
      // L(N+1,N+1)
      bincont = fChisqMargVArray->GetEntry(binNp1_X-1,binNp1_Y-1);
      if( bincont<=0.)
      bincont = CalcChisqMargPtRead(skobs_data,binNp1_X-1,binNp1_Y-1,readbinary);
      chisq[ilike_Np1][ilike_Np1] = bincont;
      
      // L(N-1,N+1)
      bincont = fChisqMargVArray->GetEntry(binNm1_X-1,binNp1_Y-1);
      if( bincont<=0.)
      bincont = CalcChisqMargPtRead(skobs_data,binNm1_X-1,binNp1_Y-1,readbinary);
      chisq[ilike_Nm1][ilike_Np1] = bincont;
      
      // L(N+1,N-1)
      bincont = fChisqMargVArray->GetEntry(binNp1_X-1,binNm1_Y-1);
      if( bincont<=0.)
      bincont = CalcChisqMargPtRead(skobs_data,binNp1_X-1,binNm1_Y-1,readbinary);
      chisq[ilike_Np1][ilike_Nm1] = bincont;
    }
    
#if defined __print_minimization_steps__
    cout << endl;
    cout << "# iteration = " << Niterations << endl;
    cout << "# iteration current = " << Niterations_Part << endl;
    cout << bin_X << " - " << binlowedge_X << " --> "
    << binNm1_X << " "
    << binN_X << " "
    << binNp1_X
    << " -> chisq_X = "
    << chisq_X[ilike_Nm1] << " "
    << chisq_X[ilike_N] << " "
    << chisq_X[ilike_Np1]
    << endl
    << chisq[ilike_N][ilike_Nm1] << " "
    << chisq[ilike_N][ilike_N] << " "
    << chisq[ilike_N][ilike_Np1] << " "
    << endl;
    if(fIsFit2D){
      cout << bin_Y << " - " << binlowedge_Y << " --> "
      << binNm1_Y << " "
      << binN_Y << " "
      << binNp1_Y
      << " -> chisq_Y = "
      << chisq_Y[ilike_Nm1] << " "
      << chisq_Y[ilike_N] << " "
      << chisq_Y[ilike_Np1]
      << endl
      << chisq[ilike_Nm1][ilike_N] << " "
      << chisq[ilike_N][ilike_N] << " "
      << chisq[ilike_Np1][ilike_N] << " "
      << endl;
    }
#endif
    
    double bincent_X = conv::kBadResult;
    double bincent_Y = conv::kBadResult;
    
    if(!fIsFit2D){
      bincent_X = fhchisq_marg_1D_Step->GetXaxis()->GetBinCenter(bin_X);
      fhchisq_marg_1D_Step->Fill(bincent_X);
    }
    else{
      bincent_X = fhchisq_marg_2D_Step->GetXaxis()->GetBinCenter(bin_X);
      bincent_Y = fhchisq_marg_2D_Step->GetYaxis()->GetBinCenter(bin_Y);
      fhchisq_marg_2D_Step->Fill(bincent_X,bincent_Y);
    }
    
    
    // Best fit bin is found
    // this condition contains the case with the boundary, when
    // z(N-1) = Z(N) or Z(N+1) = Z(N)
    
    if(
       chisq_X[ilike_N] <= chisq_X[ilike_Np1] &&
       chisq_X[ilike_N] <= chisq_X[ilike_Nm1] &&
       ( !fIsFit2D ||
        (
         chisq[ilike_N][ilike_N] <= chisq[ilike_N][ilike_Np1]   &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_N][ilike_Nm1]   &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_Nm1][ilike_N]   &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_Np1][ilike_N]   &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_Np1][ilike_Np1] &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_Nm1][ilike_Nm1] &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_Nm1][ilike_Np1] &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_Np1][ilike_Nm1]
         )
        )
       ){
      
      if( !fIsFit2D || (chisq_X[ilike_N] == chisq_Y[ilike_N]) ){
        
        // Check if it's a real minimum
        
        cout << "A minimum is found" << endl;
        cout << "Best fit Chisq = " << chisq_X[ilike_N] << endl;
        cout << "Best fit X = " << fhchisq_marg_2D->GetXaxis()->GetBinCenter(binN_X) << endl;
        cout << "Best fit Y = " << fhchisq_marg_2D->GetYaxis()->GetBinCenter(binN_Y) << endl;
        cout << "Total # of iterations = " << Niterations << endl;
        cout << "Total # of performed fits = " << Nfits << endl;
        cout << "Check if it's a real minimum" << endl;
        
        int binNm2_X = conv::kBadResult;
        double chisq_X_Nm2 = conv::kBadResult;
        int binNp2_X = conv::kBadResult;
        double chisq_X_Np2 = conv::kBadResult;
        int binNm2_Y = conv::kBadResult;
        double chisq_Y_Nm2 = conv::kBadResult;
        int binNp2_Y = conv::kBadResult;
        double chisq_Y_Np2 = conv::kBadResult;
        
        
        // L(N-2) X
        binNm2_X = binN_X - 2;
        if(binNm2_X < fFirstMinimBin_X) {
          binNm2_X = fFirstMinimBin_X;  // should happens only when Xn = Xf
        }
        if(fIsFit2D)
        bincont = fChisqMargVArray->GetEntry(binNm2_X-1,binN_Y-1);
        else
        bincont = fChisqMargVArray->GetEntry(binNm2_X-1,fmin_bin_sub_param[1]);
        if( bincont<=0.) {
          if(fIsFit2D) bincont = CalcChisqMargPtRead(skobs_data,binNm2_X-1,binN_Y-1,readbinary);
          else bincont = CalcChisqMargPtRead(skobs_data,binNm2_X-1,fmin_bin_sub_param[1],readbinary);
        }
        chisq_X_Nm2 = bincont;
        
        // L(N+2) X
        binNp2_X = binN_X + 2;
        if(binNp2_X > fLastMinimBin_X)
        binNp2_X = fLastMinimBin_X;  // should happens only when Xn = Xl
        if(fIsFit2D)
        bincont = fChisqMargVArray->GetEntry(binNp2_X-1,binN_Y-1);
        else
        bincont = fChisqMargVArray->GetEntry(binNp2_X-1,fmin_bin_sub_param[1]);
        if( bincont<=0.) {
          if(fIsFit2D) bincont = CalcChisqMargPtRead(skobs_data,binNp2_X-1,binN_Y-1,readbinary);
          else bincont = CalcChisqMargPtRead(skobs_data,binNp2_X-1,fmin_bin_sub_param[1],readbinary);
        }
        chisq_X_Np2 = bincont;
        
        if(fIsFit2D){
          
          // L(N-2) Y
          binNm2_Y = binN_Y - 2;
          if(binNm2_Y < fFirstMinimBin_Y)
          binNm2_Y = fFirstMinimBin_Y;  // should happens only when Xn = Xf
          bincont = fChisqMargVArray->GetEntry(binN_X-1,binNm2_Y-1);
          if( bincont<=0.)
          bincont = CalcChisqMargPtRead(skobs_data,binN_X-1,binNm2_Y-1,readbinary);
          chisq_Y_Nm2 = bincont;
          
          // L(N+2) Y
          binNp2_Y = binN_Y + 2;
          if(binNp2_Y > fLastMinimBin_Y) binNp2_Y = fLastMinimBin_Y;  // should happens only when Xn = Xl
          bincont = fChisqMargVArray->GetEntry(binN_X-1,binNp2_Y-1);
          if( bincont<=0.)
          bincont = CalcChisqMargPtRead(skobs_data,binN_X-1,binNp2_Y-1,readbinary);
          chisq_Y_Np2 = bincont;
        }
        
#if defined __print_minimization_steps__
        cout << endl
        << "chisq_X = "
        << chisq_X_Nm2 << " "
        << chisq_X[ilike_Nm1] << " "
        << chisq_X[ilike_N] << " "
        << chisq_X[ilike_Np1] << " "
        << chisq_X_Np2 << " "
        << endl;
        if(fIsFit2D){
          cout << "chisq_Y = "
          << chisq_Y_Nm2 << " "
          << chisq_Y[ilike_Nm1] << " "
          << chisq_Y[ilike_N] << " "
          << chisq_Y[ilike_Np1] << " "
          << chisq_Y_Np2 << " "
          << endl;
        }
#endif
        
        if( chisq_X[ilike_N] <= chisq_X_Nm2 &&
           chisq_X[ilike_N] <= chisq_X_Np2 &&
           (!fIsFit2D ||
            (chisq_Y[ilike_N] <= chisq_Y_Nm2 &&
             chisq_Y[ilike_N] <= chisq_Y_Np2)
            )
           ){
          
          Nminima++;
          
          // It's a real minimum
          cout << "It's a real minimum" << endl;
          
          // Check if another minimum was found before
          if(chisq_X[ilike_N] <= fchisq_bf){
            
            // It's a global minimum. Stop the fit.
            
            cout << "It's a global minimum" << endl;
            
            fchisq_bf = chisq_X[ilike_N];
            
            if(fIsFit2D){
              fparam_bf[0] = fhchisq_marg_2D->GetXaxis()->GetBinCenter(binN_X);
              fparam_bf[1] = fhchisq_marg_2D->GetYaxis()->GetBinCenter(binN_Y);
            }
            else{
              fparam_bf[0] = fhchisq_marg_1D->GetBinCenter(binN_X);
            }
            
          }
          
          else{
            
            // If it's a local minimum the fit has to be restarted otherwise
            // after the first minimum the fit should be stopped.
            // Indeed this case can happen only if the first minimum was on the border.
            // Restart the fit and look for another minimum not on the border
            
            cout << endl;
            cout << "A new minimum is found but it's local" << endl;
            cout << "The previous minimum is deeper:" << endl;
            cout << "Global Chisq = " << fchisq_bf << endl;
            cout << "Local Chisq = " << chisq_X[ilike_N] << endl;
            cout << "Restart the fit" << endl;
            cout << endl;
            
            binF_X = UniFirstBin_X;
            binL_X = UniLastBin_X;
            
            if(fIsFit2D){
              binF_Y = UniFirstBin_Y;
              binL_Y = UniLastBin_Y;
            }
            
            Niterations_Part = 0; // Reset # of iterations per fit
            Nfits++;
            doNewFit=true;
            
            continue;
          }
          
          
          // If on the border do the fit again
          // The reached best fit is stored and will be the global
          // one if nothing better is found
          if(fBorderDoFit &&
             ( binN_X == UniFirstBin_X ||
              binN_X == UniLastBin_X ||
              binN_Y == UniFirstBin_Y ||
              binN_Y == UniLastBin_Y )
             ){
            
            binF_X = UniFirstBin_X;
            binL_X = UniLastBin_X;
            
            if(fIsFit2D){
              binF_Y = UniFirstBin_Y;
              binL_Y = UniLastBin_Y;
            }
            
            // Set the new interval for uniform extraction
            // excluding the border where the minimum was found
            
            if(binN_X == UniFirstBin_X){ // on the first bin X
                                         // next fit will not account the first bin X for the minimization
              if( (UniFirstBin_X+1) > UniLastBin_X)
              cout << "Don't exclude the border, otherwise FirstBin_X > LastBin_X" << endl;
              else
              UniFirstBin_X = UniFirstBin_X+1;
              binF_X = UniFirstBin_X;
            }
            else if(binN_X == fLastMinimBin_X){ // on the first bin X
                                                // next fit will not account the last bin X for the minimization
              if( (UniLastBin_X-1) < UniFirstBin_X)
              cout << "Don't exclude the border, otherwise LastBin_X < FirstBin_X" << endl;
              else
              UniLastBin_X = UniLastBin_X-1;
              binL_X = UniLastBin_X;
            }
            
            if( !fIsFit2D || (binN_Y == fFirstMinimBin_Y) ){ // on the first bin Y
                                                             // next fit will not account the first bin Y for the minimization
              if( (UniFirstBin_Y+1) > UniLastBin_Y)
              cout << "Don't exclude the border, otherwise FirstBin_Y > LastBin_Y" << endl;
              else
              UniFirstBin_Y = UniFirstBin_Y+1;
              binF_Y = UniFirstBin_Y;
            }
            else if( !fIsFit2D || (binN_Y == fLastMinimBin_Y) ){ // on the first bin Y
                                                                 // next fit will not account the last bin Y for the minimization
              if( (UniLastBin_Y-1) < UniFirstBin_Y)
              cout << "Don't exclude the border, otherwise LastBin_Y < FirstBin_Y" << endl;
              else
              UniLastBin_Y = UniLastBin_Y-1;
              binL_Y = UniLastBin_Y;
            }
            
            cout << endl;
            cout << "The reached minimum is on the border" << endl;
            cout << "To be sure it's the global minimum restart the fit" << endl;
            cout << "The border is excluded from the next fit" << endl;
            cout << "Uniform [ " << binF_X << " , " << binL_X << " ] " << endl;
            if(fIsFit2D)
            cout << "Uniform [ " << binF_Y << " , " << binL_Y << " ] " << endl;
            
            fMaxFitNumber = TMath::Min(fBorderMaxFitNumber,fMaxFitNumber);
            
            cout << "Repeat the fit until " << fMaxFitNumber << " times" << endl;
            cout << "# of performed fits = " << Nfits << endl;
            cout << endl;
            
            Niterations_Part = 0; // Reset # of iterations per fit
            Nfits++; // Increment Nfits. A minimum is found
            doNewFit = true; // The next one is a new fit
            
            continue;
          } // if border fit
          
          break;
          
        } // if it's a real minimum and not a fluctuation
        
        else{
          cout << endl;
          cout << "The reached minimum is a fluctuation of the bin" << endl;
          cout << "Restart the fit" << endl;
          cout << endl;
          
          binF_X = UniFirstBin_X;
          binL_X = UniLastBin_X;
          
          if(fIsFit2D){
            binF_Y = UniFirstBin_Y;
            binL_Y = UniLastBin_Y;
          }
          
          // Don't increment Nfits but continue the fit
          
          continue;
        }
      } // if chisq_X == chisq_Y
      else {
        cout << endl;
        cout << "Problem with convergence" << endl;
        cout << "Minimum is found but chisq_X != chisq_Y" << endl;
        cout << "Restart the fit" << endl;
        cout << endl;
        
        binF_X = UniFirstBin_X;
        binL_X = UniLastBin_X;
        
        if(fIsFit2D){
          binF_Y = UniFirstBin_Y;
          binL_Y = UniLastBin_Y;
        }
        
        // Don't increment Nfits but continue the fit
        
        continue;
      }
    } // if chisq_X<=chisq_X[Np1] && chisq_Y<=chisq_Y[Np1] ...
    
    else {
      
      int offdiagbin_X = conv::kBadResult; // X off-diagonal bin in the gradient direction
      int offdiagbin_Y = conv::kBadResult; // Y off-diagonal bin in the gradient direction
      
      // Gradient along X
      
      // If minimum along one direction (X), keep it fixed
      if( fUncorrXY &&
         chisq_X[ilike_N] <= chisq_X[ilike_Np1] &&
         chisq_X[ilike_N] <= chisq_X[ilike_Nm1]
         )
      {
        binF_X = bin_X;
        binL_X = bin_X;
      }
      
      // Move toward binL (right)
      else if(
              chisq_X[ilike_N] > chisq_X[ilike_Np1] &&
              chisq_X[ilike_N] < chisq_X[ilike_Nm1]
              ){
#if defined __print_minimization_steps__
        cout << "X: Move toward last bin (right)"<< endl;
#endif
        offdiagbin_X = ilike_Np1;
        binF_X = bin_X; // the extracted bin --> use bin_X as limit which doesn't contain the entire bin values
      }
      // Move toward Xf (left)
      else if(
              chisq_X[ilike_N] < chisq_X[ilike_Np1] &&
              chisq_X[ilike_N] > chisq_X[ilike_Nm1]
              ){
#if defined __print_minimization_steps__
        cout << "X: Move toward first bin (left)"<< endl;
#endif
        offdiagbin_X = ilike_Nm1;
        binL_X = bin_X;
      }
      
      // Case chisq[ilike_N] is in the maximum
      else if(
              chisq_X[ilike_N] >= chisq_X[ilike_Np1] &&
              chisq_X[ilike_N] >= chisq_X[ilike_Nm1]
              ){
#if defined __print_minimization_steps__
        cout << "X: chisq in the maximum" << endl;
#endif
        
        // Move toward Xl (right)
        if ( chisq_X[ilike_Np1] < chisq_X[ilike_Nm1] ){
#if defined __print_minimization_steps__
          cout << "X: Move toward last bin (right)"<< endl;
#endif
          binF_X = bin_X;
        }
        
        // Move toward Xf (left)
        else if ( chisq_X[ilike_Np1] > chisq_X[ilike_Nm1] ){
#if defined __print_minimization_steps__
          cout << "X: Move toward first bin (left)"<< endl;
#endif
          binL_X = bin_X;
        }
      }
      
      
      // Gradient along Y
      
      if(fIsFit2D){
        
        // If minimum along one direction (Y), keep it fixed
        if( fUncorrXY &&
           chisq_Y[ilike_N] <= chisq_Y[ilike_Np1] &&
           chisq_Y[ilike_N] <= chisq_Y[ilike_Nm1]
           )
        {
          binF_Y = bin_Y;
          binL_Y = bin_Y;
        }
        
        // Move toward Xl (right)
        else if(
                chisq_Y[ilike_N] > chisq_Y[ilike_Np1] &&
                chisq_Y[ilike_N] < chisq_Y[ilike_Nm1]
                ){
#if defined __print_minimization_steps__
          cout << "Y: Move toward last bin (up)"<< endl;
#endif
          offdiagbin_Y = ilike_Np1;
          binF_Y = bin_Y;
        }
        // Move toward Xf (left)
        else if(
                chisq_Y[ilike_N] < chisq_Y[ilike_Np1] &&
                chisq_Y[ilike_N] > chisq_Y[ilike_Nm1]
                ){
#if defined __print_minimization_steps__
          cout << "Y: Move toward first bin (down)"<< endl;
#endif
          offdiagbin_Y = ilike_Nm1;
          binL_Y = bin_Y;
        }
        // Case chisq[ilike_N] is in the
        else if(
                chisq_Y[ilike_N] >= chisq_Y[ilike_Np1] &&
                chisq_Y[ilike_N] >= chisq_Y[ilike_Nm1]
                ){
#if defined __print_minimization_steps__
          cout << "Y: chisq in the maximum" << endl;
#endif
          // Move toward Xl (right)
          if ( chisq_Y[ilike_Np1] < chisq_Y[ilike_Nm1] ){
#if defined __print_minimization_steps__
            cout << "Y: Move toward last bin (up)"<< endl;
#endif
            binF_Y = bin_Y;
          }
          // Move toward Xf (left)
          else if ( chisq_Y[ilike_Np1] > chisq_Y[ilike_Nm1] ){
#if defined __print_minimization_steps__
            cout << "Y: Move toward first bin (down)"<< endl;
#endif
            binL_Y = bin_Y;
          }
        }
      }
      
      // Check off-diagonal bin along the gradient direction (only for 2D fit)
      // It must have a Chisq(off-diag) < Chisq(N)
      // Otherwise go to another bin
      
      if(fIsFit2D){
        
        cout << "off-diagonal bin is X = " << offdiagbin_X << " - Y = " << offdiagbin_Y << endl;
        cout << "Chisq = " << chisq[offdiagbin_X][offdiagbin_Y] << endl;
        
        // If condition not satisfied check another point in the same interval
        if( chisq[ilike_N][ilike_N] < chisq[offdiagbin_X][offdiagbin_Y] )
        continue;
      }
      
#if defined __print_minimization_steps__
      cout << "Uniform [ " << binF_X << " , " << binL_X << " ] " << endl;
      if(fIsFit2D)
      cout << "Uniform [ " << binF_Y << " , " << binL_Y << " ] " << endl;
#endif
      
    } // end move to next steps
    
    // Saddle point is reached and the fit doesn't move --> restart the fit
    //if( binF_X == binL_X && binF_Y == binL_Y ){
    if( binF_X == binL_X ){
      binF_X = UniFirstBin_X;
      binL_X = UniLastBin_X;
      
      if(fIsFit2D && binF_Y == binL_Y){
        binF_Y = UniFirstBin_Y;
        binL_Y = UniLastBin_Y;
      }
    }
    
    // Partial maximum # of iterations. If exceeded restart the fit
    if(Niterations_Part > fMaxIterationSingleFit){
      cout << endl;
      cout << "Partial # of iterations > " << fMaxIterationSingleFit << endl;
      cout << "Convergence has not found --> restart the fit" << endl;
      cout << endl;
      //
      binF_X = UniFirstBin_X;
      binL_X = UniLastBin_X;
      
      if(fIsFit2D){
        binF_Y = UniFirstBin_Y;
        binL_Y = UniLastBin_Y;
      }
      //
      Niterations_Part = 0;
      Nfits++;
      doNewFit = true; // The next one is a new fit
    }
    
    Niterations++;
    Niterations_Part++;
    
  } // end while
  
  fNFits = Nfits;
  fNiterations = Niterations;
  
  cout << "The fit is well performed" << endl;
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::MinimizeChisqStepRead(SuperkObservation *skobs_data,bool readbinary){
  
  // Descent gradient algorithm is implemented in discrete space
  
  // readbinary = true  --> minimize chisq reading a binary file at each minimization step
  // readbinary = false --> first read the binary file in a given range then minimize the chisq
  //                        in that range
  
  // MAP
  GridStartBinRead(skobs_data,readbinary); // Get the starting bin
                                           //
  
  double bigvalueUOflow = 99999999.; // NEW PRECISION
  
  // 1D
  // fhchisq_marg_1D_Step = (TH1D*)fhchisq_marg_1D->Clone("hchisq_marg_1D_Step"); // NEW PRECISION
  // fhchisq_marg_1D_Step->Reset();
  // fhchisq_marg_1D_Step->SetName("hchisq_marg_1D_Step");
  // fhchisq_marg_1D_Step->SetTitle("hchisq_marg_1D_Step");
  fhchisq_marg_1D_Step->Reset();
  for(int ibinX=0;ibinX<fhchisq_marg_1D_Step->GetNbinsX();ibinX++){
    fhchisq_marg_1D_Step->SetBinContent(ibinX+1,0.);
  }
  // Set under and over flow of X to very high value
  fhchisq_marg_1D->SetBinContent(0,bigvalueUOflow);
  fhchisq_marg_1D->SetBinContent(fhchisq_marg_1D->GetNbinsX()+1,bigvalueUOflow);
  
  // 2D
  // fhchisq_marg_2D_Step = (TH2D*)fhchisq_marg_2D->Clone("hchisq_marg_2D_Step"); // NEW PRECISION
  // fhchisq_marg_2D_Step->Reset();
  // fhchisq_marg_2D_Step->SetName("hchisq_marg_2D_Step");
  // fhchisq_marg_2D_Step->SetTitle("hchisq_marg_2D_Step");
  fhchisq_marg_2D_Step->Reset();
  for(int ibinX=0;ibinX<fhchisq_marg_2D_Step->GetNbinsX();ibinX++){
    for(int ibinY=0;ibinY<fhchisq_marg_2D_Step->GetNbinsY();ibinY++){
      fhchisq_marg_2D_Step->SetBinContent(ibinX+1,ibinY+1,0.);
    }
  }
  // Set under and over flow of X and Y to very high value
  for(int ibinX=0;ibinX<fhchisq_marg_2D->GetNbinsX();ibinX++){
    fhchisq_marg_2D->SetBinContent(ibinX+1,0,bigvalueUOflow);
    fhchisq_marg_2D->SetBinContent(ibinX+1,fhchisq_marg_2D->GetNbinsY()+1,bigvalueUOflow);
  }
  for(int ibinY=0;ibinY<fhchisq_marg_2D->GetNbinsY();ibinY++){
    fhchisq_marg_2D->SetBinContent(0,ibinY+1,bigvalueUOflow);
    fhchisq_marg_2D->SetBinContent(fhchisq_marg_2D->GetNbinsX()+1,ibinY+1,bigvalueUOflow);
  }
  
  //
  
  //TRandom3 *r = new TRandom3(0); // needed to randomly scan the distribution
  
  //double UniFirstBin_X = fFirstMinimBin_X;
  //double UniLastBin_X  = fLastMinimBin_X;
  //double UniFirstBin_Y = fFirstMinimBin_Y;
  //double UniLastBin_Y  = fLastMinimBin_Y;
  
  int binN_X = fStartBin_X;
  int binN_Y = fStartBin_Y;
  
  double bincont = 0.;
  
  const int NchisqBin = 3;
  double chisq_X[NchisqBin] = {0., 0., 0.};  // Chisq(N-1), Chisq(N) and Chisq(N+1)
  double chisq_Y[NchisqBin] = {0., 0., 0.};  // Chisq(N-1), Chisq(N) and Chisq(N+1)
  double chisq[NchisqBin][NchisqBin] = {{0.}};  // Chisq(N-1), Chisq(N) and Chisq(N+1)
  const int ilike_N = 1;
  const int ilike_Nm1 = 0;
  const int ilike_Np1 = 2;
  
  // const int NchisqBin = 5;
  // double chisq_X[NchisqBin] = {0.};  // Chisq(N-1), Chisq(N) and Chisq(N+1)
  // double chisq_Y[NchisqBin] = {0.};  // Chisq(N-1), Chisq(N) and Chisq(N+1)
  // double chisq[NchisqBin][NchisqBin] = {{0.}}; // Chisq(N-1), Chisq(N) and Chisq(N+1)
  // const int ilike_N = 2;
  // const int ilike_Nm1 = 1;
  // const int ilike_Np1 = 3;
  
  int Niterations=0;
  int Niterations_Part=0;
  int Nminima=0;
  bool doNewFit = true;
  int Nfits=1; // number of fits performed
  
  
  while(1){
    
    if(Nfits > fMaxFitNumber){
      cout << endl;
      cout << "Starting fit # " << Nfits << endl;
      cout << "Maximum allowed # of fits: " << fMaxFitNumber << endl;
      cout << "Convergence has not found" << endl;
      cout << "# of found minima: " << Nminima << endl;
      cout << endl;
      return;
    }
    else if(doNewFit){
      cout << endl;
      cout << "/////////////// Fit " << Nfits << " ///////////////" << endl;
      cout << endl;
      doNewFit = false;
    }
    
    double binlowedge_X = conv::kBadResult;
    double binlowedge_Y = conv::kBadResult;
    
    // Find direction along X axis
    
    //cout << "Find direction along X axis" << endl;
    
    
    // L(N,N)
    if(fIsFit2D)
    bincont = fChisqMargVArray->GetEntry(binN_X-1,binN_Y-1);
    else
    bincont = fChisqMargVArray->GetEntry(binN_X-1,fmin_bin_sub_param[1]);
    
    if( bincont<=0.){
      if(fIsFit2D) bincont = CalcChisqMargPtRead(skobs_data,binN_X-1,binN_Y-1,readbinary); // CalcChisqMargPt takes the point not bin
      else bincont = CalcChisqMargPtRead(skobs_data,binN_X-1,fmin_bin_sub_param[1],readbinary);
    }
    
    chisq_X[ilike_N] = bincont;
    chisq[ilike_N][ilike_N] = bincont;
    
    // L(N-1,N)
    int binNm1_X = binN_X - 1;
    // Bin is on the border, it should happens only when binN = binF
    if(binNm1_X < fFirstMinimBin_X){
      //binNm1_X = binN_X;
      chisq_X[ilike_Nm1] = conv::kALargeNum;
      chisq[ilike_Nm1][ilike_N] = conv::kALargeNum;
    }
    // normal situation
    else{
      if(fIsFit2D)
      bincont = fChisqMargVArray->GetEntry(binNm1_X-1,binN_Y-1);
      else
      bincont = fChisqMargVArray->GetEntry(binNm1_X-1,fmin_bin_sub_param[1]);
      if( bincont<=0.){
        if(fIsFit2D) bincont = CalcChisqMargPtRead(skobs_data,binNm1_X-1,binN_Y-1,readbinary);
        else bincont = CalcChisqMargPtRead(skobs_data,binNm1_X-1,fmin_bin_sub_param[1],readbinary);
      }
      chisq_X[ilike_Nm1] = bincont;
      chisq[ilike_Nm1][ilike_N] = bincont;
    }
    
    // L(N+1,N)
    int binNp1_X = binN_X + 1;
    // Bin is on the border, it should happens only when binN = binF
    if(binNp1_X > fLastMinimBin_X){
      //binNp1_X = binN_X;
      chisq_X[ilike_Np1] = conv::kALargeNum;
      chisq[ilike_Np1][ilike_N] = conv::kALargeNum;
    }
    // normal situation
    else{
      if(fIsFit2D)
      bincont = fChisqMargVArray->GetEntry(binNp1_X-1,binN_Y-1);
      else
      bincont = fChisqMargVArray->GetEntry(binNp1_X-1,fmin_bin_sub_param[1]);
      if( bincont<=0.){
        if(fIsFit2D) bincont = CalcChisqMargPtRead(skobs_data,binNp1_X-1,binN_Y-1,readbinary);
        else bincont = CalcChisqMargPtRead(skobs_data,binNp1_X-1,fmin_bin_sub_param[1],readbinary); // bin_y = 0
      }
      chisq_X[ilike_Np1] = bincont;
      chisq[ilike_Np1][ilike_N] = bincont;
    }
    
    binlowedge_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(binN_X);
    
    // Find direction along Y axis
    
    //cout << "Find direction along Y axis" << endl;
    
    int binNm1_Y = conv::kBadResult;
    int binNp1_Y = conv::kBadResult;
    
    if(fIsFit2D){
      
      // L(N,N) already computed for 2D chisq
      bincont = fChisqMargVArray->GetEntry(binN_X-1,binN_Y-1);
      if( bincont<=0.)
      bincont = CalcChisqMargPtRead(skobs_data,binN_X-1,binN_Y-1,readbinary);
      chisq_Y[ilike_N] = bincont;
      
      // L(N,N-1)
      binNm1_Y = binN_Y - 1;
      // Bin is on the border, it should happens only when binN = binF
      if(binNm1_Y < fFirstMinimBin_Y){
        //binNm1_Y = binN_Y;
        chisq_Y[ilike_Nm1] = conv::kALargeNum;
        chisq[ilike_N][ilike_Nm1] = conv::kALargeNum;
      }
      // normal situation
      else{
        bincont = fChisqMargVArray->GetEntry(binN_X-1,binNm1_Y-1);
        
        if( bincont<=0.)
        bincont = CalcChisqMargPtRead(skobs_data,binN_X-1,binNm1_Y-1,readbinary);
        chisq_Y[ilike_Nm1] = bincont;
        chisq[ilike_N][ilike_Nm1] = bincont;
      }
      
      // L(N,N+1)
      binNp1_Y = binN_Y + 1;
      // Bin is on the border, it should happens only when binN = binF
      if(binNp1_Y > fLastMinimBin_Y){
        //binNp1_Y = binN_Y;
        chisq_Y[ilike_Np1] = conv::kALargeNum;
        chisq[ilike_N][ilike_Np1] = conv::kALargeNum;
      }
      // normal situation
      else{
        bincont = fChisqMargVArray->GetEntry(binN_X-1,binNp1_Y-1);
        if( bincont<=0.)
        bincont = CalcChisqMargPtRead(skobs_data,binN_X-1,binNp1_Y-1,readbinary);
        chisq_Y[ilike_Np1] = bincont;
        chisq[ilike_N][ilike_Np1] = bincont;
      }
      
      binlowedge_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(binN_Y);
    }
    
    
    // Calculate off-diagonal bins
    
    if(fIsFit2D){
      
      // L(N-1,N-1)
      if(binNm1_X < fFirstMinimBin_X || binNm1_Y < fFirstMinimBin_Y){
        chisq[ilike_Nm1][ilike_Nm1] = conv::kALargeNum;
      }
      else{
        bincont = fChisqMargVArray->GetEntry(binNm1_X-1,binNm1_Y-1);
        if( bincont<=0.)
        bincont = CalcChisqMargPtRead(skobs_data,binNm1_X-1,binNm1_Y-1,readbinary);
        chisq[ilike_Nm1][ilike_Nm1] = bincont;
      }
      
      // L(N+1,N+1)
      if(binNp1_X > fLastMinimBin_X || binNp1_Y > fLastMinimBin_Y){
        chisq[ilike_Np1][ilike_Np1] = conv::kALargeNum;
      }
      else{
        bincont = fChisqMargVArray->GetEntry(binNp1_X-1,binNp1_Y-1);
        if( bincont<=0.)
        bincont = CalcChisqMargPtRead(skobs_data,binNp1_X-1,binNp1_Y-1,readbinary);
        chisq[ilike_Np1][ilike_Np1] = bincont;
      }
      
      // L(N-1,N+1)
      if(binNm1_X < fFirstMinimBin_X || binNp1_Y > fLastMinimBin_Y){
        chisq[ilike_Nm1][ilike_Np1] = conv::kALargeNum;
      }
      else{
        bincont = fChisqMargVArray->GetEntry(binNm1_X-1,binNp1_Y-1);
        if( bincont<=0.)
        bincont = CalcChisqMargPtRead(skobs_data,binNm1_X-1,binNp1_Y-1,readbinary);
        chisq[ilike_Nm1][ilike_Np1] = bincont;
      }
      
      // L(N+1,N-1)
      if(binNp1_X > fLastMinimBin_X || binNm1_Y < fFirstMinimBin_Y){
        chisq[ilike_Np1][ilike_Nm1] = conv::kALargeNum;
      }
      else{
        bincont = fChisqMargVArray->GetEntry(binNp1_X-1,binNm1_Y-1);
        if( bincont<=0.)
        bincont = CalcChisqMargPtRead(skobs_data,binNp1_X-1,binNm1_Y-1,readbinary);
        chisq[ilike_Np1][ilike_Nm1] = bincont;
      }
      
    } // if fIsFit2D
    
    
#if defined __print_minimization_steps__
    cout << endl;
    cout << "# iteration = " << Niterations << endl;
    cout << "# iteration current = " << Niterations_Part << endl;
    
    if(!fIsFit2D){
      cout << "Bin: " << binN_X
      << " --> " << binlowedge_X << endl;
      cout << "Chisq grid around binN: " << endl;
      for(int i=0;i<NchisqBin;i++) cout << chisq_X[i] << " ";
    }
    else{
      cout << "Bin: [ " << binN_X << " , " << binN_Y << " ] "
      << " --> [ " << binlowedge_X << " , " << binlowedge_Y << " ] "
      << fhchisq_marg_2D->GetBinContent(binN_X,binN_Y)
      << endl;
      cout << "Chisq grid around binN: " << endl;
      //for(int i=2;i>=0;i--){
      for(int i=(NchisqBin-1);i>=0;i--){
        for(int j=0;j<NchisqBin;j++){
          cout << chisq[j][i] << " ";
        }
        cout << endl;
      }
      cout << endl;
    }
#endif
    
    double bincent_X = conv::kBadResult;
    double bincent_Y = conv::kBadResult;
    
    if(!fIsFit2D){
      bincent_X = fhchisq_marg_1D_Step->GetXaxis()->GetBinCenter(binN_X);
      fhchisq_marg_1D_Step->Fill(bincent_X);
    }
    else{
      bincent_X = fhchisq_marg_2D_Step->GetXaxis()->GetBinCenter(binN_X);
      bincent_Y = fhchisq_marg_2D_Step->GetYaxis()->GetBinCenter(binN_Y);
      fhchisq_marg_2D_Step->Fill(bincent_X,bincent_Y);
    }
    
    // Best fit bin is found
    // this condition contains the case with the boundary, when
    // z(N-1) = Z(N) or Z(N+1) = Z(N)
    
    if(
       chisq_X[ilike_N] <= chisq_X[ilike_Np1] &&
       chisq_X[ilike_N] <= chisq_X[ilike_Nm1] &&
       ( !fIsFit2D ||
        (
         chisq[ilike_N][ilike_N] <= chisq[ilike_N][ilike_Np1]   &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_N][ilike_Nm1]   &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_Nm1][ilike_N]   &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_Np1][ilike_N]   &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_Np1][ilike_Np1] &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_Nm1][ilike_Nm1] &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_Nm1][ilike_Np1] &&
         chisq[ilike_N][ilike_N] <= chisq[ilike_Np1][ilike_Nm1]
         )
        )
       ){
      
      cout << "A minimum is found" << endl;
      cout << "Best fit Chisq = " << chisq_X[ilike_N] << endl;
      if(!fIsFit2D){
        cout << "Best fit X = " << fhchisq_marg_1D->GetXaxis()->GetBinCenter(binN_X) << endl;
      }
      else{
        cout << "Best fit X = " << fhchisq_marg_2D->GetXaxis()->GetBinCenter(binN_X) << endl;
        cout << "Best fit Y = " << fhchisq_marg_2D->GetYaxis()->GetBinCenter(binN_Y) << endl;
      }
      cout << "Total # of iterations = " << Niterations << endl;
      cout << "Total # of performed fits = " << Nfits << endl;
      cout << "Check if it's a real minimum" << endl;
      
      
      // Set new chisq grid to check the minimum. Larger than the other one.
      
      // +-3 bins around the minimum
      
      const int NbinCheck = NUMOFBINSAROUNDMINTOCHECK;
      const int ilike_Min_N = NbinCheck; // it would be NbinCheck+0 (0 is the 1st elem of array)
      const int NchisqBin_Min = ilike_Min_N + NbinCheck + 1; // +1 because the NbinCheck is array elem.
      
      double chisq_Min_1D[NchisqBin_Min] = {0.};  // Chisq(N-1), Chisq(N) and Chisq(N+1)
      double bin_Min_1D[NchisqBin_Min] = {0.};  // bin(N-1), bin(N) and bin(N+1)
      
      double chisq_Min[NchisqBin_Min][NchisqBin_Min] = {{0.}};  // Chisq(N-1), Chisq(N) and Chisq(N+1)
      double bin_Min_X[NchisqBin_Min][NchisqBin_Min] = {{0.}};  // bin(N-1), bin(N) and bin(N+1)
      double bin_Min_Y[NchisqBin_Min][NchisqBin_Min] = {{0.}};  // bin(N-1), bin(N) and bin(N+1)
      
      for(int ichisq=0;ichisq<NchisqBin_Min;ichisq++){
        
        chisq_Min_1D[ichisq] = conv::kBadResult;
        bin_Min_1D[ichisq] = conv::kBadResult;
        
        for(int jchisq=0;jchisq<NchisqBin_Min;jchisq++){
          chisq_Min[ichisq][jchisq] = conv::kBadResult;
          bin_Min_X[ichisq][jchisq] = conv::kBadResult;
          bin_Min_Y[ichisq][jchisq] = conv::kBadResult;
        }
      }
      
      int NchisqBin_Loop_X = NchisqBin_Min;
      int NchisqBin_Loop_Y = 1;
      if(fIsFit2D) NchisqBin_Loop_Y = NchisqBin_Min;
      
      
      
      
      if(!fIsFit2D){
        for(int i=0;i<NchisqBin_Loop_X;i++){
          int shift_binN_1D = i - NbinCheck; // shift wrt the current bin (e.g. from -2 to +2)
          int binloop_1D = binN_X + shift_binN_1D;
          bin_Min_1D[i] = binloop_1D; // store the X bins
          
          // Bin is on the border, it should happens only when binN = binF
          bincont = conv::kALargeNum;
          if( (binloop_1D < fFirstMinimBin_X) ||
             (binloop_1D > fLastMinimBin_X)
             )
          {
            chisq_Min_1D[i] = conv::kALargeNum;
          }
          else{
            bincont = fChisqMargVArray->GetEntry(binloop_1D-1,fmin_bin_sub_param[1]);
            if( bincont<=0.)
            bincont = CalcChisqMargPtRead(skobs_data,binloop_1D-1,fmin_bin_sub_param[1],readbinary);
          }
          chisq_Min_1D[i] = bincont;
        }
      }
      
      else{
        for(int i=0;i<NchisqBin_Loop_X;i++){
          for(int j=0;j<NchisqBin_Loop_Y;j++){
            
            int shift_binN_X = i - NbinCheck; // shift wrt the current bin (e.g. from -2 to +2)
            int binloop_X = binN_X + shift_binN_X;
            int shift_binN_Y = j - NbinCheck; // shift wrt the current bin (e.g. from -2 to +2)
            int binloop_Y = binN_Y + shift_binN_Y;
            
            bin_Min_X[i][j] = binloop_X; // store the X bins
            bin_Min_Y[i][j] = binloop_Y; // store the Y bins
            
            // Bin is on the border, it should happens only when binN = binF
            if( (binloop_X < fFirstMinimBin_X) ||
               (binloop_X > fLastMinimBin_X)  ||
               (binloop_Y < fFirstMinimBin_Y) ||
               (binloop_Y > fLastMinimBin_Y) )
            {
              chisq_Min[i][j] = conv::kALargeNum;
            }
            // normal situation
            else{
              if(fIsFit2D)
              bincont = fChisqMargVArray->GetEntry(binloop_X-1,binloop_Y-1);
              else
              bincont = fChisqMargVArray->GetEntry(binloop_X-1,fmin_bin_sub_param[1]);
              if( bincont<=0.){
                if(fIsFit2D) bincont = CalcChisqMargPtRead(skobs_data,binloop_X-1,binloop_Y-1,readbinary);
                else bincont = CalcChisqMargPtRead(skobs_data,binloop_X-1,fmin_bin_sub_param[1],readbinary);
              }
              chisq_Min[i][j] = bincont;
            }
          } // j loop
        } // i loop
      } // fIsFit2D
      
#if defined __print_minimization_steps__
      if(fIsFit2D){
        for(int i=(NchisqBin_Loop_X-1);i>=0;i--){
          for(int j=0;j<NchisqBin_Loop_Y;j++){
            cout << chisq_Min[j][i] << " ";
          }
          cout << endl;
        }
      }
      else{
        for(int i=0;i<NchisqBin_Loop_X;i++){
          cout << chisq_Min_1D[i] << " ";
        }
        cout << endl;
      }
#endif
      
      bool IsRealMin = true;
      
      if(!fIsFit2D){
        for(int i=0;i<NchisqBin_Loop_X;i++){
          if( chisq_Min_1D[ilike_Min_N] > chisq_Min_1D[i]){  // Store the minimum Chisq bin
            IsRealMin = false;
            binN_X = bin_Min_1D[i]; // store the X bins
          }
        }
      }
      else{
        for(int i=0;i<NchisqBin_Loop_X;i++){
          for(int j=0;j<NchisqBin_Loop_Y;j++){
            if(i!=j){
              if( chisq_Min[ilike_Min_N][ilike_Min_N] > chisq_Min[i][j]){  // Store the minimum Chisq bin
                IsRealMin = false;
                binN_X = bin_Min_X[i][j]; // store the X bins
                if(fIsFit2D) binN_Y = bin_Min_Y[i][j]; // store the Y bins
              }
            }
          }
        }
      }
      
      if(IsRealMin){
        
        Nminima++;
        
        // It's a real minimum
        cout << "It's a real minimum" << endl;
        
        fchisq_bf = chisq_X[ilike_N];
        if(fIsFit2D) fchisq_bf = chisq[ilike_N][ilike_N];
        
        if(fIsFit2D){
          fparam_bf[0] = fhchisq_marg_2D->GetXaxis()->GetBinCenter(binN_X);
          fparam_bf[1] = fhchisq_marg_2D->GetYaxis()->GetBinCenter(binN_Y);
        }
        else{
          fparam_bf[0] = fhchisq_marg_1D->GetBinCenter(binN_X);
        }
        
        Nminima++;
        
        fFitStatus = 1;
        
        break;
      }// global minimum
      
      
      // It's a bin fluctuation or fake minimum due to binning
      
      else{
        cout << endl;
        cout << "The reached minimum is a fluctuation of the bin" << endl;
        
        cout << "Continue the fit to the lowest chisq bin" << endl;
        //cout << "Restart the fit from a new starting point:" << endl;
        
        cout << endl;
        
        // // Throw the bin number
        // double r_X = r->Uniform(UniFirstBin_X,UniLastBin_X+1); // +1 because it's then truncated into a int. Extract 21th bin as well
        // double r_Y = r->Uniform(UniFirstBin_Y,UniLastBin_Y+1); // +1 because it's then truncated into a int. Extract 21th bin as well
        // binN_X = (int)r_X;
        // binN_Y = (int)r_Y;
        
        double bincontmin = conv::kALargeNum;
        if(fIsFit2D) bincontmin = CalcChisqMargPtRead(skobs_data,binN_X-1,binN_Y-1,readbinary);
        else bincontmin = CalcChisqMargPtRead(skobs_data,binN_X-1,fmin_bin_sub_param[1],readbinary);
        
        cout << "Restart the fit from a new starting point:" << endl;
        cout << "chisq = " << bincontmin << endl;
        cout << "binN_X = " << binN_X << " - binN_Y = " << binN_Y << endl;
        cout << "X = " << fhchisq_marg_2D->GetXaxis()->GetBinCenter(binN_X);
        if(fIsFit2D)
        cout << " - Y = " << fhchisq_marg_2D->GetYaxis()->GetBinCenter(binN_Y) << endl;
        cout << endl;
        
        continue;
        
      } // if it's a real minimum or a bin fluctuation
      
    } // if chisq_X<=chisq_X[Np1] && chisq_Y<=chisq_Y[Np1] ...
    
    
    
    // Get the new direction
    
    else {
      
      // Move toward the lowest adjacent Chisq bin
      
      int NewDir_X = conv::kBadResult; // It's not the bin number but an array index
      int NewDir_Y = conv::kBadResult; // It's not the bin number but an array index
      double min = conv::kALargeNum; // NEW PRECISION
      
      if(!fIsFit2D){
        for(int i=0;i<NchisqBin;i++){
          if( chisq_X[i] < min ){
            NewDir_X = i;
            min = chisq_X[i];
          }
        }
      }
      else{
        for(int i=0;i<NchisqBin;i++){
          for(int j=0;j<NchisqBin;j++){
            if( chisq[i][j] < min ){
              NewDir_X = i;
              NewDir_Y = j;
              min = chisq[i][j];
            }
          }
        }
      } // if 2D fit
      
#if defined __print_minimization_steps__
      
      cout << "New direction bin: " << endl;
      cout << "  - X: ";
      if(NewDir_X == ilike_Nm1) cout << "left";
      else if(NewDir_X == ilike_Np1) cout << "right";
      else if(NewDir_X == ilike_N) cout << "stop";
      else{
        cerr << endl << endl;
        cerr << "FitMarginalization::MinimizeChisqStepRead" << endl;
        cerr << "No new direction is found" << endl;
        cerr << "NewDir_X = " << NewDir_X << endl;
        cerr << endl;
        util::Die(1);
      }
      cout << endl;
      
      if(fIsFit2D){
        cout << "  - Y: ";
        if(NewDir_Y == ilike_Nm1) cout << "down";
        else if(NewDir_Y == ilike_Np1) cout << "up";
        else if(NewDir_Y == ilike_N) cout << "stop";
        else{
          cerr << endl << endl;
          cerr << "FitMarginalization::MinimizeChisqStepRead" << endl;
          cerr << "No new direction is found" << endl;
          cerr << "NewDir_Y = " << NewDir_Y << endl;
          cerr << endl;
          util::Die(1);
        }
        cout << endl;
      }
#endif
      
      if(
         (!fIsFit2D && NewDir_X == ilike_N) ||
         (fIsFit2D && NewDir_X == ilike_N && NewDir_Y == ilike_N)
         ){
        cerr << endl;
        cerr << "FitMarginalization::MinimizeChisqStepRead" << endl;
        cerr << "It's a minimum, you cannot find a new direction" << endl;
        cerr << "There's something wrong" << endl;
        cerr << endl;
        util::Die(1);
      }
      
      // Compute the gradient between the lowest bin and the opposite one and the new step size
      
      //double grad_X=0.;
      //double grad_Y=0.;
      //double step_X=0.;
      //double step_Y=0.;
      //double num=0.;
#if defined __print_minimization_steps__
      double den=0.;
#endif
      
      if(!fIsFit2D){
        
        int NextBin_X = binN_X + (NewDir_X-1); // NewDir=0,1,2. NewDir=1 if binN
        
#if defined __print_minimization_steps__
        double NextVal_X = fhchisq_marg_1D->GetXaxis()->GetBinCenter(NextBin_X);
        double BinVal_X = binlowedge_X;
        double NewVal_X = fhchisq_marg_1D->GetXaxis()->GetBinCenter(NextBin_X);
        
        cout << "binN_X = " << binN_X << " -  NextBin_X = " << NextBin_X << endl;
        cout << "BinVal_X = " << BinVal_X << endl;
        cout << "NextVal_X = " << NextVal_X << endl;
        cout << "Chisq_New - Chisq_N = " << den << endl;
        cout << "New point: " << endl;
        cout << "  - NewVal_X = " << NewVal_X << endl;
        cout << "New bin: " << endl;
        cout << "  - NextBin_X = " << NextBin_X << endl;
        cout << endl;
#endif
        
        binN_X = NextBin_X;
      }
      else{
        
        int NextBin_X = binN_X + (NewDir_X-1); // NewDir=0,1,2. NewDir=1 if binN
        int NextBin_Y = binN_Y + (NewDir_Y-1);
        
#if defined __print_minimization_steps__
        double NextVal_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(NextBin_X);
        double NextVal_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(NextBin_Y);
        
        double BinVal_X = binlowedge_X;
        double BinVal_Y = binlowedge_Y;
        
        double NewVal_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(NextBin_X);
        double NewVal_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(NextBin_Y);
        
        cout << "binN_X = " << binN_X << " -  NextBin_X = " << NextBin_X << endl;
        cout << "binN_Y = " << binN_Y << " -  NextBin_Y = " << NextBin_Y << endl;
        cout << "BinVal_X = " << BinVal_X << " - NextVal_X = " << NextVal_X << endl;
        cout << "BinVal_Y = " << BinVal_Y << " - NextVal_Y = " << NextVal_Y << endl;
        cout << "NextVal_X = " << NextVal_X << " - NextVal_Y = " << NextVal_Y << endl;
        cout << "Chisq_New - Chisq_N = " << den << endl;
        cout << "New point: " << endl;
        cout << "  - NewVal_X = " << NewVal_X << endl;
        cout << "  - NewVal_Y = " << NewVal_Y << endl;
        cout << "New bin: " << endl;
        cout << "  - NextBin_X = " << NextBin_X << endl;
        cout << "  - NextBin_Y = " << NextBin_Y << endl;
        cout << endl;
#endif
        
        // Get the bin of the next iteration
        
        binN_X = NextBin_X;
        binN_Y = NextBin_Y;
        
      } // if 2D fit
      
      Niterations_Part++;
      Niterations++;
      
    } // if chisq_X<=chisq_X[Np1] && chisq_Y<=chisq_Y[Np1]... or get new direction
    
  } // while(1)
  
  fNFits = Nfits;
  fNiterations = Niterations;
  
  cout << "The fit is well performed" << endl;
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetNumStartBin(int Nthrows_X,int Nthrows_Y){
  fNumStartBin_X = Nthrows_X;
  
  if(!fIsFit2D) fNumStartBin_Y = 1;
  else fNumStartBin_Y = Nthrows_Y;
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetExtractStartPt(string doextract){
  fExtractStartPt = doextract;
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetGridStartRange(double minX,double maxX,double minY,double maxY){
  
  // fGridStartMinVal_X = minX+conv::kAVerySmallNum;
  // fGridStartMaxVal_X = maxX+conv::kAVerySmallNum;
  // fGridStartMinVal_Y = minY+conv::kAVerySmallNum;
  // fGridStartMaxVal_Y = maxY+conv::kAVerySmallNum;
  fGridStartMinVal_X = minX;
  fGridStartMaxVal_X = maxX;
  fGridStartMinVal_Y = minY;
  fGridStartMaxVal_Y = maxY;
  
  if(!fIsFit2D){
    //fGridStartMinBin_X = fhchisq_marg_1D->FindBin(minX+conv::kAVerySmallNum);
    //fGridStartMaxBin_X = fhchisq_marg_1D->FindBin(maxX+conv::kAVerySmallNum);
    fGridStartMinBin_X = fhchisq_marg_1D->FindBin(minX);
    fGridStartMaxBin_X = fhchisq_marg_1D->FindBin(maxX);
    fGridStartMinBin_Y = conv::kBadResult;
    fGridStartMaxBin_Y = conv::kBadResult;
  }
  else{
    // fGridStartMinBin_X = fhchisq_marg_2D->GetXaxis()->FindBin(minX+conv::kAVerySmallNum);
    // fGridStartMaxBin_X = fhchisq_marg_2D->GetXaxis()->FindBin(maxX+conv::kAVerySmallNum);
    // fGridStartMinBin_Y = fhchisq_marg_2D->GetYaxis()->FindBin(minY+conv::kAVerySmallNum);
    // fGridStartMaxBin_Y = fhchisq_marg_2D->GetYaxis()->FindBin(maxY+conv::kAVerySmallNum);
    fGridStartMinBin_X = fhchisq_marg_2D->GetXaxis()->FindBin(minX);
    fGridStartMaxBin_X = fhchisq_marg_2D->GetXaxis()->FindBin(maxX);
    fGridStartMinBin_Y = fhchisq_marg_2D->GetYaxis()->FindBin(minY);
    fGridStartMaxBin_Y = fhchisq_marg_2D->GetYaxis()->FindBin(maxY);
  }
  return;
}
//________________________________________________________________________________________
void FitMarginalization::GridStartBinRead(SuperkObservation *skobs_data,bool readbinary)
{
  // readbinary = true  --> minimize chisq reading a binary file at each minimization step
  // readbinary = false --> first read the binary file in a given range then minimize the chisq
  //                        in that range
  
  const int NThrowsExtrBin = fNumStartBin_X*fNumStartBin_Y;
  
  int UniFirstBin_X = conv::kBadResult;
  int UniLastBin_X = conv::kBadResult;
  int UniFirstBin_Y = conv::kBadResult;
  int UniLastBin_Y = conv::kBadResult;
  
  double UniFirstVal_X = conv::kBadResult;
  double UniLastVal_X = conv::kBadResult;
  double UniFirstVal_Y = conv::kBadResult;
  double UniLastVal_Y = conv::kBadResult;
  
  // If read binary file at each grid point step, the range must be contained in full range
  // while if binary files are first stored and then the chisq is computed, since a binary
  // file is read and then the corresponding part of chisq is filled, the range must
  // be consistent with the read binary file. The same is when minimization is done without
  // without reading binary at each step: the minimization range must be contained in the
  // binary file one.
  
  // Read the binary file and then minimize all the contained parameter range
  if(!readbinary){
    UniFirstBin_X = fFirstMinimBin_X;
    UniLastBin_X  = fLastMinimBin_X;
    UniFirstBin_Y = fFirstMinimBin_Y;
    UniLastBin_Y  = fLastMinimBin_Y;
    
    UniFirstVal_X = fFirstMinimVal_X;
    UniLastVal_X  = fLastMinimVal_X;
    UniFirstVal_Y = fFirstMinimVal_Y;
    UniLastVal_Y  = fLastMinimVal_Y;
  }
  // Read binary file at each minimization step
  else{
    UniFirstBin_X = fGridStartMinBin_X;
    UniLastBin_X  = fGridStartMaxBin_X;
    UniFirstBin_Y = fGridStartMinBin_Y;
    UniLastBin_Y  = fGridStartMaxBin_Y;
    
    UniFirstVal_X = fGridStartMinVal_X;
    UniLastVal_X  = fGridStartMaxVal_X;
    UniFirstVal_Y = fGridStartMinVal_Y;
    UniLastVal_Y  = fGridStartMaxVal_Y;
    
    // Check it's in the minimization interval
    if( (fGridStartMinBin_X < fFirstMinimBin_X) ||
       (fGridStartMinBin_X > fLastMinimBin_X) ||
       (fGridStartMaxBin_X < fFirstMinimBin_X) ||
       (fGridStartMaxBin_X > fLastMinimBin_X)
       )
    {
      cerr << endl;
      cerr << "FitMarginalization::GridStartBinRead" << endl;
      cerr << "Wrong input range param0 to get the minimization starting points" << endl;
      cerr << endl;
      util::Die(1);
    }
    
    
    if(fIsFit2D){
      if( (fGridStartMinBin_Y < fFirstMinimBin_Y) ||
         (fGridStartMinBin_Y > fLastMinimBin_Y)  ||
         (fGridStartMaxBin_Y < fFirstMinimBin_Y) ||
         (fGridStartMaxBin_Y > fLastMinimBin_Y) )
      {
        cerr << endl;
        cerr << "FitMarginalization::GridStartBinRead" << endl;
        cerr << "Wrong input range param1 to get the minimization starting points" << endl;
        cerr << endl;
        util::Die(1);
      }
    }
  }
  
  cout << endl;
  cout << "X: " << UniFirstBin_X << " " << UniLastBin_X << endl;
  cout << "Y: " << UniFirstBin_Y << " " << UniLastBin_Y << endl;
  cout << endl;
  
  double bincont = conv::kBadResult;
  
  double ChisqMinStart = conv::kALargeNum;
  int binStart_X = fStartBin_X;
  int binStart_Y = fStartBin_Y;
  
  TRandom3 *r = new TRandom3(0);
  
  //double like_iter = conv::kBadResult;
  //double like_prev = conv::kBadResult;
  
  
  // Single starting point
  
  if(fExtractStartPt=="Single"){
    cout << "Single starting point:" << endl;
    ChisqMinStart = CalcChisqMargPtRead(skobs_data,fStartBin_X-1,fStartBin_Y-1,readbinary);
  }
  
  // Random extraction
  
  else if(fExtractStartPt=="GridRnd"){
    
    cout << "Find the starting point from a extracted grid of ";
    cout << NThrowsExtrBin << " points in the range:"<< endl;
    cout << "X: Uniform [ " << UniFirstBin_X << " , " << UniLastBin_X << " ] " << endl;
    if(fIsFit2D)
    cout << "Y: Uniform [ " << UniFirstBin_Y << " , " << UniLastBin_Y << " ] " << endl;
    
    for(int iStartPt=0;iStartPt<NThrowsExtrBin;iStartPt++){
      
      double r_X = r->Uniform(UniFirstBin_X,UniLastBin_X+1); // +1 because it's then truncated into a int. Extract 21th bin as well
      double r_Y = r->Uniform(UniFirstBin_Y,UniLastBin_Y+1); // +1 because it's then truncated into a int. Extract 21th bin as well
      int bs_x = (int)r_X;
      int bs_y = (int)r_Y;
      //cout << "bs_x = " << bs_x << " - bs_y = " << bs_y << endl;
      
      if(!fIsFit2D)
      bincont = CalcChisqMargPtRead(skobs_data,bs_x-1,fmin_bin_sub_param[1],readbinary);
      else
      bincont = CalcChisqMargPtRead(skobs_data,bs_x-1,bs_y-1,readbinary);
      
      if(bincont < ChisqMinStart){
        ChisqMinStart = bincont;
        binStart_X = bs_x;
        binStart_Y = bs_y;
      }
    }
    fStartBin_X = binStart_X;
    fStartBin_Y = binStart_Y;
  }
  
  // Fix a grid of NThrowsExtrBin points for extraction
  
  else if(fExtractStartPt=="GridFix"){
    
    const int Npt_X = fNumStartBin_X;
    const int Npt_Y = fNumStartBin_Y;
    
    cout << "Find the starting point from a fixed grid of " << endl;
    cout << Npt_X << " X " << Npt_Y << " points" << endl;
    cout << "in the range:"<< endl;
    cout << "X: [ " << UniFirstBin_X << " , " << UniLastBin_X << " ] " << endl;
    if(fIsFit2D)
    cout << "Y: [ " << UniFirstBin_Y << " , " << UniLastBin_Y << " ] " << endl;
    
    int Nsteps_X = 1;
    if(Npt_X>1) Nsteps_X = Npt_X-1;
    int Nsteps_Y = 1;
    if(Npt_Y>1) Nsteps_Y = Npt_Y-1;
    
    //double step_X = (double)(fLastMinimVal_X - fFirstMinimVal_X) / (double)Nsteps_X;
    //double step_Y = (double)(fLastMinimVal_Y - fFirstMinimVal_Y) / (double)Nsteps_Y;
    double step_X = (double)(UniLastVal_X - UniFirstVal_X) / (double)Nsteps_X;
    double step_Y = (double)(UniLastVal_Y - UniFirstVal_Y) / (double)Nsteps_Y;
    
    for(int iStartPt_X=0;iStartPt_X<Npt_X;iStartPt_X++){
      for(int iStartPt_Y=0;iStartPt_Y<Npt_Y;iStartPt_Y++){
        
        //double currval0 = fFirstMinimVal_X + iStartPt_X*step_X;
        //double currval1 = fFirstMinimVal_Y + iStartPt_Y*step_Y;
        double currval0 = UniFirstVal_X + iStartPt_X*step_X;
        double currval1 = UniFirstVal_Y + iStartPt_Y*step_Y;
        
        int binX = conv::kBadResult;
        int binY = conv::kBadResult;
        
        // binX = fhchisq_marg_2D->GetXaxis()->FindBin(currval0+conv::kAVerySmallNum);
        // if(fIsFit2D) binY = fhchisq_marg_2D->GetYaxis()->FindBin(currval1+conv::kAVerySmallNum);
        binX = fhchisq_marg_2D->GetXaxis()->FindBin(currval0);
        if(fIsFit2D) binY = fhchisq_marg_2D->GetYaxis()->FindBin(currval1);
        //cout << binX << " " << binY << endl;
        
        bincont = conv::kBadResult;
        
        if(!fIsFit2D) bincont = CalcChisqMargPtRead(skobs_data,binX-1,fmin_bin_sub_param[1],readbinary);
        else bincont = CalcChisqMargPtRead(skobs_data,binX-1,binY-1,readbinary);
        
        if(bincont<ChisqMinStart){
          ChisqMinStart = bincont;
          fStartBin_X = binX;
          fStartBin_Y = binY;
        }
        
      } // for iStartPt_X
    } // for iStartPt_Y
  } // if fExtractStartPt
  
  else {
    cerr << endl;
    cerr << "FitMarginalization::GridStartBinRead" << endl;
    cerr << "The chosen strategy " << fExtractStartPt;
    cerr << " for the chisq minimization starting point  is not available" << endl;
    cerr << endl;
    util::Die(1);
  }
  
  // NEW PRECISION
  double startX = fhchisq_marg_2D->GetXaxis()->GetBinCenter(fStartBin_X);
  double startY = fhchisq_marg_2D->GetYaxis()->GetBinCenter(fStartBin_Y);
  //
  
  cout << "fStartBin_X = " << fStartBin_X << " - fStartBin_Y = " << fStartBin_Y << endl;
  cout << "X = " << startX << " - Y = " << startY << endl;
  cout << "Starting Chisq = " << ChisqMinStart << endl;
  cout << endl;
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetStartBin(double start_X,double start_Y){
  
  if(!fIsFit2D){
    //fStartBin_X = fhchisq_marg_1D->GetXaxis()->FindBin(start_X+conv::kAVerySmallNum);
    fStartBin_X = fhchisq_marg_1D->GetXaxis()->FindBin(start_X);
    fStartBin_Y = 1;
  }
  else{
    // fStartBin_X = fhchisq_marg_2D->GetXaxis()->FindBin(start_X+conv::kAVerySmallNum);
    // fStartBin_Y = fhchisq_marg_2D->GetYaxis()->FindBin(start_Y+conv::kAVerySmallNum);
    fStartBin_X = fhchisq_marg_2D->GetXaxis()->FindBin(start_X);
    fStartBin_Y = fhchisq_marg_2D->GetYaxis()->FindBin(start_Y);
  }
  
  bool GoodStartBin = true;
  
  if(!fIsFit2D){
    if(fStartBin_X < 1 || fStartBin_X > fhchisq_marg_1D->GetNbinsX() ){
      cerr << endl;
      cerr << "The starting X bin is not allowed" << endl;
      cerr << "  - X = " << start_X << " --> bin = " << fStartBin_X << endl;
      cerr << endl;
      GoodStartBin = false;
    }
  }
  else{
    if(fStartBin_X < 1 || fStartBin_X > fhchisq_marg_2D->GetNbinsX() ){
      cerr << endl;
      cerr << "The starting X bin is not allowed" << endl;
      cerr << "  - X = " << start_X << " --> bin = " << fStartBin_X << endl;
      cerr << endl;
      GoodStartBin = false;
    }
    if(fStartBin_Y < 1 || fStartBin_Y > fhchisq_marg_2D->GetNbinsY() ){
      cerr << endl;
      cerr << "The starting Y bin is not allowed" << endl;
      cerr << "  - Y = " << start_Y << " --> bin = " << fStartBin_Y << endl;
      cerr << endl;
      GoodStartBin = false;
    }
  }
  if(!GoodStartBin){
    cerr << endl;
    cerr << "FitMarginalization::SetStartBin" << endl;
    cerr << "The starting bin is not good" << endl;
    util::Die(1);
  }
  
  // cout << endl;
  // cout << "Starting point: " << endl;
  // cout << "  - X = " << start_X << " --> bin = " << fStartBin_X << endl;
  // cout << "  - Y = " << start_Y << " --> bin = " << fStartBin_Y << endl;
  // cout << endl;
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetMaxIterationSingleFit(int maxiter){
  fMaxIterationSingleFit = maxiter;
  cout << "Set maximum # of iterations per fit to " << fMaxIterationSingleFit << endl;
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetMaxFitNumber(int nfit){
  fMaxFitNumber = nfit;
  cout << "Set maximum # of fits to " << fMaxFitNumber << endl;
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetBorderDoFit(bool dofit){
  fBorderDoFit = dofit;
  cout << "If minimum is found on the border, redo the fit" << endl;
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetBorderMaxFitNumber(int nfit){
  fBorderMaxFitNumber = nfit;
  cout << "Set maximum # of fits after border minimum to " << fBorderMaxFitNumber << endl;
  return;
}
//________________________________________________________________________________________
double FitMarginalization::GetChisqMargTrueRead(double valpar0,double valpar1,SuperkObservation *skobs,bool readbinary){ // NEW PRECISION
  
  int bin0 = conv::kBadResult;
  int bin1 = conv::kBadResult;
  
  int nbinsX = fnbins_param[0];
  int nbinsY = fnbins_param[1];
  
  double min_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(1);
  double max_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(nbinsX);
  double min_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(1);
  double max_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(nbinsY);
  
  
  if(valpar0<min_X || valpar0>max_X){
    cerr << endl;
    cerr << "FitMarginalization::GetChisqMargTrueRead" << endl;
    cerr << "Parameter0 " << fnamefreeparam[0] << " is outside the total allowed range" << endl;
    cerr << "Allowed range: " << min_X << " - " << max_X << endl;
    cerr << "Param0 value: " << valpar0 << endl;
    cerr << endl;
    util::Die(1);
  }
  if( fIsFit2D && (valpar1<min_Y || valpar1>max_Y) ){
    cerr << endl;
    cerr << "FitMarginalization::GetChisqMargTrueRead" << endl;
    cerr << "Parameter1 " << fnamefreeparam[1] << " is outside the total allowed range" << endl;
    cerr << "Allowed range: " << min_Y << " - " << max_Y << endl;
    cerr << "Param1 value: " << valpar1 << endl;
    cerr << endl;
    util::Die(1);
  }
  
  if(!fIsFit2D){
    //bin0 = fhchisq_marg_1D->GetXaxis()->FindBin(valpar0+conv::kAVerySmallNum);
    bin0 = fhchisq_marg_1D->GetXaxis()->FindBin(valpar0);
  }
  else{
    // bin0 = fhchisq_marg_2D->GetXaxis()->FindBin(valpar0+conv::kAVerySmallNum);
    // bin1 = fhchisq_marg_2D->GetYaxis()->FindBin(valpar1+conv::kAVerySmallNum);
    bin0 = fhchisq_marg_2D->GetXaxis()->FindBin(valpar0);
    bin1 = fhchisq_marg_2D->GetYaxis()->FindBin(valpar1);
  }
  
  // Check if the true point is included by the sub-range
  // fmin_bin_sub_param is an index while bin* is the histogram bin number
  
  // Do this check only if binary files were not read at each grid point step before.
  // Otherwise the range is not the entire one but the one corresponding to the last
  // read binary file.
  if(!readbinary){
    if(
       (bin0 < (fmin_bin_sub_param[0]+1) || bin0 > (fmax_bin_sub_param[0]+1)) ||
       ( fIsFit2D && ( bin1<(fmin_bin_sub_param[1]+1) || bin1>(fmax_bin_sub_param[1]+1) ) )
       ){
      
      cout << endl;
      cout << "True point is not included in the sub-range" << endl;
      cout << "Set the value to " << conv::kALargeNum << endl;
      return conv::kALargeNum;
    }
  }
  
  double chisqtrue = conv::kBadResult;
  if(fIsFit2D) chisqtrue = CalcChisqMargPtRead(skobs,bin0-1,bin1-1,readbinary);
  else chisqtrue = CalcChisqMargPtRead(skobs,bin0-1,fmin_bin_sub_param[1],readbinary);
  
  return chisqtrue;
}
//________________________________________________________________________________________








// NEW SPECTRA MARG
double FitMarginalization::GetLikeToyRead(double valpar0,double valpar1,int iexperiment,bool readbinary){ // Use the map of the binary files
  
  if(this->fStoreLikeWeights){
    
    int bin0 = conv::kBadResult;
    int bin1 = conv::kBadResult;
    
    int nbinsX = fnbins_param[0];
    int nbinsY = fnbins_param[1];
    
    double min_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(1);
    double max_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(nbinsX);
    double min_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(1);
    double max_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(nbinsY);
    
    if(valpar0<min_X || valpar0>max_X){
      cerr << endl;
      cerr << "FitMarginalization::GetLikeToyRead" << endl;
      cerr << "Parameter0 " << fnamefreeparam[0] << " is outside the total allowed range" << endl;
      cerr << "Allowed range: " << min_X << " - " << max_X << endl;
      cerr << "Param0 value: " << valpar0 << endl;
      cerr << endl;
      util::Die(1);
    }
    
    if( fIsFit2D && (valpar1<min_Y || valpar1>max_Y) ){
      cerr << endl;
      cerr << "FitMarginalization::GetLikeToyRead" << endl;
      cerr << "Parameter1 " << fnamefreeparam[1] << " is outside the total allowed range" << endl;
      cerr << "Allowed range: " << min_Y << " - " << max_Y << endl;
      cerr << "Param1 value: " << valpar1 << endl;
      cerr << endl;
      util::Die(1);
    }
    
    if(!fIsFit2D){
      //bin0 = fhchisq_marg_1D->GetXaxis()->FindBin(valpar0+conv::kAVerySmallNum);
      bin0 = fhchisq_marg_1D->GetXaxis()->FindBin(valpar0);
    }
    else{
      // bin0 = fhchisq_marg_2D->GetXaxis()->FindBin(valpar0+conv::kAVerySmallNum);
      // bin1 = fhchisq_marg_2D->GetYaxis()->FindBin(valpar1+conv::kAVerySmallNum);
      bin0 = fhchisq_marg_2D->GetXaxis()->FindBin(valpar0);
      bin1 = fhchisq_marg_2D->GetYaxis()->FindBin(valpar1);
    }
    
    // Check if the true point is included by the sub-range
    // fmin_bin_sub_param is an index while bin* is the histogram bin number
    
    // Do this check only if binary files were not read at each grid point step before.
    // Otherwise the range is not the entire one but the one corresponding to the last
    // read binary file.
    if(!readbinary){
      if(
         (bin0 < (fmin_bin_sub_param[0]+1) || bin0 > (fmax_bin_sub_param[0]+1)) ||
         ( fIsFit2D && ( bin1<(fmin_bin_sub_param[1]+1) || bin1>(fmax_bin_sub_param[1]+1) ) )
         ){
        
        cout << endl;
        cout << "True point is not included in the sub-range" << endl;
        cout << "Set the value to " << conv::kALargeNum << endl;
        return conv::kALargeNum;
      }
    }
    
    double liketoy = conv::kBadResult;
    if(fIsFit2D) liketoy = flikeToyVArray->GetEntry(bin0-1,bin1-1,iexperiment);
    else         liketoy = flikeToyVArray->GetEntry(bin0-1,fmin_bin_sub_param[1],iexperiment);
    
    return liketoy;
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::WriteLikeToys2Binary" << endl;
    cerr << "Storing/Reading of toy likelihood weights is not enabled" << endl;
    cerr << endl;
    return conv::kBadResult;
  }
}
//










//________________________________________________________________________________________
void FitMarginalization::SetTagFileName(const char *tag){
  ftagfile = string(tag); // Get the file tag name
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetFileID(int nfiles){
  
  cout << endl;
  cout << "Build the ma[p between the binary file ID and the free parameters sub-range" << endl;
  cout << endl;
  
  char namemh[256];
  if(fdoMargMassHierarchy) sprintf(namemh,"%s","marg");
  else{
    if( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal )
    sprintf(namemh,"%s","NH");
    else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted )
    sprintf(namemh,"%s","IH");
    else{
      cerr << endl;
      cerr << "FitMarginalization::SetFileID" << endl;
      cerr << "Mass Hierarchy not defined" << endl;
      util::Die(1);
    }
  }
  
  for(int iptX=0;iptX<fnbins_param[0];iptX++){
    for(int iptY=0;iptY<fnbins_param[1];iptY++){
      fhmap_chisqrange_fileid_2D->SetBinContent(iptX+1,iptY+1,(int)conv::kBadResult);
    }
  }
  
  const char *tag = ftagfile.c_str(); // Get the file tag name
  
  // Store file names in a vector
  fnamefile = new string[nfiles];
  
  // Loop over the sub-range binary files
  
  for(int idfile=0;idfile<nfiles;idfile++){
    
    char inname[256];
    sprintf(inname,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"",".bin");
    
    std::string sinname(inname);
    fnamefile[idfile] = sinname;
    
#ifdef __print_build_map_fileid__
    cout << "Reading free param sub-range from binary file: " << fnamefile[idfile] << endl;
#endif
    
    
    
    
    // // Check # of binary files
    // //cout << fnamefilebinary << endl;
    // //cout << tag << endl;
    // char numfiletag[256];
    // char ext[256] = ".bin";
    // sprintf(numfiletag,"%s","./");
    // //sprintf(numfiletag,"../../inputs/marginalization/toys/nueapp_2D_NormPOT_fixosc_10toys/submit_sinsqth13Vsdcp_NormPOT_postfit*");
    // //TString snumfiletag(numfiletag);
    // //TIter next(files);
    // //util::Die(1);
    // cout << numfiletag << endl;
    // TSystemDirectory dir(numfiletag, numfiletag);
    // TList *files = dir.GetListOfFiles();
    // if (files) {
    //   TSystemFile *file;
    //   TString fname;
    //   TIter next(files);
    //   while ((file=(TSystemFile*)next())) {
    //  fname = file->GetName();
    //  if (!file->IsDirectory() &&
    //      fname.BeginsWith("../../inputs/marginalization/toys/nueapp_2D_NormPOT_fixosc_10toys/submit_sinsqth13Vsdcp_NormPOT_postfit") //&&
    //      //fname.EndsWith(".bin")
    //      ){
    //    cout << fname.Data() << endl;
    //  }
    //   }
    // }
    // util::Die(1);
    
    
    
    char inname_range[256];
    sprintf(inname_range,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"_subrange",".bin");
    std::string sinname_range(inname_range);
    
    int min_bin_sub_param0=0;
    int max_bin_sub_param0=0;
    int min_bin_sub_param1=0;
    int max_bin_sub_param1=0;
    
    ifstream is(sinname_range.c_str(),ios::binary|ios::in);
    
    // Check inputs
    
    if( !is.good() ){
      cerr << endl;
      cerr << "The input file " << sinname_range << " doesn't exist" << endl;
      cerr << "Check the input parameters: " << endl;
      cerr << "  - tag = " << tag << endl;
      cerr << "  - fntoys_marg = " << fntoys_marg << endl;
      cerr << "  - # of pt param 0 = " << fnbins_param[0] << endl;
      cerr << "  - # of pt param 1 = " << fnbins_param[1] << endl;
      //cerr << "  - MH = " << namemh << endl; // NEW MARG MH 2
      //cerr << "  - # of binary files = " << << endl;
      cerr << endl;
      util::Die(1);
    }
    
    assert(is);
    is.read( (char*)&min_bin_sub_param0,sizeof(int));
    is.read( (char*)&max_bin_sub_param0,sizeof(int));
    is.read( (char*)&min_bin_sub_param1,sizeof(int));
    is.read( (char*)&max_bin_sub_param1,sizeof(int));
    is.close();
    
#ifdef __print_build_map_fileid__
    cout << "min_bin_sub_0 = " << min_bin_sub_param0 << " - "
    << "max_bin_sub_0 = " << max_bin_sub_param0 << " - "
    << "min_bin_sub_1 = " << min_bin_sub_param1 << " - "
    << "max_bin_sub_1 = " << max_bin_sub_param1 << endl;
#endif
    
    for(int iptX=0;iptX<fnbins_param[0];iptX++){
      for(int iptY=0;iptY<fnbins_param[1];iptY++){
        if( iptX >= min_bin_sub_param0 &&
           iptX <= max_bin_sub_param0 &&
           iptY >= min_bin_sub_param1 &&
           iptY <= max_bin_sub_param1 ){
          fhmap_chisqrange_fileid_2D->SetBinContent(iptX+1,iptY+1,(int)idfile);
        }
      }
    }
    
  } // loop id file
  
  // Check if all the grid points are filled with the file ID
  
  for(int iptX=0;iptX<fnbins_param[0];iptX++){
    for(int iptY=0;iptY<fnbins_param[1];iptY++){
      double bincont = (double)fhmap_chisqrange_fileid_2D->GetBinContent(iptX+1,iptY+1);
      //cout << bincont << endl;
      if(bincont==conv::kBadResult){
        cerr << endl;
        cerr << "The grid point "
        << fhmap_chisqrange_fileid_2D->GetXaxis()->GetBinCenter(iptX+1) << " - "
        << fhmap_chisqrange_fileid_2D->GetYaxis()->GetBinCenter(iptY+1)
        << " is not filled with file ID" << endl;
        cout << "bincont = " << bincont << endl;
        cerr << "Check the input total # of binary files set to " << nfiles << endl;
        cerr << endl;
        util::Die(1);
      }
    }
  }
  
  return;
}
//________________________________________________________________________________________
int FitMarginalization::GetFileIDRead(int ipar0,int ipar1){
  int fileid=conv::kBadResult;
  fileid = (int)fhmap_chisqrange_fileid_2D->GetBinContent(ipar0+1,ipar1+1);
  return fileid;
}
//________________________________________________________________________________________
void FitMarginalization::CalcBestFitErrors(SuperkObservation *skobs,string method,bool readbinary){
  
  //double error = conv::kBadResult;
  
  if(method=="MINOS"){
    this->CalcBestFitErrorsMINOS(skobs,readbinary);
  }
  else if(method=="HESSE"){
    this->CalcBestFitErrorsHESSE(skobs,readbinary);
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::CalcBestFitErrors" << endl;
    cerr << "Method " << method << " is not available for fit error computation" << endl;
    cerr << endl;
    util::Die(1);
  }
  
  return;
}
//________________________________________________________________________________________
double FitMarginalization::GetInterpolationX(double X1,double Y1,double X2, double Y2,double Yin){
  double slope = (Y2-Y1)/(X2-X1);
  double Xout = (Yin-Y1)/slope + X1;
  return Xout;
  
  // cout << "Yin = " << Yin << endl;
  // cout << "X1 = " << X1 << " - X2 = " << X2 << endl;
  // cout << "Y1 = " << Y1 << " - Y2 = " << Y2 << endl;
  // cout << "Xout = " << Xout << endl;
  // cout << endl;
  // double Xin = Xout;
  // Yin = (Y2-Y1)/(X2-X1)*(Xin-X1) + Y1;
  // cout << "Yin = " << Yin << endl;
  //TGraph *gr = new TGraph();
  //gr->SetPoint(0,X1,Y1);
  //gr->SetPoint(1,X2,Y2);
  //cout << "Eval = " << gr->Eval(Xin) << endl;
}
//________________________________________________________________________________________
void FitMarginalization::CalcBestFitErrorsMINOS(SuperkObservation *skobs_data,bool readbinary){
  
  cout << "Compute MINOS errors" << endl;
  
  //double error = conv::kBadResult;
  const double CritChisq_1D = 1.;
  const double CritChisq_2D = 2.3;
  
  // Find 68.27% CL contour (1sigma)
  
  if(!fIsFit2D){   // 1D
    
    this->FillChisqMarg(skobs_data,readbinary);
    this->MinimizeHisto();
    
    // Get Delta Chisq distribution
    this->BuildDeltaChisqMarg();
    //fhdeltachisq_marg_1D = (TH1D*) fhchisq_marg_1D->Clone("fhdeltachisq_marg_1D"); // NEW PRECISION
    
    //int bin_bf = fhchisq_marg_1D->FindBin(this->GetParamBestFit(0)+conv::kAVerySmallNum); // best-fit bin
    int bin_bf = fhchisq_marg_1D->FindBin(this->GetParamBestFit(0)); // best-fit bin
    int bin_negsigma_a = conv::kBadResult;
    int bin_negsigma_b = conv::kBadResult;
    int bin_possigma_a = conv::kBadResult;
    int bin_possigma_b = conv::kBadResult;
    
    // Find param value at DeltaChisq=1 (negative error)
    
    for(int ibin=bin_bf;ibin>=0;ibin--){ // histogram convention
      
      double bincont = fhdeltachisq_marg_1D->GetBinContent(ibin);
      if(bincont<=0.){
        bincont = this->CalcChisqMargPtRead(skobs_data,ibin-1,fmin_bin_sub_param[1],readbinary);
        bincont -= this->GetChisqBestFit();
        fhdeltachisq_marg_1D->SetBinContent(ibin,bincont);
      }
      
      // cout << "ibin = " << ibin
      //           << " - X = " << fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(ibin)
      //           << " - bincont = " << bincont << endl;
      
      if(bincont>CritChisq_1D){
        // bins before (b) and after (a) DeltaChisq=1 (following x direction)
        bin_negsigma_b = ibin;
        bin_negsigma_a = bin_negsigma_b+1;
        break;
      }
      else if(bincont<0.){
        cerr << endl;
        cerr << "FitMarginalization::CalcBestFitErrorsMINOS" << endl;
        cerr << "Negative DeltaChisq 1D:" << endl;
        cerr << "DeltaChisq_1D = " << fhdeltachisq_marg_1D->GetBinContent(ibin) << endl;
        cerr << endl;
        util::Die(1);
      }
      
    }
    // cout << "bin_negsigma_b = " << bin_negsigma_b << endl;
    // cout << "bin_negsigma_a = " << bin_negsigma_a << endl;
    
    
    // Positive error
    
    // cout << endl;
    // cout << "Positive error" << endl;
    // cout << endl;
    
    for(int ibin=bin_bf;ibin<=fhdeltachisq_marg_1D->GetNbinsX();ibin++){ // histogram convention
      double bincont = fhdeltachisq_marg_1D->GetBinContent(ibin);
      if(bincont<=0.){
        bincont = this->CalcChisqMargPtRead(skobs_data,ibin-1,fmin_bin_sub_param[1],readbinary);
        bincont -= this->GetChisqBestFit();
        fhdeltachisq_marg_1D->SetBinContent(ibin,bincont);
      }
      
      // cout << "ibin = " << ibin
      //           << " - X = " << fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(ibin)
      //           << " - bincont = " << bincont << endl;
      
      if(bincont>CritChisq_1D){
        // bins before (b) and after (a) DeltaChisq=1
        bin_possigma_a = ibin;
        bin_possigma_b = bin_possigma_a-1;
        break;
      }
    }
    //cout << "bin_possigma_b = " << bin_possigma_b << endl;
    //cout << "bin_possigma_a = " << bin_possigma_a << endl;
    
    
    // If no bin is found it means the 1sigma bin is outside the physical region
    // and bin variables are still set to default value.
    // If so set left extreme to 1 (first bin of the histogram) and right extreme
    // to Nbins (last bin of the histogram)
    
    if(bin_negsigma_a == conv::kBadResult){
      bin_negsigma_a = 1;
      bin_negsigma_b = 1;
    }
    if(bin_possigma_b == conv::kBadResult){
      bin_possigma_a = fhdeltachisq_marg_1D->GetNbinsX();
      bin_possigma_b = fhdeltachisq_marg_1D->GetNbinsX();
    }
    
    // Check if the interval bins are correctly found
    if(bin_negsigma_a > bin_possigma_b){
      cerr << endl;
      cerr << "FitMarginalization::CalcBestFitErrorsMINOS" << endl;
      cerr << "Left extreme bin is higher than right extreme bin" << endl;
      cerr << "The 1sigma error interval is messed up:" << endl;
      cerr << "bin_possigma_b = " << bin_possigma_b << endl;
      cerr << "bin_negsigma_a = " << bin_negsigma_a << endl;
      cerr << endl;
      util::Die(1);
    }
    
    double X1 = conv::kBadResult;
    double Y1 = conv::kBadResult;
    double X2 = conv::kBadResult;
    double Y2 = conv::kBadResult;
    double Xout = conv::kBadResult;
    
    // Extrapolation (left extreme). Do it only if 1sigma bin is in the physical region
    if(bin_negsigma_a>1){
      X1 = fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(bin_negsigma_b);
      Y1 = fhdeltachisq_marg_1D->GetBinContent(bin_negsigma_b);
      X2 = fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(bin_negsigma_a);
      Y2 = fhdeltachisq_marg_1D->GetBinContent(bin_negsigma_a);
      Xout = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_1D);
    }
    else{
      Xout = fhdeltachisq_marg_1D->GetBinCenter(bin_negsigma_a);
    }
    fparam_bf_errneg[0] = fabs(this->GetParamBestFit(0) - Xout);
    fparam_bf_errneg[1] = conv::kBadResult;
    //cout << "fparam_bf_errneg[0] = " << fparam_bf_errneg[0] << endl;
    
    // Extrapolation (right extreme). Do it only if 1sigma bin is in the physical region
    if(bin_possigma_b < fhdeltachisq_marg_1D->GetNbinsX()){
      X1 = fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(bin_possigma_b);
      Y1 = fhdeltachisq_marg_1D->GetBinContent(bin_possigma_b);
      X2 = fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(bin_possigma_a);
      Y2 = fhdeltachisq_marg_1D->GetBinContent(bin_possigma_a);
      Xout = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_1D);
    }
    else{
      Xout = fhdeltachisq_marg_1D->GetBinCenter(bin_possigma_b);
    }
    fparam_bf_errpos[0] = fabs(this->GetParamBestFit(0) - Xout);
    fparam_bf_errpos[1] = conv::kBadResult;
    //cout << "fparam_bf_errpos[0] = " << fparam_bf_errpos[0] << endl;
    
  }
  
  // 2D fit
  
  else{
    
    // Get the full chisq distribution, since in 2D it's much more difficult to implement
    // an algorithm which stops at DeltaChisq = 2.3, since you need the full contour
    
    this->FillChisqMarg(skobs_data,true);  // read a new binary file at each step and fill the full range
    
    // If the global minimum was not found, replace the best-fit with the global one
    this->MinimizeHisto();
    
    // Get Delta Chisq distribution
    this->BuildDeltaChisqMarg();
    //fhdeltachisq_marg_1D = (TH1D*) fhchisq_marg_1D->Clone("fhdeltachisq_marg_1D"); // NEW PRECISION
    
    // int bin_bf0 = fhdeltachisq_marg_2D->GetXaxis()->FindBin(fparam_bf[0]+conv::kAVerySmallNum);
    // int bin_bf1 = fhdeltachisq_marg_2D->GetYaxis()->FindBin(fparam_bf[1]+conv::kAVerySmallNum);
    int bin_bf0 = fhdeltachisq_marg_2D->GetXaxis()->FindBin(fparam_bf[0]);
    int bin_bf1 = fhdeltachisq_marg_2D->GetYaxis()->FindBin(fparam_bf[1]);
    
    // X error: find the X interval with extremes the bins at X^2crit = 2.3
    
    int binleftX_x=conv::kBadResult;
    int binleftX_y=conv::kBadResult;
    
    int binrightX_x=conv::kBadResult;
    int binrightX_y=conv::kBadResult;
    
    int binleftY_x=conv::kBadResult;
    int binleftY_y=conv::kBadResult;
    
    int binrightY_x=conv::kBadResult;
    int binrightY_y=conv::kBadResult;
    
    for(int ibinX=0;ibinX<fhdeltachisq_marg_2D->GetNbinsX();ibinX++){
      for(int ibinY=0;ibinY<fhdeltachisq_marg_2D->GetNbinsY();ibinY++){
        
        int bincurrX = ibinX+1;
        int bincurrY = ibinY+1;
        double bincont = fhdeltachisq_marg_2D->GetBinContent(bincurrX,bincurrY);
        
        // Find the most left extreme bin of 1sigma interval along X
        double bincont_leftprevX = fhdeltachisq_marg_2D->GetBinContent(bincurrX+1,bincurrY); // must be on the right of bincurr
        if(bincont>=CritChisq_2D && bincont_leftprevX<=CritChisq_2D){
          if(binleftX_x == conv::kBadResult || (bincurrX) < binleftX_x ){
            binleftX_x = bincurrX;
            binleftX_y = bincurrY;
          }
        }
        
        // Find the most right extreme bin of 1sigma interval along X
        double bincont_rightprevX = fhdeltachisq_marg_2D->GetBinContent(bincurrX-1,bincurrY); // must be on the left of bincurr
        if(bincont>=CritChisq_2D && bincont_rightprevX<=CritChisq_2D){
          if(binrightX_x == conv::kBadResult || (bincurrX) > binrightX_x ){
            binrightX_x = bincurrX;
            binrightX_y = bincurrY;
          }
        }
        
        // Find the most left extreme bin of 1sigma interval along Y (down)
        double bincont_leftprevY = fhdeltachisq_marg_2D->GetBinContent(bincurrX,bincurrY+1); // must be on the right of bincurr
        if(bincont>=CritChisq_2D && bincont_leftprevY<=CritChisq_2D){
          if(binleftY_y == conv::kBadResult || (bincurrY) < binleftY_y ){
            binleftY_x = bincurrX;
            binleftY_y = bincurrY;
          }
        }
        
        // Find the most right extreme bin of 1sigma interval along Y
        double bincont_rightprevY = fhdeltachisq_marg_2D->GetBinContent(bincurrX,bincurrY-1); // must be on the left of bincurr
        if(bincont>=CritChisq_2D && bincont_rightprevY<=CritChisq_2D){
          if(binrightY_y == conv::kBadResult || (bincurrY) > binrightY_y ){
            binrightY_x = bincurrX;
            binrightY_y = bincurrY;
          }
        }
        
      } // for binY
    } // for binX
    
    
    //
    // Case when +-1sigma contour correspond to the full left range:
    // check if for some values all the full range is included within +-1sigma;
    // it could happen that for certain chisq shapes, for some point the condition
    // bincont>2.3 && bincont_prev<2.3 is satisfied, but it's a open contour
    // e.g.: sin2th23 Vs deltaCP, when DeltaChisq is low and for some values
    //       of sin2th23 all the deltaCP values are inside the +-1sigma interval
    // If there is not a +-1sigma contour in the considered total parameter range
    // set +-1sigma errors to the full range (left and right parameter extremes)
    //
    
    // Fix Y to a bin and loop over X. Left side.
    for(int ibinY=0;ibinY<fhdeltachisq_marg_2D->GetNbinsY();ibinY++){
      bool NoRangeErr_curr = true;
      for(int ibinX=0;ibinX<fhdeltachisq_marg_2D->GetNbinsX();ibinX++){
        int bincurrX = ibinX+1;
        int bincurrY = ibinY+1;
        double bincont = fhdeltachisq_marg_2D->GetBinContent(bincurrX,bincurrY);
        double bincont_leftprev = fhdeltachisq_marg_2D->GetBinContent(bincurrX+1,bincurrY); // must be on the right of bincurr
        if( (bincont>=CritChisq_2D && bincont_leftprev<=CritChisq_2D) ||
           (bincont>=CritChisq_2D) ){
          NoRangeErr_curr = false;
        }
      } // for Y
      if(NoRangeErr_curr){
        binleftX_x = 1;
        break;
      }
    } // for X
    
    // Fix Y to a bin and loop over X. Right side.
    for(int ibinY=0;ibinY<fhdeltachisq_marg_2D->GetNbinsY();ibinY++){
      bool NoRangeErr_curr = true;
      for(int ibinX=0;ibinX<fhdeltachisq_marg_2D->GetNbinsX();ibinX++){
        int bincurrX = ibinX+1;
        int bincurrY = ibinY+1;
        double bincont = fhdeltachisq_marg_2D->GetBinContent(bincurrX,bincurrY);
        double bincont_rightprev = fhdeltachisq_marg_2D->GetBinContent(bincurrX-1,bincurrY); // must be on the left of bincurr
        if( (bincont>=CritChisq_2D && bincont_rightprev<=CritChisq_2D) ||
           (bincont>=CritChisq_2D) ){
          NoRangeErr_curr = false;
        }
      } // for Y
      if(NoRangeErr_curr){
        binrightX_x = fhdeltachisq_marg_2D->GetNbinsX();
        break;
      }
    } // for X
    
    // Fix X to a bin and loop over Y. Left side (down).
    for(int ibinX=0;ibinX<fhdeltachisq_marg_2D->GetNbinsX();ibinX++){
      bool NoRangeErr_curr = true;
      for(int ibinY=0;ibinY<fhdeltachisq_marg_2D->GetNbinsY();ibinY++){
        int bincurrX = ibinX+1;
        int bincurrY = ibinY+1;
        double bincont = fhdeltachisq_marg_2D->GetBinContent(bincurrX,bincurrY);
        double bincont_leftprev = fhdeltachisq_marg_2D->GetBinContent(bincurrX,bincurrY+1); // must be on the right of bincurr
        if( (bincont>=CritChisq_2D && bincont_leftprev<=CritChisq_2D) ||
           (bincont>=CritChisq_2D) ){
          NoRangeErr_curr = false;
        }
      } // for Y
      if(NoRangeErr_curr){
        binleftY_y = 1;
        break;
      }
    } // for X
    
    // Fix X to a bin and loop over Y. Right side (up).
    for(int ibinX=0;ibinX<fhdeltachisq_marg_2D->GetNbinsX();ibinX++){
      bool NoRangeErr_curr = true;
      for(int ibinY=0;ibinY<fhdeltachisq_marg_2D->GetNbinsY();ibinY++){
        int bincurrX = ibinX+1;
        int bincurrY = ibinY+1;
        double bincont = fhdeltachisq_marg_2D->GetBinContent(bincurrX,bincurrY);
        double bincont_rightprev = fhdeltachisq_marg_2D->GetBinContent(bincurrX,bincurrY-1); // must be on the left of bincurr
        if( (bincont>=CritChisq_2D && bincont_rightprev<=CritChisq_2D) ||
           (bincont>=CritChisq_2D) ){
          NoRangeErr_curr = false;
        }
      } // for Y
      if(NoRangeErr_curr){
        binrightY_y = fhdeltachisq_marg_2D->GetNbinsY();
        break;
      }
    } // for X
    
    
    // If no bin is found it means the 1sigma bin is outside the physical region
    // and bin variables are still set to default value.
    // If so set left extreme to 1 (first bin of the histogram) and right extreme
    // to Nbins (last bin of the histogram)
    
    if(binleftX_x == conv::kBadResult){
      binleftX_x = 1;
      binleftX_y = bin_bf1;
    }
    if(binrightX_x == conv::kBadResult){
      binrightX_x = fhdeltachisq_marg_2D->GetNbinsX();
      binrightX_y = bin_bf1;
    }
    if(binleftY_y == conv::kBadResult){
      binleftY_x = bin_bf0;
      binleftY_y = 1;
    }
    if(binrightY_y == conv::kBadResult){
      binrightY_x = bin_bf0;
      binrightY_y = fhdeltachisq_marg_2D->GetNbinsY();
    }
    
    // Check if the interval bins are correctly found
    if(binleftX_x > binrightX_x){
      cerr << endl;
      cerr << "FitMarginalization::CalcBestFitErrorsMINOS" << endl;
      cerr << "X left extreme bin is higher than right extreme bin" << endl;
      cerr << "The 1sigma error interval is messed up:" << endl;
      cerr << "binleftX = " << binleftX_x << endl;
      cerr << "binrightX = " << binrightX_x << endl;
      cerr << endl;
      util::Die(1);
    }
    if(binleftY_y > binrightY_y){
      cerr << endl;
      cerr << "FitMarginalization::CalcBestFitErrorsMINOS" << endl;
      cerr << "Y left extreme bin is higher than right extreme bin" << endl;
      cerr << "The 1sigma error interval is messed up:" << endl;
      cerr << "binleftY = " << binleftY_y << endl;
      cerr << "binrightY = " << binrightY_y << endl;
      cerr << endl;
      util::Die(1);
    }
    
    //
    
    double valleftX = fhdeltachisq_marg_2D->GetXaxis()->GetBinCenter(binleftX_x);
    double valrightX = fhdeltachisq_marg_2D->GetXaxis()->GetBinCenter(binrightX_x);
    double valleftY = fhdeltachisq_marg_2D->GetYaxis()->GetBinCenter(binleftY_y);
    double valrightY = fhdeltachisq_marg_2D->GetYaxis()->GetBinCenter(binrightY_y);
    
    // cout << endl;
    // cout << "binleftX = " << binleftX_x << endl;
    // cout << "valleftX = " << valleftX << endl;
    // cout << "binrightX = " << binrightX_x << endl;
    // cout << "valrightX = " << valrightX << endl;
    // cout << endl;
    // cout << "binleftY = " << binleftY_y << endl;
    // cout << "valleftY = " << valleftY << endl;
    // cout << "binrightY = " << binrightY_y << endl;
    // cout << "valrightY = " << valrightY << endl;
    // cout << endl;
    
    // Interpolate between adjacent bins with bincont>2.3 and bincont<2.3
    // Do it independently for X and Y
    
    double X1 = conv::kBadResult;
    double Y1 = conv::kBadResult;
    double X2 = conv::kBadResult;
    double Y2 = conv::kBadResult;
    double Xout_left = conv::kBadResult;
    
    // X extrapolation (left extreme). Do it only if 1sigma bin is in the physical region
    if(binleftX_x>1){
      X1 = conv::kBadResult;
      Y1 = conv::kBadResult;
      X2 = conv::kBadResult;
      Y2 = conv::kBadResult;
      
      X1 = fhdeltachisq_marg_2D->GetXaxis()->GetBinCenter(binleftX_x); // left side
      Y1 = fhdeltachisq_marg_2D->GetBinContent(binleftX_x,binleftX_y); // left side
      X2 = fhdeltachisq_marg_2D->GetXaxis()->GetBinCenter(binleftX_x+1); // right side
      Y2 = fhdeltachisq_marg_2D->GetBinContent(binleftX_x+1,binleftX_y); // right side
      Xout_left = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_2D);
    }
    else{
      Xout_left = valleftX;
    }
    
    // cout << endl;
    // cout << "binleftX_x = " << binleftX_x << endl;
    // cout << "valleftX = " << valleftX << endl;
    // cout << "X1 = " << X1 << ", Y1 = " << Y1 << ", X2 = " << X2 << ", Y2 = " << Y2;
    // cout << " --> Xout_left = " << Xout_left << endl;
    
    // X extrapolation (right extreme). Do it only if 1sigma bin is in the physical region
    double Xout_right = conv::kBadResult;
    if(binrightX_x < fhdeltachisq_marg_2D->GetNbinsX()){
      X1 = conv::kBadResult;
      Y1 = conv::kBadResult;
      X2 = conv::kBadResult;
      Y2 = conv::kBadResult;
      
      X1 = fhdeltachisq_marg_2D->GetXaxis()->GetBinCenter(binrightX_x-1); // left side
      Y1 = fhdeltachisq_marg_2D->GetBinContent(binrightX_x-1,binrightX_y); // left side
      X2 = fhdeltachisq_marg_2D->GetXaxis()->GetBinCenter(binrightX_x); // right side
      Y2 = fhdeltachisq_marg_2D->GetBinContent(binrightX_x,binrightX_y); // right side
      Xout_right = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_2D);
    }
    else{
      Xout_right = valrightX;
    }
    
    // cout << endl;
    // cout << "binrightX_x = " << binrightX_x << endl;
    // cout << "valrightX = " << valrightX << endl;
    // cout << "X1 = " << X1 << ", Y1 = " << Y1 << ", X2 = " << X2 << ", Y2 = " << Y2;
    // cout << " --> Xout_right = " << Xout_right << endl;
    
    // Y extrapolation (left extreme, down). Do it only if 1sigma bin is in the physical region
    double Yout_left = conv::kBadResult;
    if(binleftY_y > 1){
      X1 = conv::kBadResult;
      Y1 = conv::kBadResult;
      X2 = conv::kBadResult;
      Y2 = conv::kBadResult;
      
      X1 = fhdeltachisq_marg_2D->GetYaxis()->GetBinCenter(binleftY_y); // left side
      Y1 = fhdeltachisq_marg_2D->GetBinContent(binleftY_x,binleftY_y); // left side
      X2 = fhdeltachisq_marg_2D->GetYaxis()->GetBinCenter(binleftY_y+1); // right side
      Y2 = fhdeltachisq_marg_2D->GetBinContent(binleftY_x,binleftY_y+1); // right side
      Yout_left = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_2D);
    }
    else{
      Yout_left = valleftY;
    }
    
    // cout << endl;
    // cout << "binleftY_y = " << binleftY_y << endl;
    // cout << "valleftY = " << valleftY << endl;
    // cout << "X1 = " << X1 << ", Y1 = " << Y1 << ", X2 = " << X2 << ", Y2 = " << Y2;
    // cout << " --> Yout_left = " << Yout_left << endl;
    
    // Y extrapolation (right extreme, up). Do it only if 1sigma bin is in the physical region
    double Yout_right = conv::kBadResult;
    if(binrightY_y < fhdeltachisq_marg_2D->GetNbinsY()){
      X1 = conv::kBadResult;
      Y1 = conv::kBadResult;
      X2 = conv::kBadResult;
      Y2 = conv::kBadResult;
      
      X1 = fhdeltachisq_marg_2D->GetYaxis()->GetBinCenter(binrightY_y-1); // left side
      Y1 = fhdeltachisq_marg_2D->GetBinContent(binrightY_x,binrightY_y-1); // left side
      X2 = fhdeltachisq_marg_2D->GetYaxis()->GetBinCenter(binrightY_y); // right side
      Y2 = fhdeltachisq_marg_2D->GetBinContent(binrightY_x,binrightY_y); // right side
      Yout_right = GetInterpolationX(X1,Y1,X2,Y2,CritChisq_2D);
    }
    else{
      Yout_right = valrightY;
    }
    
    // cout << endl;
    // cout << "binrightY_y = " << binrightY_y << endl;
    // cout << "valrightY = " << valrightY << endl;
    // cout << "X1 = " << X1 << ", Y1 = " << Y1 << ", X2 = " << X2 << ", Y2 = " << Y2;
    // cout << " --> Yout_right = " << Yout_right << endl;
    
    // Use the bin low edge as 1sigma extremes
    // fparam_bf_errneg[0] = fabs(this->GetParamBestFit(0) - valleftX);
    // fparam_bf_errpos[0] = fabs(this->GetParamBestFit(0) - valrightX);
    // fparam_bf_errneg[1] = fabs(this->GetParamBestFit(1) - valleftY);
    // fparam_bf_errpos[1] = fabs(this->GetParamBestFit(1) - valrightY);
    
    // Use the extrapolation between bin low edges of the adjacent bins
    // which include the 1sigma extreme
    fparam_bf_errneg[0] = fabs(this->GetParamBestFit(0) - Xout_left);
    fparam_bf_errpos[0] = fabs(this->GetParamBestFit(0) - Xout_right);
    fparam_bf_errneg[1] = fabs(this->GetParamBestFit(1) - Yout_left);
    fparam_bf_errpos[1] = fabs(this->GetParamBestFit(1) - Yout_right);
    
  } // if 2D fit
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::CalcBestFitErrorsHESSE(SuperkObservation * /*skobs_data*/,bool /*readbinary*/,int Nh)
{
  cout << "Compute HESSE errors" << endl;
#ifdef __print_hesseerror__
  cout << "Bins far " << Nh << " from the best-fit bin are used" << endl;
#endif
  
  //double error = conv::kBadResult;
  const double CritChisq_1D = 1.;
  const double CritChisq_2D = 2.3;
  
  
  //const int Nh = 3;
  
  this->BuildDeltaChisqMarg();
  
  // 1D HESSE error
  
  if(!fIsFit2D){
    
    //
    // Make a 1D Taylor expansion of the chisq function X^2 around the best fit.
    // Critical value for 1sigma (68.27% CL) is dX^2 = 1
    //
    // xerr = +- sqrt( 1. / A )
    // A = 0.5 * d^2 X^2 / dx^2
    //
    
    // 1D discrete deriv: d2X^2/dx^2 = [ X^2(xbf+h) - 2X^2(xbf) + X^2(xbf-h) ] / h^2
    
    int bin_bf = conv::kBadResult;
    double chisq_bf = conv::kBadResult;
    double chisq_minus = conv::kBadResult;
    double chisq_plus = conv::kBadResult;
    
    //bin_bf = fhdeltachisq_marg_1D->FindBin(fparam_bf[0]+conv::kAVerySmallNum);
    bin_bf = fhdeltachisq_marg_1D->FindBin(fparam_bf[0]);
    chisq_bf = fhdeltachisq_marg_1D->GetBinContent(bin_bf); // X^2(xbf)
    
    // If best-fit is on the border, assume the chisq to be symmetric around the minimum
    if(bin_bf==1){
      chisq_plus = fhdeltachisq_marg_1D->GetBinContent(bin_bf+Nh); // X^2(xbf+h)
      chisq_minus = chisq_plus;
    }
    else if(bin_bf==fhdeltachisq_marg_1D->GetNbinsX()){
      chisq_minus = fhdeltachisq_marg_1D->GetBinContent(bin_bf-Nh); // X^2(xbf-h)
      chisq_plus = chisq_minus;
    }
    else{
      chisq_minus = fhdeltachisq_marg_1D->GetBinContent(bin_bf-Nh); // X^2(xbf-h)
      chisq_plus = fhdeltachisq_marg_1D->GetBinContent(bin_bf+Nh); // X^2(xbf+h)
    }
    
    double h = fhdeltachisq_marg_1D->GetBinWidth(bin_bf) * (double)Nh ;
    double deriv2nd = (chisq_plus - 2*chisq_bf + chisq_minus) / ( h*h );
    
    // cout << endl;
    // cout << "chisq_bf = " << chisq_bf << endl;
    // cout << "chisq_minus = " << chisq_minus << endl;
    // cout << "chisq_plus = " << chisq_plus << endl;
    // cout << endl;
    
    // A = 0.5 * (1D discrete derivative)
    double A = 0.5 * deriv2nd;
    
    double xerr = sqrt( CritChisq_1D / A );
    
    
    // // Check that best-fit errors are in the physical region.
    // // If not restrict them to the physical region
    // double x_p1sigma = fparam_bf[0] + xerr;
    // double x_m1sigma = fparam_bf[0] - xerr;
    // x_p1sigma = this->OscPhysicalRegion(fPhysParamFree[0],x_p1sigma);
    // x_m1sigma = this->OscPhysicalRegion(fPhysParamFree[0],x_m1sigma);
    // double xerr_p1sigma = fabs(x_p1sigma - fparam_bf[0]);
    // double xerr_m1sigma = fabs(x_m1sigma - fparam_bf[0]);
    // xerr = TMath::Max(xerr_p1sigma,xerr_m1sigma);
    
    fparam_bf_errsym[0] = xerr; // symmetric error
    fparam_bf_errsym[1] = conv::kBadResult;
  }
  
  // 2D HESSE error
  
  else{
    
#ifdef __print_hesseerror__
    cout << "2D HESSE error" << endl;
#endif
    
    
    
    // Be sure the Chisq is computed for all the points are needed for the HESSE error
    // TODO: I don't need all the grid for HESSE errors.
    // for(int i=1;i<=Nh;i++){
    //   for(int j=1;j<=Nh;j++){
    //  if( (bin_bf0+i <= fhdeltachisq_marg_2D->GetNbinsX()) &&
    //      (bin_bf1+j <= fhdeltachisq_marg_2D->GetNbinsY()) )
    //    CalcChisqMargPtRead(skobs_data,bin_bf0+i-1,bin_bf1+j-1,true);
    //  if( (bin_bf0-i >= 1) &&
    //      (bin_bf1-j >= 1) )
    //    CalcChisqMargPtRead(skobs_data,bin_bf0-i-1,bin_bf1-j-1,true);
    //  if( (bin_bf0-i >= 1) &&
    //      (bin_bf1+j <= fhdeltachisq_marg_2D->GetNbinsY()) )
    //  CalcChisqMargPtRead(skobs_data,bin_bf0-i-1,bin_bf1+j-1,true);
    //  if( (bin_bf0+i <= fhdeltachisq_marg_2D->GetNbinsX()) &&
    //      (bin_bf1-j >= 1) )
    //  CalcChisqMargPtRead(skobs_data,bin_bf0+i-1,bin_bf1-j-1,true);
    //   }
    // }
    // // Be sure the +-Nh bins are in the physical region
    // if( (bin_bf0-Nh) < 1 ||
    //  (bin_bf0+Nh) > fhdeltachisq_marg_2D->GetNbinsX() ||
    //  (bin_bf1-Nh) < 1 ||
    //  (bin_bf1+Nh) > fhdeltachisq_marg_2D->GetNbinsY()
    //  )
    //   {
    //   }
    
    
    
    //
    // Make a 2D Taylor expansion of the chisq function X^2 around the best fit.
    // Critical value for 1sigma (68.27% CL) is dX^2 = 2.3
    // Discrete derivatives at the best-fit are used
    //
    // x0err = +- sqrt( 2.3 / (A-C^2/B) )
    // x1err = +- sqrt( 2.3 / (B-C^2/A) )
    // A = 0.5 * d^2 X^2 / dx0^2
    // B = 0.5 * d^2 X^2 / dx1^2
    // 2C = d^2 X^2 / dx0dx1
    //
    
    int bin_bf0 = fhdeltachisq_marg_2D->GetXaxis()->FindBin(fparam_bf[0]);
    int bin_bf1 = fhdeltachisq_marg_2D->GetYaxis()->FindBin(fparam_bf[1]);
    
    double chisq_bf = fhdeltachisq_marg_2D->GetBinContent(bin_bf0,bin_bf1);     // X^2( x0bf, x1bf )
    double h0 = fhdeltachisq_marg_2D->GetXaxis()->GetBinWidth(bin_bf0) * (double)Nh;
    double h1 = fhdeltachisq_marg_2D->GetYaxis()->GetBinWidth(bin_bf1) * (double)Nh;
    
    double chisq_minus0 = conv::kBadResult;
    double chisq_plus0 = conv::kBadResult;
    double chisq_minus1 = conv::kBadResult;
    double chisq_plus1 = conv::kBadResult;
    double chisq_minus01 = conv::kBadResult;
    double chisq_plus01 = conv::kBadResult;
    double chisq_plus0_minus1 = conv::kBadResult;
    double chisq_minus0_plus1 = conv::kBadResult;
    
    double deriv2nd_0 = conv::kBadResult;
    double deriv2nd_1 = conv::kBadResult;
    double deriv2nd_01 = conv::kBadResult;
    
    //
    // Treat cases when the best-fit is on the physical boundary.
    // Assume the chisq to be symmetric around the minimum
    // --> symmetrically copy and paste the chisq values to the bins in the unphysical region
    //
    
    if( bin_bf0==1 &&  // Only Par0 on the lowery boundary
       (bin_bf1>1 && bin_bf1<fhdeltachisq_marg_2D->GetNbinsY()) ){
      
      chisq_plus0_minus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0+Nh, bin_bf1-Nh);  // X^2( x0bf+h0, x1bf-h1 )
      chisq_minus0_plus1 = chisq_plus0_minus1;  // symmetric                               // X^2( x0bf-h0, x1bf+h1 )
      chisq_plus0  = fhdeltachisq_marg_2D->GetBinContent   (bin_bf0+Nh, bin_bf1);          // X^2( x0bf+h0, x1bf    )
      chisq_minus0 = chisq_plus0;               // symmetric                               // X^2( x0bf-h0, x1bf    )
      chisq_plus01  = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0+Nh, bin_bf1+Nh);       // X^2( x0bf+h0, x1bf+h1 )
      chisq_minus01 = chisq_plus01;             // symmetric                               // X^2( x0bf-h0, x1bf-h1 )
      chisq_minus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0,    bin_bf1-Nh);        // X^2( x0bf   , x1bf-h1 )
      chisq_plus1 = fhdeltachisq_marg_2D->GetBinContent   (bin_bf0,    bin_bf1+Nh);        // X^2( x0bf   , x1bf+h1 )
      
      // All the bins are available
      // 2D discrete deriv 01 (2): d2X^2/dx0dx1 = [ X^2(x0bf+h0,x1bf+h1) - X^2(x0bf+h0,x0bf-h1) - X^2(x0bf-h0,x1bf+h1) + X^2(x0bf-h0,x1bf-h1) ] / (4*h0*h1)
      deriv2nd_01 = ( chisq_plus01 - chisq_plus0_minus1 - chisq_minus0_plus1 + chisq_minus01) / (4*h0*h1);
      
    }
    else if( bin_bf0==fhdeltachisq_marg_2D->GetNbinsX() && // Only Par0 on the upper boundary
            (bin_bf1>1 && bin_bf1<fhdeltachisq_marg_2D->GetNbinsY()) ){
      
      chisq_minus01 = fhdeltachisq_marg_2D->GetBinContent (bin_bf0-Nh, bin_bf1-Nh);        // X^2( x0bf-h0, x1bf-h1 )
      chisq_plus01 = chisq_minus01;            // symmetric                                // X^2( x0bf+h0, x1bf+h1 )
      chisq_minus0 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0-Nh, bin_bf1);           // X^2( x0bf-h0, x1bf    )
      chisq_plus0 = chisq_minus0;              // symmetric                                // X^2( x0bf+h0, x1bf    )
      chisq_minus0_plus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0-Nh, bin_bf1+Nh);  // X^2( x0bf-h0, x1bf+h1 )
      chisq_plus0_minus1 = chisq_minus0_plus1; // symmetric                                // X^2( x0bf+h0, x1bf-h1 )
      chisq_minus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0,    bin_bf1-Nh);        // X^2( x0bf   , x1bf-h1 )
      chisq_plus1 = fhdeltachisq_marg_2D->GetBinContent   (bin_bf0,    bin_bf1+Nh);        // X^2( x0bf   , x1bf+h1 )
      
      // All the bins are available
      // 2D discrete deriv 01 (2): d2X^2/dx0dx1 = [ X^2(x0bf+h0,x1bf+h1) - X^2(x0bf+h0,x0bf-h1) - X^2(x0bf-h0,x1bf+h1) + X^2(x0bf-h0,x1bf-h1) ] / (4*h0*h1)
      deriv2nd_01 = ( chisq_plus01 - chisq_plus0_minus1 - chisq_minus0_plus1 + chisq_minus01) / (4*h0*h1);
      
    }
    else if( (bin_bf0>1 && bin_bf0<fhdeltachisq_marg_2D->GetNbinsX()) && // Only Par1 on the lower boundary
            bin_bf1==1 ){
      chisq_plus01 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0+Nh, bin_bf1+Nh);        // X^2( x0bf+h0, x1bf+h1 )
      chisq_minus01 = chisq_plus01;            // symmetric                                // X^2( x0bf-h0, x1bf-h1 )
      chisq_plus1 = fhdeltachisq_marg_2D->GetBinContent   (bin_bf0, bin_bf1+Nh);           // X^2( x0bf, x1bf+h0    )
      chisq_minus1 = chisq_plus1;              // symmetric                                // X^2( x0bf, x1bf-h0    )
      chisq_minus0_plus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0-Nh, bin_bf1+Nh);  // X^2( x0bf-h0, x1bf+h1 )
      chisq_plus0_minus1 = chisq_minus0_plus1; // symmetric                                // X^2( x0bf+h0, x1bf-h1 )
      chisq_minus0 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0-Nh,    bin_bf1);        // X^2( x0bf-h0, x1bf    )
      chisq_plus0 = fhdeltachisq_marg_2D->GetBinContent   (bin_bf0+Nh,    bin_bf1);        // X^2( x0bf+h0, x1bf    )
      
      // All the bins are available
      // 2D discrete deriv 01 (2): d2X^2/dx0dx1 = [ X^2(x0bf+h0,x1bf+h1) - X^2(x0bf+h0,x0bf-h1) - X^2(x0bf-h0,x1bf+h1) + X^2(x0bf-h0,x1bf-h1) ] / (4*h0*h1)
      deriv2nd_01 = ( chisq_plus01 - chisq_plus0_minus1 - chisq_minus0_plus1 + chisq_minus01) / (4*h0*h1);
      deriv2nd_01 = ( chisq_plus01 - chisq_plus0 - chisq_plus1 + 2*chisq_bf - chisq_minus0 - chisq_minus1 + chisq_minus01 ) / (2*h0*h1);
    }
    else if( (bin_bf0>1 && bin_bf0<fhdeltachisq_marg_2D->GetNbinsX()) && // Only Par1 on the upper boundary
            bin_bf1==fhdeltachisq_marg_2D->GetNbinsY() ){
      chisq_minus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0,    bin_bf1-Nh);        // X^2( x0bf   , x1bf-h1 )
      chisq_plus1 = chisq_minus1;              // symmetric                                // X^2( x0bf   , x1bf+h1 )
      chisq_plus0_minus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0+Nh, bin_bf1-Nh);  // X^2( x0bf+h0, x1bf-h1 )
      chisq_minus0_plus1 = chisq_plus0_minus1; // symmetric                                // X^2( x0bf-h0, x1bf+h1 )
      chisq_minus01 = fhdeltachisq_marg_2D->GetBinContent (bin_bf0-Nh, bin_bf1-Nh);        // X^2( x0bf-h0, x1bf-h1 )
      chisq_plus01 = chisq_minus01;            // symmetric                                // X^2( x0bf+h0, x1bf+h1 )
      chisq_minus0 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0-Nh, bin_bf1);           // X^2( x0bf-h0, x1bf    )
      chisq_plus0 = fhdeltachisq_marg_2D->GetBinContent   (bin_bf0+Nh, bin_bf1);           // X^2( x0bf+h0, x1bf    )
      
      // All the bins are available
      // 2D discrete deriv 01 (2): d2X^2/dx0dx1 = [ X^2(x0bf+h0,x1bf+h1) - X^2(x0bf+h0,x0bf-h1) - X^2(x0bf-h0,x1bf+h1) + X^2(x0bf-h0,x1bf-h1) ] / (4*h0*h1)
      deriv2nd_01 = ( chisq_plus01 - chisq_plus0_minus1 - chisq_minus0_plus1 + chisq_minus01) / (4*h0*h1);
      
    }
    else if(bin_bf0==1 && bin_bf1==1){ // Par0 and Par1 on the lowery boundary
      chisq_plus01 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0+Nh, bin_bf1+Nh);        // X^2( x0bf+h0, x1bf+h1 )
      chisq_minus01 = chisq_plus01; // symmetric                                           // X^2( x0bf-h0, x1bf-h1 )
      chisq_plus0 = fhdeltachisq_marg_2D->GetBinContent   (bin_bf0+Nh, bin_bf1);           // X^2( x0bf+h0, x1bf    )
      chisq_minus0 = chisq_plus0;   // symmetric                                           // X^2( x0bf-h0, x1bf    )
      chisq_plus1 = fhdeltachisq_marg_2D->GetBinContent   (bin_bf0,    bin_bf1+Nh);        // X^2( x0bf   , x1bf+h1 )
      chisq_minus1 = chisq_plus1;   // symmetric                                           // X^2( x0bf   , x1bf-h1 )
      
      // Set to conv::kBadResult bins when also its symmetric is in the unphysical region
      chisq_plus0_minus1 = conv::kBadResult;  // X^2( x0bf+h0, x1bf-h1 )
      chisq_minus0_plus1 = conv::kBadResult;  // X^2( x0bf-h0, x1bf+h1 )
      
      // Both X^2( x0bf+h0, x0bf-h1 ) and X^2( x0bf-h0, x0bf+h1 ), which are symmetric,
      // are in the unphysical region and cannot be used
      // --> use alternative function to compute the second derivative d2X^2/dx0dx1
      // 2D discrete deriv 01 (1A): d2X^2/dx0dx1 = [ X^2(x0bf+h0,x1bf+h1) - X^2(x0bf+h0,x0bf) - X^2(x0bf,x1bf+h1) + 2*X^2(x0bf,x1bf) - X^2(x0bf-h0,x1bf)
      //                                      - X^2(x0bf,x1bf-h1) + X^2(x0bf-h0,x1bf-h1)] / (2*h0*h1)
      deriv2nd_01 = ( chisq_plus01 - chisq_plus0 - chisq_plus1 + 2*chisq_bf - chisq_minus0 - chisq_minus1 + chisq_minus01 ) / (2*h0*h1);
      
    }
    else if(bin_bf0==1 && bin_bf1==fhdeltachisq_marg_2D->GetNbinsY()){ // Par0 on the lowery boundary and Par1 on the upper boundary
      chisq_plus0_minus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0+Nh, bin_bf1-Nh);  // X^2( x0bf+h0, x1bf-h1 )
      chisq_minus0_plus1 = chisq_plus0_minus1; // symmetric                                // X^2( x0bf-h0, x1bf+h1 )
      chisq_plus0 = fhdeltachisq_marg_2D->GetBinContent   (bin_bf0+Nh, bin_bf1);           // X^2( x0bf+h0, x1bf    )
      chisq_minus0 = chisq_plus0;              // symmetric                                // X^2( x0bf-h0, x1bf    )
      chisq_minus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0,    bin_bf1-Nh);        // X^2( x0bf   , x1bf-h1 )
      chisq_plus1 = chisq_minus1;              // symmetric                                // X^2( x0bf   , x1bf+h1 )
      
      // Set to conv::kBadResult bins when also its symmetric is in the unphysical region
      chisq_minus01 = conv::kBadResult;  // X^2( x0bf-h0, x1bf-h1 )
      chisq_plus01 = conv::kBadResult;   // X^2( x0bf+h0, x1bf+h1 )
      
      // Both X^2( x0bf-h0, x0bf-h1 ) and X^2( x0bf+h0, x0bf+h1 ), which are symmetric,
      // are in the unphysical region and cannot be used
      // --> use alternative function to compute the second derivative d2X^2/dx0dx1
      // 2D discrete deriv 01 (1B): d2X^2/dx0dx1 = [ X^2(x0bf-h0,x1bf+h1) - X^2(x0bf,x1bf+h1) - X^2(x0bf-h0,x1bf) + 2*X^2(x0bf,x1bf) - X^2(x0bf,x1bf-h1)
      //                                      - X^2(x0bf+h0,x1bf) + X^2(x0bf+h0,x1bf-h1)] / (2*h0*h1)
      deriv2nd_01 = ( chisq_minus0_plus1 - chisq_plus1 - chisq_minus0 + 2*chisq_bf - chisq_minus1 - chisq_plus0 + chisq_plus0_minus1 ) / (2*h0*h1);
      
    }
    else if(bin_bf0==fhdeltachisq_marg_2D->GetNbinsX() && bin_bf1==fhdeltachisq_marg_2D->GetNbinsY()){ // Par0 and Par1 on the upper boundary
      chisq_minus01 = fhdeltachisq_marg_2D->GetBinContent (bin_bf0-Nh, bin_bf1-Nh);        // X^2( x0bf-h0, x1bf-h1 )
      chisq_plus01 = chisq_minus01; // symmetric                                           // X^2( x0bf+h0, x1bf+h1 )
      chisq_minus0 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0-Nh, bin_bf1);           // X^2( x0bf-h0, x1bf    )
      chisq_plus0 = chisq_minus0;   // symmetric                                           // X^2( x0bf+h0,
      chisq_minus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0,    bin_bf1-Nh);        // X^2( x0bf   , x1bf-h1 )
      chisq_plus1 = chisq_minus1;   // symmetric                                           // X^2( x0bf   , x1bf+h1 )
      
      // Set to conv::kBadResult bins when also its symmetric is in the unphysical region
      chisq_plus0_minus1 = conv::kBadResult;  // X^2( x0bf+h0, x1bf-h1 )
      chisq_minus0_plus1 = conv::kBadResult;  // X^2( x0bf-h0, x1bf+h1 )
      
      // Both X^2( x0bf+h0, x0bf-h1 ) and X^2( x0bf-h0, x0bf+h1 ), which are symmetric,
      // are in the unphysical region and cannot be used
      // --> use alternative function to compute the second derivative d2X^2/dx0dx1
      // 2D discrete deriv 01 (1A): d2X^2/dx0dx1 = [ X^2(x0bf+h0,x1bf+h1) - X^2(x0bf+h0,x0bf) - X^2(x0bf,x1bf+h1) + 2*X^2(x0bf,x1bf) - X^2(x0bf-h0,x1bf)
      //                                      - X^2(x0bf,x1bf-h1) + X^2(x0bf-h0,x1bf-h1)] / (2*h0*h1)
      deriv2nd_01 = ( chisq_plus01 - chisq_plus0 - chisq_plus1 + 2*chisq_bf - chisq_minus0 - chisq_minus1 + chisq_minus01 ) / (2*h0*h1);
      
    }
    else if(bin_bf0==fhdeltachisq_marg_2D->GetNbinsX() && bin_bf1==1){ // Par0 on the upper boundary and Par1 on the lower boundary
      chisq_minus0_plus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0-Nh, bin_bf1+Nh);  // X^2( x0bf-h0, x1bf+h1 )
      chisq_plus0_minus1 = chisq_minus0_plus1; // symmetric                                // X^2( x0bf+h0, x1bf-h1 )
      chisq_plus1 = fhdeltachisq_marg_2D->GetBinContent   (bin_bf0,    bin_bf1+Nh);        // X^2( x0bf   , x1bf+h1 )
      chisq_minus1 = chisq_plus1;              // symmetric                                // X^2( x0bf   , x1bf-h1 )
      chisq_minus0 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0-Nh, bin_bf1);           // X^2( x0bf-h0, x1bf    )
      chisq_plus0 = chisq_minus0;              // symmetric                                // X^2( x0bf+h0, x1bf    )
      
      // Set to conv::kBadResult bins when also its symmetric is in the unphysical region
      chisq_minus01 = conv::kBadResult;        // X^2( x0bf-h0, x1bf-h1 )
      chisq_plus01 = conv::kBadResult;        // X^2( x0bf+h0, x1bf+h1 )
      
      // Both X^2( x0bf-h0, x0bf-h1 ) and X^2( x0bf+h0, x0bf+h1 ), which are symmetric,
      // are in the unphysical region and cannot be used
      // --> use alternative function to compute the second derivative d2X^2/dx0dx1
      // 2D discrete deriv 01 (1B): d2X^2/dx0dx1 = [ X^2(x0bf-h0,x1bf+h1) - X^2(x0bf,x1bf+h1) - X^2(x0bf-h0,x1bf) + 2*X^2(x0bf,x1bf) - X^2(x0bf,x1bf-h1)
      //                                      - X^2(x0bf+h0,x1bf) + X^2(x0bf+h0,x1bf-h1)] / (2*h0*h1)
      deriv2nd_01 = ( chisq_minus0_plus1 - chisq_plus1 - chisq_minus0 + 2*chisq_bf - chisq_minus1 - chisq_plus0 + chisq_plus0_minus1 ) / (2*h0*h1);
      
    }
    else if( (bin_bf0>1 && bin_bf0<fhdeltachisq_marg_2D->GetNbinsX()) &&  // Standard case
            (bin_bf1>1 && bin_bf1<fhdeltachisq_marg_2D->GetNbinsY()) ){
      
      chisq_minus0 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0-Nh, bin_bf1);           // X^2( x0bf-h0, x1bf    )
      chisq_plus0 = fhdeltachisq_marg_2D->GetBinContent   (bin_bf0+Nh, bin_bf1);           // X^2( x0bf+h0, x1bf    )
      chisq_minus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0,    bin_bf1-Nh);        // X^2( x0bf   , x1bf-h1 )
      chisq_plus1 = fhdeltachisq_marg_2D->GetBinContent   (bin_bf0,    bin_bf1+Nh);        // X^2( x0bf   , x1bf+h1 )
      chisq_minus01 = fhdeltachisq_marg_2D->GetBinContent (bin_bf0-Nh, bin_bf1-Nh);        // X^2( x0bf-h0, x1bf-h1 )
      chisq_plus01 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0+Nh, bin_bf1+Nh);        // X^2( x0bf+h0, x1bf+h1 )
      chisq_plus0_minus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0+Nh, bin_bf1-Nh);  // X^2( x0bf+h0, x1bf-h1 )
      chisq_minus0_plus1 = fhdeltachisq_marg_2D->GetBinContent  (bin_bf0-Nh, bin_bf1+Nh);  // X^2( x0bf-h0, x1bf+h1 )
      
#ifdef __print_hesseerror_verbose__
      cout << endl;
      cout << "chisq_bf = " << chisq_bf << endl;
      cout << "chisq_minus0 = " << chisq_minus0 << endl;
      cout << "chisq_plus0 = " << chisq_plus0 << endl;
      cout << "chisq_minus1 = " << chisq_minus1 << endl;
      cout << "chisq_plus1 = " << chisq_plus1 << endl;
      cout << "chisq_minus01 = " << chisq_minus01 << endl;
      cout << "chisq_plus01 = " << chisq_plus01 << endl;
      cout << "chisq_plus0_minus1 = " << chisq_plus0_minus1 << endl;
      cout << "chisq_minus0_plus1 = " << chisq_minus0_plus1 << endl;
      cout << endl;
#endif
      
      // All the bins are available
      // 2D discrete deriv 01 (2): d2X^2/dx0dx1 = [ X^2(x0bf+h0,x1bf+h1) - X^2(x0bf+h0,x0bf-h1) - X^2(x0bf-h0,x1bf+h1) + X^2(x0bf-h0,x1bf-h1) ] / (4*h0*h1)
      deriv2nd_01 = ( chisq_plus01 - chisq_plus0_minus1 - chisq_minus0_plus1 + chisq_minus01) / (4*h0*h1);
    }
    else{
      cerr << endl;
      cerr << "FitMarginalization::CalcBestFitErrorsHESSE" << endl;
      cerr << "The best-fit in the physical boundary is not treated correctly" << endl;
      cerr << "bin_bf0 = " << bin_bf0 << endl;
      cerr << "bin_bf1 = " << bin_bf1 << endl;
      cerr << endl;
      util::Die(1);
    }
    
    // 2D discrete deriv 0: d2X^2/dx0^2 = [ X^2(x0bf+h0,x1bf) - 2X^2(x0bf,x1bf) + X^2(x0bf-h0,x1bf) ] / h0^2
    deriv2nd_0 = (chisq_plus0 - 2*chisq_bf + chisq_minus0) / ( h0*h0 );
    
    // 2D discrete deriv 1: d2X^2/dx1^2 = [ X^2(x0bf,x1bf+h1) - 2X^2(x0bf,x1bf) + X^2(x0bf,x1bf-h1) ] / h1^2
    deriv2nd_1 = (chisq_plus1 - 2*chisq_bf + chisq_minus1) / ( h1*h1 );
    
    //
    // Compute the intermediate coefficients:
    // - A = 0.5 * deriv2nd_0
    // - B = 0.5 * deriv2nd_1
    // - 2C = deriv2nd_01
    //
    
    double A = 0.5 * deriv2nd_0;
    double B = 0.5 * deriv2nd_1;
    double C = 0.5 * deriv2nd_01;
    
    double sqrtarg0 = A - C*C/B;
    double sqrtarg1 = B - C*C/A;
    
    // Check that the derivatives make sense, otherwise
    // redo iteratively again HESSE error but with points
    // farther from  the minimum
    
    
    if( sqrtarg0<=0. || sqrtarg1<=0. ){
      
#ifdef __print_hesseerror__
      cout << "Bin number wrt to best fit is " << Nh << endl;
      cout << "The second derivatives are not good:" << endl;
      cout << " - deriv2nd_0 = " << deriv2nd_0 << endl;
      cout << " - deriv2nd_1 = " << deriv2nd_1 << endl;
      cout << " - deriv2nd_01 = " << deriv2nd_01 << endl;
      cout << " - A = " << A << endl;
      cout << " - B = " << B << endl;
      cout << " - C = " << C << endl;
      cout << " - A - C*C/B = " << A - C*C/B << endl;
      cout << " - B - C*C/A = " << B - C*C/A << endl;
#endif
      
      //const int newnbinaround = Nh+5;
      //cout << "Recompute the HESSE error using bin number " << newnbinaround << endl;
      //cout << endl;
      //return CalcBestFitErrorsHESSE(skobs_data,readbinary,newnbinaround);
      
      cout << "Recompute the HESSE error treating free parameters as independent 1D fits" << endl;
      double deriv2nd_x = (chisq_plus0 - 2*chisq_bf + chisq_minus0) / ( h0*h0 );
      sqrtarg0 = 0.5 * deriv2nd_x;
      double deriv2nd_y = (chisq_plus1 - 2*chisq_bf + chisq_minus1) / ( h1*h1 );
      sqrtarg1 = 0.5 * deriv2nd_y;
    }
    
    // Get best-fit errors
    
    double bf0err = sqrt( CritChisq_2D ) / sqrt( sqrtarg0 );
    double bf1err = sqrt( CritChisq_2D ) / sqrt( sqrtarg1 );
    
    // // Check that best-fit errors are in the physical region.
    // // If not restrict them to the physical region
    // double bf0_p1sigma = fparam_bf[0] + bf0err;
    // double bf0_m1sigma = fparam_bf[0] - bf0err;
    // double bf1_p1sigma = fparam_bf[1] + bf1err;
    // double bf1_m1sigma = fparam_bf[1] - bf1err;
    // bf0_p1sigma = this->OscPhysicalRegion(fPhysParamFree[0],bf0_p1sigma);
    // bf0_m1sigma = this->OscPhysicalRegion(fPhysParamFree[0],bf0_m1sigma);
    // bf1_p1sigma = this->OscPhysicalRegion(fPhysParamFree[1],bf1_p1sigma);
    // bf1_m1sigma = this->OscPhysicalRegion(fPhysParamFree[1],bf1_m1sigma);
    // double bf0err_p1sigma = fabs(bf0_p1sigma - fparam_bf[0]);
    // double bf0err_m1sigma = fabs(bf0_m1sigma - fparam_bf[0]);
    // double bf1err_p1sigma = fabs(bf1_p1sigma - fparam_bf[1]);
    // double bf1err_m1sigma = fabs(bf1_m1sigma - fparam_bf[1]);
    // bf0err = TMath::Max(bf0err_p1sigma,bf0err_m1sigma);
    // bf1err = TMath::Max(bf1err_p1sigma,bf1err_m1sigma);
    
    
    fparam_bf_errsym[0] = bf0err; // symmetric error
    fparam_bf_errsym[1] = bf1err;
    
#ifdef __print_hesseerror__
    cout << "bf0err = " << bf0err << endl;
    cout << "bf1err = " << bf1err << endl;
#endif
  }
  
  return;
}
//________________________________________________________________________________________
double FitMarginalization::OscPhysicalRegion(conv::PhysParamId_t p,double oscval)
{
  // NuBar
  if(p == conv::kPhypMixing13Bar){
    if(oscval < mixing13Bar_min) return mixing13Bar_min;
    else if(oscval > mixing13Bar_max) return mixing13Bar_max;
    else return oscval;
  }
  else if(p == conv::kPhypMixing23Bar){
    if(oscval < mixing23Bar_min) return mixing23Bar_min;
    else if(oscval > mixing23Bar_max) return mixing23Bar_max;
    else return oscval;
  }
  else if(p == conv::kPhypMixing12Bar){
    if(oscval < mixing12Bar_min) return mixing12Bar_min;
    else if(oscval > mixing12Bar_max) return mixing12Bar_max;
    else return oscval;
  }
  else if(p == conv::kPhypDm23sqBar){
    if(oscval < dmsq32Bar_min) return dmsq32Bar_min;
    else if(oscval > dmsq32Bar_max) return dmsq32Bar_max;
    else return oscval;
  }
  else if(p == conv::kPhypDm12sqBar){
    if(oscval < dmsq12Bar_min) return dmsq12Bar_min;
    else if(oscval > dmsq12Bar_max) return dmsq12Bar_max;
    else return oscval;
  }
  else if(p == conv::kPhypCPVBar){
    // Though CPV is periodic it's better to return errors in 2 Pi-Greek period
    if(oscval < deltacpBar_min) return deltacpBar_min;
    else if(oscval > deltacpBar_max) return deltacpBar_max;
    else return oscval;
  }
  
  // Nu
  else if(p == conv::kPhypMixing13){
    if(oscval < mixing13_min) return mixing13_min;
    else if(oscval > mixing13_max) return mixing13_max;
    else return oscval;
  }
  else if(p == conv::kPhypMixing23){
    if(oscval < mixing23_min) return mixing23_min;
    else if(oscval > mixing23_max) return mixing23_max;
    else return oscval;
  }
  else if(p == conv::kPhypMixing12){
    if(oscval < mixing12_min) return mixing12_min;
    else if(oscval > mixing12_max) return mixing12_max;
    else return oscval;
  }
  else if(p == conv::kPhypDm23sq){
    if(oscval < dmsq32_min) return dmsq32_min;
    else if(oscval > dmsq32_max) return dmsq32_max;
    else return oscval;
  }
  else if(p == conv::kPhypDm12sq){
    if(oscval < dmsq12_min) return dmsq12_min;
    else if(oscval > dmsq12_max) return dmsq12_max;
    else return oscval;
  }
  else if(p == conv::kPhypCPV){
    // Though CPV is periodic it's better to return errors in 2 * Pi-Greek period
    if(oscval < deltacp_min) return deltacp_min;
    else if(oscval > deltacp_max) return deltacp_max;
    else return oscval;
  }
  
  else{
    cerr << endl;
    cerr << "FitMarginalization::OscPhysicalRegion" << endl;
    cerr << "Parameter " << util::AsString(p) << " is not in the physics list" << endl;
    cerr << endl;
    util::Die(1);
  }
  
  return conv::kBadResult;
}

//NueBarapp RAJ
double FitMarginalization::FillExpectedNuEBarAppHisto(SuperkObservation *skobs_data,int iterpar0,bool readbinary){ // NEW PRECISION
  
#ifdef __print_chisqmarg_tstopwatch__
  TStopwatch sw;
#endif
  if(fNuEAppHisto[iterpar0]->GetEntries()!=0){
    std::cerr << "FitMarginalization/FillExpectedNuEBarAppHisto >> Histogram for beta bin: " << iterpar0
    << " already filled. Reseting histograms."
    << std::endl;
    fNuEAppHisto[iterpar0]->Reset();
    fNuEAppHistoUnscaled[iterpar0]->Reset();
    
  }
  
  if(fPhysParamFree[0] != kPhypBetaBar || fPhysParamFree[1] != kPhypUndefined || fteststatistic != "rate-only")
  {
    std::cerr <<  "FitMarg::FillExpectedNuEBarAppHisto >> Par0 is not BetaBar or Par1 is not undefined. "
    <<  "Can't deal with this" << std::endl;
  }
  
  int iterpar1 = 0;
  
  // NEW SPECTRA MARG
  // Set to 0 in case it has been already calculated before, for example when fabs(chisq)<1.E-8.
  //flikeVArray->SetEntry(0.,iterpar0,iterpar1); // not used
  //
  
  fNCallsCalcChisqMargPt++;
  
  // Read binary file
  if(readbinary){
    
    // Get file ID corresponding to the input free parameters
    
    int fileid = GetFileIDRead(iterpar0,iterpar1);
    const char *tag = ftagfile.c_str();
    
    if( (fCurrentFileIDminim==conv::kBadResult) || (fCurrentFileIDminim!=fileid) ){
      //if( (fileid==conv::kBadResult) || (fCurrentFileIDminim!=fileid) ){
      
#ifdef __print_chisqmarg_changebinary__
      cout << endl;
      cout << "Change binary file from " << fCurrentFileIDminim << " to " << fileid << endl;
#endif
      
      fCurrentFileIDminim = fileid;
      this->ReadToysFromBinary(fileid,tag);
    }
  }
  
  
  for(int iexperiment = 0; iexperiment < fntoys_marg; iexperiment++) {
    //cout << "toy=" << iexperiment << " - iterpar0 = " << iterpar0 << " - iterpar1 = " << iterpar1 << endl;
    
#ifdef __print_array_toy_fillchisq__
    cout << endl;
    cout << "Toy " << iexperiment << endl;
    cout << endl;
#endif
    
    // compute the likelihood for a single marginalization toy
    double chisq = 0.;
    
    // Loop over the used SK samples
    int iscount=0;
    //for(int is = 0; is < kns; is++) {
    //SuperkSample_t s = kSKSamples[is];
    //if(s == kSKsUndefined) continue;
    double NuEAppEvents= 0;
    for(int is=0;is<fNSKsamples;is++){
      SuperkSample_t s = this->GetSKSamples(is);
      
      
      // NEW FHC+RHC JOINT
      // Skip the sk samples which shouldn't be used in the chisq computation
      bool usesample = this->UseSKSample(is);
      if(!usesample){
        iscount++; // Increase the sk-sample counter, otherwise at the next iteration it takes the wrong sample
        continue;
      }
      //
      if( fteststatistic=="rate-only"){
        
#ifdef __print_array_toy_fillchisq__
        cout << util::AsString(s) << " --> # Ereco bins = " << util::NumErecoBins(s) << endl;
#ifdef __USE_2D_TEMPLATES__
        // NEW 2D TEMPLATES
        cout << util::AsString(s) << " --> # Q2 bins = " << util::NumQ2Bins(s) << endl;
        //
#endif
        double nobs_sample = 0.;
        double chisq_sample = 0.;
#endif
        double nexp_sample = 0.;
        
        
#ifndef __USE_2D_TEMPLATES__
        
        for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){
          
          // NEW PRECISION
          
          // NEW 2D TEMPLATES (VARRAY)
          //double nexp = (double)fToyMargVArray->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],iscount,ibin);
          double nexp = (double)this->GetToyMargVArraySample(s)->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],ibin);
          //
          
          nexp *= (double)fNormPOT[is];
          double nobs = (double)skobs_data->ErecoSpectrumArray(s)[ibin];
          if(s != kSKsRHCFCFV1RingELike)
          chisq += valor::util::ChisqPoisson(nobs,nexp); // penalty term is a simple constant --> not needed
          else
          nexp_sample += nexp;                                                         //
          
#ifdef __print_array_toy_fillchisq__
          nexp_sample += nexp;
          nobs_sample += nobs;
          chisq_sample += valor::util::ChisqPoisson(nobs,nexp);
#endif
          
        } // loop over ereco bins
        
        
        
#else // __USE_2D_TEMPLATES__
        
        // NEW 2D TEMPLATES
        
        for(int ibinEreco=0;ibinEreco<util::NumErecoBins(s);ibinEreco++){
          for(int ibinQ2=0;ibinQ2<util::NumQ2Bins(s);ibinQ2++){
            
            // NEW 2D TEMPLATES (VARRAY)
            //double nexp = (double)fToyMargVArray->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],iscount,ibinEreco,ibinQ2);
            double nexp = (double)this->GetToyMargVArraySample(s)->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],ibinEreco,ibinQ2);
            //
            
            nexp *= (double)fNormPOT[is];
            double nobs = (double)skobs_data->TwoDSpectrumVArray(s)->GetEntry(ibinEreco,ibinQ2); // Maybe it's slower than the simple array --> TO CHECK
            if(s != kSKsRHCFCFV1RingELike)
            chisq += valor::util::ChisqPoisson(nobs,nexp); // penalty term is a simple constant --> not needed
            else{
              nexp_sample += nexp;
            }
#ifdef __print_array_toy_fillchisq__
            nobs_sample += nobs;
            chisq_sample += valor::util::ChisqPoisson(nobs,nexp);
#endif
          } // loop over q2 bins
        } // loop over ereco bins
          //
#endif // __USE_2D_TEMPLATES__
        
#ifdef __print_array_toy_fillchisq__
        cout << "# exp sample = " << nexp_sample << " - # obs sample = " << nobs_sample << " --> chisq sample = " << chisq_sample << " --> chisq tot = " << chisq << endl;
#endif
        if(s == kSKsRHCFCFV1RingELike)
        NuEAppEvents = nexp_sample;
        
      } // rate+shape
      else{
        cerr << endl;
        cerr << "FitMarginalization::FillExpectedNuEBarAppHisto" << endl;
        cerr << "The test statistic is not available here: " << fteststatistic << endl;
        cerr << endl;
        exit(1);
      }
      
      iscount++;
    } // loop SK samples
    
    for(int nuebarextraction_int=0;nuebarextraction_int < 10000;nuebarextraction_int++)
    {
      int NuEAppEvents_fill = Rnd::Instance()->RndToyMCstat().Poisson(NuEAppEvents);
      fNuEAppHisto[iterpar0]->Fill(NuEAppEvents_fill,exp(-chisq/2.));
      fNuEAppHistoUnscaled[iterpar0]->Fill(NuEAppEvents_fill);
    }
  } // iexperiment
  
  //Calculate data p-value
  
  SuperkSample_t s = kSKsRHCFCFV1RingELike;
  double nobs_sample = 0;
#ifndef __USE_2D_TEMPLATES__
  for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){
    
    double nobs = (double)skobs_data->ErecoSpectrumArray(s)[ibin];
    nobs_sample += nobs;
    
  } // loop over ereco bins
#else // __USE_2D_TEMPLATES__
  
  // NEW 2D TEMPLATES
  
  for(int ibinEreco=0;ibinEreco<util::NumErecoBins(s);ibinEreco++){
    for(int ibinQ2=0;ibinQ2<util::NumQ2Bins(s);ibinQ2++){
      double nobs = (double)skobs_data->TwoDSpectrumVArray(s)->GetEntry(ibinEreco,ibinQ2); // Maybe it's slower than the simple array --> TO CHECK
      nobs_sample += nobs;
    } // loop over q2 bins
  } // loop over ereco bins
    //
#endif // __USE_2D_TEMPLATES__
  
#ifdef __print_array_toy_fillchisq__
  std::cout <<"FitMarg/FillExpectedNuEBarApp>> " << util::AsString(s) <<" - # obs  = " << nobs_sample << std::endl;
#endif
  
  int obs_bin = fNuEAppHisto[iterpar0]->FindBin(nobs_sample);
  
  double p_value_tail = 0;
  if(iterpar0 == 0){//Right tailed p value for beta=0 null hypothesis
    p_value_tail =   fNuEAppHisto[iterpar0]->Integral(obs_bin,fNuEAppHisto[iterpar0]->GetNbinsX());
  }
  else{  //Use left tailed p value for beta=1 null hypothesis
    p_value_tail =   fNuEAppHisto[iterpar0]->Integral(1,obs_bin);
  }
  
  double p_value_integral = fNuEAppHisto[iterpar0]->Integral();
  double p_value = p_value_tail/p_value_integral;
  
  
#ifdef __print_rate_only_pvalue__
  std::cout << "FitMarg/FillExpectedNuEBarApp>> " << util::AsString(s) <<" Nevents: " << nobs_sample <<" P-Value: " << p_value
  << " tail: " << p_value_tail << " integral: " << p_value_integral
  << std::endl;
#endif
  return p_value;
  
}

void FitMarginalization::SaveTreeOfToyParamLists(TFile *file)
{
  
  if(!fShouldSaveToyParamLists)
  {
    std::cerr << "Marginsalistion parameter values could not be stored as they were not saved. "
    << "Use FitMarginalization::SetSaveParamLists(bool should_save_lists)" << endl;
    return;
  }
  file->cd();
  TTree* t_pars = new TTree("marg_parameter_values","Input parameters used for marginalisation");
  
  //This is the only way I could make filling a tree with an array work
  //parpoints is an array to store pointers to array objects above to allow me to fill tree :(

  //conv::knppmax: large number for marginalisation parameters
  double params[conv::knppmax];
  double* parpoints[conv::knppmax];
  for(int i=0;i<conv::knppmax;i++){
    parpoints[i] = params + i;
  }
  
  PhysParamList * iparm;
  
  //use the first list to initialise the tree
  iparm = fSavedToyParamLists[0];
  int counter = 0;
  for (int cip = 0 ; cip < conv::knppmax ; cip++) {
    PhysParamId_t pip = (PhysParamId_t)cip;
    if (iparm->Contains(pip)) {
      t_pars->Branch(util::AsStringSimple(pip).c_str(),  (parpoints[counter]));
      counter++;
    }
  }
  
  for (int i = 0 ; i < fntoys_marg ; i++) {
    iparm = fSavedToyParamLists[i];
    int counter_2 = 0;
    for (int cip = 0 ; cip < conv::knppmax ; cip++) {
      PhysParamId_t pip = (PhysParamId_t)cip;
      if (iparm->Contains(pip)) {
        params[counter_2] = iparm->Value(pip);
        counter_2++;
      }
    }
    t_pars->Fill();
    
  }
  
  t_pars->Write();
  delete t_pars;
}
