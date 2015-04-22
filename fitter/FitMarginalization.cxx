//________________________________________________________________________________________
/*
 
 \copyright (c) 2010-2015, VaLOR group.

*/
//________________________________________________________________________________________

//#define __print_build_map_fileid__
//#define __print_array_toy_maketoys__
//#define __print_array_toy_maketoys_verbose__
//#define __print_array_toy_fillchisq__
//#define __print_array_toy_merge__
//#define __print_masshierarchy__
//#define __print_systematic_variations__
#define __print_chisqmarg__
//#define __readtoys_debug__
//#define __print_minimization_steps__ //
//#define __print_hesseerror__
//#define __print_hesseerror_verbose__

#define NUMOFBINSAROUNDMINTOCHECK 3  // check bins both on the right and left 

 
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

#include "Conventions/Opt.h"
#include "Utilities/Rnd.h"
#include "Utilities/UtlHist.h"

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
  flikeVArray                        = 0;
  fChisqMargVArray                   = 0;
  fnerecobins                        = 0;
  fsksamples                         = 0;
  fToyMargErecoVArray                = 0;
  fToyMargErecoVArray_Clone          = 0;
  fToyMargErecoVArray_RHCNumu_Clone  = 0;
  fToyMargErecoVArray_RHCNue_Clone   = 0;
  fToyMargErecoVArray_Numu_Clone     = 0;
  fToyMargErecoVArray_Nue_Clone      = 0;
  fhToyMargCurr_RHCNumu               = 0;
  fToyMargErecoVArray_RHCNumu        = 0;
  fToyMargErecoVArray_RHCNue         = 0;
  fToyMargErecoVArray_Numu           = 0;
  fToyMargErecoVArray_Nue            = 0;
  fhchisq_marg_1D                     = 0;
  fhchisq_marg_2D                     = 0;
  fhchisq_marg_1D_Step                = 0;
  fhchisq_marg_2D_Step                = 0;
  fhdeltachisq_marg_1D                = 0;
  fhdeltachisq_marg_2D                = 0;
  fhmap_chisqrange_fileid_1D          = 0;
  fhmap_chisqrange_fileid_2D          = 0;
  fparm                              = 0;
  fSavedToyParamLists                = 0;
  
  // to be removed once the sk sample issue is fixed
  for(int ikns=0;ikns<conv::kns;ikns++){
    fhToyMargCurr[ikns] = 0;
  }
      
  fShouldSaveToyParamLists           = false;
  fMaxSampleErecoBins                = 1;
}
//________________________________________________________________________________________
FitMarginalization::~FitMarginalization(void)
{
  this->CleanUp();
}
//________________________________________________________________________________________
void FitMarginalization::Init(SuperkExpectation *skexp){
   
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
  
  fMergeInputFiles_Called = true; // false the first time MergeInputFiles() is calles
  fNmergedFiles = 0; // fNmergedFiles++ everytime MergeInputFiles() is called (after reading each file)

  fNCallsCalcChisqMargPt = 0; // Get track of the # of calls to CalcChisqMargPt

  const int nmh=2; // mass hierarchy  

  // Get number of sk samples used in the fit
  // Needed to initialize the Ereco arrays

  this->SetNSKsamples(); // also defines fnerecobins
  
  fToyMargErecoVArray = new VALORArrayF(1,1,1,1,1,1); 
  fToyMargErecoVArray->SetAll(-1.);   // needed for MergeInputFiles()
  fToyMargErecoVArray->SetDimensionName(0,"# toys marginalization");
  fToyMargErecoVArray->SetDimensionName(1,"# bins free parameter 0");
  fToyMargErecoVArray->SetDimensionName(2,"# bins free parameter 1");
  fToyMargErecoVArray->SetDimensionName(3,"# mass hierarchies");
  fToyMargErecoVArray->SetDimensionName(4,"# sk samples");
  fToyMargErecoVArray->SetDimensionName(5,"# bins ereco");

  // need to declare array of histogram (then set to 0) because
  // problems with ROOT when declaring double pointers as global in header file 
  for(int is=0;is<conv::kns;is++){
    if(is<fNSKsamples){
      char name[200];
      SuperkSample_t s = this->GetSKSamples(is);    
      sprintf(name,"hToyMargCurr_%s",util::AsStringCompact(s).c_str());
      fhToyMargCurr[is] = new TH1F(name,name,util::NumErecoBins(s),util::ErecoBinEdges(s));   
      cout << name << endl;
    }
    else{
      fhToyMargCurr[is] = 0;
    }
  }
  //
  
  fChisqMargVArray            = new VALORArrayF(fnbins_param[0], fnbins_param[1]);
  fChisqMargVArray->SetAll(0.);
  flikeVArray                 = new VALORArrayF(fnbins_param[0], fnbins_param[1]);
  flikeVArray->SetAll(0.);
  
  char Xtitle[200];
  char Ytitle[200];
  sprintf(Xtitle,"%s",util::AsStringPaperRooTex(fPhysParamFree[0]).c_str()); 
  sprintf(Ytitle,"%s",util::AsStringPaperRooTex(fPhysParamFree[1]).c_str()); 



  // cout << "fmin_param[0] = " << fmin_param[0] << endl;
  // cout << "fmax_param[0] = " << fmax_param[0] << endl;
  // cout << "fstep_param[0] = " << fstep_param[0] << endl;
  // cout << "fmin_param[1] = " << fmin_param[1] << endl;
  // cout << "fmax_param[1] = " << fmax_param[1] << endl;
  // cout << "fstep_param[1] = " << fstep_param[1] << endl;


  //fhchisq_marg_1D = new TH1F("hchisq_marg_1D","hchisq_marg_1D",fnbins_param[0],fmin_param[0],fmax_param[0]+fstep_param[0]);  
  fhchisq_marg_1D = new TH1F("hchisq_marg_1D","hchisq_marg_1D",fnbins_param[0],fmin_param[0]-fstep_param[0]/2.,fmax_param[0]+fstep_param[0]/2.);  
  fhchisq_marg_1D->GetXaxis()->SetTitle(Xtitle);
  fhchisq_marg_1D->GetYaxis()->SetTitle("#Delta #chi^{2}");

  cout << endl;
  cout << "Initialize histogram fchisq_marg_1D:" << endl;
  cout << "BinCenter(1) = " << fhchisq_marg_1D->GetXaxis()->GetBinCenter(1) << endl;
  cout << "BinCenter(Nbins) = " << fhchisq_marg_1D->GetXaxis()->GetBinCenter(fnbins_param[0]) << endl;
  
  //fhchisq_marg_2D = new TH2F("hchisq_marg_2D","hchisq_marg_2D",
  //fnbins_param[0],fmin_param[0],fmax_param[0]+fstep_param[0],
  //fnbins_param[1],fmin_param[1],fmax_param[1]+fstep_param[1] );
  fhchisq_marg_2D = new TH2F("hchisq_marg_2D","hchisq_marg_2D",
			     fnbins_param[0],fmin_param[0]-fstep_param[0]/2.,fmax_param[0]+fstep_param[0]/2.,
			     fnbins_param[1],fmin_param[1]-fstep_param[1]/2.,fmax_param[1]+fstep_param[1]/2. );
  fhchisq_marg_2D->GetXaxis()->SetTitle(Xtitle);
  fhchisq_marg_2D->GetYaxis()->SetTitle(Ytitle);

  cout << endl;
  cout << "Initialize histogram fchisq_marg_2D:" << endl;
  cout << "X: BinCenter(1) = " << fhchisq_marg_2D->GetXaxis()->GetBinCenter(1) << endl;
  cout << "X: BinCenter(Nbins) = " << fhchisq_marg_2D->GetXaxis()->GetBinCenter(fnbins_param[1]) << endl;
  cout << "Y: BinCenter(1) = " << fhchisq_marg_2D->GetYaxis()->GetBinCenter(1) << endl;
  cout << "Y: BinCenter(Nbins) = " << fhchisq_marg_2D->GetYaxis()->GetBinCenter(fnbins_param[1]) << endl;
  
  fhdeltachisq_marg_1D = 0;
  fhdeltachisq_marg_2D = 0; 

  cout << endl;
  cout << "Initialize VALORArray fToyMargErecoVArray with:"<< endl;
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
	 << " - # bins = " << this->fnerecobins->GetEntry(is)
	 << endl;
  }
  cout << endl;
  cout << "fnbins_sub_param[0] = " << fnbins_sub_param[0] << endl;
  cout << "fnbins_sub_param[1] = " << fnbins_sub_param[1] << endl;
  cout << "nmh = " << nmh << endl;
  cout << endl;

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
  fhmap_chisqrange_fileid_1D = (TH1I*) fhchisq_marg_1D->Clone("hmap_chisqrange_fileid_1D");
  fhmap_chisqrange_fileid_2D = (TH2I*) fhchisq_marg_2D->Clone("hmap_chisqrange_fileid_2D");

  fCurrentFileIDminim = conv::kBadResult; // current file ID during the chisq minimization process


  // Store # of POT for each used sample
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s = this->GetSKSamples(is);
    fPOT[is] = skexp->GetSuperkMC()->AnalyzedSamplePOT(s);
    fNormPOT[is] = fPOT[is] / (1E+20);
    cout << "POT ( " << util::AsString(s) << " ) = " << fPOT[is] << endl;
    cout << "Norm POT ( " << util::AsString(s) << " ) = " << fNormPOT[is] << endl;
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
  
  if (flikeVArray                       ) { delete flikeVArray                       ; flikeVArray                        = 0;}
  if (fChisqMargVArray                  ) { delete fChisqMargVArray                  ; fChisqMargVArray                   = 0;}
  if (fnerecobins                       ) { delete fnerecobins                       ; fnerecobins                        = 0;}
  if (fsksamples                        ) { delete fsksamples                        ; fsksamples                         = 0;}
  if (fToyMargErecoVArray               ) { delete fToyMargErecoVArray               ; fToyMargErecoVArray                = 0;}
  if (fToyMargErecoVArray_Clone         ) { delete fToyMargErecoVArray_Clone         ; fToyMargErecoVArray_Clone          = 0;}
  if (fToyMargErecoVArray_RHCNumu_Clone ) { delete fToyMargErecoVArray_RHCNumu_Clone ; fToyMargErecoVArray_RHCNumu_Clone  = 0;}
  if (fToyMargErecoVArray_RHCNue_Clone  ) { delete fToyMargErecoVArray_RHCNue_Clone  ; fToyMargErecoVArray_RHCNue_Clone   = 0;}
  if (fToyMargErecoVArray_Numu_Clone    ) { delete fToyMargErecoVArray_Numu_Clone    ; fToyMargErecoVArray_Numu_Clone     = 0;}
  if (fToyMargErecoVArray_Nue_Clone     ) { delete fToyMargErecoVArray_Nue_Clone     ; fToyMargErecoVArray_Nue_Clone      = 0;}
  if (fhToyMargCurr_RHCNumu             ) { delete fhToyMargCurr_RHCNumu             ; fhToyMargCurr_RHCNumu              = 0;}
  if (fToyMargErecoVArray_RHCNumu       ) { delete fToyMargErecoVArray_RHCNumu       ; fToyMargErecoVArray_RHCNumu        = 0;}
  if (fToyMargErecoVArray_RHCNue        ) { delete fToyMargErecoVArray_RHCNue        ; fToyMargErecoVArray_RHCNue         = 0;}
  if (fToyMargErecoVArray_Numu          ) { delete fToyMargErecoVArray_Numu          ; fToyMargErecoVArray_Numu           = 0;}
  if (fToyMargErecoVArray_Nue           ) { delete fToyMargErecoVArray_Nue           ; fToyMargErecoVArray_Nue            = 0;}
  if (fhchisq_marg_1D                   ) { delete fhchisq_marg_1D                   ; fhchisq_marg_1D                    = 0;}
  if (fhchisq_marg_2D                   ) { delete fhchisq_marg_2D                   ; fhchisq_marg_2D                    = 0;}
  if (fhchisq_marg_1D_Step              ) { delete fhchisq_marg_1D_Step              ; fhchisq_marg_1D_Step               = 0;}
  if (fhchisq_marg_2D_Step              ) { delete fhchisq_marg_2D_Step              ; fhchisq_marg_2D_Step               = 0;}
  if (fhdeltachisq_marg_1D              ) { delete fhdeltachisq_marg_1D              ; fhdeltachisq_marg_1D               = 0;}
  if (fhdeltachisq_marg_2D              ) { delete fhdeltachisq_marg_2D              ; fhdeltachisq_marg_2D               = 0;}
  if (fhmap_chisqrange_fileid_1D        ) { delete fhmap_chisqrange_fileid_1D        ; fhmap_chisqrange_fileid_1D         = 0;}
  if (fhmap_chisqrange_fileid_2D        ) { delete fhmap_chisqrange_fileid_2D        ; fhmap_chisqrange_fileid_2D         = 0;}
  if (fparm                             ) { delete fparm                             ; fparm                              = 0;}
  if (fSavedToyParamLists               ) { delete fSavedToyParamLists               ; fSavedToyParamLists                = 0;}
  

  // to be removed once the sk sample issue is fixed
  for(int ikns=0;ikns<conv::kns;ikns++){
    if (fhToyMargCurr[ikns]              ) { delete fhToyMargCurr[ikns]               ; fhToyMargCurr[ikns]                = 0;} 
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
    exit(1);
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
      exit(1);
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
      exit(1);
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

  // The parameter is not selected
  else{
    if(ParId==0){
      cerr << "FitMarginalization::SetFreeParam" << endl;
      cerr << endl;
      cerr << "The parameter 0 has not been chosen. Always do 1D fit with parameter 0" << endl;
      cerr << endl;
      exit(1);
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
    exit(1);
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
  fhypothesis = fparm->PhysicsHypothesis();

  if( fPhysParamFree[0] == conv::kPhypUndefined && 
      fPhysParamFree[1] == conv::kPhypUndefined ){  
    cerr << endl;
    cerr << "FitMarginalization::SetFitParamList" << endl;
    cerr << "Both free parameters are not defined" << endl;
    cerr << endl;
    exit(1);
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
    exit(1);    
  }


  // Parameter 0
  if( fparm->IsFixed(fPhysParamFree[0]) ){
    cerr << "FitMarginalization::SetFitParamList" << endl;
    cerr << "The free parameter " << fnamefreeparam[0] 
	 << " is fixed in FitParamList " << fit_param_list
	 << endl << endl;
    exit(1);
  }
  else if( fparm->IsNuisance(fPhysParamFree[0]) ){
    cerr << "FitMarginalization::SetFitParamList" << endl;
    cerr << "The free parameter " << fnamefreeparam[0] 
	 << " is nuisance in FitParamList " << fit_param_list
	 << endl << endl;
    exit(1);
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
	 << endl << endl;
    exit(1);
  }
  else if( fparm->IsNuisance(fPhysParamFree[1]) ){
    cerr << "FitMarginalization::SetFitParamList" << endl;
    cerr << "The free parameter " << fnamefreeparam[1] 
	 << " is nuisance in FitParamList " << fit_param_list
	 << endl << endl;
    exit(1);
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
  const int nmh=2; // mass hierarchy  
  fToyMargErecoVArray->Resize(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nmh, fNSKsamples, fMaxSampleErecoBins);
  
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
    
    UInt_t myseed = 99;
    Rnd::Instance()->SetSeeds(myseed+iexperiment+10,myseed+iexperiment+20,myseed+iexperiment+30,myseed+iexperiment+40,fSeed+iexperiment+50); // stat, syst, osc, fitter, other
    
#ifdef __show_each_experiment__
    cout << endl;
    cout << "iexperiment " << iexperiment << endl;
    cout << "# params 0: " << fnbins_param[0] << endl;
    cout << "# params 1: " << fnbins_param[1] << endl;
#else
    if (iexperiment%100 == 0) { cout << "iexperiment " << iexperiment << endl; }
#endif
    
    fparm->ResetCurrToInitValues();
    
    // Extract oscillation parameters to be marginalized
    this->SetMargOscParam();

    // Set seed after oscillation parameter extraction

    UInt_t seedstat = Rnd::Instance()->RndToyMCstat().GetSeed();
    UInt_t seedsyst = Rnd::Instance()->RndToyMCsyst().GetSeed();
    UInt_t seedosc = Rnd::Instance()->RndToyMCosc().GetSeed();    
    UInt_t seedfit = Rnd::Instance()->RndFitter().GetSeed();
    UInt_t seedoth = Rnd::Instance()->RndOther().GetSeed();
    Rnd::Instance()->SetSeeds(seedstat+1,seedsyst+1,seedosc+1,seedfit+1,seedoth+1); // stat, syst, osc, fitter, other
    
    //

    SuperkObservation * skobs = new SuperkObservation;
    
    bool applied_negevents_condition = true;
    bool unphysical_tweak = true;
    while(applied_negevents_condition || unphysical_tweak) {
      if(tweak_systematics) {
	
	// tweak all sources of systematic error
	if(tweak_systematics == 1) {
	  //~ cout << "- Randomizing all systematic params in pdf..." << endl;
	  skexp->TweakAllSystematics(); 
	}
	else if(tweak_systematics == 2) {
	  //~ cout << "- Randomizing BANFF systematic params in pdf..." << endl;
	  skexp->TweakAllBANFFSystematics();
	}
	else if(tweak_systematics == 3) {
	  //~ cout << "- Randomizing uncorrelated-xsec systematic params in pdf..." << endl;
	  skexp->TweakAllUncorrelatedXSecSystematics();
	}
	else if(tweak_systematics == 4) {
	  //~ cout << "- Randomizing SKDetFSI systematic params in pdf..." << endl;
	  skexp->TweakAllSKDetectorFSISystematics();
	}

#ifdef __print_systematic_variations__
	skexp->PrintSystematicVariations();
#endif

	// Re-throw the systematics if at least one in the unphysical region
	PhysParamList true_nuisance_params = skexp->Tweaks();
	if (this->fRethrowSystUnphysical) {
	  unphysical_tweak = !(valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(true_nuisance_params));
	  if (unphysical_tweak) {
	    //~ cout<<"Parameters are unphysical - rethrowing"<<std::endl;
	    continue; // repeat the loop
	  }
	}
	else {
	  unphysical_tweak = false;
	}
	
      	// If we need to and haven't already save the toy param list
	if (fShouldSaveToyParamLists && (!fSavedToyParamLists[iexperiment])) {
	  // Get Syst tweaks
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
      
      }
      else { // no toys
	skexp->ResetSystematicVariations();
	unphysical_tweak = false; // otherwise it doesn't exit the while
	//~ cout << "- Generating toy-MC without systematic variation" << endl;
      }
    

      // Set seeds after systematic paramter extraction

      seedstat = Rnd::Instance()->RndToyMCstat().GetSeed();
      seedsyst = Rnd::Instance()->RndToyMCsyst().GetSeed();
      seedosc = Rnd::Instance()->RndToyMCosc().GetSeed();    
      seedfit = Rnd::Instance()->RndFitter().GetSeed();
      seedoth = Rnd::Instance()->RndOther().GetSeed();
      Rnd::Instance()->SetSeeds(seedstat+1,seedsyst+1,seedosc+1,seedfit+1,seedoth+1); // stat, syst, osc, fitter, other
    

      // Produce toys for each osc param value
      
      for(int ipar0=fmin_bin_sub_param[0];ipar0<=fmax_bin_sub_param[0];ipar0++){
	for(int ipar1=fmin_bin_sub_param[1];ipar1<=fmax_bin_sub_param[1];ipar1++){
	
	  //cout << ipar0 << " - " << ipar1 << endl;

	  double currval0 = fmin_param[0] + ipar0*fstep_param[0];
	  double currval1 = fmin_param[1] + ipar1*fstep_param[1];
	  //cout << "currval0 = " << currval0 << ", ";
	  //cout << "currval1 = " << currval1 << ", ";
	
	  // Marginalize the Mass Hierarchy
	  
	  conv::NuMassHierarchy_t hier;
	  conv::NuMassHierarchy_t hierBar;
	  
	  int nMassHierarchy = 2;
	  for(int iMH=0;iMH<nMassHierarchy;iMH++){
	    
	    // Don't marginalize Mass Hierarchy. Only the assumed one
	    if(this->fdoMargMassHierarchy == false){
	      
	      hier    = Opt::Instance()->gMassHierarchyAssumed;
	      hierBar = Opt::Instance()->gMassHierarchyAssumed;
	      
	      // Break the loop after the first iteration
	      nMassHierarchy = 1; 
	      
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
		exit(1);
	      }
	    }
	    
	    // Loop over both MH hypotheses and assign a different MH
	    else{
	      if( iMH==0 ){
		hier    = kNmhNormal;
		hierBar = kNmhNormal;
		skexp->GetSuperkMC()->SetMassHierarchy(hier);
		skexp->GetSuperkMC()->SetMassHierarchyBar(hierBar);
	      }
	      else if( iMH==1 ){
		hier    = kNmhInverted;
		hierBar = kNmhInverted;
		skexp->GetSuperkMC()->SetMassHierarchy(hier);
		skexp->GetSuperkMC()->SetMassHierarchyBar(hierBar);
	      }
	    }
	    
#ifdef __print_masshierarchy__
	    cout << endl;
	    cout << "iMH = " << iMH << " - MH: " << util::HierarchyAsString(hier) 
		 << " - MHbar: " << util::HierarchyAsString(hierBar)
		 << endl;
#endif
	    
	    // Update the expectation and produce the mock data
	    //fparm->Update(fPhysParamFree[0],currval0); 
	    //if(fIsFit2D) fparm->Update(fPhysParamFree[1],currval1); 
	    fparm->SetInitValue(fPhysParamFree[0],currval0); 
	    if(fIsFit2D) fparm->SetInitValue(fPhysParamFree[1],currval1); 
	    //cout << *fparm << endl;

	    applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), false); // build marginalization distribution w/ syst + osc param but no stat fluctuations
	    
	    // Loop over the used SK samples (This loop can be improved. Loop only over the used SK samples, pre-stored)
	    int iscount=0;
	    //for(int is = 0; is < kns; is++) {
	    //SuperkSample_t s = kSKSamples[is];
	    //if(s == kSKsUndefined) continue;
	    for(int is=0;is<fNSKsamples;is++){
	      SuperkSample_t s = this->GetSKSamples(is);
	      double * obsErec = skobs->ErecoSpectrumArray(s);
	      
	      // Build toy Ereco distributions
	      for(int ib = 0; ib < util::NumErecoBins(s); ib++){
		
		double bincont = obsErec[ib] / fNormPOT[is];
		
		fToyMargErecoVArray->SetEntry(bincont,iexperiment,ipar0-fmin_bin_sub_param[0],ipar1-fmin_bin_sub_param[1],iMH,iscount,ib); 
	      }
	      
	      iscount++; // index of the toy array
	    } // loop SK samples
	    
#ifdef __print_array_toy_maketoys__
	    cout << endl;
	    cout << "toy=" << iexperiment << endl; 
	    cout << " - ipar0 = " << ipar0 
		 << " - currval0 = " << currval0 << endl; 
	    if(fIsFit2D){
	      cout << " - ipar1 = " << ipar1 
		   << " - currval1 = " << currval1 << endl; 
	    }	 
	    cout << " : ";
	    
	    // Build toy Ereco distributions
	    
	    iscount=0;
	    //for(int is = 0; is < kns; is++) {
	    //SuperkSample_t s = kSKSamples[is];
	    //if(s == kSKsUndefined) continue;
	    for(int is=0;is<fNSKsamples;is++){
	      SuperkSample_t s = this->GetSKSamples(is);
	      
	      //cout << util::AsStringCompact(s) << "\t"<<skobs->ErecoSpectrumVArray(s)->Integral() / fNormPOT[is]<< endl;
	      cout << util::AsStringCompact(s) << "\t"<<skobs->ErecoSpectrumVArray(s)->Integral() << endl;
	      
#ifdef __print_array_toy_maketoys_verbose__
	      for(int ib = 0; ib < util::NumErecoBins(s); ib++){
		//double bincont = skobs->ErecoSpectrumArray(s)[ib];
		double bincont = fToyMargErecoVArray->GetEntry(iexperiment,ipar0-fmin_bin_sub_param[0],ipar1-fmin_bin_sub_param[1],iMH,iscount,ib); 
		cout << bincont << " ";
	      }
	      cout << endl;
	      cout << endl;
#endif // verbose
	      
	      iscount++; // index of the toy array
	    } // loop SK samples
#endif //__print_array_toy_maketoys__
	    
	  } // loop Mass Hierarchy
	} // loop ipar1 (free parameter)
      } // loop ipar0 (free parameter)
    } // applied_negevents_condition
    
    delete skobs;
  } // iexperiment

  cout << endl;
  cout << "Toys for marginalization are extracted" << endl;
  cout << endl;

  return;
}
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

  NtoysTot = fntoys_marg * 2;

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
float FitMarginalization::CalcChisqMargPtRead(SuperkObservation *skobs_data,int iterpar0,int iterpar1,bool readbinary,bool allow_cached){
  
  // If chisq already colculated return the stored value and break
  float chisq = fChisqMargVArray->GetEntry(iterpar0,iterpar1); 
  if(allow_cached && chisq > 0.) {return chisq;}

  fNCallsCalcChisqMargPt++;
  
  // Read binary file
  if(readbinary){ 

    // Get file ID corresponding to the input free parameters
    
    int fileid = GetFileIDRead(iterpar0,iterpar1);
    const char *tag = ftagfile.c_str();
    
    if( (fCurrentFileIDminim==conv::kBadResult) || (fCurrentFileIDminim!=fileid) ){
      //if( (fileid==conv::kBadResult) || (fCurrentFileIDminim!=fileid) ){
      fCurrentFileIDminim = fileid;
      this->ReadToysFromBinary(fileid,tag);
    }
  }
  
  //cout << "iterpar0 = " << iterpar0 << " - iterpar1 = " << iterpar1 << "   ";
  //double currval0 = fmin_param[0] + iterpar0*fstep_param[0];
  //double currval1 = fmin_param[1] + iterpar1*fstep_param[1];
  //cout << "currval0 = " << currval0 << endl;
  //cout << "currval1 = " << currval1 << endl;
  
  // Marginalize the Mass Hierarchy
  int nMassHierarchy = 2;
  
  for(int iMH=0;iMH<nMassHierarchy;iMH++){
    
    // Don't marginalize Mass Hierarchy. Only the assumed one
    if(this->fdoMargMassHierarchy == false){
      
      // Break the loop after the first iteration
      nMassHierarchy = 1; 
      
      // Used convention in MakeToysMarginalization(): 
      //  - iMH=0 --> kNmhNormal
      //  - iMH=1 --> kNmhInverted	
      if     ( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal   ) iMH = 0;
      else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted ) iMH = 1; 
      
      // Mass Hierarchy has to be already defined
      else {
	cerr << endl;
	cerr << "FitMarginalization::CalcChisqMargPtRead" << endl;
	cerr << "Mass Hierarchy is not defined" << endl;
	cerr << endl;
	exit(1);
      }
    }
    
#ifdef __print_masshierarchy__
    cout << "iMH = " << iMH << endl;
#endif	  
    
    for(int iexperiment = 0; iexperiment < fntoys_marg; iexperiment++) {
      //cout << "toy=" << iexperiment << " - iterpar0 = " << iterpar0 << " - iterpar1 = " << iterpar1 << endl;;
      
      // compute the likelihood for a single marginalization toy
      chisq = 0.;
      
      // Loop over the used SK samples
      int iscount=0;
      //for(int is = 0; is < kns; is++) {
      //SuperkSample_t s = kSKSamples[is];
      //if(s == kSKsUndefined) continue;
      for(int is=0;is<fNSKsamples;is++){
	SuperkSample_t s = this->GetSKSamples(is);
	
	for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){
	  float nexp = fToyMargErecoVArray->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],iMH,iscount,ibin);

	  nexp *= (float)fNormPOT[is];

	  float nobs = skobs_data->ErecoSpectrumArray(s)[ibin];
	  chisq += valor::util::ChisqPoisson(nobs,nexp); // penalty term is a simple constant --> not needed

	  //cout << "nexp = " << nexp << " - nobs = " << nobs << " - chisq = " << chisq << endl;
	  
	}	  
	iscount++;
      } // loop SK samples
      
#ifdef __print_array_toy_fillchisq__
      cout << endl;
      cout << "toy=" << iexperiment << " - iterpar0=" << iterpar0 << " - iterpar1=" << iterpar1 << " : ";
      
      iscount=0;
      //for(int is = 0; is < kns; is++) {
      //SuperkSample_t s = kSKSamples[is];
      //if(s == kSKsUndefined) continue;
      for(int is=0;is<fNSKsamples;is++){
	SuperkSample_t s = this->GetSKSamples(is);
	
	cout << util::AsStringCompact(s) << endl;
	
	for(int ib = 0; ib < util::NumErecoBins(s); ib++){
	  float bincont = fToyMargErecoVArray->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],iMH,iscount,ib);
	  cout << bincont * (float)fNormPOT[is] << " "; 
	}
	cout << endl;
	cout << endl;
	
	iscount++;
      } // loop SK samples
      
#endif
      
      float like = exp(-chisq/2.) / fntoys_marg;
      
#ifdef __print_array_toy_fillchisq__
      cout << "chisq = " << chisq << " - like = " << like << endl; 
#endif

      if(this->fdoMargMassHierarchy) like /= 2.;  // if MH is marginalized, divide by 2
      
      flikeVArray->AddToEntry(like,iterpar0,iterpar1); // this is correct if binary files differ only for free param sub-ranges
      
    } // iexperiment
  } // iMH (mass hierarchy)
  
  // Now the marginalized likelihood is computed using toys from all the binary files
  // Compute the chisq for each free osc parameter value and fill the array and histo
  float bincont = -2*log( flikeVArray->GetEntry(iterpar0,iterpar1) );
  fChisqMargVArray->SetEntry(bincont,iterpar0,iterpar1);
  if(!fIsFit2D) fhchisq_marg_1D->SetBinContent(iterpar0+1,bincont);
  else          fhchisq_marg_2D->SetBinContent(iterpar0+1,iterpar1+1,bincont);
      
  return bincont;
}

//________________________________________________________________________________________
void FitMarginalization::WriteToys2Binary(int idfile, const char *tag){

  cout << endl;
  cout << "Write marginalization toys in a binary file:" << endl;
  
  char namemh[200];
  if(fdoMargMassHierarchy) sprintf(namemh,"%s","marg");
  else{ 
    if( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal )
      sprintf(namemh,"%s","NH");
    else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted )
      sprintf(namemh,"%s","IH");
    else{
      cerr << "FitMarginalization::WriteToys2Binary" << endl;
      cerr << "Mass Hierarchy not defined" << endl;
      exit(1);
    }
  }	    

  char outname[200];
  sprintf(outname,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"",".bin");

  std::string soutname(outname);
  fToyMargErecoVArray->WriteFile(soutname);
  
  cout << " - Toy for marginalization are stored in " << soutname << endl;

  // Write store sub-range of the free parameter

  //int min_bin_sub_param0 = fmin_bin_sub_param[0];
  //int max_bin_sub_param0 = fmax_bin_sub_param[0];
  //int min_bin_sub_param1 = fmin_bin_sub_param[1];
  //int max_bin_sub_param1 = fmax_bin_sub_param[1];

  char outname_range[200];
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
  if( this->fdoMargMassHierarchy ) ostext << " Marginalize";
  else ostext << util::AsString(Opt::Instance()->gMassHierarchyAssumed) << endl;
  ostext << endl;
  if(fRethrowSystUnphysical) ostext << "Rethrow unphysical systematic toys" << endl;
  else ostext << "Don't rethrow unphysical systematic toys" << endl;
  ostext << "fntoys_marg = " << fntoys_marg << endl;
  ostext << endl;
  ostext << "Initialize VALORArray fToyMargErecoVArray with:"<< endl;
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
  
  ostext.close();
}
//________________________________________________________________________________________
void FitMarginalization::WriteToys2BinaryMerge(const char *tag){

  cout << endl;
  cout << "Write marginalization toys in a binary file:" << endl;
  
  char namemh[200];
  if(fdoMargMassHierarchy) sprintf(namemh,"%s","marg");
  else{ 
    if( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal )
      sprintf(namemh,"%s","NH");
    else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted )
      sprintf(namemh,"%s","IH");
    else{
      cerr << "FitMarginalization::WriteToys2BinaryMerge" << endl;
      cerr << "Mass Hierarchy not defined" << endl;
      exit(1);
    }
  }	    

  char outtag[200] = "%s_%s";
  sprintf(outtag,"%s_%s",tag,"Merged");

  char outname[200];
  sprintf(outname,fnamefilebinary,outtag,0,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"",".bin");

  std::string soutname(outname);
  fToyMargErecoVArray_Clone->WriteFile(soutname);

  // Write store sub-range of the free parameter

  //int min_bin_sub_param0 = fmin_bin_sub_param[0];
  //int max_bin_sub_param0 = fmax_bin_sub_param[0];
  //int min_bin_sub_param1 = fmin_bin_sub_param[1];
  //int max_bin_sub_param1 = fmax_bin_sub_param[1];

  char outname_range[200];
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

  char outname_text[200];
  std::string soutname_text(outname_text);
  sprintf(outname_text,fnamefilebinary,outtag,0,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"",".txt");

  ofstream ostext(outname_text,ios::out);
  ostext << "fnamefreeparam0: " << fnamefreeparam[0] << endl;
  ostext << "fnamefreeparam1: " << fnamefreeparam[1] << endl;
  ostext << "Mass hierarchy: "; 
  if( this->fdoMargMassHierarchy ) ostext << " Marginalize";
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

  char namemh[200];
  if(fdoMargMassHierarchy) sprintf(namemh,"%s","marg");
  else{ 
    if( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal )
      sprintf(namemh,"%s","NH");
    else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted )
      sprintf(namemh,"%s","IH");
    else{
      cerr << "FitMarginalization::ReadToysFromBinary" << endl;
      cerr << "Mass Hierarchy not defined" << endl;
      exit(1);
    }
  }	    
 

  // Write store sub-range of the free parameter

  int min_bin_sub_param0=0;
  int max_bin_sub_param0=0;
  int min_bin_sub_param1=0;
  int max_bin_sub_param1=0;
  int nbins_param0=0;
  int nbins_param1=0;

  char inname_range[200];
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
    cerr << "  - MH = " << namemh << endl;
    //cerr << "  - # of binary files = " << << endl;
    cerr << endl;
    exit(1);
  }

  // cout << fToyMargErecoVArray->GetSize() << endl;
  // cout << fToyMargErecoVArray->GetNDim() << endl;
  // cout << fToyMargErecoVArray->GetSize(0) << endl; // # toys marginalization
  // cout << fToyMargErecoVArray->GetSize(1) << endl; // # bins free parameter 0
  // cout << fToyMargErecoVArray->GetSize(2) << endl; // # bins free parameter 1
  // cout << fToyMargErecoVArray->GetSize(3) << endl; // # mass hierarchies
  // cout << fToyMargErecoVArray->GetSize(4) << endl; // # sk samples
  // cout << fToyMargErecoVArray->GetSize(5) << endl; // # bins ereco

  // int NtoyVArray = fToyMargErecoVArray->GetSize(0); // # toys marginalization
  // int Nbins0VArray = fToyMargErecoVArray->GetSize(1);  // # bins free parameter 0
  // int Nbins1VArray = fToyMargErecoVArray->GetSize(2);  // # bins free parameter 1
  // int NmhVArray = fToyMargErecoVArray->GetSize(3); // # mass hierarchies
  // int NsksampleVArray = fToyMargErecoVArray->GetSize(4); // # sk samples
  // int NbinserecoVArray = fToyMargErecoVArray->GetSize(5); // # bins ereco

  // if(NtoyVArray != fntoys_marg){
  //   cerr << endl;
  //   cerr << "# of toys from file = " << NtoyVArray << endl;
  //   cerr << "# of toys from input = " << fntoys_marg << endl;
  //   cerr << "If the chosen binary file is correct, set # of toys to " << NtoyVArray << endl;
  //   cerr << endl;
  //   exit(1);
  // }






  char inname[200];
  sprintf(inname,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"",".bin");
  
  std::string sinname(inname);
  //cout << sinname << endl;
  
#ifdef __readtoys_debug__
  cout << "Reading toys for marginalization from binary file: " << sinname << endl;  
#endif
  
  fToyMargErecoVArray->ReadFile(sinname);
  
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
    exit(1);
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
    exit(1);
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
void FitMarginalization::MinimizeHisto(){
  
  int minbin = conv::kBadResult;
  int minbinX = conv::kBadResult;
  int minbinY = conv::kBadResult;
  int minbinZ = conv::kBadResult;

  cout << endl;
  cout << "Find minimum of marginalized chisq" << endl;
  cout << endl;

  if(!fIsFit2D){
    fchisq_bf = fhchisq_marg_1D->GetMinimum();
    minbin = fhchisq_marg_1D->GetMinimumBin();
    fparam_bf[0] = fhchisq_marg_1D->GetBinCenter(minbin);
    fparam_bf[1] = conv::kBadResult;
  }
  else{
    fchisq_bf = fhchisq_marg_2D->GetMinimum();
    fhchisq_marg_2D->GetMinimumBin(minbinX,minbinY,minbinZ);
    fparam_bf[0] = fhchisq_marg_2D->GetXaxis()->GetBinCenter(minbinX);
    fparam_bf[1] = fhchisq_marg_2D->GetYaxis()->GetBinCenter(minbinY);
  }  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::MinimizeHisto(TH1F *hInChisqMarg){
  
  int minbin = conv::kBadResult;
  
  fchisq_bf = hInChisqMarg->GetMinimum();
  minbin = hInChisqMarg->GetMinimumBin();
  fparam_bf[0] = hInChisqMarg->GetBinCenter(minbin);
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::MinimizeHisto(TH2F *hInChisqMarg){
  
  int minbin = conv::kBadResult;
  
  fchisq_bf = hInChisqMarg->GetMinimum();
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
TH1F * FitMarginalization::GetToyMargHisto(SuperkSample_t sksample,int itoy,int ifreepar0,int ifreepar1,int iMassHierarchy){
        
  for(int is=0;is<fNSKsamples;is++){
    SuperkSample_t s = this->GetSKSamples(is);    
    if(sksample == s){
      for(int ibin=0;ibin<fnerecobins->GetEntry(is);ibin++){
	float bincont = fToyMargErecoVArray->GetEntry(itoy,ifreepar0,ifreepar1,iMassHierarchy,is,ibin);
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
  exit(1);
}
//________________________________________________________________________________________
void FitMarginalization::ResetChisqMarg(int *min_bin_sub_param,int *max_bin_sub_param){
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
    }
  }
  return;
}
//________________________________________________________________________________________
void FitMarginalization::BuildDeltaChisqMarg(){

  if(fchisq_bf==conv::kALargeNum){
    cerr << endl;
    cerr << "FitMarginalization::BuildDeltaChisqMarg" << endl;
    cerr << "The best-fit has not been found" << endl;
    cerr << "Can't get the marginalized DeltaChisq distribution" << endl;
    cerr << endl;
    exit(1);
  }
  
  // 1D
  if(!fIsFit2D){
    fhdeltachisq_marg_1D = (TH1F*) fhchisq_marg_1D->Clone("fhdeltachisq_marg_1D");
    
    for(int ibin=0;ibin<fhdeltachisq_marg_1D->GetNbinsX();ibin++){
      double bincont = fhdeltachisq_marg_1D->GetBinContent(ibin+1);
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
	exit(1);
      }
      fhdeltachisq_marg_1D->SetBinContent(ibin+1,bincont);    
    }
  }

  // 2D
  else{
    fhdeltachisq_marg_2D = (TH2F*) fhchisq_marg_2D->Clone("fhdeltachisq_marg_2D");

    for(int ibinX=0;ibinX<fhdeltachisq_marg_2D->GetNbinsX();ibinX++){
      for(int ibinY=0;ibinY<fhdeltachisq_marg_2D->GetNbinsY();ibinY++){

	double bincont = fhdeltachisq_marg_2D->GetBinContent(ibinX+1,ibinY+1);
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
	  exit(1);
	}
	fhdeltachisq_marg_2D->SetBinContent(ibinX+1,ibinY+1,bincont);    
      }
    }
  }

  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetPrior(PhysParamId_t physpar,string priorname,double par1,double par2){
  fparm->SetPriorMode(physpar,priorname);

  if(priorname=="flat"){
    fparm->SetMinValueFlat(physpar,par1);
    fparm->SetMaxValueFlat(physpar,par2);
  }
  else if(priorname=="gaus"){
    fparm->SetMeanGaus(physpar,par1);
    fparm->SetSigmaGaus(physpar,par2);
  }
  else if(priorname=="user"){
    cerr << endl;
    cerr << "FitMarginalization::SetPrior" << endl;
    cerr << "You chose a user prior but you didn't give any histogram as input" << endl;
    exit(1);
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::SetPrior" << endl;
    cerr << "You chose the wrong prior option" << endl;
    exit(1);
  }
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetPrior(PhysParamId_t physpar,string priorname,TH1F *hprior){
  fparm->SetPriorMode(physpar,priorname);
  
  if(priorname=="user"){
    fparm->SetPriorHisto(physpar,hprior);
  }
  else if(priorname=="flat" || priorname=="gaus"){
    cerr << endl;
    cerr << "FitMarginalization::SetPrior" << endl;
    cerr << "You chose a " << priorname << " prior, but you gave as input an histogram" << endl;
    exit(1);
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::SetPrior" << endl;
    cerr << "You chose the wrong prior option" << endl;
    exit(1);
  }

  return;
}
//________________________________________________________________________________________
string FitMarginalization::GetPrior(PhysParamId_t physpar) const
{
  return fparm->GetPriorModeName(physpar);
}
//________________________________________________________________________________________
void FitMarginalization::SetMargOscParam(){
  
  // Loop over the fit parameter list
  int np = fparm->GetNumberParameters();
  int npfree = 0;
  for(int ip = 0; ip < np; ip++) {
    conv::PhysParamId_t p = fparm->WhichParam(ip);
    //cout << p << " " << util::AsString(p) << endl;

    //bool fdoMarginalize = true;
    
    if(p!=fPhysParamFree[0]) 
      fparm->SetDoMarg(fPhysParamFree[0],true); 
    else{
      fparm->SetDoMarg(fPhysParamFree[0],false);  // Do not marginalize if free parameter 
      continue; 
    }
    
    if(fIsFit2D){
      if(p!=fPhysParamFree[1]) 
	fparm->SetDoMarg(fPhysParamFree[1],true); 
      else{
	fparm->SetDoMarg(fPhysParamFree[1],false);  // Do not marginalize if free parameter
	continue; 
      }
    }

    if(util::IsSystematic(p)) continue; // not oscillation parameters 
    if(fparm->IsFixed(p))     continue; // do not marginalize if fixed

    // Now extract the oscillation parameter to marginalize 
    double currval = conv::kBadResult;

    bool rethrow=true;
    while(rethrow){
      currval = Randomize(p);
      rethrow = RethrowOscUnphysical(p,currval);
      if(rethrow) cout << "Rethrow " << util::AsString(p) << " " << currval << endl;
    }

    //cout << "Marginalize " << util::AsString(p) << " --> " << currval << endl;

    // Update current values of parameters 
    //fparm->Update(p,currval);
    fparm->SetInitValue(p,currval);
    
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
        this->SetPrior(p,"flat",badmin,badmax);
      }
      else {
        double sigma = fracerr*mean;
        //~ std::cout<<"Setting default prior: "<<util::AsString(p)<<"\tGaus("<<mean<<","<<sigma<<")"<<std::endl;
        this->SetPrior(p,"gaus",mean,sigma);
      }
    }
    else if (util::IsPhase(p)) {
      this->SetPrior(p,"flat",dcp_min,dcp_max);
    }
        
    else {
      this->SetPrior(p,"flat",badmin,badmax);
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
    exit(1);
  }
  
  return false;
}


//________________________________________________________________________________________
double FitMarginalization::Randomize(PhysParamId_t physpar) const
{

  if(fparm->GetPriorModeName(physpar)=="flat"){
    double min = fparm->GetMinValueFlat(physpar);
    double max = fparm->GetMaxValueFlat(physpar);
    return Rnd::Instance()->RndToyMCosc().Uniform(min,max);
  }
  else if(fparm->GetPriorModeName(physpar)=="gaus"){
    double mean = fparm->GetMeanGaus(physpar);
    double sigma = fparm->GetSigmaGaus(physpar);
    return Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
  }
  else if(fparm->GetPriorModeName(physpar)=="user"){
    return fparm->GetPriorHisto(physpar)->GetRandom();
  }
  else{
    cerr << endl;
    cerr << "FitMarginalization::Randomize" << endl;
    cerr << "Used prior is: " << fparm->GetPriorModeName(physpar) << endl;
    cerr << "Prior not well defined for " << util::AsString(physpar) << endl;
    exit(1);
  }
  
}
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
  const int nmh=2; // mass hierarchy
  if(fMergeInputFiles_Called){
    for(int is=0;is<fNSKsamples;is++)
      fToyMargErecoVArray_Clone = new VALORArrayF(fntoys_marg, fnbins_param[0],fnbins_param[1], nmh, fNSKsamples, fnerecobins->GetEntry(is)); // total range of free parameter
    
    fToyMargErecoVArray_Clone->SetAll(-1.);

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

	int nMassHierarchy = 2;
	for(int iMH=0;iMH<nMassHierarchy;iMH++){


	  // Loop over the used SK samples
	  int iscount=0;
	  //for(int is = 0; is < kns; is++) {
	  //SuperkSample_t s = kSKSamples[is];
	  //if(s == kSKsUndefined) continue;
	  for(int is=0;is<fNSKsamples;is++){
	    SuperkSample_t s = this->GetSKSamples(is);
	    
	    for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){
	      float bincont = fToyMargErecoVArray->GetEntry(iexperiment,ipar0-fmin_bin_sub_param[0],ipar1-fmin_bin_sub_param[1],iMH,iscount,ibin); // shift to sub-range // bincont = -1 by default
	      if(bincont >= 0.){
		float bincont_clone = fToyMargErecoVArray_Clone->GetEntry(iexperiment,ipar0,ipar1,iMH,iscount,ibin);  // total range
		if(bincont_clone >= 0.){
		  cerr << endl;
		  cerr << "FitMarginalization::MergeInputFiles" << endl;
		  cerr << "Following bin has been already filled:" << endl;
		  cerr << " - iexperiment = " << iexperiment << endl;
		  cerr << " - ipar0 = " << ipar0 << endl;
		  cerr << " - ipar1 = " << ipar1 << endl;
		  cerr << " - iMH = " << iMH << endl;
		  cerr << " - sample = " << util::AsStringCompact(s) << endl;
		  cerr << " - iereco = " << ibin << endl;
		  exit(1);
		}
		else{
		  fToyMargErecoVArray_Clone->SetEntry(bincont,iexperiment,ipar0,ipar1,iMH,iscount,ibin); 
		}
	      }

	    } // loop ibin
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

	    for(int ib = 0; ib < util::NumErecoBins(s); ib++){
	      double bincont = fToyMargErecoVArray_Clone->GetEntry(iexperiment,ipar0,ipar1,iMH,iscount,ib);
	      cout << bincont << " "; 
	    }
	    cout << endl;
	    cout << endl;
	    
	    iscount++;
	  } // loop SK samples
#endif
	  
	} // loop iMH
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
  this->fsksamples = new VALORArrayI(fNSKsamples);

  int iscount=0;
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     
     int nbins = util::NumErecoBins(s);
     
     this->fsksamples->SetEntry(is,iscount);
     this->fnerecobins->SetEntry(nbins,iscount);
     
     if (nbins > fMaxSampleErecoBins) { fMaxSampleErecoBins = nbins; }
     
     //cout << "is " << is << " - iscount " << iscount << endl;
     //cout << fnerecobins->GetEntry(iscount) << endl;
     
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

  float bigvalueUOflow = 99999999.; 

  // 1D
  fhchisq_marg_1D_Step = (TH1F*)fhchisq_marg_1D->Clone("hchisq_marg_1D_Step");
  for(int ibinX=0;ibinX<fhchisq_marg_1D_Step->GetNbinsX();ibinX++){
    fhchisq_marg_1D_Step->SetBinContent(ibinX+1,0.);
  }
  // Set under and over flow of X to very high value
  for(int ibinX=0;ibinX<fhchisq_marg_1D->GetNbinsX();ibinX++){    
    fhchisq_marg_1D->SetBinContent(0,bigvalueUOflow);
    fhchisq_marg_1D->SetBinContent(fhchisq_marg_2D->GetNbinsX()+1,bigvalueUOflow);
  }

  // 2D
  fhchisq_marg_2D_Step = (TH2F*)fhchisq_marg_2D->Clone("hchisq_marg_2D_Step");
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

  float bigvalueUOflow = 99999999.; 

  // 1D
  fhchisq_marg_1D_Step = (TH1F*)fhchisq_marg_1D->Clone("hchisq_marg_1D_Step");
  for(int ibinX=0;ibinX<fhchisq_marg_1D_Step->GetNbinsX();ibinX++){
    fhchisq_marg_1D_Step->SetBinContent(ibinX+1,0.);
  }
  // Set under and over flow of X to very high value
  for(int ibinX=0;ibinX<fhchisq_marg_1D->GetNbinsX();ibinX++){    
    fhchisq_marg_1D->SetBinContent(0,bigvalueUOflow);
    fhchisq_marg_1D->SetBinContent(fhchisq_marg_2D->GetNbinsX()+1,bigvalueUOflow);
  }

  // 2D
  fhchisq_marg_2D_Step = (TH2F*)fhchisq_marg_2D->Clone("hchisq_marg_2D_Step");
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
  
  TRandom3 *r = new TRandom3(0); // needed to randomly scan the distribution

  double UniFirstBin_X = fFirstMinimBin_X; 
  double UniLastBin_X  = fLastMinimBin_X; 
  double UniFirstBin_Y = fFirstMinimBin_Y; 
  double UniLastBin_Y  = fLastMinimBin_Y; 

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
	  double bincont = conv::kALargeNum;
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
      float min = conv::kALargeNum;

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
	exit(1);
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
	  exit(1);
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
	exit(1);
      }

      // Compute the gradient between the lowest bin and the opposite one and the new step size
           
      double grad_X=0.;
      double grad_Y=0.;
      double step_X=0.;
      double step_Y=0.;
      double num=0.;
      double den=0.;

      
      if(!fIsFit2D){

	int NextBin_X = binN_X + (NewDir_X-1); // NewDir=0,1,2. NewDir=1 if binN
	double NextVal_X = fhchisq_marg_1D->GetXaxis()->GetBinCenter(NextBin_X);
	double BinVal_X = binlowedge_X; 
	double NewVal_X = fhchisq_marg_1D->GetXaxis()->GetBinCenter(NextBin_X);

#if defined __print_minimization_steps__
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
	
	double NextVal_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(NextBin_X);
	double NextVal_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(NextBin_Y);
	
	double BinVal_X = binlowedge_X; 
	double BinVal_Y = binlowedge_Y;
	
	double NewVal_X = fhchisq_marg_2D->GetXaxis()->GetBinCenter(NextBin_X);
	double NewVal_Y = fhchisq_marg_2D->GetYaxis()->GetBinCenter(NextBin_Y);
		
#if defined __print_minimization_steps__
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
	exit(1);
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
	  exit(1);
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

  double like_iter = conv::kBadResult;
  double like_prev = conv::kBadResult;
  

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
    exit(1);
  }

  float startX = fhchisq_marg_2D->GetXaxis()->GetBinCenter(fStartBin_X);
  float startY = fhchisq_marg_2D->GetYaxis()->GetBinCenter(fStartBin_Y);

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
    exit(1);
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
float FitMarginalization::GetChisqMargTrueRead(double valpar0,double valpar1,SuperkObservation *skobs,bool readbinary){

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
    exit(1);
  }
  if( fIsFit2D && (valpar1<min_Y || valpar1>max_Y) ){
    cerr << endl;
    cerr << "FitMarginalization::GetChisqMargTrueRead" << endl;
    cerr << "Parameter1 " << fnamefreeparam[1] << " is outside the total allowed range" << endl;
    cerr << "Allowed range: " << min_Y << " - " << max_Y << endl;
    cerr << "Param1 value: " << valpar1 << endl;
    cerr << endl;
    exit(1);
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
void FitMarginalization::SetTagFileName(const char *tag){
  ftagfile = string(tag); // Get the file tag name 
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetFileID(int nfiles){

  cout << endl;
  cout << "Build the map between the binary file ID and the free parameters sub-range" << endl;
  cout << endl;

  char namemh[200];
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
      exit(1);
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
    
    char inname[200];
    sprintf(inname,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"",".bin");
    std::string sinname(inname);
    fnamefile[idfile] = sinname;

#ifdef __print_build_map_fileid__
    cout << "Reading free param sub-range from binary file: " << fnamefile[idfile] << endl;
#endif




    // // Check # of binary files
    // //cout << fnamefilebinary << endl;
    // //cout << tag << endl;
    // char numfiletag[200];
    // char ext[200] = ".bin";
    // sprintf(numfiletag,"%s","./");
    // //sprintf(numfiletag,"../../inputs/marginalization/toys/nueapp_2D_NormPOT_fixosc_10toys/submit_sinsqth13Vsdcp_NormPOT_postfit*");
    // //TString snumfiletag(numfiletag);
    // //TIter next(files);
    // //exit(1);
    // cout << numfiletag << endl;
    // TSystemDirectory dir(numfiletag, numfiletag);
    // TList *files = dir.GetListOfFiles();
    // if (files) {
    //   TSystemFile *file;
    //   TString fname;
    //   TIter next(files);
    //   while ((file=(TSystemFile*)next())) {
    // 	fname = file->GetName();
    // 	if (!file->IsDirectory() && 
    // 	    fname.BeginsWith("../../inputs/marginalization/toys/nueapp_2D_NormPOT_fixosc_10toys/submit_sinsqth13Vsdcp_NormPOT_postfit") //&&
    // 	    //fname.EndsWith(".bin")
    // 	    ){
    // 	  cout << fname.Data() << endl;
    // 	}
    //   }
    // }    
    // exit(1);



    char inname_range[200];
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
      cerr << "  - MH = " << namemh << endl;
      //cerr << "  - # of binary files = " << << endl;
      cerr << endl;
      exit(1);
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
  	exit(1);
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
  
  double error = conv::kBadResult;

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
    exit(1);
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
  
  double error = conv::kBadResult;
  const double CritChisq_1D = 1.;
  const double CritChisq_2D = 2.3;

  // Find 68.27% CL contour (1sigma)

  if(!fIsFit2D){   // 1D
        
    // Get Delta Chisq distribution
    this->BuildDeltaChisqMarg();
    //fhdeltachisq_marg_1D = (TH1F*) fhchisq_marg_1D->Clone("fhdeltachisq_marg_1D");
        
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
      // 	   << " - X = " << fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(ibin)
      // 	   << " - bincont = " << bincont << endl;
      
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
	exit(1);
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
      // 	   << " - X = " << fhdeltachisq_marg_1D->GetXaxis()->GetBinCenter(ibin)
      //  	   << " - bincont = " << bincont << endl;

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
      cerr << "bin_possigma_a = " << bin_possigma_a << endl;
      cerr << "bin_negsigma_a = " << bin_negsigma_a << endl;
      cerr << endl;      
      exit(1);
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
    //fhdeltachisq_marg_1D = (TH1F*) fhchisq_marg_1D->Clone("fhdeltachisq_marg_1D");

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
      exit(1);
    }
    if(binleftY_y > binrightY_y){
      cerr << endl;
      cerr << "FitMarginalization::CalcBestFitErrorsMINOS" << endl;
      cerr << "Y left extreme bin is higher than right extreme bin" << endl;
      cerr << "The 1sigma error interval is messed up:" << endl;
      cerr << "binleftY = " << binleftY_y << endl;
      cerr << "binrightY = " << binrightY_y << endl;
      cerr << endl;
      exit(1);
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
void FitMarginalization::CalcBestFitErrorsHESSE(SuperkObservation *skobs_data,bool readbinary,int Nh){

  cout << "Compute HESSE errors" << endl;
#ifdef __print_hesseerror__
  cout << "Bins far " << Nh << " from the best-fit bin are used" << endl;
#endif

  double error = conv::kBadResult;
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
    // 	if( (bin_bf0+i <= fhdeltachisq_marg_2D->GetNbinsX()) &&
    // 	    (bin_bf1+j <= fhdeltachisq_marg_2D->GetNbinsY()) ) 
    // 	  CalcChisqMargPtRead(skobs_data,bin_bf0+i-1,bin_bf1+j-1,true);	
    // 	if( (bin_bf0-i >= 1) &&
    // 	    (bin_bf1-j >= 1) ) 
    // 	  CalcChisqMargPtRead(skobs_data,bin_bf0-i-1,bin_bf1-j-1,true);	
    // 	if( (bin_bf0-i >= 1) &&
    // 	    (bin_bf1+j <= fhdeltachisq_marg_2D->GetNbinsY()) ) 
    // 	CalcChisqMargPtRead(skobs_data,bin_bf0-i-1,bin_bf1+j-1,true);	
    // 	if( (bin_bf0+i <= fhdeltachisq_marg_2D->GetNbinsX()) &&
    // 	    (bin_bf1-j >= 1) ) 
    // 	CalcChisqMargPtRead(skobs_data,bin_bf0+i-1,bin_bf1-j-1,true);
    //   }
    // }
    // // Be sure the +-Nh bins are in the physical region
    // if( (bin_bf0-Nh) < 1 ||
    // 	(bin_bf0+Nh) > fhdeltachisq_marg_2D->GetNbinsX() ||
    // 	(bin_bf1-Nh) < 1 ||
    // 	(bin_bf1+Nh) > fhdeltachisq_marg_2D->GetNbinsY() 
    // 	)
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
      exit(1);
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
    exit(1);
  }
  
  return conv::kBadResult;
}
