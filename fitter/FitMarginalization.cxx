//________________________________________________________________________________________
/*
 
 \copyright (c) 2010-2015, VaLOR group.

*/
//________________________________________________________________________________________

#define __print_array_toy_maketoys__
#define __print_minimization_steps__
//#define __print_array_toy_fillchisq__
//#define __print_array_toy_merge__
//#define __print_masshierarchy__
//#define __print_systematic_variations__
#define __print_chisqmarg__

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
  hToyMargCurr_RHCNumu               = 0;
  fToyMargErecoVArray_RHCNumu        = 0;
  fToyMargErecoVArray_RHCNue         = 0;
  fToyMargErecoVArray_Numu           = 0;
  fToyMargErecoVArray_Nue            = 0;
  hchisq_marg_1D                     = 0;
  hchisq_marg_2D                     = 0;
  hchisq_marg_1D_Step                = 0;
  hchisq_marg_2D_Step                = 0;
  fparm                              = 0;
  fSavedToyParamLists                = 0;
  
  
  fShouldSaveToyParamLists           = false;
  fMaxSampleErecoBins                = 1;
}
//________________________________________________________________________________________
FitMarginalization::~FitMarginalization(void)
{
  this->CleanUp();
}
//________________________________________________________________________________________
void FitMarginalization::Init(){
   
  // Before Init() we need SetPhysParFree()
  
  doMargMassHierarchy = false;
  fRethrowSystUnphysical = true; // rethrow toys if systematics are unphysical

  fIsFit2D = false; // 1D fit by default

  fSeed = 999;

  fchisq_bf = 9999999.;
  
  sprintf(fnamefilebinary,"%s","%s_file%d_toy%d_x_oscpar%d_oscpar%d_MH%s%s%s");  
  
  MergeInputFiles_Called = true; // false the first time MergeInputFiles() is calles
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
      hToyMargCurr[is] = new TH1F(name,name,util::NumErecoBins(s),util::ErecoBinEdges(s));   
      cout << name << endl;
    }
    else{
      hToyMargCurr[is] = 0;
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

  hchisq_marg_1D = new TH1F("hchisq_marg_1D","hchisq_marg_1D",fnbins_param[0],fmin_param[0],fmax_param[0]+fstep_param[0]);  
  hchisq_marg_1D->GetXaxis()->SetTitle(Xtitle);
  hchisq_marg_1D->GetYaxis()->SetTitle("#Delta #chi^{2}");
  
  hchisq_marg_2D = new TH2F("hchisq_marg_2D","hchisq_marg_2D",
			    fnbins_param[0],fmin_param[0],fmax_param[0]+fstep_param[0],
			    fnbins_param[1],fmin_param[1],fmax_param[1]+fstep_param[1] );
  hchisq_marg_2D->GetXaxis()->SetTitle(Xtitle);
  hchisq_marg_2D->GetYaxis()->SetTitle(Ytitle);

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

  fmin_sub_param[0] = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(fmin_bin_sub_param[0]+1);
  fmax_sub_param[0] = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(fmax_bin_sub_param[0]+1);
  fmin_sub_param[1] = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(fmin_bin_sub_param[1]+1);
  fmax_sub_param[1] = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(fmax_bin_sub_param[1]+1);

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

  // Set the interval where to find the minimum (bin numbers)
  fFirstBin_X = 1;
  fLastBin_X = hchisq_marg_2D->GetNbinsX();
  fFirstBin_Y = 1;
  fLastBin_Y = hchisq_marg_2D->GetNbinsY();
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
  if (hToyMargCurr_RHCNumu              ) { delete hToyMargCurr_RHCNumu              ; hToyMargCurr_RHCNumu               = 0;}
  if (fToyMargErecoVArray_RHCNumu       ) { delete fToyMargErecoVArray_RHCNumu       ; fToyMargErecoVArray_RHCNumu        = 0;}
  if (fToyMargErecoVArray_RHCNue        ) { delete fToyMargErecoVArray_RHCNue        ; fToyMargErecoVArray_RHCNue         = 0;}
  if (fToyMargErecoVArray_Numu          ) { delete fToyMargErecoVArray_Numu          ; fToyMargErecoVArray_Numu           = 0;}
  if (fToyMargErecoVArray_Nue           ) { delete fToyMargErecoVArray_Nue           ; fToyMargErecoVArray_Nue            = 0;}
  if (hchisq_marg_1D                    ) { delete hchisq_marg_1D                    ; hchisq_marg_1D                     = 0;}
  if (hchisq_marg_2D                    ) { delete hchisq_marg_2D                    ; hchisq_marg_2D                     = 0;}
  if (hchisq_marg_1D_Step               ) { delete hchisq_marg_1D_Step               ; hchisq_marg_1D_Step                = 0;}
  if (hchisq_marg_2D_Step               ) { delete hchisq_marg_2D_Step               ; hchisq_marg_2D_Step                = 0;}
  if (fparm                             ) { delete fparm                             ; fparm                              = 0;}
  if (fSavedToyParamLists               ) { delete fSavedToyParamLists               ; fSavedToyParamLists                = 0;}
  
}
//________________________________________________________________________________________
void FitMarginalization::SetSubRange(int *submin,int *submax){
  
  // sub-range: scan all the range by default
  fmin_bin_sub_param[0] = submin[0];
  fmax_bin_sub_param[0] = submax[0];
  fmin_bin_sub_param[1] = submin[1];
  fmax_bin_sub_param[1] = submax[1];

  // default values
  if(fmin_bin_sub_param[0] == -999) fmin_bin_sub_param[0] = 0; 
  if(fmax_bin_sub_param[0] == -999) fmax_bin_sub_param[0] = fnbins_param[0] - 1; 
  if(fmin_bin_sub_param[1] == -999) fmin_bin_sub_param[1] = 0; 
  if(fmax_bin_sub_param[1] == -999) fmax_bin_sub_param[1] = fnbins_param[1] - 1; 

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
      (fnbins_sub_param[0] > fnbins_param[0]) ||
      (fmin_bin_sub_param[0] < 0) || 
      (fmax_bin_sub_param[0] > fnbins_param[0] - 1) )
    {
      cerr << endl;
      cerr << "Param 0 sub-range outside the entire free parameter range" << endl;
      cerr << "fnbins_param[0] = " << fnbins_param[0] << endl;
      cerr << "fnbins_sub_param[0] = " << fnbins_sub_param[0] << endl;
      cerr << "fmin_bin_sub_param[0] = " << fmin_bin_sub_param[0] << endl;
      cerr << "fmax_bin_sub_param[0] = " << fmax_bin_sub_param[0] << endl;
      cerr << endl;
      exit(1);
    }

  if( (fPhysParamFree[1] != conv::kPhypUndefined) &&
      (fnbins_sub_param[1] > fnbins_param[1]) ||
      (fmin_bin_sub_param[1] < 0) || 
      (fmax_bin_sub_param[1] > fnbins_param[1] - 1) )
    {
      cerr << endl;
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
    if(minpar[ParId]==-999 && maxpar[ParId]==-999){ // if limits are not given use the default ones
      fmin_param[ParId] = mixing13Bar_min;
      fmax_param[ParId] = mixing13Bar_max;
    }    
  }
  else if(fPhysParamFree[ParId]==kPhypMixing23Bar){ // mixing23Bar
    if(minpar[ParId]==-999 && maxpar[ParId]==-999){ // if limits are not given use the default ones
      fmin_param[ParId] = mixing23Bar_min;
      fmax_param[ParId] = mixing23Bar_max;
    } 
  }
  else if(fPhysParamFree[ParId]==kPhypDm23sqBar){ // Dm23sqBar
    if(minpar[ParId]==-999 && maxpar[ParId]==-999){ // if limits are not given use the default ones
      fmin_param[ParId] = dmsq32Bar_min;
      fmax_param[ParId] = dmsq32Bar_max;
    }    
  }
  else if(fPhysParamFree[ParId]==kPhypCPVBar){ // deltaCPBar
    if(minpar[ParId]==-999 && maxpar[ParId]==-999){ // if limits are not given use the default ones
      fmin_param[ParId] = deltacpBar_min;
      fmax_param[ParId] = deltacpBar_max;
    }    
  }
  else if(fPhysParamFree[ParId]==kPhypMixing12Bar){ // mixing12Bar
    if(minpar[ParId]==-999 && maxpar[ParId]==-999){ // if limits are not given use the default ones
      fmin_param[ParId] = mixing12Bar_min;
      fmax_param[ParId] = mixing12Bar_max;
    }    
  }
  else if(fPhysParamFree[ParId]==kPhypDm12sqBar){ // Dm12sqBar
    if(minpar[ParId]==-999 && maxpar[ParId]==-999){ // if limits are not given use the default ones
      fmin_param[ParId] = dmsq12Bar_min;
      fmax_param[ParId] = dmsq12Bar_max;
    }    
  }

  // Nu Osc parameters
  else if(fPhysParamFree[ParId]==kPhypMixing13){ // mixing13
    if(minpar[ParId]==-999 && maxpar[ParId]==-999){ // if limits are not given use the default ones
      fmin_param[ParId] = mixing13_min;
      fmax_param[ParId] = mixing13_max;
    }    
  }
  else if(fPhysParamFree[ParId]==kPhypMixing23){ // mixing23
    if(minpar[ParId]==-999 && maxpar[ParId]==-999){ // if limits are not given use the default ones
      fmin_param[ParId] = mixing23_min;
      fmax_param[ParId] = mixing23_max;
    }    
  }
  else if(fPhysParamFree[ParId]==kPhypDm23sq){ // Dm23sq
    if(minpar[ParId]==-999 && maxpar[ParId]==-999){ // if limits are not given use the default ones
      fmin_param[ParId] = dmsq32_min;
      fmax_param[ParId] = dmsq32_max;
    }    
  }
  else if(fPhysParamFree[ParId]==kPhypCPV){ // deltaCP
    if(minpar[ParId]==-999 && maxpar[ParId]==-999){ // if limits are not given use the default ones
      fmin_param[ParId] = deltacp_min;
      fmax_param[ParId] = deltacp_max;
    }    
  }
  else if(fPhysParamFree[ParId]==kPhypMixing12){ // mixing12
    if(minpar[ParId]==-999 && maxpar[ParId]==-999){ // if limits are not given use the default ones
      fmin_param[ParId] = mixing12_min;
      fmax_param[ParId] = mixing12_max;
    }    
  }
  else if(fPhysParamFree[ParId]==kPhypDm12sq){ // Dm12sq
    if(minpar[ParId]==-999 && maxpar[ParId]==-999){ // if limits are not given use the default ones
      fmin_param[ParId] = dmsq12_min;
      fmax_param[ParId] = dmsq12_max;
    } 
  }

  // The parameter is not selected
  else{
    if(ParId==0){
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
    cerr << endl
	 << "ParID = " << ParId << endl
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
    cerr << "Define correctly parameter 0 for 1D fit and parameter 1 for 2D fit" << endl;
    cerr << endl;
    exit(1);    
  }


  // Parameter 0
  if( fparm->IsFixed(fPhysParamFree[0]) ){
    cerr << "The free parameter " << fnamefreeparam[0] 
	 << " is fixed in FitParamList " << fit_param_list
	 << endl << endl;
    exit(1);
  }
  else if( fparm->IsNuisance(fPhysParamFree[0]) ){
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
    cerr << "The free parameter " << fnamefreeparam[1] 
	 << " is fixed in FitParamList " << fit_param_list
	 << endl << endl;
    exit(1);
  }
  else if( fparm->IsNuisance(fPhysParamFree[1]) ){
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
  cout << hchisq_marg_2D->GetXaxis()->GetBinLowEdge(fmin_bin_sub_param[0] + 1)
       << ","
       << hchisq_marg_2D->GetXaxis()->GetBinLowEdge(fmax_bin_sub_param[0] + 1);
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
  cout << hchisq_marg_2D->GetYaxis()->GetBinLowEdge(fmin_bin_sub_param[1] + 1)
       << ","
       << hchisq_marg_2D->GetYaxis()->GetBinLowEdge(fmax_bin_sub_param[1] + 1);
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
	    if(this->doMargMassHierarchy == false){
	      
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
	    fparm->Update(fPhysParamFree[0],currval0); 
	    if(fIsFit2D) fparm->Update(fPhysParamFree[1],currval1); 

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
		double bincont = obsErec[ib];
		fToyMargErecoVArray->SetEntry(bincont,iexperiment,ipar0-fmin_bin_sub_param[0],ipar1-fmin_bin_sub_param[1],iMH,iscount,ib); 
	      }
	      
	      iscount++; // index of the toy array
	    } // loop SK samples

	    
#ifdef __print_array_toy_maketoys__
	      cout << endl;
	      cout << "toy=" << iexperiment << " - ipar0 = " << ipar0 << " - ipar1 = " << ipar1 << " : ";
	      
	      // Build toy Ereco distributions

	      iscount=0;
	      //for(int is = 0; is < kns; is++) {
	      //SuperkSample_t s = kSKSamples[is];
	      //if(s == kSKsUndefined) continue;
	      for(int is=0;is<fNSKsamples;is++){
		SuperkSample_t s = this->GetSKSamples(is);
		
		cout << util::AsStringCompact(s) << "\t"<<skobs->ErecoSpectrumVArray(s)->Integral()<< endl;

#ifdef __print_array_toy_maketoys_verbose__
		for(int ib = 0; ib < util::NumErecoBins(s); ib++){
		  double bincont = skobs->ErecoSpectrumArray(s)[ib];
		  fToyMargErecoVArray->SetEntry(bincont,iexperiment,ipar0-fmin_bin_sub_param[0],ipar1-fmin_bin_sub_param[1],iMH,iscount,ib); 
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
void FitMarginalization::FillChisqMarg(SuperkObservation *skobs_data){
  
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
       << hchisq_marg_2D->GetXaxis()->GetBinLowEdge(fmin_bin_sub_param[0] + 1)
       << ","
       << hchisq_marg_2D->GetXaxis()->GetBinLowEdge(fmax_bin_sub_param[0] + 1)
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
       << hchisq_marg_2D->GetYaxis()->GetBinLowEdge(fmin_bin_sub_param[1] + 1)
       << ","
       << hchisq_marg_2D->GetYaxis()->GetBinLowEdge(fmax_bin_sub_param[1] + 1)
       << "]"
       << endl;


  // It can marginalize sub-ranges of the free parameter
  // fmin_bin_sub_param_1 and fmax_bin_sub_param_1 read from the binary file

  NparTot0=0;
  for(int ipar0=fmin_bin_sub_param[0];ipar0<=fmax_bin_sub_param[0];ipar0++){
    NparTot1=0;
    for(int ipar1=fmin_bin_sub_param[1];ipar1<=fmax_bin_sub_param[1];ipar1++){
      
      
      CalcChisqMargPt(skobs_data,ipar0,ipar1);
      

      NparTot1++;
    } // ipar1 (free parameter)
    NparTot0++;
  } // ipar0 (free parameter)



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
  for(int ipar0=fmin_bin_sub_param[0];ipar0<=fmax_bin_sub_param[0];ipar0++){
    for(int ipar1=fmin_bin_sub_param[1];ipar1<=fmax_bin_sub_param[1];ipar1++){
      cout << fChisqMargVArray->GetEntry(ipar0,ipar1) << ", ";
    }
  }
  cout << endl;
#endif

  return;
}
//________________________________________________________________________________________
float FitMarginalization::CalcChisqMargPt(SuperkObservation *skobs_data,int iterpar0,int iterpar1){
  
  // If chisq already colculated return the stored value and break
  float chisq = fChisqMargVArray->GetEntry(iterpar0,iterpar1); 
  if(chisq > 0.) return chisq;

  fNCallsCalcChisqMargPt++;

  //cout << "iterpar0 = " << iterpar0 << " - iterpar1 = " << iterpar1 << "   ";
  
  double currval0 = fmin_param[0] + iterpar0*fstep_param[0];
  double currval1 = fmin_param[1] + iterpar1*fstep_param[1];
  //cout << "currval0 = " << currval0 << endl;
  //cout << "currval1 = " << currval1 << endl;
  
  // Marginalize the Mass Hierarchy
  int nMassHierarchy = 2;
  
  for(int iMH=0;iMH<nMassHierarchy;iMH++){
    
    // Don't marginalize Mass Hierarchy. Only the assumed one
    if(this->doMargMassHierarchy == false){
      
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
      float chisq = 0.;
      
      // Loop over the used SK samples
      int iscount=0;
      //for(int is = 0; is < kns; is++) {
      //SuperkSample_t s = kSKSamples[is];
      //if(s == kSKsUndefined) continue;
      for(int is=0;is<fNSKsamples;is++){
	SuperkSample_t s = this->GetSKSamples(is);
	
	for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){
	  float nexp = fToyMargErecoVArray->GetEntry(iexperiment,iterpar0-fmin_bin_sub_param[0],iterpar1-fmin_bin_sub_param[1],iMH,iscount,ibin);
	  float nobs = skobs_data->ErecoSpectrumArray(s)[ibin];
	  //cout << "nexp = " << nexp << " - nobs = " << nobs << endl;
	  chisq += valor::util::ChisqPoisson(nobs,nexp); // penalty term is a simple constant --> not needed
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
	  cout << bincont << " "; 
	}
	cout << endl;
	cout << endl;
	
	iscount++;
      } // loop SK samples
      
#endif
      
      float like = exp(-chisq/2.) / fntoys_marg;
      
      if(this->doMargMassHierarchy) like /= 2.;  // if MH is marginalized, divide by 2
      
      flikeVArray->AddToEntry(like,iterpar0,iterpar1); // this is correct if binary files differ only for free param sub-ranges
      
    } // iexperiment
  } // iMH (mass hierarchy)
  
  // Now the marginalized likelihood is computed using toys from all the binary files
  // Compute the chisq for each free osc parameter value and fill the array and histo
  float bincont = -2*log( flikeVArray->GetEntry(iterpar0,iterpar1) );
  fChisqMargVArray->SetEntry(bincont,iterpar0,iterpar1);
  if(!fIsFit2D) hchisq_marg_1D->SetBinContent(iterpar0+1,bincont);
  else          hchisq_marg_2D->SetBinContent(iterpar0+1,iterpar1+1,bincont);
      
  return bincont;
}
//________________________________________________________________________________________
void FitMarginalization::WriteToys2Binary(int idfile, const char *tag){

  cout << endl;
  cout << "Write marginalization toys in a binary file:" << endl;
  
  char namemh[200];
  if(doMargMassHierarchy) sprintf(namemh,"%s","marg");
  else{ 
    if( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal )
      sprintf(namemh,"%s","NH");
    else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted )
      sprintf(namemh,"%s","IH");
    else{
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

  int min_bin_sub_param0 = fmin_bin_sub_param[0];
  int max_bin_sub_param0 = fmax_bin_sub_param[0];
  int min_bin_sub_param1 = fmin_bin_sub_param[1];
  int max_bin_sub_param1 = fmax_bin_sub_param[1];

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
  if( this->doMargMassHierarchy ) ostext << " Marginalize";
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
  if(doMargMassHierarchy) sprintf(namemh,"%s","marg");
  else{ 
    if( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal )
      sprintf(namemh,"%s","NH");
    else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted )
      sprintf(namemh,"%s","IH");
    else{
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

  int min_bin_sub_param0 = fmin_bin_sub_param[0];
  int max_bin_sub_param0 = fmax_bin_sub_param[0];
  int min_bin_sub_param1 = fmin_bin_sub_param[1];
  int max_bin_sub_param1 = fmax_bin_sub_param[1];

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
  if( this->doMargMassHierarchy ) ostext << " Marginalize";
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
  if(doMargMassHierarchy) sprintf(namemh,"%s","marg");
  else{ 
    if( Opt::Instance()->gMassHierarchyAssumed == kNmhNormal )
      sprintf(namemh,"%s","NH");
    else if( Opt::Instance()->gMassHierarchyAssumed == kNmhInverted )
      sprintf(namemh,"%s","IH");
    else{
      cerr << "Mass Hierarchy not defined" << endl;
      exit(1);
    }
  }	    

  char inname[200];
  sprintf(inname,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param[0],fnbins_param[1],namemh,"",".bin");

  std::string sinname(inname);
  cout << sinname << endl;
  
  cout << "Reading toys for marginalization from binary file: " << sinname << endl;  
  fToyMargErecoVArray->ReadFile(sinname);
  cout << "Done!"<<endl;
  
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

  cout << "Reading free param sub-range from binary file:" << sinname_range << endl;

  ifstream is(sinname_range.c_str(),ios::binary|ios::in);
  assert(is);
  is.read( (char*)&min_bin_sub_param0,sizeof(int));
  is.read( (char*)&max_bin_sub_param0,sizeof(int));
  is.read( (char*)&min_bin_sub_param1,sizeof(int));
  is.read( (char*)&max_bin_sub_param1,sizeof(int));
  is.close();
  
  cout<<"Done!"<<endl;

  fmin_bin_sub_param[0] = min_bin_sub_param0;
  fmax_bin_sub_param[0] = max_bin_sub_param0;
  fmin_bin_sub_param[1] = min_bin_sub_param1;
  fmax_bin_sub_param[1] = max_bin_sub_param1;

  fmin_sub_param[0] = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(fmin_bin_sub_param[0] + 1);    
  fmax_sub_param[0] = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(fmax_bin_sub_param[0] + 1);
  fmin_sub_param[1] = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(fmin_bin_sub_param[1] + 1);    
  fmax_sub_param[1] = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(fmax_bin_sub_param[1] + 1);
    
  if(fmin_sub_param[0]<(fmin_param[0]-conv::kASmallNum) || fmax_sub_param[0]>(fmax_param[0]+conv::kASmallNum)){
    cerr << endl 
	 << "The free parameter 0 sub-range ["
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
  
  if(fmin_sub_param[1]<(fmin_param[1]-conv::kASmallNum) || fmax_sub_param[1]>(fmax_param[1]+conv::kASmallNum)){  // w/o conv::kASmallNum it enters the loop if bounds are identical. To check
    cerr << endl 
	 << "The free parameter 1 sub-range ["
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
  
  cout << "File sub-range: ["
       << fmin_sub_param[0] << "," << fmax_sub_param[0] 
       << "]" << endl;
  cout << "File bin sub-range: ["
       << fmin_bin_sub_param[0] << "," << fmax_bin_sub_param[0] 
       << "]" << endl << endl;

  cout << "File sub-range: ["
       << fmin_sub_param[1] << "," << fmax_sub_param[1] 
       << "]" << endl;
  cout << "File bin sub-range: ["
       << fmin_bin_sub_param[1] << "," << fmax_bin_sub_param[1] 
       << "]" << endl << endl;
 
  return;
}
//________________________________________________________________________________________
void FitMarginalization::MinimizeHisto(){
  
  int minbin = -999;
  int minbinX = -999;
  int minbinY = -999;
  int minbinZ = -999;

  cout << endl;
  cout << "Find minimum of marginalized chisq" << endl;
  cout << endl;

  if(!fIsFit2D){
    fchisq_bf = hchisq_marg_1D->GetMinimum();
    minbin = hchisq_marg_1D->GetMinimumBin();
    fparam_bf[0] = hchisq_marg_1D->GetBinLowEdge(minbin);
    fparam_bf[1] = -999.;
  }
  else{
    fchisq_bf = hchisq_marg_2D->GetMinimum();
    hchisq_marg_2D->GetMinimumBin(minbinX,minbinY,minbinZ);
    fparam_bf[0] = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(minbinX);
    fparam_bf[1] = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(minbinY);
  }

  return;
}
//________________________________________________________________________________________
void FitMarginalization::MinimizeHisto(TH1F *hInChisqMarg){
  
  int minbin = -999;
  
  fchisq_bf = hInChisqMarg->GetMinimum();
  minbin = hInChisqMarg->GetMinimumBin();
  fparam_bf[0] = hInChisqMarg->GetBinLowEdge(minbin);
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::MinimizeHisto(TH2F *hInChisqMarg){
  
  int minbin = -999;
  
  fchisq_bf = hInChisqMarg->GetMinimum();
  minbin = hInChisqMarg->GetMinimumBin();
  fparam_bf[0] = hInChisqMarg->GetXaxis()->GetBinLowEdge(minbin);
  fparam_bf[1] = hInChisqMarg->GetYaxis()->GetBinLowEdge(minbin);
  
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
	double bincont = fToyMargErecoVArray->GetEntry(itoy,ifreepar0,ifreepar1,iMassHierarchy,is,ibin);
	hToyMargCurr[is]->SetBinContent(ibin+1,bincont);
	//cout << bincont << " ";
      }     
      //cout << endl;
      return hToyMargCurr[is];
    }
  }
  
  cerr << endl;
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
      hchisq_marg_1D->SetBinContent(ipar0+1,bincont);
      hchisq_marg_2D->SetBinContent(ipar0+1,ipar1+1,bincont);
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
    cout << "You chose a user prior but you didn't give any histogram as input" << endl;
    exit(1);
  }
  else{
    cout << "You chose the wrong prior option" << endl;
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
    cout << "You chose a " << priorname << " prior, but you gave as input an histogram" << endl;
    exit(1);
  }
  else{
    cout << "You chose the wrong prior option" << endl;
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

    bool doMarginalize = true;
    
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
    double currval = -999;

    bool rethrow=true;
    while(rethrow){
      currval = Randomize(p);
      rethrow = RethrowOscUnphysical(p,currval);
      if(rethrow) cout << "Rethrow " << util::AsString(p) << " " << currval << endl;
    }

    // Update current values of parameters 
    fparm->Update(p,currval); 

    cout << util::AsString(p) << " " << currval << endl;
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
      
      if ( ( fabs(mean -conv::kBadResult) < conv::kASmallNum ) || 
           ( fabs(fracerr-conv::kBadResult) < conv::kASmallNum ) ) {
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
    cout << "Used prior is: " << fparm->GetPriorModeName(physpar) << endl;
    cout << "Prior not well defined for " << util::AsString(physpar) << endl;
    exit(1);
  }
  
}
//________________________________________________________________________________________
void FitMarginalization::SetMargMassHierarchy(bool domargMH){
  doMargMassHierarchy = domargMH;
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
  if(MergeInputFiles_Called){
    for(int is=0;is<fNSKsamples;is++)
      fToyMargErecoVArray_Clone = new VALORArrayF(fntoys_marg, fnbins_param[0],fnbins_param[1], nmh, fNSKsamples, fnerecobins->GetEntry(is)); // total range of free parameter
    
    fToyMargErecoVArray_Clone->SetAll(-1.);

    MergeInputFiles_Called = false;
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
       << hchisq_marg_2D->GetXaxis()->GetBinLowEdge(fmin_bin_sub_param[0] + 1)
       << ","
       << hchisq_marg_2D->GetXaxis()->GetBinLowEdge(fmax_bin_sub_param[0] + 1)
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
       << hchisq_marg_2D->GetYaxis()->GetBinLowEdge(fmin_bin_sub_param[1] + 1)
       << ","
       << hchisq_marg_2D->GetYaxis()->GetBinLowEdge(fmax_bin_sub_param[1] + 1)
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
TH1F * FitMarginalization::GetChisqMarg   (void)
{
  return hchisq_marg_1D;
};
//________________________________________________________________________________________
TH2F * FitMarginalization::GetChisqMarg2D (void)
{
  return hchisq_marg_2D;
};
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








//________________________________________________________________________________________
//________________________________________________________________________________________
//________________________________________________________________________________________
//________________________________________________________________________________________
//________________________________________________________________________________________
//________________________________________________________________________________________










//________________________________________________________________________________________
void FitMarginalization::SetMinimInterval(double XvalFirst,double XvalLast,double YvalFirst,double YvalLast){

  fFirstVal_X = XvalFirst;
  fLastVal_X = XvalLast;;
  fFirstVal_Y = YvalFirst;
  fLastVal_Y = YvalLast;

  fFirstBin_X = hchisq_marg_2D->GetXaxis()->FindBin(fFirstVal_X);
  fLastBin_X  = hchisq_marg_2D->GetXaxis()->FindBin(fLastVal_X);
  fFirstBin_Y = hchisq_marg_2D->GetYaxis()->FindBin(fFirstVal_Y);
  fLastBin_Y  = hchisq_marg_2D->GetYaxis()->FindBin(fLastVal_Y);

  int nbins_X = fnbins_param[0];
  int nbins_Y = fnbins_param[1];

  // Don't exceed the allowed bins
  if(fFirstBin_X < 1){
    cerr << endl;
    cerr << "The first bin of the X interval is outside the histogram limit" << endl;
    cerr << "Set it to 1" << endl;
    cerr << endl;
    fFirstBin_X = 1;
    fFirstVal_X = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(fFirstBin_X);
  }
  if(fFirstBin_Y < 1){
    cerr << endl;
    cerr << "The first bin of the Y interval is outside the histogram limit" << endl;
    cerr << "Set it to 1" << endl;
    cerr << endl;
    fFirstBin_Y = 1;
    fFirstVal_Y = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(fFirstBin_Y);
  }
  if(fLastBin_X > nbins_X){
    cerr << endl;
    cerr << "The last bin of the X interval is outside the histogram limit" << endl;
    cerr << "Set it to the number of bins in the histogram:" << nbins_X << endl;
    cerr << endl;
    fLastBin_X = nbins_X;
    fLastVal_X = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(fLastBin_X);
  }
  if(fLastBin_Y > nbins_Y){
    cerr << endl;
    cerr << "The last bin of the Y interval is outside the histogram limit" << endl;
    cerr << "Set it to the number of bins in the histogram:" << nbins_Y << endl;
    cerr << endl;
    fLastBin_Y = nbins_Y;
    fLastVal_Y = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(fLastBin_Y);
  }

  cout << endl;
  cout << "Minimization interval set to:" << endl;
  cout << "  - X: [" << fFirstVal_X << "," << fLastVal_X << "]" << endl;
  cout << "  - X Bin: [" << fFirstBin_X << "," << fLastBin_X << "]" << endl;
  cout << "  - Y: [" << fFirstVal_Y << "," << fLastVal_Y << "]" << endl;
  cout << "  - Y Bin: [" << fFirstBin_Y << "," << fLastBin_Y << "]" << endl;
  cout << endl;

}
//________________________________________________________________________________________
void FitMarginalization::SetMinimUncorrXY(bool flaguncorr){ 
  fUncorrXY = flaguncorr;
}
//________________________________________________________________________________________
void FitMarginalization::MinimizeChisqInterval(SuperkObservation *skobs_data){

  float bigvalueUOflow = 99999999.; 

  // 1D
  hchisq_marg_1D_Step = (TH1F*)hchisq_marg_1D->Clone("hchisq_marg_1D_Step");
  for(int ibinX=0;ibinX<hchisq_marg_1D_Step->GetNbinsX();ibinX++){
    hchisq_marg_1D_Step->SetBinContent(ibinX+1,0.);
  }
  // Set under and over flow of X to very high value
  for(int ibinX=0;ibinX<hchisq_marg_1D->GetNbinsX();ibinX++){    
    hchisq_marg_1D->SetBinContent(0,bigvalueUOflow);
    hchisq_marg_1D->SetBinContent(hchisq_marg_2D->GetNbinsX()+1,bigvalueUOflow);
  }

  // 2D
  hchisq_marg_2D_Step = (TH2F*)hchisq_marg_2D->Clone("hchisq_marg_2D_Step");
  for(int ibinX=0;ibinX<hchisq_marg_2D_Step->GetNbinsX();ibinX++){
    for(int ibinY=0;ibinY<hchisq_marg_2D_Step->GetNbinsY();ibinY++){
      hchisq_marg_2D_Step->SetBinContent(ibinX+1,ibinY+1,0.);
    }
  }
  // Set under and over flow of X and Y to very high value
  for(int ibinX=0;ibinX<hchisq_marg_2D->GetNbinsX();ibinX++){    
    hchisq_marg_2D->SetBinContent(ibinX+1,0,bigvalueUOflow);
    hchisq_marg_2D->SetBinContent(ibinX+1,hchisq_marg_2D->GetNbinsY()+1,bigvalueUOflow);
  }
  for(int ibinY=0;ibinY<hchisq_marg_2D->GetNbinsY();ibinY++){    
    hchisq_marg_2D->SetBinContent(0,ibinY+1,bigvalueUOflow);
    hchisq_marg_2D->SetBinContent(hchisq_marg_2D->GetNbinsX()+1,ibinY+1,bigvalueUOflow);
  }
  

  TRandom3 *r = new TRandom3(0); // needed to randomly scan the distribution
  
  // Get the interval
  
  double UniFirstBin_X = fFirstBin_X; 
  double UniLastBin_X  = fLastBin_X; 
  double UniFirstBin_Y = fFirstBin_Y; 
  double UniLastBin_Y  = fLastBin_Y; 
  
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
    
    int binN_X = -999;
    int binN_Y = -999;

    // same point
    binN_X = bin_X; 
    if(fIsFit2D) binN_Y = bin_Y;
    
    double binlowedge_X = -999;
    double binlowedge_Y = -999;

    // Find direction along X axis     
    
    // L(N,N)
    if(fIsFit2D) 
      bincont = fChisqMargVArray->GetEntry(binN_X-1,binN_Y-1);
    else 
      bincont = fChisqMargVArray->GetEntry(binN_X-1,fmin_bin_sub_param[1]);
    if( bincont<=0.){
      if(fIsFit2D) bincont = CalcChisqMargPt(skobs_data,binN_X-1,binN_Y-1); // CalcChisqMargPt takes the point not bin
      else bincont = CalcChisqMargPt(skobs_data,binN_X-1,fmin_bin_sub_param[1]);
    }
    chisq_X[ilike_N] = bincont;
    chisq[ilike_N][ilike_N] = bincont;

    // L(N-1,N)
    int binNm1_X = binN_X - 1;
    if(binNm1_X < fFirstBin_X) 
      binNm1_X = binN_X;  // should happens only when binN = binF
    if(fIsFit2D) 
      bincont = fChisqMargVArray->GetEntry(binNm1_X-1,binN_Y-1);
    else 
      bincont = fChisqMargVArray->GetEntry(binNm1_X-1,fmin_bin_sub_param[1]);
    if( bincont<=0.){
      if(fIsFit2D) bincont = CalcChisqMargPt(skobs_data,binNm1_X-1,binN_Y-1);
      else bincont = CalcChisqMargPt(skobs_data,binNm1_X-1,fmin_bin_sub_param[1]);
    }
    chisq_X[ilike_Nm1] = bincont;
    chisq[ilike_Nm1][ilike_N] = bincont;
    
    // L(N+1,N)
    int binNp1_X = binN_X + 1;
    if(binNp1_X > fLastBin_X) 
      binNp1_X = binN_X;  // should happens only when binN = binL
    if(fIsFit2D) 
      bincont = fChisqMargVArray->GetEntry(binNp1_X-1,binN_Y-1);
    else 
      bincont = fChisqMargVArray->GetEntry(binNp1_X-1,fmin_bin_sub_param[1]);
    if( bincont<=0.){
      if(fIsFit2D) bincont = CalcChisqMargPt(skobs_data,binNp1_X-1,binN_Y-1);
      else bincont = CalcChisqMargPt(skobs_data,binNp1_X-1,fmin_bin_sub_param[1]); // bin_y = 0
    }
    chisq_X[ilike_Np1] = bincont;
    chisq[ilike_Np1][ilike_N] = bincont;
    
    binlowedge_X = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(bin_X);

        
    // Find direction along Y axis
    
    int binNm1_Y = -999;
    int binNp1_Y = -999;
    
    if(fIsFit2D){
      // L(N,N) already computed for 2D chisq
      bincont = fChisqMargVArray->GetEntry(binN_X-1,binN_Y-1);
      if( bincont<=0.)
	bincont = CalcChisqMargPt(skobs_data,binN_X-1,binN_Y-1);
      chisq_Y[ilike_N] = bincont;
      
      // L(N,N-1)
      binNm1_Y = binN_Y - 1;
      if(binNm1_Y < fFirstBin_Y) binNm1_Y = binN_Y;  // should happens only when binN = binF
      bincont = fChisqMargVArray->GetEntry(binN_X-1,binNm1_Y-1);
      if( bincont<=0.)
	bincont = CalcChisqMargPt(skobs_data,binN_X-1,binNm1_Y-1);
      chisq_Y[ilike_Nm1] = bincont;
      chisq[ilike_N][ilike_Nm1] = bincont;
      
      // L(N,N+1)
      binNp1_Y = binN_Y + 1;
      if(binNp1_Y > fLastBin_Y) binNp1_Y = binN_Y;  // should happens only when binL = binL
      bincont = fChisqMargVArray->GetEntry(binN_X-1,binNp1_Y-1);
      if( bincont<=0.)
	bincont = CalcChisqMargPt(skobs_data,binN_X-1,binNp1_Y-1);
      chisq_Y[ilike_Np1] = bincont;
      chisq[ilike_N][ilike_Np1] = bincont;
      
      binlowedge_Y = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(bin_Y);
    }




    // Find off-diagonal directions     
    
    if(fIsFit2D){
      
      // L(N-1,N-1)
      bincont = fChisqMargVArray->GetEntry(binNm1_X-1,binNm1_Y-1);
      if( bincont<=0.)
     	bincont = CalcChisqMargPt(skobs_data,binNm1_X-1,binNm1_Y-1);
      chisq[ilike_Nm1][ilike_Nm1] = bincont;
      
      // L(N+1,N+1)
      bincont = fChisqMargVArray->GetEntry(binNp1_X-1,binNp1_Y-1);
      if( bincont<=0.)
     	bincont = CalcChisqMargPt(skobs_data,binNp1_X-1,binNp1_Y-1);
      chisq[ilike_Np1][ilike_Np1] = bincont;

      // L(N-1,N+1)
      bincont = fChisqMargVArray->GetEntry(binNm1_X-1,binNp1_Y-1);
      if( bincont<=0.)
     	bincont = CalcChisqMargPt(skobs_data,binNm1_X-1,binNp1_Y-1);
      chisq[ilike_Nm1][ilike_Np1] = bincont;

      // L(N+1,N-1)
      bincont = fChisqMargVArray->GetEntry(binNp1_X-1,binNm1_Y-1);
      if( bincont<=0.)
     	bincont = CalcChisqMargPt(skobs_data,binNp1_X-1,binNm1_Y-1);
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

    double bincent_X = -999;
    double bincent_Y = -999;

    if(!fIsFit2D){
      double bincent_X = hchisq_marg_1D_Step->GetXaxis()->GetBinCenter(bin_X);
      hchisq_marg_1D_Step->Fill(bincent_X);
    }
    else{
      double bincent_X = hchisq_marg_2D_Step->GetXaxis()->GetBinCenter(bin_X);
      double bincent_Y = hchisq_marg_2D_Step->GetYaxis()->GetBinCenter(bin_Y);
      hchisq_marg_2D_Step->Fill(bincent_X,bincent_Y);
    }
    


    // Best fit bin is found
    // this condition contains the case with the boundary, when
    // z(N-1) = Z(N) or Z(N+1) = Z(N)
 
    // if( 
    //    chisq_X[ilike_N] <= chisq_X[ilike_Np1] && 
    //    chisq_X[ilike_N] <= chisq_X[ilike_Nm1] && 
    //    ( !fIsFit2D ||
    // 	 ( chisq_Y[ilike_N] <= chisq_Y[ilike_Np1] && 
    // 	   chisq_Y[ilike_N] <= chisq_Y[ilike_Nm1] )
    // 	 )	
    // 	){


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
	cout << "Best fit X = " << hchisq_marg_2D->GetXaxis()->GetBinLowEdge(binN_X) << endl;
	cout << "Best fit Y = " << hchisq_marg_2D->GetYaxis()->GetBinLowEdge(binN_Y) << endl;
	cout << "Total # of iterations = " << Niterations << endl;
	cout << "Total # of performed fits = " << Nfits << endl;
	cout << "Check if it's a real minimum" << endl;

	int binNm2_X = -999;
	double chisq_X_Nm2 = -999;
	int binNp2_X = -999;
	double chisq_X_Np2 = -999;
	int binNm2_Y = -999;
	double chisq_Y_Nm2 = -999;
	int binNp2_Y = -999;
	double chisq_Y_Np2 = -999;


	// L(N-2) X
	binNm2_X = binN_X - 2;
	if(binNm2_X < fFirstBin_X) 
	  binNm2_X = fFirstBin_X;  // should happens only when Xn = Xf
	if(fIsFit2D) 
	  bincont = fChisqMargVArray->GetEntry(binNm2_X-1,binN_Y-1);
	else 
	  bincont = fChisqMargVArray->GetEntry(binNm2_X-1,fmin_bin_sub_param[1]);
	if( bincont<=0.)
	  if(fIsFit2D) bincont = CalcChisqMargPt(skobs_data,binNm2_X-1,binN_Y-1);
	  else bincont = CalcChisqMargPt(skobs_data,binNm2_X-1,fmin_bin_sub_param[1]);
	chisq_X_Nm2 = bincont;
	
	// L(N+2) X
	binNp2_X = binN_X + 2;
	if(binNp2_X > fLastBin_X) 
	  binNp2_X = fLastBin_X;  // should happens only when Xn = Xl
	if(fIsFit2D)
	  bincont = fChisqMargVArray->GetEntry(binNp2_X-1,binN_Y-1);
	else 
	  bincont = fChisqMargVArray->GetEntry(binNp2_X-1,fmin_bin_sub_param[1]);
	if( bincont<=0.)
	  if(fIsFit2D) bincont = CalcChisqMargPt(skobs_data,binNp2_X-1,binN_Y-1);
	  else bincont = CalcChisqMargPt(skobs_data,binNp2_X-1,fmin_bin_sub_param[1]);
	chisq_X_Np2 = bincont;

	if(fIsFit2D){
	  
	  // L(N-2) Y
	  binNm2_Y = binN_Y - 2;
	  if(binNm2_Y < fFirstBin_Y) 
	    binNm2_Y = fFirstBin_Y;  // should happens only when Xn = Xf
	  bincont = fChisqMargVArray->GetEntry(binN_X-1,binNm2_Y-1);
	  if( bincont<=0.)
	    bincont = CalcChisqMargPt(skobs_data,binN_X-1,binNm2_Y-1);
	  chisq_Y_Nm2 = bincont;
	  
	  // L(N+2) Y
	  binNp2_Y = binN_Y + 2;
	  if(binNp2_Y > fLastBin_Y) binNp2_Y = fLastBin_Y;  // should happens only when Xn = Xl
	  bincont = fChisqMargVArray->GetEntry(binN_X-1,binNp2_Y-1);
	  if( bincont<=0.)
	    bincont = CalcChisqMargPt(skobs_data,binN_X-1,binNp2_Y-1);
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
	      fparam_bf[0] = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(binN_X);
	      fparam_bf[1] = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(binN_Y);
	    }
	    else{
	      fparam_bf[0] = hchisq_marg_1D->GetBinLowEdge(binN_X);
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
	    else if(binN_X == fLastBin_X){ // on the first bin X
	      // next fit will not account the last bin X for the minimization
	      if( (UniLastBin_X-1) < UniFirstBin_X)
		cout << "Don't exclude the border, otherwise LastBin_X < FirstBin_X" << endl;
	      else
		UniLastBin_X = UniLastBin_X-1;
	      binL_X = UniLastBin_X;
	    }

	    if( !fIsFit2D || (binN_Y == fFirstBin_Y) ){ // on the first bin Y	    
	      // next fit will not account the first bin Y for the minimization
	      if( (UniFirstBin_Y+1) > UniLastBin_Y)
		cout << "Don't exclude the border, otherwise FirstBin_Y > LastBin_Y" << endl;
	      else
		UniFirstBin_Y = UniFirstBin_Y+1;
	      binF_Y = UniFirstBin_Y;
	    }
	    else if( !fIsFit2D || (binN_Y == fLastBin_Y) ){ // on the first bin Y
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
      
      int offdiagbin_X = -999; // X off-diagonal bin in the gradient direction
      int offdiagbin_Y = -999; // Y off-diagonal bin in the gradient direction

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
  
  cout << "The fit is well performed" << endl;

  return;
}
//________________________________________________________________________________________
void FitMarginalization::MinimizeChisqStep(SuperkObservation *skobs_data){

  GridStartBin(skobs_data); // Get the starting bin
  
  float bigvalueUOflow = 99999999.; 

  // 1D
  hchisq_marg_1D_Step = (TH1F*)hchisq_marg_1D->Clone("hchisq_marg_1D_Step");
  for(int ibinX=0;ibinX<hchisq_marg_1D_Step->GetNbinsX();ibinX++){
    hchisq_marg_1D_Step->SetBinContent(ibinX+1,0.);
  }
  // Set under and over flow of X to very high value
  for(int ibinX=0;ibinX<hchisq_marg_1D->GetNbinsX();ibinX++){    
    hchisq_marg_1D->SetBinContent(0,bigvalueUOflow);
    hchisq_marg_1D->SetBinContent(hchisq_marg_2D->GetNbinsX()+1,bigvalueUOflow);
  }

  // 2D
  hchisq_marg_2D_Step = (TH2F*)hchisq_marg_2D->Clone("hchisq_marg_2D_Step");
  for(int ibinX=0;ibinX<hchisq_marg_2D_Step->GetNbinsX();ibinX++){
    for(int ibinY=0;ibinY<hchisq_marg_2D_Step->GetNbinsY();ibinY++){
      hchisq_marg_2D_Step->SetBinContent(ibinX+1,ibinY+1,0.);
    }
  }
  // Set under and over flow of X and Y to very high value
  for(int ibinX=0;ibinX<hchisq_marg_2D->GetNbinsX();ibinX++){    
    hchisq_marg_2D->SetBinContent(ibinX+1,0,bigvalueUOflow);
    hchisq_marg_2D->SetBinContent(ibinX+1,hchisq_marg_2D->GetNbinsY()+1,bigvalueUOflow);
  }
  for(int ibinY=0;ibinY<hchisq_marg_2D->GetNbinsY();ibinY++){    
    hchisq_marg_2D->SetBinContent(0,ibinY+1,bigvalueUOflow);
    hchisq_marg_2D->SetBinContent(hchisq_marg_2D->GetNbinsX()+1,ibinY+1,bigvalueUOflow);
  }

  //
  
  TRandom3 *r = new TRandom3(0); // needed to randomly scan the distribution

  double UniFirstBin_X = fFirstBin_X; 
  double UniLastBin_X  = fLastBin_X; 
  double UniFirstBin_Y = fFirstBin_Y; 
  double UniLastBin_Y  = fLastBin_Y; 

  int binN_X = fStartBin_X;
  int binN_Y = fStartBin_Y;

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

    double binlowedge_X = -999;
    double binlowedge_Y = -999;
    
    // L(N,N)
    if(fIsFit2D) 
      bincont = fChisqMargVArray->GetEntry(binN_X-1,binN_Y-1);
    else 
      bincont = fChisqMargVArray->GetEntry(binN_X-1,fmin_bin_sub_param[1]);
    if( bincont<=0.){
      if(fIsFit2D) bincont = CalcChisqMargPt(skobs_data,binN_X-1,binN_Y-1); // CalcChisqMargPt takes the point not bin
      else bincont = CalcChisqMargPt(skobs_data,binN_X-1,fmin_bin_sub_param[1]);
    }
    chisq_X[ilike_N] = bincont;
    chisq[ilike_N][ilike_N] = bincont;

    // L(N-1,N)
    int binNm1_X = binN_X - 1;
    // Bin is on the border, it should happens only when binN = binF
    if(binNm1_X < fFirstBin_X){ 
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
	if(fIsFit2D) bincont = CalcChisqMargPt(skobs_data,binNm1_X-1,binN_Y-1);
	else bincont = CalcChisqMargPt(skobs_data,binNm1_X-1,fmin_bin_sub_param[1]);
      }
      chisq_X[ilike_Nm1] = bincont;
      chisq[ilike_Nm1][ilike_N] = bincont;
    }

    // L(N+1,N)
    int binNp1_X = binN_X + 1;
    // Bin is on the border, it should happens only when binN = binF
    if(binNp1_X > fLastBin_X){
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
	if(fIsFit2D) bincont = CalcChisqMargPt(skobs_data,binNp1_X-1,binN_Y-1);
	else bincont = CalcChisqMargPt(skobs_data,binNp1_X-1,fmin_bin_sub_param[1]); // bin_y = 0
      }
      chisq_X[ilike_Np1] = bincont;
      chisq[ilike_Np1][ilike_N] = bincont;
    }
    
    binlowedge_X = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(binN_X);
  
        
    // Find direction along Y axis
    
    int binNm1_Y = -999;
    int binNp1_Y = -999;
    
    if(fIsFit2D){
      
      // L(N,N) already computed for 2D chisq
      bincont = fChisqMargVArray->GetEntry(binN_X-1,binN_Y-1);
      if( bincont<=0.)
	bincont = CalcChisqMargPt(skobs_data,binN_X-1,binN_Y-1);
      chisq_Y[ilike_N] = bincont;
      
      // L(N,N-1)
      binNm1_Y = binN_Y - 1;
      // Bin is on the border, it should happens only when binN = binF
      if(binNm1_Y < fFirstBin_Y){
	//binNm1_Y = binN_Y;  
 	chisq_Y[ilike_Nm1] = conv::kALargeNum;
	chisq[ilike_N][ilike_Nm1] = conv::kALargeNum;
      }
      // normal situation
      else{
	bincont = fChisqMargVArray->GetEntry(binN_X-1,binNm1_Y-1);
	if( bincont<=0.)
	  bincont = CalcChisqMargPt(skobs_data,binN_X-1,binNm1_Y-1);
	chisq_Y[ilike_Nm1] = bincont;
	chisq[ilike_N][ilike_Nm1] = bincont;
      }
      
      // L(N,N+1)
      binNp1_Y = binN_Y + 1;
      // Bin is on the border, it should happens only when binN = binF
      if(binNp1_Y > fLastBin_Y){
	//binNp1_Y = binN_Y;  
	chisq_Y[ilike_Np1] = conv::kALargeNum;
	chisq[ilike_N][ilike_Np1] = conv::kALargeNum;
       }
      // normal situation
      else{
	bincont = fChisqMargVArray->GetEntry(binN_X-1,binNp1_Y-1);
	if( bincont<=0.)
	  bincont = CalcChisqMargPt(skobs_data,binN_X-1,binNp1_Y-1);
	chisq_Y[ilike_Np1] = bincont;
	chisq[ilike_N][ilike_Np1] = bincont;
      }

      binlowedge_Y = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(binN_Y);
    }


    // Calculate off-diagonal bins     
    
    if(fIsFit2D){
      
      // L(N-1,N-1)
      if(binNm1_X < fFirstBin_X || binNm1_Y < fFirstBin_Y){ 
	chisq[ilike_Nm1][ilike_Nm1] = conv::kALargeNum;
      }
      else{
	bincont = fChisqMargVArray->GetEntry(binNm1_X-1,binNm1_Y-1);
	if( bincont<=0.)
	  bincont = CalcChisqMargPt(skobs_data,binNm1_X-1,binNm1_Y-1);
	chisq[ilike_Nm1][ilike_Nm1] = bincont;
      }

      // L(N+1,N+1)
      if(binNp1_X > fLastBin_X || binNp1_Y > fLastBin_Y){ 
	chisq[ilike_Np1][ilike_Np1] = conv::kALargeNum;
      }
      else{
	bincont = fChisqMargVArray->GetEntry(binNp1_X-1,binNp1_Y-1);
	if( bincont<=0.)
	  bincont = CalcChisqMargPt(skobs_data,binNp1_X-1,binNp1_Y-1);
	chisq[ilike_Np1][ilike_Np1] = bincont;
      }

      // L(N-1,N+1)
      if(binNm1_X < fFirstBin_X || binNp1_Y > fLastBin_Y){ 
	chisq[ilike_Nm1][ilike_Np1] = conv::kALargeNum;
      }
      else{
	bincont = fChisqMargVArray->GetEntry(binNm1_X-1,binNp1_Y-1);
	if( bincont<=0.)
	  bincont = CalcChisqMargPt(skobs_data,binNm1_X-1,binNp1_Y-1);
	chisq[ilike_Nm1][ilike_Np1] = bincont;
      }

      // L(N+1,N-1)
      if(binNp1_X > fLastBin_X || binNm1_Y < fFirstBin_Y){ 
	chisq[ilike_Np1][ilike_Nm1] = conv::kALargeNum;
      }
      else{
	bincont = fChisqMargVArray->GetEntry(binNp1_X-1,binNm1_Y-1);
	if( bincont<=0.)
	  bincont = CalcChisqMargPt(skobs_data,binNp1_X-1,binNm1_Y-1);
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
      for(int i=0;i<3;i++) cout << chisq_X[i] << " ";      
    }
    else{
      cout << "Bin: [ " << binN_X << " , " << binN_Y << " ] " 
	   << " --> [ " << binlowedge_X << " , " << binlowedge_Y << " ] " 
	   << hchisq_marg_2D->GetBinContent(binN_X,binN_Y)
	   << endl;
      cout << "Chisq grid around binN: " << endl; 
      for(int i=2;i>=0;i--){
	for(int j=0;j<3;j++){
	  cout << chisq[j][i] << " ";
	}
	cout << endl;
      }
      cout << endl;
    }
#endif
    
    double bincent_X = -999;
    double bincent_Y = -999;

    if(!fIsFit2D){
      bincent_X = hchisq_marg_1D_Step->GetXaxis()->GetBinCenter(binN_X);
      hchisq_marg_1D_Step->Fill(bincent_X);
    }
    else{
      bincent_X = hchisq_marg_2D_Step->GetXaxis()->GetBinCenter(binN_X);
      bincent_Y = hchisq_marg_2D_Step->GetYaxis()->GetBinCenter(binN_Y);
      hchisq_marg_2D_Step->Fill(bincent_X,bincent_Y);
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
      cout << "Best fit X = " << hchisq_marg_2D->GetXaxis()->GetBinLowEdge(binN_X) << endl;
      cout << "Best fit Y = " << hchisq_marg_2D->GetYaxis()->GetBinLowEdge(binN_Y) << endl;
      cout << "Total # of iterations = " << Niterations << endl;
      cout << "Total # of performed fits = " << Nfits << endl;
      cout << "Check if it's a real minimum" << endl;
      
      int binNm2_X = -999;
      double chisq_X_Nm2 = -999;
      int binNp2_X = -999;
      double chisq_X_Np2 = -999;
      int binNm2_Y = -999;
      double chisq_Y_Nm2 = -999;
      int binNp2_Y = -999;
      double chisq_Y_Np2 = -999;
      
      // L(N-2) X
      binNm2_X = binN_X - 2;
      // Bin is on the border, it should happens only when binN = binF
      if(binNm2_X < fFirstBin_X){
	//binNm2_X = fFirstBin_X;  
	chisq_X_Nm2 = conv::kALargeNum;
      }
      // normal situation
      else{
	if(fIsFit2D) 
	  bincont = fChisqMargVArray->GetEntry(binNm2_X-1,binN_Y-1);
	else 
	  bincont = fChisqMargVArray->GetEntry(binNm2_X-1,fmin_bin_sub_param[1]);
	if( bincont<=0.)
	  if(fIsFit2D) bincont = CalcChisqMargPt(skobs_data,binNm2_X-1,binN_Y-1);
	  else bincont = CalcChisqMargPt(skobs_data,binNm2_X-1,fmin_bin_sub_param[1]);
	chisq_X_Nm2 = bincont;
      }
      
      // L(N+2) X
      binNp2_X = binN_X + 2;
      // Bin is on the border, it should happens only when binN = binF
      if(binNp2_X > fLastBin_X){
	//binNp2_X = fLastBin_X; 
	chisq_X_Np2 = conv::kALargeNum;
      }
      // normal situation
      else{
	if(fIsFit2D)
	  bincont = fChisqMargVArray->GetEntry(binNp2_X-1,binN_Y-1);
	else 
	  bincont = fChisqMargVArray->GetEntry(binNp2_X-1,fmin_bin_sub_param[1]);
	if( bincont<=0.)
	  if(fIsFit2D) bincont = CalcChisqMargPt(skobs_data,binNp2_X-1,binN_Y-1);
	  else bincont = CalcChisqMargPt(skobs_data,binNp2_X-1,fmin_bin_sub_param[1]);
	chisq_X_Np2 = bincont;
      }

      if(fIsFit2D){
	
	// L(N-2) Y
	binNm2_Y = binN_Y - 2;
	// Bin is on the border, it should happens only when binN = binF
	if(binNm2_Y < fFirstBin_Y){ 
	  //binNm2_Y = fFirstBin_Y;  
	  chisq_Y_Nm2 = conv::kALargeNum;
	}
	// normal situation
	else{
	  bincont = fChisqMargVArray->GetEntry(binN_X-1,binNm2_Y-1);
	  if( bincont<=0.)
	    bincont = CalcChisqMargPt(skobs_data,binN_X-1,binNm2_Y-1);
	  chisq_Y_Nm2 = bincont;
	}
	
	// L(N+2) Y
	binNp2_Y = binN_Y + 2;
	// Bin is on the border, it should happens only when binN = binF
	if(binNp2_Y > fLastBin_Y){
	  //binNp2_Y = fLastBin_Y;  
	  chisq_Y_Np2 = conv::kALargeNum; 
	}
	else{
	  bincont = fChisqMargVArray->GetEntry(binN_X-1,binNp2_Y-1);
	  if( bincont<=0.)
	    bincont = CalcChisqMargPt(skobs_data,binN_X-1,binNp2_Y-1);
	  chisq_Y_Np2 = bincont;
	}
      
      } // if fIsFit2D
   
#if defined __print_minimization_steps__
      cout << "chisq_X = " 
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

      // Check if it's a real minimum
      
      if( chisq_X[ilike_N] <= chisq_X_Nm2 && 
      	  chisq_X[ilike_N] <= chisq_X_Np2 &&
      	  (!fIsFit2D ||
      	   (chisq_Y[ilike_N] <= chisq_Y_Nm2 && 
      	    chisq_Y[ilike_N] <= chisq_Y_Np2)
      	   )
      	  )
      	{

	  Nminima++;
	  
	  // It's a real minimum
	  cout << "It's a real minimum" << endl;
	  
	  fchisq_bf = chisq_X[ilike_N];
	  if(fIsFit2D) fchisq_bf = chisq[ilike_N][ilike_N];
	  
	  if(fIsFit2D){
	    fparam_bf[0] = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(binN_X);
	    fparam_bf[1] = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(binN_Y);
	  }
	  else{
	    fparam_bf[0] = hchisq_marg_1D->GetBinLowEdge(binN_X);
	  }

	  Nminima++;

	  break;
	} // global minimum
      

      // It's a bin fluctuation

      else{
	cout << endl;
	cout << "The reached minimum is a fluctuation of the bin" << endl; 
	cout << "Restart the fit from a new starting point:" << endl;
	cout << endl;
	
	// Throw the bin number
	double r_X = r->Uniform(UniFirstBin_X,UniLastBin_X+1); // +1 because it's then truncated into a int. Extract 21th bin as well
	double r_Y = r->Uniform(UniFirstBin_Y,UniLastBin_Y+1); // +1 because it's then truncated into a int. Extract 21th bin as well
	binN_X = (int)r_X; 
	binN_Y = (int)r_Y;
	
	cout << "Restart the fit from a new starting point:" << endl;
	cout << "binN_X = " << binN_X << " - binN_Y = " << binN_Y << endl;
	cout << "X = " << hchisq_marg_2D->GetXaxis()->GetBinLowEdge(binN_X);
	if(fIsFit2D)
	  cout << " - Y = " << hchisq_marg_2D->GetYaxis()->GetBinLowEdge(binN_Y) << endl;
	cout << endl;
	
       	continue;
	
      } // if it's a real minimum or a bin fluctuation
    
    } // if chisq_X<=chisq_X[Np1] && chisq_Y<=chisq_Y[Np1] ... 

    
    
    // Get the new direction 

    else {
      
      int offdiagbin_X = -999; // X off-diagonal bin in the gradient direction
      int offdiagbin_Y = -999; // Y off-diagonal bin in the gradient direction
      
      // Move toward the lowest adjacent Chisq bin
      
      int NewDir_X = -999; // It's not the bin number but an array index
      int NewDir_Y = -999; // It's not the bin number but an array index
      float min = 999999;

      if(!fIsFit2D){
	for(int i=0;i<3;i++){
	  if( chisq_X[i] < min ){
	    NewDir_X = i;
	    min = chisq_X[i];
	  }
	}
      }
      else{
	for(int i=0;i<3;i++){
	  for(int j=0;j<3;j++){
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
      cout << endl;
      
      if(fIsFit2D){
	cout << "  - Y: ";
	if(NewDir_Y == ilike_Nm1) cout << "down"; 
	else if(NewDir_Y == ilike_Np1) cout << "up"; 
	else if(NewDir_Y == ilike_N) cout << "stop";
	cout << endl;
      }
#endif
      
      if(
	 (!fIsFit2D && NewDir_X == ilike_N) ||
	 (fIsFit2D && NewDir_X == ilike_N && NewDir_Y == ilike_N)
	 ){
	cerr << endl;
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
	double NextVal_X = hchisq_marg_1D->GetXaxis()->GetBinLowEdge(NextBin_X);
	double BinVal_X = binlowedge_X; 
	double NewVal_X = hchisq_marg_1D->GetXaxis()->GetBinLowEdge(NextBin_X);

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
	
	double NextVal_X = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(NextBin_X);
	double NextVal_Y = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(NextBin_Y);
	
	double BinVal_X = binlowedge_X; 
	double BinVal_Y = binlowedge_Y;
	
	double NewVal_X = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(NextBin_X);
	double NewVal_Y = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(NextBin_Y);
		
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
void FitMarginalization::GridStartBin(SuperkObservation *skobs_data)
{
  
  // If MCMC is used we need to set the starting point before
  
  const int NThrowsExtrBin = fNumStartBin_X*fNumStartBin_Y;
  
  double UniFirstBin_X = fFirstBin_X; 
  double UniLastBin_X  = fLastBin_X; 
  double UniFirstBin_Y = fFirstBin_Y; 
  double UniLastBin_Y  = fLastBin_Y; 
  
  double bincont = -999;
  
  double ChisqMinStart = conv::kALargeNum;
  int binStart_X = fStartBin_X;
  int binStart_Y = fStartBin_Y;

  TRandom3 *r = new TRandom3(0); 

  double like_iter = -999.;
  double like_prev = -999.;
  

  // Single starting point

  if(fExtractStartPt=="Single"){
    cout << "Single starting point:" << endl;
    ChisqMinStart = CalcChisqMargPt(skobs_data,fStartBin_X-1,fStartBin_Y-1);
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

      if(!fIsFit2D){
	//bincont = fChisqMargVArray->GetEntry(bs_x-1,fmin_bin_sub_param[1]);
	//if(bincont<=0.) bincont = CalcChisqMargPt(skobs_data,bs_x-1,fmin_bin_sub_param[1]);
	bincont = CalcChisqMargPt(skobs_data,bs_x-1,fmin_bin_sub_param[1]);
      }
      else{
	//bincont = fChisqMargVArray->GetEntry(bs_x-1,bs_y-1);
	//if(bincont<=0.) bincont = CalcChisqMargPt(skobs_data,bs_x-1,bs_y-1);
	bincont = CalcChisqMargPt(skobs_data,bs_x-1,bs_y-1);
      }
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
    
    double step_X = (double)(fLastVal_X - fFirstVal_X) / (double)Nsteps_X;
    double step_Y = (double)(fLastVal_Y - fFirstVal_Y) / (double)Nsteps_Y;
    
    for(int iStartPt_X=0;iStartPt_X<Npt_X;iStartPt_X++){
      for(int iStartPt_Y=0;iStartPt_Y<Npt_Y;iStartPt_Y++){
	
	double currval0 = fFirstVal_X + iStartPt_X*step_X;
	double currval1 = fFirstVal_Y + iStartPt_Y*step_Y;
	//cout << currval0 << " " << currval1 << endl;

	int binX = -999;
	int binY = -999;
	
	binX = hchisq_marg_2D->GetXaxis()->FindBin(currval0+conv::kASmallNum);
	if(fIsFit2D) binY = hchisq_marg_2D->GetYaxis()->FindBin(currval1+conv::kASmallNum);
	//cout << binX << " " << binY << endl;

	float bincont = -999;

	if(!fIsFit2D) bincont = CalcChisqMargPt(skobs_data,binX-1,fmin_bin_sub_param[1]);
	else bincont = CalcChisqMargPt(skobs_data,binX-1,binY-1);

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
    cerr << "The chosen strategy " << fExtractStartPt;
    cerr << " for the chisq minimization starting point  is not available" << endl;
    cerr << endl;
    exit(1);
  }

  float startX = hchisq_marg_2D->GetXaxis()->GetBinLowEdge(fStartBin_X);
  float startY = hchisq_marg_2D->GetYaxis()->GetBinLowEdge(fStartBin_Y);

  cout << "fStartBin_X = " << fStartBin_X << " - fStartBin_Y = " << fStartBin_Y << endl;
  cout << "X = " << startX << " - Y = " << startY << endl;
  cout << "Starting Chisq = " << ChisqMinStart << endl;
  cout << endl;
    

  return;



  // //////////// MCMC extraction (Metropolis-Hasting) //////////
  // int binStartArray_X[NThrowsExtrBin];
  // int binStartArray_Y[NThrowsExtrBin];
  // double ChisqStartArray[NThrowsExtrBin];
  // // Compute chisq and likelihood in the starting point
  // if(!fIsFit2D) bincont = CalcChisqMargPt(skobs_data,binStart_X-1,fmin_bin_sub_param[1]);
  // else bincont = CalcChisqMargPt(skobs_data,binStart_X-1,binStart_Y-1);
  // like_prev = exp(-bincont/2.);
  // binStartArray_X[0] = binStart_X;
  // binStartArray_Y[0] = binStart_Y;
  // ChisqStartArray[0] = bincont;
  // for(int iStartPt=1;iStartPt<NThrowsExtrBin;iStartPt++){    
  //   // Extract iteration likelihood
  //   //double r_X = r->Uniform(UniFirstBin_X,UniLastBin_X+1); // +1 because it's then truncated
  //   //double r_Y = r->Uniform(UniFirstBin_Y,UniLastBin_Y+1); // +1 because it's then truncated
  //   //int bs_x = (int)r_X; 
  //   //int bs_y = (int)r_Y;
  //   double step_X = 10;
  //   double step_Y = 10;
  //   double r_X = r->Uniform(-step_X,step_X+1); // shift of the bin number. +1 for rounding
  //   double r_Y = r->Uniform(-step_Y,step_Y+1); // shift of the bin number. +1 for rounding
  //   int bs_x = binStart_X + (int)r_X; 
  //   int bs_y = binStart_Y + (int)r_Y;
  //   if(bs_x < 1) bs_x = 1;
  //   if(bs_y < 1) bs_y = 1;
  //   if(bs_x > hchisq_marg_2D->GetNbinsX()) bs_x = hchisq_marg_2D->GetNbinsX();
  //   if(bs_y > hchisq_marg_2D->GetNbinsY()) bs_y = hchisq_marg_2D->GetNbinsY();
  //   if(!fIsFit2D){
  //     bincont = fChisqMargVArray->GetEntry(bs_x-1,fmin_bin_sub_param[1]);
  //     if( bincont<=0.) bincont = CalcChisqMargPt(skobs_data,bs_x-1,fmin_bin_sub_param[1]);
  //   }
  //   else{
  //     bincont = fChisqMargVArray->GetEntry(bs_x-1,bs_y-1);
  //     if( bincont<=0.) bincont = CalcChisqMargPt(skobs_data,bs_x-1,bs_y-1);
  //   }
  //   like_iter = exp(-bincont/2.);
  //   // Store all the extracted chisq values
  //   binStartArray_X[iStartPt] = bs_x;
  //   binStartArray_Y[iStartPt] = bs_y;
  //   ChisqStartArray[iStartPt] = bincont;
  //   // Decide next step
  //   double alpha = like_iter / like_prev;
  //   if(alpha>=1){
  //     cout << "alpha > 1" << endl;
  //     cout << "Move to the next point" << endl;
  //     like_prev = like_iter;
  //     binStart_X = bs_x;
  //     binStart_Y = bs_y;
  //     ChisqMinStart = bincont;
  //     continue;
  //   }
  //   double unif = r->Uniform(0,1);
  //   cout << endl;
  //   cout << "like_prev = " << like_prev << endl;
  //   cout << "like_iter = " << like_iter << endl;
  //   cout << "bs_x = " << bs_x << " - bs_y = " << bs_y << endl;
  //   cout << "val_x = " << hchisq_marg_2D->GetXaxis()->GetBinLowEdge(bs_x) << " - ";
  //   cout << "val_y = " << hchisq_marg_2D->GetYaxis()->GetBinLowEdge(bs_y) << endl;
  //   cout << "alpha = " << alpha << " - unif = " << unif << endl;
  //   if(alpha>=unif){ // move to the extracted point
  //     cout << "alpha > unif" << endl;
  //     cout << "Move to the next point" << endl;
  //     like_prev = like_iter;
  //     binStart_X = bs_x;
  //     binStart_Y = bs_y;
  //     ChisqMinStart = bincont;
  //   }
  //   else{ // stay in the extracted point
  //     // Don't do anything
  //     cout << "Keep the previous point" << endl;
  //   }
  // }

}
//________________________________________________________________________________________
void FitMarginalization::SetStartBin(double start_X,double start_Y){
  
  if(!fIsFit2D){
    fStartBin_X = hchisq_marg_1D->GetXaxis()->FindBin(start_X);
    fStartBin_Y = 1;
  }
  else{
    fStartBin_X = hchisq_marg_2D->GetXaxis()->FindBin(start_X);
    fStartBin_Y = hchisq_marg_2D->GetYaxis()->FindBin(start_Y);
  }

  bool GoodStartBin = true;
  
  if(!fIsFit2D){
    if(fStartBin_X < 1 || fStartBin_X > hchisq_marg_1D->GetNbinsX() ){ 
      cerr << endl;
      cerr << "The starting X bin is not allowed" << endl;
      cerr << "  - X = " << start_X << " --> bin = " << fStartBin_X << endl;
      cerr << endl;
      GoodStartBin = false;
    }
  }
  else{
    if(fStartBin_X < 1 || fStartBin_X > hchisq_marg_2D->GetNbinsX() ){ 
      cerr << endl;
      cerr << "The starting X bin is not allowed" << endl;
      cerr << "  - X = " << start_X << " --> bin = " << fStartBin_X << endl;
      cerr << endl;
      GoodStartBin = false;
    }
    if(fStartBin_Y < 1 || fStartBin_Y > hchisq_marg_2D->GetNbinsY() ){ 
      cerr << endl;
      cerr << "The starting Y bin is not allowed" << endl;
      cerr << "  - Y = " << start_Y << " --> bin = " << fStartBin_Y << endl;
      cerr << endl;
      GoodStartBin = false;
    }
  }
  if(!GoodStartBin) exit(1);

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







void FitMarginalization::CalcBestFitErrors(){
  
  // Calculate the best-fit errors as the interval between the points
  // ChisqMarg = 
  if(fIsFit2D){
  }

  else{
  }

  return;
}
//________________________________________________________________________________________
