//________________________________________________________________________________________
/*
 
 \copyright (c) 2010-2015, VaLOR group.

*/
//________________________________________________________________________________________

//#define __print_array_toy_maketoys__
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

#include "Conventions/Opt.h"
#include "Utilities/Rnd.h"

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

  sprintf(fnamefilebinary,"%s","%s_file%d_toy%d_x_oscpar%d_oscpar%d_MH%s%s%s");  
  
  MergeInputFiles_Called = true; // false the first time MergeInputFiles() is calles
  fNmergedFiles = 0; // fNmergedFiles++ everytime MergeInputFiles() is called (after reading each file)

  const int nmh=2; // mass hierarchy  


  // Get number of sk samples used in the fit
  // Needed to initialize the Ereco arrays

  this->SetNSKsamples(); // also defines fnerecobins
  for(int is=0;is<fNSKsamples;is++) fToyMargErecoVArray = new VALORArrayF(fntoys_marg, fnbins_sub_param[0], fnbins_sub_param[1], nmh, fNSKsamples, fnerecobins->GetEntry(is)); 
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
}
//________________________________________________________________________________________
void FitMarginalization::CleanUp(void){

  if(fToyMargErecoVArray) delete fToyMargErecoVArray; fToyMargErecoVArray=0;
  if(fToyMargErecoVArray_RHCNumu) delete fToyMargErecoVArray_RHCNumu; fToyMargErecoVArray_RHCNumu=0;  
  if(fChisqMargVArray) delete fChisqMargVArray; fChisqMargVArray=0;
  if(flikeVArray) delete flikeVArray; flikeVArray=0;
  if(hchisq_marg_1D) delete hchisq_marg_1D; hchisq_marg_1D=0;
  if(hchisq_marg_2D) delete hchisq_marg_2D; hchisq_marg_2D=0;
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
  
  if(!fIsFit2D){
    fmin_bin_sub_param[1] = 0.;
    fmax_bin_sub_param[1] = 0.;
  }

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
	else unphysical_tweak = false;
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
	      
	      // Build toy Ereco distributions
	      for(int ib = 0; ib < util::NumErecoBins(s); ib++){
		double bincont = skobs->ErecoSpectrumArray(s)[ib];
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
      
      //cout << ipar0 << "-" << ipar1 << "   ";

      double currval0 = fmin_param[0] + ipar0*fstep_param[0];
      double currval1 = fmin_param[1] + ipar1*fstep_param[1];
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
	  //cout << "toy=" << iexperiment << " - ipar0 = " << ipar0 << " - ipar1 = " << ipar1 << endl;;

	  // compute the likelihood for a single marginalization toy
	  double chisq = 0.;
	  
	  // Loop over the used SK samples
	    int iscount=0;
	    //for(int is = 0; is < kns; is++) {
	    //SuperkSample_t s = kSKSamples[is];
	    //if(s == kSKsUndefined) continue;
	    for(int is=0;is<fNSKsamples;is++){
	      SuperkSample_t s = this->GetSKSamples(is);
	      
	      for(int ibin=0;ibin<util::NumErecoBins(s);ibin++){
		double nexp = fToyMargErecoVArray->GetEntry(iexperiment,ipar0-fmin_bin_sub_param[0],ipar1-fmin_bin_sub_param[1],iMH,iscount,ibin);
		double nobs = skobs_data->ErecoSpectrumArray(s)[ibin];
		chisq += valor::util::ChisqPoisson(nobs,nexp); // penalty term is a simple constant --> not needed
	      }	  
	      iscount++;
	    } // loop SK samples
	      
#ifdef __print_array_toy_fillchisq__
	  cout << endl;
	  cout << "toy=" << iexperiment << " - ipar0=" << ipar0 << " - ipar1=" << ipar1 << " : ";
	  
	  iscount=0;
	  //for(int is = 0; is < kns; is++) {
	  //SuperkSample_t s = kSKSamples[is];
	  //if(s == kSKsUndefined) continue;
	  for(int is=0;is<fNSKsamples;is++){
	    SuperkSample_t s = this->GetSKSamples(is);
	    
	    cout << util::AsStringCompact(s) << endl;

	    for(int ib = 0; ib < util::NumErecoBins(s); ib++){
	      double bincont = fToyMargErecoVArray->GetEntry(iexperiment,ipar0-fmin_bin_sub_param[0],ipar1-fmin_bin_sub_param[1],iMH,iscount,ib);
	      cout << bincont << " "; 
	    }
	    cout << endl;
	    cout << endl;
	    
	    iscount++;
	  } // loop SK samples

#endif
	  
	  double like = exp(-chisq/2.) / fntoys_marg;

	  if(this->doMargMassHierarchy) like /= 2.;  // if MH is marginalized, divide by 2
	  
	  flikeVArray->AddToEntry(like,ipar0,ipar1); // this is correct if binary files differ only for free param sub-ranges
	  
	  NtoysTot++;
	} // iexperiment
	
      } // iMH (mass hierarchy)
      
      NparTot1++;
    } // ipar1 (free parameter)
    NparTot0++;
  } // ipar0 (free parameter)

  cout << endl;
  cout << "# of scanned free parameters 0: " << NparTot0 << endl; 
  cout << "# of scanned free parameters 1: " << NparTot1 << endl; 
  if(doMargMassHierarchy) cout << "# of toys used: " << (double)NtoysTot / 2. / (double)NparTot0 / (double)NparTot1 << endl;
  else cout << "# of toys used: " << (double)NtoysTot / (double)NparTot1 /  (double)NparTot0 << endl;
  cout << endl;

  // Now the marginalized likelihood is computed using toys from all the binary files
  // Compute the chisq for each free osc parameter value
  
  for(int ipar0=fmin_bin_sub_param[0];ipar0<=fmax_bin_sub_param[0];ipar0++){
    for(int ipar1=fmin_bin_sub_param[1];ipar1<=fmax_bin_sub_param[1];ipar1++){
      double bincont = -2*log( flikeVArray->GetEntry(ipar0,ipar1) );
      fChisqMargVArray->SetEntry(bincont,ipar0,ipar1);
      if(!fIsFit2D) hchisq_marg_1D->SetBinContent(ipar0+1,bincont);
      else          hchisq_marg_2D->SetBinContent(ipar0+1,ipar1+1,bincont);
    }
  }

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
void FitMarginalization::WriteToys2Binary(int idfile,const char *tag){

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
void FitMarginalization::Minimize(){
  
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
void FitMarginalization::Minimize(TH1F *hInChisqMarg){
  
  int minbin = -999;
  
  fchisq_bf = hInChisqMarg->GetMinimum();
  minbin = hInChisqMarg->GetMinimumBin();
  fparam_bf[0] = hInChisqMarg->GetBinLowEdge(minbin);
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::Minimize(TH2F *hInChisqMarg){
  
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
string FitMarginalization::GetPrior(PhysParamId_t physpar){
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
bool FitMarginalization::RethrowOscUnphysical(conv::PhysParamId_t p,double oscval){
  
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
double FitMarginalization::Randomize(PhysParamId_t physpar){

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
     
     this->fsksamples->SetEntry(is,iscount);
     this->fnerecobins->SetEntry(util::NumErecoBins(s),iscount); 
     
     //cout << "is " << is << " - iscount " << iscount << endl;
     //cout << fnerecobins->GetEntry(iscount) << endl;
     
     iscount++;
  }
  return;
}
//________________________________________________________________________________________
conv::SuperkSample_t FitMarginalization::GetSKSamples(int isksample){
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
