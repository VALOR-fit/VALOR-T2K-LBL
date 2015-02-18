//________________________________________________________________________________________
/*
 
 \copyright (c) 2010-2015, VaLOR group.

*/
//________________________________________________________________________________________


#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <TSystem.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TString.h>
#include <TLegend.h>
#include <TSpline.h>
#include <TRandom3.h>

#include "FitMarginalization.h"
#include "SuperkMC.h"
#include "SuperkDetectorFSI.h"
#include "NuFlux.h"
#include "NuXSec.h"
#include "BANFF.h"
#include "Utilities.h"
#include "Conventions/Opt.h"
#include "Utilities/Rnd.h"

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
FitMarginalization::FitMarginalization(conv::PhysParamId_t parfree, int ntoys_marg, int nbins_param)
{
  this->Init(parfree,ntoys_marg,nbins_param);
}
//________________________________________________________________________________________
FitMarginalization::~FitMarginalization(void)
{
  this->CleanUp();
}
//________________________________________________________________________________________
void FitMarginalization::Init(conv::PhysParamId_t parfree, int ntoys_marg, int nbins_param){
   
  fntoys_marg = ntoys_marg;
  fnbins_param = nbins_param;  

  fmin_bin_sub_param = 0; // scan all the range by default
  fmax_bin_sub_param = 0;

  this->SetFreeParam(parfree);
 
  cout << endl;
  cout << endl;
  cout << "Minimize " << fnamefreeparam << endl;
  cout << "# toys for marginalization: " << fntoys_marg << endl;
  cout << "# bins of free parameter to scan: " << fnbins_param << endl;
  cout << "# free parameter step size: " << fstep_param << endl;
  cout << endl;

  hchisq_marg = new TH1F("hchisq_marg","hchisq_marg",fnbins_param,fmin_param,fmax_param+fstep_param);

  fToyMargErecoVArray_RHCNumu = new VALORArrayF(ntoys_marg, nbins_param, util::NumErecoBins(kSKsRHCFCFV1RingMuLike));
  fChisqMargVArray            = new VALORArrayF(nbins_param);
  flikeVArray                 = new VALORArrayF(nbins_param);
  
  fToyMargErecoVArray_RHCNumu->SetDimensionName(0, "# toys marginalization");
  fToyMargErecoVArray_RHCNumu->SetDimensionName(1, "# bins free parameter");
  fToyMargErecoVArray_RHCNumu->SetDimensionName(2, "# bins ereco");
}
//________________________________________________________________________________________
void FitMarginalization::CleanUp(void){

  if(fToyMargErecoVArray_RHCNumu) delete fToyMargErecoVArray_RHCNumu; fToyMargErecoVArray_RHCNumu=0;
  if(fChisqMargVArray) delete fChisqMargVArray; fChisqMargVArray=0;
  if(flikeVArray) delete flikeVArray; flikeVArray=0;
  if(hchisq_marg) delete hchisq_marg; hchisq_marg=0;
}
//________________________________________________________________________________________
void FitMarginalization::SetFreeParam(conv::PhysParamId_t parfree){
  
  fPhysParamFree = parfree;
  fnamefreeparam=util::AsString(fPhysParamFree);    

  if(fPhysParamFree==kPhypMixing23Bar){ // mixing23Bar
    fmin_param = mixing23Bar_min;
    fmax_param = mixing23Bar_max;
    fstep_param = (fmax_param-fmin_param)/(fnbins_param-1);
  }    
  else{
    cout << "Free parameter has not well set" << endl;
    exit(1);
  }
  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetFitParamList(char * fit_param_list){
  fparm  = util::NamedParameterList(fit_param_list);
  fhypothesis = fparm->PhysicsHypothesis();
}
//________________________________________________________________________________________
void FitMarginalization::PrintPhysHypothesis(){
  cout << endl;
  cout << "- The following physics hypothesis is used in the fit: " << util::AsString(fhypothesis) << endl;
  cout << endl;
}
//________________________________________________________________________________________
void FitMarginalization::MakeToysMarginalization(int tweak_systematics,SuperkExpectation *skexp,int min_bin_sub_param,int max_bin_sub_param){

  cout << "////////////////////////////////////////////////////////" << endl;
  cout << endl;
  cout << "Make toys for marginalization and store in a binary file" << endl;
  cout << endl;

  if(min_bin_sub_param==0 && max_bin_sub_param==0){ // scan all the range
    max_bin_sub_param = fnbins_param-1;
  }
  
  // Get free param sub-range
  fmin_bin_sub_param = min_bin_sub_param;
  fmax_bin_sub_param = max_bin_sub_param;

  cout << endl;
  cout << "Scan the parameter range:" << endl
       << "bin: ["  
       << fmin_bin_sub_param
       << ","
       << fmax_bin_sub_param
       << "]"
       << endl
       << fnamefreeparam << ": " 
       << "[" 
       << hchisq_marg->GetBinLowEdge(fmin_bin_sub_param+1)
       << ","
       << hchisq_marg->GetBinLowEdge(fmax_bin_sub_param+1)
       << "]"
       << endl;

  // Initialize the free param sub-range to scan
  fmin_bin_sub_param = min_bin_sub_param;
  fmax_bin_sub_param = max_bin_sub_param;

  // It can extract toys for sub-ranges of the free parameter

  for(int iexperiment = 0; iexperiment < fntoys_marg; iexperiment++) {
    
    cout << endl;
    cout << "iexperiment " << iexperiment << endl;
    cout << "# params " << fnbins_param << endl;
    
    fparm->ResetCurrToInitValues();
    
    SuperkObservation * skobs = new SuperkObservation;
    
    bool applied_negevents_condition = true;
    while(applied_negevents_condition) {
      if(tweak_systematics) {
	
	// tweak all sources of systematic error
	if(tweak_systematics == 1) {
	  cout << "- Randomizing all systematic params in pdf..." << endl;
	  skexp->TweakAllSystematics(); 
	}
	else if(tweak_systematics == 2) {
	  cout << "- Randomizing BANFF systematic params in pdf..." << endl;
	  skexp->TweakAllBANFFSystematics();
	}
	else if(tweak_systematics == 3) {
	  cout << "- Randomizing uncorrelated-xsec systematic params in pdf..." << endl;
	  skexp->TweakAllUncorrelatedXSecSystematics();
	}
	else if(tweak_systematics == 4) {
	  cout << "- Randomizing SKDetFSI systematic params in pdf..." << endl;
	  skexp->TweakAllSKDetectorFSISystematics();
	}
	//skexp->PrintSystematicVariations();
	//cout << "- Reseting all systematic params in pdf..." << endl;
	//skexp->ResetSystematicVariations();
      }
      else { // no toys
	skexp->ResetSystematicVariations();
	cout << "- Generating toy-MC without systematic variation" << endl;
      }

      // Produce toys for each osc param value
      
      //for(int ipar=0; ipar<fnbins_param; ipar++){ 
      for(int ipar=min_bin_sub_param;ipar<=max_bin_sub_param;ipar++){

	double currval = fmin_param + ipar*fstep_param;
	//cout << "currval = " << currval << endl;
	
	fparm->Update(fPhysParamFree,currval);
	
	applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), false); // build marginalization distribution w/ syst + osc param but no stat fluctuations
	
	for(int ib = 0; ib < util::NumErecoBins(kSKsRHCFCFV1RingMuLike); ib++){
	  double bincont = skobs->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ib]; 
	  fToyMargErecoVArray_RHCNumu->SetEntry(bincont,iexperiment,ipar,ib); 
	}
	
      } // fnbins_param
    } // applied_negevents_condition
  } // iexperiment
  
  cout << endl;
  cout << "Toys for marginalization are extracted" << endl;
  cout << endl;

  return;
}
//________________________________________________________________________________________
void FitMarginalization::FillChisqMarg(SuperkObservation *skobs_data){
  
  int NtoysTot=0;
  int NparTot=0;
 
  cout << endl;
  cout << "Fill marginalized Chisq in the parameter range:" << endl
       << "bin: ["  
       << fmin_bin_sub_param
       << ","
       << fmax_bin_sub_param
       << "]"
       << endl
       << fnamefreeparam << ": " 
       << "[" 
       << hchisq_marg->GetBinLowEdge(fmin_bin_sub_param+1)
       << ","
       << hchisq_marg->GetBinLowEdge(fmax_bin_sub_param+1)
       << "]"
       << endl;

  // It can marginalized sub-ranges of the free parameter
  // fmin_bin_sub_param and fmax_bin_sub_param read from the binary file

  for(int ipar=fmin_bin_sub_param;ipar<=fmax_bin_sub_param;ipar++){
    for(int iexperiment = 0; iexperiment < fntoys_marg; iexperiment++) {
      
      // compute the likelihood for a single marginalization toy
      double chisq = 0.;
      for(int ibin=0;ibin<util::NumErecoBins(kSKsRHCFCFV1RingMuLike);ibin++){
	double nexp = fToyMargErecoVArray_RHCNumu->GetEntry(iexperiment,ipar,ibin);;
	double nobs = skobs_data->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ibin];
	chisq += valor::util::ChisqPoisson(nobs,nexp); // penalty term is a simple constant --> not needed
      }

      double like = exp(-chisq/2.) / fntoys_marg;
      
      flikeVArray->AddToEntry(like,ipar); // this is correct if binary files differ only for free param sub-ranges
      
      NtoysTot++;

    } // iexperiment

    NparTot++;

  } // ipar
  
  cout << endl;
  cout << "# of scanned free parameters: " << NparTot << endl; 
  cout << "# of toys used: " << (double)NtoysTot / (double)NparTot << endl;
  cout << endl;

  // Now the marginalized likelihood is computed using toys from all the binary files
  // Compute the chisq for each free osc parameter value
  
  for(int ipar=fmin_bin_sub_param;ipar<=fmax_bin_sub_param;ipar++){
    double bincont = -2*log( flikeVArray->GetEntry(ipar) );
    fChisqMargVArray->SetEntry(bincont,ipar);
    hchisq_marg->SetBinContent(ipar+1,bincont);
  }

  return;
}
//________________________________________________________________________________________
void FitMarginalization::WriteToys2Binary(int idfile,char *tag){

  cout << endl;
  cout << "Write marginalization toys in a binary file:" << endl;
  
  char outname[200];
  sprintf(outname,"%s_file%d_toy%d_x_oscpar%d.bin",tag,idfile,fntoys_marg,fnbins_param);
  
  std::string soutname(outname);
  fToyMargErecoVArray_RHCNumu->WriteFile(soutname);
  
  cout << " - Toy for marginalization are stored in " << soutname << endl;

  // Write store sub-range of the free parameter

  int min_bin_sub_param = fmin_bin_sub_param;
  int max_bin_sub_param = fmax_bin_sub_param;

  char outname_range[200];
  sprintf(outname_range,"%s_file%d_toy%d_x_oscpar%d_subrange.bin",tag,idfile,fntoys_marg,fnbins_param);
  std::string soutname_range(outname_range);
  
  ofstream os(soutname_range.c_str(),ios::binary|ios::out);
  assert(os);
  os.write( (char*)&fmin_bin_sub_param,sizeof(int));
  os.write( (char*)&fmax_bin_sub_param,sizeof(int));
  os.close();

  cout << " - Free param sub-range for marginalization is stored in " << soutname_range << endl;
  cout << endl;

  cout << "fmin_bin_sub_param = " << fmin_bin_sub_param << endl;
  cout << "fmax_bin_sub_param = " << fmax_bin_sub_param << endl;

}
//________________________________________________________________________________________
void FitMarginalization::ReadToysFromBinary(int idfile, char *tag){ 

  char inname[200];
  sprintf(inname,"%s_file%d_toy%d_x_oscpar%d.bin",tag,idfile,fntoys_marg,fnbins_param);
  std::string sinname(inname);

  fToyMargErecoVArray_RHCNumu->ReadFile(sinname);

  cout << endl;
  cout << "Read toys for marginalization from binary file:" << sinname << endl;  
  
  // Write store sub-range of the free parameter

  int min_bin_sub_param=0;
  int max_bin_sub_param=0;

  char inname_range[200];
  sprintf(inname_range,"%s_file%d_toy%d_x_oscpar%d_subrange.bin",tag,idfile,fntoys_marg,fnbins_param);
  std::string sinname_range(inname_range);

  cout << "Read free param sub-range from binary file:" << sinname_range << endl;

  ifstream is(sinname_range.c_str(),ios::binary|ios::in);
  assert(is);
  is.read( (char*)&min_bin_sub_param,sizeof(int));
  is.read( (char*)&max_bin_sub_param,sizeof(int));
  is.close();

  fmin_bin_sub_param = min_bin_sub_param;
  fmax_bin_sub_param = max_bin_sub_param;

  fmin_sub_param = hchisq_marg->GetBinLowEdge(fmin_bin_sub_param+1);    
  fmax_sub_param = hchisq_marg->GetBinLowEdge(fmax_bin_sub_param+1);
    
  if(fmin_sub_param<fmin_param || fmax_sub_param>fmax_param){
    cerr << endl 
	 << "The free parameter sub-range ["
	 << fmin_sub_param << "," << fmax_sub_param << "] " 
	 << "provided by the files"
	 << inname << endl
	 << "and" 
	 << inname_range << endl
	 << "exceeds the total range set for" 
	 << fnamefreeparam
	 << " --> [" << fmin_param << "," << fmax_param << "]"
	 << endl;
    exit(1);
  }
  
  cout << "File sub-range: ["
       << fmin_sub_param << "," << fmax_sub_param 
       << "]" << endl;
  cout << "File bin sub-range: ["
       << fmin_bin_sub_param << "," << fmax_bin_sub_param 
       << "]" << endl << endl;
 
  return;
}
//________________________________________________________________________________________
void FitMarginalization::Minimize(){
  
  int minbin = -999;

  fchisq_bf = hchisq_marg->GetMinimum();
  minbin = hchisq_marg->GetMinimumBin();
  fparam_bf = hchisq_marg->GetBinLowEdge(minbin);;

  return;
}
//________________________________________________________________________________________
void FitMarginalization::Minimize(TH1F *hInChisqMarg){
  
  int minbin = -999;
  
  fchisq_bf = hInChisqMarg->GetMinimum();
  minbin = hInChisqMarg->GetMinimumBin();
  fparam_bf = hInChisqMarg->GetBinLowEdge(minbin);
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::PrintBestFit(void){

  cout << endl;
  cout << "Best fit ChisqMarg = " << this->GetChisqBestFit() << endl;
  cout << "Best fit " << fnamefreeparam << " = " << this->GetParamBestFit() << endl;
  cout << endl;

  return;
}
//________________________________________________________________________________________
VALORArrayF * FitMarginalization::GetToyMargErecoVArray(SuperkSample_t s){

  if(s==kSKsRHCFCFV1RingMuLike) return fToyMargErecoVArray_RHCNumu;
  else{
    cout << "You have chosen a wrong Superk Sample" << endl;
    exit(1);
  }
}







// // TO COMPLETE

// //________________________________________________________________________________________
// void FitMarginalization::SetMargOscParam(PhysParamId_t physpar,string randmod,TH1F *hprior){
//   if(physpar==kPhypMixing23Bar){
//   }
// }
// //________________________________________________________________________________________
// void FitMarginalization::RandomizeGaus(PhysParamId_t physpar,double mean,double sigma){
//   double dm    = Rnd::Instance()->RndToyMCosc().Gaus(mean,sigma);
//   if(dm >= 0.)
//     return dm;
//   //or try again
//   return RandomizeGaus();
// }
// //________________________________________________________________________________________
// void FitMarginalization::RandomizeFlat(PhysParamId_t physpar,double min,double max){
//   return Rnd::Instance()->RndToyMCosc().Uniform(min,max);
// }
// //________________________________________________________________________________________
// void FitMarginalization::RandomizePrior(PhysParamId_t physpar,TH1F *hprior){
//   return hprior->GetRandom();
// }
// //________________________________________________________________________________________
