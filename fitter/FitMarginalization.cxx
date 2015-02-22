//________________________________________________________________________________________
/*
 
 \copyright (c) 2010-2015, VaLOR group.

*/
//________________________________________________________________________________________

//#define __print_array_toy__

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
FitMarginalization::FitMarginalization(conv::PhysParamId_t parfree, int ntoys_marg, int nbins_param,double minpar,double maxpar)
{
  this->Init(parfree,ntoys_marg,nbins_param,minpar,maxpar);
}
//________________________________________________________________________________________
FitMarginalization::~FitMarginalization(void)
{
  this->CleanUp();
}
//________________________________________________________________________________________
void FitMarginalization::Init(conv::PhysParamId_t parfree, int ntoys_marg, int nbins_param,double minpar,double maxpar){
   
  fntoys_marg = ntoys_marg;
  fnbins_param = nbins_param;  

  fmin_bin_sub_param = -999; // scan all the range by default
  fmax_bin_sub_param = -999;

  fSeed = 999;

  //sprintf(fnamefilebinary,"%s","%s_file%d_toy%d_x_oscpar%d.bin");
  //sprintf(fnamefilebinary,"%s","%s_file%d_toy%d_x_oscpar%d%s%s");   // OK
  sprintf(fnamefilebinary,"%s","%s_file%d_toy%d_x_oscpar%d%s%s");  
  
  this->SetFreeParam(parfree,minpar,maxpar);
 
  cout << endl;
  cout << endl;
  cout << "Minimize " << fnamefreeparam << endl;
  cout << "# toys for marginalization: " << fntoys_marg << endl;
  cout << "# bins of free parameter to scan: " << fnbins_param << endl;
  cout << "# free parameter step size: " << fstep_param << endl;
  cout << "Interval: [" << fmin_param << "," << fmax_param << "]" << endl;
  cout << endl;

  fnerecobins_RHCNumu = util::NumErecoBins(kSKsRHCFCFV1RingMuLike),util::ErecoBinEdges(kSKsRHCFCFV1RingMuLike);

  hchisq_marg = new TH1F("hchisq_marg","hchisq_marg",fnbins_param,fmin_param,fmax_param+fstep_param);
  hToyMargCurr = new TH1F("hToyMargCurr","hToyMargCurr",util::NumErecoBins(kSKsRHCFCFV1RingMuLike),util::ErecoBinEdges(kSKsRHCFCFV1RingMuLike));

  fToyMargErecoVArray_RHCNumu = new VALORArrayF(ntoys_marg, nbins_param, util::NumErecoBins(kSKsRHCFCFV1RingMuLike));
  fToyMargErecoVArray_RHCNumu->SetAll(0.);
  fChisqMargVArray            = new VALORArrayF(nbins_param);
  fChisqMargVArray->SetAll(0.);
  flikeVArray                 = new VALORArrayF(nbins_param);
  flikeVArray->SetAll(0.);

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
void FitMarginalization::SetFreeParam(conv::PhysParamId_t parfree,double minpar,double maxpar){
  
  fPhysParamFree = parfree;
  fnamefreeparam=util::AsString(fPhysParamFree);    

  if(minpar==-999 && maxpar==-999){ // if limits are not given use the default ones
 
    // NuBar Osc parameters
    if(fPhysParamFree==kPhypMixing13Bar){ // mixing13Bar
      fmin_param = mixing13Bar_min;
      fmax_param = mixing13Bar_max;
    }    
    else if(fPhysParamFree==kPhypMixing23Bar){ // mixing23Bar
      fmin_param = mixing23Bar_min;
      fmax_param = mixing23Bar_max;
    }    
    else if(fPhysParamFree==kPhypDm23sqBar){ // Dm23sqBar
      fmin_param = dmsq32Bar_min;
      fmax_param = dmsq32Bar_max;
    }    
    else if(fPhysParamFree==kPhypCPVBar){ // deltaCPBar
      fmin_param = deltacpBar_min;
      fmax_param = deltacpBar_max;
    }    
    else if(fPhysParamFree==kPhypMixing12Bar){ // mixing12Bar
      fmin_param = mixing12Bar_min;
      fmax_param = mixing12Bar_max;
    }    
    else if(fPhysParamFree==kPhypDm12sqBar){ // Dm12sqBar
      fmin_param = dmsq12Bar_min;
      fmax_param = dmsq12Bar_max;
    }    
    // Nu Osc parameters
    else if(fPhysParamFree==kPhypMixing13){ // mixing13
      fmin_param = mixing13_min;
      fmax_param = mixing13_max;
    }    
    else if(fPhysParamFree==kPhypMixing23){ // mixing23
      fmin_param = mixing23_min;
      fmax_param = mixing23_max;
    }    
    else if(fPhysParamFree==kPhypDm23sq){ // Dm23sq
      fmin_param = dmsq32_min;
      fmax_param = dmsq32_max;
    }    
    else if(fPhysParamFree==kPhypCPV){ // deltaCP
      fmin_param = deltacp_min;
      fmax_param = deltacp_max;
    }    
    else if(fPhysParamFree==kPhypMixing12){ // mixing12
      fmin_param = mixing12_min;
      fmax_param = mixing12_max;
    }    
    else if(fPhysParamFree==kPhypDm12sq){ // Dm12sq
      fmin_param = dmsq12_min;
      fmax_param = dmsq12_max;
    }    
    else{
      cout << "Free parameter has not well set" << endl;
      exit(1);
    }
  }
  
  else if(minpar > maxpar){
    cerr << endl
	 << "minpar >= maxpar && minpar!=-999 && maxpar!=-999" << endl
	 << " - minpar = " << minpar << endl
	 << " - maxpar = " << maxpar << endl
	 << endl;
    exit(1);
  }
  
  else{
    fmin_param = minpar;
    fmax_param = maxpar;
  }

  fstep_param = (fmax_param-fmin_param)/(double)(fnbins_param-1);

  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetFitParamList(char * fit_param_list){
  fparm  = util::NamedParameterList(fit_param_list);
  fhypothesis = fparm->PhysicsHypothesis();

  if( fparm->IsFixed(fPhysParamFree) ){
    cerr << "The free parameter " << fnamefreeparam 
	 << " is fixed in FitParamList " << fit_param_list
	 << endl << endl;
    exit(1);
  }
  else if( fparm->IsNuisance(fPhysParamFree) ){
    cerr << "The free parameter " << fnamefreeparam 
	 << " is nuisance in FitParamList " << fit_param_list
	 << endl << endl;
    exit(1);
  }

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

  if(min_bin_sub_param==-999 && max_bin_sub_param==-999){ // scan all the range
    min_bin_sub_param = 0;
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
    
    // Initialize the seed to the experiment number. fSeed is 999

    UInt_t myseed = 99;
    Rnd::Instance()->SetSeeds(myseed+iexperiment+10,myseed+iexperiment+20,myseed+iexperiment+30,myseed+iexperiment+40,fSeed+iexperiment+50); // stat, syst, osc, fitter, other

    cout << endl;
    cout << "iexperiment " << iexperiment << endl;
    cout << "# params " << fnbins_param << endl;
    
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
      }
      else { // no toys
	skexp->ResetSystematicVariations();
	cout << "- Generating toy-MC without systematic variation" << endl;
      }


      // Set seeds after systematic paramter extraction

      seedstat = Rnd::Instance()->RndToyMCstat().GetSeed();
      seedsyst = Rnd::Instance()->RndToyMCsyst().GetSeed();
      seedosc = Rnd::Instance()->RndToyMCosc().GetSeed();    
      seedfit = Rnd::Instance()->RndFitter().GetSeed();
      seedoth = Rnd::Instance()->RndOther().GetSeed();
      Rnd::Instance()->SetSeeds(seedstat+1,seedsyst+1,seedosc+1,seedfit+1,seedoth+1); // stat, syst, osc, fitter, other
    
      // cout << endl;
      // cout << "seedstat = " << seedstat+1 << endl;
      // cout << "seedsyst = " << seedsyst+1 << endl;
      // cout << "seedosc = " << seedosc+1 << endl;
      // cout << "seedfit = " << seedfit+1 << endl;
      // cout << "seedoth = " << seedoth+1 << endl;
      // Rnd::Instance()->Print(std::cout);
      
      //


      // Produce toys for each osc param value
      
      for(int ipar=min_bin_sub_param;ipar<=max_bin_sub_param;ipar++){

	double currval = fmin_param + ipar*fstep_param;
	//cout << "currval = " << currval << ", ";
	
	fparm->Update(fPhysParamFree,currval); // NOT THE SAME FITPARAMLIST AS THE ONE USED FOR THE TOYS
	applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), false); // build marginalization distribution w/ syst + osc param but no stat fluctuations
	//cout << *fparm << endl;

	for(int ib = 0; ib < util::NumErecoBins(kSKsRHCFCFV1RingMuLike); ib++){
	  double bincont = skobs->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ib]; 
	  fToyMargErecoVArray_RHCNumu->SetEntry(bincont,iexperiment,ipar,ib); 
	}
	

#ifdef __print_array_toy__
	cout << endl;
	cout << "toy=" << iexperiment << " - ipar=" << ipar << " : ";
	for(int ib = 0; ib < util::NumErecoBins(kSKsRHCFCFV1RingMuLike); ib++){
	  double bincont = skobs->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ib]; 
	  cout << bincont << " "; 
	}
	cout << endl;	
#endif

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
    
    double currval = fmin_param + ipar*fstep_param;
    //cout << "currval = " << currval << endl;

    for(int iexperiment = 0; iexperiment < fntoys_marg; iexperiment++) {
      
      // compute the likelihood for a single marginalization toy
      double chisq = 0.;
      for(int ibin=0;ibin<util::NumErecoBins(kSKsRHCFCFV1RingMuLike);ibin++){
	double nexp = fToyMargErecoVArray_RHCNumu->GetEntry(iexperiment,ipar,ibin);
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
  //sprintf(outname,"%s_file%d_toy%d_x_oscpar%d.bin",tag,idfile,fntoys_marg,fnbins_param);
  sprintf(outname,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param,"",".bin");
  
  std::string soutname(outname);
  fToyMargErecoVArray_RHCNumu->WriteFile(soutname);
  
  cout << " - Toy for marginalization are stored in " << soutname << endl;

  // Write store sub-range of the free parameter

  int min_bin_sub_param = fmin_bin_sub_param;
  int max_bin_sub_param = fmax_bin_sub_param;

  char outname_range[200];
  //sprintf(outname_range,"%s_file%d_toy%d_x_oscpar%d_subrange.bin",tag,idfile,fntoys_marg,fnbins_param);
  sprintf(outname_range,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param,"_subrange",".bin");
  std::string soutname_range(outname_range);
  
  ofstream os(soutname_range.c_str(),ios::binary|ios::out);
  assert(os);
  os.write( (char*)&fmin_bin_sub_param,sizeof(int));
  os.write( (char*)&fmax_bin_sub_param,sizeof(int));
  os.close();

  cout << " - Free param sub-range for marginalization is stored in " << soutname_range << endl;
  cout << endl;

  //sprintf(outname_range,"%s_file%d_toy%d_x_oscpar%d_subrange.txt",tag,idfile,fntoys_marg,fnbins_param);
  sprintf(outname_range,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param,"_subrange",".txt");


  fmin_sub_param = hchisq_marg->GetBinLowEdge(fmin_bin_sub_param+1);    
  fmax_sub_param = hchisq_marg->GetBinLowEdge(fmax_bin_sub_param+1);    

  ofstream ostext(outname_range,ios::out);
  ostext << "Used sub-range of the free parameter: " << fnamefreeparam << endl;
  ostext << "First bin of free param: " << fmin_bin_sub_param << endl;
  ostext << "Last bin of free param: " << fmax_bin_sub_param << endl;  
  ostext << endl;
  ostext << "fnamefreeparam: " << endl
	 << " - " << soutname << endl
	 << " - " << soutname_range << endl;
  ostext << "fnbins_param = " << fnbins_param << endl;
  ostext << "fntoys_marg = " << fntoys_marg << endl;
  ostext << "fnerecobins_RHCNumu = " << fnerecobins_RHCNumu << endl;
  ostext << "fmin_bin_sub_param = " << fmin_bin_sub_param << endl;
  ostext << "fmax_bin_sub_param = " << fmax_bin_sub_param << endl;
  ostext << "fmin_sub_param = " << fmin_sub_param << endl;
  ostext << "fmax_sub_param = " << fmax_sub_param << endl;
  ostext << "fmin_param = " << fmin_param << endl;
  ostext << "fmax_param = " << fmax_param << endl;
  ostext << "fstep_param = " << fstep_param << endl;
  ostext << endl;
  ostext << "fhypothesis = " << util::AsString(fhypothesis) << endl;
  ostext << endl;
  ostext << "FitParamList:" << endl;
  ostext << *fparm << endl;
  
  ostext.close();
}
//________________________________________________________________________________________
void FitMarginalization::ReadToysFromBinary(int idfile, char *tag){ 

  char inname[200];
  //sprintf(inname,"%s_file%d_toy%d_x_oscpar%d.bin",tag,idfile,fntoys_marg,fnbins_param);
  sprintf(inname,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param,"",".bin");

  std::string sinname(inname);

  cout << sinname << endl;

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
//________________________________________________________________________________________
TH1F * FitMarginalization::GetToyMargHisto(SuperkSample_t s,int itoy,int ifreepar){

  if(s==kSKsRHCFCFV1RingMuLike){ 
    for(int ibin=0;ibin<util::NumErecoBins(kSKsRHCFCFV1RingMuLike);ibin++){
      double bincont = fToyMargErecoVArray_RHCNumu->GetEntry(itoy,ifreepar,ibin);
      hToyMargCurr->SetBinContent(ibin+1,bincont);
      //cout << bincont << " ";
    }
    //cout << endl;
    return hToyMargCurr;
  }
  else{
    cout << "You have chosen a wrong Superk Sample" << endl;
    exit(1);
  }
}
//________________________________________________________________________________________
void FitMarginalization::ResetChisqMarg(int min_bin_sub_param,int max_bin_sub_param){
  for(int ipar=fmin_bin_sub_param;ipar<=fmax_bin_sub_param;ipar++){
    double bincont = 0;
    fChisqMargVArray->SetEntry(bincont,ipar);
    hchisq_marg->SetBinContent(ipar+1,bincont);
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
    
    bool doMarginalize = true;
    
    if(p!=fPhysParamFree) 
      fparm->SetDoMarg(fPhysParamFree,true); 
    else{
      fparm->SetDoMarg(fPhysParamFree,false);  // Do not marginalize if free parameter
      continue; 
    }

    if(util::IsSystematic(p)) continue; // not oscillation parameters 
    if(fparm->IsFixed(p))     continue; // do not marginalize if fixed

    // Now extract the oscillation parameter to marginalize 
    double currval = Randomize(p);
    
    fparm->Update(p,currval); // Update current values of parameters 

    cout << util::AsString(p) << " " << currval << endl;
   }

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
