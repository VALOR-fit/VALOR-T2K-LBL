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
   
  doMargMassHierarchy = false;
  fRethrowSystUnphysical = true; // rethrow toys if systematics are unphysical

  fntoys_marg = ntoys_marg;
  fnbins_param_1 = nbins_param;  

  fmin_bin_sub_param_1 = -999; // scan all the range by default
  fmax_bin_sub_param_1 = -999;

  fSeed = 999;

  sprintf(fnamefilebinary,"%s","%s_file%d_toy%d_x_oscpar%d_MH%s%s%s");  
  
  MergeInputFiles_Called = true; // false the first time MergeInputFiles() is calles
  fNmergedFiles = 0; // ++ everytime MergeInputFiles() is called (after reading each file)

  this->SetFreeParam(parfree,minpar,maxpar);
 
  cout << endl;
  cout << endl;
  cout << "Minimize " << fnamefreeparam_1 << endl;
  cout << "# toys for marginalization: " << fntoys_marg << endl;
  cout << "# bins of free parameter to scan: " << fnbins_param_1 << endl;
  cout << "# free parameter step size: " << fstep_param_1 << endl;
  cout << "Interval: [" << fmin_param_1 << "," << fmax_param_1 << "]" << endl;


  fnerecobins_RHCNumu = util::NumErecoBins(kSKsRHCFCFV1RingMuLike),util::ErecoBinEdges(kSKsRHCFCFV1RingMuLike);

  hchisq_marg_1D = new TH1F("hchisq_marg_1D","hchisq_marg_1D",fnbins_param_1,fmin_param_1,fmax_param_1+fstep_param_1);
  hToyMargCurr_1D = new TH1F("hToyMargCurr_1D","hToyMargCurr_1D",util::NumErecoBins(kSKsRHCFCFV1RingMuLike),util::ErecoBinEdges(kSKsRHCFCFV1RingMuLike));

  const int nmh=2; // mass hierarchy
  fToyMargErecoVArray_1D_RHCNumu = new VALORArrayF(ntoys_marg, nbins_param, nmh, util::NumErecoBins(kSKsRHCFCFV1RingMuLike));
  //fToyMargErecoVArray_1D_RHCNumu->SetAll(0.);
  fToyMargErecoVArray_1D_RHCNumu->SetAll(-1.);   // needed for MergeInputFiles()

  fChisqMargVArray_1D            = new VALORArrayF(nbins_param);
  fChisqMargVArray_1D->SetAll(0.);
  flikeVArray_1D                 = new VALORArrayF(nbins_param);
  flikeVArray_1D->SetAll(0.);

  fToyMargErecoVArray_1D_RHCNumu->SetDimensionName(0, "# toys marginalization");
  fToyMargErecoVArray_1D_RHCNumu->SetDimensionName(1, "# bins free parameter");
  fToyMargErecoVArray_1D_RHCNumu->SetDimensionName(2, "# mass hierarchies");
  fToyMargErecoVArray_1D_RHCNumu->SetDimensionName(3, "# bins ereco");
}
//________________________________________________________________________________________
void FitMarginalization::CleanUp(void){

  if(fToyMargErecoVArray_1D_RHCNumu) delete fToyMargErecoVArray_1D_RHCNumu; fToyMargErecoVArray_1D_RHCNumu=0;
  if(fChisqMargVArray_1D) delete fChisqMargVArray_1D; fChisqMargVArray_1D=0;
  if(flikeVArray_1D) delete flikeVArray_1D; flikeVArray_1D=0;
  if(hchisq_marg_1D) delete hchisq_marg_1D; hchisq_marg_1D=0;
}
//________________________________________________________________________________________
void FitMarginalization::SetFreeParam(conv::PhysParamId_t parfree,double minpar,double maxpar){
  
  fPhysParamFree_1 = parfree;
  fnamefreeparam_1=util::AsString(fPhysParamFree_1);    

  if(minpar==-999 && maxpar==-999){ // if limits are not given use the default ones
 
    // NuBar Osc parameters
    if(fPhysParamFree_1==kPhypMixing13Bar){ // mixing13Bar
      fmin_param_1 = mixing13Bar_min;
      fmax_param_1 = mixing13Bar_max;
    }    
    else if(fPhysParamFree_1==kPhypMixing23Bar){ // mixing23Bar
      fmin_param_1 = mixing23Bar_min;
      fmax_param_1 = mixing23Bar_max;
    }    
    else if(fPhysParamFree_1==kPhypDm23sqBar){ // Dm23sqBar
      fmin_param_1 = dmsq32Bar_min;
      fmax_param_1 = dmsq32Bar_max;
    }    
    else if(fPhysParamFree_1==kPhypCPVBar){ // deltaCPBar
      fmin_param_1 = deltacpBar_min;
      fmax_param_1 = deltacpBar_max;
    }    
    else if(fPhysParamFree_1==kPhypMixing12Bar){ // mixing12Bar
      fmin_param_1 = mixing12Bar_min;
      fmax_param_1 = mixing12Bar_max;
    }    
    else if(fPhysParamFree_1==kPhypDm12sqBar){ // Dm12sqBar
      fmin_param_1 = dmsq12Bar_min;
      fmax_param_1 = dmsq12Bar_max;
    }    
    // Nu Osc parameters
    else if(fPhysParamFree_1==kPhypMixing13){ // mixing13
      fmin_param_1 = mixing13_min;
      fmax_param_1 = mixing13_max;
    }    
    else if(fPhysParamFree_1==kPhypMixing23){ // mixing23
      fmin_param_1 = mixing23_min;
      fmax_param_1 = mixing23_max;
    }    
    else if(fPhysParamFree_1==kPhypDm23sq){ // Dm23sq
      fmin_param_1 = dmsq32_min;
      fmax_param_1 = dmsq32_max;
    }    
    else if(fPhysParamFree_1==kPhypCPV){ // deltaCP
      fmin_param_1 = deltacp_min;
      fmax_param_1 = deltacp_max;
    }    
    else if(fPhysParamFree_1==kPhypMixing12){ // mixing12
      fmin_param_1 = mixing12_min;
      fmax_param_1 = mixing12_max;
    }    
    else if(fPhysParamFree_1==kPhypDm12sq){ // Dm12sq
      fmin_param_1 = dmsq12_min;
      fmax_param_1 = dmsq12_max;
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
    fmin_param_1 = minpar;
    fmax_param_1 = maxpar;
  }

  fstep_param_1 = (fmax_param_1-fmin_param_1)/(double)(fnbins_param_1 - 1);

  return;
}
//________________________________________________________________________________________
void FitMarginalization::SetFitParamList(char * fit_param_list){
  fparm  = util::NamedParameterList(fit_param_list);
  fhypothesis = fparm->PhysicsHypothesis();

  if( fparm->IsFixed(fPhysParamFree_1) ){
    cerr << "The free parameter " << fnamefreeparam_1 
	 << " is fixed in FitParamList " << fit_param_list
	 << endl << endl;
    exit(1);
  }
  else if( fparm->IsNuisance(fPhysParamFree_1) ){
    cerr << "The free parameter " << fnamefreeparam_1 
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
    max_bin_sub_param = fnbins_param_1 - 1;
  }
  
  // Get free param sub-range
  fmin_bin_sub_param_1 = min_bin_sub_param;
  fmax_bin_sub_param_1 = max_bin_sub_param;

  cout << endl;
  cout << "Scan the parameter range:" << endl
       << "bin: ["  
       << fmin_bin_sub_param_1
       << ","
       << fmax_bin_sub_param_1
       << "]"
       << endl
       << fnamefreeparam_1 << ": " 
       << "[" 
       << hchisq_marg_1D->GetBinLowEdge(fmin_bin_sub_param_1 + 1)
       << ","
       << hchisq_marg_1D->GetBinLowEdge(fmax_bin_sub_param_1 + 1)
       << "]"
       << endl;

  // Initialize the free param sub-range to scan
  fmin_bin_sub_param_1 = min_bin_sub_param;
  fmax_bin_sub_param_1 = max_bin_sub_param;

  // It can extract toys for sub-ranges of the free parameter

  for(int iexperiment = 0; iexperiment < fntoys_marg; iexperiment++) {
    
    // Initialize the seed to the experiment number. fSeed is 999

    UInt_t myseed = 99;
    Rnd::Instance()->SetSeeds(myseed+iexperiment+10,myseed+iexperiment+20,myseed+iexperiment+30,myseed+iexperiment+40,fSeed+iexperiment+50); // stat, syst, osc, fitter, other

    cout << endl;
    cout << "iexperiment " << iexperiment << endl;
    cout << "# params " << fnbins_param_1 << endl;
    
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

#ifdef __print_systematic_variations__
	skexp->PrintSystematicVariations();
#endif

	// Re-throw the systematics if at least one in the unphysical region
	PhysParamList true_nuisance_params = skexp->Tweaks();
	if (this->fRethrowSystUnphysical) {
	  unphysical_tweak = !(valor::T2KPhysicalParameterLimits::Instance()->AllWithinLimits(true_nuisance_params));
	  if (unphysical_tweak) {
	    cout<<"Parameters are unphysical - rethrowing"<<std::endl;
	    continue; // repeat the loop
	  }
	}
	else {
	  unphysical_tweak = false;
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

	double currval = fmin_param_1 + ipar*fstep_param_1;
	//cout << "currval = " << currval << ", ";
	

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
	  fparm->Update(fPhysParamFree_1,currval); 
	  applied_negevents_condition = skobs->CreateMockData(*skexp, PhysHypothesis(*fparm), false); // build marginalization distribution w/ syst + osc param but no stat fluctuations


	  // Build toy Ereco distributions

	  for(int ib = 0; ib < util::NumErecoBins(kSKsRHCFCFV1RingMuLike); ib++){
	    double bincont = skobs->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ib]; 
	    fToyMargErecoVArray_1D_RHCNumu->SetEntry(bincont,iexperiment,ipar,iMH,ib); 
	  }
	 
	  

#ifdef __print_array_toy_maketoys__
	cout << endl;
	cout << "toy=" << iexperiment << " - ipar=" << ipar << " : ";
	for(int ib = 0; ib < util::NumErecoBins(kSKsRHCFCFV1RingMuLike); ib++){
	  double bincont = skobs->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ib]; 
	  cout << bincont << " "; 
	}
	cout << endl;	
#endif

	} // loop Mass Hierarchy
      } // loop ipar (free parameter)
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
       << fmin_bin_sub_param_1
       << ","
       << fmax_bin_sub_param_1
       << "]"
       << endl
       << fnamefreeparam_1 << ": " 
       << "[" 
       << hchisq_marg_1D->GetBinLowEdge(fmin_bin_sub_param_1 + 1)
       << ","
       << hchisq_marg_1D->GetBinLowEdge(fmax_bin_sub_param_1 + 1)
       << "]"
       << endl;

  // It can marginalize sub-ranges of the free parameter
  // fmin_bin_sub_param_1 and fmax_bin_sub_param_1 read from the binary file

  for(int ipar=fmin_bin_sub_param_1;ipar<=fmax_bin_sub_param_1;ipar++){
    
    double currval = fmin_param_1 + ipar*fstep_param_1;
    //cout << "currval = " << currval << endl;
    
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
	
	// compute the likelihood for a single marginalization toy
	double chisq = 0.;
	for(int ibin=0;ibin<util::NumErecoBins(kSKsRHCFCFV1RingMuLike);ibin++){
	  double nexp = fToyMargErecoVArray_1D_RHCNumu->GetEntry(iexperiment,ipar,iMH,ibin);
	  double nobs = skobs_data->ErecoSpectrumArray(kSKsRHCFCFV1RingMuLike)[ibin];
	  chisq += valor::util::ChisqPoisson(nobs,nexp); // penalty term is a simple constant --> not needed
	}

      
#ifdef __print_array_toy_fillchisq__
	cout << endl;
	cout << "toy=" << iexperiment << " - ipar=" << ipar << " : ";
	for(int ib = 0; ib < util::NumErecoBins(kSKsRHCFCFV1RingMuLike); ib++){
	  double bincont = fToyMargErecoVArray_1D_RHCNumu->GetEntry(iexperiment,ipar,iMH,ib);
	  cout << bincont << " "; 
	}
	cout << endl;	
#endif

	
	double like = exp(-chisq/2.) / fntoys_marg;
	
	if(this->doMargMassHierarchy) like /= 2.;  // if MH is marginalized, divide by 2
	
	flikeVArray_1D->AddToEntry(like,ipar); // this is correct if binary files differ only for free param sub-ranges
	
	NtoysTot++;
	
      } // iexperiment
      
    } // iMH (mass hierarchy)
    
    NparTot++;

  } // ipar (free parameter)
  
  cout << endl;
  cout << "# of scanned free parameters: " << NparTot << endl; 
  if(doMargMassHierarchy) cout << "# of toys used: " << (double)NtoysTot / (double)NparTot / 2. << endl;
  else cout << "# of toys used: " << (double)NtoysTot / (double)NparTot << endl;
  cout << endl;

  // Now the marginalized likelihood is computed using toys from all the binary files
  // Compute the chisq for each free osc parameter value
  
  for(int ipar=fmin_bin_sub_param_1;ipar<=fmax_bin_sub_param_1;ipar++){
    double bincont = -2*log( flikeVArray_1D->GetEntry(ipar) );
    fChisqMargVArray_1D->SetEntry(bincont,ipar);
    hchisq_marg_1D->SetBinContent(ipar+1,bincont);
  }

#ifdef __print_chisqmarg__
  cout << endl;
  cout << "Marginalized Chisq distribution:" << endl; 
  cout << "fmin_bin_sub_param_1 = " << fmin_bin_sub_param_1 << endl;
  cout << "fmax_bin_sub_param_1 = " << fmax_bin_sub_param_1 << endl;
  for(int ipar=fmin_bin_sub_param_1;ipar<=fmax_bin_sub_param_1;ipar++){
    cout << fChisqMargVArray_1D->GetEntry(ipar) << ", ";
  }
  cout << endl;
#endif

  return;
}
//________________________________________________________________________________________
void FitMarginalization::WriteToys2Binary(int idfile,char *tag){

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
  sprintf(outname,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param_1,namemh,"",".bin");

  std::string soutname(outname);
  fToyMargErecoVArray_1D_RHCNumu->WriteFile(soutname);
  
  cout << " - Toy for marginalization are stored in " << soutname << endl;

  // Write store sub-range of the free parameter

  int min_bin_sub_param = fmin_bin_sub_param_1;
  int max_bin_sub_param = fmax_bin_sub_param_1;

  char outname_range[200];
  sprintf(outname_range,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param_1,namemh,"_subrange",".bin");
  
  std::string soutname_range(outname_range);
  
  ofstream os(soutname_range.c_str(),ios::binary|ios::out);
  assert(os);
  os.write( (char*)&fmin_bin_sub_param_1,sizeof(int));
  os.write( (char*)&fmax_bin_sub_param_1,sizeof(int));
  os.close();

  cout << " - Free param sub-range for marginalization is stored in " << soutname_range << endl;
  cout << endl;

  //

  sprintf(outname_range,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param_1,namemh,"_subrange",".txt");

  fmin_sub_param_1 = hchisq_marg_1D->GetBinLowEdge(fmin_bin_sub_param_1 + 1);    
  fmax_sub_param_1 = hchisq_marg_1D->GetBinLowEdge(fmax_bin_sub_param_1 + 1);    

  ofstream ostext(outname_range,ios::out);
  ostext << "Used sub-range of the free parameter: " << fnamefreeparam_1 << endl;
  ostext << "First bin of free param: " << fmin_bin_sub_param_1 << endl;
  ostext << "Last bin of free param: " << fmax_bin_sub_param_1 << endl;  
  ostext << endl;
  ostext << "fnamefreeparam_1: " << endl
	 << " - " << soutname << endl
	 << " - " << soutname_range << endl;
  ostext << "Mass hierarchy: "; 
  if( this->doMargMassHierarchy ) ostext << " Marginalize";
  else ostext << util::AsString(Opt::Instance()->gMassHierarchyAssumed) << endl;
  ostext << endl;
  if(fRethrowSystUnphysical) ostext << "Rethrow unphysical systematic toys" << endl;
  else ostext << "Don't rethrow unphysical systematic toys" << endl;
  ostext << "fnbins_param_1 = " << fnbins_param_1 << endl;
  ostext << "fntoys_marg = " << fntoys_marg << endl;
  ostext << "fnerecobins_RHCNumu = " << fnerecobins_RHCNumu << endl;
  ostext << "fmin_bin_sub_param_1 = " << fmin_bin_sub_param_1 << endl;
  ostext << "fmax_bin_sub_param_1 = " << fmax_bin_sub_param_1 << endl;
  ostext << "fmin_sub_param_1 = " << fmin_sub_param_1 << endl;
  ostext << "fmax_sub_param_1 = " << fmax_sub_param_1 << endl;
  ostext << "fmin_param_1 = " << fmin_param_1 << endl;
  ostext << "fmax_param_1 = " << fmax_param_1 << endl;
  ostext << "fstep_param_1 = " << fstep_param_1 << endl;
  ostext << endl;
  ostext << "fhypothesis = " << util::AsString(fhypothesis) << endl;
  ostext << endl;
  ostext << "FitParamList:" << endl;
  ostext << *fparm << endl;
  
  ostext.close();
}



//________________________________________________________________________________________
void FitMarginalization::WriteToys2BinaryMerge(char *tag){

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
  sprintf(outname,fnamefilebinary,outtag,fNmergedFiles,fntoys_marg,fnbins_param_1,namemh,"",".bin");

  std::string soutname(outname);
  fToyMargErecoVArray_1D_RHCNumu_Clone->WriteFile(soutname);

  // Write store sub-range of the free parameter

  int min_bin_sub_param = fmin_bin_sub_param_1;
  int max_bin_sub_param = fmax_bin_sub_param_1;

  char outname_range[200];
  sprintf(outname_range,fnamefilebinary,outtag,fNmergedFiles,fntoys_marg,fnbins_param_1,namemh,"_subrange",".bin");
  
  std::string soutname_range(outname_range);
  
  int min_bin = 0;
  int max_bin = fnbins_param_1 - 1;

  ofstream os(soutname_range.c_str(),ios::binary|ios::out);
  assert(os);
  os.write( (char*)&min_bin,sizeof(int));
  os.write( (char*)&max_bin,sizeof(int));
  os.close();

  cout << " - Free param sub-range for marginalization is stored in " << soutname_range << endl;
  cout << endl;

  //

  char outname_text[200];
  std::string soutname_text(outname_text);
  sprintf(outname_text,fnamefilebinary,outtag,fNmergedFiles,fntoys_marg,fnbins_param_1,namemh,"",".txt");

  ofstream ostext(outname_text,ios::out);
  ostext << "fnamefreeparam_1: " << fnamefreeparam_1 << endl;
  ostext << "Mass hierarchy: "; 
  if( this->doMargMassHierarchy ) ostext << " Marginalize";
  else ostext << util::AsString(Opt::Instance()->gMassHierarchyAssumed) << endl;
  ostext << endl;
  if(fRethrowSystUnphysical) ostext << "Rethrow unphysical systematic toys" << endl;
  else ostext << "Don't rethrow unphysical systematic toys" << endl;
  ostext << "# of merged files: " << fNmergedFiles << endl;
  ostext << "fnbins_param_1 = " << fnbins_param_1 << endl;
  ostext << "fntoys_marg = " << fntoys_marg << endl;
  ostext << "fnerecobins_RHCNumu = " << fnerecobins_RHCNumu << endl;
  ostext << "min_bin = " << min_bin << endl;
  ostext << "max_bin = " << max_bin << endl;
  ostext << "fmin_param_1 = " << fmin_param_1 << endl;
  ostext << "fmax_param_1 = " << fmax_param_1 << endl;
  ostext << "fstep_param_1 = " << fstep_param_1 << endl;
  ostext << endl;
  ostext << "fhypothesis = " << util::AsString(fhypothesis) << endl;
  ostext << endl;
  ostext << "FitParamList:" << endl;
  ostext << *fparm << endl;
  
  ostext.close();
}
//________________________________________________________________________________________
void FitMarginalization::ReadToysFromBinary(int idfile, char *tag){ 

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
  sprintf(inname,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param_1,namemh,"",".bin");

  std::string sinname(inname);
  cout << sinname << endl;

  fToyMargErecoVArray_1D_RHCNumu->ReadFile(sinname);

  cout << endl;
  cout << "Read toys for marginalization from binary file:" << sinname << endl;  
  
  // Write store sub-range of the free parameter

  int min_bin_sub_param=0;
  int max_bin_sub_param=0;

  char inname_range[200];
  sprintf(inname_range,fnamefilebinary,tag,idfile,fntoys_marg,fnbins_param_1,namemh,"_subrange",".bin");

  std::string sinname_range(inname_range);

  cout << "Read free param sub-range from binary file:" << sinname_range << endl;

  ifstream is(sinname_range.c_str(),ios::binary|ios::in);
  assert(is);
  is.read( (char*)&min_bin_sub_param,sizeof(int));
  is.read( (char*)&max_bin_sub_param,sizeof(int));
  is.close();

  fmin_bin_sub_param_1 = min_bin_sub_param;
  fmax_bin_sub_param_1 = max_bin_sub_param;

  fmin_sub_param_1 = hchisq_marg_1D->GetBinLowEdge(fmin_bin_sub_param_1 + 1);    
  fmax_sub_param_1 = hchisq_marg_1D->GetBinLowEdge(fmax_bin_sub_param_1 + 1);
    
  if(fmin_sub_param_1<fmin_param_1 || fmax_sub_param_1>fmax_param_1){
    cerr << endl 
	 << "The free parameter sub-range ["
	 << fmin_sub_param_1 << "," << fmax_sub_param_1 << "] " 
	 << "provided by the files"
	 << inname << endl
	 << "and" 
	 << inname_range << endl
	 << "exceeds the total range set for" 
	 << fnamefreeparam_1
	 << " --> [" << fmin_param_1 << "," << fmax_param_1 << "]"
	 << endl;
    exit(1);
  }
  
  cout << "File sub-range: ["
       << fmin_sub_param_1 << "," << fmax_sub_param_1 
       << "]" << endl;
  cout << "File bin sub-range: ["
       << fmin_bin_sub_param_1 << "," << fmax_bin_sub_param_1 
       << "]" << endl << endl;
 
  return;
}
//________________________________________________________________________________________
void FitMarginalization::Minimize(){
  
  int minbin = -999;

  fchisq_bf = hchisq_marg_1D->GetMinimum();
  minbin = hchisq_marg_1D->GetMinimumBin();
  fparam_1_bf = hchisq_marg_1D->GetBinLowEdge(minbin);;

  return;
}
//________________________________________________________________________________________
void FitMarginalization::Minimize(TH1F *hInChisqMarg){
  
  int minbin = -999;
  
  fchisq_bf = hInChisqMarg->GetMinimum();
  minbin = hInChisqMarg->GetMinimumBin();
  fparam_1_bf = hInChisqMarg->GetBinLowEdge(minbin);
  
  return;
}
//________________________________________________________________________________________
void FitMarginalization::PrintBestFit(void){

  cout << endl;
  cout << "Best fit ChisqMarg = " << this->GetChisqBestFit() << endl;
  cout << "Best fit " << fnamefreeparam_1 << " = " << this->GetParamBestFit() << endl;
  cout << endl;

  return;
}
//________________________________________________________________________________________
VALORArrayF * FitMarginalization::GetToyMargErecoVArray(SuperkSample_t s){

  if(s==kSKsRHCFCFV1RingMuLike) return fToyMargErecoVArray_1D_RHCNumu;
  else{
    cout << "You have chosen a wrong Superk Sample" << endl;
    exit(1);
  }
}
//________________________________________________________________________________________
TH1F * FitMarginalization::GetToyMargHisto(SuperkSample_t s,int itoy,int ifreepar,int iMassHierarchy){

  if(s==kSKsRHCFCFV1RingMuLike){ 
    for(int ibin=0;ibin<util::NumErecoBins(kSKsRHCFCFV1RingMuLike);ibin++){
      double bincont = fToyMargErecoVArray_1D_RHCNumu->GetEntry(itoy,ifreepar,iMassHierarchy,ibin);
      hToyMargCurr_1D->SetBinContent(ibin+1,bincont);
      //cout << bincont << " ";
    }
    //cout << endl;
    return hToyMargCurr_1D;
  }
  else{
    cout << "You have chosen a wrong Superk Sample" << endl;
    exit(1);
  }
}
//________________________________________________________________________________________
void FitMarginalization::ResetChisqMarg(int min_bin_sub_param,int max_bin_sub_param){
  for(int ipar=fmin_bin_sub_param_1;ipar<=fmax_bin_sub_param_1;ipar++){
    double bincont = 0;
    fChisqMargVArray_1D->SetEntry(bincont,ipar);
    hchisq_marg_1D->SetBinContent(ipar+1,bincont);
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
    
    if(p!=fPhysParamFree_1) 
      fparm->SetDoMarg(fPhysParamFree_1,true); 
    else{
      fparm->SetDoMarg(fPhysParamFree_1,false);  // Do not marginalize if free parameter
      continue; 
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
    fToyMargErecoVArray_1D_RHCNumu_Clone = new VALORArrayF(fntoys_marg, fnbins_param_1, nmh, util::NumErecoBins(kSKsRHCFCFV1RingMuLike));
    
    fToyMargErecoVArray_1D_RHCNumu_Clone->SetAll(-1.);

    MergeInputFiles_Called = false;
  }

  cout << endl;
  cout << "Merge toys of input file:" << endl
       << "bin: ["  
       << fmin_bin_sub_param_1
       << ","
       << fmax_bin_sub_param_1
       << "]"
       << endl
       << fnamefreeparam_1 << ": " 
       << "[" 
       << hchisq_marg_1D->GetBinLowEdge(fmin_bin_sub_param_1 + 1)
       << ","
       << hchisq_marg_1D->GetBinLowEdge(fmax_bin_sub_param_1 + 1)
       << "]"
       << endl;

  // Copy read array elements to the clone one. Same loop as in FillChisqMarg()
  for(int iexperiment = 0; iexperiment < fntoys_marg; iexperiment++) {
    for(int ipar=fmin_bin_sub_param_1;ipar<=fmax_bin_sub_param_1;ipar++){
      int nMassHierarchy = 2;
      for(int iMH=0;iMH<nMassHierarchy;iMH++){
	for(int ibin=0;ibin<util::NumErecoBins(kSKsRHCFCFV1RingMuLike);ibin++){
	  
	  float bincont = fToyMargErecoVArray_1D_RHCNumu->GetEntry(iexperiment,ipar,iMH,ibin); // bincont = -1 by default
	  if(bincont >= 0.){
	    float bincont_clone = fToyMargErecoVArray_1D_RHCNumu_Clone->GetEntry(iexperiment,ipar,iMH,ibin); 
	    if(bincont_clone >= 0.){
	      cerr << endl;
	      cerr << "Following bin has been already filled:" << endl;
	      cerr << " - iexperiment = " << iexperiment << endl;
	      cerr << " - ipar = " << ipar << endl;
	      cerr << " - iMH = " << iMH << endl;
	      cerr << " - iereco = " << ibin << endl;
	      exit(1);
	    }
	    else{
	      fToyMargErecoVArray_1D_RHCNumu_Clone->SetEntry(bincont,iexperiment,ipar,iMH,ibin); 
	    }

	    //cout << fToyMargErecoVArray_1D_RHCNumu_Clone->GetEntry(iexperiment,ipar,iMH,ibin) << ", ";
	  }
	  
	} // loop ibin


#ifdef __print_array_toy_merge__
	cout << endl;
	cout << "toy=" << iexperiment << " - ipar=" << ipar << " : ";
	for(int ib = 0; ib < util::NumErecoBins(kSKsRHCFCFV1RingMuLike); ib++){
	  double bincont = fToyMargErecoVArray_1D_RHCNumu_Clone->GetEntry(iexperiment,ipar,iMH,ib);
	  cout << bincont << " "; 
	}
	cout << endl;	
#endif


      } // loop iMH
    } // loop sub-range
  } // loop iexperiment
    
  fNmergedFiles++;

  return;
}
//________________________________________________________________________________________

