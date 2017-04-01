//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis) 
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 5-1-2014  Nick
 Added reverse horn current fluxes.

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__  
//#define __show_flux_twk_mesg_info__   
//#define __show_flux_twk_mesg_debug__   
//#define __show_flux_twk_mesg_debug_level2__  //even more verbose
//#define __show_banffsystconv_check_nuflux_mesg_debug__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <TVirtualPad.h>
#include <TLegend.h>

#include "NuFlux.h" 
#include "BANFF.h" 
#include "Utilities.h" 

using std::cout;
using std::cerr;
using std::setw;
using std::endl;

using namespace valor;
using namespace valor::util;
using namespace valor::conv;

ClassImp(NuFlux)

//________________________________________________________________________________________
NuFlux::NuFlux()
{
  this->Init();

  //check that BANFFSystConv.h is setup correctly
  CheckListOfBANFFParams(kNumOfParamsInBANFFfit, kNumOfFluxParamsInBANFFfit, kListOfBANFFParams,
			 "kListOfBANFFParams", false);
  CheckListOfBANFFParams(kNumOfRHCParamsInBANFFfit, kNumOfRHCFluxParamsInBANFFfit, kListOfRHCBANFFParams,
			 "kListOfRHCBANFFParams", false);
  CheckListOfBANFFParams(kNumOfFHCpRHCParamsInBANFFfit, kNumOfFHCpRHCFluxParamsInBANFFfit, kListOfFHCpRHCBANFFParams,
			 "kListOfFHCpRHCBANFFParams", true);
  
  CheckIsFluxParam(kNumOfParamsInBANFFfit, kNumOfFluxParamsInBANFFfit, kListOfBANFFParams,
		   kParamInBANFFfitIsSuperkFluxParam, "kParamInBANFFfitIsSuperkFluxParam");
  CheckIsFluxParam(kNumOfRHCParamsInBANFFfit, kNumOfRHCFluxParamsInBANFFfit, kListOfRHCBANFFParams,
		   kRHCParamInBANFFfitIsSuperkFluxParam, "kRHCParamInBANFFfitIsSuperkFluxParam");
  CheckIsFluxParam(kNumOfFHCpRHCParamsInBANFFfit, kNumOfFHCpRHCFluxParamsInBANFFfit, kListOfFHCpRHCBANFFParams,
		   kFHCpRHCParamInBANFFfitIsSuperkFluxParam, "kFHCpRHCParamInBANFFfitIsSuperkFluxParam");

  CheckFluxPdg(kNumOfParamsInBANFFfit, kNumOfFluxParamsInBANFFfit, kListOfBANFFParams,
	       kNuPdgCodeOfFluxParamInBANFFfit, "kNuPdgCodeOfFluxParamInBANFFfit");
  CheckFluxPdg(kNumOfRHCParamsInBANFFfit, kNumOfRHCFluxParamsInBANFFfit, kListOfRHCBANFFParams,
	       kNuPdgCodeOfRHCFluxParamInBANFFfit, "kNuPdgCodeOfRHCFluxParamInBANFFfit");
  CheckFluxPdg(kNumOfFHCpRHCParamsInBANFFfit, kNumOfFHCpRHCFluxParamsInBANFFfit, kListOfFHCpRHCBANFFParams,
	       kNuPdgCodeOfFHCpRHCFluxParamInBANFFfit, "kNuPdgCodeOfFHCpRHCFluxParamInBANFFfit");

  CheckFluxERange(kNumOfParamsInBANFFfit, kNumOfFluxParamsInBANFFfit, kListOfBANFFParams,
		  (const double *)kEtrueRangeOfFluxParamInBANFFfit, "kEtrueRangeOfFluxParamInBANFFfit");
  CheckFluxERange(kNumOfRHCParamsInBANFFfit, kNumOfRHCFluxParamsInBANFFfit, kListOfRHCBANFFParams,
		  (const double *)kEtrueRangeOfRHCFluxParamInBANFFfit, "kEtrueRangeOfRHCFluxParamInBANFFfit");
  CheckFluxERange(kNumOfFHCpRHCParamsInBANFFfit, kNumOfFHCpRHCFluxParamsInBANFFfit, kListOfFHCpRHCBANFFParams,
		  (const double *)kEtrueRangeOfFHCpRHCFluxParamInBANFFfit, "kEtrueRangeOfFHCpRHCFluxParamInBANFFfit");

}
//________________________________________________________________________________________
NuFlux::~NuFlux()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void NuFlux::CheckListOfBANFFParams(int nparams, int nfluxparams, const conv::PhysParamId_t * params, const char * vname, bool allow_gap)
{
  int ib = 0;
  PhysParamId_t curr, prev;

  //check that BANFF flux parameters are monotonically increasing (a single gap is allowed in FHCpRHC)
  for(ib = 1; ib < nparams; ib++) {
    prev = params[ib-1];
    curr = params[ib];
#ifdef __show_banffsystconv_check_nuflux_mesg_debug__
    if(ib == 1)
      cout << "NuFlux/Init >> " << vname << "[" << ib-1 << "] = " << prev << " " << util::AsString(prev) << endl;
    cout << "NuFlux/Init >> " << vname << "[" << ib   << "] = " << curr << " " << util::AsString(curr) << endl;
#endif
    if(ib < nfluxparams) {
      if((int)curr != ((int)prev + 1)) {
	//this is to allow jumps between FHC + RHC flux parameters
	if(allow_gap) {
	  allow_gap = false;
	  continue;
	}
	cerr << "NuFlux/Init >> " << vname << "[" << ib   << "] = " << curr << " " << util::AsString(curr) << endl
	     << " comes before " << vname << "[" << ib-1 << "] = " << prev << " " << util::AsString(prev)
	     << endl;
	exit(-1);
      }
    }
  }//ib
}
//________________________________________________________________________________________
void NuFlux::CheckIsFluxParam(int nparams, int nfluxparams, const PhysParamId_t * params, const bool * flux, const char * vname)
{
  PhysParamId_t curr;
  //check that BANFF parameters mark themselves as flux/not correctly
  for(int ib = 0; ib < nparams; ib++) {
    curr = params[ib];
    bool is_flux_expect = (ib < nfluxparams ? true : false);
    bool is_flux  = flux[ib];
#ifdef __show_banffsystconv_check_nuflux_mesg_debug__
    cout << "NuFlux/Init >> " << vname << "[" << ib   << "] = " << is_flux << " " << util::AsString(curr) << endl;
#endif
    if(is_flux != is_flux_expect) {
      cerr << "NuFlux/Init >> " << vname << "[" << ib   << "] = " << is_flux << " " << util::AsString(curr) << endl
	   << "INCORRECT!" << endl;
      exit(-1);
    }
  }//ib
}
//________________________________________________________________________________________
void NuFlux::CheckFluxPdg(int nparams, int nfluxparams, const PhysParamId_t * params, const int * pdg, const char * vname)
{
  PhysParamId_t curr, prev;
  //check that BANFF parameters effect the correct PDG
  for(int ib = 1; ib < nparams; ib++) {
    prev = params[ib-1];
    curr = params[ib];
    int this_flux_pdg = pdg[ib];
    int prev_flux_pdg = pdg[ib-1];
#ifdef __show_banffsystconv_check_nuflux_mesg_debug__
    if(ib == 1)
      cout << "NuFlux/Init >> " << vname << "[" << ib-1 << "] = " << prev_flux_pdg << " " << util::AsString(prev) << endl;
    cout << "NuFlux/Init >> " << vname << "[" << ib   << "] = " << this_flux_pdg << " " << util::AsString(curr) << endl;
#endif
    if(ib < nfluxparams) {
      int this_flux_position = std::distance(kOrderOfBANFFFluxPDG, std::find(kOrderOfBANFFFluxPDG, kOrderOfBANFFFluxPDG + kNBANFFluxPDGs, this_flux_pdg));
      int prev_flux_position = std::distance(kOrderOfBANFFFluxPDG, std::find(kOrderOfBANFFFluxPDG, kOrderOfBANFFFluxPDG + kNBANFFluxPDGs, prev_flux_pdg));
      int next_flux_position = (prev_flux_position + 1) % 4;
      if(this_flux_position > kNBANFFluxPDGs) {
	cerr << "NuFlux/Init >> Unknown flux PDG code" << endl
	     << vname << "[" << ib   << "] = " << this_flux_pdg << " " << util::AsString(curr) << endl;
	exit(-1);
      }
      if(prev_flux_position > kNBANFFluxPDGs) {
	cerr << "NuFlux/Init >> Unknown flux PDG code" << endl
	     << vname << "[" << ib-1 << "] = " << prev_flux_pdg << " " << util::AsString(prev) << endl;
	exit(-1);
      }
      if(this_flux_position != prev_flux_position && this_flux_position != next_flux_position){
	cerr << "NuFlux/Init >> Wrong flux PDG code" << endl
	     << vname << "[" << ib   << "] = " << this_flux_pdg << " " << util::AsString(curr) << endl
	     << "Should be " << pdg[prev_flux_position]
	     << " or "       << pdg[next_flux_position] << endl;
	exit(-1);
      }
    }//flux param
    else {
      if(this_flux_pdg) {
	cerr << "NuFlux/Init >> Wrong flux PDG code" << endl
	     << vname << "[" << ib   << "] = " << this_flux_pdg << " " << util::AsString(curr) << endl
	     << "Should be 0" << endl;
	exit(-1);
      }
    }//not flux param
  }//ib
}
//________________________________________________________________________________________
void NuFlux::CheckFluxERange(int nparams, int nfluxparams, const PhysParamId_t * params, const double * erange, const char * vname)
{
  PhysParamId_t curr, prev;
  //check that BANFF parameters effect the correct PDG
  for(int ib = 1; ib < nparams; ib++) {
    prev = params[ib-1];
    curr = params[ib];
    double this_lo = erange[ib*2];
    double this_hi = erange[ib*2 + 1];
    double prev_lo = erange[(ib-1)*2];
    double prev_hi = erange[(ib-1)*2 + 1];
#ifdef __show_banffsystconv_check_nuflux_mesg_debug__
    if(ib == 1)
      cout << "NuFlux/Init >> " << vname << "[" << ib-1 << "] = " << prev_lo << "," << prev_hi << " " << util::AsString(prev) << endl;
    cout << "NuFlux/Init >> " << vname << "[" << ib << "] = " << this_lo << "," << this_hi << " " << util::AsString(curr) << endl;
#endif
    if(ib < nfluxparams) {
      if(TMath::Abs(this_lo - prev_hi) > 1E-6) {
	if(TMath::Abs(prev_hi - 30) < 1E-6 && TMath::Abs(this_lo) < 1E-6) {} //it's ok!
	else {
	  cerr << "NuFlux/Init >> Wrong energy range" << endl
	       << vname << "[" << ib-1 << "] = " << prev_lo << "," << prev_hi << " " << util::AsString(prev) << endl
	       << vname << "[" << ib << "] = " << this_lo << "," << this_hi << " " << util::AsString(curr) << endl
	       << "Should be previous high and this low should be identical (or 30 and 0)" << endl;
	  exit(-1);
	}
      }
      if(TMath::Abs(this_lo + 1) < 1E-6 || TMath::Abs(this_hi + 1) < 1E-6) {
	cerr << "NuFlux/Init >> Wrong energy range" << endl
	     << vname << "[" << ib << "] = " << this_lo << "," << this_hi << " " << util::AsString(curr) << endl
	     << "Flux parameters shouldn't have -1 here" << endl;	
      }
    }
    else {
      if(TMath::Abs(this_lo + 1) > 1E-6 || TMath::Abs(this_hi + 1) > 1E-6) {
	cerr << "NuFlux/Init >> Wrong energy range" << endl
	     << vname << "[" << ib << "] = " << this_lo << "," << this_hi << " " << util::AsString(curr) << endl
	     << "Should be -1,-1" << endl;
	exit(-1);
      }
    }
  }//ib
}
//________________________________________________________________________________________
void NuFlux::Init(void)
{
  fSKNuMuFluxNominal          = 0; 
  fSKNuMuBarFluxNominal       = 0; 
  fSKNuEFluxNominal           = 0; 
  fSKNuEBarFluxNominal        = 0; 
  fSKNuMuFluxTunedExt         = 0; 
  fSKNuMuBarFluxTunedExt      = 0; 
  fSKNuEFluxTunedExt          = 0; 
  fSKNuEBarFluxTunedExt       = 0; 
  fSKNuMuFluxTunedFinal       = 0; 
  fSKNuMuBarFluxTunedFinal    = 0; 
  fSKNuEFluxTunedFinal        = 0; 
  fSKNuEBarFluxTunedFinal     = 0; 
  fSKNuMuFluxCurrent          = 0; 
  fSKNuMuBarFluxCurrent       = 0; 
  fSKNuEFluxCurrent           = 0; 
  fSKNuEBarFluxCurrent        = 0; 
  fSKNuMuFluxExtTuningWght    = 0;
  fSKNuMuBarFluxExtTuningWght = 0;
  fSKNuEFluxExtTuningWght     = 0;
  fSKNuEBarFluxExtTuningWght  = 0;
  fBANFF                      = 0;
  fSKRHCNuMuFluxNominal          = 0; 
  fSKRHCNuMuBarFluxNominal       = 0; 
  fSKRHCNuEFluxNominal           = 0; 
  fSKRHCNuEBarFluxNominal        = 0; 
  fSKRHCNuMuFluxTunedExt         = 0; 
  fSKRHCNuMuBarFluxTunedExt      = 0; 
  fSKRHCNuEFluxTunedExt          = 0; 
  fSKRHCNuEBarFluxTunedExt       = 0; 
  fSKRHCNuMuFluxTunedFinal       = 0; 
  fSKRHCNuMuBarFluxTunedFinal    = 0; 
  fSKRHCNuEFluxTunedFinal        = 0; 
  fSKRHCNuEBarFluxTunedFinal     = 0; 
  fSKRHCNuMuFluxCurrent          = 0; 
  fSKRHCNuMuBarFluxCurrent       = 0; 
  fSKRHCNuEFluxCurrent           = 0; 
  fSKRHCNuEBarFluxCurrent        = 0; 
  fSKRHCNuMuFluxExtTuningWght    = 0;
  fSKRHCNuMuBarFluxExtTuningWght = 0;
  fSKRHCNuEFluxExtTuningWght     = 0;
  fSKRHCNuEBarFluxExtTuningWght  = 0;
}
//________________________________________________________________________________________
void NuFlux::CleanUp(void)
{
  if ( fSKNuMuFluxNominal          ) delete fSKNuMuFluxNominal;
  if ( fSKNuMuBarFluxNominal       ) delete fSKNuMuBarFluxNominal;
  if ( fSKNuEFluxNominal           ) delete fSKNuEFluxNominal;
  if ( fSKNuEBarFluxNominal        ) delete fSKNuEBarFluxNominal;
  if ( fSKNuMuFluxTunedExt         ) delete fSKNuMuFluxTunedExt;
  if ( fSKNuMuBarFluxTunedExt      ) delete fSKNuMuBarFluxTunedExt;
  if ( fSKNuEFluxTunedExt          ) delete fSKNuEFluxTunedExt;
  if ( fSKNuEBarFluxTunedExt       ) delete fSKNuEBarFluxTunedExt;
  if ( fSKNuMuFluxTunedFinal       ) delete fSKNuMuFluxTunedFinal;
  if ( fSKNuMuBarFluxTunedFinal    ) delete fSKNuMuBarFluxTunedFinal;
  if ( fSKNuEFluxTunedFinal        ) delete fSKNuEFluxTunedFinal;
  if ( fSKNuEBarFluxTunedFinal     ) delete fSKNuEBarFluxTunedFinal;
  if ( fSKNuMuFluxCurrent          ) delete fSKNuMuFluxCurrent;
  if ( fSKNuMuBarFluxCurrent       ) delete fSKNuMuBarFluxCurrent;
  if ( fSKNuEFluxCurrent           ) delete fSKNuEFluxCurrent;
  if ( fSKNuEBarFluxCurrent        ) delete fSKNuEBarFluxCurrent;
  if ( fSKNuMuFluxExtTuningWght    ) delete fSKNuMuFluxExtTuningWght;
  if ( fSKNuMuBarFluxExtTuningWght ) delete fSKNuMuBarFluxExtTuningWght;
  if ( fSKNuEFluxExtTuningWght     ) delete fSKNuEFluxExtTuningWght;
  if ( fSKNuEBarFluxExtTuningWght  ) delete fSKNuEBarFluxExtTuningWght;
  if ( fSKRHCNuMuFluxNominal          ) delete fSKRHCNuMuFluxNominal;
  if ( fSKRHCNuMuBarFluxNominal       ) delete fSKRHCNuMuBarFluxNominal;
  if ( fSKRHCNuEFluxNominal           ) delete fSKRHCNuEFluxNominal;
  if ( fSKRHCNuEBarFluxNominal        ) delete fSKRHCNuEBarFluxNominal;
  if ( fSKRHCNuMuFluxTunedExt         ) delete fSKRHCNuMuFluxTunedExt;
  if ( fSKRHCNuMuBarFluxTunedExt      ) delete fSKRHCNuMuBarFluxTunedExt;
  if ( fSKRHCNuEFluxTunedExt          ) delete fSKRHCNuEFluxTunedExt;
  if ( fSKRHCNuEBarFluxTunedExt       ) delete fSKRHCNuEBarFluxTunedExt;
  if ( fSKRHCNuMuFluxTunedFinal       ) delete fSKRHCNuMuFluxTunedFinal;
  if ( fSKRHCNuMuBarFluxTunedFinal    ) delete fSKRHCNuMuBarFluxTunedFinal;
  if ( fSKRHCNuEFluxTunedFinal        ) delete fSKRHCNuEFluxTunedFinal;
  if ( fSKRHCNuEBarFluxTunedFinal     ) delete fSKRHCNuEBarFluxTunedFinal;
  if ( fSKRHCNuMuFluxCurrent          ) delete fSKRHCNuMuFluxCurrent;
  if ( fSKRHCNuMuBarFluxCurrent       ) delete fSKRHCNuMuBarFluxCurrent;
  if ( fSKRHCNuEFluxCurrent           ) delete fSKRHCNuEFluxCurrent;
  if ( fSKRHCNuEBarFluxCurrent        ) delete fSKRHCNuEBarFluxCurrent;
  if ( fSKRHCNuMuFluxExtTuningWght    ) delete fSKRHCNuMuFluxExtTuningWght;
  if ( fSKRHCNuMuBarFluxExtTuningWght ) delete fSKRHCNuMuBarFluxExtTuningWght;
  if ( fSKRHCNuEFluxExtTuningWght     ) delete fSKRHCNuEFluxExtTuningWght;
  if ( fSKRHCNuEBarFluxExtTuningWght  ) delete fSKRHCNuEBarFluxExtTuningWght;
}
//________________________________________________________________________________________
bool NuFlux::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read flux inputs to the oscillation fit
//
#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Reading flux inputs from input XML file." << endl;
#endif

  //
  // Load fluxes tuned based on external data only
  //

  // SuperK numu flux
  fSKNuMuFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned/numu");
  if(!fSKNuMuFluxTunedExt) {
    cerr << "Exp/Flux >> No input SuperK numu ext/tuned flux histogram!" << endl;
    exit(1);
  }
  // SuperK numubar flux
  fSKNuMuBarFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned/numubar");
  if(!fSKNuMuBarFluxTunedExt) {
    cerr << "Exp/Flux >> No input SuperK numubar ext/tuned flux histogram!" << endl;
    exit(1);
  }
  // SuperK nue flux
  fSKNuEFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned/nue");
  if(!fSKNuEFluxTunedExt) {
    cerr << "Exp/Flux >> No input SuperK nue ext/tuned flux histogram!" << endl;
    exit(1);
  }
  // SuperK nuebar flux
  fSKNuEBarFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned/nuebar");
  if(!fSKNuEBarFluxTunedExt) {
    cerr << "Exp/Flux >> No input SuperK nuebar ext/tuned flux histogram!" << endl;
    exit(1);
  }

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Done loading ext/tuned flux histograms" << endl;
#endif

  //
  // Load nominal -> ext.-tuned flux weights
  //

  // SuperK numu flux tuning weights
  fSKNuMuFluxExtTuningWght = 
     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/numu");
  if(!fSKNuMuFluxExtTuningWght) {
    cerr << "Exp/Flux >> No input SuperK numu flux ext/tuning weights!" << endl;
    exit(1);
  }
  // SuperK numubar flux tuning weights
  fSKNuMuBarFluxExtTuningWght = 
     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/numubar");
  if(!fSKNuMuBarFluxExtTuningWght) {
    cerr << "Exp/Flux >> No input SuperK numubar flux ext/tuning weights!" << endl;
    exit(1);
  }
  // SuperK nue flux tuning weights
  fSKNuEFluxExtTuningWght = 
     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/nue");
  if(!fSKNuEFluxExtTuningWght) {
    cerr << "Exp/Flux >> No input SuperK nue flux ext/tuning weights!" << endl;
    exit(1);
  }
  // SuperK nuebar flux tuning weights
  fSKNuEBarFluxExtTuningWght = 
     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/nuebar");
  if(!fSKNuEBarFluxExtTuningWght) {
    cerr << "Exp/Flux >> No input SuperK nuebar flux ext/tuning weights!" << endl;
    exit(1);
  }

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Done loading flux tuning weights" << endl;
#endif

  //
  // Book and initialize other auxiliary flux histograms not loaded from files
  //
  this->BookAndInitAuxSpectra();

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Done reading flux inputs!" << endl;
#endif

  return true;
}
//________________________________________________________________________________________
bool NuFlux::ReadRHC(TXMLEngine* xml, XMLNodePointer_t root_node)
{
  // Read flux inputs to the oscillation fit
  //
#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Reading RHC flux inputs from input XML file." << endl;
#endif

  //
  // Load fluxes tuned based on external data only
  //

  // SuperK numu flux (reverse horn current)
  fSKRHCNuMuFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned_rhc/numu");
  if(!fSKRHCNuMuFluxTunedExt) {
    cerr << "Exp/Flux >> No input SuperK RHC numu ext/tuned flux histogram!" << endl;
    exit(1);
  }
  // SuperK numubar flux (reverse horn current)
  fSKRHCNuMuBarFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned_rhc/numubar");
  if(!fSKRHCNuMuBarFluxTunedExt) {
    cerr << "Exp/Flux >> No input SuperK RHC numubar ext/tuned flux histogram!" << endl;
    exit(1);
  }
  // SuperK nue flux (reverse horn current)
  fSKRHCNuEFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned_rhc/nue");
  if(!fSKRHCNuEFluxTunedExt) {
    cerr << "Exp/Flux >> No input SuperK RHC nue ext/tuned flux histogram!" << endl;
    exit(1);
  }
  // SuperK nuebar flux (reverse horn current)
  fSKRHCNuEBarFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned_rhc/nuebar");
  if(!fSKRHCNuEBarFluxTunedExt) {
    cerr << "Exp/Flux >> No input SuperK RHC nuebar ext/tuned flux histogram!" << endl;
    exit(1);
  }

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Done loading RHC ext/tuned flux histograms" << endl;
#endif

  //
  // Load nominal -> ext.-tuned flux weights
  //

  // SuperK numu flux tuning weights (reverse horn current) 
  fSKRHCNuMuFluxExtTuningWght = 
    util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights_rhc/numu");
  if(!fSKRHCNuMuFluxExtTuningWght) {
    cerr << "Exp/Flux >> No input SuperK RHC numu flux ext/tuning weights!" << endl;
    exit(1);
  }
  // SuperK numubar flux tuning weights (reverse horn current) 
  fSKRHCNuMuBarFluxExtTuningWght = 
    util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights_rhc/numubar");
  if(!fSKRHCNuMuBarFluxExtTuningWght) {
    cerr << "Exp/Flux >> No input SuperK RHC numubar flux ext/tuning weights!" << endl;
    exit(1);
  }
  // SuperK nue flux tuning weights (reverse horn current) 
  fSKRHCNuEFluxExtTuningWght = 
    util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights_rhc/nue");
  if(!fSKRHCNuEFluxExtTuningWght) {
    cerr << "Exp/Flux >> No input SuperK RHC nue flux ext/tuning weights!" << endl;
    exit(1);
  }
  // SuperK nuebar flux tuning weights (reverse horn current) 
  fSKRHCNuEBarFluxExtTuningWght = 
    util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights_rhc/nuebar");
  if(!fSKRHCNuEBarFluxExtTuningWght) {
    cerr << "Exp/Flux >> No input SuperK RHC nuebar flux ext/tuning weights!" << endl;
    exit(1);
  }

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Done loading RHC flux tuning weights" << endl;
#endif

  //
  // Book and initialize other auxiliary flux histograms not loaded from files
  //
  this->BookAndInitRHCAuxSpectra();

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Done reading RHC flux inputs!" << endl;
#endif

  return true;
}
//________________________________________________________________________________________











// NEW FLUX FHC+RHC

// bool NuFlux::ReadFHCpRHC(TXMLEngine* xml, XMLNodePointer_t root_node)
// {
//   // Read flux inputs to the oscillation fit
//   //
// #ifdef __show_init_mesg_info__
//   cout << "Exp/Flux >> Reading FHC+RHC flux inputs from input XML file." << endl;
// #endif

//   //
//   // Load fluxes tuned based on external data only
//   //

//   // SuperK numu flux
//   fSKNuMuFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned/numu");
//   if(!fSKNuMuFluxTunedExt) {
//     cerr << "Exp/Flux >> No input SuperK numu ext/tuned flux histogram!" << endl;
//     exit(1);
//   }
//   // SuperK numubar flux
//   fSKNuMuBarFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned/numubar");
//   if(!fSKNuMuBarFluxTunedExt) {
//     cerr << "Exp/Flux >> No input SuperK numubar ext/tuned flux histogram!" << endl;
//     exit(1);
//   }
//   // SuperK nue flux
//   fSKNuEFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned/nue");
//   if(!fSKNuEFluxTunedExt) {
//     cerr << "Exp/Flux >> No input SuperK nue ext/tuned flux histogram!" << endl;
//     exit(1);
//   }
//   // SuperK nuebar flux
//   fSKNuEBarFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned/nuebar");
//   if(!fSKNuEBarFluxTunedExt) {
//     cerr << "Exp/Flux >> No input SuperK nuebar ext/tuned flux histogram!" << endl;
//     exit(1);
//   }

//   // SuperK numu flux (reverse horn current)
//   fSKRHCNuMuFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned_rhc/numu");
//   if(!fSKRHCNuMuFluxTunedExt) {
//     cerr << "Exp/Flux >> No input SuperK RHC numu ext/tuned flux histogram!" << endl;
//     exit(1);
//   }
//   // SuperK numubar flux (reverse horn current)
//   fSKRHCNuMuBarFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned_rhc/numubar");
//   if(!fSKRHCNuMuBarFluxTunedExt) {
//     cerr << "Exp/Flux >> No input SuperK RHC numubar ext/tuned flux histogram!" << endl;
//     exit(1);
//   }
//   // SuperK nue flux (reverse horn current)
//   fSKRHCNuEFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned_rhc/nue");
//   if(!fSKRHCNuEFluxTunedExt) {
//     cerr << "Exp/Flux >> No input SuperK RHC nue ext/tuned flux histogram!" << endl;
//     exit(1);
//   }
//   // SuperK nuebar flux (reverse horn current)
//   fSKRHCNuEBarFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned_rhc/nuebar");
//   if(!fSKRHCNuEBarFluxTunedExt) {
//     cerr << "Exp/Flux >> No input SuperK RHC nuebar ext/tuned flux histogram!" << endl;
//     exit(1);
//   }

// #ifdef __show_init_mesg_info__
//   cout << "Exp/Flux >> Done loading RHC ext/tuned flux histograms" << endl;
// #endif

//   //
//   // Load nominal -> ext.-tuned flux weights
//   //

//   // SuperK numu flux tuning weights
//   fSKNuMuFluxExtTuningWght = 
//      util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/numu");
//   if(!fSKNuMuFluxExtTuningWght) {
//     cerr << "Exp/Flux >> No input SuperK numu flux ext/tuning weights!" << endl;
//     exit(1);
//   }
//   // SuperK numubar flux tuning weights
//   fSKNuMuBarFluxExtTuningWght = 
//      util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/numubar");
//   if(!fSKNuMuBarFluxExtTuningWght) {
//     cerr << "Exp/Flux >> No input SuperK numubar flux ext/tuning weights!" << endl;
//     exit(1);
//   }
//   // SuperK nue flux tuning weights
//   fSKNuEFluxExtTuningWght = 
//      util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/nue");
//   if(!fSKNuEFluxExtTuningWght) {
//     cerr << "Exp/Flux >> No input SuperK nue flux ext/tuning weights!" << endl;
//     exit(1);
//   }
//   // SuperK nuebar flux tuning weights
//   fSKNuEBarFluxExtTuningWght = 
//      util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/nuebar");
//   if(!fSKNuEBarFluxExtTuningWght) {
//     cerr << "Exp/Flux >> No input SuperK nuebar flux ext/tuning weights!" << endl;
//     exit(1);
//   }

//   // SuperK numu flux tuning weights (reverse horn current) 
//   fSKRHCNuMuFluxExtTuningWght = 
//     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights_rhc/numu");
//   if(!fSKRHCNuMuFluxExtTuningWght) {
//     cerr << "Exp/Flux >> No input SuperK RHC numu flux ext/tuning weights!" << endl;
//     exit(1);
//   }
//   // SuperK numubar flux tuning weights (reverse horn current) 
//   fSKRHCNuMuBarFluxExtTuningWght = 
//     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights_rhc/numubar");
//   if(!fSKRHCNuMuBarFluxExtTuningWght) {
//     cerr << "Exp/Flux >> No input SuperK RHC numubar flux ext/tuning weights!" << endl;
//     exit(1);
//   }
//   // SuperK nue flux tuning weights (reverse horn current) 
//   fSKRHCNuEFluxExtTuningWght = 
//     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights_rhc/nue");
//   if(!fSKRHCNuEFluxExtTuningWght) {
//     cerr << "Exp/Flux >> No input SuperK RHC nue flux ext/tuning weights!" << endl;
//     exit(1);
//   }
//   // SuperK nuebar flux tuning weights (reverse horn current) 
//   fSKRHCNuEBarFluxExtTuningWght = 
//     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights_rhc/nuebar");
//   if(!fSKRHCNuEBarFluxExtTuningWght) {
//     cerr << "Exp/Flux >> No input SuperK RHC nuebar flux ext/tuning weights!" << endl;
//     exit(1);
//   }

// #ifdef __show_init_mesg_info__
//   cout << "Exp/Flux >> Done loading FHC+RHC flux tuning weights" << endl;
// #endif

//   //
//   // Book and initialize other auxiliary flux histograms not loaded from files
//   //
//   this->BookAndInitFHCpRHCAuxSpectra();

// #ifdef __show_init_mesg_info__
//   cout << "Exp/Flux >> Done reading FHC+RHC flux inputs!" << endl;
// #endif

//   return true;
// }
// //






//________________________________________________________________________________________
void NuFlux::BookAndInitAuxSpectra(void)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Booking and initializing/computing auxiliary flux spectra" << endl;
#endif

  //
  // Create nominal flux histograms and initialize them using the ext/tuned flux.
  // Then de-weight them with input ext/tuning weights to recreate the nominal flux.
  //
#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Booking and re-creating nominal flux histograms" << endl;
#endif

  fSKNuMuFluxNominal    = new TH1D(*fSKNuMuFluxTunedExt   );
  fSKNuMuBarFluxNominal = new TH1D(*fSKNuMuBarFluxTunedExt);
  fSKNuEFluxNominal     = new TH1D(*fSKNuEFluxTunedExt    );
  fSKNuEBarFluxNominal  = new TH1D(*fSKNuEBarFluxTunedExt );

  fSKNuMuFluxNominal    -> SetDirectory(0);
  fSKNuMuBarFluxNominal -> SetDirectory(0);
  fSKNuEFluxNominal     -> SetDirectory(0);
  fSKNuEBarFluxNominal  -> SetDirectory(0);

  const int n_nu = 4;
  int nu_pdg[n_nu] = { conv::kPdgNuE, conv::kPdgNuMu, conv::kPdgNuEBar, conv::kPdgNuMuBar };

  for(int inu = 0; inu < n_nu; inu++) {
     TH1D * hnom = this->SuperkFluxSpectrum(nu_pdg[inu], conv::kFlxNominal  );
     for(int ie = 1; ie <= hnom->GetNbinsX(); ie++) {
        double E = hnom->GetBinCenter (ie);
        double f = hnom->GetBinContent(ie);
        double w = this->SuperkFluxExtTuningWght(nu_pdg[inu],E);
        if(w >=0) { f /= w; }
        hnom->SetBinContent(ie,f);
     }//ie
  }//inu

  //
  // Create the fin/tuned flux histograms and initialize them using the ext/tuned flux
  // Then re-weight them using the weights that correspond to the BANFF best-fit params.
  //
  fSKNuMuFluxTunedFinal    = new TH1D(*fSKNuMuFluxTunedExt   );
  fSKNuMuBarFluxTunedFinal = new TH1D(*fSKNuMuBarFluxTunedExt);
  fSKNuEFluxTunedFinal     = new TH1D(*fSKNuEFluxTunedExt    );
  fSKNuEBarFluxTunedFinal  = new TH1D(*fSKNuEBarFluxTunedExt );

  fSKNuMuFluxTunedFinal    -> SetDirectory(0);
  fSKNuMuBarFluxTunedFinal -> SetDirectory(0);
  fSKNuEFluxTunedFinal     -> SetDirectory(0);
  fSKNuEBarFluxTunedFinal  -> SetDirectory(0);

  this->ApplyBANFFWeights(false); // reweight ext/tuned to fin/tuned flux (corresponds to BANFF best-fit)

  //
  // Create current flux histograms and initialize them using the fin/tuned flux.
  //
#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Booking and initializing current flux histograms" << endl;
#endif

  fSKNuMuFluxCurrent    = new TH1D(*fSKNuMuFluxTunedFinal   );
  fSKNuMuBarFluxCurrent = new TH1D(*fSKNuMuBarFluxTunedFinal);
  fSKNuEFluxCurrent     = new TH1D(*fSKNuEFluxTunedFinal    );
  fSKNuEBarFluxCurrent  = new TH1D(*fSKNuEBarFluxTunedFinal );

  fSKNuMuFluxCurrent    -> SetDirectory(0);
  fSKNuMuBarFluxCurrent -> SetDirectory(0);
  fSKNuEFluxCurrent     -> SetDirectory(0);
  fSKNuEBarFluxCurrent  -> SetDirectory(0);
}
//________________________________________________________________________________________
void NuFlux::BookAndInitRHCAuxSpectra(void)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Booking and initializing/computing RHC auxiliary flux spectra" << endl;
#endif

  //
  // Create nominal flux histograms and initialize them using the ext/tuned flux.
  // Then de-weight them with input ext/tuning weights to recreate the nominal flux.
  //

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Booking and recreating RHC nominal flux histograms" << endl;
#endif

  fSKRHCNuMuFluxNominal    = new TH1D(*fSKRHCNuMuFluxTunedExt   );
  fSKRHCNuMuBarFluxNominal = new TH1D(*fSKRHCNuMuBarFluxTunedExt);
  fSKRHCNuEFluxNominal     = new TH1D(*fSKRHCNuEFluxTunedExt    );
  fSKRHCNuEBarFluxNominal  = new TH1D(*fSKRHCNuEBarFluxTunedExt );

  fSKRHCNuMuFluxNominal    -> SetDirectory(0);
  fSKRHCNuMuBarFluxNominal -> SetDirectory(0);
  fSKRHCNuEFluxNominal     -> SetDirectory(0);
  fSKRHCNuEBarFluxNominal  -> SetDirectory(0);

  const int n_nu = 4;
  int nu_pdg[n_nu] = { conv::kPdgNuE, conv::kPdgNuMu, conv::kPdgNuEBar, conv::kPdgNuMuBar };

  for(int inu = 0; inu < n_nu; inu++) {
    TH1D * hnomrhc = this->SuperkFluxSpectrum(nu_pdg[inu], conv::kRHCFlxNominal  );
    for(int ie = 1; ie <= hnomrhc->GetNbinsX(); ie++) {
      double E = hnomrhc->GetBinCenter (ie);
      double f = hnomrhc->GetBinContent(ie);
      double w = this->SuperkRHCFluxExtTuningWght(nu_pdg[inu],E);
      if(w >=0) { f /= w; }
      hnomrhc->SetBinContent(ie,f);
    }//ie
  }//inu

  //
  // Create the fin/tuned flux histograms and initialize them using the ext/tuned flux
  // Then re-weight them using the weights that correspond to the BANFF best-fit params.
  //

  fSKRHCNuMuFluxTunedFinal    = new TH1D(*fSKRHCNuMuFluxTunedExt   );
  fSKRHCNuMuBarFluxTunedFinal = new TH1D(*fSKRHCNuMuBarFluxTunedExt);
  fSKRHCNuEFluxTunedFinal     = new TH1D(*fSKRHCNuEFluxTunedExt    );
  fSKRHCNuEBarFluxTunedFinal  = new TH1D(*fSKRHCNuEBarFluxTunedExt );

  fSKRHCNuMuFluxTunedFinal    -> SetDirectory(0);
  fSKRHCNuMuBarFluxTunedFinal -> SetDirectory(0);
  fSKRHCNuEFluxTunedFinal     -> SetDirectory(0);
  fSKRHCNuEBarFluxTunedFinal  -> SetDirectory(0);

  this->ApplyRHCBANFFWeights(false); // reweight ext/tuned to fin/tuned flux (corresponds to BANFF best-fit)

  //
  // Create current flux histograms and initialize them using the fin/tuned flux.
  //

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Booking and initializing RHC current flux histograms" << endl;
#endif

  fSKRHCNuMuFluxCurrent    = new TH1D(*fSKRHCNuMuFluxTunedFinal   );
  fSKRHCNuMuBarFluxCurrent = new TH1D(*fSKRHCNuMuBarFluxTunedFinal);
  fSKRHCNuEFluxCurrent     = new TH1D(*fSKRHCNuEFluxTunedFinal    );
  fSKRHCNuEBarFluxCurrent  = new TH1D(*fSKRHCNuEBarFluxTunedFinal );

  fSKRHCNuMuFluxCurrent    -> SetDirectory(0);
  fSKRHCNuMuBarFluxCurrent -> SetDirectory(0);
  fSKRHCNuEFluxCurrent     -> SetDirectory(0);
  fSKRHCNuEBarFluxCurrent  -> SetDirectory(0);
}
//________________________________________________________________________________________













// NEW FLUX FHC+RHC

// void NuFlux::BookAndInitFHCpRHCAuxSpectra(void)
// {
// #ifdef __show_init_mesg_info__
//   cout << "Exp/Flux >> Booking and initializing/computing FHC+RHC auxiliary flux spectra" << endl;
// #endif

//   //
//   // Create nominal flux histograms and initialize them using the ext/tuned flux.
//   // Then de-weight them with input ext/tuning weights to recreate the nominal flux.
//   //

// #ifdef __show_init_mesg_info__
//   cout << "Exp/Flux >> Booking and recreating RHC nominal flux histograms" << endl;
// #endif

//   fSKNuMuFluxNominal    = new TH1D(*fSKNuMuFluxTunedExt   );
//   fSKNuMuBarFluxNominal = new TH1D(*fSKNuMuBarFluxTunedExt);
//   fSKNuEFluxNominal     = new TH1D(*fSKNuEFluxTunedExt    );
//   fSKNuEBarFluxNominal  = new TH1D(*fSKNuEBarFluxTunedExt );

//   fSKNuMuFluxNominal    -> SetDirectory(0);
//   fSKNuMuBarFluxNominal -> SetDirectory(0);
//   fSKNuEFluxNominal     -> SetDirectory(0);
//   fSKNuEBarFluxNominal  -> SetDirectory(0);

//   fSKRHCNuMuFluxNominal    = new TH1D(*fSKRHCNuMuFluxTunedExt   );
//   fSKRHCNuMuBarFluxNominal = new TH1D(*fSKRHCNuMuBarFluxTunedExt);
//   fSKRHCNuEFluxNominal     = new TH1D(*fSKRHCNuEFluxTunedExt    );
//   fSKRHCNuEBarFluxNominal  = new TH1D(*fSKRHCNuEBarFluxTunedExt );

//   fSKRHCNuMuFluxNominal    -> SetDirectory(0);
//   fSKRHCNuMuBarFluxNominal -> SetDirectory(0);
//   fSKRHCNuEFluxNominal     -> SetDirectory(0);
//   fSKRHCNuEBarFluxNominal  -> SetDirectory(0);

//   const int n_nu = 4;
//   int nu_pdg[n_nu] = { conv::kPdgNuE, conv::kPdgNuMu, conv::kPdgNuEBar, conv::kPdgNuMuBar };

//   for(int inu = 0; inu < n_nu; inu++) {
//      TH1D * hnom = this->SuperkFluxSpectrum(nu_pdg[inu], conv::kFlxNominal  );
//      for(int ie = 1; ie <= hnom->GetNbinsX(); ie++) {
//         double E = hnom->GetBinCenter (ie);
//         double f = hnom->GetBinContent(ie);
//         double w = this->SuperkFluxExtTuningWght(nu_pdg[inu],E);
//         if(w >=0) { f /= w; }
//         hnom->SetBinContent(ie,f);
//      }//ie
//   }//inu

//   for(int inu = 0; inu < n_nu; inu++) {
//     TH1D * hnomrhc = this->SuperkFluxSpectrum(nu_pdg[inu], conv::kRHCFlxNominal  );
//     for(int ie = 1; ie <= hnomrhc->GetNbinsX(); ie++) {
//       double E = hnomrhc->GetBinCenter (ie);
//       double f = hnomrhc->GetBinContent(ie);
//       double w = this->SuperkRHCFluxExtTuningWght(nu_pdg[inu],E);
//       if(w >=0) { f /= w; }
//       hnomrhc->SetBinContent(ie,f);
//     }//ie
//   }//inu

//   //
//   // Create the fin/tuned flux histograms and initialize them using the ext/tuned flux
//   // Then re-weight them using the weights that correspond to the BANFF best-fit params.
//   //

//   fSKNuMuFluxTunedFinal    = new TH1D(*fSKNuMuFluxTunedExt   );
//   fSKNuMuBarFluxTunedFinal = new TH1D(*fSKNuMuBarFluxTunedExt);
//   fSKNuEFluxTunedFinal     = new TH1D(*fSKNuEFluxTunedExt    );
//   fSKNuEBarFluxTunedFinal  = new TH1D(*fSKNuEBarFluxTunedExt );

//   fSKNuMuFluxTunedFinal    -> SetDirectory(0);
//   fSKNuMuBarFluxTunedFinal -> SetDirectory(0);
//   fSKNuEFluxTunedFinal     -> SetDirectory(0);
//   fSKNuEBarFluxTunedFinal  -> SetDirectory(0);

//   fSKRHCNuMuFluxTunedFinal    = new TH1D(*fSKRHCNuMuFluxTunedExt   );
//   fSKRHCNuMuBarFluxTunedFinal = new TH1D(*fSKRHCNuMuBarFluxTunedExt);
//   fSKRHCNuEFluxTunedFinal     = new TH1D(*fSKRHCNuEFluxTunedExt    );
//   fSKRHCNuEBarFluxTunedFinal  = new TH1D(*fSKRHCNuEBarFluxTunedExt );

//   fSKRHCNuMuFluxTunedFinal    -> SetDirectory(0);
//   fSKRHCNuMuBarFluxTunedFinal -> SetDirectory(0);
//   fSKRHCNuEFluxTunedFinal     -> SetDirectory(0);
//   fSKRHCNuEBarFluxTunedFinal  -> SetDirectory(0);

//   this->ApplyFHCpRHCBANFFWeights(false); // reweight ext/tuned to fin/tuned flux (corresponds to BANFF best-fit)

//   //
//   // Create current flux histograms and initialize them using the fin/tuned flux.
//   //

// #ifdef __show_init_mesg_info__
//   cout << "Exp/Flux >> Booking and initializing FHC+RHC current flux histograms" << endl;
// #endif

//   fSKNuMuFluxCurrent    = new TH1D(*fSKNuMuFluxTunedFinal   );
//   fSKNuMuBarFluxCurrent = new TH1D(*fSKNuMuBarFluxTunedFinal);
//   fSKNuEFluxCurrent     = new TH1D(*fSKNuEFluxTunedFinal    );
//   fSKNuEBarFluxCurrent  = new TH1D(*fSKNuEBarFluxTunedFinal );

//   fSKNuMuFluxCurrent    -> SetDirectory(0);
//   fSKNuMuBarFluxCurrent -> SetDirectory(0);
//   fSKNuEFluxCurrent     -> SetDirectory(0);
//   fSKNuEBarFluxCurrent  -> SetDirectory(0);

//   fSKRHCNuMuFluxCurrent    = new TH1D(*fSKRHCNuMuFluxTunedFinal   );
//   fSKRHCNuMuBarFluxCurrent = new TH1D(*fSKRHCNuMuBarFluxTunedFinal);
//   fSKRHCNuEFluxCurrent     = new TH1D(*fSKRHCNuEFluxTunedFinal    );
//   fSKRHCNuEBarFluxCurrent  = new TH1D(*fSKRHCNuEBarFluxTunedFinal );

//   fSKRHCNuMuFluxCurrent    -> SetDirectory(0);
//   fSKRHCNuMuBarFluxCurrent -> SetDirectory(0);
//   fSKRHCNuEFluxCurrent     -> SetDirectory(0);
//   fSKRHCNuEBarFluxCurrent  -> SetDirectory(0);
// }
// //







//________________________________________________________________________________________
void NuFlux::ApplyBANFFWeights(bool get_current) const
{
// Apply BANFF weights to ext/tuned flux.
//
// Input param     | Action    
// ---------------------------------------------------
// false           | apply BANFF best-fit params to fin/tuned flux
// true            | apply BANFF currrent params to current   flux 

#ifdef __show_flux_twk_mesg_info__
   cout << "Exp/Flux >> Applying BANFF weights to SuperK flux " 
        << ( (get_current) ? "[ext/tuned --> current]" : "[ext/tuned --> fin/tuned]" )
        << endl;   
#endif
#ifdef __show_flux_twk_mesg_debug_level2__
   fBANFF->PrintSystematicVariations();
#endif

   const int np = conv::kNumOfParamsInBANFFfit;
   
  for(int ip = 0; ip < np; ip++) {
     conv::PhysParamId_t p = conv::PhysParamId_t ( conv::kListOfBANFFParams[ip]);
     // flux param?
     bool is_sk_flux_param = conv::kParamInBANFFfitIsSuperkFluxParam[ip];
     if(is_sk_flux_param) 
     {
       
        // Get the current or best-fit BANFF Superk flux param value
        double banff_param_value = (get_current) ? 
             fBANFF->CurrentValue(p) :
             fBANFF->InitialValue(p) ; 

#ifdef __show_flux_twk_mesg_debug__
        cout << "Exp/Flux >> Taking into account the effect of BANFF flux param: "
             << util::AsString(p) << " (value = " << banff_param_value << ")" << endl;
#endif

        // BANFF param -> neutrino PDG code
        int nu_pdg = conv::kNuPdgCodeOfFluxParamInBANFFfit[ip];

        // BANFF param -> true energy range [variable size bin]
        double Emin = conv::kEtrueRangeOfFluxParamInBANFFfit[ip][0];
        double Emax = conv::kEtrueRangeOfFluxParamInBANFFfit[ip][1];
        assert(Emin >= 0. && Emax > 0);

        // Access flux histograms

        // original spectrum:
        TH1D * flux_spectrum_orig = this->SuperkFluxSpectrum(nu_pdg, conv::kFlxExtTuning); 
        // original spectrum + BANFF weights
        TH1D * flux_spectrum_wght = (get_current) ?
           this->SuperkFluxSpectrum(nu_pdg, conv::kFlxCurrent) :
           this->SuperkFluxSpectrum(nu_pdg, conv::kFlxFinalTuning);

        // Note: 
        // The flux_spectrum histograms used in this code, use a binning scheme
        // (typically fixed-size 50-MeV bins, although other fixed or variable binning
        // schemes are possible) which is, in principle, different from the coarser
        // binning scheme used by BANFF.
        // Below, we will refer to the two binning schemes as `own' and `BANFF'

        // Find range of `own' bins which corresponds to the current `BANFF' bin ([Emin,Emax] range)
        int bin_min = flux_spectrum_orig->FindBin(Emin + conv::kASmallNum); 
        int bin_max = flux_spectrum_orig->FindBin(Emax - conv::kASmallNum);

#ifdef __show_flux_twk_mesg_debug__
        cout << "Exp/Flux >> BANFF parameter tweaks Superk flux for neutrino PDG = " << nu_pdg 
             << " in true energy range = [" << Emin << ", " << Emax << "]"
             << endl;
        cout << "Exp/Flux >> In our binning scheme, this energy range corresponds to bins = ["
             << bin_min << ", " << bin_max << "]" 
             << endl;
#endif	
        // Loop over `own' bins contained in the current `BANFF' bin and
        // "distribute" the `BANFF' bin flux variation to our `own' bins.
        for ( int bin = bin_min; bin <= bin_max; bin++ ) {
           // Un-tweaked flux in current bin 
           double flux = flux_spectrum_orig->GetBinContent(bin);
           // Tweaked flux in current bin
           double flux_twk = flux * banff_param_value;
#ifdef __show_flux_twk_mesg_debug_level2__
           cout << "Exp/Flux >> - Our bin = " << bin
                << ", flux = " << flux << " --> " << flux_twk << endl;
#endif
           if(flux_twk < 0) {
                cerr << "Exp/Flux >> Tweaked flux for neutrino PDG = " << nu_pdg 
                     << " in bin = " << bin << " is = " << flux_twk << " !!! " << endl;
                cerr << "Exp/Flux >> Negative value was the result of tweaking BANFF param: " 
                     << util::AsString(p) << " (value = " << banff_param_value << ")" << endl;
                exit(1);
           }
           flux_spectrum_wght->SetBinContent(bin,flux_twk);
        } // bin
     } 
     else 
     {
#ifdef __show_flux_twk_mesg_debug__
        cout << "Exp/Flux >> BANFF param: " << util::AsString(p) 
             << " is not a flux param and it is skipped..." << endl;
#endif
     }//sk flux param?
  }//ip loop  
}
//________________________________________________________________________________________
void NuFlux::ApplyRHCBANFFWeights(bool get_current) const
{
  // Apply RHC BANFF weights to ext/tuned flux.
  //
  // Input param     | Action    
  // ---------------------------------------------------
  // false           | apply RHC BANFF best-fit params to fin/tuned flux
  // true            | apply RHC BANFF currrent params to current   flux 

#ifdef __show_flux_twk_mesg_info__
  cout << "Exp/Flux >> Applying RHC BANFF weights to SuperK flux " 
       << ( (get_current) ? "[ext/tuned --> current]" : "[ext/tuned --> fin/tuned]" )
       << endl;   
#endif
#ifdef __show_flux_twk_mesg_debug_level2__
  fBANFF->PrintRHCSystematicVariations();
#endif

  const int np = conv::kNumOfRHCParamsInBANFFfit;

  for(int ip = 0; ip < np; ip++) {
    conv::PhysParamId_t p = conv::PhysParamId_t ( (int)conv::kRHCNuisBANFF00 + ip );
    // flux param?
    bool is_sk_flux_param = conv::kRHCParamInBANFFfitIsSuperkFluxParam[ip];

    if(is_sk_flux_param) 
      {
        // Get the current or best-fit BANFF Superk flux param value
        double banff_param_value = (get_current) ? 
          fBANFF->CurrentValue(p) :
          fBANFF->InitialValue(p) ;

#ifdef __show_flux_twk_mesg_debug__
        cout << "Exp/Flux >> Taking into account the effect of BANFF flux param: "
             << util::AsString(p) << " (value = " << banff_param_value << ")" << endl;
#endif

        // BANFF param -> neutrino PDG code
        int nu_pdg = conv::kNuPdgCodeOfRHCFluxParamInBANFFfit[ip];

        // BANFF param -> true energy range [variable size bin]
        double Emin = conv::kEtrueRangeOfRHCFluxParamInBANFFfit[ip][0];
        double Emax = conv::kEtrueRangeOfRHCFluxParamInBANFFfit[ip][1];
        assert(Emin >= 0. && Emax > 0);

        // Access flux histograms

        // original spectrum:
        TH1D * flux_spectrum_orig = this->SuperkFluxSpectrum(nu_pdg, conv::kRHCFlxExtTuning); 
        // original spectrum + BANFF weights
        TH1D * flux_spectrum_wght = (get_current) ?
	  this->SuperkFluxSpectrum(nu_pdg, conv::kRHCFlxCurrent) :
	  this->SuperkFluxSpectrum(nu_pdg, conv::kRHCFlxFinalTuning);

        // Note: 
        // The flux_spectrum histograms used in this code, use a binning scheme
        // (typically fixed-size 50-MeV bins, although other fixed or variable binning
        // schemes are possible) which is, in principle, different from the coarser
        // binning scheme used by BANFF.
        // Below, we will refer to the two binning schemes as `own' and `BANFF'

        // Find range of `own' bins which corresponds to the current `BANFF' bin ([Emin,Emax] range)
        int bin_min = flux_spectrum_orig->FindBin(Emin + conv::kASmallNum); 
        int bin_max = flux_spectrum_orig->FindBin(Emax - conv::kASmallNum);

#ifdef __show_flux_twk_mesg_debug__
        cout << "Exp/Flux >> RHC BANFF parameter tweaks Superk flux for neutrino PDG = " << nu_pdg 
             << " in true energy range = [" << Emin << ", " << Emax << "]"
             << endl;
        cout << "Exp/Flux >> In our binning scheme, this energy range corresponds to bins = ["
             << bin_min << ", " << bin_max << "]" 
             << endl;
#endif

        // Loop over `own' bins contained in the current `BANFF' bin and
        // "distribute" the `BANFF' bin flux variation to our `own' bins.
        for ( int bin = bin_min; bin <= bin_max; bin++ ) {
	  // Un-tweaked flux in current bin 
	  double flux = flux_spectrum_orig->GetBinContent(bin);
	  // Tweaked flux in current bin
	  double flux_twk = flux * banff_param_value;
#ifdef __show_flux_twk_mesg_debug_level2__
	  cout << "Exp/Flux >> - Our bin = " << bin
	       << ", RHC flux = " << flux << " --> " << flux_twk << endl;
#endif
	  if(flux_twk < 0) {
	    cerr << "Exp/Flux >> Tweaked RHC flux for neutrino PDG = " << nu_pdg 
		 << " in bin = " << bin << " is = " << flux_twk << " !!! " << endl;
	    cerr << "Exp/Flux >> Negative value was the result of tweaking RHC BANFF param: " 
		 << util::AsString(p) << " (value = " << banff_param_value << ")" << endl;
	    exit(1);
	  }
	  flux_spectrum_wght->SetBinContent(bin,flux_twk);
        } // bin
      } 
    else 
      {
#ifdef __show_flux_twk_mesg_debug__
        cout << "Exp/Flux >> RHC BANFF param: " << util::AsString(p) 
             << " is not a flux param and it is skipped..." << endl;
#endif
      }//sk flux param?
  }//ip loop  
}
//________________________________________________________________________________________
//Todo: change these methods so that they return the correct BANFF fit
//________________________________________________________________________________________
void NuFlux::TellMeAboutTheBANFFfit(BANFF * banff)
{
  fBANFF = banff;
}
//________________________________________________________________________________________
void NuFlux::TellMeAboutTheRHCBANFFfit(BANFF * banff)
{
  fBANFF = banff;
}
//________________________________________________________________________________________


// // NEW FLUX FHC+RHC
// void NuFlux::TellMeAboutTheFHCpRHCBANFFfit(BANFF * banff)
// {
//   fBANFF = banff;
// }
// //


//________________________________________________________________________________________
double NuFlux::SuperkFlux(conv::SuperkSample_t s, int nu_pdg, double Etrue, conv::FluxType_t flxtype) const
{
  
    if(IsRHCSample(s))
    {
        return this->RHCSuperkFlux(nu_pdg, Etrue, flxtype);
    }
    else
        return this->SuperkFlux(nu_pdg, Etrue, flxtype);
}
//________________________________________________________________________________________
double NuFlux::SuperkFlux(int nu_pdg, double Etrue, conv::FluxType_t flxtype) const
{
    // update current flux (if not already) as a result of the latest BANFF param twk
    if( flxtype == conv::kFlxCurrent) {
        this->ApplyBANFFWeights(true);
    }
    
    TH1D * flux_spectrum = this->SuperkFluxSpectrum(nu_pdg,flxtype);
    if(!flux_spectrum) {
        cerr << "Exp/Flux >> No " << util::AsString(flxtype)
        << " flux spectrum for neutrino PDG = " << nu_pdg
        << endl;
    }
    assert(flux_spectrum);
    int bin = flux_spectrum->FindBin(Etrue);
    double flux = flux_spectrum->GetBinContent(bin);
    return flux;
}
//________________________________________________________________________________________

double NuFlux::RHCSuperkFlux( int nu_pdg, double Etrue, conv::FluxType_t flxtype) const
{
  // update current flux (if not already) as a result of the latest BANFF param twk
  if( flxtype == conv::kRHCFlxCurrent) {
    this->ApplyRHCBANFFWeights(true);
  }


  TH1D * flux_spectrum = this->SuperkFluxSpectrum(nu_pdg,flxtype);
  if(!flux_spectrum) {
    cerr << "Exp/Flux >> No RHC " << util::AsString(flxtype)
         << " flux spectrum for neutrino PDG = " << nu_pdg
         << endl;
  }
  assert(flux_spectrum);
  int bin = flux_spectrum->FindBin(Etrue);
  double flux = flux_spectrum->GetBinContent(bin);
  return flux;
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxExtTuningWght(int nu_pdg, double Etrue) const
{
// Weights required to re-weight nominal flux to ext/tuned flux
//
//  if(Etrue > conv::kEtrueMaxForFluxTuneWeights) 
//  {
//    return 1;
//  }

  TH1D * wght_spectrum = this->SuperkFluxExtTuningWghts(nu_pdg);
  assert(wght_spectrum);

  int bin = wght_spectrum->FindBin(Etrue);
  double wght = wght_spectrum->GetBinContent(bin);
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkRHCFluxExtTuningWght(int nu_pdg, double Etrue) const
{
  // Weights required to re-weight nominal flux to ext/tuned flux
  //
  //  if(Etrue > conv::kEtrueMaxForFluxTuneWeights) 
  //  {
  //    return 1;
  //  }

  TH1D * wght_spectrum = this->SuperkRHCFluxExtTuningWghts(nu_pdg);
  assert(wght_spectrum);

  int bin = wght_spectrum->FindBin(Etrue);
  double wght = wght_spectrum->GetBinContent(bin);
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxFinTuningWght(int nu_pdg, double Etrue) const
{
// Weights required to re-weight ext/tuned flux to fin/tuned flux
//
  double flux_ext = this->SuperkFlux(nu_pdg, Etrue, conv::kFlxExtTuning  );
  double flux_fin = this->SuperkFlux(nu_pdg, Etrue, conv::kFlxFinalTuning);

  double wght = 1.0;
  if(flux_ext >= 0) {
    wght = flux_fin / flux_ext;
  }
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkRHCFluxFinTuningWght(int nu_pdg, double Etrue) const
{
  // Weights required to re-weight ext/tuned flux to fin/tuned flux
  //
  double flux_ext = this->RHCSuperkFlux(nu_pdg, Etrue, conv::kRHCFlxExtTuning  );
  double flux_fin = this->RHCSuperkFlux(nu_pdg, Etrue, conv::kRHCFlxFinalTuning);

  double wght = 1.0;
  if(flux_ext >= 0) {
    wght = flux_fin / flux_ext;
  }
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxSystWght(int nu_pdg, double Etrue) const
{
// Weights required to re-weight fin/tuned flux to current (fin/tuned + systematic
// variarions) flux
//
  double flux_twk = this->SuperkFlux(nu_pdg, Etrue, conv::kFlxCurrent    );
  double flux     = this->SuperkFlux(nu_pdg, Etrue, conv::kFlxFinalTuning);

  double wght = 1.0;
  if(flux >= 0) {
    wght = flux_twk / flux;
  }
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkRHCFluxSystWght(int nu_pdg, double Etrue) const
{
  // Weights required to re-weight fin/tuned flux to current (fin/tuned + systematic
  // variarions) flux
  //
  double flux_twk = this->RHCSuperkFlux(nu_pdg, Etrue, conv::kRHCFlxCurrent    );
  double flux     = this->RHCSuperkFlux(nu_pdg, Etrue, conv::kRHCFlxFinalTuning);

  double wght = 1.0;
  if(flux >= 0) {
    wght = flux_twk / flux;
  }
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkFlux2p2hSystWght(int nu_pdg, double Etrue) const
{
// Weights required to re-weight ext/tuned flux to current (fin/tuned + systematic
// variarions) flux
//
  double flux_twk = this->SuperkFlux(nu_pdg, Etrue, conv::kFlxCurrent    );
  double flux     = this->SuperkFlux(nu_pdg, Etrue, conv::kFlxExtTuning  );

  double wght = 1.0;
  if(flux >= 0) {
    wght = flux_twk / flux;
  }
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkRHCFlux2p2hSystWght(int nu_pdg, double Etrue) const
{
  // Weights required to re-weight ext/tuned flux to current (fin/tuned + systematic
  // variarions) flux
  //
  double flux_twk = this->RHCSuperkFlux(nu_pdg, Etrue, conv::kRHCFlxCurrent    );
  double flux     = this->RHCSuperkFlux(nu_pdg, Etrue, conv::kRHCFlxExtTuning  );

  double wght = 1.0;
  if(flux >= 0) {
    wght = flux_twk / flux;
  }
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxExtTuningWght(
  conv::SuperkSample_t s, conv::SuperkMode_t m, double Etrue) const
{
  int flux_nu_pdg = util::SuperkMode2FluxNuPdg(m);
  if(IsRHCSample(s))
    return this->SuperkRHCFluxExtTuningWght(flux_nu_pdg,Etrue);
  else
    return this->SuperkFluxExtTuningWght(flux_nu_pdg,Etrue);
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxFinTuningWght(
  conv::SuperkSample_t s, conv::SuperkMode_t m, double Etrue) const
{
  int flux_nu_pdg = util::SuperkMode2FluxNuPdg(m);
  if(IsRHCSample(s))
    return this->SuperkRHCFluxFinTuningWght(flux_nu_pdg,Etrue);
  else
    return this->SuperkFluxFinTuningWght(flux_nu_pdg,Etrue);
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxSystWght(
  conv::SuperkSample_t s, conv::SuperkMode_t m, double Etrue) const
{
  int flux_nu_pdg = util::SuperkMode2FluxNuPdg(m);
  if(IsRHCSample(s))
    {
      return this->SuperkRHCFluxSystWght(flux_nu_pdg,Etrue);
    }
  else
    {
      return this->SuperkFluxSystWght(flux_nu_pdg,Etrue);
    }
}
//________________________________________________________________________________________
TH1D * NuFlux::SuperkFluxSpectrum(int nu_pdg, conv::FluxType_t flxtype) const
{
  switch(flxtype) {

     // 
     case (conv::kFlxNominal) :
     {
        if      ( nu_pdg == conv::kPdgNuE     ) return fSKNuEFluxNominal;
        else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKNuEBarFluxNominal;
        else if ( nu_pdg == conv::kPdgNuMu    ) return fSKNuMuFluxNominal;
        else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKNuMuBarFluxNominal;
        else
        {
           cerr << "Exp/Flux >> No nominal flux spectrum for particle with PDG = " 
             << nu_pdg << endl;
           exit(1);
        }
        break;
     }

     //
     case (conv::kFlxExtTuning) : 
     {
        if      ( nu_pdg == conv::kPdgNuE     ) return fSKNuEFluxTunedExt;
        else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKNuEBarFluxTunedExt;
        else if ( nu_pdg == conv::kPdgNuMu    ) return fSKNuMuFluxTunedExt;
        else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKNuMuBarFluxTunedExt;
        else
        {
           cerr << "Exp/Flux >> No ext/tuned flux spectrum for particle with PDG = " 
             << nu_pdg << endl;
           exit(1);
        }
        break;
     }

     //
     case (conv::kFlxFinalTuning) :
     {
        if      ( nu_pdg == conv::kPdgNuE     ) return fSKNuEFluxTunedFinal;
        else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKNuEBarFluxTunedFinal;
        else if ( nu_pdg == conv::kPdgNuMu    ) return fSKNuMuFluxTunedFinal;
        else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKNuMuBarFluxTunedFinal;
        else
        {
           cerr << "Exp/Flux >> No fin/tuned flux spectrum for particle with PDG = " 
             << nu_pdg << endl;
           exit(1);
        }
        break;
     }

     //
     case (conv::kFlxCurrent) :
     {
        if      ( nu_pdg == conv::kPdgNuE     ) return fSKNuEFluxCurrent;
        else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKNuEBarFluxCurrent;
        else if ( nu_pdg == conv::kPdgNuMu    ) return fSKNuMuFluxCurrent;
        else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKNuMuBarFluxCurrent;
        else
        {
           cerr << "Exp/Flux >> No current flux spectrum for particle with PDG = " 
             << nu_pdg << endl;
           exit(1);
        }
        break;
     }
  
  //
  case (conv::kRHCFlxNominal) :
    {
      if      ( nu_pdg == conv::kPdgNuE     ) return fSKRHCNuEFluxNominal;
      else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKRHCNuEBarFluxNominal;
      else if ( nu_pdg == conv::kPdgNuMu    ) return fSKRHCNuMuFluxNominal;
      else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKRHCNuMuBarFluxNominal;
      else
        {
	  cerr << "Exp/Flux >> No RHC nominal flux spectrum for particle with PDG = " 
	       << nu_pdg << endl;
	  exit(1);
        }
      break;
    }

    //
  case (conv::kRHCFlxExtTuning) : 
    {
      if      ( nu_pdg == conv::kPdgNuE     ) return fSKRHCNuEFluxTunedExt;
      else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKRHCNuEBarFluxTunedExt;
      else if ( nu_pdg == conv::kPdgNuMu    ) return fSKRHCNuMuFluxTunedExt;
      else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKRHCNuMuBarFluxTunedExt;
      else
        {
	  cerr << "Exp/Flux >> No RHC ext/tuned flux spectrum for particle with PDG = " 
	       << nu_pdg << endl;
	  exit(1);
        }
      break;
    }

    //
  case (conv::kRHCFlxFinalTuning) :
    {
      if      ( nu_pdg == conv::kPdgNuE     ) return fSKRHCNuEFluxTunedFinal;
      else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKRHCNuEBarFluxTunedFinal;
      else if ( nu_pdg == conv::kPdgNuMu    ) return fSKRHCNuMuFluxTunedFinal;
      else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKRHCNuMuBarFluxTunedFinal;
      else
        {
	  cerr << "Exp/Flux >> No RHC fin/tuned flux spectrum for particle with PDG = " 
	       << nu_pdg << endl;
	  exit(1);
        }
      break;
    }

    //
  case (conv::kRHCFlxCurrent) :
    {
      if      ( nu_pdg == conv::kPdgNuE     ) return fSKRHCNuEFluxCurrent;
      else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKRHCNuEBarFluxCurrent;
      else if ( nu_pdg == conv::kPdgNuMu    ) return fSKRHCNuMuFluxCurrent;
      else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKRHCNuMuBarFluxCurrent;
      else
        {
	  cerr << "Exp/Flux >> No RHC current flux spectrum for particle with PDG = "
	       << nu_pdg << endl;
	  exit(1);
        }
      break;
    }

     //
     default:
        return 0;
  }

  return 0;
}
//________________________________________________________________________________________
TH1D * NuFlux::SuperkFluxExtTuningWghts(int nu_pdg) const
{
  if      ( nu_pdg == conv::kPdgNuE     ) return fSKNuEFluxExtTuningWght;
  else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKNuEBarFluxExtTuningWght;
  else if ( nu_pdg == conv::kPdgNuMu    ) return fSKNuMuFluxExtTuningWght;
  else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKNuMuBarFluxExtTuningWght;
  else
    {
      cerr << "Exp/Flux >> Unknown neutrino PDG value: " << nu_pdg << endl;
      exit(1);
    }
  return 0;
}
//________________________________________________________________________________________
TH1D * NuFlux::SuperkRHCFluxExtTuningWghts(int nu_pdg) const
{
  if      ( nu_pdg == conv::kPdgNuE     ) return fSKRHCNuEFluxExtTuningWght;
  else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKRHCNuEBarFluxExtTuningWght;
  else if ( nu_pdg == conv::kPdgNuMu    ) return fSKRHCNuMuFluxExtTuningWght;
  else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKRHCNuMuBarFluxExtTuningWght;
  else
    {
      cerr << "Exp/Flux >> Unknown neutrino PDG value: " << nu_pdg << endl;
      exit(1);
    }
  return 0;
}
//________________________________________________________________________________________
void NuFlux::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  if(opt.Contains("Flux/SK/NuMu/TuneWghts")) { 
    TH1D * h = new TH1D(*fSKNuMuFluxExtTuningWght);
    h->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK numu flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMuBar/TuneWghts")) { 
    TH1D * h = new TH1D(*fSKNuMuBarFluxExtTuningWght);
    h->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK numubar flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuE/TuneWghts")) { 
    TH1D * h = new TH1D(*fSKNuEFluxExtTuningWght);
    h->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK nue flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuEBar/TuneWghts")) { 
    TH1D * h = new TH1D(*fSKNuEBarFluxExtTuningWght);
    h->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK nuebar flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMu")) { 
    TH1D * hnom  = new TH1D(*fSKNuMuFluxNominal);
    TH1D * htune = new TH1D(*fSKNuMuFluxTunedExt);
    TH1D * htunf = new TH1D(*fSKNuMuFluxTunedFinal);
    hnom ->SetDirectory(0);
    htune->SetDirectory(0);
    htunf->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK numu flux (#neutrinos /cm2 /x POT)");
    hnom->Draw();
    htune->Draw("same");
    htunf->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMuBar")) { 
    TH1D * hnom  = new TH1D(*fSKNuMuBarFluxNominal);
    TH1D * htune = new TH1D(*fSKNuMuBarFluxTunedExt);
    TH1D * htunf = new TH1D(*fSKNuMuBarFluxTunedFinal);
    hnom ->SetDirectory(0);
    htune->SetDirectory(0);
    htunf->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK numubar flux (#neutrinos /cm2 /x POT)");
    hnom->Draw();
    htune->Draw("same");
    htunf->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuE")) { 
    TH1D * hnom  = new TH1D(*fSKNuEFluxNominal);
    TH1D * htune = new TH1D(*fSKNuEFluxTunedExt);
    TH1D * htunf = new TH1D(*fSKNuEFluxTunedFinal);
    hnom ->SetDirectory(0);
    htune->SetDirectory(0);
    htunf->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK nue flux (#neutrinos /cm2 /x POT)");
    hnom->Draw();
    htune->Draw("same");
    htunf->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuEBar")) { 
    TH1D * hnom  = new TH1D(*fSKNuEBarFluxNominal);
    TH1D * htune = new TH1D(*fSKNuEBarFluxTunedExt);
    TH1D * htunf = new TH1D(*fSKNuEBarFluxTunedFinal);
    hnom->SetDirectory(0);
    htune->SetDirectory(0);
    htunf->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK nuebar flux (#neutrinos /cm2 /x POT)");
    hnom->Draw();
    htune->Draw("same");
    htunf->Draw("same");
    c->Update();
  }

}
//________________________________________________________________________________________
void NuFlux::DrawRHC(TVirtualPad * c, Option_t * option) const
  {
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  if(opt.Contains("Flux/SK/NuMuRHC/TuneWghts")) { 
    TH1D * hrhc = new TH1D(*fSKRHCNuMuFluxExtTuningWght);
    hrhc->SetDirectory(0);
    hrhc->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hrhc->GetYaxis()->SetTitle("SuperK RHC numu flux weight");
    hrhc->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMuBarRHC/TuneWghts")) { 
    TH1D * hrhc = new TH1D(*fSKRHCNuMuBarFluxExtTuningWght);
    hrhc->SetDirectory(0);
    hrhc->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hrhc->GetYaxis()->SetTitle("SuperK RHC numubar flux weight");
    hrhc->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuERHC/TuneWghts")) { 
    TH1D * hrhc = new TH1D(*fSKRHCNuEFluxExtTuningWght);
    hrhc->SetDirectory(0);
    hrhc->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hrhc->GetYaxis()->SetTitle("SuperK RHC nue flux weight");
    hrhc->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuEBarRHC/TuneWghts")) { 
    TH1D * hrhc = new TH1D(*fSKRHCNuEBarFluxExtTuningWght);
    hrhc->SetDirectory(0);
    hrhc->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hrhc->GetYaxis()->SetTitle("SuperK RHC nuebar flux weight");
    hrhc->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMuRHC")) { 
    TH1D * hnomrhc  = new TH1D(*fSKRHCNuMuFluxNominal);
    TH1D * htunerhc = new TH1D(*fSKRHCNuMuFluxTunedExt);
    TH1D * htunfrhc = new TH1D(*fSKRHCNuMuFluxTunedFinal);
    hnomrhc ->SetDirectory(0);
    htunerhc->SetDirectory(0);
    htunfrhc->SetDirectory(0);
    hnomrhc->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnomrhc->GetYaxis()->SetTitle("SuperK RHC numu flux (#antineutrinos /cm2 /x POT)");
    hnomrhc->Draw();
    htunerhc->Draw("same");
    htunfrhc->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMuBarRHC")) { 
    TH1D * hnomrhc  = new TH1D(*fSKRHCNuMuBarFluxNominal);
    TH1D * htunerhc = new TH1D(*fSKRHCNuMuBarFluxTunedExt);
    TH1D * htunfrhc = new TH1D(*fSKRHCNuMuBarFluxTunedFinal);
    hnomrhc ->SetDirectory(0);
    htunerhc->SetDirectory(0);
    htunfrhc->SetDirectory(0);
    hnomrhc->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnomrhc->GetYaxis()->SetTitle("SuperK RHC numubar flux (#antineutrinos /cm2 /x POT)");
    hnomrhc->Draw();
    htunerhc->Draw("same");
    htunfrhc->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuERHC")) { 
    TH1D * hnomrhc  = new TH1D(*fSKRHCNuEFluxNominal);
    TH1D * htunerhc = new TH1D(*fSKRHCNuEFluxTunedExt);
    TH1D * htunfrhc = new TH1D(*fSKRHCNuEFluxTunedFinal);
    hnomrhc ->SetDirectory(0);
    htunerhc->SetDirectory(0);
    htunfrhc->SetDirectory(0);
    hnomrhc->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnomrhc->GetYaxis()->SetTitle("SuperK RHC nue flux (#antineutrinos /cm2 /x POT)");
    hnomrhc->Draw();
    htunerhc->Draw("same");
    htunfrhc->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuEBarRHC")) { 
    TH1D * hnomrhc  = new TH1D(*fSKRHCNuEBarFluxNominal);
    TH1D * htunerhc = new TH1D(*fSKRHCNuEBarFluxTunedExt);
    TH1D * htunfrhc = new TH1D(*fSKRHCNuEBarFluxTunedFinal);
    hnomrhc->SetDirectory(0);
    htunerhc->SetDirectory(0);
    htunfrhc->SetDirectory(0);
    hnomrhc->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnomrhc->GetYaxis()->SetTitle("SuperK RHC nuebar flux (#antineutrinos /cm2 /x POT)");
    hnomrhc->Draw();
    htunerhc->Draw("same");
    htunfrhc->Draw("same");
    c->Update();
  }

}
//________________________________________________________________________________________
