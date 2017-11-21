//________________________________________________________________________________________
/*
 
 \copyright (c) 2010-2013, VaLOR group.
 
 30-12-2013  Nick
 Added samples, modes and systematics for reverse horn current running
 
 */
//________________________________________________________________________________________

#define __show_init_mesg_info__
//#define __show_init_mesg_debug__
//#define __show_init_mesg_debug2__
//#define __show_xsecsystconv_check_mesg_info__
//#define __show_xsecsystconv_check_mesg_debug__
//#define __show_banff_tuning_wght_mesg_info__
//#define __show_banff_tuning_wght_mesg_debug__
//#define __show_syst_var_wght_mesg_debug__
//#define __show_build_xsbinning_template_mesg_info__
//#define __show_build_xsbinning_template_mesg_debug__
#define __show_wght_spline_idx_mesg_debug__
//#define __show_isanevent_mesg_debug__
//#define __show_find_wght_spline_mesg_debug__
//#define __save_ori_templates_init_debug__
#define __skip_unexpected_spline_test_in_readsplines__ // skip looking for splines we don't expect (due to our empty bins) - speedy!
//#define __check_splines_not_flat_on_read__ // Make sure the input splines aren't flat - something has gone wrong. needs __use_flat_splines_for_missing_splines__
//#define check_spline_limits       //slow, but should be used when testing
//#define __use_1Rsplines_for_MultiR_sample__

//#define __use_flat_splines_for_missing_splines__ // DON'T USE! will use flat splines for missing splines. Only used between SK MC versions (new SKMC, old splines)

#define WARN_ON_FLAT_SPLINE

//
#define DISABLE_SCC // we don't have the SCC splines yet (waiting on the NIWG) so we'll turn it off.
#define __ALLOW_SPLINES_OUTSIDE_LIMIT__  // allow splines to go outside the min and max limit where defined. To use only for syst effect cross checks

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <TMath.h>

#include "valor/Utilities/UtlSystem.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Utilities/VThreadPool.h"

#include "NuXSec.h"
#include "BANFF.h"
#include "Utilities.h"
#include "Build.h"

using std::ostringstream;
using std::setprecision;

using namespace valor;
using namespace valor::conv;
using namespace valor::util;

ClassImp(NuXSec)

//________________________________________________________________________________________
NuXSec::NuXSec()
{
    this->Init();

    //check that XSecSystConv.h is setup correctly
    CheckListOfXSecSyst(kNumOfXSecSyst,    kListOfXSecSyst,    "kListOfXSecSyst",    false);
    CheckListOfXSecSyst(kNumOfRHCXSecSyst, kListOfRHCXSecSyst, "kListOfRHCXSecSyst", false);
    
    CheckIsSplineParam(kNumOfXSecSyst, kListOfXSecSyst, kXSecSystHasWghtSplines,    "kXSecSystHasWghtSplines");
    CheckIsSplineParam(kNumOfXSecSyst, kListOfXSecSyst, kRHCXSecSystHasWghtSplines, "kRHCXSecSystHasWghtSplines");
    
    CheckIsBANFFParam(kNumOfXSecSyst, kNumOfParamsInBANFFfit - kNumOfFluxParamsInBANFFfit,
                      kListOfXSecSyst, kXSecSystIncludedInBANFF, "kXSecSystIncludedInBANFF");
    CheckIsBANFFParam(kNumOfRHCXSecSyst, kNumOfParamsInBANFFfit - kNumOfFluxParamsInBANFFfit,
                      kListOfRHCXSecSyst, kRHCXSecSystIncludedInBANFF, "kRHCXSecSystIncludedInBANFF");

    CheckERange(kNumOfXSecSyst,    kListOfXSecSyst,    (const double *)kXSecSystErecoCut,    "kXSecSystErecoCut");
    CheckERange(kNumOfRHCXSecSyst, kListOfRHCXSecSyst, (const double *)kRHCXSecSystErecoCut, "kRHCXSecSystErecoCut");
    
    CheckERange(kNumOfXSecSyst,    kListOfXSecSyst,    (const double *)kXSecSystEtrueCut,    "kXSecSystEtrueCut");
    CheckERange(kNumOfRHCXSecSyst, kListOfRHCXSecSyst, (const double *)kRHCXSecSystEtrueCut, "kRHCXSecSystEtrueCut");
    
    CheckSKMode(kNumOfXSecSyst, kListOfXSecSyst, (const bool *)kXsApplyXSecSystToSKMode, "kXsApplyXSecSystToSKMode");
#ifdef __use_1Rsplines_for_MultiR_sample__
    LOG("NuXSec",pWARN) << "Changing spline name to read 1Re splines for the MultiR samples";
#endif

}
//________________________________________________________________________________________
NuXSec::~NuXSec()
{
    this->CleanUp();
}
//________________________________________________________________________________________
void NuXSec::CheckListOfXSecSyst(int nparams, const conv::PhysParamId_t * params, const char * vname, bool is_fatal)
{
  int ix;
  PhysParamId_t curr, prev;

  //check that XSec systematics are monotonically increasing
  for(ix = 1; ix < nparams; ix++) {
    prev = params[ix-1];
    curr = params[ix];
#ifdef __show_xsecsystconv_check_mesg_debug__
    if(ix == 1) {
      LOG("NuXSec",pINFO) << vname << "[" << ix-1 << "] = " << prev << " " << util::AsString(prev);
    }
    LOG("NuXSec",pINFO)   << vname << "[" << ix   << "] = " << curr << " " << util::AsString(curr);
#endif
    if((int)curr != ((int)prev + 1)) {
      LOG("NuXSec",pDEBUG) << "NuXSec/Init >> " << vname << "[" << ix   << "] = " << curr << " " << util::AsString(curr)
           << " comes before " << vname << "[" << ix-1 << "] = " << prev << " " << util::AsString(prev);
      if(is_fatal) {
        util::Die(-1);
      }
    }
  }//ix
}
//________________________________________________________________________________________
void NuXSec::CheckIsSplineParam(int nparams, const conv::PhysParamId_t * params, const bool * spline, const char * vname)
{
  int ix;
  PhysParamId_t curr;
  //Check the spline parameters all don't have "norm" in the AsString (and vice versa)
  for(ix = 0; ix < nparams; ix++) {
    curr     = params[ix];
    bool isspline = spline[ix];
    string name = util::AsString(curr);
#ifdef __show_xsecsystconv_check_mesg_debug__
    LOG("NuXSec",pINFO) << vname << "[" << ix << "] = " << isspline << " " << name;
#endif
    bool expected_spline = (name.find("Norm") == string::npos && name.find("Tuning") == string::npos);
    if(isspline != expected_spline) {
      LOG("NuXSec",pDEBUG) << "NuXSec/Init >> " << vname << (isspline ? " is" : " isn't") << " a spline. It "
           << (expected_spline ? "should" : "shouldn't") << " be " << util::AsString(curr);
      //~ util::Die(-1);
    }
  }//ix
}
//________________________________________________________________________________________
void NuXSec::CheckIsBANFFParam(int nparams, int nbanffxsecparams, const conv::PhysParamId_t * params, const bool * banff, const char * vname)
{
  int ix, nbanff = 0;
  bool isbanff;
  //Check the total number of BANFF XSec parameters are what's expected
  for(ix = 0; ix < nparams; ix++) {
    isbanff = banff[ix];
#ifdef __show_xsecsystconv_check_mesg_debug__
    PhysParamId_t curr    = params[ix];
     LOG("NuXSec",pINFO) << vname << "[" << ix << "] = " << isbanff << " " << util::AsString(curr);
#endif
    if(isbanff)
      nbanff++;
  }//ix
  if(nbanff != nbanffxsecparams) {
    LOG("NuXSec",pDEBUG) << vname << " has " << nbanff
         << " parameters. Expected " << nbanffxsecparams;
    util::Die(-1);
  }
}
//________________________________________________________________________________________
void NuXSec::CheckERange(int nparams, const PhysParamId_t * params, const double * erange, const char * vname)
{
  PhysParamId_t curr, prev;
  //check that XSec parameters effect the correct energy range
  for(int ix = 1; ix < nparams; ix++) {
    prev = params[ix-1];
    curr = params[ix];
    double this_lo = erange[ix*2];
    double this_hi = erange[ix*2 + 1];
    double prev_lo = erange[(ix-1)*2];
    double prev_hi = erange[(ix-1)*2 + 1];
#ifdef __show_xsecsystconv_check_mesg_debug__
    if(ix == 1)
       LOG("NuXSec",pINFO) << vname << "[" << ix-1 << "] = " << prev_lo << "," << prev_hi << " " << util::AsString(prev);
     LOG("NuXSec",pINFO) << vname << "[" << ix << "] = " << this_lo << "," << this_hi << " " << util::AsString(curr);
#endif
    if(TMath::Abs(this_lo - prev_hi) > 1E-6) {
      if(TMath::Abs(prev_hi - 100) < 1E-6 && TMath::Abs(this_lo) < 1E-6) {} //it's ok!
      else {
        LOG("NuXSec",pDEBUG) << "Wrong energy range\n"
             << vname << "[" << ix-1 << "] = " << prev_lo << "," << prev_hi << " " << util::AsString(prev) << "\n"
             << vname << "[" << ix << "] = " << this_lo << "," << this_hi << " " << util::AsString(curr) << "\n"
             << "Should be previous high and this low should be identical (or 100 and 0)";
        util::Die(-1);
      }
    }
    if(TMath::Abs(this_lo + 1) < 1E-6 || TMath::Abs(this_hi + 1) < 1E-6) {
      LOG("NuXSec",pDEBUG) << "Wrong energy range\n" 
           << vname << "[" << ix << "] = " << this_lo << "," << this_hi << " " << util::AsString(curr) << "\n"
           << "XSec parameters shouldn't have -1 here";
    }
  }//ix
}
//________________________________________________________________________________________
void NuXSec::CheckSKMode(int nparams, const conv::PhysParamId_t * params, const bool * modes, const char * vname)
{
  PhysParamId_t curr;
  //check that XSec parameters effect the correct energy range
  for(int ix = 0; ix < nparams; ix++) {
    curr = params[ix];
#if defined(__show_xsecsystconv_check_mesg_debug__) || defined(__show_xsecsystconv_check_mesg_info__)
     LOG("NuXSec",pINFO) << vname << " Parameter " << ix << " " << util::AsString(curr);
#endif
    vector<SuperkMode_t> first_modes_set;
    bool found_good_sample = false;
    SuperkSample_t first_sample = conv::kSKsUndefined;
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      if(s == kSKsUndefined) continue;
#ifdef __show_xsecsystconv_check_mesg_debug__
      LOG("NuXSec",pDEBUG) << util::AsString(s);
#endif
      vector<SuperkMode_t> modes_set;
      for(int im = 0; im < knm; im++) {
        SuperkMode_t m = kSKReactionModes[is][im];
        if(m == kSKmUndefined) continue;
        bool applies_to_mode = modes[is*knm*nparams + im*nparams + ix];
        if(applies_to_mode) {
          modes_set.push_back(m);
#ifdef __show_xsecsystconv_check_mesg_debug__
          LOG("NuXSec",pDEBUG) << util::AsString(m);
#endif
        }
      }//im
      if(!found_good_sample) {
        found_good_sample = true;
        first_modes_set = modes_set;
        first_sample = s;
      }
      else {
        //check that they're equal
        if(first_modes_set != modes_set) {
          LOG("NuXSec",pDEBUG) << vname << " Parameter " << ix << " " << util::AsString(curr)
               << "\tDifference between samples " <<util::AsString(s)<<" and good sample "<<util::AsString(first_sample);
        }
      }
    }//is
    if(!first_modes_set.size()) {
      LOG("NuXSec",pDEBUG) << vname << " Parameter " << ix << " " << util::AsString(curr)
           << "\tHas no effect on any mode!";
#if defined(DISABLE_SCC)
      //      if (!util::IsSCC(curr)) {util::Die(-1);}
#else
      //      util::Die(-1);
#endif
    }
    else {
      //check that they apply to all e.g. CCQE
      // this must be done by printing to screen & analyser
#ifdef __show_xsecsystconv_check_mesg_info__
      PrintAppliesType(util::AppliesToCCQE    (first_modes_set), "CCQE"     );
      PrintAppliesType(util::AppliesToCC2p2h  (first_modes_set), "CC2p2h"   );
      PrintAppliesType(util::AppliesToCCcoh   (first_modes_set), "CCcoh"    );
      PrintAppliesType(util::AppliesToCC1pi   (first_modes_set), "CC1pi"    );
      PrintAppliesType(util::AppliesToCCoth   (first_modes_set), "CCoth"    );
      PrintAppliesType(util::AppliesToNC1pi0  (first_modes_set), "NC1pi0"   );
      PrintAppliesType(util::AppliesToNC1piC  (first_modes_set), "NC1piC"   );
      PrintAppliesType(util::AppliesToNCcoh   (first_modes_set), "NCcoh"    );
      PrintAppliesType(util::AppliesToNC1Gamma(first_modes_set), "NC1Gamma" );
      PrintAppliesType(util::AppliesToNCoth   (first_modes_set), "NCoth"    );
#endif
    }
  }//ix
}
//________________________________________________________________________________________
void NuXSec::Init(void)
{
    fBANFF = 0;
   
    for(int ispl = 0; ispl < knwspl1Rmu; ispl++)    { fXSecWghtSpl1Rmu[ispl]    = 0; }
    for(int ispl = 0; ispl < knwspl1Re; ispl++)     { fXSecWghtSpl1Re[ispl]     = 0; }
    for(int ispl = 0; ispl < kRHCnwspl1Rmu; ispl++) { fRHCXSecWghtSpl1Rmu[ispl] = 0; }
    for(int ispl = 0; ispl < kRHCnwspl1Re; ispl++)  { fRHCXSecWghtSpl1Re[ispl]  = 0; }
    
    /// NEW MULTIRING
    for(int ispl = 0; ispl < knwsplMultiRmu; ispl++) { fXSecWghtSplMultiRmu[ispl] = 0; }
    for(int ispl = 0; ispl < knwsplMultiRe; ispl++)  { fXSecWghtSplMultiRe[ispl]  = 0; }
    ////
    
    for(int is=0; is<kns; is++) {
        for(int im=0; im<knm; im++) {
            fXsBinningTempl   [is][im] = 0;
            fXsBinningTempl2D [is][im] = 0;
        }
    }
    ////
    
    fRPATwkDial            = 1.; // relRPA is the default parameter
    fFermiMomTwkDial       = 0;
    fSpecFuncTwkDial       = 0;
    fWshapeTwkDial         = 0;
    fShapeCCothTwkDial     = 0;
    fPiLessDeltaTwkDial    = 0;
    fBindETwkDial          = 0;
    fNormCCcohTwkDial      = 0;
    fNormNC1piCTwkDial     = 0;
    fNormNCcohTwkDial      = 0;
    fNormNCothTwkDial      = 0;
    fNormNuEToNuMuTwkDial  = 0;
    fNormAntiNuToNuTwkDial = 0;
    
    fRHCRPATwkDial            = 1; // NEW XSEC //relRPA is the default parameter
    fRHCFermiMomTwkDial       = 0;
    fRHCSpecFuncTwkDial       = 0;
    fRHCWshapeTwkDial         = 0;
    fRHCShapeCCothTwkDial     = 0;
    fRHCPiLessDeltaTwkDial    = 0;
    fRHCBindETwkDial          = 0;
    fRHCNormCCcohTwkDial      = 0;
    fRHCNormNC1piCTwkDial     = 0;
    fRHCNormNCcohTwkDial      = 0;
    fRHCNormNCothTwkDial      = 0;
    fRHCNormNuEToNuMuTwkDial  = 0;
    fRHCNormAntiNuToNuTwkDial = 0;
    fNormRHCToFHCTwkDial = 0;
    
    fUntunedCentralValues = 0;
    fRHCUntunedCentralValues = 0;
    
    arr_index = util::ArrIdx::Instance();
    
}
//________________________________________________________________________________________
void NuXSec::CleanUp(void)
{
    for(int ispl = 0; ispl < knwspl1Rmu; ispl++) {
        if (fXSecWghtSpl1Rmu[ispl]) delete fXSecWghtSpl1Rmu[ispl];
    }
    for(int ispl = 0; ispl < knwspl1Re; ispl++) {
        if (fXSecWghtSpl1Re[ispl]) delete fXSecWghtSpl1Re[ispl];
    }
    for(int ispl = 0; ispl < kRHCnwspl1Rmu; ispl++) {
        if (fRHCXSecWghtSpl1Rmu[ispl]) delete fRHCXSecWghtSpl1Rmu[ispl];
    }
    for(int ispl = 0; ispl < kRHCnwspl1Re; ispl++) {
        if (fRHCXSecWghtSpl1Re[ispl]) delete fRHCXSecWghtSpl1Re[ispl];
    }
    
    /// NEW MULTIRING
    for(int ispl = 0; ispl < knwsplMultiRmu; ispl++) {
        if (fXSecWghtSplMultiRmu[ispl]) delete fXSecWghtSplMultiRmu[ispl];
    }
    for(int ispl = 0; ispl < knwsplMultiRe; ispl++) {
        if (fXSecWghtSplMultiRe[ispl]) delete fXSecWghtSplMultiRe[ispl];
    }
    ////
    
    for(int is=0; is<kns; is++) {
        for(int im=0; im<knm; im++) {
            if (fXsBinningTempl   [is][im]) delete fXsBinningTempl   [is][im];
            if (fXsBinningTempl2D [is][im]) delete fXsBinningTempl2D [is][im];

        }
    }
    
    if (fUntunedCentralValues)    delete fUntunedCentralValues;
    if (fRHCUntunedCentralValues) delete fRHCUntunedCentralValues;
}
//________________________________________________________________________________________
bool NuXSec::ReadSplines(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  LOG("NuXSec", pINFO) << "Exp/NuXSec >> Reading neutrino generator inputs (splines only) from input XML file.";
#endif
    
    //
    // Load cross-section weight splines
    //
    
#if SPLINE_TYPE == -1
    LOG("NuXSec",pDEBUG) << "Special mode turned on in which splines are not read, and weights involving splines are set to 1";
    return true;
#endif
    
    std::vector<TFile*> spline_files(conv::kns);
    
    spline_files[kSKsFCFV1RingMuLike    ] = util::XML2TFile(xml,root_node,"xsec/weight_splines/one_mu_ring");
    spline_files[kSKsFCFV1RingELike     ] = util::XML2TFile(xml,root_node,"xsec/weight_splines/one_e_ring" );
    
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )
    spline_files[kSKsRHCFCFV1RingMuLike ] = util::XML2TFile(xml,root_node,"xsec/weight_splines/one_mu_ring_rhc");
    spline_files[kSKsRHCFCFV1RingELike  ] = util::XML2TFile(xml,root_node,"xsec/weight_splines/one_e_ring_rhc" );
#endif
    
    /// NEW MULTIRING (multiring-only splines)
    //~ spline_files[kSKsFCFVMultiRingMuLike]= util::XML2TFile(xml,root_node,"xsec/weight_splines/one_mu_ring");
    //~ spline_files[kSKsFCFVMultiRingELike ]= util::XML2TFile(xml,root_node,"xsec/weight_splines/one_e_ring");
    spline_files[kSKsFCFVMultiRingMuLike] = util::XML2TFile(xml,root_node,"xsec/weight_splines/multi_mu_ring");
    spline_files[kSKsFCFVMultiRingELike ] = util::XML2TFile(xml,root_node,"xsec/weight_splines/multi_e_ring");
    
    assert(spline_files[kSKsFCFV1RingMuLike    ]);
    assert(spline_files[kSKsFCFV1RingELike     ]);
    assert(spline_files[kSKsRHCFCFV1RingMuLike ]);
    assert(spline_files[kSKsRHCFCFV1RingELike  ]);
    assert(spline_files[kSKsFCFVMultiRingMuLike]);
    assert(spline_files[kSKsFCFVMultiRingELike ]);
    
    // Make a flat spline
#ifdef __use_flat_splines_for_missing_splines__
    double x[5] = {-5,-2,0,+2,+5};
    double y[5] = {1,1,1,1,1};
    int n = 5;
    TSpline3 *flatspline = new TSpline3("Missing spline",x,y,n);
#endif //__use_flat_splines_for_missing_splines__
    
    PhysParamId_t p;
    bool need_splines;
    int nXSecSyst;
    
    // loop over samples
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        
#ifdef __show_init_mesg_info__
        LOG("NuXSec",pINFO) << "Reading splines for " << util::AsString(s);
#endif
        
        bool is_rhc = util::IsRHCSample(s);
        
        if(is_rhc)
            nXSecSyst = kNumOfRHCXSecSyst;
        else
            nXSecSyst = kNumOfXSecSyst;
        
        // loop over modes for current sample
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined)    continue;
            if(util::IsOscNuMu(m))    continue; //using NuMu splines
            if(util::IsOscNuMuBar(m)) continue; //using NuMuBar splines
            
            // NEW XSEC
#ifdef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
            if(util::IsOscNuEBar(m))  continue; //using NuEBar splines
#endif
            
            // loop over cross-section systematic parameters
            for(int ip = 0; ip < nXSecSyst; ip++) {
                if(is_rhc) {
                    p = kListOfRHCXSecSyst[ip];
                    if(p == kPhypUndefined) continue;
                    //do we need weight splines for the current systematic parameter?
                    need_splines = kRHCXSecSystHasWghtSplines[ip];
                    if(!need_splines) continue;
                }
                else {
                    p = kListOfXSecSyst[ip];
                    if(p == kPhypUndefined) continue;
                    // do we need weight splines for the current systematic parameter?
                    need_splines = kXSecSystHasWghtSplines[ip];
                    if(!need_splines) continue;
                }
#if defined(DISABLE_SCC) && BANFF_BINS != 2
                if (p == kNuisBANFFSCCA          || p == kRHCNuisBANFFSCCA       ||
                    p == kNuisBANFFSCCV          || p == kRHCNuisBANFFSCCV       )
                    continue;
#endif
                // is the current systematic parameter going to be applied to current sample/mode?
                //~ bool is_applied = this->ApplySystToSKMode(p,s,m);
                bool is_applied = kXsApplyXSecSystToSKMode[is][im][ip];
                if(!is_applied) continue;
                // loop over the var.-size  Ereco and Etrue bins for which weight splines were pre-computed
                for(int ir = 0; ir < kNXsErecoBins[is]; ir++) {
                    double Ereco     = util::XsecErecoBinEdges(s,ir) + kOneMeV;
                    
                    for(int it = 0; it < kXsNEtrueBins; it++) {
                        double Etrue     = kXsEtrueBinEdges[is][im][it] + kOneMeV;
                        
                        bool bin_is_beyond_end_of_range = (kXsEtrueBinEdges[is][im][it+1] < 0.);
                        if(bin_is_beyond_end_of_range) continue;
                        
                        bool bin_is_not_empty = this->IsThereAnEvent(s, m, Ereco, Etrue);
                        
#ifdef __skip_unexpected_spline_test_in_readsplines__
                        if(!bin_is_not_empty)
                            continue;
#endif
                        std::string splname = this->BuildSplineName(s,m,p,ir,it);
        		
#ifdef __show_init_mesg_debug__
                        LOG("NuXSec",pINFO) << "Loading cross-section weight spline: "
                        << splname << " for Ereco: " << Ereco << ", Etrue: " << Etrue;
#endif
                        
                        // Copy spline from input file
                        TSpline3 * spl = NULL;
                        //Todo: add separate splines for RHC running if necessary
                        spl = (TSpline3 *) spline_files[s]->Get(splname.c_str());
                        
                        if(!spl) {
                            if(!bin_is_not_empty) {
#ifdef __show_init_mesg_debug__
                                LOG("NuXSec",pDEBUG) << "No events or spline for: "
                                << splname;
#endif
                                continue;
                            }
                            else {
#ifndef __use_flat_splines_for_missing_splines__
                              LOG("NuXSec",pERROR) << "Reading from file: "
                                                   <<  spline_files[s]->GetName() << "\n"
                                                   << ">> Null cross-section weight spline: "
                                                   << splname
#if SPLINE_TYPE == 0
                                << ", " << AsString(m)
#endif
                                << " Ereco: " << Ereco << ", Etrue: " << Etrue;
                                LOG("NuXSec",pDEBUG) << AsString(s);
                                util::Die(1);
#else
                                spl = flatspline;
#ifdef WARN_ON_FLAT_SPLINE
                                LOG("NuXSec",pDEBUG) << "Using flat spline: " <<splname
                                << ", " << AsString(m)
                                << " Ereco: " << Ereco << ", Etrue: " << Etrue;
#endif // WARN_ON_FLAT_SPLINE
#endif //__use_flat_splines_for_missing_splines__
                            } //bin_is_not_empty
                        } // !spl
                        else {
                            if(!bin_is_not_empty) {
                                util::Die(1);
                            }
                        }
                        
#ifndef __skip_unexpected_spline_test_in_readsplines__
#ifndef __USE_NIWG2015__
                        string spln(spl->GetTitle());
                        spln.erase(spln.end()-7,spln.end());
                        int nevents_spl = atoi(spln.c_str());
                        int nevents_hist = this->HowManyEvents(s, m, Ereco, Etrue);
                        if(nevents_spl != nevents_hist) {
                            LOG("NuXSec",pDEBUG) << "Spline and histogram disagree on number of events! "
                            << "Spline thinks there are: " << nevents_spl
                            << ", Histogram thinks there are: " << nevents_hist
                            << "; " << splname
                            << "; Ereco: " << Ereco - kOneMeV << " to " << util::XsecErecoBinEdges(s,ir+1)
                            << ", Etrue: " << Etrue - kOneMeV << " to " << kXsEtrueBinEdges[is][im][it+1];
                        }
#endif
#endif
                        // Get index for current weight spline
                        int idx = this->WghtSplineIdx(p,s,m,Ereco,Etrue);
                        if(idx < 0) {
                            LOG("NuXSec",pDEBUG) << "Exp/NuXSec >> WghtSplineIdx returned invalid index for spline: "
                            << splname;
                            util::Die(1);
                        }
                        
#ifdef __show_init_mesg_debug__
                        LOG("NuXSec",pDEBUG) << "Storing cross-section weight spline at slot: "
                        << idx << "  for spline: " << splname;
                        LOG("NuXSec",pDEBUG) << "Spline thinks there are "
                        << spl->GetTitle();
#endif
                        
#ifdef __check_splines_not_flat_on_read__
                        if (spl!=flatspline && bin_is_not_empty) {
                            int n_points_to_check = 11;
                            double zero_point = spl->Eval(0);
                            bool flat=true;
                            for (int ipoint = 0 ; ipoint < n_points_to_check ; ipoint++) {
                                double val = spl->GetXmin()+(spl->GetXmax()-spl->GetXmin()*ipoint/n_points_to_check);
                                if (TMath::Abs(spl->Eval(val)-zero_point) > 1E-3){
                                    flat=false;
                                    break;
                                }
                            }
                            if (flat) {
			      LOG("NuXSec",pERROR) << "Warning, input spline from file is flat: "<<std::setw(40)<<std::left<<splname
                                <<std::fixed
                                << "( Ereco: " << setprecision(2) << Ereco_min << "-" << setprecision(2) << Ereco_max
                                << ", Etrue: " << setprecision(2) << Etrue_min << "-" << setprecision(2) << Etrue_max <<")";
                            }
                        }
#endif
                        
                        if     (s == kSKsFCFV1RingMuLike)    fXSecWghtSpl1Rmu   [idx] = new TSpline3(*spl);
                        else if(s == kSKsFCFV1RingELike)     fXSecWghtSpl1Re    [idx] = new TSpline3(*spl);
                        else if(s == kSKsRHCFCFV1RingMuLike) fRHCXSecWghtSpl1Rmu[idx] = new TSpline3(*spl);
                        else if(s == kSKsRHCFCFV1RingELike)  fRHCXSecWghtSpl1Re [idx] = new TSpline3(*spl);
                        
                        
                        /// NEW MULTIRING (update with splines for MultiRing sample)
                        else if(s == kSKsFCFVMultiRingMuLike) fXSecWghtSplMultiRmu[idx] = new TSpline3(*spl);
                        else if(s == kSKsFCFVMultiRingELike)  fXSecWghtSplMultiRe [idx] = new TSpline3(*spl);
                        ///
                        
#ifdef __show_init_mesg_debug2__
                        //check that the spline is not flat
                        LOG("NuXSec",pDEBUG) << "Spline weights at -1,0,+1 are: " << spl->Eval(-1) << "\t"
                        << spl->Eval(0) << "\t"
                        << spl->Eval(+1);
                        
#endif
                    }//Etrue
                }//Ereco
            }//p
        }//m
    }//s
    
#ifdef __show_init_mesg_info__
    LOG("NuXSec",pDEBUG) << "Done reading neutrino cross-section inputs (splines only)!";
#endif
    delete spline_files[kSKsFCFV1RingMuLike    ];
    delete spline_files[kSKsFCFV1RingELike     ];
    delete spline_files[kSKsRHCFCFV1RingMuLike ];
    delete spline_files[kSKsRHCFCFV1RingELike  ];
    delete spline_files[kSKsFCFVMultiRingMuLike];
    delete spline_files[kSKsFCFVMultiRingELike ];
    
    return true;
}
//________________________________________________________________________________________
bool NuXSec::ReadSplines2D(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  LOG("NuXSec", pINFO) << "Exp/NuXSec >> Reading neutrino generator inputs (splines only) from input XML file.";
#endif
    
    //
    // Load cross-section weight splines
    //
    
#if SPLINE_TYPE == -1
  LOG("NuXSec",pDEBUG) << "Exp/NuXSec >> Special mode turned on in which splines are not read, and weights involving splines are set to 1";
  return true;
#endif
    
  std::vector<TFile*>  spline_files(conv::kns);
    
  spline_files[kSKsFCFV1RingMuLike    ] = util::XML2TFile(xml,root_node,"xsec/weight_splines_2d/one_mu_ring");
  spline_files[kSKsFCFV1RingELike     ] = util::XML2TFile(xml,root_node,"xsec/weight_splines_2d/one_e_ring" );
    
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )
  spline_files[kSKsRHCFCFV1RingMuLike ] = util::XML2TFile(xml,root_node,"xsec/weight_splines_2d/one_mu_ring_rhc");
  spline_files[kSKsRHCFCFV1RingELike  ] = util::XML2TFile(xml,root_node,"xsec/weight_splines_2d/one_e_ring_rhc" );
#endif
    
  /// NEW MULTIRING (multiring-only splines)
  //~ spline_files[kSKsFCFVMultiRingMuLike]= util::XML2TFile(xml,root_node,"xsec/weight_splines_2d/one_mu_ring");
  //~ spline_files[kSKsFCFVMultiRingELike ]= util::XML2TFile(xml,root_node,"xsec/weight_splines_2d/one_e_ring");
  spline_files[kSKsFCFVMultiRingMuLike] = util::XML2TFile(xml,root_node,"xsec/weight_splines_2d/multi_mu_ring");
  spline_files[kSKsFCFVMultiRingELike ] = util::XML2TFile(xml,root_node,"xsec/weight_splines_2d/multi_e_ring");
  
  assert(spline_files[kSKsFCFV1RingMuLike    ]);
  assert(spline_files[kSKsFCFV1RingELike     ]);
  assert(spline_files[kSKsRHCFCFV1RingMuLike ]);
  assert(spline_files[kSKsRHCFCFV1RingELike  ]);
  assert(spline_files[kSKsFCFVMultiRingMuLike]);
  assert(spline_files[kSKsFCFVMultiRingELike ]);
    
    // Make a flat spline
#ifdef __use_flat_splines_for_missing_splines__
  double x[5] = {-5,-2,0,+2,+5};
  double y[5] = {1,1,1,1,1};
  int n = 5;
  TSpline3 *flatspline = new TSpline3("Missing spline",x,y,n);
#endif //__use_flat_splines_for_missing_splines__
  
  PhysParamId_t p;
  bool need_splines;
  int nXSecSyst;

  //  std::atomic_bool splines_status(false);
  
  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    LOG("NuXSec",pINFO) << "Looping for sample: " << util::AsString(s);

    
    //std::function<void(void)> read_splines_2D = [&, s, is]{

#ifdef __show_init_mesg_info__
    LOG("NuXSec",pINFO) << "Reading splines for " << util::AsString(s) << "\n"
                        << "From file: " << spline_files[s]->GetName();
#endif
        
    bool is_rhc = util::IsRHCSample(s);
    
    if(is_rhc)
      nXSecSyst = kNumOfRHCXSecSyst;
    else
      nXSecSyst = kNumOfXSecSyst;
    
    // loop over modes for current sample
    for(int im = 0; im < knm; im++) {
      SuperkMode_t m = kSKReactionModes[is][im];
      if(m == kSKmUndefined)    continue;
      if(util::IsOscNuMu(m))    continue; //using NuMu splines
      if(util::IsOscNuMuBar(m)) continue; //using NuMuBar splines
            
      // NEW XSEC
#ifdef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
      if(util::IsOscNuEBar(m))  continue; //using NuEBar splines
#endif
            
      // loop over cross-section systematic parameters
      for(int ip = 0; ip < nXSecSyst; ip++) {
        if(is_rhc) {
          p = kListOfRHCXSecSyst[ip];
          if(p == kPhypUndefined) continue;
          //do we need weight splines for the current systematic parameter?
          need_splines = kRHCXSecSystHasWghtSplines[ip];
          if(!need_splines) continue;
        }
        else {
          p = kListOfXSecSyst[ip];
          if(p == kPhypUndefined) continue;
          // do we need weight splines for the current systematic parameter?
          need_splines = kXSecSystHasWghtSplines[ip];
          if(!need_splines) continue;
        }
#if defined(DISABLE_SCC) && BANFF_BINS != 2
        if (p == kNuisBANFFSCCA          || p == kRHCNuisBANFFSCCA       ||
            p == kNuisBANFFSCCV          || p == kRHCNuisBANFFSCCV       )
          continue;
#endif
        // is the current systematic parameter going to be applied to current sample/mode?
        //~ bool is_applied = this->ApplySystToSKMode(p,s,m);
        bool is_applied = kXsApplyXSecSystToSKMode[is][im][ip];       
                
        if(!is_applied) continue;

        std::string splname = this->BuildSplineName(s, m, p);
        TObjArray* splineArray = static_cast<TObjArray*>(spline_files[s]->Get(splname.c_str()));

        // loop over the var.-size  Ereco and Etrue bins for which weight splines were pre-computed
        for(int ir = 0; ir < kNXsErecoBins[is]; ir++) {
          double Ereco     = util::XsecErecoBinEdges(s,ir) + kOneMeV;
                    
          for(int it = 0; it < kXsNEtrueBins; it++) {
            double Etrue     = kXsEtrueBinEdges[is][im][it] + kOneMeV;
                        
            for(int iv = 0; iv < kNXsQ2Bins[is]; iv++) {
              double Q2     = util::XsecQ2BinEdges(s,iv) + kOneMeV; // Edit this for diff Q2 bins per mode
                            
              bool bin_is_beyond_end_of_etrue_range = (kXsEtrueBinEdges[is][im][it+1] < 0.);
              bool bin_is_beyond_end_of_q2_range = (util::XsecQ2BinEdges(s,iv+1) < 0.);
                            
              if(bin_is_beyond_end_of_etrue_range || bin_is_beyond_end_of_q2_range) continue;
                            
              bool bin_is_not_empty = this->IsThereAnEvent(s, m, Ereco, Etrue,Q2);
                            
#ifdef __skip_unexpected_spline_test_in_readsplines__
              if(!bin_is_not_empty)
                continue;
#endif
              int inIdx = it * kNXsErecoBins[is] * kNXsQ2Bins[is] + ir * kNXsQ2Bins[is] + iv;
                            
//              std::string splname = this->BuildSplineName(s,m,p,ir,it,iv);
                            
#ifdef __show_init_mesg_debug__
              LOG("NuXSec",pINFO) << "Loading cross-section weight spline: "
                                  << splname << " for Ereco: " << Ereco << ", Etrue: " << Etrue << ", Q2: " << Q2;
#endif
                            
              // Copy spline from input file
              TSpline3 * spl = NULL;
              spl = static_cast<TSpline3*>(splineArray->At(inIdx));

              if(!spl) {
                if(!bin_is_not_empty) {
#ifdef __show_init_mesg_debug__
                  LOG("NuXSec",pDEBUG) << "No events or spline for: "
                       << splname;
#endif
                  continue;
                }
                else {
#ifndef __use_flat_splines_for_missing_splines__
                  LOG("NuXSec",pERROR) 
                    << "Reading from file: "
                    << spline_files[s]->GetName() << "\n"
                    << ">> Null cross-section weight spline: "
                    << splname
#if SPLINE_TYPE == 0
                    << ", " << AsString(m)
#endif
                    << " Ereco: " << Ereco << ", Etrue: " << Etrue << ", Q2: " << Q2;
                  LOG("NuXSec",pDEBUG) << AsString(s);
                  //util::Die(1);
#else
                  spl = flatspline;
#ifdef WARN_ON_FLAT_SPLINE
                  LOG("NuXSec", pWARN) << "Using flat spline: "<< splname
                                << ", " << AsString(m)
                                << " Ereco: " << Ereco << ", Etrue: " << Etrue << ", Q2: " << Q2
                                << " Events in MC: " << this->HowManyEvents(s, m, Ereco, Etrue,Q2);
#endif // WARN_ON_FLAT_SPLINE
#endif //__use_flat_splines_for_missing_splines__
                } //bin_is_not_empty
              } // !spl
              else {
                if(!bin_is_not_empty) {
                  util::Die(1);
                }
              }
                            
#ifndef __skip_unexpected_spline_test_in_readsplines__
#ifndef __USE_NIWG2015__
              string spln(spl->GetTitle());
              spln.erase(spln.end()-7,spln.end());
              int nevents_spl = atoi(spln.c_str());
              int nevents_hist = this->HowManyEvents(s, m, Ereco, Etrue, Q2);
              if(nevents_spl != nevents_hist) {
                LOG("NuXSec",pDEBUG) << "Exp/NuXSec >> Spline and histogram disagree on number of events! "
                     << "Spline thinks there are: " << nevents_spl
                     << ", Histogram thinks there are: " << nevents_hist
                     << "; " << splname
                     << "; Ereco: " << Ereco - kOneMeV << " to " << util::XsecErecoBinEdges(s,ir+1)
                     << ", Etrue: " << Etrue - kOneMeV << " to " << kXsEtrueBinEdges[is][im][it+1]
                     << ", Q2   : " << Q2    - kOneMeV << " to " << util::XsecQ2BinEdges(s,iv+1);
              }
#endif
#endif
              // Get index for current weight spline
              int idx = this->WghtSplineIdx(p,s,m,Ereco,Etrue,Q2);
              if(idx < 0) {
                LOG("NuXSec",pERROR) << "WghtSplineIdx returned invalid index for spline: "
                                     << splname;
                util::Die(1);
              }
                            
#ifdef __show_init_mesg_debug__
              LOG("NuXSec",pDEBUG) << "Storing cross-section weight spline at slot: "
                   << idx << "  for spline: " << splname;
              LOG("NuXSec",pDEBUG) << "Spline thinks there are "
                   << spl->GetTitle();
#endif
                            
#ifdef __check_splines_not_flat_on_read__
              if (spl!=flatspline && bin_is_not_empty) {
                int n_points_to_check = 11;
                double zero_point = spl->Eval(0);
                bool flat=true;
                for (int ipoint = 0 ; ipoint < n_points_to_check ; ipoint++) {
                  double val = spl->GetXmin()+(spl->GetXmax()-spl->GetXmin()*ipoint/n_points_to_check);
                  if (TMath::Abs(spl->Eval(val)-zero_point) > 1E-3){
                    flat=false;
                    break;
                  }
                }
                if (flat) {
                  LOG("NuXSec",pERROR) << "Warning, input spline from file is flat: "<<std::setw(40)<<std::left<<splname
                     <<std::fixed
                     << "( Ereco: " << setprecision(2) << Ereco_min << "-" << setprecision(2) << Ereco_max
                     << ", Etrue: " << setprecision(2) << Etrue_min << "-" << setprecision(2) << Etrue_max <<")"
                     << ", Q2   : " << setprecision(2) << Q2_min << "-" << setprecision(2) << Q2_max <<")";
                }
              }
#endif
                            
              if     (s == kSKsFCFV1RingMuLike)    fXSecWghtSpl1Rmu   [idx] = new TSpline3(*spl);
              else if(s == kSKsFCFV1RingELike)     fXSecWghtSpl1Re    [idx] = new TSpline3(*spl);
              else if(s == kSKsRHCFCFV1RingMuLike) fRHCXSecWghtSpl1Rmu[idx] = new TSpline3(*spl);
              else if(s == kSKsRHCFCFV1RingELike)  fRHCXSecWghtSpl1Re [idx] = new TSpline3(*spl);
                            
              /// NEW MULTIRING (update with splines for MultiRing sample)
              else if(s == kSKsFCFVMultiRingMuLike) fXSecWghtSplMultiRmu[idx] = new TSpline3(*spl);
              else if(s == kSKsFCFVMultiRingELike)  fXSecWghtSplMultiRe [idx] = new TSpline3(*spl);
              ///
                            
#ifdef __show_init_mesg_debug2__
              //check that the spline is not flat
              LOG("NuXSec",pDEBUG) << "Spline: " << spl->GetName()
                   << "\n Spline weights at -2,-1,0,+1,+2 are: " ;
              LOG("NuXSec",pDEBUG) << spl->Eval(-2) << "\t"
                   << spl->Eval(-1) << "\t"
                   << spl->Eval(0)  << "\t"
                   << spl->Eval(+1) << "\t"
                   << spl->Eval(+2);
#endif
            }//Q2
          }//Etrue
        }//Ereco
      }//p
    }//m
    //    };
    //    util::VThreadPool::Instance()->AddTask(read_splines_2D);        
  } //s
//  LOG("NuXSec", pDEBUG) << util::VThreadPool::Instance()->GetNTasksInQueue();
//util::VThreadPool::Instance()->Wait();

#ifdef __show_init_mesg_info__
    LOG("NuXSec",pDEBUG) << "Done reading neutrino cross-section inputs (splines only)!";
#endif

  delete spline_files[kSKsFCFV1RingMuLike    ];
  delete spline_files[kSKsFCFV1RingELike     ];
  delete spline_files[kSKsRHCFCFV1RingMuLike ];
  delete spline_files[kSKsRHCFCFV1RingELike  ];
  delete spline_files[kSKsFCFVMultiRingMuLike];
  delete spline_files[kSKsFCFVMultiRingELike ];

    // return splines_status(true);
  return true;
}

//________________________________________________________________________________________
bool NuXSec::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
    LOG("NuXSec",pDEBUG) << "Reading neutrino generator inputs (minus splines) from input XML file.";
#endif
    
    //
    //
    //
   
    //Old 1D
    
    //book the bool templ histograms
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        const int      n_ereco_bins = kNXsErecoBins[is];
        const double * ereco_bin_edges = util::XsecErecoBinEdges(s);
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            const int      n_etrue_bins    = util::NumVariableXsEtrueBins (s,m);
            const double * etrue_bin_edges = util::VariableXsEtrueBinEdges(s,m);
            fXsBinningTempl[is][im] = new TH2D(
                                               Form("bool_templ_%d_%d",is,im), util::AsString(m).c_str(),
                                               n_ereco_bins, ereco_bin_edges, n_etrue_bins, etrue_bin_edges);
            fXsBinningTempl[is][im]->SetDirectory(0);
            fXsBinningTempl[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
            fXsBinningTempl[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
        }
    }
    
    //New 2D
    //book the bool templ histograms
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        const int      n_ereco_bins = kNXsErecoBins[is];
        const double * ereco_bin_edges = util::XsecErecoBinEdges(s);
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            const int      n_etrue_bins    = util::NumVariableXsEtrueBins (s,m);
            const double * etrue_bin_edges = util::VariableXsEtrueBinEdges(s,m);
            const int      n_Q2_bins       = util::NumXsecQ2Bins(s);
            const double * q2_bin_edges    = util::XsecQ2BinEdges(s);
            fXsBinningTempl2D[is][im] = new TH3D(
                                               Form("bool_templ_%d_%d",is,im), util::AsString(m).c_str(),n_Q2_bins,q2_bin_edges,
                                               n_ereco_bins, ereco_bin_edges, n_etrue_bins, etrue_bin_edges);
            fXsBinningTempl2D[is][im]->SetDirectory(0);
            fXsBinningTempl2D[is][im]->GetYaxis()->SetTitle("E_{reco} (GeV)");
            fXsBinningTempl2D[is][im]->GetZaxis()->SetTitle("E_{true} (GeV)");
            fXsBinningTempl2D[is][im]->GetXaxis()->SetTitle("Q^{2} (GeV)");

        }
    }
    
#if  ( defined (__USE_NIWG2015__) || defined (__USE_NIWG2017__) )
    if (fUntunedCentralValues    ) delete fUntunedCentralValues   ;
    if (fRHCUntunedCentralValues ) delete fRHCUntunedCentralValues;
    
    fUntunedCentralValues    = util::XML2TVectorD(xml,root_node,"xsec/untuned_parameter_central_values");
    fRHCUntunedCentralValues = util::XML2TVectorD(xml,root_node,"xsec/untuned_parameter_central_values");
    
    if (not (fUntunedCentralValues && fRHCUntunedCentralValues) ) {
        LOG("NuXSec",pDEBUG)<<"Failed to read <untuned_parameter_central_values>";
        util::Die(-1);
    }
    if (fUntunedCentralValues->GetNoElements()!=kNumOfParamsInBANFFfit-kNumOfFluxParamsInBANFFfit) {
      LOG("NuXSec",pERROR) << __func__ << " >> Expected "<<kNumOfParamsInBANFFfit-kNumOfFluxParamsInBANFFfit
        << " elements in <untuned_parameter_central_values>, got "<<fUntunedCentralValues->GetNoElements();
        util::Die(-1);
    }
    if (fRHCUntunedCentralValues->GetNoElements()!=kNumOfRHCParamsInBANFFfit-kNumOfRHCFluxParamsInBANFFfit) {
         LOG("NuXSec",pERROR) << __func__ << " >> Expected "<<kNumOfRHCParamsInBANFFfit-kNumOfRHCFluxParamsInBANFFfit
        << " elements in <untuned_parameter_central_values>, got "<<fRHCUntunedCentralValues->GetNoElements();
        util::Die(-1);
    }
#endif
    
    // NEW XSEC --> comment old xsec systematics --> not commented any longer because used
    fNormCCcohRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/cccoh");
    fNormNC1piCRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/nc1piC");
    fNormNCcohRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/nccoh");
    fNormNCothRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/ncoth");
    fNormNuEToNuMuRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/numu_to_nue");
    fNormAntiNuToNuRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/nu_to_antinu");
    fRHCNormCCcohRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/cccoh");
    fRHCNormNC1piCRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/nc1piC");
    fRHCNormNCcohRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/nccoh");
    fRHCNormNCothRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/ncoth");
    fRHCNormNuEToNuMuRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/numu_to_nue");
    fRHCNormAntiNuToNuRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/nu_to_antinu");
    fNormRHCToFHCRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/rhc_to_fhc");
#ifdef __show_init_mesg_info__
    LOG("NuXSec",pDEBUG) << "Done reading neutrino cross-section inputs (not splines)!";
#endif

    
    return true;
}

void NuXSec::PrintXSecSystApplication()
{
    cout << "Systematics:" << std::endl;
    int nXSecSyst = 0;
    for(int is = 0; is < kns; ++is)
    {
       SuperkSample_t s = kSKSamples[is];
       if(s == kSKsUndefined)
       {
          continue;
       }
       cout << util::AsString(s) << std::endl;
       bool is_rhc = util::IsRHCSample(s);
       if(is_rhc)
          nXSecSyst = kNumOfRHCXSecSyst;
       else
          nXSecSyst = kNumOfXSecSyst;
       for(int im = 0; im < knm; ++im)
       {
          SuperkMode_t m = kSKReactionModes[is][im];
          if(m == kSKmUndefined)
          {
             continue;
          }
          cout << "   " << util::AsString(m) << std::endl;
          for(int ip = 0; ip < nXSecSyst; ++ip)
          {
             PhysParamId_t p = kPhypUndefined;
             if(!is_rhc)
             {
                p = kListOfXSecSyst[ip];
             }
             else
             {
                p = kListOfRHCXSecSyst[ip];
             }
             if(p == kPhypUndefined)
             {
                continue;
             }
             cout << "      " << util::AsString(p) << ": " << this->ApplySystToSKMode(p, s, m) << std::endl;
          }
       }
    }
}

//Todo: change these methods so that they return the correct BANFF fit
//________________________________________________________________________________________
void NuXSec::TellMeAboutTheBANFFfit(BANFF * banff)
{
    fBANFF = banff;
}
//________________________________________________________________________________________
void NuXSec::TellMeAboutTheRHCBANFFfit(BANFF * banff)
{
    fBANFF = banff;
}
//________________________________________________________________________________________
double NuXSec::FinTuningWght(SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
    double tot_wght = 1.;
    
    int nXSecSyst;
    PhysParamId_t p;
    bool is_banff_xsec_param;
    bool is_rhc = util::IsRHCSample(s);
    
    if(is_rhc)
        nXSecSyst = kNumOfRHCXSecSyst;
    else
        nXSecSyst = kNumOfXSecSyst;
    
    for(int isyst = 0; isyst < nXSecSyst; isyst++) {
        if(is_rhc) {
            p = kListOfRHCXSecSyst[isyst];
            // "Tuning" involves reweighting to BANFF best-fit values.
            // Consider only the cross-section systematics that were included in the BANFF fit.
            is_banff_xsec_param = kRHCXSecSystIncludedInBANFF[isyst];
        }
        else {
            p = kListOfXSecSyst[isyst];
            is_banff_xsec_param = kXSecSystIncludedInBANFF[isyst];
        }
        if(!is_banff_xsec_param) continue;
        // Is the current systematic parameter applied to current sample/mode end energy bin?
        bool apply = this->ApplySyst(p,s,m,Ereco,Etrue);
        if(!apply) continue;
        // Calculate weight
        double wght = this->FinTuningWght_1p(p,s,m,Ereco,Etrue);
        tot_wght *= wght;
    }
    
    return tot_wght;
}
//________________________________________________________________________________________
double NuXSec::FinTuningWght_1p(PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
#if defined(DISABLE_SCC) && BANFF_BINS != 2
    if (p == kNuisBANFFSCCA          || p == kRHCNuisBANFFSCCA       ||
        p == kNuisBANFFSCCV          || p == kRHCNuisBANFFSCCV       )
        return 1;
#endif
    
    // NEW XSEC
    
    //
    // This method provides weights to reweight MC to the BANFF best-fit parameters.
    // It handles the BANFF cross-section parameters which are correlated with the BANFF
    // flux parameters. In the current incarnation of the BANFF fit, these parameters are:
    //
    // - NIWG2015 BANFF param, 2p2h        normalization
    // - NIWG2015 BANFF param, CA5
    // - NIWG2015 BANFF param, BgRES       normalization
    // - NIWG2015 BANFF param, M_{A} QE
    // - NIWG2015 BANFF param, M_{A} RES
    // - NIWG2015 BANFF param, SCCA 
    // - NIWG2015 BANFF param, SCCV 
    // - NIWG2015 BANFF param, Pf
    // - NIWG2015 BANFF param, CCothShape
    // - NIWG2015 BANFF param, Eb  
    // - NIWG2015 BANFF param, CCcoh      normalization
    // - NIWG2015 BANFF param, NCcoh      normalization
    // - NIWG2015 BANFF param, NCoth      normalization
    // - NIWG2015 BANFF param, NuEToNuMu  normalization
    //
    
    // if( ! this->ApplySyst(p,s,m,Ereco,Etrue) ) return 1.; // check already made in function that calls this one
    TSpline3 * wght_spl = 0;
    
    double wght = 1.;
    double tb;
    
    // RESPONSE FUNCTION PARAMETERS
    
    // NEW XSEC NIWG
    
#if ( defined(__USE_NIWG2015__) || defined (__USE_NIWG2017__) )
    
    //
    // BANFF param, M_{A} QE
    // BANFF param, M_{A} RES
    // BANFF param, CA5
    // BANFF param, SCCA (not for numu disappearance)
    // BANFF param, SCCV (not for numu disappearance)
    // BANFF param, Pf
    // BANFF param, CCothShape
    // BANFF param, Eb
    //
    
    if( util::IsUsingSplines(p) )
    {
        //
#if SPLINE_TYPE == -1
        return wght;
#endif

        wght_spl = this->FindWghtSpline(p,s,m,Ereco,Etrue);
        assert(wght_spl);
        double untuned_value = 0;
        // NEW XSEC
        if(util::IsRHCBANFF(p))
        {
            tb = fBANFF->InitialValue(p);
            untuned_value = (*fRHCUntunedCentralValues)[arr_index->RHCBANFFXSecSystIdx(p)];
        }
        else
        {
            tb = fBANFF->InitialValue(p);
            untuned_value = (*fUntunedCentralValues)[arr_index->BANFFXSecSystIdx(p)];
        }
        
        double spline_point = tb-untuned_value;
        
#ifdef check_spline_limits
        if (spline_point > wght_spl->GetXmax() || spline_point < wght_spl->GetXmin()) {
            LOG("NuXSec",pDEBUG) << "EVALUATING SPLINE OUTSIDE LIMITS!: "<<util::AsString(p)<<"\t"<<spline_point
                     <<"\tlimits: "<<wght_spl->GetXmin()<<" , "<<wght_spl->GetXmax();
        }
#endif
        
        
#ifndef __ALLOW_SPLINES_OUTSIDE_LIMIT__
        if (spline_point > wght_spl->GetXmax()) spline_point = wght_spl->GetXmax();
        if (spline_point < wght_spl->GetXmin()) spline_point = wght_spl->GetXmin();
#endif
        
        wght = wght_spl->Eval(spline_point);
    }
#endif
    
    
    // NORMALIZATION PARAMETERS
#if ( defined(__USE_NIWG2015__) || defined (__USE_NIWG2017__) )
    // - BANFF param, 2p2h       normalization
    // - BANFF param, BgRES      normalization
    // - BANFF param, CCcoh      normalization
    // - BANFF param, NCcoh      normalization
    // - BANFF param, NCoth      normalization
    // - BANFF param, NuEToNuMu  normalization
    // - BANFF param, CCNu       normalization
    // - BANFF param, CCNuBar    normalization
    else if( util::IsNormXSecParam(p) )
    {
      if (util::IsBeRPA(p))
      {
         wght = 1.0;
      }
      else if(util::IsNormCC(p) && (Ereco < 0.4 || Ereco > 0.6))
      {  // Apply spline for CC nu/nubar only in 400-600 MeV neutrino energy
         wght = 1.0;
      }
      else
      {
         wght = fBANFF->InitialValue(p);
      }
    }

#endif
    
#ifdef __show_banff_tuning_wght_mesg_debug__
    LOG("NuXSec",pDEBUG)
    << "\nSample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
    << ", Ereco = " << Ereco << " GeV"
    << ", Etrue = " << Etrue << " GeV => Weight ["
    << util::AsString(p) << "] = " << wght;
#endif
        // Hack to get around some small negative values in some of the fermi momentum response functions -1e-5 level
    if (wght < 0){
      wght = 2 * conv::kASmallNum;
      LOG("NuXSec", pWARN) << "Setting negative weight to zero "
                           << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
                           << " Current value = " << this->CurrentValue(p);
      if (wght_spl) {
        LOG("NuXSec", pWARN) << "Spline name, title: "<<wght_spl->GetName()
                             <<" , "<<wght_spl->GetTitle();
      }
    }

    if(fabs(wght) > conv::kMaxAllowedWght || wght < conv::kASmallNum) {
        //get the ereco & etrue bin indexes (to be used to check whether the spline is sensible)
        int ireco = -1;
        for (int i = 0; i < kNXsErecoBins[s]; i++) {
            if(util::XsecErecoBinEdges(s,i)   <= Ereco &&
               util::XsecErecoBinEdges(s,i+1) >  Ereco)
            {
                ireco = i;
                break;
            }
        }
        int itrue = -1;
        int is = arr_index->SuperkSampleIdx(s);
        int im = arr_index->SuperkModeIdx(s, m);
        for (int i = 0; i < kXsNEtrueBins; i++) {
            if(kXsEtrueBinEdges[is][im][i]   <= Etrue &&
               kXsEtrueBinEdges[is][im][i+1] >  Etrue)
            {
                itrue = i;
                break;
            }
        }
        LOG("NuXSec",pDEBUG) << "Got an absurd weight! "
        << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
        << ", Ereco = " << Ereco << " GeV"
        << ", Etrue = " << Etrue << " GeV =>"
        << " Bins are Ereco,Etrue = " << ireco << "," << itrue << " =>"
        << " Weight [" << util::AsString(p) << "] = " << wght
        << " Current value = " << this->CurrentValue(p);
        if (wght_spl) {
            LOG("NuXSec",pDEBUG) << "Spline name, title: "<<wght_spl->GetName()<<" , "<<wght_spl->GetTitle();
        }
        else {
            LOG("NuXSec",pDEBUG) <<"Not using spline";
        }
        util::Die(1);
    }
    return wght;

    if ( s == kSKsFCFVMultiRingELike ) LOG("NuXSec",pINFO) << "Weight for CC1" << wght;
}
//________________________________________________________________________________________
double NuXSec::FinTuningWght(SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue, double Q2) const
{
    double tot_wght = 1.;
    
    int nXSecSyst;
    PhysParamId_t p;
    bool is_banff_xsec_param;
    bool is_rhc = util::IsRHCSample(s);
    
    if(is_rhc)
        nXSecSyst = kNumOfRHCXSecSyst;
    else
        nXSecSyst = kNumOfXSecSyst;
    
    for(int isyst = 0; isyst < nXSecSyst; isyst++) {
        if(is_rhc) {
            p = kListOfRHCXSecSyst[isyst];
            // "Tuning" involves reweighting to BANFF best-fit values.
            // Consider only the cross-section systematics that were included in the BANFF fit.
            is_banff_xsec_param = kRHCXSecSystIncludedInBANFF[isyst];
        }
        else {
            p = kListOfXSecSyst[isyst];
            is_banff_xsec_param = kXSecSystIncludedInBANFF[isyst];
        }
        if(!is_banff_xsec_param) continue;
        // Is the current systematic parameter applied to current sample/mode end energy bin?
        bool apply = this->ApplySyst(p,s,m,Ereco,Etrue);
        if(!apply) continue;
        // Calculate weight
        double wght = this->FinTuningWght_1p(p,s,m,Ereco,Etrue,Q2);
        // LOG("NuXSec", pDEBUG) << wght << " Fin Tun Weight for " << util::AsString(s)
        //                       << " " << util::AsString(p)  << " " <<  util::AsString(m)
        //                       << " Ereco: " << Ereco << " Etrue: " << Etrue << " Theta: " << Q2;
        tot_wght *= wght;
    }
    
    return tot_wght;
}
//________________________________________________________________________________________
double NuXSec::FinTuningWght_1p(PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue, double Q2) const
{
#if defined(DISABLE_SCC) && BANFF_BINS != 2
    if (p == kNuisBANFFSCCA          || p == kRHCNuisBANFFSCCA       ||
        p == kNuisBANFFSCCV          || p == kRHCNuisBANFFSCCV       )
        return 1;
#endif
    
    // NEW XSEC
    
    //
    // This method provides weights to reweight MC to the BANFF best-fit parameters.
    // It handles the 4 BANFF cross-section parameters which are correlated with the BANFF
    // flux parameters. In the current incarnation of the BANFF fit, these parameters are:
    //
    // - NIWG2015 BANFF param, 2p2h        normalization
    // - NIWG2015 BANFF param, CA5
    // - NIWG2015 BANFF param, BgRES      normalization
    // - NIWG2015 BANFF param, M_{A} QE
    // - NIWG2015 BANFF param, M_{A} RES
    // - NIWG2015 BANFF param, SCCA (not for numu disappearance)
    // - NIWG2015 BANFF param, SCCV (not for numu disappearance)
    // - NIWG2015 BANFF param, Pf
    // - NIWG2015 BANFF param, CCothShape
    // - NIWG2015 BANFF param, Eb
    // - NIWG2015 BANFF param, CCcoh      normalization
    // - NIWG2015 BANFF param, NCcoh      normalization
    // - NIWG2015 BANFF param, NCoth      normalization
    // - NIWG2015 BANFF param, NuEToNuMu  normalization
    //
    
    TSpline3 * wght_spl = 0;
    
    double wght = 1.;
    double tb;
    double spline_point = 0;

    // RESPONSE FUNCTION PARAMETERS
    
    // NEW XSEC NIWG   
#if ( defined(__USE_NIWG2015__) || defined (__USE_NIWG2017__) )
    
    // NEW XSEC
    //
    // BANFF param, M_{A} QE
    // BANFF param, M_{A} RES
    // BANFF param, CA5
    // BANFF param, SCCA (not for numu disappearance)
    // BANFF param, SCCV (not for numu disappearance)
    // BANFF param, Pf
    // BANFF param, CCothShape
    // BANFF param, Eb
    //
    
    if(util::IsUsingSplines(p))
    {
#if SPLINE_TYPE == -1
        return wght;
#endif
        wght_spl = this->FindWghtSpline(p,s,m,Ereco,Etrue,Q2);
        assert(wght_spl);
        double untuned_value = 0;
        // NEW XSEC
        if(util::IsRHCBANFF(p))
        {
            tb = fBANFF->InitialValue(p);
            untuned_value = (*fRHCUntunedCentralValues)[arr_index->RHCBANFFXSecSystIdx(p)];
        }
        else
        {
            tb = fBANFF->InitialValue(p);
            untuned_value = (*fUntunedCentralValues)[arr_index->BANFFXSecSystIdx(p)];
        }

        spline_point = tb-untuned_value;

#ifdef check_spline_limits
        if (spline_point > wght_spl->GetXmax() || spline_point < wght_spl->GetXmin()) {
                    LOG("NuXSec",pERROR) <<"EVALUATING SPLINE OUTSIDE LIMITS!: "<<util::AsString(p)<<"\t"<<spline_point
                               <<"\tlimits: "<<wght_spl->GetXmin()<<" , "<<wght_spl->GetXmax() ;
          if (util::IsRHCBANFF(p)) {
            LOG("NuXSec",pERROR) << "tb = " << fBANFF->InitialValue(p) 
                                 << " untuned_value = " << (*fRHCUntunedCentralValues)[arr_index->RHCBANFFXSecSystIdx(p)];
          }
          if (util::IsBANFF(p)) {
            LOG("NuXSec",pERROR) << "tb = " << fBANFF->InitialValue(p) 
                                 << " untuned_value = " << (*fUntunedCentralValues)[arr_index->BANFFXSecSystIdx(p)];
          }
        }
#endif
        
#ifndef __ALLOW_SPLINES_OUTSIDE_LIMIT__
        if (spline_point > wght_spl->GetXmax()) spline_point = wght_spl->GetXmax();
        if (spline_point < wght_spl->GetXmin()) spline_point = wght_spl->GetXmin();
#endif
        wght = wght_spl->Eval(spline_point);
    }
#endif
    // NORMALIZATION PARAMETERS
#if ( defined(__USE_NIWG2015__) || defined (__USE_NIWG2017__) )
    else if( util::IsNormXSecParam(p) ) {
      if (util::IsBeRPA(p))
      {
         wght = 1.0;
      }
      else if(util::IsNormCC(p) && (Ereco < 0.4 || Ereco > 0.6))
      {  // Apply spline for CC nu/nubar only in 400-600 MeV neutrino energy
         wght = 1.0;
      }
      else
      {
         wght = fBANFF->InitialValue(p);
      }
    }
#endif
    
#ifdef __show_banff_tuning_wght_mesg_debug__
    LOG("NuXSec", pDEBUG) 
    << "\nSample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
    << ", Ereco = " << Ereco << " GeV"
    << ", Etrue = " << Etrue << " GeV"
    << ", Q2    = " << Q2 << " GeV=> Weight ["
    << util::AsString(p) << "] = " << wght;
#endif
    // Hack to get around some small negative values in some of the fermi momentum response functions -1e-5 level
    if (wght < 0){
      wght = 2 * conv::kASmallNum;
      LOG("NuXSec", pWARN) << "Setting negative weight to zero "
                           << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
                           << " Current value = " << this->CurrentValue(p);
      if (wght_spl) {
        LOG("NuXSec", pWARN) << "Spline name, title: "<<wght_spl->GetName()
                             <<" , "<<wght_spl->GetTitle();
      }
    }

    if(fabs(wght) > conv::kMaxAllowedWght || wght < conv::kASmallNum) {
        //get the ereco & etrue bin indexes (to be used to check whether the spline is sensible)
        int ireco = -1;
        for (int i = 0; i < kNXsErecoBins[s]; i++) {
            if(util::XsecErecoBinEdges(s,i)   <= Ereco &&
               util::XsecErecoBinEdges(s,i+1) >  Ereco)
            {
                ireco = i;
                break;
            }
        }
        int itrue = -1;
        int is = arr_index->SuperkSampleIdx(s);
        int im = arr_index->SuperkModeIdx(s, m);
        for (int i = 0; i < kXsNEtrueBins; i++) {
            if(kXsEtrueBinEdges[is][im][i]   <= Etrue &&
               kXsEtrueBinEdges[is][im][i+1] >  Etrue)
            {
                itrue = i;
                break;
            }
        }
        int iq2 = -1;
        for (int i = 0; i < kNXsQ2Bins[is]; i++) {
            if(util::XsecQ2BinEdges(s,i)   <= Q2 &&
               util::XsecQ2BinEdges(s,i+1) >  Q2)
            {
                iq2 = i;
                break;
            }
        }
        LOG("NuXSec", pERROR) << "Got an absurd weight! "
        << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
        << ", Ereco = " << Ereco << " GeV"
        << ", Etrue = " << Etrue << " GeV =>"
        << ", Q2    = " << Q2    << " GeV =>"
        << " Bins are Ereco,Etrue,Q2 = " << ireco << "," << itrue << "," << iq2 << " =>"
        << " Weight [" << util::AsString(p) << "] = " << wght
        << " At point " <<  spline_point 
        << " Current value = " << this->CurrentValue(p);
        if (wght_spl) {
            LOG("NuXSec",pDEBUG) << "Spline name, title: "<<wght_spl->GetName()<<" , "<<wght_spl->GetTitle();
        }
        else {
            LOG("NuXSec",pDEBUG) <<"Not using spline";
        }
        util::Die(1);
    }
    return wght;
}
//________________________________________________________________________________________
double NuXSec::SystWght(SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
    double tot_wght = 1.;
    
    int nXSecSyst;
    PhysParamId_t p;
    double dt;
    bool is_rhc = util::IsRHCSample(s);
    
    if(is_rhc)
        nXSecSyst = kNumOfRHCXSecSyst;
    else
        nXSecSyst = kNumOfXSecSyst;
    
    for(int isyst = 0; isyst < nXSecSyst; isyst++) {
        if(is_rhc) {
            p = kListOfRHCXSecSyst[isyst];
            // Was tweaked?
            if(kRHCXSecSystIncludedInBANFF[isyst])
                dt = fabs(fBANFF->CurrentValue(p) - fBANFF->InitialValue(p)); // banff parameters
            else
                dt = fabs(this->CurrentValue(p)); // non-banff normalisations
        }//RHC
        else {
            p = kListOfXSecSyst[isyst];
            // Was tweaked?
            if(kXSecSystIncludedInBANFF[isyst])
                dt = fabs(fBANFF->CurrentValue(p) - fBANFF->InitialValue(p)); // banff parameters
            else
                dt = fabs(this->CurrentValue(p)); // non-banff normalisations
        }//FHC
        // Was tweaked?
        if(dt < kASmallNum) continue;
        // Is the current systematic parameter applied to current sample/mode end energy bin?
        bool apply = this->ApplySyst(p,s,m,Ereco,Etrue);
        if(!apply) continue;
        // Calculate weight
        double wght = this->SystWght_1p(p,s,m,Ereco,Etrue);
        tot_wght *= wght;
    }//isyst
    return tot_wght;
}
//________________________________________________________________________________________
double NuXSec::SystWght_1p(PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
#if defined(DISABLE_SCC) && BANFF_BINS != 2
    if (p == kNuisBANFFSCCA          || p == kRHCNuisBANFFSCCA       ||
        p == kNuisBANFFSCCV          || p == kRHCNuisBANFFSCCV       )
        return 1;
#endif
    //
    // This method provides weights to account for systematic variations to the XSec params.
    // The weights are relative to MC reweighted to the BANFF best-fit point.
    // It handles both the BANFF and uncorelated XSec params.
    //
    
    // if( ! this->ApplySyst(p,s,m,Ereco,Etrue) ) return 1.; // check already made in function that calls this one
    
    TSpline3 * wght_spl = 0;
    
    double wght = 1.;
    double tb, tc;


    //
    // BANFF param w/ response functions:
    //
    
    // NEW XSEC NIWG
    
#if ( defined(__USE_NIWG2015__) || defined (__USE_NIWG2017__) )
    //
    // NEW XSEC
    //
    // BANFF param, M_{A} QE
    // BANFF param, M_{A} RES
    // BANFF param, CA5
    // BANFF param, SCCA (not for numu disappearance)
    // BANFF param, SCCV (not for numu disappearance)
    // BANFF param, Pf
    // BANFF param, CCothShape
    // BANFF param, Eb


    // Axial masses
    if (util::IsUsingSplines(p))
    {
#if SPLINE_TYPE == -1
        return wght;
#endif
        double untuned_value = 0;
        wght_spl = this->FindWghtSpline(p,s,m,Ereco,Etrue);
        assert(wght_spl);
        
        if(util::IsRHCBANFF(p))
        {
            tb = fBANFF->InitialValue(p);
            tc = fBANFF->CurrentValue(p);
            untuned_value = (*fRHCUntunedCentralValues)[arr_index->RHCBANFFXSecSystIdx(p)];
        }
        else {
            tb = fBANFF->InitialValue(p);
            tc = fBANFF->CurrentValue(p);
            untuned_value = (*fUntunedCentralValues)[arr_index->BANFFXSecSystIdx(p)];
        }
        double spl_tb = tb-untuned_value;
        double spl_tc = tc-untuned_value;
#ifdef check_spline_limits
        if (spl_tc > wght_spl->GetXmax() || spl_tc < wght_spl->GetXmin()) {
            LOG("NuXSec",pDEBUG) << "EVALUATING SPLINE OUTSIDE LIMITS!: "<<util::AsString(p)<<"\t"<<spl_tc<<"\tlimits: "<<wght_spl->GetXmin()<<" , "<<wght_spl->GetXmax();
        }
#endif
        
#ifndef __ALLOW_SPLINES_OUTSIDE_LIMIT__
        if (spl_tc > wght_spl->GetXmax()) spl_tc = wght_spl->GetXmax();
        if (spl_tc < wght_spl->GetXmin()) spl_tc = wght_spl->GetXmin();
#endif
        
        wght = wght_spl->Eval(spl_tc) / wght_spl->Eval(spl_tb);
    }    
    else if( util::IsNormXSecParam(p) )
    {
       if(util::IsBeRPA(p))
       {
          wght = 1.0;
       }
       else if(util::IsNormCC(p) && (Ereco < 0.4 || Ereco > 0.6))
       {  // Apply spline for CC nu/nubar only in 400-600 MeV neutrino energy
          wght = 1.0;
       }
       else
       {
          wght = fBANFF->CurrentValue(p) / fBANFF->InitialValue(p);
       }
    }
    else if(p == kNuisXSecRPA || p == kRHCNuisXSecRPA)
    {
        return 1; // This is implented in SKMC for now
        //~ #if SPLINE_TYPE == -1
        //~ return wght;
        //~ #endif
        //~ wght_spl = this->FindWghtSpline(p,s,m,Ereco,Etrue);
        //~ assert(wght_spl);
        //~ tc = this->CurrentValue(p);
        //~ wght = wght_spl->Eval(tc);
    }
#endif
    
#ifdef __show_syst_var_wght_mesg_debug__
    if (wght != 1.)
    {
        LOG("NuXSec",pDEBUG)
        << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
        << ", Ereco = " << Ereco << " GeV"
        << ", Etrue = " << Etrue << " GeV => Weight ["
        << util::AsString(p) << "] = " << wght;
    }
#endif
    
    if(fabs(wght) > conv::kMaxAllowedWght) {
        //get the ereco & etrue bin indexes (to be used to check whether the spline is sensible)
        int ireco = -1;
        for (int i = 0; i < kNXsErecoBins[s]; i++) {
            if(util::XsecErecoBinEdges(s,i)   <= Ereco &&
               util::XsecErecoBinEdges(s,i+1) >  Ereco)
            {
                ireco = i;
                break;
            }
        }
        int itrue = -1;
        int is = arr_index->SuperkSampleIdx(s);
        int im = arr_index->SuperkModeIdx(s, m);
        for (int i = 0; i < kXsNEtrueBins; i++) {
            if(kXsEtrueBinEdges[is][im][i]   <= Etrue &&
               kXsEtrueBinEdges[is][im][i+1] >  Etrue)
            {
                itrue = i;
                break;
            }
        }
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec/SystWght_1p >> Got an absurd weight! "
        << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
        << ", Ereco = " << Ereco << " GeV"
        << ", Etrue = " << Etrue << " GeV =>"
        << " Bins are Ereco,Etrue = " << ireco << "," << itrue << " =>"
        << " Weight [" << util::AsString(p) << "] = " << wght
        << " Current value = " << this->CurrentValue(p);
        if (wght_spl) {
            LOG("NuXSec",pDEBUG) << "Spline name, title: "<<wght_spl->GetName()<<" , "<<wght_spl->GetTitle();
        }
        else {
            LOG("NuXSec",pDEBUG) <<"Not using spline";
        }
        util::Die(1);
    }
    
    return wght;
}
//________________________________________________________________________________________
double NuXSec::SystWght(SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue, double Q2) const
{
    double tot_wght = 1.;
    
    int nXSecSyst;
    PhysParamId_t p;
    double dt;
    bool is_rhc = util::IsRHCSample(s);
    
    if(is_rhc)
        nXSecSyst = kNumOfRHCXSecSyst;
    else
        nXSecSyst = kNumOfXSecSyst;
    
    for(int isyst = 0; isyst < nXSecSyst; isyst++) {
        if(is_rhc) {
            p = kListOfRHCXSecSyst[isyst];
            // Was tweaked?
            if(kRHCXSecSystIncludedInBANFF[isyst])
                dt = 1.0;//fabs(fBANFF->CurrentValue(p) - fBANFF->InitialValue(p)); // banff parameters
            else
                dt = fabs(this->CurrentValue(p)); // non-banff normalisations
        }//RHC
        else {
            p = kListOfXSecSyst[isyst];
            // Was tweaked?
            if(kXSecSystIncludedInBANFF[isyst])
            {
                dt = 1.0;//fabs(fBANFF->CurrentValue(p) - fBANFF->InitialValue(p)); // banff parameters
            }
            else
                dt = fabs(this->CurrentValue(p)); // non-banff normalisations
        }//FHC
        // Was tweaked?
        if(dt < kASmallNum) continue;
        // Is the current systematic parameter applied to current sample/mode end energy bin?
        bool apply = this->ApplySyst(p,s,m,Ereco,Etrue);// Q2 update for different Q2 range
        if(!apply) continue;
        // Calculate weight
        double wght = this->SystWght_1p(p,s,m,Ereco,Etrue,Q2);
        LOG("NuXSec", pDEBUG) << wght << " Syst Weight for " << util::AsString(s) << " "
                              << util::AsString(p)  << " "   <<  util::AsString(m)
                              << " Ereco: " << Ereco << " Etrue: " << Etrue << "Theta: " << Q2 ;
        tot_wght *= wght;
    }//isyst
    return tot_wght;
}
//________________________________________________________________________________________
double NuXSec::SystWght_1p(PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue, double Q2) const
{
#if defined(DISABLE_SCC) && BANFF_BINS != 2
    if (p == kNuisBANFFSCCA          || p == kRHCNuisBANFFSCCA       ||
        p == kNuisBANFFSCCV          || p == kRHCNuisBANFFSCCV       )
        return 1;
#endif
    //
    // This method provides weights to account for systematic variations to the XSec params.
    // It handles both the BANFF and uncorelated XSec params.
    //
    
    TSpline3 * wght_spl = 0;
    
    double wght = 1.;
    double tc;
    
    double spl_tc(0);
    
    //
    // BANFF param w/ response functions:
    //
#if ( defined(__USE_NIWG2015__) || defined(__USE_NIWG2017__) )
    if (util::IsUsingSplines(p))
    {
#if SPLINE_TYPE == -1
        return wght;
#endif
        double untuned_value = 0;
        wght_spl = this->FindWghtSpline(p,s,m,Ereco,Etrue,Q2);
        assert(wght_spl);
        
        if(util::IsRHCBANFF(p))
        {
            tc = fBANFF->CurrentValue(p);
            untuned_value = (*fRHCUntunedCentralValues)[arr_index->RHCBANFFXSecSystIdx(p)];
        }
        else {
            tc = fBANFF->CurrentValue(p);
            untuned_value = (*fUntunedCentralValues)[arr_index->BANFFXSecSystIdx(p)];
        }
        spl_tc = tc-untuned_value;

#ifdef check_spline_limits
        if (spl_tc > wght_spl->GetXmax() || spl_tc < wght_spl->GetXmin()) {
          LOG("NuXSec", pERROR) <<"EVALUATING SPLINE OUTSIDE LIMITS!: "<<util::AsString(p)<<"\t"<<spl_tc
                                <<"\tlimits: "<<wght_spl->GetXmin()<<" , "<<wght_spl->GetXmax();
        }
#endif
        
#ifndef __ALLOW_SPLINES_OUTSIDE_LIMIT__
        if (spl_tc > wght_spl->GetXmax()) spl_tc = wght_spl->GetXmax();
        if (spl_tc < wght_spl->GetXmin()) spl_tc = wght_spl->GetXmin();
#endif
        
        wght = wght_spl->Eval(spl_tc);
        LOG("NuXSec", pDEBUG) << "Parameter: " << util::AsString(p) 
                              <<" fBANFF->CurrentValue(p): "<< wght_spl->Eval(spl_tc) << " @ " <<  spl_tc;
    }
    else if( util::IsNormXSecParam(p))
    {
      if (util::IsBeRPA(p))
      {
         wght = 1.0; 
      }
      else if(util::IsNormCC(p) && (Ereco < 0.4 || Ereco > 0.6))
      {  // Apply spline for CC nu/nubar only in 400-600 MeV neutrino energy
         wght = 1.0;
      }
      else
      {
         wght = fBANFF->CurrentValue(p);
      }
    }
    else if(p == kNuisXSecRPA || p == kRHCNuisXSecRPA) {
        return 1; // This is implented in SKMC for now
    }
#endif
    
#ifdef __show_syst_var_wght_mesg_debug__
    if (wght != 1.)
    {
        LOG("NuXSec",pDEBUG)
        << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
        << ", Ereco = " << Ereco << " GeV"
        << ", Etrue = " << Etrue << " GeV"
        << ", Q2    = " << Q2 << " GeV=> Weight ["
        << util::AsString(p) << "] = " << wght;
    }
#endif
    // Hack to get around some small negative values in some of the fermi momentum response functions -1e-5 level
    if (wght < 0){
      wght = 2 * conv::kASmallNum;
      LOG("NuXSec", pWARN) << "Setting negative weight to zero "
                           << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
                           << " Current value = " << this->CurrentValue(p);
      if (wght_spl) {
        LOG("NuXSec", pWARN) << "Spline name, title: "<<wght_spl->GetName()
                             <<" , "<<wght_spl->GetTitle();
      }
    }

    if(fabs(wght) > conv::kMaxAllowedWght) {
        //get the ereco & etrue bin indexes (to be used to check whether the spline is sensible)
        int ireco = -1;
        for (int i = 0; i < kNXsErecoBins[s]; i++) {
            if(util::XsecErecoBinEdges(s,i)   <= Ereco &&
               util::XsecErecoBinEdges(s,i+1) >  Ereco)
            {
                ireco = i;
                break;
            }
        }
        int itrue = -1;
        int is = arr_index->SuperkSampleIdx(s);
        int im = arr_index->SuperkModeIdx(s, m);
        for (int i = 0; i < kXsNEtrueBins; i++) {
            if(kXsEtrueBinEdges[is][im][i]   <= Etrue &&
               kXsEtrueBinEdges[is][im][i+1] >  Etrue)
            {
                itrue = i;
                break;
            }
        }
        
        int iq2 = -1;
        for (int i = 0; i < kNXsQ2Bins[is]; i++) {
            if(util::XsecQ2BinEdges(s,i)   <= Q2 &&
               util::XsecQ2BinEdges(s,i+1) >  Q2)
            {
                iq2 = i;
                break;
            }
        }
        
        LOG("NuXSec", pERROR) << "Got an absurd weight! "
        << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
        << ", Ereco = " << Ereco << " GeV"
        << ", Etrue = " << Etrue << " GeV =>"
        << ", Q2    = " << Q2    << " GeV =>"
        << " Bins are Ereco,Etrue,Q2 = " << ireco << "," << itrue << "," << iq2 << " =>"
        << " Weight [" << util::AsString(p) << "] = " << wght
        << " Current value = " << this->CurrentValue(p)
        << "\nEvaluated at: tc "<< spl_tc ;
        if (wght_spl) {
            LOG("NuXSec",pDEBUG) << "Spline name, title: "<<wght_spl->GetName()<<" , "<<wght_spl->GetTitle();
        }
        else {
            LOG("NuXSec",pDEBUG) <<"Not using spline";
        }
        util::Die(1);
    }
    //    wght = 1.0;    
    return wght;
}
//________________________________________________________________________________________
double NuXSec::CurrentValue(PhysParamId_t p) const
{
  if (util::IsBANFF(p)) {
    return fBANFF->CurrentValue(p);
  }
    // NEW XSEC NIWG
    
#if defined(__USE_NIWG2013__)
    switch(p) {
        case kNuisXSecFermiMom:
            return fFermiMomTwkDial;
            break;
        case kNuisXSecSpecFunc:
            return fSpecFuncTwkDial;
            break;
        case kNuisXSecWshape:
            return fWshapeTwkDial;
            break;
        case kNuisXSecShapeCCoth:
            return fShapeCCothTwkDial;
            break;
        case kNuisXSecPiLessDelta:
            return fPiLessDeltaTwkDial;
            break;
        case kNuisXSecBindE:
            return fBindETwkDial;
            break;
        case kNuisXSecNormCCcoh:
            return fNormCCcohTwkDial;
            break;
        case kNuisXSecNormNC1piC:
            return fNormNC1piCTwkDial;
            break;
        case kNuisXSecNormNCoth:
            return fNormNCothTwkDial;
            break;
        case kNuisXSecNormNCcoh:
            return fNormNCcohTwkDial;
            break;
        case kNuisXSecNormNuEToNuMu:
            return fNormNuEToNuMuTwkDial;
            break;
        case kNuisXSecNormAntiNuToNu:
            return fNormAntiNuToNuTwkDial;
            break;
        case kRHCNuisXSecFermiMom:
            return fRHCFermiMomTwkDial;
            break;
        case kRHCNuisXSecSpecFunc:
            return fRHCSpecFuncTwkDial;
            break;
        case kRHCNuisXSecWshape:
            return fRHCWshapeTwkDial;
            break;
        case kRHCNuisXSecShapeCCoth:
            return fRHCShapeCCothTwkDial;
            break;
        case kRHCNuisXSecPiLessDelta:
            return fRHCPiLessDeltaTwkDial;
            break;
        case kRHCNuisXSecBindE:
            return fRHCBindETwkDial;
            break;
        case kRHCNuisXSecNormCCcoh:
            return fRHCNormCCcohTwkDial;
            break;
        case kRHCNuisXSecNormNC1piC:
            return fRHCNormNC1piCTwkDial;
            break;
        case kRHCNuisXSecNormNCoth:
            return fRHCNormNCothTwkDial;
            break;
        case kRHCNuisXSecNormNCcoh:
            return fRHCNormNCcohTwkDial;
            break;
        case kRHCNuisXSecNormNuEToNuMu:
            return fRHCNormNuEToNuMuTwkDial;
            break;
        case kRHCNuisXSecNormAntiNuToNu:
            return fRHCNormAntiNuToNuTwkDial;
            break;
        case kNuisRHCToFHC:
            return fNormRHCToFHCTwkDial;
            break;
        default:
            return conv::kBadResult;
    }
    
#elif ( defined(__USE_NIWG2015__) || defined (__USE_NIWG2017__) )
    
    // NEW XSEC
    
    switch(p) {
        case kNuisXSecRPA:
        case kRHCNuisXSecRPA:
            return fRPATwkDial;
            break;
        default:
            return conv::kBadResult;
    }
#endif
    
    return conv::kBadResult;
}
//________________________________________________________________________________________
//
// Cross-section weight splines
//

//________________________________________________________________________________________
TSpline3 * NuXSec::FindWghtSpline(PhysParamId_t p, SuperkSample_t s, SuperkMode_t mp, double Ereco, double Etrue) const
{
    SuperkMode_t m = mp; // default
    // HACK: for `osc numu' templates, use `numu' response functions
    if      ( mp == kSKmOscNuMuCCQE       ) { m = kSKmNuMuCCQE;   } //both use the same binning
    else if ( mp == kSKmOscNuMuCC1pi      ) { m = kSKmNuMuCC1pi;  }
    else if ( mp == kSKmOscNuMuCCcoh      ) { m = kSKmNuMuCCcoh;  }
    else if ( mp == kSKmOscNuMuCCoth      ) { m = kSKmNuMuCCoth;  }
    
    else if ( mp == kSKmOscNuMuBarCCQE    ) { m = kSKmNuMuBarCCQE;  }//both use the same binning
    else if ( mp == kSKmOscNuMuBarCC1pi   ) { m = kSKmNuMuBarCC1pi; }
    else if ( mp == kSKmOscNuMuBarCCcoh   ) { m = kSKmNuMuBarCCcoh; }
    else if ( mp == kSKmOscNuMuBarCCoth   ) { m = kSKmNuMuBarCCoth; }
#ifdef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
    else if ( mp == kSKmOscNuEBarCCQE     ) { m = kSKmNuEBarCCQE;  }//problem if they dont use the same binning?
    else if ( mp == kSKmOscNuEBarCC1pi    ) { m = kSKmNuEBarCC1pi; }//fewer bins then numubar
    else if ( mp == kSKmOscNuEBarCCcoh    ) { m = kSKmNuEBarCCcoh; }
    else if ( mp == kSKmOscNuEBarCCoth    ) { m = kSKmNuEBarCCoth; }
#endif
    
    
    int idx = this->WghtSplineIdx(p,s,m,Ereco,Etrue);
    
    if(idx == -1 || idx == -2)
    {
      LOG("NuXSec",pERROR) << "Didn't find weight spline for "
			   << "sample = "      << util::AsString(s)
			   << ", mode = "      << util::AsString(m)
			   << ", parameter = " << util::AsString(p)
			   << ", Ereco = "     << Ereco << " GeV"
			   << ", Etrue = "     << Etrue << " GeV!"
			   << " idx = " << idx;
        util::Die(1); return 0;
    }
    
#ifdef __show_find_wght_spline_mesg_debug__
    LOG("NuXSec",pDEBUG) << "Looking for spline: "
			 << "sample = "      << util::AsString(s)
			 << ", mode = "      << util::AsString(m)
			 << ", parameter = " << util::AsString(p)
			 << ", Ereco = "     << Ereco << " GeV"
			 << ", Etrue = "     << Etrue << " GeV!"
			 << " idx = ";
#endif
    
    //Todo: need to cover case when sample is none of kSKsFCFV1RingMuLike,
    //kSKsFCFV1RingELike, kSKsRHCFCFV1RingMuLike or kSKsRHCFCFV1RingELike
    
    /// NEW MULTIRING

    switch(s){
      case( kSKsRHCFCFV1RingMuLike ) :  return fRHCXSecWghtSpl1Rmu[idx]; break; 
      case( kSKsRHCFCFV1RingELike  ) :  return fRHCXSecWghtSpl1Re[idx];  break; 
      case( kSKsFCFV1RingELike     ) :  return fXSecWghtSpl1Re[idx];     break;  
      case( kSKsFCFV1RingMuLike    ) :  return fXSecWghtSpl1Rmu[idx];    break; 
      case( kSKsFCFVMultiRingMuLike) :
        assert(fXSecWghtSplMultiRmu[idx]);
        return fXSecWghtSplMultiRmu[idx];
        break;
      case( kSKsFCFVMultiRingELike ) :
        assert(fXSecWghtSplMultiRe[idx]);
        return fXSecWghtSplMultiRe[idx];
        break;
    default                        :
      LOG("NuXSec", pERROR) << "Couldn't find spline for: " << AsString(s)
                            << ", mode = " << AsString(m) << ", for parameter = " << AsString(p);
      util::Die(1); return 0;
    }    
    
    ///
}
//________________________________________________________________________________________
TSpline3 * NuXSec::FindWghtSpline(PhysParamId_t p, SuperkSample_t s, SuperkMode_t mp, double Ereco, double Etrue, double Q2) const
{
#ifdef __use_1Rsplines_for_MultiR_sample__
  // ~#{ Hack for testing set multiring splines to single ring splines
  if (util::IsMultiRMuSample(s)) s = kSKsFCFV1RingMuLike;
  if (util::IsMultiReSample(s))  s = kSKsFCFV1RingELike;
  // }#~
#endif 		        

    SuperkMode_t m = mp; // default
    // HACK: for `osc numu' templates, use `numu' response functions
    if      ( mp == kSKmOscNuMuCCQE       ) { m = kSKmNuMuCCQE;    } //both use the same binning
    else if ( mp == kSKmOscNuMuCC1pi      ) { m = kSKmNuMuCC1pi;   }
    else if ( mp == kSKmOscNuMuCCcoh      ) { m = kSKmNuMuCCcoh;   }
    else if ( mp == kSKmOscNuMuCCoth      ) { m = kSKmNuMuCCoth;   }
    
    else if ( mp == kSKmOscNuMuBarCCQE    ) { m = kSKmNuMuBarCCQE;   }//both use the same binning
    else if ( mp == kSKmOscNuMuBarCC1pi   ) { m = kSKmNuMuBarCC1pi;  }
    else if ( mp == kSKmOscNuMuBarCCcoh   ) { m = kSKmNuMuBarCCcoh;  }
    else if ( mp == kSKmOscNuMuBarCCoth   ) { m = kSKmNuMuBarCCoth;  }
#ifdef __USE_OLD_OSC_NUEBAR_MC_TEMPLATES__
    else if ( mp == kSKmOscNuEBarCCQE     ) { m = kSKmNuEBarCCQE;  }//problem if they dont use the same binning?
    else if ( mp == kSKmOscNuEBarCC1pi    ) { m = kSKmNuEBarCC1pi; }//fewer bins then numubar
    else if ( mp == kSKmOscNuEBarCCcoh    ) { m = kSKmNuEBarCCcoh; }
    else if ( mp == kSKmOscNuEBarCCoth    ) { m = kSKmNuEBarCCoth; }
#endif
        
    int idx = this->WghtSplineIdx(p,s,m,Ereco,Etrue,Q2);
    
    if(idx == -1 || idx == -2)
    {
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec/FindWghtSpline >> Didn't find weight spline for "
        << "sample = "      << util::AsString(s)
        << ", mode = "      << util::AsString(m)
        << ", parameter = " << util::AsString(p)
        << ", Ereco = "     << Ereco << " GeV"
        << ", Etrue = "     << Etrue << " GeV!"
        << ", Q2    = "     << Q2    << " GeV!"
        << " idx = " << idx;
        util::Die(1); return 0;
    }
    
#ifdef __show_find_wght_spline_mesg_debug__
  LOG("NuXSec", pDEBUG) << "Looking for spline: "
    << "sample = "      << util::AsString(s)
    << ", mode = "      << util::AsString(m)
    << ", parameter = " << util::AsString(p)
    << ", Ereco = "     << Ereco << " GeV"
    << ", Etrue = "     << Etrue << " GeV!"
    << ", Q2    = "     << Q2    << " GeV!"
    << " idx = " << idx;
#endif
    
    //Todo: need to cover case when sample is none of kSKsFCFV1RingMuLike,
    //kSKsFCFV1RingELike, kSKsRHCFCFV1RingMuLike or kSKsRHCFCFV1RingELike
    
    // if(s == kSKsRHCFCFV1RingMuLike)
    //   return fRHCXSecWghtSpl1Rmu[idx];
    // else if(s == kSKsRHCFCFV1RingELike)
    //   return fRHCXSecWghtSpl1Re[idx];
    // else if(s == kSKsFCFV1RingELike)
    //   return fXSecWghtSpl1Re[idx];
    // else
    //   return fXSecWghtSpl1Rmu[idx];
    
    /// NEW MULTIRING (use 1 ring splines for multiring sample)
    if(s == kSKsRHCFCFV1RingMuLike)
    {
        return fRHCXSecWghtSpl1Rmu[idx];
    }
    else if(s == kSKsRHCFCFV1RingELike)
    {
        return fRHCXSecWghtSpl1Re[idx];
    }
    else if(s == kSKsFCFV1RingELike)
    {
        return fXSecWghtSpl1Re[idx];
    }
    else if(s == kSKsFCFV1RingMuLike)
    {
        return fXSecWghtSpl1Rmu[idx];
    }
    else if(s == kSKsFCFVMultiRingELike){ // multiring
        assert(fXSecWghtSplMultiRe[idx]);
        return fXSecWghtSplMultiRe[idx];
    }
    else {
        assert(fXSecWghtSplMultiRmu[idx]);
        return fXSecWghtSplMultiRmu[idx];
    }
    ///
    
}
//________________________________________________________________________________________
int NuXSec::WghtSplineIdx(PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
    // 1-D array of weight splines is basically a linearized version of a 5-D array
    // with the following indices:
    //
    //  0:  0 --> kNSuperkSamples
    //  1:  0 --> kNSuperkModes
    //  2:  0 --> kNumOfXSecSyst
    //  3:  0 --> kXsNErecoBins
    //  4:  0 --> kXsNEtrueBins
    //
    
    int is = arr_index->SuperkSampleIdx(s);
    if(is == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec/WghtSplineIdx >> Didn't find Idx for sample "
        << util::AsString(s)
        << ". util::ArrIdx returned "
        << arr_index->SuperkSampleIdx(s);
#endif //__show_wght_spline_idx_mesg_debug__
        return -1;
    }
    int im = arr_index->SuperkModeIdx(s,m);
    if(im == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec/WghtSplineIdx >> Didn't find Idx for mode "
        << util::AsString(m)
        << " in sample "
        << util::AsString(s)
        << ". util::ArrIdx returned "
        << arr_index->SuperkModeIdx(s,m);
#endif //__show_wght_spline_idx_mesg_debug__
        return -1;
    }
    int ip;
    if(util::IsRHCSample(s))
        ip = arr_index->RHCXSecSystIdx(p);
    else
        ip = arr_index->XSecSystIdx(p);
    if(ip == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec/WghtSplineIdx >> Didn't find Idx for xsec syst "
        << util::AsString(p)
        << ". util::ArrIdx returned "
        << arr_index->XSecSystIdx(p);
#endif//__show_wght_spline_idx_mesg_debug__
        return -1;
    }
    int ireco = -1;
    int nErecoBins = kNXsErecoBins[is];
    const double * ErecoBinEdges = util::XsecErecoBinEdges(s);
    for (int i = 0; i < nErecoBins ; i++) {
        if(ErecoBinEdges[i]   <= Ereco &&
           ErecoBinEdges[i+1] >  Ereco)
        {
            ireco = i;
            break;
        }
    }
    if(ireco == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec/WghtSplineIdx >> Didn't find ERecoBin for ereco value "
        << Ereco;
#endif//__show_wght_spline_idx_mesg_debug__
        return -1;
    }
    
    int itrue = -1;
    
    for (int i = 0; i < kXsNEtrueBins; i++) {
        if(kXsEtrueBinEdges[is][im][i]   <= Etrue &&
           kXsEtrueBinEdges[is][im][i+1] >  Etrue)
        {
            itrue = i;
            break;
        }
    }
    if(itrue == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec/WghtSplineIdx >> Didn't find ETrueBin for etrue value "
        << Etrue
        << " in mode " << m << " "
        << util::AsString(m)
        << " in sample " << s << " "
        << util::AsString(s);
#endif //__show_wght_spline_idx_mesg_debug__
        return -2;
    }
    
    int idx = itrue
    + ireco * kXsNEtrueBins
    + ip    * kXsNEtrueBins * nErecoBins
    + im    * kXsNEtrueBins * nErecoBins * kNumOfXSecSyst;
    //    + is    * kXsNEtrueBins * util::NumXsecErecoBins(s) * kNumOfXSecSyst * kNSuperkModes;
    
    
#ifdef __show_wght_spline_idx_mesg_debug__
    LOG("NuXSec",pDEBUG) << "Exp/NuXSec/WghtSplineIdx >> Spline information:"
    << "\tEtrue: "<<Etrue<<" ("<<itrue<<") , Ereco: "<<Ereco<<" ("<<ireco<<") \n"
    << " in mode " << m << " "
    << util::AsString(m)
    << " in sample " << s << " "
    << util::AsString(s)
    << " has index "<<idx<<"\n";
#endif //__show_wght_spline_idx_mesg_debug__
    
    return idx;
}

int NuXSec::WghtSplineIdx( PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue, double Q2) const
{
    // 1-D array of weight splines is basically a linearized version of a 6-D array
    // with the following indices:
    //
    //  0:  0 --> kNSuperkSamples
    //  1:  0 --> kNSuperkModes
    //  2:  0 --> kNumOfXSecSyst
    //  3:  0 --> kXsNErecoBins
    //  4:  0 --> kNXsQ2Bins
    //  5:  0 --> kXsNEtrueBins
    //
    
    int is = arr_index->SuperkSampleIdx(s);
    if(is == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec/WghtSplineIdx >> Didn't find Idx for sample "
        << util::AsString(s)
        << ". util::ArrIdx returned "
        << arr_index->SuperkSampleIdx(s);
#endif //__show_wght_spline_idx_mesg_debug__
        return -1;
    }
    int im = arr_index->SuperkModeIdx(s,m);
    if(im == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec/WghtSplineIdx >> Didn't find Idx for mode "
        << util::AsString(m)
        << " in sample "
        << util::AsString(s)
        << ". util::ArrIdx returned "
        << arr_index->SuperkModeIdx(s,m);
#endif //__show_wght_spline_idx_mesg_debug__
        return -1;
    }
    int ip;
    if(util::IsRHCSample(s))
        ip = arr_index->RHCXSecSystIdx(p);
    else
        ip = arr_index->XSecSystIdx(p);
    if(ip == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec/WghtSplineIdx >> Didn't find Idx for xsec syst "
        << util::AsString(p)
        << ". util::ArrIdx returned "
        << arr_index->XSecSystIdx(p);
#endif//__show_wght_spline_idx_mesg_debug__
        return -1;
    }
    int ireco = -1;
    int nErecoBins = kNXsErecoBins[is];
    const double * ErecoBinEdges = util::XsecErecoBinEdges(s);
    for (int i = 0; i < nErecoBins ; i++) {
        if(ErecoBinEdges[i]   <= Ereco &&
           ErecoBinEdges[i+1] >  Ereco)
        {
            ireco = i;
            break;
        }
    }
    if(ireco == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec/WghtSplineIdx >> Didn't find ERecoBin for ereco value "
        << Ereco;
#endif//__show_wght_spline_idx_mesg_debug__
        return -1;
    }
    
    int itrue = -1;
    
    for (int i = 0; i < kXsNEtrueBins; i++) {
        if(kXsEtrueBinEdges[is][im][i]   <= Etrue &&
           kXsEtrueBinEdges[is][im][i+1] >  Etrue)
        {
            itrue = i;
            break;
        }
    }
    if(itrue == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec/WghtSplineIdx >> Didn't find ETrueBin for etrue value "
        << Etrue
        << " in mode " << m << " "
        << util::AsString(m)
        << " in sample " << s << " "
        << util::AsString(s);
#endif //__show_wght_spline_idx_mesg_debug__
        return -2;
    }
    
    int iq2 = -1;
    
    for (int i = 0; i < kNXsQ2Bins[is]; i++) {
        if(util::XsecQ2BinEdges(s,i)   <= Q2 &&
           util::XsecQ2BinEdges(s,i+1) >  Q2)
        {
            iq2 = i;
            break;
        }
    }
    if(iq2 == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec/WghtSplineIdx >> Didn't find Q2 bin for Q2 value "
        << Q2
        << " in mode " << m << " "
        << util::AsString(m)
        << " in sample " << s << " "
        << util::AsString(s);
#endif //__show_wght_spline_idx_mesg_debug__
        return -2;
    }
    
    int idx = itrue
    + ireco * kXsNEtrueBins
    + iq2   * kXsNEtrueBins * nErecoBins
    + ip    * kXsNEtrueBins * nErecoBins * kNXsQ2Bins[is]
    + im    * kXsNEtrueBins * nErecoBins * kNXsQ2Bins[is] * kNumOfXSecSyst;
    //    + is    * kXsNEtrueBins * util::NumXsecErecoBins(s) * kNumOfXSecSyst * kNSuperkModes;
    
    
#ifdef __show_wght_spline_idx_mesg_debug__
    LOG("NuXSec",pDEBUG) << "Exp/NuXSec/WghtSplineIdx >> Spline information:"
    << "\tEtrue: "<<Etrue<<" ("<<itrue<<") , Ereco: "<<Ereco<<" ("<<ireco<<"), Q2: "<<Q2<<" ("<<iq2<<") \n"
    << " in mode " << m << " "
    << util::AsString(m)
    << " in sample " << s << " "
    << util::AsString(s)
    << " has index "<<idx<<"\n";
#endif //__show_wght_spline_idx_mesg_debug__
    
    return idx;
}
//________________________________________________________________________________________
bool NuXSec::ApplySyst( PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
    bool apply2mode = this->ApplySystToSKMode(p,s,m);
    if (not apply2mode) return false;
    bool apply2bin  = this->ApplySystToEnergyBin(p,Ereco,Etrue);
    return apply2bin;
}
//________________________________________________________________________________________
bool NuXSec::ApplySystToSKMode(PhysParamId_t p, SuperkSample_t s, SuperkMode_t m) const
{
    int im = arr_index->SuperkModeIdx(s,m);
    if (im == -1) return false;
    
    int isyst;
    if(util::IsRHCSample(s)) {
        isyst = arr_index->RHCXSecSystIdx(p);
        if (isyst == -1) return false;
    }
    else {
        isyst = arr_index->XSecSystIdx(p);
        if (isyst == -1) return false;
    }
    
    int is = arr_index->SuperkSampleIdx(s);
    
    bool apply = kXsApplyXSecSystToSKMode[is][im][isyst];
    
    return apply;
}
//________________________________________________________________________________________
bool NuXSec::ApplySystToEnergyBin(PhysParamId_t p, double Ereco, double Etrue) const
{
    int isyst;
    bool apply;
    
    if(IsRHCSystematic(p)) {
        isyst = arr_index->RHCXSecSystIdx(p);
        if (isyst == -1) return false;
        apply = (Ereco >= kRHCXSecSystErecoCut[isyst][0] && Ereco < kRHCXSecSystErecoCut[isyst][1]) &&
        (Etrue >= kRHCXSecSystEtrueCut[isyst][0] && Etrue < kRHCXSecSystEtrueCut[isyst][1]);
    }
    else {
        isyst = arr_index->XSecSystIdx(p);
        if (isyst == -1) return false;
        apply = (Ereco >= kXSecSystErecoCut[isyst][0] && Ereco < kXSecSystErecoCut[isyst][1]) &&
        (Etrue >= kXSecSystEtrueCut[isyst][0] && Etrue < kXSecSystEtrueCut[isyst][1]);
    }
    
    return apply;
}
//________________________________________________________________________________________
void NuXSec::BuildXsBinningMCTemplates(TH2D * input_mctempl_array[kns][knm])
{
#ifdef __show_init_mesg_info__
    LOG("NuXSec",pDEBUG) << "BuildXsBinningMCTemplates ";
#endif
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
#ifdef __show_build_xsbinning_template_mesg_info__
            LOG("NuXSec",pDEBUG) << "Building Xsec binning MC template for "
            << "sample = " << util::AsString(s)
            << ", mode = " << util::AsString(m);
#endif
            const TH2D * mctempl = input_mctempl_array[is][im];
            assert (mctempl);
            TH2D * mctempl_xsbinning = fXsBinningTempl[is][im];
            
            for(int ir = 0; ir <= mctempl->GetNbinsX(); ir++) {
                for(int it = 0; it <= mctempl->GetNbinsY(); it++) {
                    if(mctempl->GetBinContent(ir,it) > kASmallNum) {
                        mctempl_xsbinning->Fill(mctempl->GetXaxis()->GetBinCenter(ir), mctempl->GetYaxis()->GetBinCenter(it), mctempl->GetBinContent(ir,it));
#ifdef __show_build_xsbinning_template_mesg_debug__
                        int bx = 0, by = 0, bz = 0;
                        int bg = mctempl_xsbinning->FindBin(mctempl->GetXaxis()->GetBinCenter(ir), mctempl->GetYaxis()->GetBinCenter(it));
                        mctempl_xsbinning->GetBinXYZ(bg, bx, by, bz);
                        LOG("NuXSec",pDEBUG) << std::fixed
                        << "Exp/NuXSec >> Filled Xsec binning MC template for "
                        << "ereco = " << mctempl->GetXaxis()->GetBinCenter(ir)
                        << " etrue = " << mctempl->GetYaxis()->GetBinCenter(it)
                        << "\n Put into bins: ereco "
                        << mctempl_xsbinning->GetXaxis()->GetBinLowEdge(bx) << " - "
                        << mctempl_xsbinning->GetXaxis()->GetBinUpEdge(bx) << "\t etrue "
                        << mctempl_xsbinning->GetYaxis()->GetBinLowEdge(by) << " - "
                        << mctempl_xsbinning->GetYaxis()->GetBinUpEdge(by);
#endif
                    }
                }//it
            }//ir
            
        }//im
    }//is
    
#ifdef __save_ori_templates_init_debug__
    TFile out("XsBinningTemplates.root","RECREATE");
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            fXsBinningTempl[is][im]->Write();
        }
    }
    out.Close();
#endif
}
//________________________________________________________________________________________
void NuXSec::BuildXsBinningMCTemplates(TH3D * input_mctempl_array[kns][knm])
{
#ifdef __show_init_mesg_info__
  LOG("NuXSec", pINFO) << "BuildXsBinningMCTemplates(2D) ";
#endif
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
#ifdef __show_build_xsbinning_template_mesg_info__
            LOG("NuXSec",pDEBUG) << "Building Xsec binning MC template for "
            << "sample = " << util::AsString(s)
            << ", mode = " << util::AsString(m);
#endif
            const TH3D * mctempl = input_mctempl_array[is][im];
            assert (mctempl);
            TH3D * mctempl_xsbinning = fXsBinningTempl2D[is][im];
            
            for(int ir = 0; ir <= mctempl->GetNbinsX(); ir++) {
                for(int it = 0; it <= mctempl->GetNbinsY(); it++) {
                    for(int iv = 0; iv <= mctempl->GetNbinsZ(); iv++) {
                        
                        if(mctempl->GetBinContent(ir,it,iv) > kASmallNum) {
                            mctempl_xsbinning->Fill(mctempl->GetXaxis()->GetBinCenter(ir), mctempl->GetYaxis()->GetBinCenter(it),mctempl->GetZaxis()->GetBinCenter(iv), mctempl->GetBinContent(ir,it,iv));
#ifdef __show_build_xsbinning_template_mesg_debug__
                            int bx = 0, by = 0, bz = 0;
                            int bg = mctempl_xsbinning->FindBin(mctempl->GetXaxis()->GetBinCenter(ir), mctempl->GetYaxis()->GetBinCenter(it),mctempl->GetZaxis()->GetBinCenter(iv));
                            mctempl_xsbinning->GetBinXYZ(bg, bx, by, bz);
                            LOG("NuXSec",pDEBUG)
                            << "Filled Xsec binning MC template for "
                            << " Q2    = " << mctempl->GetXaxis()->GetBinCenter(ir)
                            << " ereco = " << mctempl->GetYaxis()->GetBinCenter(it)
                            << " etrue = " << mctempl->GetZaxis()->GetBinCenter(iv);
                            << "\n Put into bins: ereco "
                            << mctempl_xsbinning->GetXaxis()->GetBinLowEdge(bx) << " - "
                            << mctempl_xsbinning->GetXaxis()->GetBinUpEdge(bx) << "\t Q2 "
                            << mctempl_xsbinning->GetYaxis()->GetBinLowEdge(by) << " - "
                            << mctempl_xsbinning->GetYaxis()->GetBinUpEdge(by)  << "\t ereco "
                            << mctempl_xsbinning->GetZaxis()->GetBinLowEdge(bz) << " - "
                            << mctempl_xsbinning->GetZaxis()->GetBinUpEdge(bz) << "\t etrue ";
#endif
                        }
                    }//iv
                }//it
            }//ir
        }//im
    }//is
    
#ifdef __save_ori_templates_init_debug__
    TFile out("XsBinningTemplates.root","RECREATE");
    for(int is = 0; is < kns; is++) {
        SuperkSample_t s = kSKSamples[is];
        if(s == kSKsUndefined) continue;
        for(int im = 0; im < knm; im++) {
            SuperkMode_t m = kSKReactionModes[is][im];
            if(m == kSKmUndefined) continue;
            fXsBinningTempl2D[is][im]->Write();
        }
    }
    out.Close();
#endif
}

//________________________________________________________________________________________
bool NuXSec::IsThereAnEvent(SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
    int is = arr_index->SuperkSampleIdx(s);
    int im = arr_index->SuperkModeIdx(s,m);
    assert(is != -1 && im != -1);
    
    TH2D * mctempl_xsbinning = fXsBinningTempl[is][im];
    assert(mctempl_xsbinning);
    
    int ir = mctempl_xsbinning->GetXaxis()->FindBin(Ereco);
    int it = mctempl_xsbinning->GetYaxis()->FindBin(Etrue);
    
    bool ev = mctempl_xsbinning->GetBinContent(ir,it) > kASmallNum;
    
#ifdef __show_isanevent_mesg_debug__
    if(ev)
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec >> There are " << mctempl_xsbinning->GetBinContent(ir,it)
        << " events for sample = " << util::AsString(s)
        << ", mode = " << util::AsString(m)
        << ", Etrue = " << Etrue
        << ", Ereco = " << Ereco;
#endif
    return ev;
}
//________________________________________________________________________________________
bool NuXSec::IsThereAnEvent(SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue, double Q2) const
{
    int is = arr_index->SuperkSampleIdx(s);
    int im = arr_index->SuperkModeIdx(s,m);
    assert(is != -1 && im != -1);
    
    TH3D * mctempl_xsbinning = fXsBinningTempl2D[is][im];
    assert(mctempl_xsbinning);
    
    int ir = mctempl_xsbinning->GetXaxis()->FindBin(Q2);
    int it = mctempl_xsbinning->GetYaxis()->FindBin(Ereco);
    int iv = mctempl_xsbinning->GetZaxis()->FindBin(Etrue);
    
    double ev = mctempl_xsbinning->GetBinContent(ir,it,iv);
    
#ifdef __show_isanevent_mesg_debug__
    if(ev > kASmallNum)
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec >> There are " << mctempl_xsbinning->GetBinContent(ir,it,iv)
        << " events for sample = " << util::AsString(s)
        << ", mode = " << util::AsString(m)
        << ", Etrue = " << Etrue
        << ", Ereco = " << Ereco
        << ", Q2    = " << Q2;
#endif
    return ev;
}
//________________________________________________________________________________________
int NuXSec::HowManyEvents(SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
    int is = arr_index->SuperkSampleIdx(s);
    int im = arr_index->SuperkModeIdx(s,m);
    assert(is != -1 && im != -1);
    
    TH2D * mctempl_xsbinning = fXsBinningTempl[is][im];
    assert(mctempl_xsbinning);
    
    int ir = mctempl_xsbinning->GetXaxis()->FindBin(Ereco);
    int it = mctempl_xsbinning->GetYaxis()->FindBin(Etrue);
    
    int ev = (int)mctempl_xsbinning->GetBinContent(ir,it);
    
#ifdef __show_isanevent_mesg_debug__
    if(ev)
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec >> There are " << ev
        << " events for sample = " << util::AsString(s)
        << ", mode = " << util::AsString(m)
        << ", Etrue = " << Etrue
        << ", Ereco = " << Ereco;
#endif
    
    return ev;
}
//________________________________________________________________________________________

int NuXSec::HowManyEvents(SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue, double Q2) const
{
    int is = arr_index->SuperkSampleIdx(s);
    int im = arr_index->SuperkModeIdx(s,m);
    assert(is != -1 && im != -1);
    
    TH3D * mctempl_xsbinning = fXsBinningTempl2D[is][im];
    assert(mctempl_xsbinning);
    
    int ir = mctempl_xsbinning->GetXaxis()->FindBin(Q2);
    int it = mctempl_xsbinning->GetYaxis()->FindBin(Ereco);
    int iv = mctempl_xsbinning->GetZaxis()->FindBin(Etrue);
    
    double ev = mctempl_xsbinning->GetBinContent(ir,it,iv);
    
#ifdef __show_isanevent_mesg_debug__
    if((int)ev)
        LOG("NuXSec",pDEBUG) << "Exp/NuXSec >> There are " << ev
        << " events for sample = " << util::AsString(s)
        << ", mode = " << util::AsString(m)
        << ", Etrue = " << Etrue
        << ", Ereco = " << Ereco
        << ", Q2    = " << Q2;
#endif
    
    return ev;
}

//________________________________________________________________________________________
TH2D * NuXSec::Template(SuperkSample_t s, SuperkMode_t m) const
{
    int sidx = arr_index->SuperkSampleIdx(s);
    if(sidx == -1) {
        return 0;
    }
    int midx = arr_index->SuperkModeIdx(s,m);
    if(midx == -1) {
        return 0;
    }
    
    TH2D * h = 0;
    h = fXsBinningTempl[sidx][midx];
    assert(h);
    
    return h;
}
//________________________________________________________________________________________
TH3D * NuXSec::Template2D(SuperkSample_t s, SuperkMode_t m) const
{
    int sidx = arr_index->SuperkSampleIdx(s);
    if(sidx == -1) {
        return 0;
    }
    int midx = arr_index->SuperkModeIdx(s,m);
    if(midx == -1) {
        return 0;
    }
    
    TH3D * h = 0;
    h = fXsBinningTempl2D[sidx][midx];
    assert(h);
    
    return h;
}
//________________________________________________________________________________________
void NuXSec::Tweak(conv::PhysParamId_t p , double twkdial)
{
    // NEW XSEC --> comment all because many xsec syst are now in the BANFF and the others are not used any longer
    // -->  not commented any longer
#if defined(__USE_NIWG2015__)
    // RPA is a discrete value (0 or 1)
    if      (p==kNuisXSecRPA      ){
        if( fabs(twkdial)>conv::kASmallNum && fabs(twkdial-1)>conv::kASmallNum ){
            if(twkdial < 0.5) twkdial = 0.;
            else              twkdial = 1.;
            LOG("NuXSec",pDEBUG) << "Tweak Dial RPA must be 0 or 1 ";
            LOG("NuXSec",pDEBUG) << "Tweak Dial RPA restored to " << twkdial;
        }
        fRPATwkDial      = twkdial;
    }
    else if (p==kRHCNuisXSecRPA   ){
        if( fabs(twkdial)>conv::kASmallNum && fabs(twkdial-1)>conv::kASmallNum ){
            if(twkdial < 0.5) twkdial = 0.;
            else              twkdial = 1.;
            LOG("NuXSec",pDEBUG) << "Tweak Dial RHC RPA must be 0 or 1 ";
            LOG("NuXSec",pDEBUG) << "Tweak Dial RHC RPA restored to " << twkdial;
        }
        fRHCRPATwkDial   = twkdial;
    }

    else {
        LOG("NuXSec",pDEBUG)<<"NuXSec/Tweak >> Failed to tweak parameter "<<util::AsString(p);
        util::Die(1);
    }
#endif    
}
//________________________________________________________________________________________
void NuXSec::SyncRHCToFHCSystematics()
{
    fRHCRPATwkDial            = fRPATwkDial           ; // NEW XSEC
    fRHCFermiMomTwkDial       = fFermiMomTwkDial      ;
    fRHCSpecFuncTwkDial       = fSpecFuncTwkDial      ;
    fRHCWshapeTwkDial         = fWshapeTwkDial        ;
    fRHCShapeCCothTwkDial     = fShapeCCothTwkDial    ;
    fRHCPiLessDeltaTwkDial    = fPiLessDeltaTwkDial   ;
    fRHCBindETwkDial          = fBindETwkDial         ;
    fRHCNormCCcohTwkDial      = fNormCCcohTwkDial     ;
    fRHCNormNC1piCTwkDial     = fNormNC1piCTwkDial    ;
    fRHCNormNCcohTwkDial      = fNormNCcohTwkDial     ;
    fRHCNormNCothTwkDial      = fNormNCothTwkDial     ;
    fRHCNormNuEToNuMuTwkDial  = fNormNuEToNuMuTwkDial ;
    fRHCNormAntiNuToNuTwkDial = fNormAntiNuToNuTwkDial;
}
//________________________________________________________________________________________
std::string NuXSec::BuildSplineName(conv::SuperkSample_t s, conv::SuperkMode_t m, conv::PhysParamId_t p) const
{
   ostringstream splname;
   splname << "spl_" << util::AsStringCompact(s) << "_" <<
      util::AsStringSpline(p)  << "_" << util::AsStringCompact(m);

   return splname.str();
}

//________________________________________________________________________________________
std::string NuXSec::BuildSplineName(conv::SuperkSample_t s, conv::SuperkMode_t m, conv::PhysParamId_t p, int reco_bin, int true_bin) const
{
#ifdef __use_1Rsplines_for_MultiR_sample__
  // ~#{ Hack for testing set multiring splines to single ring splines
  if (util::IsMultiRMuSample(s)) s = kSKsFCFV1RingMuLike;
  if (util::IsMultiReSample(s))  s = kSKsFCFV1RingELike;
  // }#~
#endif
  
  ostringstream splname;
  splname << "spl_";
#if ( defined(__USE_NIWG2015__) || defined (__USE_NIWG2017__) )
    splname << util::AsStringCompact(s) <<"_"
    << util::AsStringSpline(p)  <<"_"
    << util::AsStringCompact(m) <<"_"
    << true_bin << "_"
    << reco_bin;
#endif
    return splname.str();
}
//________________________________________________________________________________________
std::string NuXSec::BuildSplineName(conv::SuperkSample_t s, conv::SuperkMode_t m, conv::PhysParamId_t p, int reco_bin, int true_bin, int q2bin) const
{
#ifdef __use_1Rsplines_for_MultiR_sample__
  // ~#{ Hack for testing set multiring splines to single ring splines
  if (util::IsMultiRMuSample(s)) s = kSKsFCFV1RingMuLike;
  if (util::IsMultiReSample(s))  s = kSKsFCFV1RingELike;
  // }#~
#endif 		      
    ostringstream splname;
    splname << "spl_";
#if ( defined(__USE_NIWG2015__) || defined (__USE_NIWG2017__) )
    splname << util::AsStringCompact(s) <<"_"
    << util::AsStringSpline(p)  <<"_"
    << util::AsStringCompact(m) <<"_"
    << true_bin << "_"
    << reco_bin  << "_"
    << q2bin;
#endif
    return splname.str();
}
//________________________________________________________________________________________
