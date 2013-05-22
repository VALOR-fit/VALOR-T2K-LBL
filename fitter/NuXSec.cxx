//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @ 18/2/13 LE
 Added parameters for NC1pi0 (BANFF) and NCcoh normalization, added splines for 1Re as well.

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__
//#define __show_init_mesg_debug__
//#define __show_init_mesg_debug2__
//#define __show_banff_tuning_wght_mesg_info__
//#define __show_banff_tuning_wght_mesg_debug__
//#define __show_syst_var_wght_mesg_debug__
//#define __show_syst_var_wght_mesg_debug2__
//#define __show_build_xsbinning_template_mesg_info__
//#define __show_build_xsbinning_template_mesg_debug__
//#define __show_wght_spline_idx_mesg_debug__
//#define __show_isanevent_mesg_debug__
//#define __save_ori_templates_init_debug__
#define __skip_unexpected_spline_test_in_readsplines__ // skip looking for splines we don't expect (due to our empty bins) - speedy!
#define __use_flat_splines_for_missing_splines__ // DON'T USE! will use flat splines for missing splines. Only used between SK MC versions (new SKMC, old splines)

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <TMath.h>

#include "NuXSec.h"
#include "BANFF.h"
#include "Utilities.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;

using namespace t2kana;
using namespace t2kana::conv;
using namespace t2kana::util;

ClassImp(NuXSec)

//________________________________________________________________________________________
NuXSec::NuXSec()
{
  this->Init();
}
//________________________________________________________________________________________
NuXSec::~NuXSec()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void NuXSec::Init(void)
{
  fBANFF = 0;

  for(int ispl = 0; ispl < knwspl1Rmu; ispl++) {
     fXSecWghtSpl1Rmu[ispl] = 0;
  }
  for(int ispl = 0; ispl < knwspl1Re; ispl++) {
     fXSecWghtSpl1Re[ispl] = 0;
  }

  for(int is=0; is<kns; is++) {
    for(int im=0; im<knm; im++) {
      fXsBinningTempl [is][im] = 0;
    }
  }

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

  for(int is=0; is<kns; is++) {
    for(int im=0; im<knm; im++) {
      if (fXsBinningTempl [is][im]) delete fXsBinningTempl [is][im]; 
    }
  }
}
//________________________________________________________________________________________
bool NuXSec::ReadSplines(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/NuXSec >> Reading neutrino generator inputs (splines only) from input XML file." << endl;
#endif

  //
  // Load cross-section weight splines
  //

#if SPLINE_TYPE == -1
  cout << "Exp/NuXSec >> Special mode turned on in which splines are not read, and weights involving splines are set to 1" << endl;
  return true;
#elif SPLINE_TYPE == 1
  TFile * spline_file_1Rmu = util::XML2TFile(xml,root_node,"xsec/weight_splines/one_mu_ring");
  assert(spline_file_1Rmu);
  TFile * spline_file_1Re = util::XML2TFile(xml,root_node,"xsec/weight_splines/one_e_ring");
  assert(spline_file_1Re);
#endif

  // Make a flat spline
#ifdef __use_flat_splines_for_missing_splines__
  double x[5] = {-5,-2,0,+2,+5};
  double y[5] = {1,1,1,1,1};
  int n = 5;
  TSpline3 *flatspline = new TSpline3("",x,y,n);
#endif //__use_flat_splines_for_missing_splines__

  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    // loop over modes for current sample
    for(int im = 0; im < knm; im++) {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;

       // loop over cross-section systematic parameters
       for(int ip = 0; ip < kNumOfXSecSyst; ip++) {
          PhysParamId_t p = kListOfXSecSyst[ip];          
          if(p == kPhypUndefined) continue;
          // do we need weight splines for the current systematic parameter?
          bool need_splines = kXSecSystHasWghtSplines[ip];
          if(!need_splines) continue;
          // is the current systematic parameter going to be applied to current sample/mode?
          bool is_applied = kXsApplyXSecSystToSKMode[is][im][ip];
          if(!is_applied) continue;
          // loop over the var.-size  Ereco and Etrue bins for which weight splines were pre-computed
          for(int ir = 0; ir < util::NumXsecErecoBins(s); ir++) {
	    double Ereco = util::XsecErecoBinEdges(s,ir) + kOneMeV;
             for(int it = 0; it < kXsNEtrueBins; it++) {
                double Etrue = kXsEtrueBinEdges[is][im][it] + kOneMeV;

                bool bin_is_beyond_end_of_range = (kXsEtrueBinEdges[is][im][it+1] < 0.);
                if(bin_is_beyond_end_of_range) continue;

		bool bin_is_not_empty = this->IsThereAnEvent(s, m, Ereco, Etrue);

#ifdef __skip_unexpected_spline_test_in_readsplines__
		if(!bin_is_not_empty)                 
		  continue;
#endif
                // Build spline name
                ostringstream splname;
                splname << "spl_"
		  //<< util::AsStringCompact(kSKsFCFV1RingMuLike) << "_"
			<< util::AsStringCompact(s) << "_"
                        << util::AsStringCompact(m) << "_"
                        << util::AsStringCompact(p) << "_"
                        << ir << "_" 
                        << it;

#ifdef __show_init_mesg_debug__
                cout << "Exp/NuXSec >> Loading cross-section weight spline: " 
                     << splname.str() << " for Ereco: " << Ereco << ", Etrue: " << Etrue << endl;
#endif
                // Copy spline from input file
		TSpline3 * spl;
		if (s==kSKsFCFV1RingMuLike)   spl = (TSpline3 *) spline_file_1Rmu->Get(splname.str().c_str());
		else                          spl = (TSpline3 *) spline_file_1Re ->Get(splname.str().c_str()); 
	     
		if(!spl) {
		  if(!bin_is_not_empty) {
#ifdef __show_init_mesg_debug__
		    cout << "Exp/NuXSec >> No events or spline for: "
			 << splname.str()
			 << endl;
#endif
		    continue;
		  }
		  else {
#ifndef __use_flat_splines_for_missing_splines__
		    cerr << "Exp/NuXSec >> Null cross-section weight spline: " 
			 << splname.str() 
#if SPLINE_TYPE == 0
			 << ", " << mode[0]
#endif
			 << " Ereco: " << Ereco << ", Etrue: " << Etrue << endl;
		    exit(1);
#else
		    spl = flatspline;
#endif //__use_flat_splines_for_missing_splines__
		  }
		}
		else {
		  if(!bin_is_not_empty) {
		    cerr << "Exp/NuXSec >> Spline exists for empty bin: "
			 << splname.str() 
			 << " Ereco: " << Ereco << ", Etrue: " << Etrue
			 << ". Spline thinks there are " << spl->GetTitle()
			 << endl;
		    exit(1);
		  }
		}

#ifndef __skip_unexpected_spline_test_in_readsplines__
		string spln(spl->GetTitle());
		spln.erase(spln.end()-7,spln.end());
		int nevents_spl = atoi(spln.c_str());
		int nevents_hist = this->HowManyEvents(s, m, Ereco, Etrue);
		if(nevents_spl != nevents_hist) {
		  cerr << "Exp/NuXSec >> Spline and histogram disagree on number of events! "
		       << "Spline thinks there are: " << nevents_spl
		       << ", Histogram thinks there are: " << nevents_hist
		       << "; " << splname.str()
		       << "; Ereco: " << Ereco << ", Etrue: " << Etrue << endl;
		  exit(1);
		}
#endif	

                // Get index for current weight spline
                int idx = this->WghtSplineIdx(p,s,m,Ereco,Etrue);     
		if(idx < 0) {
		  cerr << "Exp/NuXSec >> WghtSplineIdx returned invalid index for spline: " 
                       << splname.str()
			 << endl;
		  exit(1);
		}

#ifdef __show_init_mesg_debug__
                cout << "Exp/NuXSec >> Storing cross-section weight spline at slot: " 
                     << idx << endl;
		cout << "Exp/NuXSec >> Spline thinks there are " 
		     << spl->GetTitle() << endl;
#endif
                if (s == kSKsFCFV1RingMuLike) fXSecWghtSpl1Rmu[idx] = new TSpline3(*spl);
		else                          fXSecWghtSpl1Re[idx] = new TSpline3(*spl);

#ifdef __show_init_mesg_debug2__
		//check that the spline is not flat
		cout << "Exp/NuXSec >> Spline weights at -1,0,+1 are: " ;
		  if (s == kSKsFCFV1RingMuLike){
		    cout << fXSecWghtSpl1Rmu[idx]->Eval(-1) << "\t"
			 << fXSecWghtSpl1Rmu[idx]->Eval(0) << "\t"
			 << fXSecWghtSpl1Rmu[idx]->Eval(+1) << endl;
		  }
		  else {
		    cout << fXSecWghtSpl1Re[idx]->Eval(-1) << "\t"
			 << fXSecWghtSpl1Re[idx]->Eval(0) << "\t"
			 << fXSecWghtSpl1Re[idx]->Eval(+1) << endl;
		  }
#endif
             }//Etrue
          }//Ereco
       }//p
    }//m
  }//s

#ifdef __show_init_mesg_info__
  cout << "Exp/NuXSec >> Done reading neutrino cross-section inputs (splines only)!" << endl;
#endif

  return true;  
}

//________________________________________________________________________________________
bool NuXSec::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/NuXSec >> Reading neutrino generator inputs (minus splines) from input XML file." << endl;
#endif

  //
  //
  //

  //book the bool templ histograms
  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     const int      n_ereco_bins    = util::NumXsecErecoBins(s);
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
  fSplineMAQERelErr =
    util::XML2Dbl(xml,root_node,"xsec/splines/maqe");
  fSpline1overMAQE2RelErr =
    util::XML2Dbl(xml,root_node,"xsec/splines/oneovermaqesq");
  fSplineMARESRelErr =
    util::XML2Dbl(xml,root_node,"xsec/splines/mares");

#ifdef __show_init_mesg_info__
  cout << "Exp/NuXSec >> Done reading neutrino cross-section inputs (not splines)!" << endl;
#endif

  return true;  
}
//________________________________________________________________________________________
void NuXSec::TellMeAboutTheBANFFfit(BANFF * banff)
{
  fBANFF = banff;
}
//________________________________________________________________________________________
double NuXSec::FinTuningWght(
  SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
  double tot_wght = 1.;

  for(int isyst = 0; isyst < kNumOfXSecSyst; isyst++) {
     PhysParamId_t p = kListOfXSecSyst[isyst];
     // "Tuning" involves reweighting to BANFF best-fit values.
     // Consider only the cross-section systematics that were included in the BANFF fit.
     bool is_banff_xsec_param = kXSecSystIncludedInBANFF[isyst];
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
double NuXSec::FinTuningWght_1p(
  PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
//
// This method provides weights to reweight MC to the BANFF best-fit parameters.
// It handles the 4 BANFF cross-section parameters which are correlated with the BANFF 
// flux parameters. In the current incarnation of the BANFF fit, these parameters are:
// - BANFF param, M_{A}  QE
// - BANFF param, M_{A}  RES
// - BANFF param, CCQE   normalization, 0 - 1.5 GeV (Etrue)
// - BANFF param, CCQE   normalization, 1.5 - 3.5 GeV (Etrue)
// - BANFF param, CCQE   normalization, 3.5 - 30.0 GeV (Etrue)
// - BANFF param, CC1pi  normalization, 0 - 2.5 GeV (Etrue)
// - BANFF param, CC1pi  normalization, 2.5 - 30.0 GeV (Etrue)
//

  // if( ! this->ApplySyst(p,s,m,Ereco,Etrue) ) return 1.; // check already made in function that calls this one

  double wght = 1.;

  //
  // BANFF param, M_{A} QE
  // BANFF param, M_{A} RES
  //
  if(p == kNuisBANFFMAQE || p == kNuisBANFFMARES)
  {
#if SPLINE_TYPE == -1
    return wght;
#endif
    // bool ev = this->IsThereAnEvent(s, m, Ereco, Etrue); //already done in SKMC
    //if(ev) {
    TSpline3 * wght_spl = this->FindWghtSpline(p,s,m,Ereco,Etrue);
    assert(wght_spl);
    double tb = fBANFF->BestFitValue(p);
    double relerr = 0;       
    if(p == kNuisBANFFMAQE)
#ifndef ONEOVERMAQE
      relerr = fSplineMAQERelErr;
#else
    relerr = fSpline1overMAQE2RelErr;
#endif
    else if(p == kNuisBANFFMARES)
      relerr = fSplineMARESRelErr;
    double numerator = (tb-1) / relerr;
    wght = wght_spl->Eval(numerator) / 1.;
    // }
  }
  
  //
  // BANFF param, CCQE   normalization, 0 - 1.5 GeV
  // BANFF param, CCQE   normalization, 1.5 - 3.5 GeV
  // BANFF param, CCQE   normalization, 3.5 - 30.0 GeV
  // BANFF param, CC1pi  normalization, 0 - 2.5 GeV
  // BANFF param, CC1pi  normalization, 2.5 - 30.0 GeV
  // BANFF param, NC1pi0 normalization
  //  
  else
  if(p == kNuisBANFFCCQEE0
     || p == kNuisBANFFCCQEE1
     || p == kNuisBANFFCCQEE2
     || p == kNuisBANFFCC1piE0
     || p == kNuisBANFFCC1piE1
     || p == kNuisBANFFNC1pi0
     )
  {
     wght = fBANFF->BestFitValue(p);
  }

#ifdef __show_banff_tuning_wght_mesg_debug__
  if (wght != 1.)
  {
     cout << "Exp/NuXSec/FinTuningWght >> "
          << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
          << ", Ereco = " << Ereco << " GeV" 
          << ", Etrue = " << Etrue << " GeV => Weight [" 
          << util::AsString(p) << "] = " << wght 
          << endl;   
  }
#endif

  if(TMath::Abs(wght) > conv::kMaxAllowedWght) {
    //get the ereco & etrue bin indexes (to be used to check whether the spline is sensible)
    int ireco = -1; 
    for (int i = 0; i < util::NumXsecErecoBins(s); i++) {
      if(util::XsecErecoBinEdges(s,i)   <= Ereco && 
	 util::XsecErecoBinEdges(s,i+1) >  Ereco) 
	{
	  ireco = i;
	  break;
	}
    }
    int itrue = -1; 
    int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
    int im = util::ArrIdx::Instance()->SuperkModeIdx(s, m);
    for (int i = 0; i < kXsNEtrueBins; i++) {
      if(kXsEtrueBinEdges[is][im][i]   <= Etrue && 
	 kXsEtrueBinEdges[is][im][i+1] >  Etrue) 
	{
	  itrue = i;
	  break;
	}
    }

    cout << "Exp/NuXSec/FinTuningWght >> Got an absurd weight! "
	 << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
	 << ", Ereco = " << Ereco << " GeV" 
	 << ", Etrue = " << Etrue << " GeV =>"
	 << " Bins are Ereco,Etrue = " << ireco << "," << itrue << " =>"
	 << " Weight [" << util::AsString(p) << "] = " << wght
	 << " Current value = " << this->CurrentValue(p)
	 << endl;   
    exit(1);
  }

  return wght;
}
//________________________________________________________________________________________
double NuXSec::SystWght(
  SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
  double tot_wght = 1.;

  for(int isyst = 0; isyst < kNumOfXSecSyst; isyst++) {
     PhysParamId_t p = kListOfXSecSyst[isyst];
     // Was tweaked?
     double dt = 0.;
     if(kXSecSystIncludedInBANFF[isyst])
       dt = TMath::Abs(fBANFF->CurrentValue(p) - fBANFF->BestFitValue(p)); // banff parameters
     else
       dt = TMath::Abs(this->CurrentValue(p)); // non-banff normalisations
     if(dt < kASmallNum) continue;
     
     // Is the current systematic parameter applied to current sample/mode end energy bin?
     bool apply = this->ApplySyst(p,s,m,Ereco,Etrue);
     if(!apply) continue;
     // Calculate weight
     double wght = this->SystWght_1p(p,s,m,Ereco,Etrue);
     tot_wght *= wght;
  }

  return tot_wght;
}
//________________________________________________________________________________________
double NuXSec::SystWght_1p(
  PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
//
// This method provides weights to account for systematic variations to the XSec params.
// The weights are relative to MC reweighted to the BANFF best-fit point.
// It handles both the BANFF and uncorelated XSec params.
//

  // if( ! this->ApplySyst(p,s,m,Ereco,Etrue) ) return 1.; // check already made in function that calls this one

  double wght = 1.;

  //
  // BANFF param # 31, M_{A} QE
  // BANFF param # 32, M_{A} RES
  //
  if(p == kNuisBANFFMAQE || p == kNuisBANFFMARES)
  {
#if SPLINE_TYPE == -1
    return wght;
#endif
     //bool ev = this->IsThereAnEvent(s, m, Ereco, Etrue);//already done in SKMC
     //if(ev) {
       TSpline3 * wght_spl = this->FindWghtSpline(p,s,m,Ereco,Etrue);
       assert(wght_spl);
       double tb = fBANFF->BestFitValue(p);
       double tc = fBANFF->CurrentValue(p);   
       double relerr = 0;
       if(p == kNuisBANFFMAQE)
#ifndef ONEOVERMAQE
	 relerr = fSplineMAQERelErr;
#else
	 relerr = fSpline1overMAQE2RelErr;
#endif
       else if(p == kNuisBANFFMARES)
	 relerr = fSplineMARESRelErr;
       double denominator = (tb-1) / relerr;
       double numerator = (tc-1) / relerr;
       wght = wght_spl->Eval(numerator) / wght_spl->Eval(denominator);
       //}
  }

  //
  //
  // BANFF param, CCQE   normalization, 0 - 1.5 GeV
  // BANFF param, CCQE   normalization, 1.5 - 3.5 GeV
  // BANFF param, CCQE   normalization, 3.5 - 30.0 GeV
  // BANFF param, CC1pi  normalization, 0 - 2.5 GeV
  // BANFF param, CC1pi  normalization, 2.5 - 30.0 GeV
  // BANFF param, NC1pi0 normalization
  // 
  else
  if(p == kNuisBANFFCCQEE0
     || p == kNuisBANFFCCQEE1
     || p == kNuisBANFFCCQEE2
     || p == kNuisBANFFCC1piE0
     || p == kNuisBANFFCC1piE1
     || p == kNuisBANFFNC1pi0
     )
  {
     wght = fBANFF->CurrentValue(p) / fBANFF->BestFitValue(p);
  }

  //
  // Fermi momentum
  // spectral function on/off
  // W-shape
  // CC other shape
  //
  else
  if(p == kNuisXSecFermiMom || p == kNuisXSecSpecFunc || p == kNuisXSecWshape || p == kNuisXSecShapeCCoth ||
     p == kNuisXSecPiLessDelta || p == kNuisXSecBindE)
  {
#if SPLINE_TYPE == -1
    return wght;
#endif
    //    bool ev = this->IsThereAnEvent(s, m, Ereco, Etrue);//already done in SKMC
    //if(ev) {
      TSpline3 * wght_spl = this->FindWghtSpline(p,s,m,Ereco,Etrue);
      assert(wght_spl);
      double tc = this->CurrentValue(p);
      wght = wght_spl->Eval(tc) / 1.; //wght_spl->Eval(0); 
      //}
  }

  //
  // CCQE norm 1.5 - 3.5 GeV
  // CCQE norm > 3.5 GeV
  // CC1pi norm > 2.5 GeV   
  // CC coherent norm
  // NC 1piC norm
  // NC coherent norm
  // NC other norm
  // nue / numu cross section
  //
  else
  if(
     p == kNuisXSecNormCCcoh || p == kNuisXSecNormNC1piC || p == kNuisXSecNormNCcoh ||
     p == kNuisXSecNormNCoth || p == kNuisXSecNormNuEToNuMu || p == kNuisXSecNormAntiNuToNu)
  {
    wght = 1 + (this->CurrentValue(p) * this->RelativeError(p));
  }


#ifdef __show_syst_var_wght_mesg_debug__
  if (wght != 1.)
  {
     cout << "Exp/NuXSec/SystWght >> "
          << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
          << ", Ereco = " << Ereco << " GeV" 
          << ", Etrue = " << Etrue << " GeV => Weight [" 
          << util::AsString(p) << "] = " << wght 
          << endl;   
  }
#endif

  if(TMath::Abs(wght) > conv::kMaxAllowedWght) {
    //get the ereco & etrue bin indexes (to be used to check whether the spline is sensible)
    int ireco = -1; 
    for (int i = 0; i < util::NumXsecErecoBins(s); i++) {
      if(util::XsecErecoBinEdges(s,i)   <= Ereco && 
	 util::XsecErecoBinEdges(s,i+1) >  Ereco) 
	{
	  ireco = i;
	  break;
	}
    }
    int itrue = -1; 
    int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
    int im = util::ArrIdx::Instance()->SuperkModeIdx(s, m);
    for (int i = 0; i < kXsNEtrueBins; i++) {
      if(kXsEtrueBinEdges[is][im][i]   <= Etrue && 
	 kXsEtrueBinEdges[is][im][i+1] >  Etrue) 
	{
	  itrue = i;
	  break;
	}
    }
    cout << "Exp/NuXSec/SystWght >> Got an absurd weight! "
	 << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
	 << ", Ereco = " << Ereco << " GeV" 
	 << ", Etrue = " << Etrue << " GeV =>"
	 << " Bins are Ereco,Etrue = " << ireco << "," << itrue << " =>"
	 << " Weight [" << util::AsString(p) << "] = " << wght
	 << " Current value = " << this->CurrentValue(p)
	 << endl;
    exit(1);
  }
  
  return wght;
}
//________________________________________________________________________________________
double NuXSec::CurrentValue(PhysParamId_t p) const
{
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
    case kNuisBANFFMAQE:
    case kNuisBANFFMARES:
    case kNuisBANFFCCQEE0:
    case kNuisBANFFCCQEE1:
    case kNuisBANFFCCQEE2:
    case kNuisBANFFCC1piE0:
    case kNuisBANFFCC1piE1:
    case kNuisBANFFNC1pi0:
      return fBANFF->CurrentValue(p);
      break;
    default:
      return 0;
  }
  return 0;
}
//________________________________________________________________________________________
double NuXSec::RelativeError(PhysParamId_t p) const
{
  switch(p) {
    case kNuisXSecNormCCcoh:
      return fNormCCcohRelErr;
      break;
    case kNuisXSecNormNC1piC:
      return fNormNC1piCRelErr;
      break;
    case kNuisXSecNormNCcoh:
      return fNormNCcohRelErr;
      break;
    case kNuisXSecNormNCoth:
      return fNormNCothRelErr;
      break;
    case kNuisXSecNormNuEToNuMu:
      return fNormNuEToNuMuRelErr;
      break;
    case kNuisXSecNormAntiNuToNu:
      return fNormAntiNuToNuRelErr;
      break;
    case kNuisBANFFMAQE:
#ifndef ONEOVERMAQE
	 return fSplineMAQERelErr;
#else
	 return fSpline1overMAQE2RelErr;
#endif
      break;
    case kNuisBANFFMARES:
      return fSplineMARESRelErr;
      break;
    default:
      cout << "Relative error is not defined for xsec parameter: " << AsString(p) << endl;
      return 0;
  }
  return 0;
}
//________________________________________________________________________________________
void NuXSec::PrintSystematicVariations(ostream & stream) const
{
  stream << "Exp/XSec >> XSec parameters:"                              << endl; 
  stream << "         >> Fermi momentum      " << fFermiMomTwkDial      << endl;
  stream << "         >> Spectral function   " << fSpecFuncTwkDial      << endl;
  stream << "         >> W-shape             " << fWshapeTwkDial        << endl;
  stream << "         >> CC-oth shape        " << fShapeCCothTwkDial    << endl;
  stream << "         >> Pi-less delta decay " << fPiLessDeltaTwkDial   << endl;
  stream << "         >> Binding energy      " << fBindETwkDial         << endl;
  stream << "         >> CC-coh norm         " << fNormCCcohTwkDial     << endl;
  stream << "         >> NC-1piC norm        " << fNormNC1piCTwkDial    << endl;
  stream << "         >> NC-coh norm         " << fNormNCcohTwkDial     << endl;
  stream << "         >> NC-oth norm         " << fNormNCothTwkDial     << endl;
  stream << "         >> NuE to NuMu norm    " << fNormNuEToNuMuTwkDial << endl;
  stream << "         >> AntiNu to Nu norm   " << fNormAntiNuToNuTwkDial<< endl;
}
//________________________________________________________________________________________

//
// Cross-section weight splines
//

//________________________________________________________________________________________
TSpline3 * NuXSec::FindWghtSpline(
   PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
  int idx = this->WghtSplineIdx(p,s,m,Ereco,Etrue);
  if(idx == -1 || idx == -2) 
  {
    cerr << "Exp/NuXSec/FindWghtSpline >> Didn't find weight spline for "
         << "sample = "      << util::AsString(s) 
         << ", mode = "      << util::AsString(m)
         << ", parameter = " << util::AsString(p)
         << ", Ereco = "     << Ereco << " GeV"
         << ", Etrue = "     << Etrue << " GeV!"
         << endl;
    exit(1);
  }

  if (s == kSKsFCFV1RingMuLike) return fXSecWghtSpl1Rmu[idx];
  else return fXSecWghtSpl1Re[idx];
}
//________________________________________________________________________________________
int NuXSec::WghtSplineIdx(
   PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
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

  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(is == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
    cout << "Exp/NuXSec/WghtSplineIdx >> Didn't find Idx for sample "
	 << util::AsString(s)
	 << ". util::ArrIdx returned "
	 << util::ArrIdx::Instance()->SuperkSampleIdx(s)
	 << endl;
#endif
    return -1;
  }
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(im == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
    cout << "Exp/NuXSec/WghtSplineIdx >> Didn't find Idx for mode "
	 << util::AsString(m)
	 << " in sample "
	 << util::AsString(s)
	 << ". util::ArrIdx returned "
	 << util::ArrIdx::Instance()->SuperkModeIdx(s,m)
	 << endl;
#endif
    return -1;
  }
  int ip = util::ArrIdx::Instance()->XSecSystIdx(p);
  if(ip == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
    cout << "Exp/NuXSec/WghtSplineIdx >> Didn't find Idx for xsec sysy "
	 << util::AsString(p)
	 << ". util::ArrIdx returned "
	 << util::ArrIdx::Instance()->XSecSystIdx(p)
	 << endl;
#endif
    return -1;
  }
  int ireco = -1; 
  for (int i = 0; i < util::NumXsecErecoBins(s); i++) {
    if(util::XsecErecoBinEdges(s,i)   <= Ereco && 
       util::XsecErecoBinEdges(s,i+1) >  Ereco) 
    {
        ireco = i;
        break;
    }
  }
  if(ireco == -1) {
#ifdef __show_wght_spline_idx_mesg_debug__
    cout << "Exp/NuXSec/WghtSplineIdx >> Didn't find ERecoBin for ereco value "
	 << Ereco << endl;
#endif
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
    cout << "Exp/NuXSec/WghtSplineIdx >> Didn't find ETrueBin for etrue value "
	 << Ereco 
	 << " in mode "
	 << util::AsString(m)
	 << " in sample "
	 << util::AsString(s)
	 << endl;
      }
#endif
    return -2;
  }

  int idx = itrue 
          + ireco * kXsNEtrueBins 
	  + ip    * kXsNEtrueBins * util::NumXsecErecoBins(s)
	  + im    * kXsNEtrueBins * util::NumXsecErecoBins(s) * kNumOfXSecSyst;
	    //	  + is    * kXsNEtrueBins * util::NumXsecErecoBins(s) * kNumOfXSecSyst * kNSuperkModes;

  return idx;
}
//________________________________________________________________________________________
bool NuXSec::ApplySyst(
  PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
  bool apply2mode = this->ApplySystToSKMode(p,s,m);
  bool apply2bin  = this->ApplySystToEnergyBin(p,Ereco,Etrue);
  bool apply = apply2mode && apply2bin;
  return apply;
}
//________________________________________________________________________________________
bool NuXSec::ApplySystToSKMode(PhysParamId_t p, SuperkSample_t s, SuperkMode_t m) const
{
  int isyst = util::ArrIdx::Instance()->XSecSystIdx(p);
  if (isyst == -1) return false;
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if (is == -1) return false;
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if (im == -1) return false;

  bool apply = kXsApplyXSecSystToSKMode[is][im][isyst];

  return apply;
}
//________________________________________________________________________________________
bool NuXSec::ApplySystToEnergyBin(PhysParamId_t p, double Ereco, double Etrue) const
{
  int isyst = util::ArrIdx::Instance()->XSecSystIdx(p);
  if (isyst == -1) return false;

  bool apply = 
     (Ereco >= kXSecSystErecoCut[isyst][0] && Ereco < kXSecSystErecoCut[isyst][1]) &&
     (Etrue >= kXSecSystEtrueCut[isyst][0] && Etrue < kXSecSystEtrueCut[isyst][1]);

  return apply;
}
//________________________________________________________________________________________
void NuXSec::BuildXsBinningMCTemplates(
  TH2D * input_mctempl_array[kns][knm])
{
  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     for(int im = 0; im < knm; im++) {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue; 
#ifdef __show_build_xsbinning_template_mesg_info__
       cout << "Exp/NuXSec >> Building Xsec binning MC template for "
            << "sample = " << util::AsString(s) 
            << ", mode = " << util::AsString(m)
            << endl;
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
	     cout << std::fixed
		  << "Exp/NuXSec >> Filled Xsec binning MC template for "
		  << "ereco = " << mctempl->GetBinCenter(ir)
		  << " etrue = " << mctempl->GetBinCenter(it)
		  << endl
		  << "\t Put into bins: ereco " 
		  << mctempl_xsbinning->GetXaxis()->GetBinLowEdge(bx) << " - " 
		  << mctempl_xsbinning->GetXaxis()->GetBinUpEdge(bx) << "\t etrue " 
		  << mctempl_xsbinning->GetYaxis()->GetBinLowEdge(by) << " - " 
		  << mctempl_xsbinning->GetYaxis()->GetBinUpEdge(by)
		  << endl;
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
bool NuXSec::IsThereAnEvent(SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  assert(is != -1 && im != -1);

  TH2D * mctempl_xsbinning = fXsBinningTempl[is][im];
  assert(mctempl_xsbinning);

  int ir = mctempl_xsbinning->GetXaxis()->FindBin(Ereco);
  int it = mctempl_xsbinning->GetYaxis()->FindBin(Etrue);

  bool ev = mctempl_xsbinning->GetBinContent(ir,it) > kASmallNum;

#ifdef __show_isanevent_mesg_debug__
  if(ev)//&& (s == conv::kSKsFCFV1RingELike))
    cout << "Exp/NuXSec >> There are " << mctempl_xsbinning->GetBinContent(ir,it)
	 << " events for sample = " << util::AsString(s) 
	 << ", mode = " << util::AsString(m)
	 << ", Etrue = " << Etrue
	 << ", Ereco = " << Ereco
	 << endl;
#endif

  return ev;
} 
//________________________________________________________________________________________
int NuXSec::HowManyEvents(SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  int im = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  assert(is != -1 && im != -1);

  TH2D * mctempl_xsbinning = fXsBinningTempl[is][im];
  assert(mctempl_xsbinning);

  int ir = mctempl_xsbinning->GetXaxis()->FindBin(Ereco);
  int it = mctempl_xsbinning->GetYaxis()->FindBin(Etrue);

  int ev = (int)mctempl_xsbinning->GetBinContent(ir,it);

#ifdef __show_isanevent_mesg_debug__
  if(ev)
    cout << "Exp/NuXSec >> There are " << ev
	 << " events for sample = " << util::AsString(s) 
	 << ", mode = " << util::AsString(m)
	 << ", Etrue = " << Etrue
	 << ", Ereco = " << Ereco
	 << endl;
#endif

  return ev;
} 
//________________________________________________________________________________________
TH2D * NuXSec::Template(SuperkSample_t s, SuperkMode_t m) const
{
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }

  TH2D * h = 0;
  h = fXsBinningTempl[sidx][midx];
  assert(h);

  return h;
}
//________________________________________________________________________________________
double NuXSec::ChisqPenalty(bool all_syst) const
{
  //only contains the contributions from the non-BANFF xsec parameters

  double chisq_penalty = 0;

  // when only fitting for most important systematics
  if(!all_syst) {
    if (fPiLessDeltaTwkDial)  chisq_penalty += fPiLessDeltaTwkDial*fPiLessDeltaTwkDial;
    if (fSpecFuncTwkDial)     chisq_penalty += fSpecFuncTwkDial*fSpecFuncTwkDial;
    if (fBindETwkDial)        chisq_penalty += fBindETwkDial*fBindETwkDial;
    if (fWshapeTwkDial)       chisq_penalty += fWshapeTwkDial*fWshapeTwkDial;
    return chisq_penalty;
  }

  if (fFermiMomTwkDial)       chisq_penalty += fFermiMomTwkDial*fFermiMomTwkDial;
  if (fSpecFuncTwkDial)       chisq_penalty += fSpecFuncTwkDial*fSpecFuncTwkDial;
  if (fWshapeTwkDial)         chisq_penalty += fWshapeTwkDial*fWshapeTwkDial;
  if (fShapeCCothTwkDial)     chisq_penalty += fShapeCCothTwkDial*fShapeCCothTwkDial;
  if (fPiLessDeltaTwkDial)    chisq_penalty += fPiLessDeltaTwkDial*fPiLessDeltaTwkDial;
  if (fBindETwkDial)          chisq_penalty += fBindETwkDial*fBindETwkDial;
  if (fNormCCcohTwkDial)      chisq_penalty += fNormCCcohTwkDial*fNormCCcohTwkDial;
  if (fNormNC1piCTwkDial)     chisq_penalty += fNormNC1piCTwkDial*fNormNC1piCTwkDial;
  if (fNormNCcohTwkDial)      chisq_penalty += fNormNCcohTwkDial*fNormNCcohTwkDial;
  if (fNormNCothTwkDial)      chisq_penalty += fNormNCothTwkDial*fNormNCothTwkDial;
  if (fNormNuEToNuMuTwkDial)  chisq_penalty += fNormNuEToNuMuTwkDial*fNormNuEToNuMuTwkDial;
  if (fNormAntiNuToNuTwkDial) chisq_penalty += fNormAntiNuToNuTwkDial*fNormAntiNuToNuTwkDial;

  return chisq_penalty;
}
//________________________________________________________________________________________
