//________________________________________________________________________________________
/*

 >> ChangeLog:

 For changes up to version 2010_v3.0.0 (last implementation of the 2010-era analysis)
 see code in the corresponding SVN branch.

 Changes since 2011_v0.9.0:

 @

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__
//#define __show_init_mesg_debug__
//#define __show_init_mesg_debug2__
//#define __show_banff_tuning_wght_mesg_info__
//#define __show_banff_tuning_wght_mesg_debug__
//#define __show_syst_var_wght_mesg_debug__
//#define __show_build_bool_template_mesg_info__
//#define __show_wght_spline_idx_mesg_debug__
//#define __show_isanevent_mesg_debug__
//#define __save_ori_templates_init_debug__

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

  for(int ispl = 0; ispl < knwspl; ispl++) {
     fXSecWghtSpl[ispl] = 0;
  }

  for(int is=0; is<kns; is++) {
    for(int im=0; im<knm; im++) {
      fXsBinningTempl [is][im] = 0;
    }
  }

  fFermiMomTwkDial      = 0;
  fSpecFuncTwkDial      = 0;
  fWshapeTwkDial        = 0;
  fShapeCCothTwkDial    = 0;
  fPiLessDeltaTwkDial   = 0;
  fENuShapeTwkDial      = 0;
  fNormCCQEE1TwkDial    = 0;
  fNormCCQEE2TwkDial    = 0;
  fNormCC1piE1TwkDial   = 0;
  fNormCCcohTwkDial     = 0;
  fNormNCcohTwkDial     = 0;
  fNormNCothTwkDial     = 0;
  fNormNuEToNuMuTwkDial = 0;

}
//________________________________________________________________________________________
void NuXSec::CleanUp(void)
{
  for(int ispl = 0; ispl < knwspl; ispl++) {
     if (fXSecWghtSpl[ispl]) delete fXSecWghtSpl[ispl];
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

  TFile * spline_file = util::XML2TFile(xml,root_node,"xsec/weight_splines");
  assert(spline_file);

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
          for(int ir = 0; ir < kXsNErecoBins; ir++) {
             double Ereco = kXsErecoBinEdges[ir] + kOneMeV;
             for(int it = 0; it < kXsNEtrueBins; it++) {
                double Etrue = kXsEtrueBinEdges[is][im][it] + kOneMeV;

                bool bin_is_beyond_end_of_range = (kXsEtrueBinEdges[is][im][it+1] < 0.);
                if(bin_is_beyond_end_of_range) continue;

                // Build spline name
                ostringstream splname;
                splname << "spl_"
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
                TSpline3 * spl = (TSpline3 *) spline_file->Get(splname.str().c_str());
		bool bin_is_not_empty = this->IsThereAnEvent(s, m, Ereco, Etrue);
	     
		if(!spl) {
		  if(!bin_is_not_empty) {
#ifdef __show_init_mesg_debug__
		    cout << "Exp/NuXSec >> No events or spline for: "
			 << splname.str() << endl;
#endif
		    continue;
		  }
		  else {
		    cerr << "Exp/NuXSec >> Null cross-sectionweight spline: " 
			 << splname.str() << " Ereco: " << Ereco << ", Etrue: " << Etrue << endl;
		    exit(1); 
		  }
		}
		else {
		  if(!bin_is_not_empty) {
		    cerr << "Exp/NuXSec >> Spline exists for empty bin: "
			 << splname.str() << " Ereco: " << Ereco << ", Etrue: " << Etrue
			 << ". Spline thinks there are " << spl->GetTitle() << endl;
		    exit(1);
		  }
		}

                // Get index for current weight spline
                int idx = this->WghtSplineIdx(p,s,m,Ereco,Etrue);     
//		if(idx == -2) {
//#ifdef __show_init_mesg_debug__
//                cout << "Exp/NuXSec >> Skipping spline due to invalid Etrue bin for Etrue: " 
//                     << Etrue << endl;
//#endif
//		  continue;
//		}
		if(idx < 0) {
		  cerr << "Exp/NuXSec >> WghtSplineIdx returned invalid index for spline: " 
                       << splname.str();
		  exit(1);
		}

#ifdef __show_init_mesg_debug__
                cout << "Exp/NuXSec >> Storing cross-section weight spline at slot: " 
                     << idx << endl;
#endif
                fXSecWghtSpl[idx] = new TSpline3(*spl);

#ifdef __show_init_mesg_debug2__
		//check that the spline is not flat
		cout << "Exp/NuXSec >> Spline weights at -1,0,+1 are: "
		     << fXSecWghtSpl[idx]->Eval(-1) << "\t"
		     << fXSecWghtSpl[idx]->Eval(0) << "\t"
		     << fXSecWghtSpl[idx]->Eval(+1) << endl;
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
     const int      n_ereco_bins    = util::NumErecoBins (s);
     const double * ereco_bin_edges = util::ErecoBinEdges(s);
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

  fNormCCQEE1RelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/ccqee1");
  fNormCCQEE2RelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/ccqee2");
  fNormCC1piE1RelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/cc1pie1");
  fNormCCcohRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/cccoh");
  fNormNCcohRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/nccoh");
  fNormNCothRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/ncoth");
  fNormNuEToNuMuRelErr =
    util::XML2Dbl(xml,root_node,"xsec/normalisations/numu_to_nue");

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
// - BANFF param # 22, M_{A}  QE
// - BANFF param # 23, M_{A}  RES
// - BANFF param # 24, CCQE   normalization, 0 - 1.5 GeV (Etrue)
// - BANFF param # 25, CC1pi  normalization, 0 - 2.5 GeV (Etrue)
// - BANFF param # 26, NC1pi0 normalization
//

  // if( ! this->ApplySyst(p,s,m,Ereco,Etrue) ) return 1.; // check already made in function that calls this one

  double wght = 1.;

  //
  // BANFF param # 22, M_{A} QE
  // BANFF param # 23, M_{A} RES
  //
  if(p == kNuisBANFF22 || p == kNuisBANFF23)
  {
     bool ev = this->IsThereAnEvent(s, m, Ereco, Etrue);
     if(ev) {
       TSpline3 * wght_spl = this->FindWghtSpline(p,s,m,Ereco,Etrue);
       assert(wght_spl);
       double t = fBANFF->BestFitValue(p);
       wght = wght_spl->Eval(t);
       //wght = wght_spl->Eval(0);
       wght = TMath::Max(0.,wght);
     }
  }

  //
  // BANFF param # 24, CCQE  normalization, 0 - 1.5 GeV
  // BANFF param # 25, CCRES normalization, 0 - 2.5 GeV
  //  
  else
  if(p == kNuisBANFF24 || p == kNuisBANFF25 || p == kNuisBANFF26)
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

  if(wght < 0 || wght > conv::kMaxAllowedWght) {
     cout << "Exp/NuXSec/FinTuningWght >> Got an absurd weight! "
          << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
          << ", Ereco = " << Ereco << " GeV" 
          << ", Etrue = " << Etrue << " GeV => Weight [" 
          << util::AsString(p) << "] = " << wght 
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
// This method provides weights to account for systematic variations to the BANFF params.
// The weights are relative to MC reweighted to the BANFF best-fit point.
// It handles the 4 BANFF cross-section parameters which are correlated with the BANFF 
// flux parameters. In the current incarnation of the BANFF fit, these parameters are:
// - BANFF param # 22, M_{A} QE
// - BANFF param # 23, M_{A} RES
// - BANFF param # 24, CCQE  normalization, 0 - 1.5 GeV (Etrue)
// - BANFF param # 25, CCRES normalization, 0 - 2.5 GeV (Etrue)
//

  // if( ! this->ApplySyst(p,s,m,Ereco,Etrue) ) return 1.; // check already made in function that calls this one

  double wght = 1.;

  //
  // BANFF param # 22, M_{A} QE
  // BANFF param # 23, M_{A} RES
  //
  if(p == kNuisBANFF22 || p == kNuisBANFF23)
  {
     bool ev = this->IsThereAnEvent(s, m, Ereco, Etrue);
     if(ev) {
       TSpline3 * wght_spl = this->FindWghtSpline(p,s,m,Ereco,Etrue);
       assert(wght_spl);
       double tb = fBANFF->BestFitValue(p);
       double tc = fBANFF->CurrentValue(p);
       double diag = 0;
       if(p == kNuisBANFF22)
	 diag = fBANFF->CovMatrix()[22][22];
       else if(p == kNuisBANFF23)
	 diag = fBANFF->CovMatrix()[23][23];
       wght = wght_spl->Eval(((tc / tb) - 1) / diag) / wght_spl->Eval(tb);
       wght = TMath::Max(0.,wght);
     }
  }

  //
  // BANFF param # 24, CCQE  normalization, 0 - 1.5 GeV
  // BANFF param # 25, CCRES normalization, 0 - 2.5 GeV
  //  
  else
  if(p == kNuisBANFF24 || p == kNuisBANFF25 || p == kNuisBANFF26)
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
     p == kNuisXSecPiLessDelta || p == kNuisXSecENuShape)
  {
    bool ev = this->IsThereAnEvent(s, m, Ereco, Etrue);
    if(ev) {
      TSpline3 * wght_spl = this->FindWghtSpline(p,s,m,Ereco,Etrue);
      assert(wght_spl);
      double tc = this->CurrentValue(p);
      wght = wght_spl->Eval(tc) / 1.; //wght_spl->Eval(0);
      wght = TMath::Max(0.,wght);
    }
  }

  //
  // CCQE norm 1.5 - 3.5 GeV
  // CCQE norm > 3.5 GeV
  // CC1pi norm > 2.5 GeV   
  // CC coherent norm
  // NC coherent norm
  // NC other norm
  // nue / numu cross section
  //
  else
  if(p == kNuisXSecNormCCQEE1 || p == kNuisXSecNormCCQEE2 || p == kNuisXSecNormCC1piE1 || p == kNuisXSecNormCCcoh || 
     p == kNuisXSecNormNCcoh || p == kNuisXSecNormNCoth || p == kNuisXSecNormNuEToNuMu)
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

  if(wght < 0 || wght > conv::kMaxAllowedWght) {
     cout << "Exp/NuXSec/SystWght >> Got an absurd weight! "
          << "Sample/Mode = " << util::AsString(s) << " / " << util::AsString(m)
          << ", Ereco = " << Ereco << " GeV" 
          << ", Etrue = " << Etrue << " GeV => Weight [" 
          << util::AsString(p) << "] = " << wght 
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
    case kNuisXSecENuShape:
      return fENuShapeTwkDial;
      break;
    case kNuisXSecNormCCQEE1:
      return fNormCCQEE1TwkDial;
      break;
    case kNuisXSecNormCCQEE2:
      return fNormCCQEE2TwkDial;
      break;
    case kNuisXSecNormCC1piE1:
      return fNormCC1piE1TwkDial;
      break;
    case kNuisXSecNormCCcoh:
      return fNormCCcohTwkDial;
      break;
    case kNuisXSecNormNCcoh:
      return fNormNCcohTwkDial;
      break;
    case kNuisXSecNormNCoth:
      return fNormNCothTwkDial;
      break;
    case kNuisXSecNormNuEToNuMu:
      return fNormNuEToNuMuTwkDial;
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
    case kNuisXSecNormCCQEE1:
      return fNormCCQEE1RelErr;
      break;
    case kNuisXSecNormCCQEE2:
      return fNormCCQEE2RelErr;
      break;
    case kNuisXSecNormCC1piE1:
      return fNormCC1piE1RelErr;
      break;
    case kNuisXSecNormCCcoh:
      return fNormCCcohRelErr;
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
    default:
      cout << "Relative error is not defined for xsec parameter: " << AsString(p) << endl;
      return 0;
  }
  return 0;
}
//________________________________________________________________________________________
void NuXSec::PrintSystematicVariations(void) const
{
  cout << "Exp/XSec >> XSec parameters:"                            << endl; 
  cout << "         >> Fermi momentum      " << fFermiMomTwkDial      << endl;
  cout << "         >> Spectral function   " << fSpecFuncTwkDial      << endl;
  cout << "         >> W-shape             " << fWshapeTwkDial        << endl;
  cout << "         >> CC-oth shape        " << fShapeCCothTwkDial    << endl;
  cout << "         >> Pi-less delta decay " << fPiLessDeltaTwkDial   << endl;
  cout << "         >> CC1pi ENu Shape     " << fENuShapeTwkDial      << endl;
  cout << "         >> CCQE bin1 norm      " << fNormCCQEE1TwkDial    << endl;
  cout << "         >> CCQE bin2 norm      " << fNormCCQEE2TwkDial    << endl;
  cout << "         >> CC1pi bin1 norm     " << fNormCC1piE1TwkDial   << endl;
  cout << "         >> CC-coh norm         " << fNormCCcohTwkDial     << endl;
  cout << "         >> NC-coh norm         " << fNormNCcohTwkDial     << endl;
  cout << "         >> NC-oth norm         " << fNormNCothTwkDial     << endl;
  cout << "         >> NuE to NuMu norm    " << fNormNuEToNuMuTwkDial << endl;
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

  return fXSecWghtSpl[idx];
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
  for (int i = 0; i < kXsNErecoBins; i++) {
    if(kXsErecoBinEdges[i]   <= Ereco && 
       kXsErecoBinEdges[i+1] >  Ereco) 
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
          + ip    * kXsNEtrueBins * kXsNErecoBins
          + im    * kXsNEtrueBins * kXsNErecoBins * kNumOfXSecSyst
          + is    * kXsNEtrueBins * kXsNErecoBins * kNumOfXSecSyst * kNSuperkModes;

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

       for(int ir = 0; ir <= kNErecoBins1Rmu; ir++) {
	 for(int it = 0; it <= kNEtrueBins; it++) {
	   if(mctempl->GetBinContent(ir,it) > kASmallNum)
	     mctempl_xsbinning->Fill(mctempl->GetBinCenter(ir), mctempl->GetBinCenter(it), mctempl->GetBinContent(ir,it));
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
  if(ev)
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
