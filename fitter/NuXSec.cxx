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
//#define __show_banff_tuning_wght_mesg_debug__
//#define __show_syst_var_wght_mesg_debug__
#define __show_fsi_cholesky_mesg_info__
#define __show_fsi_cholesky_mesg_debug__
#define __show_fsi_tweak_mesg_debug__

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

  fFSICurrentValues = 0;
  fFSICovMtrx       = 0;
  fFSICovMtrxChoL   = 0;

  for(int ispl = 0; ispl < knwspl; ispl++) {
     fXSecWghtSpl[ispl] = 0;
  }
}
//________________________________________________________________________________________
void NuXSec::CleanUp(void)
{
  if ( fFSICurrentValues ) delete fFSICurrentValues;
  if ( fFSICovMtrx       ) delete fFSICovMtrx;
  if ( fFSICovMtrxChoL   ) delete fFSICovMtrxChoL;

  for(int ispl = 0; ispl < knwspl; ispl++) {
     if (fXSecWghtSpl[ispl]) delete fXSecWghtSpl[ispl];
  }
}
//________________________________________________________________________________________
bool NuXSec::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/NuXSec >> Reading neutrino generator inputs from input XML file." << endl;
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
                     << splname.str() << endl;
#endif
                // Copy spline from input file
                TSpline3 * spl = (TSpline3 *) spline_file->Get(splname.str().c_str());
                if(!spl) {
                   cerr << "Exp/NuXSec >> Null cross-sectionweight spline: " 
                        << splname.str() << endl;
                   exit(1); 
                }
                // Get index for current weight spline
                int idx = this->WghtSplineIdx(p,s,m,Ereco,Etrue);     
#ifdef __show_init_mesg_debug__
                cout << "Exp/NuXSec >> Storing cross-section weight spline at slot: " 
                     << idx << endl;
#endif
                fXSecWghtSpl[idx] = new TSpline3(*spl);
             }//Etrue
          }//Ereco
       }//p
    }//m
  }//s

/*
  //
  // Load info needed for handling FSI systematics
  //

  const int np = conv::kNumOfFSIParams;

  // Get BANFF fit covariance matrix
  fFSICovMtrx = util::XML2TMatrixD(xml,root_node,"xsec/fsi/covariance_matrix");
  if(!fFSICovMtrx) {
    cerr << "Exp/NuXSec >> No BANFF covariance matrix was found." << endl;
    exit(1);
  }
#ifdef __show_init_mesg_info__
  cout << "Exp/NuXSec >> Input FSI covariance matrix: " << endl;
  util::PrintCompact(*fFSICovMtrx,7,4);
#endif

  // check cov matrix dimension against no
  // TODO

  // Create/init current FSI params
  fFSICurrentValues = new TVectorD(np);
  for(int ip = 0; ip < np; ip++) 
  {
     (*fFSICurrentValues)[ip] = 1.;
  }
*/

#ifdef __show_init_mesg_info__
  cout << "Exp/NuXSec >> Done reading neutrino cross-section inputs!" << endl;
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
  PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue) const
{
//
// This method provides weights to reweight MC to the BANFF best-fit parameters.
// It handles the 4 BANFF cross-section parameters which are correlated with the BANFF 
// flux parameters. In the current incarnation of the BANFF fit, these parameters are:
// - BANFF param # 22, M_{A} QE
// - BANFF param # 23, M_{A} RES
// - BANFF param # 24, CCQE  normalization, 0 - 1.5 GeV (Etrue)
// - BANFF param # 25, CCRES normalization, 0 - 2.5 GeV (Etrue)
//

  if( ! this->ApplySyst(p,s,m,Ereco,Etrue) ) return 1.;

  double wght = 1.;

  //
  // BANFF param # 22, M_{A} QE
  // BANFF param # 23, M_{A} RES
  //
  if(p == kNuisBANFF22 || p == kNuisBANFF23)
  {
     TSpline3 * wght_spl = this->FindWghtSpline(p,s,m,Ereco,Etrue);
     assert(wght_spl);
     double t = fBANFF->BestFitValue(p);
     wght = wght_spl->Eval(t);
     wght = TMath::Max(0.,wght); 
  }

  //
  // BANFF param # 24, CCQE  normalization, 0 - 1.5 GeV
  // BANFF param # 25, CCRES normalization, 0 - 2.5 GeV
  //  
  else
  if(p == kNuisBANFF24 || p == kNuisBANFF25)
  {
     wght = fBANFF->BestFitValue(p);
  }
  
#ifdef __show_banff_tuning_wght_mesg_debug__
  if (wght != 1.)
  {
     cout << "Exp/NuXSec/FinTuningWght >> "
          << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m)
          << ", Ereco = " << Ereco << " GeV" 
          << ", Etrue = " << Etrue << " GeV => Weight [" 
          << util::AsString(p) << "] = " << wght 
          << endl;   
  }
#endif

  return wght;
}
//________________________________________________________________________________________
double NuXSec::FinTuningWght(
  SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue, const TH2D * mctempl) const
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
     double wght = this->FinTuningWght(p,s,m,Ereco,Etrue);
     // If the weight was calculated for a large var-size Ereco,Etrue bin and needs to
     // shared between multiple MC template bins, calculate the weight fraction that
     // belongs to the current MC template bin.
     bool share = kXSecSystWghtUnequallySharedBetweenBins[isyst];
     double frac = (share) ? 
        this->FractionOfEvents (p,s,m,Ereco,Etrue,mctempl) : 1.;
     // Add to total weight
     tot_wght *= (frac*wght);
  }

  return tot_wght;
}
//________________________________________________________________________________________
double NuXSec::SystWght(
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

  if( ! this->ApplySyst(p,s,m,Ereco,Etrue) ) return 1.;

  double wght = 1.;

  //
  // BANFF param # 22, M_{A} QE
  // BANFF param # 23, M_{A} RES
  //
  if(p == kNuisBANFF22 || p == kNuisBANFF23)
  {
     TSpline3 * wght_spl = this->FindWghtSpline(p,s,m,Ereco,Etrue);
     assert(wght_spl);
     double tb = fBANFF->BestFitValue(p);
     double tc = fBANFF->CurrentValue(p);
     wght = wght_spl->Eval(tc) / wght_spl->Eval(tb);
     wght = TMath::Max(0.,wght); 
  }

  //
  // BANFF param # 24, CCQE  normalization, 0 - 1.5 GeV
  // BANFF param # 25, CCRES normalization, 0 - 2.5 GeV
  //  
  else
  if(p == kNuisBANFF24 || p == kNuisBANFF25)
  {
     wght = fBANFF->CurrentValue(p) / fBANFF->BestFitValue(p);
  }

#ifdef __show_syst_var_wght_mesg_debug__
  if (wght != 1.)
  {
     cout << "Exp/NuXSec/SystWght >> "
          << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m)
          << ", Ereco = " << Ereco << " GeV" 
          << ", Etrue = " << Etrue << " GeV => Weight [" 
          << util::AsString(p) << "] = " << wght 
          << endl;   
  }
#endif
  
  return wght;
}
//________________________________________________________________________________________
double NuXSec::SystWght(
  SuperkSample_t s, SuperkMode_t m, double Ereco, double Etrue, const TH2D * mctempl) const
{
  double tot_wght = 1.;

  for(int isyst = 0; isyst < kNumOfXSecSyst; isyst++) {
     PhysParamId_t p = kListOfXSecSyst[isyst];
     // Is the current systematic parameter applied to current sample/mode end energy bin?
     bool apply = this->ApplySyst(p,s,m,Ereco,Etrue);
     if(!apply) continue;
     // Calculate weight
     double wght = this->SystWght(p,s,m,Ereco,Etrue);
     // If the weight was calculated for a large var-size Ereco,Etrue bin and needs to
     // shared between multiple MC template bins, calculate the weight fraction that
     // belongs to the current MC template bin.
     bool share = kXSecSystWghtUnequallySharedBetweenBins[isyst];
     double frac = (share) ? 
         this->FractionOfEvents (p,s,m,Ereco,Etrue,mctempl) : 1.;
     // Add to total weight
     tot_wght *= (frac*wght);
  }

  return tot_wght;
}
//________________________________________________________________________________________

double NuXSec::CurrentValue(PhysParamId_t p) const
{
  return 0;
}
//________________________________________________________________________________________
void NuXSec::PrintSystematicVariations(void) const
{
  TVectorD & fsi_current_values = *fFSICurrentValues;
  cout << "Exp/NuXSec >> FSI parameters: " << endl;
  util::PrintCompact(fsi_current_values,12,6,"FSI params");
}
//________________________________________________________________________________________

//
// FSI 
//

//________________________________________________________________________________________
void NuXSec::TweakFSI(
  bool       force_inputs, /* force the input params or generate random ones internally  */
  TVectorD * input_vec     /* input vector of params - Set to 0 to reset to nominal ones */
)
{
#ifdef __show_fsi_cholesky_mesg_debug__
  cout << "Exp/NuXSec/FSI/Tweak >> Tweaking BANFF fit parameters" << endl;
#endif

  //
  // If this is the first time this function is called, then run the Cholesky
  // decomposition of the flux covariance matrix
  //

  if(fFSICovMtrxChoL == 0)
  {
#ifdef __show_fsi_cholesky_mesg_info__
    cout << "Exp/NuXSec/FSI/Tweak >> "
         << "Running Cholesky's decomposition method on the FSI covariance matrix." << endl;
#endif
    
    TMatrixD cov(this->FSICovMatrix());
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/NuXSec/FSI/Tweak >> FSI covariance matrix: " << endl;
    util::PrintCompact(cov,7,4);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fFSICovMtrxChoL = new TMatrixD(L);
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/NuXSec/FSI/Tweak >> Cholesky decomposition of FSI covariance matrix: " << endl;
    util::PrintCompact(*fFSICovMtrxChoL,7,4);
#endif
  }

  TVectorD & current_values = *fFSICurrentValues;

  //
  // Force input FSI parameters
  //
  if(force_inputs)
  {
     // Input vector is NULL.
     // This resets the FSI params to their best-fit values.
     if(input_vec==0) {
#ifdef __show_fsi_tweak_mesg_debug__
        cout << "Exp/NuXSec/FSI/Tweak >> Reset FSI params to nominal values." << endl;
#endif
        const int np = conv::kNumOfFSIParams;
        for(int ip = 0; ip < np; ip++) {
           current_values[ip] = 1.;
        }
        return;
     }
     else
     {
#ifdef __show_fsi_tweak_mesg_debug__
        cout << "Exp/NuXSec/FSI/Tweak >> Forcing input BANFF parameters. " << endl;
#endif
        TVectorD & vec = *input_vec;
        current_values = vec;
        return;
     }
  }//force_inputs?

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  else
  {   
#ifdef __show_fsi_tweak_mesg_debug__
     cout << "Exp/NuXSec/FSI/Tweak >> Generating random FSI params (correlated)! " << endl;
#endif
     TMatrixD & L = *fFSICovMtrxChoL; 
     TVectorD vec = util::CholeskyGenerateCorrelatedParamVariations(L);
     current_values = vec;
     return;
  }

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
  if(idx == -1) 
  {
    cerr << "Exp/NuXSec/FindWghtSpline >> Didn't find weight spline for "
         << "sample = "      << util::AsString(s) 
         << ", mode = "      << util::AsString(m)
         << ", parameter = " << util::AsString(p)
         << ", Ereco = "     << Ereco << " GeV"
         << ", Etrue = "     << Etrue << " GeV!"
         << endl;
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

  int is = util::FindSampleSlot(s);
  if(is == -1) return -1;

  int im = util::FindModeSlot(s,m);
  if(im == -1) return -1;

  int ip = util::FindSlotInListOfXSecSyst(p);
  if(ip == -1) return -1;

  int ireco = -1; 
  for (int i = 0; i < kXsNErecoBins; i++) {
    if(kXsErecoBinEdges[i]   <= Ereco && 
       kXsErecoBinEdges[i+1] >  Ereco) 
    {
        ireco = i;
        break;
    }
  }
  if(ireco == -1) return -1;

  int itrue = -1; 
  for (int i = 0; i < kXsNEtrueBins; i++) {
    if(kXsEtrueBinEdges[is][im][i]   <= Etrue && 
       kXsEtrueBinEdges[is][im][i+1] >  Etrue) 
    {
        itrue = i;
        break;
    }
  }
  if(itrue == -1) return -1;

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
  int isyst = util::FindSlotInListOfXSecSyst(p);
  if (isyst == -1) return false;

  int is = util::FindSampleSlot(s);
  if (is == -1) return false;
  int im = util::FindModeSlot(s,m);
  if (im == -1) return false;

  bool apply = kXsApplyXSecSystToSKMode[is][im][isyst];
  return apply;
}
//________________________________________________________________________________________
bool NuXSec::ApplySystToEnergyBin(PhysParamId_t p, double Ereco, double Etrue) const
{
  int isyst = util::FindSlotInListOfXSecSyst(p);
  if (isyst == -1) return false;

  bool apply = 
     (Ereco >= kXSecSystErecoCut[isyst][0] && Ereco < kXSecSystErecoCut[isyst][1]) &&
     (Etrue >= kXSecSystEtrueCut[isyst][0] && Etrue < kXSecSystEtrueCut[isyst][1]);

  return apply;
}
//________________________________________________________________________________________
double NuXSec::FractionOfEvents(
  PhysParamId_t p, SuperkSample_t s, SuperkMode_t m, 
  double Ereco, double Etrue, const TH2D * mctempl) const
{
  assert(mctempl);

  int is = util::FindSampleSlot(s);
  int im = util::FindModeSlot(s,m);
  assert(is != -1 && im != -1);

  int ip = util::FindSlotInListOfXSecSyst(p);
  if(ip == -1) {
     cerr << "Exp/NuXSec/FractionOfEvents >> Can't find index for cross-section "
          << "systematic param: " << util::AsString(p) << endl;
     exit(1);
  }

  // Look-up the Ereco and Etrue energy range (in the variable-size binning
  // scheme used when pre-computing weight splines)
  double Ereco_min = -1;
  double Ereco_max = -1;
  double Etrue_min = -1;
  double Etrue_max = -1;
  for(int ir = 0; ir < kXsNErecoBins; ir++) {
    if(Ereco >= kXsErecoBinEdges[ir] && Ereco < kXsErecoBinEdges[ir+1]) {
        Ereco_min = kXsErecoBinEdges[ir];
        Ereco_max = kXsErecoBinEdges[ir+1];
        break;
    }
  }
  for(int it = 0; it < kXsNEtrueBins; it++) {
    if(Etrue >= kXsEtrueBinEdges[is][im][it] && Etrue < kXsEtrueBinEdges[is][im][it+1]) {
        Etrue_min = kXsEtrueBinEdges[is][im][it];
        Etrue_max = kXsEtrueBinEdges[is][im][it+1];
        break;
    }
  }

  if(Ereco < Ereco_min || Ereco > Ereco_max) return 0.;
  if(Etrue < Etrue_min || Etrue > Etrue_max) return 0.;

  int ir0 = mctempl->GetXaxis()->FindBin(Ereco);
  int it0 = mctempl->GetYaxis()->FindBin(Etrue);

  double n0 = mctempl->GetBinContent(ir0,it0);

  int ir_min = mctempl->GetXaxis()->FindBin(Ereco_min + kOneMeV);
  int ir_max = mctempl->GetXaxis()->FindBin(Ereco_max - kOneMeV);
  int it_min = mctempl->GetYaxis()->FindBin(Etrue_min + kOneMeV);
  int it_max = mctempl->GetYaxis()->FindBin(Etrue_max - kOneMeV);

  double ntot_range = 0;
  for(int ir = ir_min; ir <= ir_max; ir++) {
    for(int it = it_min; it <= it_max; it++) {
       ntot_range += mctempl->GetBinContent(ir,it);
    }
  }

  assert(ntot_range >= n0);   
      
  double frac = 0;
  if(ntot_range > 0) {
    frac = n0 / ntot_range;
  }
          
  return frac;
}
//________________________________________________________________________________________
