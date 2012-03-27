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
#define __show_fsi_cholesky_mesg_info__
#define __show_fsi_cholesky_mesg_debug__
#define __show_fsi_tweak_mesg_debug__
//#define __show_build_frac_template_mesg_info__

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include <TMath.h>

#include "FSI.h"
#include "Utilities.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;

using namespace t2kana;

ClassImp(FSI)

//________________________________________________________________________________________
FSI::FSI()
{
  this->Init();
}
//________________________________________________________________________________________
FSI::~FSI()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void FSI::Init(void)
{
  fCurrentValues = 0;
  fNominalValues = 0;
  fCovMtrx       = 0;
  fCovMtrxChoL   = 0;

  for(int is = 0; is < kns; is++) {
     fFracTun[is] = 0;
  }

}
//________________________________________________________________________________________
void FSI::CleanUp(void)
{
  if ( fCurrentValues ) delete fCurrentValues;
  if ( fNominalValues ) delete fNominalValues;
  if ( fCovMtrx       ) delete fCovMtrx;
  if ( fCovMtrxChoL   ) delete fCovMtrxChoL;

  for(int is = 0; is < kns; is++) {
     if (fFracTun[is]) delete fFracTun[is];
  }
}
//________________________________________________________________________________________
bool FSI::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/FSI >> Reading FSI inputs from input XML file." << endl;
#endif

  const int np = conv::kNumOfFSIParams;

/*
  // Get FSI covariance matrix
  fCovMtrx = util::XML2TMatrixD(xml,root_node,"fsi/covariance_matrix");
  if(!fCovMtrx) {
    cerr << "Exp/FSI >> No FSI covariance matrix was found." << endl;
    exit(1);
  }

  const TMatrixD & cov = *fCovMtrx;
  if(cov.GetNrows() != np || cov.GetNcols() != np)
  {
    cerr << "Exp/FSI >> BANFF covariance matrix has unexpected size!" << endl;
    cerr << "Exp/FSI >> Expecting: " <<  np << " x " << np << " matrix" << endl;
    cerr << "Exp/FSI >> Found: " << cov.GetNrows() << " x " <<  cov.GetNcols() << " matrix" << endl;
  }
*/

// Temporary code: 
// As I don't have an FSI cov matrix, create a fake one so as to have something to
// test the code with. 
// ------------------------------------------------------------
fCovMtrx = new TMatrixD(np,np);
TMatrixD & cov = *fCovMtrx;
for(int i=0; i<np; i++) {
  cov(i,i) = 0.5/(1.+i);
}
for(int i=0; i<np; i++) {
  for(int j=0; j<np; j++) {
    if(i >= j) continue;
    cov(i,j) = 0;
    cov(j,i) = 0;
  }
}
// ------------------------------------------------------------

#ifdef __show_init_mesg_info__
  cout << "Exp/FSI >> Input FSI covariance matrix: " << endl;
  util::PrintCompact(cov,7,4);
#endif

  // Create/init current FSI params
  fCurrentValues = new TVectorD(np);
  fNominalValues = new TVectorD(np);
  for(int ip = 0; ip < np; ip++) 
  {
     (*fCurrentValues)[ip] = 1.;
     (*fNominalValues)[ip] = 1.;
  }

  //
  //
  //

  for(int is = 0; is < kns; is++) {
     conv::SuperkSample_t s = conv::kSKSamples[is];
     if(s == conv::kSKsUndefined) continue;
     const int      n_ereco_bins    = util::NumErecoBins (s);
     const double * ereco_bin_edges = util::ErecoBinEdges(s);
     fFracTun[is] = new TH1D(
        Form("nsk_tun_%d",is), util::AsStringCompact(s).c_str(), n_ereco_bins, ereco_bin_edges);
     fFracTun[is]->SetDirectory(0);
     fFracTun[is]->GetXaxis()->SetTitle("E_{reco} (GeV)");
     fFracTun[is]->GetYaxis()->SetTitle("Fraction");
  }//is

#ifdef __show_init_mesg_info__
  cout << "Exp/FSI >> Done reading neutrino cross-section inputs!" << endl;
#endif

  return true;  
}
//________________________________________________________________________________________
double FSI::SystWght(
  conv::SuperkSample_t s, conv::SuperkMode_t /*m*/, double Ereco, double /*Etrue*/) const
{
  // find FSI param which corresponds to the current Ereco param
  const int np = conv::kNumOfFSIParams;
  int ip = TMath::BinarySearch(np+1, conv::kFSIParamsErecoBinEdges, Ereco);
  if(ip < 0 || ip >= np) {
     cerr << "Exp/FSI >> Weird FSI Ereco index for Ereco value = " << Ereco << endl;
     exit(1);
  }
 
  // get current value of ip^th FSI param
  double fsi_param_value = (*fCurrentValues)[ip];

  // param tweaked
  bool tweaked = TMath::Abs(fsi_param_value - 1.) > conv::kASmallNum;
  if(!tweaked) {
    return 1.;
  }

  // get fraction of events in current MC template Ereco bin relative to the 
  // number of events in the current (wider) FSI Ereco bin
  int is =  util::ArrIdx::Instance()->SuperkSampleIdx(s); 
  if(is < 0) {
     cerr << "Exp/FSI >> Negative idx for sample: " << util::AsString(s) << endl;
     exit(1);
  }
  int ir =  fFracTun[is]->FindBin(Ereco); 
  double frac = fFracTun[is]->GetBinContent(ir);

  // Calculate FSI weight for current MC template bin
  double fsi_wght = frac * fsi_param_value;
  return fsi_wght;
}
//________________________________________________________________________________________
double FSI::CurrentValue(conv::PhysParamId_t p) const
{
  const TVectorD & current_values = *fCurrentValues;
  int idx = this->Param2Idx(p);
  return current_values[idx];
}
//________________________________________________________________________________________
void FSI::Tweak(
  bool       force_inputs,    /* force the input params or generate random ones internally  */
  TVectorD * input_fsi_params /* input vector of FSI params - Set to NULL to reset to nominal ones */
)
{
#ifdef __show_fsi_cholesky_mesg_debug__
  cout << "Exp/FSI/Tweak >> Tweaking FSI parameters" << endl;
#endif

  //
  // If this is the first time this function is called, then run the Cholesky
  // decomposition of the flux covariance matrix
  //

  if(fCovMtrxChoL == 0)
  {
#ifdef __show_fsi_cholesky_mesg_info__
    cout << "Exp/FSI/Tweak >> "
         << "Running Cholesky's decomposition method on the FSI covariance matrix." << endl;
#endif
    
    TMatrixD cov(this->CovMatrix());
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/FSI/Tweak >> FSI covariance matrix: " << endl;
    util::PrintCompact(cov,7,4);
#endif

    TMatrixD L = util::CholeskyDecomposition(cov);
    fCovMtrxChoL = new TMatrixD(L);
#ifdef __show_fsi_cholesky_mesg_debug__
    cout << "Exp/FSI/Tweak >> Cholesky decomposition of FSI covariance matrix: " << endl;
    util::PrintCompact(*fCovMtrxChoL,7,4);
#endif
  }

  TVectorD & current_values = *fCurrentValues;

  //
  // Force input FSI parameters
  //
  if(force_inputs)
  {
     // Input vector is NULL.
     // This resets the FSI params to their best-fit values.
     if(input_fsi_params==0) {
#ifdef __show_fsi_tweak_mesg_debug__
        cout << "Exp/FSI/Tweak >> Reset FSI params to nominal values." << endl;
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
        cout << "Exp/FSI/Tweak >> Forcing input FSI parameters. " << endl;
#endif
        TVectorD & vec = *input_fsi_params;
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
     cout << "Exp/FSI/Tweak >> Generating random FSI params (correlated)! " << endl;
#endif
     TMatrixD & L = *fCovMtrxChoL; 
     TVectorD & nominal_values = *fNominalValues;
     TVectorD vec = util::CholeskyGenerateCorrelatedParams(L,nominal_values);
     current_values = vec;
     return;
  }
}
//________________________________________________________________________________________
double FSI::ChisqPenalty(void) const
{
  return 0;
}
//________________________________________________________________________________________
void FSI::PrintSystematicVariations(void) const
{
  const TVectorD & curr_fsi_params = *fCurrentValues;
  cout << "Exp/FSI >> FSI parameters: " << endl;
  util::PrintCompact(curr_fsi_params,15,6,"FSI params");
}
//________________________________________________________________________________________
int FSI::Param2Idx(conv::PhysParamId_t p) const
{
  int iparam = (int) p;
  int first  = (int) conv::kNuisFSI00;
  int last   = (int) conv::kNuisFSI09;
  assert(iparam >= first && iparam <= last);
  int idx = iparam-first;
  return idx;
}
//________________________________________________________________________________________
void FSI::BuildFracMCTemplates(TH2D* input_mctempl_array[kns][knm])
{
#ifdef __show_build_frac_template_mesg_info__
       cout << "Exp/FSI >> Building MC template fractions... " << endl;
#endif  

  const int np = conv::kNumOfFSIParams;

  for(int is = 0; is < kns; is++) {
     conv::SuperkSample_t s = conv::kSKSamples[is];
     if(s == conv::kSKsUndefined) continue;
     TH1D * mctempl_frac = fFracTun[is];
     assert(mctempl_frac);   
     // Loop over Ereco bins (MC template binning scheme)
     int nr = mctempl_frac->GetXaxis()->GetNbins();
     for(int ir = 1; ir <= nr; ir++) {
        // Loop over all modes and all Etrue bins and calculate the number of events
        double nev = 0;
        for(int im = 0; im < knm; im++) {
           conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
           if(m == conv::kSKmUndefined) continue;
           const TH2D * mctempl = input_mctempl_array[is][im];
           assert(mctempl);
           int nt = mctempl->GetYaxis()->GetNbins();
           for(int it = 1; it <= nt; it++) {
              nev += mctempl->GetBinContent(ir,it);
           }//it
        }//im
        mctempl_frac->SetBinContent(ir,nev);
     }//ir

     // Loop over FSI params and look-up the Ereco energy range for each 
     for(int irv = 0; irv < np; irv++) {
        double Ereco_min = conv::kFSIParamsErecoBinEdges[irv];
        double Ereco_max = conv::kFSIParamsErecoBinEdges[irv+1];
        if(Ereco_min < 0 || Ereco_max < 0) continue;
        Ereco_min += conv::kOneMeV;
        Ereco_max -= conv::kOneMeV;
        int ir_min = mctempl_frac->GetXaxis()->FindBin(Ereco_min);
        int ir_max = mctempl_frac->GetXaxis()->FindBin(Ereco_max);
        double sum = 0;
        for(int ir = ir_min; ir <= ir_max; ir++) {
          double n = mctempl_frac->GetBinContent(ir);
          sum += n;
        }//ir
        for(int ir = ir_min; ir <= ir_max; ir++) {
          double n = mctempl_frac->GetBinContent(ir);
          double frac = (sum>0) ? n/sum : 0;
          mctempl_frac->SetBinContent(ir,frac);
        }//ir
     }//irv
  }//is
     
}
//________________________________________________________________________________________
  

