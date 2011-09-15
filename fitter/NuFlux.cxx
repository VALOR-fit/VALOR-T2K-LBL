//________________________________________________________________________________________
/*
 Changes since v1.0.0:

 @ 16/03/2011 - CA
   Created from the old SuperkExpectation::BeamInputs nested class. Rm numu10a namespace.
 @ 10/06/2011 - CA
   Modified code so as to use directly the numu flux fractional covariance matrix provided 
   by the Beam WG. The flux covariance matrix needed for applying the Cholesky method is
   now constructed on the fly during initialisation. 
 @ 05/07/2011 - CA
   Modified TweakShape() method to allow fitting for flux shape systematics. Using the
   same strategy as with the SK CCQE ring-counting err. The TweakShape() method can 
   generate random variations (as we were doing before in order to include errors in the
   contour using the Cousins-Highland methd) but one can also force externally specified 
   variations (as we would need to do if we are fitting for flux shape systematics).
   The forced tweak dials can be the correlated or de-correlated ones. Forcing the values
   of de-correlated tweak dials is employed in the fit. Forcing the values of the 
   correlated tweak dials is employed in certain validation tests (hence this capability
   was retained). Added method to IgnoreShapeTweaking() so as to inhibit effect and
   calculate the appropriate normalisation factor in SuperkExpectation (before modifying
   the TweakShape() method this was accomplished by running TweakShape(false) but now
   the TweakShape() arguments are used for different purpose). In order to simplify the
   TweakShape() method, I added the FillInTweakedNuMuFluxHisto().
*/
//________________________________________________________________________________________

//#define __debug_level_1__
//#define __debug_level_2__
//#define __debug_level_3__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TVectorD.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TSystem.h>

#include "NuFlux.h" 
#include "Conventions.h" 
#include "Utilities.h" 
#include "HadErrCalc.h" ///<-- code by Kubo-san

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;
using namespace t2kana::util;

ClassImp(NuFlux)

//________________________________________________________________________________________
NuFlux::NuFlux()
{
  this->Init();
}
//________________________________________________________________________________________
NuFlux::~NuFlux()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void NuFlux::Init(void)
{
  fSuperkNuMuFlux10a           = 0; 
  fSuperkNuMuFlux10aVarBin     = 0; 
  fSuperkNuMuFluxPOT           = 0;
  fSuperkNuMuFluxCov           = 0;
  fSuperkNuMuWght10aTo10d      = 0;
  fSuperkNuMuBarWght10aTo10d   = 0;
  fSuperkNuEWght10aTo10d       = 0;
  fSuperkNuMuFluxNormTwkDial   = 0;
  fSuperkNuMuFluxCovCholeskyL  = 0;
  fSuperkNuMuFluxHaveCholeskyL = false;
  fSuperkNuMuFluxTweaked       = 0;
  fFluxShapeWasTweaked         = false;
  fIgnoreFluxShapeTweaks       = false;
  fHEC                         = 0;
  fSuperkNormRelErr_simple     = 0;
}
//________________________________________________________________________________________
void NuFlux::CleanUp(void)
{
  if ( fSuperkNuMuFlux10a          ) delete fSuperkNuMuFlux10a;
  if ( fSuperkNuMuFlux10aVarBin    ) delete fSuperkNuMuFlux10aVarBin;
  if ( fSuperkNuMuFluxCov          ) delete fSuperkNuMuFluxCov;
  if ( fSuperkNuMuWght10aTo10d     ) delete fSuperkNuMuWght10aTo10d;
  if ( fSuperkNuMuBarWght10aTo10d  ) delete fSuperkNuMuBarWght10aTo10d;
  if ( fSuperkNuEWght10aTo10d      ) delete fSuperkNuEWght10aTo10d;
  if ( fSuperkNuMuFluxTweaked      ) delete fSuperkNuMuFluxTweaked;
  if ( fSuperkNuMuFluxCovCholeskyL ) delete fSuperkNuMuFluxCovCholeskyL;
  if ( fHEC                        ) delete fHEC;
}
//________________________________________________________________________________________
bool NuFlux::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read beam inputs to disappearance analysis
//
  cout << "Exp/Beam >> Reading beam inputs from input XML file." << endl;

  //
  // 10a flux
  //

  // Get SuperK numu flux
  fSuperkNuMuFlux10a = util::XML2TH1D(xml,root_node,"beam/superk_flux/v10a/numu");
  if(!fSuperkNuMuFlux10a) {
    cerr << "Exp/Beam >> No input SuperK numu flux histogram!" << endl;
    exit(1);
  }
  // Get SuperK numu flux (var binning)
  fSuperkNuMuFlux10aVarBin = util::XML2TH1D(xml,root_node,"beam/superk_flux/v10a/numu_varbin");
  if(!fSuperkNuMuFlux10aVarBin) {
    cerr << "Exp/Beam >> No input SuperK numu (var binning) flux histogram!" << endl;
    exit(1);
  }
  // Integrated POT for numu flux simulation
  fSuperkNuMuFluxPOT = util::XML2Dbl(xml,root_node,"beam/superk_flux/pot_mc");
  cout << "Exp/Beam >> SK flux simulation POT = " << fSuperkNuMuFluxPOT << endl;

  //
  // 10a->10d weights
  //

  fSuperkNuMuWght10aTo10d = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/weight_v10a_to_v10d/numu");
  if(!fSuperkNuMuWght10aTo10d) {
    cerr << "Exp/Beam >> No input 10a->10d SuperK numu flux weights!" << endl;
    exit(1);
  }
  fSuperkNuMuBarWght10aTo10d = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/weight_v10a_to_v10d/numubar");
  if(!fSuperkNuMuBarWght10aTo10d) {
    cerr << "Exp/Beam >> No input 10a->10d SuperK numubar flux weights!" << endl;
    exit(1);
  }
  fSuperkNuEWght10aTo10d = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/weight_v10a_to_v10d/nue");
  if(!fSuperkNuEWght10aTo10d) {
    cerr << "Exp/Beam >> No input 10a->10d SuperK nue flux weights!" << endl;
    exit(1);
  }

  //
  // flux normalisation err
  //

  cout << "Exp/Beam >> Loading Kubo-san's flux norm err (Nsk/Nnd) inputs!" << endl;

  fSuperkNormRelErr_simple = 0;

  string kubo_data_path =
      util::XML2String(xml,root_node,"beam/superk_flux/systematics/norm_error_inputs");
  if(kubo_data_path != "none") {
     string full_kubo_data_path = gSystem->Getenv("T2KOSCANA") + string("/") + kubo_data_path;
     cout << "Exp/Beam >> Using input files from: " << full_kubo_data_path << endl;
     fHEC = new HadErrCalc(full_kubo_data_path.c_str());
  } else {
     cout << "Exp/Beam >> Not using Kubo's code / Use a fixed Nsk/Nd280 err instead." << endl;
     fSuperkNormRelErr_simple = util::XML2Dbl(
         xml,root_node,"beam/superk_flux/systematics/norm_error_simple_version");
     cout << "Exp/Beam >> Norm relative err (simple) : " << fSuperkNormRelErr_simple << endl;
  }

  //
  // flux shape err
  //

  cout << "Exp/Beam >> Loading flux shape err inputs!" << endl;

  TH2D * sk_numu_flux_fractional_covar = util::XML2TH2D(
       xml,root_node,"beam/superk_flux/systematics/fractional_covariance/numu");
  this->BuildNuMuFluxCovMatrix(sk_numu_flux_fractional_covar);
  delete sk_numu_flux_fractional_covar;

  cout << "Exp/Beam >> Done reading beam inputs!" << endl;

  return true;
}
//________________________________________________________________________________________
void NuFlux::BuildNuMuFluxCovMatrix(TH2D * fractional_covariance)
{
  // Build 10d numu flux (50 MeV binning)
  int    nbins = fSuperkNuMuFlux10a->GetXaxis()->GetNbins();
  double emin  = fSuperkNuMuFlux10a->GetXaxis()->GetXmin();
  double emax  = fSuperkNuMuFlux10a->GetXaxis()->GetXmax();
  TH1D * flux10d = new TH1D("flux10d","",nbins,emin,emax);
  bool wght_to_10d = true;
  for(int ibin = 1; ibin <= nbins; ibin++) {
     double Etrue = flux10d->GetBinCenter(ibin);
     double flux  = this->SuperkNuMuFlux(Etrue,wght_to_10d);
     flux10d->SetBinContent(ibin,flux);
  }
  
  // Re-bin using the Beam WG variable binning scheme
  TH1D * flux10d_varbin = (TH1D*) flux10d->Rebin(
        conv::kNEtrueBinsBeamWG,"10dvarbin",conv::kEtrueBinEdgesBeamWG);
       
 // Now build the flux covariance matrix from this histogram 
 fSuperkNuMuFluxCov = new TMatrixD(conv::kNEtrueBinsBeamWG,conv::kNEtrueBinsBeamWG);   
 TMatrixD & cov = *fSuperkNuMuFluxCov;
 for(int ibin = 1; ibin <= conv::kNEtrueBinsBeamWG; ibin++) {
   double flux_i = flux10d_varbin->GetBinContent(ibin);
   for(int jbin = 1; jbin <= conv::kNEtrueBinsBeamWG; jbin++) {
      double flux_j = flux10d_varbin->GetBinContent(jbin);
      double frac_cov = fractional_covariance->GetBinContent(ibin,jbin);
      cov[ibin-1][jbin-1] = (flux_i*flux_j) * frac_cov;
#ifdef __debug_level_3__
      cout << "Exp/Beam/BuildCovarianceMatrix >> " 
           << "CovFrac(" << ibin << ", " << jbin << ") = " << frac_cov
           << ", F(" << ibin << ") = " << flux_i 
           << ", F(" << jbin << ") = " << flux_j
           << " => Cov(" << ibin << ", " << jbin << ") = " << cov[ibin-1][jbin-1]
           << endl;
#endif
    }
  }

#ifdef __debug_level_3__
  fSuperkNuMuFluxCov->Print();
#endif

  delete flux10d;
  delete flux10d_varbin;
}
//________________________________________________________________________________________
void NuFlux::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  if(opt.Contains("SKFlux/Weights/10aTo10d/NuMu")) { 
    fSuperkNuMuWght10aTo10d->Draw();
    fSuperkNuMuWght10aTo10d->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSuperkNuMuWght10aTo10d->GetYaxis()->SetTitle("SuperK numu flux weight");
    c->Update();
  }
  else
  if(opt.Contains("SKFlux/10a/NuMu")) { 
    fSuperkNuMuFlux10a->Draw();
    fSuperkNuMuFlux10a->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSuperkNuMuFlux10a->GetYaxis()->SetTitle("SuperK numu flux (#neutrinos /cm2 /POT)");
    c->Update();
  } 
  else
  if(opt.Contains("SKFlux/10d/NuMu")) { 
    TH1D * weighted_flux = new TH1D(*fSuperkNuMuFlux10a);
    weighted_flux->SetDirectory(0);
    for(int bin=1; bin <= weighted_flux->GetNbinsX(); bin++) {
	double Et = weighted_flux->GetBinCenter(bin);
	double f  = this->SuperkNuMuFlux(Et,true);
        weighted_flux->SetBinContent(bin,f);
    }
    weighted_flux->Draw();
    weighted_flux->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    weighted_flux->GetYaxis()->SetTitle("SuperK numu flux (#neutrinos /cm2 /POT)");
    c->Update();
  }
  else
  if(opt.Contains("SKFlux/Both10a10d/NuMu")) { 
    TH1D * weighted_flux = new TH1D(*fSuperkNuMuFlux10a);
    weighted_flux->SetDirectory(0);
    for(int bin=1; bin <= weighted_flux->GetNbinsX(); bin++) {
	double Et = weighted_flux->GetBinCenter(bin);
	double f  = this->SuperkNuMuFlux(Et,true);
        weighted_flux->SetBinContent(bin,f);
    }
    weighted_flux->SetLineStyle(kDashed);
    fSuperkNuMuFlux10a->Draw();
    weighted_flux->Draw("same");
    fSuperkNuMuFlux10a->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fSuperkNuMuFlux10a->GetYaxis()->SetTitle("SuperK numu flux (#neutrinos /cm2 /POT)");
    c->Update();
  }
  else
  if(opt.Contains("SKFlux/CovarianceMatrix/NuMu")) { 
    fSuperkNuMuFluxCov->Draw("colz");
    c->Update();
  }
}
//________________________________________________________________________________________
double NuFlux::SuperkNuMuFlux10aTo10dWght(double Etrue) const
{
  double wght = 1.0;
  if(Etrue < 10/*GeV*/) {
    int wght_bin = fSuperkNuMuWght10aTo10d->FindBin(Etrue);
    wght = fSuperkNuMuWght10aTo10d->GetBinContent(wght_bin);
  }
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkNuMuBarFlux10aTo10dWght(double Etrue) const
{
  double wght = 1.0;
  if(Etrue < 10/*GeV*/) {
    int wght_bin = fSuperkNuMuBarWght10aTo10d->FindBin(Etrue);
    wght = fSuperkNuMuBarWght10aTo10d->GetBinContent(wght_bin);
  }
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkNuEFlux10aTo10dWght(double Etrue) const
{
  double wght = 1.0;
  if(Etrue < 10/*GeV*/) {
    int wght_bin = fSuperkNuEWght10aTo10d->FindBin(Etrue);
    wght = fSuperkNuEWght10aTo10d->GetBinContent(wght_bin);
  }
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkNuMuFlux(double Etrue, bool wght_to_10d) const
{
  assert(fSuperkNuMuFlux10a);

  int flux_bin = fSuperkNuMuFlux10a->FindBin(Etrue);
  double flux  = fSuperkNuMuFlux10a->GetBinContent(flux_bin);
 
  if(wght_to_10d) {
    double wght = this->SuperkNuMuFlux10aTo10dWght(Etrue);
    flux *= wght;
  }

  return flux;
}
//________________________________________________________________________________________
double NuFlux::SuperkNuMuFluxPOT (void) const
{
  return fSuperkNuMuFluxPOT;
}
//________________________________________________________________________________________
const TMatrixD &  NuFlux::SuperkNuMuFluxCov(void) const 
{ 
  return *fSuperkNuMuFluxCov; 
} 
//________________________________________________________________________________________
double NuFlux::SuperkFluxNormErr(
  TH1D * hSkNuMu_1RingMu, TH1D * hSkNuMuBar_1RingMu, TH1D * hSkNuE_1RingMu,
  TH1D * hNdNuMu,         TH1D * hNdNuMuBar,         TH1D * hNdNuE) const
{
// normalization err Nsk/Nd280

  // don't  use Kubo-san's code. 
  if(!fHEC) {
    return fSuperkNormRelErr_simple;
  }

#ifdef __debug_level_3__
  TFile f("./inputs_to_beam_wg_norm_err_calc_code.root","recreate");
  TH1D * h1 = new TH1D(*hSkNuMu_1RingMu);
  TH1D * h2 = new TH1D(*hSkNuMuBar_1RingMu);
  TH1D * h3 = new TH1D(*hSkNuE_1RingMu);
  TH1D * h4 = new TH1D(*hNdNuMu);
  TH1D * h5 = new TH1D(*hNdNuMuBar);
  TH1D * h6 = new TH1D(*hNdNuE);
  h1 -> Write("sk1rmu_numu");
  h2 -> Write("sk1rmu_numubar");
  h3 -> Write("sk1rmu_nue");
  h4 -> Write("nd_numu");
  h5 -> Write("nd_numubar");
  h6 -> Write("nd_nue");
  f.Close();
#endif

  // use code contributed by Kubo-san's on behalf of the Beam WG
  fHEC->Set1RmuHist(hSkNuMu_1RingMu, hSkNuMuBar_1RingMu, hSkNuE_1RingMu);
  fHEC->SetNDHist  (hNdNuMu,         hNdNuMuBar,         hNdNuE);
  fHEC->Calc();

#ifdef __debug_level_3__
  fHEC->ErrorComponents();
  fHEC->ShowDetail();
#endif

  // get relative err on Nsk/Nnd280
  double relative_err = fHEC->GetSumErr(0);

#ifdef __debug_level_1__
  cout << "Exp/Beam >> Nsk/Nnd280 relative err calculated by Beam WG code = " 
       << relative_err << endl;
#endif

  return relative_err;
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxNorm(
  TH1D * hSkNuMu_1RingMu, TH1D * hSkNuMuBar_1RingMu, TH1D * hSkNuE_1RingMu,
  TH1D * hNdNuMu,         TH1D * hNdNuMuBar,         TH1D * hNdNuE) const
{
// normalization

  double twkdial = fSuperkNuMuFluxNormTwkDial;

  if(TMath::Abs(twkdial) < kASmallNum) return 1; // untweaked - skip calculation

  double relerr  = this->SuperkFluxNormErr(
     hSkNuMu_1RingMu,hSkNuMuBar_1RingMu,hSkNuE_1RingMu,hNdNuMu,hNdNuMuBar,hNdNuE);

#ifdef __debug_level_3__
  cout << "Exp/Beam/SuperkFluxNorm() >> "
       << "Normalization relative error = " << relerr << endl;
#endif
	
  double norm = 1 + twkdial * relerr;
  return norm;
}
//________________________________________________________________________________________
double NuFlux::SuperkNuMuFluxShapeWeight(double Etrue)
{
  // was forced to ignore shape tweaks?
  if(fIgnoreFluxShapeTweaks) {
     return 1.;
  }

  if(!fFluxShapeWasTweaked) {
     return 1;
  }

  if(!fSuperkNuMuFluxHaveCholeskyL) {
     exit(1);
  }
  if(!fSuperkNuMuFluxTweaked || !fSuperkNuMuFlux10a) {
     exit(1);
  }

  int bin = fSuperkNuMuFlux10a->FindBin(Etrue);

  double flux_nom = this->SuperkNuMuFlux(Etrue, true);            //10d nominal
  double flux_twk = fSuperkNuMuFluxTweaked ->GetBinContent(bin);  //10d tweaked

  double wght = 1.;
  if(flux_nom > 0) {
     wght = flux_twk/flux_nom;
  }

#ifdef __debug_level_3__
  cout << "Exp/Beam/SuperkNuMuFluxShapeWeight >> weight(Etrue = " 
       << Etrue << ") = " << wght << endl;
#endif

  return wght;
}
//________________________________________________________________________________________
void NuFlux::TweakNorm (double twkdial) 
{ 
  fSuperkNuMuFluxNormTwkDial = twkdial; 
}
//________________________________________________________________________________________
void NuFlux::TweakShape(
  bool       force_input_twkdials,  /* force the input twk dials or generate twk dials internally */
  TVectorD * input_twkv,            /* input twk dials (if forced) - set 0 to reset */
  bool       correlated,            /* twk dials are the correlated or uncorrelated ones? */
  bool       fractional             /* twk dials expressing a fractional or absolute flux change? */
)    
{ 
#ifdef __debug_level_3__
  cout << "Exp/Beam/TweakShape >> Tweaking SuperK nu_mu flux shape" << endl;
#endif

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the flux covariance matrix
  //
  if(!fSuperkNuMuFluxHaveCholeskyL) {
    cout << "Exp/Beam/TweakShape >> "
         << "Running Cholesky's decomposition method on the flux covariance matrix." << endl;
    TMatrixD  sk_flux_cov(this->SuperkNuMuFluxCov());
    TMatrixD L = util::CholeskyDecomposition(sk_flux_cov);
#ifdef __debug_level_3__
    L.Print();
#endif

    fSuperkNuMuFluxCovCholeskyL = new TMatrixD(L);
#ifdef __debug_level_3__
    fSuperkNuMuFluxCovCholeskyL->Print();
#endif
    fSuperkNuMuFluxHaveCholeskyL = true;
    fSuperkNuMuFluxTweaked = new TH1D(*fSuperkNuMuFlux10a);
  }

  //
  // Forcing null tweaking dials.
  // Reset tweaked flux to nominal flux.
  //
  if(force_input_twkdials && input_twkv==0) {
#ifdef __debug_level_3__
     cout << "Exp/Beam/TweakShape >> RESET! " << endl;
#endif
     this->FillInTweakedNuMuFluxHisto(0);
     return;
  }

  //
  // Check whether the values of the tweaking dials are forced.
  // If the forced tweaking dials are supposed to be the correlated ones then we don't
  // need to use the Cholesky decomposition of the flux covariance matrix.
  //
  if(force_input_twkdials && correlated) {
#ifdef __debug_level_3__
     cout << "Exp/Beam/TweakShape >> Forcing input twk dials (correlated)! " << endl;
#endif
     this->FillInTweakedNuMuFluxHisto(input_twkv,fractional);
     return;
  }

  //
  // Check whether the values of the tweaking dials are forced.
  // If the forced tweaking dials are supposed to be the un-correlated ones then use
  // the Cholesky decomposition of the flux covariance matrix to generate
  // correlated tweaking dials.
  //
  if(force_input_twkdials && !correlated) {
#ifdef __debug_level_3__
     cout << "Exp/Beam/TweakShape >> Forcing input twk dials (uncorrelated)! " << endl;
#endif
     TMatrixD & L = *fSuperkNuMuFluxCovCholeskyL;
     TVectorD twkv = util::CholeskyCalculateCorrelatedParamVariations(L,*input_twkv);
     this->FillInTweakedNuMuFluxHisto(&twkv,false);
     return;
  }

  //
  // If the tweaking dials were not forced, generate correlated tweaking dials
  // using the Cholesky method
  //
  // Use Cholseky method to generate random fluxes
  //
  if(!force_input_twkdials) {
#ifdef __debug_level_3__
    cout << "Exp/Beam/TweakShape >> Generating random twk dials (correlated)! " << endl;
#endif
    TMatrixD & L = *fSuperkNuMuFluxCovCholeskyL;
    TVectorD twkv = util::CholeskyGenerateCorrelatedParamVariations(L);
    this->FillInTweakedNuMuFluxHisto(&twkv,false);
    return;
  }

}
//________________________________________________________________________________________
void NuFlux::FillInTweakedNuMuFluxHisto(TVectorD * ptwkv, bool fractional)
{
// Fill-in the tweaked flux histogram
//

 // reset to nominal
 if(ptwkv==0) {
    fFluxShapeWasTweaked = false;
    for(int j=1; j <= 200; j++){ ///TODO
       double Etrue = fSuperkNuMuFlux10a->GetBinCenter(j);
       double nominal_flux = this->SuperkNuMuFlux(Etrue,true); // include 10a->10d weight
       fSuperkNuMuFluxTweaked->SetBinContent(j, nominal_flux);
    }//end of bin loop
    return;
 } 

 // apply tweak 
 TVectorD & twkv = *ptwkv;
 assert(twkv.GetNrows() == kNEtrueBinsBeamWG);

 fFluxShapeWasTweaked = true;

 //int nbins = fSuperkNuMuFlux10a->GetXaxis()->GetNbins();
 //for(int j=1; j <= nbins; j++){
 for(int j=1; j <= 200; j++){ ///TODO
     // true energy and nominal flux in the current bin
     double Etrue = fSuperkNuMuFlux10a->GetBinCenter(j);
     double nominal_flux = this->SuperkNuMuFlux(Etrue,true); // include 10a->10d weight

     // which bin (in the BeamWG binning scheme) contains the j^th 50-MeV bin?
     int varBin = util::FluxFix2VarBinNumber(j);

     // calculate contribution of the j^th 50-MeV bin in the BeamWG bin which contains it
     double rto = util::FluxFractionsFixedBinRelToVarBin(
                         j, fSuperkNuMuFlux10a, fSuperkNuMuFlux10aVarBin);

     // calculate the tweaked flux
     double tweaked_flux = 0;
     if(fractional) {
        // input tweak vector describes fractional flux change in each bin 
        // BeamWG binning scheme)       
        double fractional_change_beamwg_bin = twkv[varBin - 1];
        tweaked_flux = nominal_flux * (1 + fractional_change_beamwg_bin);
     } else {
        // input tweak vector describes absolute flux change in each bin 
        // (in the BeamWG binning scheme)
        double absolute_change_beamwg_bin = twkv[varBin - 1];
        tweaked_flux = nominal_flux + absolute_change_beamwg_bin * rto;
     }

     // update tweaked flux histogram
     fSuperkNuMuFluxTweaked->SetBinContent(j, tweaked_flux);
  }//end of bin loop
}
//________________________________________________________________________________________

