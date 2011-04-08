//________________________________________________________________________________________
/*
 Changes since v1.0.0:

 @ 16/03/2011 - CA
   Created from the old SuperkExpectation::BeamInputs nested class. Rm numu10a namespace.

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
  fEnableFluxShapeTweaking     = false;
  fHEC                         = 0;
  fSuperkNormRelErr_simple     = 0;
  fTweakShapeVldMode           = false; 
  fTweakShapeVldModeBin        = 0; 
  fTweakShapeVldModeDial       = 0; 
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

  // Get SuperK numu flux covariance matrix
  fSuperkNuMuFluxCov = util::XML2TMatrixD(
      xml,root_node,"beam/superk_flux/systematics/covariance/numu");
  if(!fSuperkNuMuFluxCov) {
    cerr << "Exp/Beam >> No input SuperK numu flux covariance matrix!" << endl;
    exit(1);
  }
#ifdef __debug_level_3__
  fSuperkNuMuFluxCov->Print();
#endif

  // 10a->10d weights
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

  cout << "Exp/Beam >> Loading Kubo-san's norm err inputs!" << endl;

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

  cout << "Exp/Beam >> Done reading beam inputs!" << endl;

  return true;
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
  if(!fEnableFluxShapeTweaking) return 1;

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
void NuFlux::TweakShape(bool enable_tweaking)    
{ 
  fEnableFluxShapeTweaking = enable_tweaking;
  if(!fEnableFluxShapeTweaking) return;

  if(!fSuperkNuMuFluxHaveCholeskyL) {
    cout << "Exp/Beam/TweakShape >> Running Cholesky's method." << endl;
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

#ifdef __debug_level_3__
  cout << "Exp/Beam >> Tweaking SuperK flux shape" << endl;
#endif
 
  //
  //
  //

  if(fTweakShapeVldMode && fTweakShapeVldModeBin>0) {

    // fill in the tweaked histo
    //int nbins = fSuperkNuMuFlux10a->GetXaxis()->GetNbins();
    //for(int j=1; j <= nbins; j++){
      for(int j=1; j <= 200; j++){ ///TODO
         double Etrue = fSuperkNuMuFlux10a->GetBinCenter(j);
         double fixBinFlux10d = this->SuperkNuMuFlux(Etrue,true); // include 10a->10d weight
         int varbin_curr = util::FluxFix2VarBinNumber(j);
         if(varbin_curr == fTweakShapeVldModeBin) {
           fSuperkNuMuFluxTweaked->SetBinContent(j, fixBinFlux10d*(1+fTweakShapeVldModeDial));
         } else {
           fSuperkNuMuFluxTweaked->SetBinContent(j, fixBinFlux10d);
         }
      }//end of bin loop

      return;
  }

  //
  // Use Cholseky method to generate random fluxes
  //

  //get the variations
  TVectorD tweek_vec = 
     util::CholeskyGenerateCorrelatedParamVariations(*fSuperkNuMuFluxCovCholeskyL);

  // fill in the tweaked histo
//int nbins = fSuperkNuMuFlux10a->GetXaxis()->GetNbins();
//for(int j=1; j <= nbins; j++){
  for(int j=1; j <= 200; j++){ ///TODO
     double rto = util::FluxFractionsFixedBinRelToVarBin(
         j, fSuperkNuMuFlux10a, fSuperkNuMuFlux10aVarBin);
     double Etrue = fSuperkNuMuFlux10a->GetBinCenter(j);
     double fixBinFlux10d = this->SuperkNuMuFlux(Etrue,true); // include 10a->10d weight
     int varBin = util::FluxFix2VarBinNumber(j);
     fSuperkNuMuFluxTweaked->SetBinContent(j, fixBinFlux10d + tweek_vec[varBin - 1]*rto );
  }//end of bin loop

}
//________________________________________________________________________________________


