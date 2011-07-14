//________________________________________________________________________________________
/*
 Changes since v1.0.0:

 @ 16/03/2011 - CA
   Created from the old SuperkExpectation::ND280Inputs nested class. Removed the numu10a 
   namespace.
 @ 13/07/11 - CA
   Take into account the statistical error on the measured number of CC-like events.
   Update the efficiency error from +/- 5% to +5.6% -5.2%.
*/
//________________________________________________________________________________________

//#define __debug_level_1__
//#define __debug_level_2__
//#define __debug_level_3__

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TH1D.h>
#include <TH2D.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TMath.h>

#include "ND280Measurement.h"
#include "Utilities.h"
#include "NuXSec.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;

ClassImp(ND280Measurement)

//________________________________________________________________________________________
ND280Measurement::ND280Measurement()
{
  this->Init();
}
//________________________________________________________________________________________
ND280Measurement::~ND280Measurement()
{  
  this->CleanUp();
}
//________________________________________________________________________________________
void ND280Measurement::Init(void)
{
  fPOTData       = 0;
  fNumObsCCLike  = 0;
  fPOTMC         = 0;
  for(int i = 0; i < knmnd; i++) {
    fNumSelEt[i] = 0;
  }
  fSelEffTwkDial     = 0.;
  fSelEffRelErrPlus  = 0.;
  fSelEffRelErrMinus = 0.;
  fNObsCCLikeTwkDial = 0.;
  fRelStatErr        = 0.;
  fGenInputs         = 0;
     
  fNuMuEt    = 0;  
  fNuMuBarEt = 0;
  fNuEEt     = 0;
     
  fNomCorr   = 0;
}    
//________________________________________________________________________________________
void ND280Measurement::CleanUp(void)
{
  for(int i = 0; i < knmnd; i++) {
    if (fNumSelEt[i]) delete fNumSelEt[i];
  }
     
  if ( fNuMuEt    ) delete fNuMuEt;
  if ( fNuMuBarEt ) delete fNuMuBarEt;
  if ( fNuEEt     ) delete fNuEEt;
}
//________________________________________________________________________________________
bool ND280Measurement::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
  cout << "Exp/nd280 >> Reading nd280 inputs from input XML file." << endl;
    
  // Get Data POT
  fPOTData = util::XML2Dbl(xml,root_node,"nd280/pot_data");
  cout << "Exp/nd280 >> " << "N_{POT;data} : " << fPOTData << endl;
  
  // Get observed number of CC-like events
  fNumObsCCLike = util::XML2Dbl(xml,root_node,"nd280/saclay_selection/nobs_cclike");
  cout << "Exp/nd280 >> " << "N^{obs}_{CClike} : " << fNumObsCCLike << endl;
 
  // Get MC POT
  fPOTMC = util::XML2Dbl(xml,root_node,"nd280/saclay_selection/pot_mc");
  cout << "Exp/nd280 >> " << "N_{POT;mc} : " << fPOTMC << endl;

  // Get NMC_{mode}(Etrue) after cuts
  cout << "Exp/nd280 >> " << "Getting NMC_{mode}(Etrue) histograms" << endl;
  fNumSelEt [kNDmNuMuCCQE]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_ccqe");
  fNumSelEt [kNDmNuMuCCRES]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_ccres");
  fNumSelEt [kNDmNuMuCCCOH]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_cccoh");
  fNumSelEt [kNDmNuMuCCOth]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_ccoth");
  fNumSelEt [kNDmNuMuNCRES]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_ncres");
  fNumSelEt [kNDmNuMuNCCOH]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_nccoh");
  fNumSelEt [kNDmNuMuNCOth]     =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numu_ncoth");
  fNumSelEt [kNDmNuMuBarCCQE]   =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_ccqe");
  fNumSelEt [kNDmNuMuBarCCRES]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_ccres");
  fNumSelEt [kNDmNuMuBarCCCOH]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_cccoh");
  fNumSelEt [kNDmNuMuBarCCOth]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_ccoth");
  fNumSelEt [kNDmNuMuBarNCRES]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_ncres");
  fNumSelEt [kNDmNuMuBarNCCOH]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_nccoh");
  fNumSelEt [kNDmNuMuBarNCOth]  =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/numubar_ncoth");
  fNumSelEt [kNDmNuECCQE]       =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_ccqe");
  fNumSelEt [kNDmNuECCRES]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_ccres");
  fNumSelEt [kNDmNuECCCOH]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_cccoh");
  fNumSelEt [kNDmNuECCOth]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_ccoth");
  fNumSelEt [kNDmNuENCRES]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_ncres");
  fNumSelEt [kNDmNuENCCOH]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_nccoh");
  fNumSelEt [kNDmNuENCOth]      =
    util::XML2TH1F(xml,root_node,"nd280/saclay_selection/true_energy_spectrum/nue_ncoth");

  cout << "Exp/nd280 >> Done loading NMC_{mode}(Etrue) histograms" << endl;
 
  double data_to_mc_pot = fPOTData/fPOTMC;
  cout << "Exp/nd280 >> " << "MC POT scaling factor : " << data_to_mc_pot << endl;
   
  fNumSelEt [kNDmNuMuCCQE]     -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuCCRES]    -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuCCCOH]    -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuCCOth]    -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuNCRES]    -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuNCCOH]    -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuNCOth]    -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarCCQE]  -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarCCRES] -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarCCCOH] -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarCCOth] -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarNCRES] -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarNCCOH] -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuMuBarNCOth] -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuECCQE]      -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuECCRES]     -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuECCCOH]     -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuECCOth]     -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuENCRES]     -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuENCCOH]     -> Scale (data_to_mc_pot);
  fNumSelEt [kNDmNuENCOth]     -> Scale (data_to_mc_pot);

  cout << "Exp/nd280 >> Done rescaling NMC_{mode} histograms to the data POT" << endl;
  
  // Get error on efficiency
  fSelEffRelErrPlus  = util::XML2Dbl(xml,root_node,"nd280/saclay_selection/systematics/efficiency_relative_error/plus");
  fSelEffRelErrMinus = util::XML2Dbl(xml,root_node,"nd280/saclay_selection/systematics/efficiency_relative_error/minus");
  cout << "Exp/nd280 >> Efficiency error : "
       << "+ " << 100.*fSelEffRelErrPlus  << " %, "
       << "- " << 100.*fSelEffRelErrMinus << " %" << endl;

  // Get statistical error
  fRelStatErr = TMath::Sqrt(fNumObsCCLike)/fNumObsCCLike;
  cout << "Exp/nd280 >> Statistical error : " << 100.*fRelStatErr << " %" << endl;
    
  // Set colors and line styles for drawing
    
  fNumSelEt [kNDmNuMuCCQE]     -> SetLineColor (kBlack);
  fNumSelEt [kNDmNuMuCCRES]    -> SetLineColor (kBlack);
  fNumSelEt [kNDmNuMuCCCOH]    -> SetLineColor (kBlack);
  fNumSelEt [kNDmNuMuCCOth]    -> SetLineColor (kBlack);
  fNumSelEt [kNDmNuMuNCRES]    -> SetLineColor (kRed);
  fNumSelEt [kNDmNuMuNCCOH]    -> SetLineColor (kRed);
  fNumSelEt [kNDmNuMuNCOth]    -> SetLineColor (kRed);
  fNumSelEt [kNDmNuMuBarCCQE]  -> SetLineColor (kBlue);
  fNumSelEt [kNDmNuMuBarCCRES] -> SetLineColor (kBlue);
  fNumSelEt [kNDmNuMuBarCCCOH] -> SetLineColor (kBlue);
  fNumSelEt [kNDmNuMuBarCCOth] -> SetLineColor (kBlue);
  fNumSelEt [kNDmNuMuBarNCRES] -> SetLineColor (kMagenta);
  fNumSelEt [kNDmNuMuBarNCCOH] -> SetLineColor (kMagenta);
  fNumSelEt [kNDmNuMuBarNCOth] -> SetLineColor (kMagenta);
  fNumSelEt [kNDmNuECCQE]      -> SetLineColor (kGreen);
  fNumSelEt [kNDmNuECCRES]     -> SetLineColor (kGreen);
  fNumSelEt [kNDmNuECCCOH]     -> SetLineColor (kGreen);
  fNumSelEt [kNDmNuECCOth]     -> SetLineColor (kGreen);
  fNumSelEt [kNDmNuENCRES]     -> SetLineColor (kPink);
  fNumSelEt [kNDmNuENCCOH]     -> SetLineColor (kPink);
  fNumSelEt [kNDmNuENCOth]     -> SetLineColor (kPink);

  fNumSelEt [kNDmNuMuCCQE]     -> SetLineStyle (kSolid);
  fNumSelEt [kNDmNuMuCCRES]    -> SetLineStyle (kDashed);
  fNumSelEt [kNDmNuMuCCCOH]    -> SetLineStyle (kDotted);
  fNumSelEt [kNDmNuMuCCOth]    -> SetLineStyle (kDashDotted);
  fNumSelEt [kNDmNuMuNCRES]    -> SetLineStyle (kDashed);
  fNumSelEt [kNDmNuMuNCCOH]    -> SetLineStyle (kDotted);
  fNumSelEt [kNDmNuMuNCOth]    -> SetLineStyle (kDashDotted);
  fNumSelEt [kNDmNuMuBarCCQE]  -> SetLineStyle (kSolid);
  fNumSelEt [kNDmNuMuBarCCRES] -> SetLineStyle (kDashed);
  fNumSelEt [kNDmNuMuBarCCCOH] -> SetLineStyle (kDotted);
  fNumSelEt [kNDmNuMuBarCCOth] -> SetLineStyle (kDashDotted);
  fNumSelEt [kNDmNuMuBarNCRES] -> SetLineStyle (kDashed);
  fNumSelEt [kNDmNuMuBarNCCOH] -> SetLineStyle (kDotted);
  fNumSelEt [kNDmNuMuBarNCOth] -> SetLineStyle (kDashDotted);
  fNumSelEt [kNDmNuECCQE]      -> SetLineStyle (kSolid);
  fNumSelEt [kNDmNuECCRES]     -> SetLineStyle (kDashed);
  fNumSelEt [kNDmNuECCCOH]     -> SetLineStyle (kDotted);
  fNumSelEt [kNDmNuECCOth]     -> SetLineStyle (kDashDotted);
  fNumSelEt [kNDmNuENCRES]     -> SetLineStyle (kDashed);
  fNumSelEt [kNDmNuENCCOH]     -> SetLineStyle (kDotted);
  fNumSelEt [kNDmNuENCOth]     -> SetLineStyle (kDashDotted);
  
  //
  // histograms needed for beam norm err calculation
  //
  
  fNuMuEt    = new TH1D("fNuMuEt",    "Etrue for numu",    200, 0, 10); /// binning expected by Kubo's code
  fNuMuBarEt = new TH1D("fNuMuBarEt", "Etrue for numubar", 200, 0, 10);
  fNuEEt     = new TH1D("fNuMuEEt",   "Etrue for nue",     200, 0, 10);
  
  cout << "Exp/nd280 >> Done reading nd280 inputs!" << endl;
  
  //
  // store nominal correction factor
  //
  
  bool nominal = true;
  double nobscclike_nd280  = this -> NObsCCLike(nominal);
  double nexpncbkg_nd280   = this -> NExpNCBkg (nominal);
  double nexpcc_nd280      = this -> NExpCC    (nominal);
  fNomCorr = nobscclike_nd280 / (nexpncbkg_nd280 + nexpcc_nd280);
  cout << "Exp/nd280 >> Nominal correction factor from nd280: " << fNomCorr << endl;
  
  return true;
}
//________________________________________________________________________________________
void ND280Measurement::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  util::SetDefaultStyle();

  TString opt = option;

  double Emax = 10000;
  double Nmax = 100;
  double Nmin = 0.001;

  if(opt.Contains("Saclay/EtrueSpectrum/NuMu/CC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,Nmin,Emax,Nmax);
    hframe -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hframe -> GetYaxis()->SetTitle(Form("events / %E POT", fPOTData));
    hframe->Draw();
    fNumSelEt [kNDmNuMuCCQE]     -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuCCRES]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuCCCOH]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuCCOth]    -> Draw("LPERRSAME");
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fNumSelEt [kNDmNuMuCCQE] , "#nu_{#mu} CCQE",    "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuCCRES], "#nu_{#mu} CCRES",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuCCCOH], "#nu_{#mu} CCCOH",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuCCOth], "#nu_{#mu} CCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/EtrueSpectrum/NuMu/NC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,Nmin,Emax,Nmax);
    hframe->Draw();
    fNumSelEt [kNDmNuMuNCRES]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuNCCOH]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuNCOth]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuNCRES]    -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fNumSelEt [kNDmNuMuNCRES]    -> GetYaxis()->SetTitle(Form("events / %E POT", fPOTData));
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fNumSelEt [kNDmNuMuNCRES], "#nu_{#mu} NCRES",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuNCCOH], "#nu_{#mu} NCCOH",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuNCOth], "#nu_{#mu} NCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/EtrueSpectrum/NuMuBar/CC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,Nmin,Emax,Nmax);
    hframe->Draw();
    fNumSelEt [kNDmNuMuBarCCQE]     -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarCCRES]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarCCCOH]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarCCOth]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarCCQE]     -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fNumSelEt [kNDmNuMuBarCCQE]     -> GetYaxis()->SetTitle(Form("events / %E POT", fPOTData));
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fNumSelEt [kNDmNuMuBarCCQE] , "#bar{#nu_{#mu}} CCQE",    "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuBarCCRES], "#bar{#nu_{#mu}} CCRES",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuBarCCCOH], "#bar{#nu_{#mu}} CCCOH",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuBarCCOth], "#bar{#nu_{#mu}} CCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/EtrueSpectrum/NuMuBar/NC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,Nmin,Emax,Nmax);
    hframe->Draw();
    fNumSelEt [kNDmNuMuBarNCRES]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarNCCOH]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarNCOth]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuMuBarNCRES]    -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fNumSelEt [kNDmNuMuBarNCRES]    -> GetYaxis()->SetTitle(Form("events / %E POT", fPOTData));
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fNumSelEt [kNDmNuMuBarNCRES], "#bar{#nu_{#mu}} NCRES",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuBarNCCOH], "#bar{#nu_{#mu}} NCCOH",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuMuBarNCOth], "#bar{#nu_{#mu}} NCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/EtrueSpectrum/NuE/CC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,Nmin,Emax,Nmax);
    hframe->Draw();
    fNumSelEt [kNDmNuECCQE]     -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuECCRES]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuECCCOH]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuECCOth]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuECCQE]     -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fNumSelEt [kNDmNuECCQE]     -> GetYaxis()->SetTitle(Form("events / %E POT", fPOTData));
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fNumSelEt [kNDmNuECCQE] , "#nu_{e} CCQE",    "LP");
    l->AddEntry(fNumSelEt [kNDmNuECCRES], "#nu_{e} CCRES",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuECCCOH], "#nu_{e} CCCOH",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuECCOth], "#nu_{e} CCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
  else
  if(opt.Contains("Saclay/EtrueSpectrum/NuE/NC")) { 
    TH1F * hframe = (TH1F*)c->DrawFrame(0,Nmin,Emax,Nmax);
    hframe->Draw();
    fNumSelEt [kNDmNuENCRES]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuENCCOH]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuENCOth]    -> Draw("LPERRSAME");
    fNumSelEt [kNDmNuENCRES]    -> GetXaxis()->SetTitle("E_{#nu} (GeV)");
    fNumSelEt [kNDmNuENCRES]    -> GetYaxis()->SetTitle(Form("events / %E POT", fPOTData));
    TLegend * l = new TLegend(0.7,0.6,0.9,0.9);
    l->SetFillColor(0);
    l->AddEntry(fNumSelEt [kNDmNuENCRES], "#nu_{e} NCRES",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuENCCOH], "#nu_{e} NCCOH",   "LP");
    l->AddEntry(fNumSelEt [kNDmNuENCOth], "#nu_{e} NCOther", "LP");
    l->Draw();
    c->SetLogy();
  }
	 
  c->Update();
}
//________________________________________________________________________________________
double ND280Measurement::CorrectionFactor(void) const
{
  //
  // avoid looping over all histograms and histogram bins if not needed
  //
  bool xsec_tweaked = fGenInputs->TweakAffectsND280();
  bool eff_tweaked  = (TMath::Abs(fSelEffTwkDial) > kASmallNum);
  bool norm_tweaked = (TMath::Abs(fNObsCCLikeTwkDial) > kASmallNum);
  bool tweaked      = eff_tweaked || xsec_tweaked || norm_tweaked;

  if(!tweaked) {
#ifdef __debug_level_2__
    cout << "Exp / nd280 >> nd280 event spectra not tweaked."
         << "Returning nominal value = " << fNomCorr << endl;
#endif
    return fNomCorr;
  }

  //
  // calculate correction factor
  //
  double nobscclike_nd280  = this -> NObsCCLike();
  double nexpncbkg_nd280   = this -> NExpNCBkg();
  double nexpcc_nd280      = this -> NExpCC();

  double correction = nobscclike_nd280 / (nexpncbkg_nd280 + nexpcc_nd280);

#ifdef __debug_level_2__
  cout << "Exp / nd280 >> "
       << "nobscclike_nd280 = "  << nobscclike_nd280
       << ", nexpncbkg_nd280 = " << nexpncbkg_nd280
       << ", nexpcc_nd280 = " << nexpcc_nd280
       << " ==> correction factor = " << correction << endl;
#endif

  return correction;
}
//________________________________________________________________________________________
double ND280Measurement::NObsCCLike(bool nominal) const 
{
// observed number of CC-like events

  if(nominal) return fNumObsCCLike;

  double n = fNumObsCCLike * (1. + fNObsCCLikeTwkDial * fRelStatErr);
  return n;
}
//________________________________________________________________________________________
double ND280Measurement::NExpNCBkg(bool nominal) const
{
// expected NC contamination of CC-like sample

  double nexpncbkg = 
       this->NumOfEvents( kNDmNuMuNCRES,    nominal) +
       this->NumOfEvents( kNDmNuMuNCCOH,    nominal) +
       this->NumOfEvents( kNDmNuMuNCOth,    nominal) +
       this->NumOfEvents( kNDmNuMuBarNCRES, nominal) +
       this->NumOfEvents( kNDmNuMuBarNCCOH, nominal) +
       this->NumOfEvents( kNDmNuMuBarNCOth, nominal) +
       this->NumOfEvents( kNDmNuENCRES,     nominal) +
       this->NumOfEvents( kNDmNuENCCOH,     nominal) +
       this->NumOfEvents( kNDmNuENCOth,     nominal);

  return nexpncbkg;
}
//________________________________________________________________________________________
double ND280Measurement::NExpCC(bool nominal) const
{
// expected number of CC events

  double nexpcc = 
       this->NumOfEvents( kNDmNuMuCCQE,     nominal) +
       this->NumOfEvents( kNDmNuMuCCRES,    nominal) +
       this->NumOfEvents( kNDmNuMuCCCOH,    nominal) +
       this->NumOfEvents( kNDmNuMuCCOth,    nominal) +
       this->NumOfEvents( kNDmNuMuBarCCQE,  nominal) +
       this->NumOfEvents( kNDmNuMuBarCCRES, nominal) +
       this->NumOfEvents( kNDmNuMuBarCCCOH, nominal) +
       this->NumOfEvents( kNDmNuMuBarCCOth, nominal) +
       this->NumOfEvents( kNDmNuECCQE,      nominal) +
       this->NumOfEvents( kNDmNuECCRES,     nominal) +
       this->NumOfEvents( kNDmNuECCCOH,     nominal) +
       this->NumOfEvents( kNDmNuECCOth,     nominal);

  return nexpcc;
}
//________________________________________________________________________________________
double ND280Measurement::NumOfEvents(ND280Mode_t mode, bool nominal) const
{
  TH1F * spectrum = fNumSelEt[mode];
  assert(spectrum);
  double sum_nev = 0;
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
     sum_nev += this->NumOfEvents(mode,ibin, nominal);
  }
  return sum_nev;
}
//________________________________________________________________________________________
double ND280Measurement::NumOfEvents(
       t2kana::conv::ND280Mode_t mode, int ibin, bool nominal) const
{
  TH1F * spectrum = fNumSelEt[mode];
  assert(spectrum);
  double mev2gev = 1E-3;
  double Etrue = spectrum->GetBinCenter(ibin) * mev2gev; // GeV

  double nev = spectrum->GetBinContent(ibin);
  if(nominal) return nev;

  double efftwkdial = fSelEffTwkDial;
  double efferr     = this->EfficiencyRelErr(Etrue, util::Sign(efftwkdial));
  double effwght    = 1 + efftwkdial * efferr;
  double sigwght    = fGenInputs->ND280XSecRatioWght(Etrue, mode);

  nev *= (effwght*sigwght);
  return nev;
}
//________________________________________________________________________________________
double ND280Measurement::NumOfEvents(
       t2kana::conv::ND280Mode_t mode, double Etrue, bool nominal) const
{
  TH1F * spectrum   = fNumSelEt[mode];
  assert(spectrum);
  double mev2gev = 1E-3;
  int ibin = spectrum->FindBin(Etrue/mev2gev);
  return this->NumOfEvents(mode,ibin, nominal);
}
//________________________________________________________________________________________
void ND280Measurement::UpdateEtrueSpectra(void)
{
  TH1D * spectrum = 0;

  // update numu
  spectrum = fNuMuEt;
  assert(spectrum);
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
    double nev = 
         this->NumOfEvents( kNDmNuMuCCQE,  ibin) +
         this->NumOfEvents( kNDmNuMuCCRES, ibin) +
         this->NumOfEvents( kNDmNuMuCCCOH, ibin) +
         this->NumOfEvents( kNDmNuMuCCOth, ibin) +
         this->NumOfEvents( kNDmNuMuNCRES, ibin) +
         this->NumOfEvents( kNDmNuMuNCCOH, ibin) +
         this->NumOfEvents( kNDmNuMuNCOth, ibin);
    spectrum->SetBinContent(ibin,nev);
  }//bins

  // update numubar
  spectrum = fNuMuBarEt;
  assert(spectrum);
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
    double nev = 
         this->NumOfEvents( kNDmNuMuBarCCQE,  ibin) +
         this->NumOfEvents( kNDmNuMuBarCCRES, ibin) +
         this->NumOfEvents( kNDmNuMuBarCCCOH, ibin) +
         this->NumOfEvents( kNDmNuMuBarCCOth, ibin) +
         this->NumOfEvents( kNDmNuMuBarNCRES, ibin) +
         this->NumOfEvents( kNDmNuMuBarNCCOH, ibin) +
         this->NumOfEvents( kNDmNuMuBarNCOth, ibin);
    spectrum->SetBinContent(ibin,nev);
  }//bins

  // update nue
  spectrum = fNuEEt;
  assert(spectrum);
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
    double nev = 
         this->NumOfEvents( kNDmNuECCQE,  ibin) +
         this->NumOfEvents( kNDmNuECCRES, ibin) +
         this->NumOfEvents( kNDmNuECCCOH, ibin) +
         this->NumOfEvents( kNDmNuECCOth, ibin) +
         this->NumOfEvents( kNDmNuENCRES, ibin) +
         this->NumOfEvents( kNDmNuENCCOH, ibin) +
         this->NumOfEvents( kNDmNuENCOth, ibin);
    spectrum->SetBinContent(ibin,nev);
  }//bins

}
//________________________________________________________________________________________
double ND280Measurement::EfficiencyRelErr(double /*Et*/, int sign) const
{
// relative uncertainty on efficiency 
// (can be asymmetric)

  if(sign > 0) return fSelEffRelErrPlus;
  if(sign < 0) return fSelEffRelErrMinus;
  return 0;
}
//________________________________________________________________________________________
void ND280Measurement::TweakEfficiency(double twkdial)
{
// tweak efficiency of all modes by the same amount
//
  fSelEffTwkDial = twkdial;
}
//________________________________________________________________________________________
void ND280Measurement::FluctuateNObsCCLike(double twkdial)
{
  fNObsCCLikeTwkDial = twkdial;
}
//________________________________________________________________________________________

