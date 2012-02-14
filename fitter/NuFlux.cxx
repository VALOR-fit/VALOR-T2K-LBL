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
#define __show_init_mesg_debug__   

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TVectorD.h>
#include <TVirtualPad.h>
#include <TLegend.h>

#include "NuFlux.h" 
#include "Conventions.h" 
#include "Utilities.h" 
#include "BANFF.h" 

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
  fSuperkFluxPOT               = 0;
  fSuperkNuMuFluxNom           = 0; 
  fSuperkNuMuBarFluxNom        = 0; 
  fSuperkNuEFluxNom            = 0; 
  fSuperkNuMuFluxTuneWght      = 0;
  fSuperkNuMuBarFluxTuneWght   = 0;
  fSuperkNuEFluxTuneWght       = 0;
  fBANFF                       = 0;
}
//________________________________________________________________________________________
void NuFlux::CleanUp(void)
{
  if ( fSuperkNuMuFluxNom          ) delete fSuperkNuMuFluxNom;
  if ( fSuperkNuMuBarFluxNom       ) delete fSuperkNuMuBarFluxNom;
  if ( fSuperkNuEFluxNom           ) delete fSuperkNuEFluxNom;
  if ( fSuperkNuMuFluxTuneWght     ) delete fSuperkNuMuFluxTuneWght;
  if ( fSuperkNuMuBarFluxTuneWght  ) delete fSuperkNuMuBarFluxTuneWght;
  if ( fSuperkNuEFluxTuneWght      ) delete fSuperkNuEFluxTuneWght;
}
//________________________________________________________________________________________
bool NuFlux::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read beam inputs to the oscillation fit
//
#ifdef __show_init_mesg_info__
  cout << "Exp/Beam >> Reading beam inputs from input XML file." << endl;
#endif

  //
  // Integrated POT for flux simulation
  //

  fSuperkFluxPOT = util::XML2Dbl(xml,root_node,"beam/superk_flux/pot_mc");
#ifdef __show_init_mesg_info__
  cout << "Exp/Beam >> POT_{SK flux simulation} = " << fSuperkFluxPOT << endl;
#endif

  //
  // Nominal fluxes
  //

  // SuperK numu flux
  fSuperkNuMuFluxNom = util::XML2TH1D(xml,root_node,"beam/superk_flux/nominal/numu");
  if(!fSuperkNuMuFluxNom) {
    cerr << "Exp/Beam >> No input SuperK numu flux histogram!" << endl;
    exit(1);
  }
  // SuperK numubar flux
  fSuperkNuMuBarFluxNom = util::XML2TH1D(xml,root_node,"beam/superk_flux/nominal/numubar");
  if(!fSuperkNuMuBarFluxNom) {
    cerr << "Exp/Beam >> No input SuperK numubar flux histogram!" << endl;
    exit(1);
  }
  // SuperK nue flux
  fSuperkNuEFluxNom = util::XML2TH1D(xml,root_node,"beam/superk_flux/nominal/nue");
  if(!fSuperkNuEFluxNom) {
    cerr << "Exp/Beam >> No input SuperK nue flux histogram!" << endl;
    exit(1);
  }

  //
  // Nominal -> tuned flux weights
  //

  // SuperK numu flux tuning weights
  fSuperkNuMuFluxTuneWght = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/tuning_weights/numu");
  if(!fSuperkNuMuFluxTuneWght) {
    cerr << "Exp/Beam >> No input SuperK numu flux tuning weights!" << endl;
    exit(1);
  }
  // SuperK numubar flux tuning weights
  fSuperkNuMuBarFluxTuneWght = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/tuning_weights/numubar");
  if(!fSuperkNuMuBarFluxTuneWght) {
    cerr << "Exp/Beam >> No input SuperK numubar flux tuning weights!" << endl;
    exit(1);
  }
  // SuperK nue flux tuning weights
  fSuperkNuEFluxTuneWght = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/tuning_weights/nue");
  if(!fSuperkNuEFluxTuneWght) {
    cerr << "Exp/Beam >> No input SuperK nue flux tuning weights!" << endl;
    exit(1);
  }

#ifdef __show_init_mesg_info__
  cout << "Exp/Beam >> Done reading beam inputs!" << endl;
#endif

  return true;
}
//________________________________________________________________________________________
void NuFlux::TellMeAboutTheBANFFfit(BANFF * banff)
{
  fBANFF = banff;
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxPOT (void) const
{
  return fSuperkFluxPOT;
}
//________________________________________________________________________________________
double NuFlux::SuperkNuMuFlux(double Etrue, bool tuned) const
{
  return this->Flux(Etrue,tuned,fSuperkNuMuFluxNom,fSuperkNuMuFluxTuneWght);
}
//________________________________________________________________________________________
double NuFlux::SuperkNuMuBarFlux(double Etrue, bool tuned) const
{
  return this->Flux(Etrue,tuned,fSuperkNuMuBarFluxNom,fSuperkNuMuBarFluxTuneWght);
}
//________________________________________________________________________________________
double NuFlux::SuperkNuEFlux(double Etrue, bool tuned) const
{
  return this->Flux(Etrue,tuned,fSuperkNuEFluxNom,fSuperkNuEFluxTuneWght);
}
//________________________________________________________________________________________
double NuFlux::SuperkNuMuFluxTuneWght(double Etrue) const
{
  return this->FluxTuneWght(Etrue,fSuperkNuMuFluxTuneWght);
}
//________________________________________________________________________________________
double NuFlux::SuperkNuMuBarFluxTuneWght(double Etrue) const
{
  return this->FluxTuneWght(Etrue,fSuperkNuMuBarFluxTuneWght);
}
//________________________________________________________________________________________
double NuFlux::SuperkNuEFluxTuneWght(double Etrue) const
{
  return this->FluxTuneWght(Etrue,fSuperkNuEFluxTuneWght);
}
//________________________________________________________________________________________
double NuFlux::SuperkNuMuFluxWght(double Etrue) const
{
  return 1.;
}
//________________________________________________________________________________________
double NuFlux::SuperkNuMuBarFluxWght(double Etrue) const
{
  return 1.;
}
//________________________________________________________________________________________
double NuFlux::SuperkNuEFluxWght(double Etrue) const
{
  return 1.;
}
//________________________________________________________________________________________
double NuFlux::Flux(
  double Etrue, bool tuned, TH1D * nominal_flux, TH1D * tuning_weights) const
{
  assert(nominal_flux);
  int flux_bin = nominal_flux->FindBin(Etrue);
  double flux  = nominal_flux->GetBinContent(flux_bin);
  if(tuned) {
    assert(tuning_weights);
    double wght = this->FluxTuneWght(Etrue,tuning_weights);
    flux *= wght;
  }
  return flux;
}
//________________________________________________________________________________________
double NuFlux::FluxTuneWght(
  double Etrue, TH1D * tuning_weights) const
{
  double wght = 1.0;
  if(Etrue < kEtrueMaxForFluxTuneWeights) {
    assert(tuning_weights);
    int wght_bin = tuning_weights->FindBin(Etrue);
    wght = tuning_weights->GetBinContent(wght_bin);
  }
  return wght;
}
//________________________________________________________________________________________
void NuFlux::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  if(opt.Contains("Flux/SK/NuMu/TuneWghts")) { 
    TH1D * h = new TH1D(*fSuperkNuMuFluxTuneWght);
    hnom->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK numu flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMuBar/TuneWghts")) { 
    TH1D * h = new TH1D(*fSuperkNuMuBarFluxTuneWght);
    hnom->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK numubar flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuE/TuneWghts")) { 
    TH1D * h = new TH1D(*fSuperkNuEFluxTuneWght);
    hnom->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK nue flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMu")) { 
    TH1D * hnom = new TH1D(*fSuperkNuMuFlux);
    hnom->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK numu flux (#neutrinos /cm2 /POT)");
    hnom->Divide(fSuperkFluxPOT);
    hnom->Draw();
    TH1D * htun = new TH1D(*hnom);
    htun->SetDirectory(0);
    for(int i = 1; i <= htun->GetNbinsX(); i++) {
	double E = htun->GetBinCenter(i);
	double f = this->SuperkNuMuFlux(E,true);
        htun->SetBinContent(i,f);
    }
    htun->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMuBar")) { 
    TH1D * hnom = new TH1D(*fSuperkNuMuBarFlux);
    hnom->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK numubar flux (#neutrinos /cm2 /POT)");
    hnom->Divide(fSuperkFluxPOT);
    hnom->Draw();
    TH1D * htun = new TH1D(*hnom);
    htun->SetDirectory(0);
    for(int i = 1; i <= htun->GetNbinsX(); i++) {
	double E = htun->GetBinCenter(i);
	double f = this->SuperkNuMuBarFlux(E,true);
        htun->SetBinContent(i,f);
    }
    htun->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuE")) { 
    TH1D * hnom = new TH1D(*fSuperkNuEFlux);
    hnom->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK nue flux (#neutrinos /cm2 /POT)");
    hnom->Divide(fSuperkFluxPOT);
    hnom->Draw();
    TH1D * htun = new TH1D(*hnom);
    htun->SetDirectory(0);
    for(int i = 1; i <= htun->GetNbinsX(); i++) {
	double E = htun->GetBinCenter(i);
	double f = this->SuperkNuEFlux(E,true);
        htun->SetBinContent(i,f);
    }
    htun->Draw("same");
    c->Update();
  }
}
//________________________________________________________________________________________
