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
  fSuperkNuBarFluxNom          = 0; 
  fSuperkNuMuFluxTuneWght      = 0;
  fSuperkNuMuBarFluxTuneWght   = 0;
  fSuperkNuEFluxTuneWght       = 0;
  fSuperkNuEBarFluxTuneWght    = 0;
  fBANFF                       = 0;
}
//________________________________________________________________________________________
void NuFlux::CleanUp(void)
{
  if ( fSuperkNuMuFluxNom          ) delete fSuperkNuMuFluxNom;
  if ( fSuperkNuMuBarFluxNom       ) delete fSuperkNuMuBarFluxNom;
  if ( fSuperkNuEFluxNom           ) delete fSuperkNuEFluxNom;
  if ( fSuperkNuEBarFluxNom        ) delete fSuperkNuEBarFluxNom;
  if ( fSuperkNuMuFluxTuneWght     ) delete fSuperkNuMuFluxTuneWght;
  if ( fSuperkNuMuBarFluxTuneWght  ) delete fSuperkNuMuBarFluxTuneWght;
  if ( fSuperkNuEFluxTuneWght      ) delete fSuperkNuEFluxTuneWght;
  if ( fSuperkNuEBarFluxTuneWght   ) delete fSuperkNuEBarFluxTuneWght;
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
  // SuperK nuebar flux
  fSuperkNuEBarFluxNom = util::XML2TH1D(xml,root_node,"beam/superk_flux/nominal/nuebar");
  if(!fSuperkNuEBarFluxNom) {
    cerr << "Exp/Beam >> No input SuperK nuebar flux histogram!" << endl;
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
  // SuperK nuebar flux tuning weights
  fSuperkNuEBarFluxTuneWght = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/tuning_weights/nuebar");
  if(!fSuperkNuEBarFluxTuneWght) {
    cerr << "Exp/Beam >> No input SuperK nuebar flux tuning weights!" << endl;
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
double NuFlux::SuperkFlux(int nu_pdg, double Etrue, bool tuned) const
{
  const TH1D * nom_flux_spectrum = this->SuperkNominalFluxSpectrum(nu_pdg);
  assert(nom_flux_spectrum);

  int bin = nom_flux_spectrum->FindBin(Etrue);
  double flux = nom_flux_spectrum->GetBinContent(bin);

  if(tuned) {
     double wght = this->SuperkFluxTuneWght(nu_pdg,Etrue);
     flux *= wght;
  }

  return flux;
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxTuneWght(int nu_pdg, double Etrue) const
{
  if(Etrue > kEtrueMaxForFluxTuneWeights) 
  {
    return 1;
  }

  const TH1D * wght_spectrum = this->SuperkFluxTuningWghtSpectrum(nu_pdg);
  assert(wght_spectrum);

  int bin = wght_spectrum->FindBin(Etrue);
  double wght = wght_spectrum->GetBinContent(bin);
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxWght(int nu_pdg, double Etrue) const
{
  return 1.;
}
//________________________________________________________________________________________
const TH1D * NuFlux::SuperkNominalFluxSpectrum(nu_pdg) const
{
  if      ( nu_pdg == kPdgNuE     ) return fSuperkNuEFluxNom;
  else if ( nu_pdg == kPdgNuEBar  ) return fSuperkNuEBarFluxNom;
  else if ( nu_pdg == kPdgNuMu    ) return fSuperkNuMuFluxNom;
  else if ( nu_pdg == kPdgNuMuBar ) return fSuperkNuMuBarFluxNom;
  else
  {
    exit(1);
  }
  return 0;
}
//________________________________________________________________________________________
const TH1D * NuFlux::SuperkFluxTuningWghtSpectrum(nu_pdg) const
{
  if      ( nu_pdg == kPdgNuE     ) return fSuperkNuEFluxTuneWght;
  else if ( nu_pdg == kPdgNuEBar  ) return fSuperkNuEBarFluxTuneWght;
  else if ( nu_pdg == kPdgNuMu    ) return fSuperkNuMuFluxTuneWght;
  else if ( nu_pdg == kPdgNuMuBar ) return fSuperkNuMuBarFluxTuneWght;
  else
  {
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
  if(opt.Contains("Flux/SK/NuEBar/TuneWghts")) { 
    TH1D * h = new TH1D(*fSuperkNuEBarFluxTuneWght);
    hnom->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK nuebar flux weight");
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
	double f = this->SuperFlux(kPdgNuMu,E,true);
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
	double f = this->SuperkFlux(kPdgNuMuBar,E,true);
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
	double f = this->SuperkFlux(kPdgNuE,E,true);
        htun->SetBinContent(i,f);
    }
    htun->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuEBar")) { 
    TH1D * hnom = new TH1D(*fSuperkNuEBarFlux);
    hnom->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK nuebar flux (#neutrinos /cm2 /POT)");
    hnom->Divide(fSuperkFluxPOT);
    hnom->Draw();
    TH1D * htun = new TH1D(*hnom);
    htun->SetDirectory(0);
    for(int i = 1; i <= htun->GetNbinsX(); i++) {
	double E = htun->GetBinCenter(i);
	double f = this->SuperkFlux(kPdgNuEBar,E,true);
        htun->SetBinContent(i,f);
    }
    htun->Draw("same");
    c->Update();
  }
}
//________________________________________________________________________________________
