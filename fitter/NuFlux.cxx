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

#include <TVirtualPad.h>
#include <TLegend.h>

#include "NuFlux.h" 
#include "BANFF.h" 
#include "Conventions.h" 
#include "Utilities.h" 

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
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
  fSuperkNuEBarFluxNom         = 0; 
  fSuperkNuMuFluxTun           = 0; 
  fSuperkNuMuBarFluxTun        = 0; 
  fSuperkNuEFluxTun            = 0; 
  fSuperkNuEBarFluxTun         = 0; 
  fSuperkNuMuFluxTunWght       = 0;
  fSuperkNuMuBarFluxTunWght    = 0;
  fSuperkNuEFluxTunWght        = 0;
  fSuperkNuEBarFluxTunWght     = 0;
  fBANFF                       = 0;
}
//________________________________________________________________________________________
void NuFlux::CleanUp(void)
{
  if ( fSuperkNuMuFluxNom          ) delete fSuperkNuMuFluxNom;
  if ( fSuperkNuMuBarFluxNom       ) delete fSuperkNuMuBarFluxNom;
  if ( fSuperkNuEFluxNom           ) delete fSuperkNuEFluxNom;
  if ( fSuperkNuEBarFluxNom        ) delete fSuperkNuEBarFluxNom;
  if ( fSuperkNuMuFluxTun          ) delete fSuperkNuMuFluxTun;
  if ( fSuperkNuMuBarFluxTun       ) delete fSuperkNuMuBarFluxTun;
  if ( fSuperkNuEFluxTun           ) delete fSuperkNuEFluxTun;
  if ( fSuperkNuEBarFluxTun        ) delete fSuperkNuEBarFluxTun;
  if ( fSuperkNuMuFluxTunWght      ) delete fSuperkNuMuFluxTunWght;
  if ( fSuperkNuMuBarFluxTunWght   ) delete fSuperkNuMuBarFluxTunWght;
  if ( fSuperkNuEFluxTunWght       ) delete fSuperkNuEFluxTunWght;
  if ( fSuperkNuEBarFluxTunWght    ) delete fSuperkNuEBarFluxTunWght;
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
  // Load nominal fluxes
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
  // Load nominal -> tuned flux weights
  //

  // SuperK numu flux tuning weights
  fSuperkNuMuFluxTunWght = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/tuning_weights/numu");
  if(!fSuperkNuMuFluxTunWght) {
    cerr << "Exp/Beam >> No input SuperK numu flux tuning weights!" << endl;
    exit(1);
  }
  // SuperK numubar flux tuning weights
  fSuperkNuMuBarFluxTunWght = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/tuning_weights/numubar");
  if(!fSuperkNuMuBarFluxTunWght) {
    cerr << "Exp/Beam >> No input SuperK numubar flux tuning weights!" << endl;
    exit(1);
  }
  // SuperK nue flux tuning weights
  fSuperkNuEFluxTunWght = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/tuning_weights/nue");
  if(!fSuperkNuEFluxTunWght) {
    cerr << "Exp/Beam >> No input SuperK nue flux tuning weights!" << endl;
    exit(1);
  }
  // SuperK nuebar flux tuning weights
  fSuperkNuEBarFluxTunWght = 
     util::XML2TH1D(xml,root_node,"beam/superk_flux/tuning_weights/nuebar");
  if(!fSuperkNuEBarFluxTunWght) {
    cerr << "Exp/Beam >> No input SuperK nuebar flux tuning weights!" << endl;
    exit(1);
  }

  //
  // Create / init (to nominal) tuned flux histograms & then
  // apply tuning weights
  //
/*
  fSuperkNuMuFluxTun    = new TH1D(*fSuperkNuMuFluxNom);  
  fSuperkNuMuBarFluxTun = new TH1D(*fSuperkNuMuBarFluxNom); 
  fSuperkNuEFluxTun     = new TH1D(*fSuperkNuEFluxNom);     
  fSuperkNuEBarFluxTun  = new TH1D(*fSuperkNuEBarFluxNom);          

  fSuperkNuMuFluxTun    -> SetDirectory(0);
  fSuperkNuMuBarFluxTun -> SetDirectory(0); 
  fSuperkNuEFluxTun     -> SetDirectory(0);
  fSuperkNuEBarFluxTun  -> SetDirectory(0);
*/
  this->ApplyFluxTuneWghts(conv::kPdgNuE    );
  this->ApplyFluxTuneWghts(conv::kPdgNuEBar );
  this->ApplyFluxTuneWghts(conv::kPdgNuMu   );
  this->ApplyFluxTuneWghts(conv::kPdgNuMuBar);

#ifdef __show_init_mesg_info__
  cout << "Exp/Beam >> Done reading beam inputs!" << endl;
#endif

  return true;
}
//________________________________________________________________________________________
void NuFlux::ApplyFluxTuneWghts(int nu_pdg)
{
  TH1D * hnom = this->SuperkFluxSpectrum(nu_pdg, false);
  TH1D * htun = this->SuperkFluxSpectrum(nu_pdg, true );
  htun = new TH1D(*hnom);
  htun -> SetDirectory(0);
  for(int i = 1; i <= htun->GetNbinsX(); i++) {
     double E = htun->GetBinCenter (i);
     double f = htun->GetBinContent(i);
     double w = this->SuperkFluxTuneWght(nu_pdg,E);
     f *= w;
     htun->SetBinContent(i,f);
  }
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
  TH1D * flux_spectrum = this->SuperkFluxSpectrum(nu_pdg,tuned);
  assert(flux_spectrum);

  int bin = flux_spectrum->FindBin(Etrue);
  double flux = flux_spectrum->GetBinContent(bin);

  return flux;
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxTuneWght(int nu_pdg, double Etrue) const
{
  if(Etrue > conv::kEtrueMaxForFluxTuneWeights) 
  {
    return 1;
  }

  TH1D * wght_spectrum = this->SuperkFluxTuningWghts(nu_pdg);
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
TH1D * NuFlux::SuperkFluxSpectrum(int nu_pdg, bool tuned) const
{
  if(tuned) 
  {
     if      ( nu_pdg == conv::kPdgNuE     ) return fSuperkNuEFluxTun;
     else if ( nu_pdg == conv::kPdgNuEBar  ) return fSuperkNuEBarFluxTun;
     else if ( nu_pdg == conv::kPdgNuMu    ) return fSuperkNuMuFluxTun;
     else if ( nu_pdg == conv::kPdgNuMuBar ) return fSuperkNuMuBarFluxTun;
     else
     {
       exit(1);
     }
  }
  else
  {
     if      ( nu_pdg == conv::kPdgNuE     ) return fSuperkNuEFluxNom;
     else if ( nu_pdg == conv::kPdgNuEBar  ) return fSuperkNuEBarFluxNom;
     else if ( nu_pdg == conv::kPdgNuMu    ) return fSuperkNuMuFluxNom;
     else if ( nu_pdg == conv::kPdgNuMuBar ) return fSuperkNuMuBarFluxNom;
     else
     {
       exit(1);
     }
  }
  return 0;
}
//________________________________________________________________________________________
TH1D * NuFlux::SuperkFluxTuningWghts(int nu_pdg) const
{
  if      ( nu_pdg == conv::kPdgNuE     ) return fSuperkNuEFluxTunWght;
  else if ( nu_pdg == conv::kPdgNuEBar  ) return fSuperkNuEBarFluxTunWght;
  else if ( nu_pdg == conv::kPdgNuMu    ) return fSuperkNuMuFluxTunWght;
  else if ( nu_pdg == conv::kPdgNuMuBar ) return fSuperkNuMuBarFluxTunWght;
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
    TH1D * h = new TH1D(*fSuperkNuMuFluxTunWght);
    h->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK numu flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMuBar/TuneWghts")) { 
    TH1D * h = new TH1D(*fSuperkNuMuBarFluxTunWght);
    h->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK numubar flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuE/TuneWghts")) { 
    TH1D * h = new TH1D(*fSuperkNuEFluxTunWght);
    h->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK nue flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuEBar/TuneWghts")) { 
    TH1D * h = new TH1D(*fSuperkNuEBarFluxTunWght);
    h->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK nuebar flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMu")) { 
    TH1D * hnom = new TH1D(*fSuperkNuMuFluxNom);
    TH1D * htun = new TH1D(*fSuperkNuMuFluxTun);
    hnom->SetDirectory(0);
    htun->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK numu flux (#neutrinos /cm2 /x POT)");
    hnom->Draw();
    htun->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMuBar")) { 
    TH1D * hnom = new TH1D(*fSuperkNuMuBarFluxNom);
    TH1D * htun = new TH1D(*fSuperkNuMuBarFluxTun);
    hnom->SetDirectory(0);
    htun->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK numubar flux (#neutrinos /cm2 /x POT)");
    hnom->Draw();
    htun->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuE")) { 
    TH1D * hnom = new TH1D(*fSuperkNuEFluxNom);
    TH1D * htun = new TH1D(*fSuperkNuEFluxTun);
    hnom->SetDirectory(0);
    htun->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK nue flux (#neutrinos /cm2 /x POT)");
    hnom->Draw();
    htun->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuEBar")) { 
    TH1D * hnom = new TH1D(*fSuperkNuEBarFluxNom);
    TH1D * htun = new TH1D(*fSuperkNuEBarFluxTun);
    hnom->SetDirectory(0);
    htun->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK nuebar flux (#neutrinos /cm2 /x POT)");
    hnom->Draw();
    htun->Draw("same");
    c->Update();
  }
}
//________________________________________________________________________________________
