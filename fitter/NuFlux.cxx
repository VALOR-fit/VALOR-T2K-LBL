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
#define __show_flux_twk_mesg_info__   
#define __show_flux_twk_mesg_debug__   
#define __show_flux_twk_mesg_debug_level2__  #even more verbose

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <TVirtualPad.h>
#include <TLegend.h>

#include "NuFlux.h" 
#include "BANFF.h" 
#include "Conventions.h" 
#include "Utilities.h" 

using std::cout;
using std::cerr;
using std::setw;
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
  fSuperkNuMuFluxNom           = 0; 
  fSuperkNuMuBarFluxNom        = 0; 
  fSuperkNuEFluxNom            = 0; 
  fSuperkNuEBarFluxNom         = 0; 
  fSuperkNuMuFluxTun           = 0; 
  fSuperkNuMuBarFluxTun        = 0; 
  fSuperkNuEFluxTun            = 0; 
  fSuperkNuEBarFluxTun         = 0; 
  fSuperkNuMuFluxTwk           = 0; 
  fSuperkNuMuBarFluxTwk        = 0; 
  fSuperkNuEFluxTwk            = 0; 
  fSuperkNuEBarFluxTwk         = 0; 
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
  if ( fSuperkNuMuFluxTwk          ) delete fSuperkNuMuFluxTwk;
  if ( fSuperkNuMuBarFluxTwk       ) delete fSuperkNuMuBarFluxTwk;
  if ( fSuperkNuEFluxTwk           ) delete fSuperkNuEFluxTwk;
  if ( fSuperkNuEBarFluxTwk        ) delete fSuperkNuEBarFluxTwk;
  if ( fSuperkNuMuFluxTunWght      ) delete fSuperkNuMuFluxTunWght;
  if ( fSuperkNuMuBarFluxTunWght   ) delete fSuperkNuMuBarFluxTunWght;
  if ( fSuperkNuEFluxTunWght       ) delete fSuperkNuEFluxTunWght;
  if ( fSuperkNuEBarFluxTunWght    ) delete fSuperkNuEBarFluxTunWght;
}
//________________________________________________________________________________________
bool NuFlux::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read flux inputs to the oscillation fit
//
#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Reading flux inputs from input XML file." << endl;
#endif

  //
  // Load tuned fluxes
  //

  // SuperK numu flux
  fSuperkNuMuFluxTun = util::XML2TH1D(xml,root_node,"flux/superk/tuned/numu");
  if(!fSuperkNuMuFluxTun) {
    cerr << "Exp/Flux >> No input SuperK numu tuned flux histogram!" << endl;
    exit(1);
  }
  // SuperK numubar flux
  fSuperkNuMuBarFluxTun = util::XML2TH1D(xml,root_node,"flux/superk/tuned/numubar");
  if(!fSuperkNuMuBarFluxTun) {
    cerr << "Exp/Flux >> No input SuperK numubar tuned flux histogram!" << endl;
    exit(1);
  }
  // SuperK nue flux
  fSuperkNuEFluxTun = util::XML2TH1D(xml,root_node,"flux/superk/tuned/nue");
  if(!fSuperkNuEFluxTun) {
    cerr << "Exp/Flux >> No input SuperK nue tuned flux histogram!" << endl;
    exit(1);
  }
  // SuperK nuebar flux
  fSuperkNuEBarFluxTun = util::XML2TH1D(xml,root_node,"flux/superk/tuned/nuebar");
  if(!fSuperkNuEBarFluxTun) {
    cerr << "Exp/Flux >> No input SuperK nuebar tuned flux histogram!" << endl;
    exit(1);
  }

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Done loading tuned flux histograms" << endl;
#endif

  //
  // Load nominal -> tuned flux weights
  //

  // SuperK numu flux tuning weights
  fSuperkNuMuFluxTunWght = 
     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/numu");
  if(!fSuperkNuMuFluxTunWght) {
    cerr << "Exp/Flux >> No input SuperK numu flux tuning weights!" << endl;
    exit(1);
  }
  // SuperK numubar flux tuning weights
  fSuperkNuMuBarFluxTunWght = 
     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/numubar");
  if(!fSuperkNuMuBarFluxTunWght) {
    cerr << "Exp/Flux >> No input SuperK numubar flux tuning weights!" << endl;
    exit(1);
  }
  // SuperK nue flux tuning weights
  fSuperkNuEFluxTunWght = 
     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/nue");
  if(!fSuperkNuEFluxTunWght) {
    cerr << "Exp/Flux >> No input SuperK nue flux tuning weights!" << endl;
    exit(1);
  }
  // SuperK nuebar flux tuning weights
  fSuperkNuEBarFluxTunWght = 
     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/nuebar");
  if(!fSuperkNuEBarFluxTunWght) {
    cerr << "Exp/Flux >> No input SuperK nuebar flux tuning weights!" << endl;
    exit(1);
  }

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Done loading flux tuning weights" << endl;
#endif

  //
  // Book and initialize other auxiliary histograms (nominal flux, tweaked flux)
  //

  this->BookAndInitNomFlux ();
  this->BookAndInitTwkFlux ();

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Done reading flux inputs!" << endl;
#endif

  return true;
}
//________________________________________________________________________________________
void NuFlux::BookAndInitNomFlux(void)
{
// Create nominal flux histograms and initialize them using the tuned flux.
// Then de-weight with flux tuning weights to recreate the nominal flux.

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Booking and re-creating nominal flux histograms" << endl;
#endif

  const int n_nu = 4;
  int nu_pdg[n_nu] = { conv::kPdgNuE, conv::kPdgNuMu, conv::kPdgNuEBar, conv::kPdgNuMuBar };
  for(int inu = 0; inu < n_nu; inu++) {
     TH1D * hnom = this->SuperkFluxSpectrum(nu_pdg[inu], false, false);
     TH1D * htun = this->SuperkFluxSpectrum(nu_pdg[inu], true,  false);
     hnom = new TH1D(*htun);
     hnom -> SetDirectory(0);
     for(int ie = 1; ie <= htun->GetNbinsX(); ie++) {
        double E = htun->GetBinCenter (ie);
        double f = htun->GetBinContent(ie);
        double w = this->SuperkFluxTuneWght(nu_pdg[inu],E);
        if(w >=0) { f /= w; }
        hnom->SetBinContent(ie,f);
     }//ie
  }//inu
}
//________________________________________________________________________________________
void NuFlux::BookAndInitTwkFlux(void)
{
// Create tweaked flux histograms and initialize them using the tuned flux.

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Booking and initializing tuned+tweaked flux histograms" << endl;
#endif

  fSuperkNuMuFluxTwk    = new TH1D(*fSuperkNuMuFluxTun   );
  fSuperkNuMuBarFluxTwk = new TH1D(*fSuperkNuMuBarFluxTun);
  fSuperkNuEFluxTwk     = new TH1D(*fSuperkNuEFluxTun    );
  fSuperkNuEBarFluxTwk  = new TH1D(*fSuperkNuEBarFluxTun );

  fSuperkNuMuFluxTwk    -> SetDirectory(0);
  fSuperkNuMuBarFluxTwk -> SetDirectory(0);
  fSuperkNuEFluxTwk     -> SetDirectory(0);
  fSuperkNuEBarFluxTwk  -> SetDirectory(0);
}
//________________________________________________________________________________________
void NuFlux::CalcTwkFlux(void) const
{
#ifdef __show_flux_twk_mesg_info__
   cout << "Exp/Flux >> Calculating tweaked flux ..." << endl;
#endif
#ifdef __show_flux_twk_mesg_debug_level2__
   fBANFF->PrintSystematicVariations();
#endif

  const int np = conv::kNumOfParamsInBANFFfit;

  for(int ip = 0; ip < np; ip++) {
     conv::PhysParamId_t p = conv::PhysParamId_t ( (int)conv::kNuisBANFF00 + ip );
     // flux param?
     bool is_sk_flux_param = util::IsBANFFfitSuperkFluxParam(p);
     if(is_sk_flux_param) 
     {
        // Get current & best-fit values of current BANFF Superk flux param
        double banff_current = fBANFF->CurrentValue(p);
      //double banff_bestfit = fBANFF->BestFitValue(p);       

#ifdef __show_flux_twk_mesg_debug__
        cout << "Exp/Flux >> Taking into account the effect of tweaking BANFF flux param: "
             << util::AsString(p) << " (current value = " << banff_current << ")" << endl;
#endif

        // BANFF param -> neutrino PDG code
        int nu_pdg = util::BANFFfitSuperkFluxParam2NuPdg(p);

        // BANFF param -> true energy range [variable size bin]
        double Emin = -1;
        double Emax = -1;
        util::BANFFfitSuperkFluxParamEtrueRange(p,Emin,Emax);

        // Access tweaked & untweaked flux histograms
        TH1D * flux_spectrum_twk = this->SuperkFluxSpectrum(nu_pdg,true,true );
        TH1D * flux_spectrum     = this->SuperkFluxSpectrum(nu_pdg,true,false);

        // Note: 
        // The flux_spectrum histograms used in this code, use a binning scheme
        // (typically fixed-size 50-MeV bins, although other fixed or variable binning
        // schemes are possible) which is, in principle, different from the coarser
        // binning scheme used by BANFF.
        // Below, we will refer to the two binning schemes as `own' and `BANFF'

        // Find range of `own' bins which corresponds to the current `BANFF' bin ([Emin,Emax] range)
        int bin_min = flux_spectrum->FindBin(Emin + conv::kASmallNum); 
        int bin_max = flux_spectrum->FindBin(Emax - conv::kASmallNum);

#ifdef __show_flux_twk_mesg_debug__
        cout << "Exp/Flux >> BANFF parameter tweaks Superk flux for neutrino PDG = " << nu_pdg 
             << " in true energy range = [" << Emin << ", " << Emax << "]"
             << endl;
        cout << "Exp/Flux >> In our binning scheme, this energy range corresponds to bins = ["
             << bin_min << ", " << bin_max << "]" 
             << endl;
#endif
        // Loop over `own' bins contained in the current `BANFF' bin and
        // calculate the "tuned" and "tuned+tweaked" flux in that `BANFF' bin.
        double flux_var_twk = 0;
        double flux_var     = 0;
        for ( int bin = bin_min; bin <= bin_max; bin++ ) {
           // Un-tweaked flux in current bin 
           double flux = flux_spectrum->GetBinContent(bin);
           // Add to total flux in larger variable-size bin
           flux_var += flux;
        }
        flux_var_twk = flux_var * banff_current;
#ifdef __show_flux_twk_mesg_debug_level2__
        cout << "Exp/Flux >> Flux in current BANFF bin: " 
             << flux_var << " --> " << flux_var_twk << endl;
#endif

        // Loop over `own' bins contained in the current `BANFF' bin and
        // "distribute" the `BANFF' bin flux variation to our `own' bins.
        for ( int bin = bin_min; bin <= bin_max; bin++ ) {
           // Un-tweaked flux in current bin 
           double flux = flux_spectrum->GetBinContent(bin);
           // Fraction of flux of the larger variable-size bin which is contained 
           // in the current bin
           double frac_flux = flux / flux_var;
           // Tweaked flux in current bin
           double flux_twk = flux + frac_flux * (flux_var_twk - flux_var);
           flux_spectrum_twk->SetBinContent(bin,flux_twk);
#ifdef __show_flux_twk_mesg_debug_level2__
        cout << "Exp/Flux >> - Our bin = " << bin 
             << ", contains " << frac_flux*100 << "% of flux of current BANFF bin"
             << ", flux = " << flux << " --> " << flux_twk << endl;
#endif
        } // bin
     } 
     else 
     {
#ifdef __show_flux_twk_mesg_debug__
        cout << "Exp/Flux >> BANFF param: " << util::AsString(p) 
             << " is not a flux param and it is skipped..." << endl;
#endif
     }//sk flux param?
  }//ip loop  
}
//________________________________________________________________________________________
void NuFlux::TellMeAboutTheBANFFfit(BANFF * banff)
{
  fBANFF = banff;
}
//________________________________________________________________________________________
double NuFlux::SuperkFlux(int nu_pdg, double Etrue, bool tuned, bool tweaked) const
{
  // update tweaked flux (if not already) as a result of the latest BANFF param twk
  if( ! fBANFF->IsFluxWghtUpdated() ) {
    this->CalcTwkFlux();
    fBANFF->MarkFluxWghtUpdated();
  }

  TH1D * flux_spectrum = this->SuperkFluxSpectrum(nu_pdg,tuned,tweaked);
  if(!flux_spectrum) {
    cerr << "Exp/Flux >> No " 
         << ((tuned) ? "`tuned'" : "`nominal'")
         << " + "
         << ((tweaked) ? "`tweaked'" : "`untweaked'")
         << " flux spectrum for neutrino PDG = " << nu_pdg
         << endl;
  }
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
  double flux_twk = this->SuperkFlux(nu_pdg, Etrue, true, true );
  double flux     = this->SuperkFlux(nu_pdg, Etrue, true, false);

  double wght = 1.0;
  if(flux >= 0) {
    wght = flux_twk / flux;
  }
  return wght;
}
//________________________________________________________________________________________
TH1D * NuFlux::SuperkFluxSpectrum(int nu_pdg, bool tuned, bool tweaked) const
{
  if(tuned) 
  {
     if(tweaked) {
        if      ( nu_pdg == conv::kPdgNuE     ) return fSuperkNuEFluxTwk;
        else if ( nu_pdg == conv::kPdgNuEBar  ) return fSuperkNuEBarFluxTwk;
        else if ( nu_pdg == conv::kPdgNuMu    ) return fSuperkNuMuFluxTwk;
        else if ( nu_pdg == conv::kPdgNuMuBar ) return fSuperkNuMuBarFluxTwk;
        else
        {
           cerr << "Exp/Flux >> No tweaked flux for particle with PDG = " 
                << nu_pdg << endl;
           exit(1);
        }
     }//tweaked?
     else {
        if      ( nu_pdg == conv::kPdgNuE     ) return fSuperkNuEFluxTun;
        else if ( nu_pdg == conv::kPdgNuEBar  ) return fSuperkNuEBarFluxTun;
        else if ( nu_pdg == conv::kPdgNuMu    ) return fSuperkNuMuFluxTun;
        else if ( nu_pdg == conv::kPdgNuMuBar ) return fSuperkNuMuBarFluxTun;
        else
        {
           cerr << "Exp/Flux >> No tuned flux for particle with PDG = " 
                << nu_pdg << endl;
           exit(1);
        }
     }
  }//tuned?

  else
  {
     if(tweaked) {
       cerr << "Exp/Flux >> BANFF tweaks are applied on the tuned flux, not the nominal!" 
            << endl;
       exit(1);
     }

     if      ( nu_pdg == conv::kPdgNuE     ) return fSuperkNuEFluxNom;
     else if ( nu_pdg == conv::kPdgNuEBar  ) return fSuperkNuEBarFluxNom;
     else if ( nu_pdg == conv::kPdgNuMu    ) return fSuperkNuMuFluxNom;
     else if ( nu_pdg == conv::kPdgNuMuBar ) return fSuperkNuMuBarFluxNom;
     else
     {
       cerr << "Exp/Flux >> No nominal flux for particle with PDG = " 
            << nu_pdg << endl;
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
