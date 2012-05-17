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
//#define __show_flux_twk_mesg_debug__   
//#define __show_flux_twk_mesg_debug_level2__  #even more verbose

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <TVirtualPad.h>
#include <TLegend.h>

#include "NuFlux.h" 
#include "BANFF.h" 
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
  fSKNuMuFluxNominal          = 0; 
  fSKNuMuBarFluxNominal       = 0; 
  fSKNuEFluxNominal           = 0; 
  fSKNuEBarFluxNominal        = 0; 
  fSKNuMuFluxTunedExt         = 0; 
  fSKNuMuBarFluxTunedExt      = 0; 
  fSKNuEFluxTunedExt          = 0; 
  fSKNuEBarFluxTunedExt       = 0; 
  fSKNuMuFluxTunedFinal       = 0; 
  fSKNuMuBarFluxTunedFinal    = 0; 
  fSKNuEFluxTunedFinal        = 0; 
  fSKNuEBarFluxTunedFinal     = 0; 
  fSKNuMuFluxCurrent          = 0; 
  fSKNuMuBarFluxCurrent       = 0; 
  fSKNuEFluxCurrent           = 0; 
  fSKNuEBarFluxCurrent        = 0; 
  fSKNuMuFluxExtTuningWght    = 0;
  fSKNuMuBarFluxExtTuningWght = 0;
  fSKNuEFluxExtTuningWght     = 0;
  fSKNuEBarFluxExtTuningWght  = 0;
  fBANFF                      = 0;
}
//________________________________________________________________________________________
void NuFlux::CleanUp(void)
{
  if ( fSKNuMuFluxNominal          ) delete fSKNuMuFluxNominal;
  if ( fSKNuMuBarFluxNominal       ) delete fSKNuMuBarFluxNominal;
  if ( fSKNuEFluxNominal           ) delete fSKNuEFluxNominal;
  if ( fSKNuEBarFluxNominal        ) delete fSKNuEBarFluxNominal;
  if ( fSKNuMuFluxTunedExt         ) delete fSKNuMuFluxTunedExt;
  if ( fSKNuMuBarFluxTunedExt      ) delete fSKNuMuBarFluxTunedExt;
  if ( fSKNuEFluxTunedExt          ) delete fSKNuEFluxTunedExt;
  if ( fSKNuEBarFluxTunedExt       ) delete fSKNuEBarFluxTunedExt;
  if ( fSKNuMuFluxTunedFinal       ) delete fSKNuMuFluxTunedFinal;
  if ( fSKNuMuBarFluxTunedFinal    ) delete fSKNuMuBarFluxTunedFinal;
  if ( fSKNuEFluxTunedFinal        ) delete fSKNuEFluxTunedFinal;
  if ( fSKNuEBarFluxTunedFinal     ) delete fSKNuEBarFluxTunedFinal;
  if ( fSKNuMuFluxCurrent          ) delete fSKNuMuFluxCurrent;
  if ( fSKNuMuBarFluxCurrent       ) delete fSKNuMuBarFluxCurrent;
  if ( fSKNuEFluxCurrent           ) delete fSKNuEFluxCurrent;
  if ( fSKNuEBarFluxCurrent        ) delete fSKNuEBarFluxCurrent;
  if ( fSKNuMuFluxExtTuningWght    ) delete fSKNuMuFluxExtTuningWght;
  if ( fSKNuMuBarFluxExtTuningWght ) delete fSKNuMuBarFluxExtTuningWght;
  if ( fSKNuEFluxExtTuningWght     ) delete fSKNuEFluxExtTuningWght;
  if ( fSKNuEBarFluxExtTuningWght  ) delete fSKNuEBarFluxExtTuningWght;
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
  // Load fluxes tuned based on external data only
  //

  // SuperK numu flux
  fSKNuMuFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned/numu");
  if(!fSKNuMuFluxTunedExt) {
    cerr << "Exp/Flux >> No input SuperK numu ext/tuned flux histogram!" << endl;
    exit(1);
  }
  // SuperK numubar flux
  fSKNuMuBarFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned/numubar");
  if(!fSKNuMuBarFluxTunedExt) {
    cerr << "Exp/Flux >> No input SuperK numubar ext/tuned flux histogram!" << endl;
    exit(1);
  }
  // SuperK nue flux
  fSKNuEFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned/nue");
  if(!fSKNuEFluxTunedExt) {
    cerr << "Exp/Flux >> No input SuperK nue ext/tuned flux histogram!" << endl;
    exit(1);
  }
  // SuperK nuebar flux
  fSKNuEBarFluxTunedExt = util::XML2TH1D(xml,root_node,"flux/superk/tuned/nuebar");
  if(!fSKNuEBarFluxTunedExt) {
    cerr << "Exp/Flux >> No input SuperK nuebar ext/tuned flux histogram!" << endl;
    exit(1);
  }

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Done loading ext/tuned flux histograms" << endl;
#endif

  //
  // Load nominal -> ext.-tuned flux weights
  //

  // SuperK numu flux tuning weights
  fSKNuMuFluxExtTuningWght = 
     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/numu");
  if(!fSKNuMuFluxExtTuningWght) {
    cerr << "Exp/Flux >> No input SuperK numu flux ext/tuning weights!" << endl;
    exit(1);
  }
  // SuperK numubar flux tuning weights
  fSKNuMuBarFluxExtTuningWght = 
     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/numubar");
  if(!fSKNuMuBarFluxExtTuningWght) {
    cerr << "Exp/Flux >> No input SuperK numubar flux ext/tuning weights!" << endl;
    exit(1);
  }
  // SuperK nue flux tuning weights
  fSKNuEFluxExtTuningWght = 
     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/nue");
  if(!fSKNuEFluxExtTuningWght) {
    cerr << "Exp/Flux >> No input SuperK nue flux ext/tuning weights!" << endl;
    exit(1);
  }
  // SuperK nuebar flux tuning weights
  fSKNuEBarFluxExtTuningWght = 
     util::XML2TH1D(xml,root_node,"flux/superk/tuning_weights/nuebar");
  if(!fSKNuEBarFluxExtTuningWght) {
    cerr << "Exp/Flux >> No input SuperK nuebar flux ext/tuning weights!" << endl;
    exit(1);
  }

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Done loading flux tuning weights" << endl;
#endif

  //
  // Book and initialize other auxiliary flux histograms not loaded from files
  //
  this->BookAndInitAuxSpectra();

#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Done reading flux inputs!" << endl;
#endif

  return true;
}
//________________________________________________________________________________________
void NuFlux::BookAndInitAuxSpectra(void)
{
#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Booking and initializing/computing auxiliary flux spectra" << endl;
#endif

  //
  // Create nominal flux histograms and initialize them using the ext/tuned flux.
  // Then de-weight them with input ext/tuning weights to recreate the nominal flux.
  //
#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Booking and re-creating nominal flux histograms" << endl;
#endif

  fSKNuMuFluxNominal    = new TH1D(*fSKNuMuFluxTunedExt   );
  fSKNuMuBarFluxNominal = new TH1D(*fSKNuMuBarFluxTunedExt);
  fSKNuEFluxNominal     = new TH1D(*fSKNuEFluxTunedExt    );
  fSKNuEBarFluxNominal  = new TH1D(*fSKNuEBarFluxTunedExt );

  fSKNuMuFluxNominal    -> SetDirectory(0);
  fSKNuMuBarFluxNominal -> SetDirectory(0);
  fSKNuEFluxNominal     -> SetDirectory(0);
  fSKNuEBarFluxNominal  -> SetDirectory(0);

  const int n_nu = 4;
  int nu_pdg[n_nu] = { conv::kPdgNuE, conv::kPdgNuMu, conv::kPdgNuEBar, conv::kPdgNuMuBar };
  for(int inu = 0; inu < n_nu; inu++) {
     TH1D * hnom = this->SuperkFluxSpectrum(nu_pdg[inu], conv::kFlxNominal  );
     for(int ie = 1; ie <= hnom->GetNbinsX(); ie++) {
        double E = hnom->GetBinCenter (ie);
        double f = hnom->GetBinContent(ie);
        double w = this->SuperkFluxExtTuningWght(nu_pdg[inu],E);
        if(w >=0) { f /= w; }
        hnom->SetBinContent(ie,f);
     }//ie
  }//inu

  //
  // Create the fin/tuned flux histograms and initialize them using the ext/tuned flux
  // Then re-weight them using the weights that correspond to the BANFF best-fit params.
  //
  fSKNuMuFluxTunedFinal    = new TH1D(*fSKNuMuFluxTunedExt   );
  fSKNuMuBarFluxTunedFinal = new TH1D(*fSKNuMuBarFluxTunedExt);
  fSKNuEFluxTunedFinal     = new TH1D(*fSKNuEFluxTunedExt    );
  fSKNuEBarFluxTunedFinal  = new TH1D(*fSKNuEBarFluxTunedExt );

  fSKNuMuFluxTunedFinal    -> SetDirectory(0);
  fSKNuMuBarFluxTunedFinal -> SetDirectory(0);
  fSKNuEFluxTunedFinal     -> SetDirectory(0);
  fSKNuEBarFluxTunedFinal  -> SetDirectory(0);

  this->ApplyBANFFWeights(false); // reweight ext/tuned to fin/tuned flux (corresponds to BANFF best-fit)

  //
  // Create current flux histograms and initialize them using the fin/tuned flux.
  //
#ifdef __show_init_mesg_info__
  cout << "Exp/Flux >> Booking and initializing current flux histograms" << endl;
#endif

  fSKNuMuFluxCurrent    = new TH1D(*fSKNuMuFluxTunedFinal   );
  fSKNuMuBarFluxCurrent = new TH1D(*fSKNuMuBarFluxTunedFinal);
  fSKNuEFluxCurrent     = new TH1D(*fSKNuEFluxTunedFinal    );
  fSKNuEBarFluxCurrent  = new TH1D(*fSKNuEBarFluxTunedFinal );

  fSKNuMuFluxCurrent    -> SetDirectory(0);
  fSKNuMuBarFluxCurrent -> SetDirectory(0);
  fSKNuEFluxCurrent     -> SetDirectory(0);
  fSKNuEBarFluxCurrent  -> SetDirectory(0);

}
//________________________________________________________________________________________
void NuFlux::ApplyBANFFWeights(bool get_current) const
{
// Apply BANFF weights to ext/tuned flux.
//
// Input param     | Action    
// ---------------------------------------------------
// false           | apply BANFF best-fit params to fin/tuned flux
// true            | apply BANFF currrent params to current   flux 

#ifdef __show_flux_twk_mesg_info__
   cout << "Exp/Flux >> Applying BANFF weights to SuperK flux " 
        << ( (get_current) ? "[ext/tuned --> current]" : "[ext/tuned --> fin/tuned]" )
        << endl;   
#endif
#ifdef __show_flux_twk_mesg_debug_level2__
   fBANFF->PrintSystematicVariations();
#endif

  const int np = conv::kNumOfParamsInBANFFfit;

  for(int ip = 0; ip < np; ip++) {
     conv::PhysParamId_t p = conv::PhysParamId_t ( (int)conv::kNuisBANFF00 + ip );
     // flux param?
     bool is_sk_flux_param = conv::kParamInBANFFfitIsSuperkFluxParam[ip];
     if(is_sk_flux_param) 
     {
        // Get the current or best-fit BANFF Superk flux param value
        double banff_param_value = (get_current) ? 
             fBANFF->CurrentValue(p) :
             fBANFF->BestFitValue(p) ; 

#ifdef __show_flux_twk_mesg_debug__
        cout << "Exp/Flux >> Taking into account the effect of BANFF flux param: "
             << util::AsString(p) << " (value = " << banff_param_value << ")" << endl;
#endif

        // BANFF param -> neutrino PDG code
        int nu_pdg = conv::kNuPdgCodeOfFluxParamInBANFFfit[ip];

        // BANFF param -> true energy range [variable size bin]
        double Emin = conv::kEtrueRangeOfFluxParamInBANFFfit[ip][0];
        double Emax = conv::kEtrueRangeOfFluxParamInBANFFfit[ip][1];
        assert(Emin >= 0. && Emax > 0);

        // Access flux histograms

        // original spectrum:
        TH1D * flux_spectrum_orig = this->SuperkFluxSpectrum(nu_pdg, conv::kFlxExtTuning); 
        // original spectrum + BANFF weights
        TH1D * flux_spectrum_wght = (get_current) ?
           this->SuperkFluxSpectrum(nu_pdg, conv::kFlxCurrent) :
           this->SuperkFluxSpectrum(nu_pdg, conv::kFlxFinalTuning);

        // Note: 
        // The flux_spectrum histograms used in this code, use a binning scheme
        // (typically fixed-size 50-MeV bins, although other fixed or variable binning
        // schemes are possible) which is, in principle, different from the coarser
        // binning scheme used by BANFF.
        // Below, we will refer to the two binning schemes as `own' and `BANFF'

        // Find range of `own' bins which corresponds to the current `BANFF' bin ([Emin,Emax] range)
        int bin_min = flux_spectrum_orig->FindBin(Emin + conv::kASmallNum); 
        int bin_max = flux_spectrum_orig->FindBin(Emax - conv::kASmallNum);

#ifdef __show_flux_twk_mesg_debug__
        cout << "Exp/Flux >> BANFF parameter tweaks Superk flux for neutrino PDG = " << nu_pdg 
             << " in true energy range = [" << Emin << ", " << Emax << "]"
             << endl;
        cout << "Exp/Flux >> In our binning scheme, this energy range corresponds to bins = ["
             << bin_min << ", " << bin_max << "]" 
             << endl;
#endif

        // Loop over `own' bins contained in the current `BANFF' bin and
        // "distribute" the `BANFF' bin flux variation to our `own' bins.
        for ( int bin = bin_min; bin <= bin_max; bin++ ) {
           // Un-tweaked flux in current bin 
           double flux = flux_spectrum_orig->GetBinContent(bin);
           // Tweaked flux in current bin
           double flux_twk = flux * banff_param_value;
#ifdef __show_flux_twk_mesg_debug_level2__
           cout << "Exp/Flux >> - Our bin = " << bin
                << ", flux = " << flux << " --> " << flux_twk << endl;
#endif
           if(flux_twk < 0) {
                cerr << "Exp/Flux >> Tweaked flux for neutrino PDG = " << nu_pdg 
                     << " in bin = " << bin << " is = " << flux_twk << " !!! " << endl;
                cerr << "Exp/Flux >> Negative value was the result of tweaking BANFF param: " 
                     << util::AsString(p) << " (value = " << banff_param_value << ")" << endl;
                exit(1);
           }
           flux_spectrum_wght->SetBinContent(bin,flux_twk);
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
double NuFlux::SuperkFlux(int nu_pdg, double Etrue, conv::FluxType_t flxtype) const
{
  // update current flux (if not already) as a result of the latest BANFF param twk
  if( flxtype == conv::kFlxCurrent && !fBANFF->IsFluxWghtUpdated() ) {
    this->ApplyBANFFWeights(true);
    fBANFF->MarkFluxWghtUpdated();
  }

  TH1D * flux_spectrum = this->SuperkFluxSpectrum(nu_pdg,flxtype);
  if(!flux_spectrum) {
    cerr << "Exp/Flux >> No " << util::AsString(flxtype)
         << " flux spectrum for neutrino PDG = " << nu_pdg
         << endl;
  }
  assert(flux_spectrum);
  int bin = flux_spectrum->FindBin(Etrue);
  double flux = flux_spectrum->GetBinContent(bin);
  return flux;
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxExtTuningWght(int nu_pdg, double Etrue) const
{
// Weights required to re-weight nominal flux to ext/tuned flux
//
//  if(Etrue > conv::kEtrueMaxForFluxTuneWeights) 
//  {
//    return 1;
//  }

  TH1D * wght_spectrum = this->SuperkFluxExtTuningWghts(nu_pdg);
  assert(wght_spectrum);

  int bin = wght_spectrum->FindBin(Etrue);
  double wght = wght_spectrum->GetBinContent(bin);
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxFinTuningWght(int nu_pdg, double Etrue) const
{
// Weights required to re-weight ext/tuned flux to fin/tuned flux
//
  double flux_ext = this->SuperkFlux(nu_pdg, Etrue, conv::kFlxExtTuning  );
  double flux_fin = this->SuperkFlux(nu_pdg, Etrue, conv::kFlxFinalTuning);

  double wght = 1.0;
  if(flux_ext >= 0) {
    wght = flux_fin / flux_ext;
  }
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxSystWght(int nu_pdg, double Etrue) const
{
// Weights required to re-weight fin/tuned flux to current (fin/tuned + systematic
// variarions) flux
//
  double flux_twk = this->SuperkFlux(nu_pdg, Etrue, conv::kFlxCurrent    );
  double flux     = this->SuperkFlux(nu_pdg, Etrue, conv::kFlxFinalTuning);

  double wght = 1.0;
  if(flux >= 0) {
    wght = flux_twk / flux;
  }
  return wght;
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxExtTuningWght(
  conv::SuperkSample_t /*s*/, conv::SuperkMode_t m, double Etrue) const
{
  int flux_nu_pdg = util::SuperkMode2FluxNuPdg(m);
  return this->SuperkFluxExtTuningWght(flux_nu_pdg,Etrue);
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxFinTuningWght(
  conv::SuperkSample_t /*s*/, conv::SuperkMode_t m, double Etrue) const
{
  int flux_nu_pdg = util::SuperkMode2FluxNuPdg(m);
  return this->SuperkFluxFinTuningWght(flux_nu_pdg,Etrue);
}
//________________________________________________________________________________________
double NuFlux::SuperkFluxSystWght(
  conv::SuperkSample_t /*s*/, conv::SuperkMode_t m, double Etrue) const
{
  int flux_nu_pdg = util::SuperkMode2FluxNuPdg(m);
  return this->SuperkFluxSystWght(flux_nu_pdg,Etrue);
}
//________________________________________________________________________________________
TH1D * NuFlux::SuperkFluxSpectrum(int nu_pdg, conv::FluxType_t flxtype) const
{
  switch(flxtype) {

     // 
     case (conv::kFlxNominal) :
     {
        if      ( nu_pdg == conv::kPdgNuE     ) return fSKNuEFluxNominal;
        else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKNuEBarFluxNominal;
        else if ( nu_pdg == conv::kPdgNuMu    ) return fSKNuMuFluxNominal;
        else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKNuMuBarFluxNominal;
        else
        {
           cerr << "Exp/Flux >> No nominal flux spectrum for particle with PDG = " 
             << nu_pdg << endl;
           exit(1);
        }
        break;
     }

     //
     case (conv::kFlxExtTuning) : 
     {
        if      ( nu_pdg == conv::kPdgNuE     ) return fSKNuEFluxTunedExt;
        else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKNuEBarFluxTunedExt;
        else if ( nu_pdg == conv::kPdgNuMu    ) return fSKNuMuFluxTunedExt;
        else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKNuMuBarFluxTunedExt;
        else
        {
           cerr << "Exp/Flux >> No ext/tuned flux spectrum for particle with PDG = " 
             << nu_pdg << endl;
           exit(1);
        }
        break;
     }

     //
     case (conv::kFlxFinalTuning) :
     {
        if      ( nu_pdg == conv::kPdgNuE     ) return fSKNuEFluxTunedFinal;
        else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKNuEBarFluxTunedFinal;
        else if ( nu_pdg == conv::kPdgNuMu    ) return fSKNuMuFluxTunedFinal;
        else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKNuMuBarFluxTunedFinal;
        else
        {
           cerr << "Exp/Flux >> No fin/tuned flux spectrum for particle with PDG = " 
             << nu_pdg << endl;
           exit(1);
        }
        break;
     }

     //
     case (conv::kFlxCurrent) :
     {
        if      ( nu_pdg == conv::kPdgNuE     ) return fSKNuEFluxCurrent;
        else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKNuEBarFluxCurrent;
        else if ( nu_pdg == conv::kPdgNuMu    ) return fSKNuMuFluxCurrent;
        else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKNuMuBarFluxCurrent;
        else
        {
           cerr << "Exp/Flux >> No current flux spectrum for particle with PDG = " 
             << nu_pdg << endl;
           exit(1);
        }
        break;
     }

     //
     default:
        return 0;
  }

  return 0;
}
//________________________________________________________________________________________
TH1D * NuFlux::SuperkFluxExtTuningWghts(int nu_pdg) const
{
  if      ( nu_pdg == conv::kPdgNuE     ) return fSKNuEFluxExtTuningWght;
  else if ( nu_pdg == conv::kPdgNuEBar  ) return fSKNuEBarFluxExtTuningWght;
  else if ( nu_pdg == conv::kPdgNuMu    ) return fSKNuMuFluxExtTuningWght;
  else if ( nu_pdg == conv::kPdgNuMuBar ) return fSKNuMuBarFluxExtTuningWght;
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
    TH1D * h = new TH1D(*fSKNuMuFluxExtTuningWght);
    h->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK numu flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMuBar/TuneWghts")) { 
    TH1D * h = new TH1D(*fSKNuMuBarFluxExtTuningWght);
    h->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK numubar flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuE/TuneWghts")) { 
    TH1D * h = new TH1D(*fSKNuEFluxExtTuningWght);
    h->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK nue flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuEBar/TuneWghts")) { 
    TH1D * h = new TH1D(*fSKNuEBarFluxExtTuningWght);
    h->SetDirectory(0);
    h->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    h->GetYaxis()->SetTitle("SuperK nuebar flux weight");
    h->Draw();
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMu")) { 
    TH1D * hnom  = new TH1D(*fSKNuMuFluxNominal);
    TH1D * htune = new TH1D(*fSKNuMuFluxTunedExt);
    TH1D * htunf = new TH1D(*fSKNuMuFluxTunedFinal);
    hnom ->SetDirectory(0);
    htune->SetDirectory(0);
    htunf->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK numu flux (#neutrinos /cm2 /x POT)");
    hnom->Draw();
    htune->Draw("same");
    htunf->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuMuBar")) { 
    TH1D * hnom  = new TH1D(*fSKNuMuBarFluxNominal);
    TH1D * htune = new TH1D(*fSKNuMuBarFluxTunedExt);
    TH1D * htunf = new TH1D(*fSKNuMuBarFluxTunedFinal);
    hnom ->SetDirectory(0);
    htune->SetDirectory(0);
    htunf->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK numubar flux (#neutrinos /cm2 /x POT)");
    hnom->Draw();
    htune->Draw("same");
    htunf->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuE")) { 
    TH1D * hnom  = new TH1D(*fSKNuEFluxNominal);
    TH1D * htune = new TH1D(*fSKNuEFluxTunedExt);
    TH1D * htunf = new TH1D(*fSKNuEFluxTunedFinal);
    hnom ->SetDirectory(0);
    htune->SetDirectory(0);
    htunf->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK nue flux (#neutrinos /cm2 /x POT)");
    hnom->Draw();
    htune->Draw("same");
    htunf->Draw("same");
    c->Update();
  }
  else
  if(opt.Contains("Flux/SK/NuEBar")) { 
    TH1D * hnom  = new TH1D(*fSKNuEBarFluxNominal);
    TH1D * htune = new TH1D(*fSKNuEBarFluxTunedExt);
    TH1D * htunf = new TH1D(*fSKNuEBarFluxTunedFinal);
    hnom->SetDirectory(0);
    htune->SetDirectory(0);
    htunf->SetDirectory(0);
    hnom->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hnom->GetYaxis()->SetTitle("SuperK nuebar flux (#neutrinos /cm2 /x POT)");
    hnom->Draw();
    htune->Draw("same");
    htunf->Draw("same");
    c->Update();
  }
}
//________________________________________________________________________________________
