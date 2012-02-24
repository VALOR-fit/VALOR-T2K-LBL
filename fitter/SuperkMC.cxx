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
//#define __show_hypothesis__ 
//#define __show_hypothesis_probabilities__ 
//#define __show_syst_var_wght_mesg_info__ 
//#define __show_syst_var_wght_mesg_debug__ 
#define __show_rc_cholesky_mesg_info__
//#define __show_rc_cholesky_mesg_debug__ 
//#define __show_spectra_update_mesg_debug__ 
//#define __show_skipped_templates__ 

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TTree.h>
#include <TH2D.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TMath.h>

#include "SuperkMC.h"
#include "NuXSec.h"
#include "NuFlux.h"
#include "Utilities.h"
#include "Opt.h"

using std::cout;
using std::cerr;
using std::endl;

using namespace t2kana;
using namespace t2kana::conv;
using namespace t2kana::util;

ClassImp(SuperkMC)

//________________________________________________________________________________________
SuperkMC::SuperkMC()
{
  this->Init();
}
//________________________________________________________________________________________
SuperkMC::~SuperkMC()
{
  this->CleanUp();
}
//________________________________________________________________________________________
void SuperkMC::Init(void)
{
  fPOT = 0;

  for(int is=0; is<kns; is++) {
    fMCErecoSpectrumTot[is] = 0;
    for(int im=0; im<knm; im++) {
      fMCTemplNom      [is][im] = 0;
      fMCTemplTun      [is][im] = 0;
      fMCTemplSys      [is][im] = 0;
      fMCTemplFin      [is][im] = 0;
      fMCErecoSpectrum [is][im] = 0;
    }
  }
  
  fEnergyScale1RmuRelErr      = 0.;
  fEff1RmuCCQEOthRelErr       = 0.;
  fEff1RmuCCQERCCovarRelErr   = 0;
  fEff1RmuCCQERCCovarE1RelErr = 0.;
  fEff1RmuCCQERCCovarE2RelErr = 0.;      
  fEff1RmuCCnonQERelErr       = 0.;
  fEff1RmuNCRelErr            = 0.;
  fEff1RmuNuECCRelErr         = 0.;
  
  fEnergyScaleTwkDial         = 0.;
  fEff1RmuCCQEOthTwkDial      = 0.;
  fEff1RmuCCQERC0TwkDial      = 0.;
  fEff1RmuCCQERC1TwkDial      = 0.;
  fEff1RmuCCQERC2TwkDial      = 0.;
  fEff1RmuCCnonQETwkDial      = 0.;
  fEff1RmuNCTwkDial           = 0.;
  fEff1RmuNuECCTwkDial        = 0.;

  fEff1RmuCCQERCCovarChoLRelErr = 0;
  
  fMCTreeNuMu1Rmu      = 0;
  fMCTreeNuMuBar1Rmu   = 0;
  fMCTreeNuE1Rmu       = 0;
  fMCTreeNuESig1Rmu    = 0;
  fMCTreeNuMu1Re       = 0;
  fMCTreeNuMuBar1Re    = 0;
  fMCTreeNuE1Re        = 0;
  fMCTreeNuESig1Re     = 0;

  fMCPOTNuMu           = 0;
  fMCPOTNuMuBar        = 0;
  fMCPOTNuE            = 0;
  fMCPOTNuESig         = 0;

  fBrEreco        = 0;
  fBrEtrue        = 0;
  fBrPdg          = 0;
  fBrGenCode      = 0;
         
  fNuFlux         = 0;
  fNuXSec         = 0;  
}
//________________________________________________________________________________________
void SuperkMC::CleanUp(void)
{
  if(fMCTreeNuMu1Rmu    ) delete fMCTreeNuMu1Rmu;
  if(fMCTreeNuMuBar1Rmu ) delete fMCTreeNuMuBar1Rmu;
  if(fMCTreeNuE1Rmu     ) delete fMCTreeNuE1Rmu;
  if(fMCTreeNuESig1Rmu  ) delete fMCTreeNuESig1Rmu;
  if(fMCTreeNuMu1Re     ) delete fMCTreeNuMu1Re;
  if(fMCTreeNuMuBar1Re  ) delete fMCTreeNuMuBar1Re;
  if(fMCTreeNuE1Re      ) delete fMCTreeNuE1Re;
  if(fMCTreeNuESig1Re   ) delete fMCTreeNuESig1Re;

  for(int is = 0; is < kns; is++) {
     for(int im = 0; im < knm; im++) {
       if( fMCTemplNom     [is][im] ) delete fMCTemplNom     [is][im];
       if( fMCTemplTun     [is][im] ) delete fMCTemplTun     [is][im];
       if( fMCTemplSys     [is][im] ) delete fMCTemplSys     [is][im];
       if( fMCTemplFin     [is][im] ) delete fMCTemplFin     [is][im];
       if( fMCErecoSpectrum[is][im] ) delete fMCErecoSpectrum[is][im];
     }
  }
}
//________________________________________________________________________________________
bool SuperkMC::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
  this->CleanUp();

  // instantiate all necessary histograms
  this->Book();

#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> Reading SuperK inputs from input XML file." << endl;
#endif

  //
  // Get data POT
  // 
  fPOT = util::XML2Dbl(xml,root_node,"data_pot/superk");
#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> " << "POT_{data} = " << fPOT << endl;
#endif
  
  //
  // Get MC trees and MC POT
  //
  fMCPOTNuMu    = util::XML2Dbl  (xml,root_node,"mc/superk/nominal/numu/pot");
  fMCPOTNuMuBar = util::XML2Dbl  (xml,root_node,"mc/superk/nominal/numubar/pot");
  fMCPOTNuE     = util::XML2Dbl  (xml,root_node,"mc/superk/nominal/nuebeam/pot");
  fMCPOTNuESig  = util::XML2Dbl  (xml,root_node,"mc/superk/nominal/nuesignal/pot");

  fMCTreeNuMu1Rmu    = util::XML2TTree (xml,root_node,"mc/superk/nominal/numu/one_mu_ring");
  fMCTreeNuMuBar1Rmu = util::XML2TTree (xml,root_node,"mc/superk/nominal/numubar/one_mu_ring");
  fMCTreeNuE1Rmu     = util::XML2TTree (xml,root_node,"mc/superk/nominal/nuebeam/one_mu_ring");
  fMCTreeNuESig1Rmu  = util::XML2TTree (xml,root_node,"mc/superk/nominal/nuesignal/one_mu_ring");
  fMCTreeNuMu1Re     = util::XML2TTree (xml,root_node,"mc/superk/nominal/numu/one_e_ring");
  fMCTreeNuMuBar1Re  = util::XML2TTree (xml,root_node,"mc/superk/nominal/numubar/one_e_ring");
  fMCTreeNuE1Re      = util::XML2TTree (xml,root_node,"mc/superk/nominal/nuebeam/one_e_ring");
  fMCTreeNuESig1Re   = util::XML2TTree (xml,root_node,"mc/superk/nominal/nuesignal/one_e_ring");

  assert (fMCTreeNuMu1Rmu   );
  assert (fMCTreeNuMuBar1Rmu);
  assert (fMCTreeNuE1Rmu    );
  assert (fMCTreeNuESig1Rmu );
  assert (fMCTreeNuMu1Re    );
  assert (fMCTreeNuMuBar1Re );
  assert (fMCTreeNuE1Re     );
  assert (fMCTreeNuESig1Re  );

#ifdef __show_init_mesg_info__

  cout << "Exp/SKMC >> " << "POT {MC; numu}    = " << fMCPOTNuMu    << endl;
  cout << "Exp/SKMC >> " << "POT {MC; numubar} = " << fMCPOTNuMuBar << endl;
  cout << "Exp/SKMC >> " << "POT {MC; nue_bkg} = " << fMCPOTNuE     << endl;
  cout << "Exp/SKMC >> " << "POT {MC; nue_sig} = " << fMCPOTNuESig  << endl;

  cout << "Exp/SKMC >> " << "Nev {MC; numu;    1Rmu} = " << fMCTreeNuMu1Rmu    -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; numubar; 1Rmu} = " << fMCTreeNuMuBar1Rmu -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; nue_bkg; 1Rmu} = " << fMCTreeNuE1Rmu     -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; nue_sig; 1Rmu} = " << fMCTreeNuESig1Rmu  -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; numu;    1Re } = " << fMCTreeNuMu1Re     -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; numubar; 1Re } = " << fMCTreeNuMuBar1Re  -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; nue_bkg; 1Re } = " << fMCTreeNuE1Re      -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; nue_sig; 1Re } = " << fMCTreeNuESig1Re   -> GetEntries()  << endl;
#endif
  
  //
  // Systematics
  //

  // Get relative error on energy scale
  fEnergyScale1RmuRelErr = 
    util::XML2Dbl(xml,root_node,
     "detector_systematics/superk/energy_scale_relative_error/one_mu_ring");
  // Get relative error on various efficiencies
  fEff1RmuCCQEOthRelErr =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/other_error");
  fEff1RmuCCQERCCovarRelErr = 
    util::XML2TMatrixD(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/rc_error");
  fEff1RmuCCQERCCovarE1RelErr =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/rc_ereco_high_bin1");
  fEff1RmuCCQERCCovarE2RelErr =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccqe/rc_ereco_high_bin2");
  fEff1RmuCCnonQERelErr   =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/ccnonqe");
  fEff1RmuNCRelErr        =
    util::XML2Dbl(xml,root_node,
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/nc");
  fEff1RmuNuECCRelErr =
    util::XML2Dbl(xml,root_node, 
      "detector_systematics/superk/efficiency_relative_error/one_mu_ring/nuecc");

  if(!fEff1RmuCCQERCCovarRelErr) {
    cerr << "Exp/SK >> No input SuperK CCQE ring-counting error covariance matrix!" << endl;
    exit(1);
  }
      
#ifdef __show_init_mesg_info__      
  cout << "Exp/SKMC >> " << "Error {EnergyScale;1Rmu}    : " << fEnergyScale1RmuRelErr  << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency;1Rmu/CCQE/Oth} : " << fEff1RmuCCQEOthRelErr << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency;1Rmu/CCQE/RC}  : " << endl;
  cout << " - The ring-counting efficiency error is given for the following 3 Ereco bins : ";
  cout << "(0, " << fEff1RmuCCQERCCovarE1RelErr << " GeV), "
       << "[" << fEff1RmuCCQERCCovarE1RelErr << " GeV, " << fEff1RmuCCQERCCovarE2RelErr << " GeV], "
       << "(" << fEff1RmuCCQERCCovarE2RelErr << " GeV, oo)" << endl;
  cout << " - The ring-counting error matrix is : " << endl;
  util::PrintCompact(*fEff1RmuCCQERCCovarRelErr,10,6);
  cout << "Exp/SKMC >> " << "Error {Efficiency;1Rmu/CCnonQE} : " << fEff1RmuCCnonQERelErr  << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency;1Rmu/NC} : "      << fEff1RmuNCRelErr       << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency;1Rmu/NuECC} : "   << fEff1RmuNuECCRelErr    << endl;
#endif
  
  // initialise physics hypothesis
  fPhysHypothesis.Init(kHypUndefined);

  // fill-in MC templates from the input event tree
  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     for(int im = 0; im < knm; im++) {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       this->FillTemplate(s,m);
     }
  }

  // apply flux tuning (based on external data) to templates
  this->ApplyFluxExtTuningToTemplates();
     
  // apply flux and cross-section tuning (based on nd280 data fit by BANFF) to templates
  this->ApplyBANFFTuningToTemplates();

#ifdef __show_init_mesg_info__      
  cout << "Exp/SKMC >> Done reading SuperK inputs!" << endl;
#endif
  
  return true;
}
//________________________________________________________________________________________
void SuperkMC::Book(void)
{
  //
  // Ereco histograms filled-in by SK spectrum prediction code
  // Histograms corresponding to different samples can have different binning
  //

  for(int is = 0; is < kns; is++) 
  {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     const int      n_ereco_bins    = util::NumErecoBins (s);
     const double * ereco_bin_edges = util::ErecoBinEdges(s);
     fMCErecoSpectrumTot[is] = new TH1D("","", n_ereco_bins, ereco_bin_edges);
     fMCErecoSpectrumTot[is] -> SetDirectory(0);
     fMCErecoSpectrumTot[is]->SetLineWidth(2);
     fMCErecoSpectrumTot[is]->SetLineColor(kBlack);
     for(int im=0; im < knm; im++) 
     {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       fMCErecoSpectrum[is][im] = new TH1D("","", n_ereco_bins, ereco_bin_edges);
       fMCErecoSpectrum[is][im] -> SetDirectory(0);
       util::SetDefaultSKSpectrumStyle(s,m,fMCErecoSpectrum[is][im]);
     }//m
  }//s

  //
  // 2-D(Ereco,Etrue) MC templates
  // Histograms corresponding to different samples can have different binning
  //

  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     const int      n_ereco_bins    = util::NumErecoBins (s);
     const double * ereco_bin_edges = util::ErecoBinEdges(s);
     for(int im = 0; im < knm; im++) {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       fMCTemplNom[is][im] = new TH2D(
          Form("nsk_nom_%d_%d",is,im), util::AsString(m).c_str(),
               n_ereco_bins, ereco_bin_edges, kNEtrueBins, kEtrueMin, kEtrueMax);
       fMCTemplTun[is][im] = new TH2D(
          Form("nsk_tun_%d_%d",is,im), util::AsString(m).c_str(),
               n_ereco_bins, ereco_bin_edges, kNEtrueBins, kEtrueMin, kEtrueMax);
       fMCTemplSys[is][im] = new TH2D(
          Form("nsk_sys_%d_%d",is,im), util::AsString(m).c_str(),
               n_ereco_bins, ereco_bin_edges, kNEtrueBins, kEtrueMin, kEtrueMax);
       fMCTemplFin[is][im] = new TH2D(
          Form("nsk_fin_%d_%d",is,im), util::AsString(m).c_str(),
               n_ereco_bins, ereco_bin_edges, kNEtrueBins, kEtrueMin, kEtrueMax);
       fMCTemplNom[is][im]->SetDirectory(0);
       fMCTemplNom[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fMCTemplNom[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
       fMCTemplTun[is][im]->SetDirectory(0);
       fMCTemplTun[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fMCTemplTun[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
       fMCTemplSys[is][im]->SetDirectory(0);
       fMCTemplSys[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fMCTemplSys[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
       fMCTemplFin[is][im]->SetDirectory(0);
       fMCTemplFin[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fMCTemplFin[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
     }
  }

}
//________________________________________________________________________________________
void SuperkMC::FillTemplate(SuperkSample_t s, SuperkMode_t m)
{
#ifdef __show_init_mesg_debug__      
  cout << "Exp/SKMC/CreateTemplate >> Creating MC template for: "
       << util::AsString(s) << " / " << util::AsString(m) << endl;
#endif
    
  TTree * mctree = this->MCTree(s,m);
  double  mcpot  = this->MCPOT (s,m);

  if(!mctree) {
     cerr << "Exp/SKMC/CreateTemplate >> Couldn't find event tree for template: " 
          << util::AsString(s) << "/" << util::AsString(m) << endl;
     exit(1);
  }

  double pot_wght = this->AnalyzedSamplePOT() / mcpot;

#ifdef __show_init_mesg_debug__        
  cout << "Exp/SKMC/CreateTemplate >>      "
       << "Weight (data POT / MC POT) = " << pot_wght << endl;
#endif
  if(pot_wght <= 0) {
     cerr << "Exp/SKMC/CreateTemplate >> Weight (data POT / MC POT) for " 
          << util::AsString(s) << "/" << util::AsString(m)
          << " is " << pot_wght << endl;
     exit(1);
  }

  TBranch * branchEreco   = mctree->GetBranch("EnuReco");    
  TBranch * branchEtrue   = mctree->GetBranch("Enu");
  TBranch * branchPdg     = mctree->GetBranch("Pdg");
  TBranch * branchGenCode = mctree->GetBranch("Mode");

  branchEreco   -> SetAddress(&fBrEreco);
  branchEtrue   -> SetAddress(&fBrEtrue);
  branchPdg     -> SetAddress(&fBrPdg);
  branchGenCode -> SetAddress(&fBrGenCode);
  
  for(int iev = 0; iev < mctree->GetEntries(); iev++) 
  {
     mctree->GetEntry(iev);

     // use current entry for mode m template?
     bool use_it = util::UseEventToFillInTemplate(m,fBrPdg,fBrGenCode);
     if(!use_it) continue;
     
     this -> Template(s,m,kSKMCtNom) -> Fill (fBrEreco, fBrEtrue, pot_wght);
     this -> Template(s,m,kSKMCtTun) -> Fill (fBrEreco, fBrEtrue, pot_wght);
     this -> Template(s,m,kSKMCtSys) -> Fill (fBrEreco, fBrEtrue, pot_wght);
     this -> Template(s,m,kSKMCtFin) -> Fill (fBrEreco, fBrEtrue, pot_wght);

  }//iev

#ifdef __show_init_mesg_debug__        
  cout << "Exp/SKMC/CreateTemplate >>      "
       << "Number of template entries = " 
       << this -> Template(s,m,kSKMCtNom) -> GetEntries() << endl;
#endif
}
//________________________________________________________________________________________
TH1D* SuperkMC::CurrErecoSpectrum(SuperkSample_t s) const                               
{ 
  int sidx = util::FindSampleSlot(s);
  if(sidx == -1) {
    return 0;
  }
  return fMCErecoSpectrumTot[sidx]; 
}
//________________________________________________________________________________________
TH1D* SuperkMC::CurrErecoSpectrum(SuperkSample_t s, SuperkMode_t m) const 
{ 
  int sidx = util::FindSampleSlot(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::FindModeSlot(s,m);
  if(midx == -1) {
    return 0;
  }            
  return fMCErecoSpectrum[sidx][midx]; 
}
//________________________________________________________________________________________
TH2D * SuperkMC::Template(SuperkSample_t s, SuperkMode_t m, SKMCTemplateType_t type) const
{
  int sidx = util::FindSampleSlot(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::FindModeSlot(s,m);
  if(midx == -1) {
    return 0;
  }            

  TH2D * h = 0;

  switch(type) {
   case (kSKMCtNom) : 
     h = fMCTemplNom[sidx][midx]; 
     break;
   case (kSKMCtTun) : 
     h = fMCTemplTun[sidx][midx];
     break;
   case (kSKMCtSys) : 
     h = fMCTemplSys[sidx][midx];
     break;
   case (kSKMCtFin) : 
     h = fMCTemplFin[sidx][midx];
     break;
   default:
     return 0;
  }
  assert(h);

  return h;
}
//________________________________________________________________________________________
void SuperkMC::ApplyFluxExtTuningToTemplates(void) 
{
// Take into account the flux tuning based on external (NA61,...) data.
// This method reweights/overwrites the nominal templates & it is run just once after
// filling-in the MC templates from the input mini-DSTs.

  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    // loop over sample modes
    for(int im = 0; im < knm; im++) { 
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       // consider template?
       bool skip = this->SkipTemplate(m);
       if(skip) continue;
       // get corresponding template
       TH2D * mctempl = this -> Template(s,m,kSKMCtNom);
       if(!mctempl) continue; 
       // loop over Ereco
       for(int ir = 1; ir <= mctempl->GetXaxis()->GetNbins(); ir++) {
         //double Ereco = mctempl->GetXaxis()->GetBinCenter(ir);
         // loop over Etrue
         for(int it = 1; it <= mctempl->GetYaxis()->GetNbins(); it++) {
           double Etrue = mctempl->GetYaxis()->GetBinCenter(it);
           // get number of events
           double nev = mctempl->GetBinContent(ir,it);
           // get tuning weights
           int flux_nu_pdg = util::SuperkMode2FluxNuPdg(m);
           double wght = fNuFlux->SuperkFluxExtTuningWght(flux_nu_pdg,Etrue);
           // update nominal MC template
           nev *= wght;
           mctempl->SetBinContent(ir,it,nev);
           // initialize "tun", "sys" & "fin" templates to "nom" values
           this->Template(s,m,kSKMCtTun)->SetBinContent(ir,it,nev);
           this->Template(s,m,kSKMCtSys)->SetBinContent(ir,it,nev);
           this->Template(s,m,kSKMCtFin)->SetBinContent(ir,it,nev);
         }//it
       }//ir
    }//im
  }//is
}
//________________________________________________________________________________________
void SuperkMC::ApplyBANFFTuningToTemplates(void)
{
// Take BANFF tuning into account (reweight MC to BANFF best-fit parameters)
// The weights returned by NuFlux::SuperkFluxFinTuningWght() are relative to the flux
// tuned based on external (NA61,...) data, so make sure that you run 
// SuperkMC::ApplyFluxExtTuningToTemplates() before running this method.
// Source template type: kSKMCtNom
// Target template type: kSKMCtTun
//

  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    // loop over sample modes
    for(int im = 0; im < knm; im++) { 
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       // consider template?
       bool skip = this->SkipTemplate(m);
       if(skip) continue;
       // get corresponding template
       TH2D * mctempl_src = this -> Template(s,m,kSKMCtNom);
       if(!mctempl_src) continue; 
       TH2D * mctempl_tgt = this -> Template(s,m,kSKMCtTun);
       if(!mctempl_tgt) continue; 
       // loop over Ereco
       for(int ir = 1; ir <= mctempl_src->GetXaxis()->GetNbins(); ir++) {
         //double Ereco = mctempl_src->GetXaxis()->GetBinCenter(ir);
         // loop over Etrue
         for(int it = 1; it <= mctempl_src->GetYaxis()->GetNbins(); it++) {
           double Etrue = mctempl_src->GetYaxis()->GetBinCenter(it);
           // get number of events
           double nev = mctempl_src->GetBinContent(ir,it);
           // get tuning weights
           int flux_nu_pdg = util::SuperkMode2FluxNuPdg(m);
           double flux_wght = fNuFlux->SuperkFluxFinTuningWght(flux_nu_pdg,Etrue);
           double xsec_wght = 1.; /// TODO            
           // update tuned MC template
           double nev_tuned = nev * flux_wght * xsec_wght;
           mctempl_tgt->SetBinContent(ir,it,nev_tuned);
           // initialize "sys" & "fin" templates to "tun" values
           this->Template(s,m,kSKMCtSys)->SetBinContent(ir,it,nev);
           this->Template(s,m,kSKMCtFin)->SetBinContent(ir,it,nev);
         }//it
       }//ir
    }//im
  }//is
}
//________________________________________________________________________________________
void SuperkMC::UpdateTemplates(void)
{
  this->ApplySystematicsToTemplates   ();
  this->ApplyPhysHypothesisToTemplates();
}
//________________________________________________________________________________________
void SuperkMC::ApplySystematicsToTemplates(void) 
{
// Apply systematic variations to MC templates
//
  this->ApplySystematicsToTemplates_Multiplicative   (); 
  this->ApplySystematicsToTemplates_NonMultiplicative();
}
//________________________________________________________________________________________
void SuperkMC::ApplySystematicsToTemplates_Multiplicative(void) 
{
// All multiplicative systematics are handled in the following loop.
// Basically all systematic variations for each MC template bin is computed 
// in SuperkMC::NumOfEvents(Ereco,Etrue,sample,mode)
//

  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    // loop over sample modes
    for(int im = 0; im < knm; im++) { 
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       // consider template?
       bool skip = this->SkipTemplate(m);
       if(skip) continue;
       // get corresponding templates
       TH2D * mctempl_src = this -> Template(s,m,kSKMCtTun);
       if(!mctempl_src) continue; 
       TH2D * mctempl_tgt = this -> Template(s,m,kSKMCtSys);
       if(!mctempl_tgt) continue; 
       // loop over Ereco
       for(int ir = 1; ir <= mctempl_src->GetXaxis()->GetNbins(); ir++) {
         double Ereco = mctempl_src->GetXaxis()->GetBinCenter(ir);
         // loop over Etrue
         for(int it = 1; it <= mctempl_src->GetYaxis()->GetNbins(); it++) {
           double Etrue = mctempl_src->GetYaxis()->GetBinCenter(it);
           // get number of events in current bin of source MC template
           double nev = mctempl_src->GetBinContent(ir,it);
           // get combined weight due to all multiplicative systematic variations
           double wght = this->SystVariationWeight (Ereco, Etrue, s,m);
           // update target MC template
           mctempl_tgt->SetBinContent(ir,it,nev*wght);
         }//it
       }//ir
    }//im
  }//is
}
//________________________________________________________________________________________
void SuperkMC::ApplySystematicsToTemplates_NonMultiplicative(void) 
{
//
// Now that all simple multiplicative systematics were applied, apply the reconstructed 
// energy scale systematic. This systematic is applied by scaling the Ereco bin edges
// and, by assuming uniform distribution of events within each bin, calculating how many
// events a lost to / gained from neighbouring bins.
//

  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    // loop over sample modes
    for(int im = 0; im < knm; im++) { 
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       // consider template?
       bool skip = this->SkipTemplate(m);
       if(skip) continue;
       // get corresponding template
       TH2D * mctempl = this -> Template(s,m,kSKMCtSys);
       if(!mctempl) continue; 
       // clone
       TH2D * mctempl_no_ereco_scale = new TH2D(*mctempl);
       // get energy scale
       double escale = this->EnergyScale(s,m);

       // energy scale > 1 : 
       // For a given bin, events are gained from bins on the left (lower Ereco bins) and 
       // are lost to bins on the right (higher Ereco bins)
       if(escale > 1. + kASmallNum) {
           // loop over Ereco
           for(int ir = 1; ir <= mctempl->GetXaxis()->GetNbins(); ir++) {
              // nominal reco bin edges
              double bin_max = mctempl->GetXaxis()->GetBinUpEdge (ir);
              double bin_min = mctempl->GetXaxis()->GetBinLowEdge(ir);
              // scaled reco bin edges
              bin_max *= escale;
              bin_min *= escale;
              // loop over Ereco' (> Ereco)
              for(int irp = ir+1; irp <= mctempl->GetXaxis()->GetNbins(); irp++) {
                 // bin edges for reco bin ir'
                 double bin_max_p = mctempl->GetXaxis()->GetBinUpEdge (irp);
                 double bin_min_p = mctempl->GetXaxis()->GetBinLowEdge(irp);
                 // do bins ir and ir' overlap?
                 if(bin_max <= bin_min_p || bin_min >= bin_max_p) continue; // no overlap
                 // calculate overlap
                 double overlap = 0;
                 if(bin_min_p < bin_max && bin_min_p > bin_min) overlap = bin_max   - bin_min_p;
                 if(bin_max_p < bin_max && bin_max_p > bin_min) overlap = bin_max_p - bin_min;
                 if(bin_max_p > bin_max && bin_min_p < bin_min) overlap = bin_max   - bin_min;
                 if(bin_max_p < bin_max && bin_min_p > bin_min) overlap = bin_max_p - bin_min_p;                 
                 // loop over Etrue
                 for(int it = 1; it <= mctempl->GetYaxis()->GetNbins(); it++) {
                      // calculate how many events from nominal bin ir get lost to bin ir'
                      double n  = mctempl_no_ereco_scale->GetBinContent(ir, it);
                      double dn = n * overlap / (bin_max-bin_min); ///< (orerlap) * (nominal # of events) / (bin width)
                      dn = TMath::Min(dn,n);
                      // update the systematically varied templates
                      mctempl->SetBinContent(ir, it, mctempl->GetBinContent(ir,  it) - dn);
                      mctempl->SetBinContent(irp,it, mctempl->GetBinContent(irp, it) + dn);
                 }//it
              }//irp
           }//ir
       }// scale > 1

       // energy scale < 1 : 
       // For a given bin, events are gained from bins on the right (higher Ereco bins) and 
       // are lost to bins on the left (lower Ereco bins)
       if(escale < 1. - kASmallNum) {
           // loop over Ereco
           for(int ir = 1; ir <= mctempl->GetXaxis()->GetNbins(); ir++) {
              // nominal reco bin edges
              double bin_max = mctempl->GetXaxis()->GetBinUpEdge (ir);
              double bin_min = mctempl->GetXaxis()->GetBinLowEdge(ir);
              // scaled reco bin edges
              bin_max *= escale;
              bin_min *= escale;
              // loop over Ereco' (< Ereco)
              for(int irp = 1; irp <= ir-1; irp++) {
                 // bin edges for reco bin ir'
                 double bin_max_p = mctempl->GetXaxis()->GetBinUpEdge (irp);
                 double bin_min_p = mctempl->GetXaxis()->GetBinLowEdge(irp);
                 // do bins ir and ir' overlap?
                 if(bin_max <= bin_min_p || bin_min >= bin_max_p) continue; // no overlap
                 // calculate overlap
                 double overlap = 0;
                 if(bin_min_p < bin_max && bin_min_p > bin_min) overlap = bin_max   - bin_min_p;
                 if(bin_max_p < bin_max && bin_max_p > bin_min) overlap = bin_max_p - bin_min;
                 if(bin_max_p > bin_max && bin_min_p < bin_min) overlap = bin_max   - bin_min;
                 if(bin_max_p < bin_max && bin_min_p > bin_min) overlap = bin_max_p - bin_min_p;
                 // loop over Etrue
                 for(int it = 1; it <= mctempl->GetYaxis()->GetNbins(); it++) {
                      // calculate how many events from nominal bin ir get lost to bin ir'
                      double n  = mctempl_no_ereco_scale->GetBinContent(ir, it);
                      double dn = n * overlap / (bin_max-bin_min); ///< (orerlap) * (nominal # of events) / (bin width)
                      dn = TMath::Min(dn,n);
                      // update the systematically varied templates
                      mctempl->SetBinContent(ir, it, mctempl->GetBinContent(ir,  it) - dn);
                      mctempl->SetBinContent(irp,it, mctempl->GetBinContent(irp, it) + dn);
                 }//it
              }//irp
           }//ir
       }// scale < 1
       delete mctempl_no_ereco_scale;
    }//im
  }//is

}
//________________________________________________________________________________________
void SuperkMC::ApplyPhysHypothesisToTemplates(void)
{
// Apply physics hypothesis {oscillation (2-flavour / 3-flavour), decay or 
// decoherence} to MC templates
//

  PhysHypothesis_t hyp = fPhysHypothesis.Hypothesis();

  //
  // 3-flavour oscillations
  //
  if(hyp == kHypOsc3F) 
  {
    // set oscillation parameters to the oscillation probability calculation code
    double mixing12 = fPhysHypothesis.CurrValue(kPhypMixing12);
    double mixing13 = fPhysHypothesis.CurrValue(kPhypMixing13);
    double mixing23 = fPhysHypothesis.CurrValue(kPhypMixing23);
    double dm12sq   = fPhysHypothesis.CurrValue(kPhypDm12sq);
    double dm23sq   = fPhysHypothesis.CurrValue(kPhypDm23sq);
    double delta_cp = fPhysHypothesis.CurrValue(kPhypDeltaCP);
#ifdef __show_hypothesis__
    cout << "**  Exp/SKMC >> " << fPhysHypothesis;
#endif
    NuOscMixParam_t omp = Opt::Instance()->gOscMixingParamChoice;
    switch(omp) {
       // mixing params are sin^2(2*theta)
       case(kNuMixSinsq2Theta) : 
         fNuOsc3Flavour.SetOscParameters(mixing12,mixing13,mixing23,dm12sq,dm23sq,delta_cp);
         break;
       // mixing params are sin^2(theta)
       case(kNuMixSinsqTheta) : 
         fNuOsc3Flavour.SetOscParametersSinSqTAll(mixing12,mixing13,mixing23,dm12sq,dm23sq,delta_cp);
         break;
       default : 
         cerr << "**  Exp/SKMC >> Can't handle choice of oscillation mixing params (" 
              << util::AsString(omp) << ")" << endl;
         exit(1);
    }
    
    // loop over samples
    for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;

     assert(this -> Template(s,kSKmNuMuCCQE,kSKMCtSys));
     TAxis * ereco_axis = this->Template(s,kSKmNuMuCCQE,kSKMCtSys)->GetXaxis(); 
     TAxis * etrue_axis = this->Template(s,kSKmNuMuCCQE,kSKMCtSys)->GetYaxis(); 

     // loop over Etrue
     for(int it = 1; it <= etrue_axis->GetNbins(); it++) {

        // calculate oscillation probabilities (incl. matter effects) for the current Etrue value
        double Etrue = etrue_axis->GetBinCenter(it);

        double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuMu,    kPdgNuMu,    Etrue);
      //double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuE,     kPdgNuMu,    Etrue);
        double p_numubar_numubar = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuMuBar, kPdgNuMuBar, Etrue);
        double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuE,     kPdgNuE,     Etrue);
        double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuMu,    kPdgNuE,     Etrue);
#ifdef __show_hypothesis_probabilities__
        cout << "Exp/SKMC >>  E = " << Etrue << " GeV" << endl;
        cout << "Exp/SKMC >>  * P (numu     -> numu)    = " << p_numu_numu       << endl;
      //cout << "Exp/SKMC >>  * P (nue      -> numu)    = " << p_nue_numu        << endl;
        cout << "Exp/SKMC >>  * P (numubar  -> numubar) = " << p_numubar_numubar << endl;
        cout << "Exp/SKMC >>  * P (nue      -> nue)     = " << p_nue_nue         << endl;
        cout << "Exp/SKMC >>  * P (numu     -> nue)     = " << p_numu_nue        << endl;
#endif

        // loop over Ereco
        for(int ir = 1; ir <= ereco_axis->GetNbins(); ir++) {

           // loop over modes
           for(int im=0; im < knm; im++) {
              SuperkMode_t m = kSKReactionModes[is][im];
              if(m == kSKmUndefined) continue;
              // number of events in (Ereco,Etrue) bin (systematics applied)
              double nev = this -> Template(s,m,kSKMCtSys) -> GetBinContent(ir,it);
              // osc. probability
              double posc = 1.;
              if( util::IsCC(m) ) 
              {
                if( util::IsNuMu(m)      ) { posc = p_numu_numu;       }
                if( util::IsNuMuBar(m)   ) { posc = p_numubar_numubar; }
                if( util::IsBeamNuE(m)   ) { posc = p_nue_nue;         }
                if( util::IsSignalNuE(m) ) { posc = p_numu_nue;        }
              }
              // number of events in (Ereco,Etrue) bin (oscillations + systematics applied)
              double nev_osc = nev * posc;
              this -> Template(s,m,kSKMCtFin) -> SetBinContent(ir,it, nev_osc);
           }//im
       }//ir
     }//it
   }//is
  }

  //
  // 2-flavour oscillations / numu disappearance
  //
  else
  if(hyp == kHypNuMuDisapp2F) {
    double mixing23 = fPhysHypothesis.CurrValue(kPhypMixing23);
    double dm23sq   = fPhysHypothesis.CurrValue(kPhypDm23sq);
    // loop over samples
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      if(s != kSKsFCFV1RingMuLike) continue;
      // loop over sample modes
      for(int im = 0; im < knm; im++) { 
         SuperkMode_t m = (SuperkMode_t) im;
         // consider template?
         bool skip = this->SkipTemplate(m);
         if(skip) continue;
         // get source (systematics) and target (systematics + phys.hypothesis) template
         TH2D * mctempl_src = this -> Template(s,m,kSKMCtSys);
         if(!mctempl_src) continue; 
         TH2D * mctempl_tgt = this -> Template(s,m,kSKMCtFin);
         if(!mctempl_tgt) continue; 
         // loop over Etrue
         for(int it = 1; it <= mctempl_src->GetYaxis()->GetNbins(); it++) {
           double Etrue = mctempl_src->GetYaxis()->GetBinCenter(it);
           double psurv = util::NuMuSurvProbOsc2Flavour(mixing23, dm23sq, Etrue, m);
           // loop over Ereco
           for(int ir = 1; ir <= mctempl_src->GetXaxis()->GetNbins(); ir++) {
              double nev = mctempl_src->GetBinContent(ir,it);
              nev *= psurv;
              mctempl_tgt->SetBinContent(ir,it,nev);
           }//ir
         }//it
      }//im
    }//is
  }

  //
  // decay
  //
  else
  if(hyp == kHypDecay) {
    double sinsq2theta = fPhysHypothesis.CurrValue(kPhypSinsq2ThetaDecay);
    double alpha       = fPhysHypothesis.CurrValue(kPhypAlphaDecay);
    // loop over samples
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      if(s != kSKsFCFV1RingMuLike) continue;
      // loop over sample modes
      for(int im = 0; im < knm; im++) { 
         SuperkMode_t m = (SuperkMode_t) im;
         // consider template?
         bool skip = this->SkipTemplate(m);
         if(skip) continue;
         // get source (systematics) and target (systematics + phys.hypothesis) template
         TH2D * mctempl_src = this -> Template(s,m,kSKMCtSys);
         if(!mctempl_src) continue; 
         TH2D * mctempl_tgt = this -> Template(s,m,kSKMCtFin);
         if(!mctempl_tgt) continue; 
         // loop over Etrue
         for(int it = 1; it <= mctempl_src->GetYaxis()->GetNbins(); it++) {
           double Etrue = mctempl_src->GetYaxis()->GetBinCenter(it);
           double psurv = util::NuMuSurvProbDecay(sinsq2theta, alpha, Etrue, m);
           // loop over Ereco
           for(int ir = 1; ir <= mctempl_src->GetXaxis()->GetNbins(); ir++) {
              double nev = mctempl_src->GetBinContent(ir,it);
              nev *= psurv;
              mctempl_tgt->SetBinContent(ir,it,nev);
           }//ir
         }//it
      }//im
    }//is
  }


  //
  // decoherence
  //
  else
  if(hyp == kHypDecoh) {
    double sinsq2theta = fPhysHypothesis.CurrValue(kPhypSinsq2ThetaDecoh);
    double musq        = fPhysHypothesis.CurrValue(kPhypMusqDecoh);
    // loop over samples
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = kSKSamples[is];
      if(s != kSKsFCFV1RingMuLike) continue;
      // loop over sample modes
      for(int im = 0; im < knm; im++) { 
         SuperkMode_t m = (SuperkMode_t) im;
         // consider template?
         bool skip = this->SkipTemplate(m);
         if(skip) continue;
         // get source (systematics) and target (systematics + phys.hypothesis) template
         TH2D * mctempl_src = this -> Template(s,m,kSKMCtSys);
         if(!mctempl_src) continue; 
         TH2D * mctempl_tgt = this -> Template(s,m,kSKMCtFin);
         if(!mctempl_tgt) continue; 
         // loop over Etrue
         for(int it = 1; it <= mctempl_src->GetYaxis()->GetNbins(); it++) {
           double Etrue = mctempl_src->GetYaxis()->GetBinCenter(it);
           double psurv = util::NuMuSurvProbDecoh(sinsq2theta, musq, Etrue, m);
           // loop over Ereco
           for(int ir = 1; ir <= mctempl_src->GetXaxis()->GetNbins(); ir++) {
              double nev = mctempl_src->GetBinContent(ir,it);
              nev *= psurv;
              mctempl_tgt->SetBinContent(ir,it,nev);
           }//ir
         }//it
      }//im
    }//is
  }

}
//________________________________________________________________________________________
double SuperkMC::SystVariationWeight(
    double Ereco, double Etrue, SuperkSample_t s, SuperkMode_t m) const
{
// The return value is the combined wght for the (Ereco,Etrue) bin of the (s,m) template
// due to all multiplicative systematic param variation
// This method considers:
// - SuperK efficiencies variations
// - BANFF parameter variations
// - Cross-section variations

#ifdef __show_syst_var_wght_mesg_info__ 
  cout << "Exp/SKMC/SystVariationWeight >> "
       << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m) 
       << ", Ereco = " << Ereco << " GeV (bin = " << reco_bin << ")" 
       << ", Etrue = " << Etrue << " GeV (bin = " << true_bin << ")" << endl;
#endif

  //
  // Plug-in BANFF flux systematics
  //
  int flux_nu_pdg = util::SuperkMode2FluxNuPdg(m);
  double banff_flux_wght = fNuFlux->SuperkFluxSystWght(flux_nu_pdg, Etrue);

  //
  // Plug-in BANFF xsec systematics 
  //
  double banff_xsec_wght = 1;

  //
  // ... other xsec, FSI etc
  //

  //
  // Plug-in SK efficiency uncertainty due to fiducial volume, reduction chain, OD cut,
  // ring-counting, PID, momentum cut and decay e- cut.
  // The ring-counting error for CCQE has Ereco-dependence.
  //
  double skeff_wght = this->EfficiencyWeight(s,m,Ereco);
#ifdef __show_syst_var_wght_mesg_info__
  if(skeff_wght != 1.) {
    cout << "Exp/SKMC/SystVariationWeight >> Applied efficiency-systematic weight = " 
         << skeff_wght << endl;
  }
#endif

  double wght = 
    banff_flux_wght *
    banff_xsec_wght *
    skeff_wght;

  return wght;
}
//________________________________________________________________________________________
void SuperkMC::UpdateErecoSpectra(void)
{
  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    TH1D * reco_spectrum_tot = this->CurrErecoSpectrum(s);

    // loop over Ereco
    for(int ir = 1; ir <= reco_spectrum_tot->GetNbinsX(); ir++) {
       double nev_tot = 0;

       // loop over modes
       for(int im = 0; im < knm; im++) { 
          SuperkMode_t m = kSKReactionModes[is][im];
          if(m == kSKmUndefined) continue;
          // consider template?
          bool skip = this->SkipTemplate(m);
          if(skip) continue;

          TH2D * mctempl = this -> Template(s,m,kSKMCtFin);
          TH1D * reco_spectrum = this->CurrErecoSpectrum(s,m);

          // loop over Etrue
          double nev = 0;
          for(int it = 1; it <= mctempl->GetYaxis()->GetNbins(); it++) {
             nev += mctempl->GetBinContent(ir,it);
          }//it
          reco_spectrum->SetBinContent(ir,nev);
          nev_tot += nev;
       }//im
       reco_spectrum_tot->SetBinContent(ir,nev_tot);
    }//ir
  }//is
}
//________________________________________________________________________________________
double SuperkMC::AnalyzedSamplePOT(void) const
{
  return fPOT;  
}
//________________________________________________________________________________________
TTree * SuperkMC::MCTree(SuperkSample_t s, SuperkMode_t m)
{
  if (s==kSKsFCFV1RingMuLike)
  {
     if ( util::IsNuMu(m)      ) return fMCTreeNuMu1Rmu;
     if ( util::IsNuMuBar(m)   ) return fMCTreeNuMuBar1Rmu;
     if ( util::IsBeamNuE(m)   ) return fMCTreeNuE1Rmu;
     if ( util::IsSignalNuE(m) ) return fMCTreeNuESig1Rmu;
  }//1Rmu

  if (s==kSKsFCFV1RingELike)
  {
     if ( util::IsNuMu(m)      ) return fMCTreeNuMu1Re;
     if ( util::IsNuMuBar(m)   ) return fMCTreeNuMuBar1Re;
     if ( util::IsBeamNuE(m)   ) return fMCTreeNuE1Re;
     if ( util::IsSignalNuE(m) ) return fMCTreeNuESig1Re;
  }//1Re

  return 0;
}
//________________________________________________________________________________________
double SuperkMC::MCPOT(SuperkSample_t s, SuperkMode_t m)
{
  if(s == kSKsFCFV1RingMuLike || s == kSKsFCFV1RingELike)
  {
     if ( util::IsNuMu(m)      ) return fMCPOTNuMu;
     if ( util::IsNuMuBar(m)   ) return fMCPOTNuMuBar;
     if ( util::IsBeamNuE(m)   ) return fMCPOTNuE;
     if ( util::IsSignalNuE(m) ) return fMCPOTNuESig;
  }
  return 0;
}
//________________________________________________________________________________________
void SuperkMC::TweakEnergyScale(double twkdial) 
{ 
  fEnergyScaleTwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkMC::TweakEfficiency1RmuCCQEOther(double twkdial) 
{ 
  fEff1RmuCCQEOthTwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkMC::TweakEfficiency1RmuCCQERingCounting(
   bool   force_input_twkdials, /* Force input twk dials? If not, generate using the Cholesky method */
   double twkdial0,             /* Value of 1st twk dial, needed if values are forced */
   double twkdial1,             /* Value of 2nd twk dial, needed if values are forced */
   double twkdial2,             /* Value of 3rd twk dial, needed if values are forced */
   bool   corellated            /* Are forced twk dials the correlated or uncorrelated ones */
)
{
#ifdef __show_rc_cholesky_mesg_debug__
  cout << "Exp/SKMC >> Tweaking CCQE ring-counting efficiency" << endl;
#endif

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the ring-counting error matrix (covariance matrix)
  //
  if(!fEff1RmuCCQERCCovarChoLRelErr) 
  {
#ifdef __show_rc_cholesky_mesg_info__
    cout << "Exp/SKMC/TweakRingCountingEffCCQE >> "
         << "Run Cholesky's decomposition method on the CCQE RC error matrix." << endl;
#endif

    TMatrixD ccqe_rc_err_matrix(*fEff1RmuCCQERCCovarRelErr);
#ifdef __show_rc_cholesky_mesg_debug__
    cout << "Exp/SKMC/TweakRingCountingEffCCQE >> RC error matrix:" << endl;
    ccqe_rc_err_matrix.Print();
#endif

    TMatrixD L = util::CholeskyDecomposition(ccqe_rc_err_matrix);
    fEff1RmuCCQERCCovarChoLRelErr = new TMatrixD(L);
#ifdef __show_rc_cholesky_mesg_info__
    cout << "Exp/SKMC/TweakRingCountingEffCCQE >> "
         << "Cholesky decomposition of RC error matrix:" << endl;
    fEff1RmuCCQERCCovarChoLRelErr->Print();
#endif
  }


  //
  // Check whether the values of the tweaking dials are forced. 
  // If the forced tweaking dials are supposed to be the correlated ones then we don't
  // need to use the Cholesky decomposition of the  ring-counting error matrix.
  //
  if(force_input_twkdials && corellated) {
    fEff1RmuCCQERC0TwkDial = twkdial0;
    fEff1RmuCCQERC1TwkDial = twkdial1;
    fEff1RmuCCQERC2TwkDial = twkdial2;
    return;
  }

  TMatrixD & L = *fEff1RmuCCQERCCovarChoLRelErr;
  TVectorD twkv(3);

  //
  // Check whether the values of the tweaking dials are forced. 
  // If the forced tweaking dials are supposed to be the un-correlated ones then use
  // the Cholesky decomposition of the ring-counting error matrix to generate
  // correlated tweaking dials.
  //

  if(force_input_twkdials && !corellated) {
    TVectorD twkv_uncorrelated(3);
    twkv_uncorrelated[0] = twkdial0;
    twkv_uncorrelated[1] = twkdial1;
    twkv_uncorrelated[2] = twkdial2;
    twkv = util::CholeskyCalculateCorrelatedParamVariations(L,twkv_uncorrelated);
  } 

  //
  // If the tweaking dials were not forced, generate gaussian correlated tweaking dials.
  // Then use the Cholesky decomposition of the ring-counting error matrix to generate
  // correlated tweaking dials.
  //

  else {
    twkv = util::CholeskyGenerateCorrelatedParamVariations(L);
  }

  //
  // These generated tweaking dials are distributed around 0 and have the spread and 
  // correlation described by the input fEffCCQERelErrCovarRC error matrix.
  // Divide these tweaking dials with the sqrt(diagonal error matrix elements) so that
  // they have a `canonical' definition (similar with all other systematic tweaking dials 
  // used in GENIE: mean = 0, sigma = 1.
  //

  TMatrixD & rc_err_matrix = *fEff1RmuCCQERCCovarRelErr;
  fEff1RmuCCQERC0TwkDial = twkv(0) / TMath::Sqrt(rc_err_matrix(0,0)); 
  fEff1RmuCCQERC1TwkDial = twkv(1) / TMath::Sqrt(rc_err_matrix(1,1));
  fEff1RmuCCQERC2TwkDial = twkv(2) / TMath::Sqrt(rc_err_matrix(2,2));
}
//________________________________________________________________________________________
void SuperkMC::TweakEfficiency1RmuCCnonQE(double twkdial) 
{ 
  fEff1RmuCCnonQETwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkMC::TweakEfficiency1RmuNC(double twkdial) 
{ 
  fEff1RmuNCTwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkMC::TweakEfficiency1RmuNuECC(double twkdial) 
{ 
  fEff1RmuNuECCTwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkMC::PrintSystematicVariations(void) const
{
  cout << "Exp/SKMC >> Twk 1Rmu/1Re energy scales = "     << fEnergyScaleTwkDial      << endl;
  cout << "Exp/SKMC >> Twk 1Rmu CCQE RC0 efficiency = "   << fEff1RmuCCQERC0TwkDial   << endl;
  cout << "Exp/SKMC >> Twk 1Rmu CCQE RC1 efficiency = "   << fEff1RmuCCQERC1TwkDial   << endl;
  cout << "Exp/SKMC >> Twk 1Rmu CCQE RC2 efficiency = "   << fEff1RmuCCQERC2TwkDial   << endl;
  cout << "Exp/SKMC >> Twk 1Rmu CCQE Oth efficiency = "   << fEff1RmuCCQEOthTwkDial   << endl;
  cout << "Exp/SKMC >> Twk 1Rmu CCnonQE efficiency = "    << fEff1RmuCCnonQETwkDial   << endl;
  cout << "Exp/SKMC >> Twk 1Rmu NC efficiency = "         << fEff1RmuNCTwkDial        << endl;
  cout << "Exp/SKMC >> Twk 1Rmu nu_e CC efficiency = "    << fEff1RmuNuECCTwkDial     << endl;
}
//________________________________________________________________________________________
double SuperkMC::EfficiencyWeight(
    SuperkSample_t s, SuperkMode_t m, double Ereco) const
{
  //
  // 1-ring mu-like sample
  //
  if(s == kSKsFCFV1RingMuLike) {
    switch(m) {

       // numu CCQE, numubar CCQE:
       // (1) Apply the CCQE other efficiency error. No energy dependence.
       // (2) Apply the CCQE ring-counting efficiency error. Ereco-dependence.
       //
       case (kSKmNuMuCCQE)    : 
       case (kSKmNuMuBarCCQE) : 
       {
          // other efficiency err
          double relerr_oth  = fEff1RmuCCQEOthRelErr;
          double twkdial_oth = fEff1RmuCCQEOthTwkDial;
          double weight_oth  = (1. + twkdial_oth*relerr_oth);
          // ring-counting efficiency err
          TMatrixD & rc_err_matrix = *fEff1RmuCCQERCCovarRelErr;
          double weight_rc0 = 1 + fEff1RmuCCQERC0TwkDial * TMath::Sqrt(rc_err_matrix(0,0));
          double weight_rc1 = 1 + fEff1RmuCCQERC1TwkDial * TMath::Sqrt(rc_err_matrix(1,1));
          double weight_rc2 = 1 + fEff1RmuCCQERC2TwkDial * TMath::Sqrt(rc_err_matrix(2,2));
          double weight_rc  = 1.;
          if      (Ereco < fEff1RmuCCQERCCovarE1RelErr) { weight_rc = weight_rc0; }
          else if (Ereco > fEff1RmuCCQERCCovarE2RelErr) { weight_rc = weight_rc2; }
          else                            { weight_rc = weight_rc1; }
          // combined weight
          double weight     = weight_oth * weight_rc;
          return weight;
          break;
       }

       //
       // numu CC1pi, numu CCoth, numubar CC1pi, numubar CCoth:
       // Apply the CCnonQE efficiency error. No energy dependence.
       //
       case (kSKmNuMuCC1pi)    :
       case (kSKmNuMuCCoth_mu)    : 
       case (kSKmNuMuBarCC1pi) : 
       case (kSKmNuMuBarCCoth_mu) : 
       {
          double relerr  = fEff1RmuCCnonQERelErr;
          double twkdial = fEff1RmuCCnonQETwkDial;
          double weight  = (1. + twkdial*relerr);
          return weight;
          break;
       }

       //
       // numu/nutau NC, numubar/nutaubar NC, nue NC:
       // Apply the NC efficiency error. No energy dependence.
       //
       case (kSKmNuMuNC)        : 
       case (kSKmNuMuBarNC) : 
       case (kSKmNuENC)     : 
       {
          double relerr  = fEff1RmuNCRelErr;
          double twkdial = fEff1RmuNCTwkDial;
          double weight  = (1. + twkdial*relerr);
          return weight;
          break;
       }

       //
       // nue CCQE, nue CC1pi, nue CCoth:
       // Apply the nue CC efficiency error. No energy dependence.
       //
       case (kSKmNuECCQE)     : 
       case (kSKmNuECC1pi)    : 
       case (kSKmNuECCoth_mu)    : 
       case (kSKmSigNuECCQE)  : 
       case (kSKmSigNuECC1pi) : 
       case (kSKmSigNuECCoth) : 
       {
          double relerr  = fEff1RmuNuECCRelErr;
          double twkdial = fEff1RmuNuECCTwkDial;
          double weight  = (1. + twkdial*relerr);
          return weight;
          break;
       }

       default : 
       {
          return 0.;                 
          break;
       }
    }//mode
  }//1-ring mu-like sample
  //
  // 1-ring e-like sample
  //
  else if ( s == kSKsFCFV1RingELike ) {
    switch(m) {
       case ( kSKmNuMuCCQE        ) : 
       case ( kSKmNuMuCC1pi       ) : 
       case ( kSKmNuMuCCcoh       ) : 
       case ( kSKmNuMuCCoth_e     ) : 
       case ( kSKmNuMuNC1pi0      ) : 
       case ( kSKmNuMuNCcoh       ) : 
       case ( kSKmNuMuNCoth_e     ) : 
       case ( kSKmNuMuBarCCQE     ) : 
       case ( kSKmNuMuBarCC1pi    ) : 
       case ( kSKmNuMuBarCCcoh    ) : 
       case ( kSKmNuMuBarCCoth_e  ) : 
       case ( kSKmNuMuBarNC1pi0   ) : 
       case ( kSKmNuMuBarNCcoh    ) : 
       case ( kSKmNuMuBarNCoth_e  ) : 
       case ( kSKmNuECCQE         ) : 
       case ( kSKmNuECC1pi        ) : 
       case ( kSKmNuECCcoh        ) : 
       case ( kSKmNuECCoth_e      ) : 
       case ( kSKmNuENC1pi0       ) : 
       case ( kSKmNuENCcoh        ) : 
       case ( kSKmNuENCoth_e      ) : 
/////////	 {
/////////	   double relerr  = fEff1ReBkgRelErr;
/////////	   double twkdial = fEff1ReBkgTwkDial;
/////////	   double weight  = (1. + twkdial*relerr);
/////////	   return weight;
/////////	   break;
/////////	 }
       case ( kSKmSigNuECCQE      ) : 
       case ( kSKmSigNuECC1pi     ) : 
       case ( kSKmSigNuECCoth     ) : 
/////////	 {
/////////	   double relerr  = fEff1ReSigRelErr;
/////////	   double twkdial = fEff1ReSigTwkDial;
/////////	   double weight  = (1. + twkdial*relerr);
/////////	   return weight;
/////////	   break;
/////////	 }
       default:
	 return 0.;
    }//mode
  }//1-ring e-like sample

  return 0.;
}
//________________________________________________________________________________________
double SuperkMC::EnergyScaleRelErr(SuperkSample_t s, SuperkMode_t m) const
{
// energy scale relative error

  if(s == kSKsFCFV1RingMuLike) 
  {
     return fEnergyScale1RmuRelErr;
  }
  else 
  if (s == kSKsFCFV1RingELike) {
    if( util::IsSignalNuE(m)) {
/////////      return fEnergyScale1ReSigRelErr;
    }
    else {
/////////      return fEnergyScale1ReBkgRelErr;
    }
  }
  return 0; 
}
//________________________________________________________________________________________
double SuperkMC::EnergyScale(SuperkSample_t s, SuperkMode_t m) const
{
// energy scale

  double twkdial = fEnergyScaleTwkDial;
  double err     = this->EnergyScaleRelErr(s,m);
  double scale   = 1. + twkdial * err;
  return scale;
}
//________________________________________________________________________________________
PhysHypothesis & SuperkMC::Hypothesis(void) 
{ 
  return fPhysHypothesis; 
}
//________________________________________________________________________________________
void SuperkMC::Draw(TVirtualPad * c, Option_t * option) const
{
  if(!c) return;

  c->cd();

  TString opt = option;
  util::SetDefaultStyle();

  if(opt.Contains("EtrueSpectra/NuMu")) { 
      fMCTreeNuMu1Rmu->Draw("Enu");
      fMCTreeNuMu1Rmu->Draw("Enu","Mode==1","same");
      fMCTreeNuMu1Rmu->Draw("Enu","Mode==11||Mode==12||Mode==13","same");
      c->Update();
  }
/*
  if(opt.Contains("1RingMuLike/ErecoEtrueTemplates")) { 
    double emin=0;
    double emax=6;
    const char * draw_opt = "COLZ";
    c->Divide(2,6);
    c->cd(1);
    c->GetPad(1)->SetLogz();
    c->GetPad(1)->SetBottomMargin(0.32);
    c->GetPad(1)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->Draw(draw_opt);
    c->cd(2);
    c->GetPad(2)->SetLogz();
    c->GetPad(2)->SetBottomMargin(0.32);
    c->GetPad(2)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->Draw(draw_opt);
    c->cd(3);
    c->GetPad(3)->SetLogz();
    c->GetPad(3)->SetBottomMargin(0.32);
    c->GetPad(3)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth_mu]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth_mu]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth_mu]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth_mu]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth_mu]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth_mu]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth_mu]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth_mu]->Draw(draw_opt);
    c->cd(4);
    c->GetPad(4)->SetLogz();
    c->GetPad(4)->SetBottomMargin(0.32);
    c->GetPad(4)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuNC]->Draw(draw_opt);
    c->cd(5);
    c->GetPad(5)->SetLogz();
    c->GetPad(5)->SetBottomMargin(0.32);
    c->GetPad(5)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->Draw(draw_opt);
    c->cd(6);
    c->GetPad(6)->SetLogz();
    c->GetPad(6)->SetBottomMargin(0.32);
    c->GetPad(6)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->Draw(draw_opt);
    c->cd(7);
    c->GetPad(7)->SetLogz();
    c->GetPad(7)->SetBottomMargin(0.32);
    c->GetPad(7)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth_mu]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth_mu]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth_mu]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth_mu]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth_mu]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth_mu]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth_mu]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth_mu]->Draw(draw_opt);
    c->cd(8);
    c->GetPad(8)->SetLogz();
    c->GetPad(8)->SetBottomMargin(0.32);
    c->GetPad(8)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->Draw(draw_opt);
    c->cd(9);
    c->GetPad(9)->SetLogz();
    c->GetPad(9)->SetBottomMargin(0.32);
    c->GetPad(9)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCQE]->Draw(draw_opt);
    c->cd(10);
    c->GetPad(10)->SetLogz();
    c->GetPad(10)->SetBottomMargin(0.32);
    c->GetPad(10)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->Draw(draw_opt);
    c->cd(11);
    c->GetPad(11)->SetLogz();
    c->GetPad(11)->SetBottomMargin(0.32);
    c->GetPad(11)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth_mu]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth_mu]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth_mu]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth_mu]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth_mu]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth_mu]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth_mu]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth_mu]->Draw(draw_opt);
    c->Update();
    c->cd(12);
    c->GetPad(12)->SetLogz();
    c->GetPad(12)->SetBottomMargin(0.32);
    c->GetPad(12)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuENC]->Draw(draw_opt);
    c->Update();
  }
*/
  c->Update();
}
//________________________________________________________________________________________
bool SuperkMC::SkipTemplate(SuperkMode_t m) 
{
// ignore nue_signal templates for all but the 3-flavour oscillation hypothesis

  bool mode_is_nuesig = util::IsSignalNuE(m);

  PhysHypothesis_t hyp = fPhysHypothesis.Hypothesis();
  bool hyp_is_not_3flavour_osc = (hyp != kHypOsc3F);

  bool skip = (mode_is_nuesig && hyp_is_not_3flavour_osc);

#ifdef __show_skipped_templates__
  if(skip) {
    cout << "Exp/SKMC >> Skipping template: " << util::AsString(m) << endl;
  }
#endif

  return skip;
}
//________________________________________________________________________________________


