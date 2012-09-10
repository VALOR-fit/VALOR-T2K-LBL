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
//#define __show_init_mesg_debug2__
#define __show_update_mesg_info__
//#define __show_hypothesis__ 
//#define __show_hypothesis_probabilities__ 
//#define __show_syst_var_wght_mesg_debug__ 
//#define __show_spectra_update_mesg_debug__ 
//#define __show_skipped_templates__ 
//#define __show_isanevent_mesg_debug__
//#define __save_templates_debug__   // saves to files, and cout's number of events, at various stages of applying tunes to the templates:
// 1 number of events
// 2 number of events, weighted by pot & mc normalisation
// 3 number of events, weighted by pot & mc normalisation, flux tuning
// 4 number of events, weighted by pot & mc normalisation, flux tuning, banff tuning
// 5 number of events, weighted by pot & mc normalisation, flux tuning, banff tuning, multiplicative systematic variations
// 6 number of events, weighted by pot & mc normalisation, flux tuning, banff tuning, multiplicative systematic variations, nonmultiplicative systematic variations
// 7 number of events, weighted by pot & mc normalisation, flux tuning, banff tuning, multiplicative systematic variations, nonmultiplicative systematic variations, physics hypothesis

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <TTree.h>
#include <TH2D.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TMath.h>

#include "SuperkMC.h"
#include "SuperkDetector.h"
#include "FSI.h"
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
    fMCEtrueSpectrumTot[is] = 0;
    for(int im=0; im<knm; im++) {
      fMCTemplOri      [is][im] = 0;
      fMCTemplNom      [is][im] = 0;
      fMCTemplTun      [is][im] = 0;
      fMCTemplSys      [is][im] = 0;
      fMCTemplFin      [is][im] = 0;
      fMCErecoSpectrum [is][im] = 0;
    }
  }
    
  fMCTreeNuMu1Rmu      = 0;
  fMCTreeNuMuBar1Rmu   = 0;
  fMCTreeNuE1Rmu       = 0;
  fMCTreeNuEBar1Rmu    = 0;
  fMCTreeSigNuE1Rmu    = 0;
  fMCTreeNuMu1Re       = 0;
  fMCTreeNuMuBar1Re    = 0;
  fMCTreeNuE1Re        = 0;
  fMCTreeNuEBar1Re     = 0;
  fMCTreeSigNuE1Re     = 0;

  fMCPOTNuMu           = 0;
  fMCPOTNuMuBar        = 0;
  fMCPOTNuE            = 0;
  fMCPOTNuEBar         = 0;
  fMCPOTSigNuE         = 0;

  fBrEreco             = 0;
  fBrEtrue             = 0;
  fBrGenCode           = 0;
         
  fNuFlux              = 0;
  fNuXSec              = 0;  
  fFSI                 = 0;  
  fSKDet               = 0;
}
//________________________________________________________________________________________
void SuperkMC::CleanUp(void)
{
  if(fMCTreeNuMu1Rmu    ) delete fMCTreeNuMu1Rmu;
  if(fMCTreeNuMuBar1Rmu ) delete fMCTreeNuMuBar1Rmu;
  if(fMCTreeNuE1Rmu     ) delete fMCTreeNuE1Rmu;
  if(fMCTreeNuEBar1Rmu  ) delete fMCTreeNuEBar1Rmu;
  if(fMCTreeSigNuE1Rmu  ) delete fMCTreeSigNuE1Rmu;
  if(fMCTreeNuMu1Re     ) delete fMCTreeNuMu1Re;
  if(fMCTreeNuMuBar1Re  ) delete fMCTreeNuMuBar1Re;
  if(fMCTreeNuE1Re      ) delete fMCTreeNuE1Re;
  if(fMCTreeNuEBar1Re   ) delete fMCTreeNuEBar1Re;
  if(fMCTreeSigNuE1Re   ) delete fMCTreeSigNuE1Re;

  for(int is = 0; is < kns; is++) {
     for(int im = 0; im < knm; im++) {
       if( fMCTemplOri     [is][im] ) delete fMCTemplOri     [is][im];
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
  cout << "Exp/SKMC >> Reading SuperK MC inputs from input XML file." << endl;
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
  fMCPOTHistNuMu    = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numu/pot");
  fMCPOTHistNuMuBar = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/numubar/pot");
  fMCPOTHistNuE     = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nue/pot");
  fMCPOTHistNuEBar  = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/nuebar/pot");
  fMCPOTHistSigNuE  = util::XML2TH1D  (xml,root_node,"mc/superk/nominal/signal_nue/pot");

  fMCPOTNuMu    = fMCPOTHistNuMu->GetBinContent(2)    / fMCPOTHistNuMu->GetBinContent(1);
  fMCPOTNuMuBar = fMCPOTHistNuMuBar->GetBinContent(2) / fMCPOTHistNuMuBar->GetBinContent(1);
  fMCPOTNuE     = fMCPOTHistNuE->GetBinContent(2)     / fMCPOTHistNuE->GetBinContent(1);
  fMCPOTNuEBar  = fMCPOTHistNuEBar->GetBinContent(2)  / fMCPOTHistNuEBar->GetBinContent(1);
  fMCPOTSigNuE  = fMCPOTHistSigNuE->GetBinContent(2)  / fMCPOTHistSigNuE->GetBinContent(1);

  fMCTreeNuMu1Rmu    = util::XML2TTree (xml,root_node,"mc/superk/nominal/numu/one_mu_ring");
  fMCTreeNuMuBar1Rmu = util::XML2TTree (xml,root_node,"mc/superk/nominal/numubar/one_mu_ring");
  fMCTreeNuE1Rmu     = util::XML2TTree (xml,root_node,"mc/superk/nominal/nue/one_mu_ring");
  fMCTreeNuEBar1Rmu  = util::XML2TTree (xml,root_node,"mc/superk/nominal/nuebar/one_mu_ring");
  fMCTreeSigNuE1Rmu  = util::XML2TTree (xml,root_node,"mc/superk/nominal/signal_nue/one_mu_ring");
  fMCTreeNuMu1Re     = util::XML2TTree (xml,root_node,"mc/superk/nominal/numu/one_e_ring");
  fMCTreeNuMuBar1Re  = util::XML2TTree (xml,root_node,"mc/superk/nominal/numubar/one_e_ring");
  fMCTreeNuE1Re      = util::XML2TTree (xml,root_node,"mc/superk/nominal/nue/one_e_ring");
  fMCTreeNuEBar1Re   = util::XML2TTree (xml,root_node,"mc/superk/nominal/nuebar/one_e_ring");
  fMCTreeSigNuE1Re   = util::XML2TTree (xml,root_node,"mc/superk/nominal/signal_nue/one_e_ring");

  assert (fMCTreeNuMu1Rmu   );
  assert (fMCTreeNuMuBar1Rmu);
  assert (fMCTreeNuE1Rmu    );
  assert (fMCTreeNuEBar1Rmu );
  assert (fMCTreeSigNuE1Rmu );
  assert (fMCTreeNuMu1Re    );
  assert (fMCTreeNuMuBar1Re );
  assert (fMCTreeNuE1Re     );
  assert (fMCTreeNuEBar1Re  );
  assert (fMCTreeSigNuE1Re  );

#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> " << "POT {MC; numu}    = " << fMCPOTNuMu    << endl;
  cout << "Exp/SKMC >> " << "POT {MC; numubar} = " << fMCPOTNuMuBar << endl;
  cout << "Exp/SKMC >> " << "POT {MC; nue_bkg} = " << fMCPOTNuE     << endl;
  cout << "Exp/SKMC >> " << "POT {MC; nue_sig} = " << fMCPOTSigNuE  << endl;

  cout << "Exp/SKMC >> " << "Nev {MC; numu;     1Rmu} = " << fMCTreeNuMu1Rmu    -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; numubar;  1Rmu} = " << fMCTreeNuMuBar1Rmu -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; nue;      1Rmu} = " << fMCTreeNuE1Rmu     -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; nuebar;   1Rmu} = " << fMCTreeNuEBar1Rmu  -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; nue[sig]; 1Rmu} = " << fMCTreeSigNuE1Rmu  -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; numu;     1Re } = " << fMCTreeNuMu1Re     -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; numubar;  1Re } = " << fMCTreeNuMuBar1Re  -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; nue;      1Re } = " << fMCTreeNuE1Re      -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; nuebar;   1Re } = " << fMCTreeNuEBar1Re   -> GetEntries()  << endl;
  cout << "Exp/SKMC >> " << "Nev {MC; nue[sig]; 1Re } = " << fMCTreeSigNuE1Re   -> GetEntries()  << endl;
#endif
    
  // initialise physics hypothesis
  // need to initialise as 3 flavour oscillation, otherwise ExtTuning (NA61) and FinTuning (BANFF) is not applied to nuesig modes
  fPhysHypothesis.Init(kHypOsc3F);

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

#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtOri, "nevents");
  this->NevDebug(kSKMCtNom, "mcnorm");
#endif

  // need to read in the xsec splines, before applying any BANFF tweaks
  fNuXSec -> BuildXsBinningMCTemplates(fMCTemplOri);
  if( !fNuXSec -> ReadSplines(xml,root_node) ) return false;

  // apply flux tuning (based on external data) to templates
  //this->ApplyFluxExtTuningToTemplates();

#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtNom, "mcnorm_fluxtuned");
#endif

  // apply flux and cross-section tuning (based on nd280 data fit by BANFF) to templates
  this->ApplyBANFFTuningToTemplates();

  // fill the etrue spectra
  this->FillEtrueSpectra();

#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtTun, "mcnorm_fluxtuned_banfftuned");
#endif

#ifdef __show_init_mesg_info__      
  cout << "Exp/SKMC >> Done reading SuperK MC inputs!" << endl;
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
     fMCEtrueSpectrumTot[is] = new TH1D("","", kNEtrueBins, kEtrueBinEdges);
     fMCEtrueSpectrumTot[is] -> SetDirectory(0);
     fMCEtrueSpectrumTot[is]->SetLineWidth(2);
     fMCEtrueSpectrumTot[is]->SetLineColor(kBlack);
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
       fMCTemplOri[is][im] = new TH2D(
          Form("nsk_ori_%d_%d",is,im), util::AsString(m).c_str(),
	       n_ereco_bins, ereco_bin_edges, kNEtrueBins, kEtrueBinEdges);
       fMCTemplNom[is][im] = new TH2D(
          Form("nsk_nom_%d_%d",is,im), util::AsString(m).c_str(),
               n_ereco_bins, ereco_bin_edges, kNEtrueBins, kEtrueBinEdges);
       fMCTemplTun[is][im] = new TH2D(
          Form("nsk_tun_%d_%d",is,im), util::AsString(m).c_str(),
               n_ereco_bins, ereco_bin_edges, kNEtrueBins, kEtrueBinEdges);
       fMCTemplSys[is][im] = new TH2D(
          Form("nsk_sys_%d_%d",is,im), util::AsString(m).c_str(),
               n_ereco_bins, ereco_bin_edges, kNEtrueBins, kEtrueBinEdges);
       fMCTemplFin[is][im] = new TH2D(
          Form("nsk_fin_%d_%d",is,im), util::AsString(m).c_str(),
               n_ereco_bins, ereco_bin_edges, kNEtrueBins, kEtrueBinEdges);
       fMCTemplOri[is][im]->SetDirectory(0);
       fMCTemplOri[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
       fMCTemplOri[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
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

  TBranch * branchEreco   = mctree->GetBranch("erec");
  TBranch * branchEtrue   = mctree->GetBranch("enu");
  TBranch * branchGenCode = mctree->GetBranch("mode");
  TBranch * branchWght    = mctree->GetBranch("wgt");

  branchEreco   -> SetAddress(&fBrEreco);
  branchEtrue   -> SetAddress(&fBrEtrue);
  branchGenCode -> SetAddress(&fBrGenCode);
  branchWght    -> SetAddress(&fBrWght);

  mctree->SetBranchStatus("*",0);
  mctree->SetBranchStatus("erec",1);
  mctree->SetBranchStatus("enu",1);
  mctree->SetBranchStatus("mode",1);
  mctree->SetBranchStatus("wgt",1);
  
  for(int iev = 0; iev < mctree->GetEntries(); iev++) 
  {
     mctree->GetEntry(iev);

     // use current entry for mode m template?
     bool use_it = util::UseEventToFillInTemplate(m,util::SuperkMode2SKMCPdg(m),fBrGenCode);
     if(!use_it) continue;
  
     this -> Template(s,m,kSKMCtOri) -> Fill (fBrEreco/1000., fBrEtrue);
     this -> Template(s,m,kSKMCtNom) -> Fill (fBrEreco/1000., fBrEtrue, fBrWght * pot_wght);
     this -> Template(s,m,kSKMCtTun) -> Fill (fBrEreco/1000., fBrEtrue, fBrWght * pot_wght);
     this -> Template(s,m,kSKMCtSys) -> Fill (fBrEreco/1000., fBrEtrue, fBrWght * pot_wght);
     this -> Template(s,m,kSKMCtFin) -> Fill (fBrEreco/1000., fBrEtrue, fBrWght * pot_wght);

#ifdef __show_init_mesg_debug__
     if(fBrEtrue <= 0 || fBrEtrue >= this->Template(s,m,kSKMCtOri)->GetYaxis()->GetBinUpEdge(this->Template(s,m,kSKMCtOri)->GetNbinsY()))
       cout << "Exp/SKMC/CreateTemplate >> Etrue is out of template range " << AsString(s) << " " << AsString(m) << " " << fBrEtrue << endl;
#ifdef __show_init_mesg_debug2__
     if(fBrEreco/1000 <= 0)
       cout << "Exp/SKMC/CreateTemplate >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << fBrEreco/1000 << endl;
#endif
     if(fBrEreco/1000 >= this->Template(s,m,kSKMCtOri)->GetXaxis()->GetBinUpEdge(this->Template(s,m,kSKMCtOri)->GetNbinsX()))
       cout << "Exp/SKMC/CreateTemplate >> Ereco is out of template range " << AsString(s) << " " << AsString(m) << " " << fBrEreco/1000 << endl;
#endif

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
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  return fMCErecoSpectrumTot[sidx]; 
}
//________________________________________________________________________________________
TH1D* SuperkMC::CurrErecoSpectrum(SuperkSample_t s, SuperkMode_t m) const 
{ 
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }            
  return fMCErecoSpectrum[sidx][midx]; 
}
//________________________________________________________________________________________
TH2D * SuperkMC::Template(SuperkSample_t s, SuperkMode_t m, SKMCTemplateType_t type) const
{
  int sidx = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  if(sidx == -1) {
    return 0;
  }
  int midx = util::ArrIdx::Instance()->SuperkModeIdx(s,m);
  if(midx == -1) {
    return 0;
  }            

  TH2D * h = 0;

  switch(type) {
   case (kSKMCtOri) : 
     h = fMCTemplOri[sidx][midx]; 
     break;
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

#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> Applying flux tuning based on external data (NA61,...)." << endl;
#endif

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
           double wght = fNuFlux->SuperkFluxExtTuningWght(s,m,Etrue);
 	   if(wght < kMinAllowedWght) {
	     cout << "set weight " << wght << " to kMinAllowedWght " << kMinAllowedWght << endl;
	     wght = kMinAllowedWght;
	   }
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

#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> Applying BANFF flux and cross-section tuning." << endl;
#endif

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
         double Ereco = mctempl_src->GetXaxis()->GetBinCenter(ir);

	 //
	 // There are 600 50MeV Etrue energy bins in the templates, but
	 // the flux/xsec systematics only have up to ~10 Etrue bins (same as flux binning), so
	 // we can call NuXSec::FinTuningWght() and NuFlux::SuperkFluxFinTuningWght()  
	 // a factor of 10 less often by looping over the kXsEtrueBinEdges.
	 // Note that the Ereco binning is identical between XSec and SuperkMC templates
	 //
	 // Loop over variable flux Etrue bins
	 for(int itv = 1; itv <= util::NumVariableMCEtrueBins(s,m); itv ++) {

	   // Get the min/max energy of the current variable flux Etrue bin
	   double EtrueVarMin = util::VariableMCEtrueBinEdges(s,m)[itv-1] + kOneMeV;
	   double EtrueVarMax = util::VariableMCEtrueBinEdges(s,m)[itv] - kOneMeV;

           // Plug-in BANFF flux systematics
           double flux_wght = fNuFlux->SuperkFluxFinTuningWght(s,m,EtrueVarMin);
	   
           // Plug-in BANFF & other xsec systematics 
           double xsec_wght = fNuXSec->FinTuningWght(s,m,Ereco,EtrueVarMin); 

           // Combined weight
           double wght = flux_wght * xsec_wght;
	   if(wght < kMinAllowedWght) {
	     cout << "set weight " << wght << " to kMinAllowedWght " << kMinAllowedWght << endl;
	     wght = kMinAllowedWght;
	   }

	   // Get the 50MeV Etrue bin number range
	   int firstbin = mctempl_src->GetYaxis()->FindBin(EtrueVarMin);
	   int lastbin = mctempl_src->GetYaxis()->FindBin(EtrueVarMax);

	   // Loop over 50MeV Etrue bins
	   for(int it = firstbin; it <= lastbin; it++) {
	     // get number of events in current bin of source MC template
	     double nev = mctempl_src->GetBinContent(ir,it);
	     nev = TMath::Max(0.,nev);
	     double nev_tuned = nev * wght;
	     mctempl_tgt->SetBinContent(ir,it,nev_tuned);
	     // initialize "sys" & "fin" templates to "tun" values
	     this->Template(s,m,kSKMCtSys)->SetBinContent(ir,it,nev_tuned);
	     this->Template(s,m,kSKMCtFin)->SetBinContent(ir,it,nev_tuned);
	   }//it
	 }//itv
       }//ir
    }//im
  }//is

  fFSI    -> BuildFracMCTemplates(fMCTemplTun);
}
//________________________________________________________________________________________
void SuperkMC::UpdateTemplates(void)
{
#ifdef __show_update_mesg_info__
  cout << "Exp/SKMC >> Updating MC templates...." << endl;
#endif

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
//
#ifdef __show_update_mesg_info__
  cout << "Exp/SKMC >> Applying multiplicative systematics" << endl;
#endif

  // Loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    // Loop over sample modes
    for(int im = 0; im < knm; im++) { 
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       // Consider template?
       bool skip = this->SkipTemplate(m);
       if(skip) continue;
       // Get corresponding templates
       TH2D * mctempl_src = this -> Template(s,m,kSKMCtTun);
       if(!mctempl_src) continue;
       TH2D * mctempl_tgt = this -> Template(s,m,kSKMCtSys);
       if(!mctempl_tgt) continue;

       // Loop over Ereco
       for(int ir = 1; ir <= mctempl_src->GetXaxis()->GetNbins(); ir++) {
         double Ereco = mctempl_src->GetXaxis()->GetBinCenter(ir);
	 //
	 // Get combined weight due to all multiplicative systematic variations.
	 // This method considers:
	 // - SuperK efficiencies variations
	 // - BANFF parameter variations
	 // - Cross-section variations
	 //
	 // Plug-in SK efficiency uncertainty due to fiducial volume, reduction chain, OD cut,
	 // ring-counting, PID, momentum cut and decay e- cut.
	 // The ring-counting error for CCQE has Ereco-dependence.
	 double skeff_wght = fSKDet->EfficiencyWeight(s,m,Ereco);

	 //
	 // There are 600 50MeV Etrue energy bins in the templates, but
	 // the flux/xsec systematics only have up to ~10 Etrue bins (same as flux binning), so
	 // we can call NuXSec::SystWght() and FSI::SystWght() and NuFlux::SuperkFluxSystWght()  
	 // a factor of 10 less often by looping over the kXsEtrueBinEdges.
	 // Note that the Ereco binning is identical between XSec and SuperkMC templates
	 //
	 // Loop over variable flux Etrue bins
	 for(int itv = 1; itv <= util::NumVariableMCEtrueBins(s,m); itv ++) {

	   // Get the min/max energy of the current variable flux Etrue bin
	   double EtrueVarMin = util::VariableMCEtrueBinEdges(s,m)[itv-1] + kOneMeV;
	   double EtrueVarMax = util::VariableMCEtrueBinEdges(s,m)[itv] - kOneMeV;

           // Plug-in BANFF flux systematics
           double flux_wght = fNuFlux->SuperkFluxSystWght(s,m,EtrueVarMin);
	   
           // Plug-in BANFF & other xsec systematics 
           double xsec_wght = fNuXSec->SystWght(s,m,Ereco,EtrueVarMin); 

           // Plug-in FSI weights
           double fsi_wght = fFSI->SystWght(s,m,Ereco,EtrueVarMin); 

           // Combined weight
           double wght = skeff_wght * flux_wght * xsec_wght * fsi_wght;
	   if(wght < kMinAllowedWght) {
	     cout << "set weight " << wght << " to kMinAllowedWght " << kMinAllowedWght << endl;
	     wght = kMinAllowedWght;
	   }

	   // Get the 50MeV Etrue bin number range
	   int firstbin = mctempl_src->GetYaxis()->FindBin(EtrueVarMin);
	   int lastbin = mctempl_src->GetYaxis()->FindBin(EtrueVarMax);

	   // Loop over 50MeV Etrue bins
	   for(int it = firstbin; it <= lastbin; it++) {
	     // get number of events in current bin of source MC template
	     double nev = mctempl_src->GetBinContent(ir,it);
	     nev = TMath::Max(0.,nev);

#ifdef __show_syst_var_wght_mesg_debug__ 
	     if( skeff_wght != 1. || 
		 flux_wght  != 1. || 
		 xsec_wght  != 1. )
	       {
		 // Get the energy of the current 50MeV Etrue bin
		 double Etrue = mctempl_src->GetYaxis()->GetBinCenter(it);
		 cout << "Exp/SKMC/SystVarWght >> "
		      << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m) 
		      << ", Ereco = " << Ereco << " GeV" 
		      << ", Etrue = " << Etrue << " GeV => Weights: " 
		      << "SK eff = "  << skeff_wght
		      << ", Flux = "  << flux_wght
		      << ", XSec = "  << xsec_wght
		      << ", FSI  = "  << fsi_wght
		      << " (total = " << wght << ")"
		      << endl;
	       }
#endif

	     // Update target MC template
	     mctempl_tgt->SetBinContent(ir,it,nev*wght);
	   }//it
	 }//itv
       }//ir
    }//im
  }//is

#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtSys, "mcnorm_fluxtuned_banfftuned_syst_multi");
#endif
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

#ifdef __show_update_mesg_info__
  cout << "Exp/SKMC >> Applying non-multiplicative systematics" << endl;
#endif

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
       double escale = fSKDet->EnergyScale(s,m);

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
#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtSys, "mcnorm_fluxtuned_banfftuned_syst_nonmulti");
#endif
}
//________________________________________________________________________________________
void SuperkMC::ApplyPhysHypothesisToTemplates(void)
{
// Apply physics hypothesis {oscillation (2-flavour / 3-flavour), decay or 
// decoherence} to MC templates
//

#ifdef __show_update_mesg_info__
  cout << "Exp/SKMC >> Applying physics hypothesis" << endl;
#endif

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
	// skip unnecessary calculations
	bool bin_is_not_empty = this->IsThereAnEvent(s, Etrue);
	if(!bin_is_not_empty) continue;

        double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuMu,    kPdgNuMu,    Etrue);
        double p_numubar_numubar = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuMuBar, kPdgNuMuBar, Etrue);
        double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuE,     kPdgNuE,     Etrue);
        double p_nuebar_nuebar   = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuEBar,  kPdgNuEBar,  Etrue);
      //double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuE,     kPdgNuMu,    Etrue);
        double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuMu,    kPdgNuE,     Etrue);
#ifdef __show_hypothesis_probabilities__
        cout << "Exp/SKMC >>  E = " << Etrue << " GeV" << endl;
        cout << "Exp/SKMC >>  * P (numu     -> numu   ) = " << p_numu_numu       << endl;
        cout << "Exp/SKMC >>  * P (numubar  -> numubar) = " << p_numubar_numubar << endl;
        cout << "Exp/SKMC >>  * P (nue      -> nue    ) = " << p_nue_nue         << endl;
        cout << "Exp/SKMC >>  * P (nuebar   -> nuebar ) = " << p_nuebar_nuebar   << endl;
      //cout << "Exp/SKMC >>  * P (nue      -> numu   ) = " << p_nue_numu        << endl;
        cout << "Exp/SKMC >>  * P (numu     -> nue    ) = " << p_numu_nue        << endl;
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
                if( util::IsNuE(m)       ) { posc = p_nue_nue;         }
                if( util::IsNuEBar(m)    ) { posc = p_nuebar_nuebar;   }
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
#ifdef __save_templates_debug__
  this->NevDebug(kSKMCtFin, "mcnorm_fluxtuned_banfftuned_syst_phys");
#endif
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
void SuperkMC::FillEtrueSpectra(void)
{
  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    TH1D * true_spectrum_tot = fMCEtrueSpectrumTot[is];

    // loop over Etrue
    for(int it = 1; it <= true_spectrum_tot->GetNbinsX(); it++) {
       double nev_tot = 0;

       // loop over modes
       for(int im = 0; im < knm; im++) { 
          SuperkMode_t m = kSKReactionModes[is][im];
          if(m == kSKmUndefined) continue;
          // consider template?
          bool skip = this->SkipTemplate(m);
          if(skip) continue;

          TH2D * mctempl = this -> Template(s,m,kSKMCtOri);
 
          // loop over Ereco
          double nev = 0;
          for(int ir = 1; ir <= mctempl->GetXaxis()->GetNbins(); ir++) {
             nev += mctempl->GetBinContent(ir,it);
          }//ir
          nev_tot += nev;
       }//im
       true_spectrum_tot->SetBinContent(it,nev_tot);
    }//it
  }//is
}
//________________________________________________________________________________________
bool SuperkMC::IsThereAnEvent(SuperkSample_t s, double Etrue) const
{
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  assert(is != -1);

  TH1D * true_spectrum_tot = fMCEtrueSpectrumTot[is];
  assert(true_spectrum_tot);

  int it = true_spectrum_tot->GetXaxis()->FindBin(Etrue);

  bool ev = true_spectrum_tot->GetBinContent(it) > kASmallNum;

#ifdef __show_isanevent_mesg_debug__
  if(ev)
    cout << "Exp/SuperkMC >> There are " << true_spectrum_tot->GetBinContent(it)
	 << " events for sample = " << util::AsString(s) 
	 << ", Etrue = " << Etrue
	 << endl;
#endif

  return ev;
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
     if ( util::IsNuE(m)       ) return fMCTreeNuE1Rmu;
     if ( util::IsNuEBar(m)    ) return fMCTreeNuEBar1Rmu;
     if ( util::IsSignalNuE(m) ) return fMCTreeSigNuE1Rmu;
  }//1Rmu

  if (s==kSKsFCFV1RingELike)
  {
     if ( util::IsNuMu(m)      ) return fMCTreeNuMu1Re;
     if ( util::IsNuMuBar(m)   ) return fMCTreeNuMuBar1Re;
     if ( util::IsNuE(m)       ) return fMCTreeNuE1Re;
     if ( util::IsNuEBar(m)    ) return fMCTreeNuEBar1Re;
     if ( util::IsSignalNuE(m) ) return fMCTreeSigNuE1Re;
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
     if ( util::IsNuE(m)       ) return fMCPOTNuE;
     if ( util::IsNuEBar(m)    ) return fMCPOTNuEBar;
     if ( util::IsSignalNuE(m) ) return fMCPOTSigNuE;
  }
  return 0;
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
    c->Divide(4,4);
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
    c->cd(13);
    c->GetPad(13)->SetLogz();
    c->GetPad(13)->SetBottomMargin(0.32);
    c->GetPad(13)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCQE]->Draw(draw_opt);
    c->cd(14);
    c->GetPad(14)->SetLogz();
    c->GetPad(14)->SetBottomMargin(0.32);
    c->GetPad(14)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCC1pi]->Draw(draw_opt);
    c->cd(15);
    c->GetPad(15)->SetLogz();
    c->GetPad(15)->SetBottomMargin(0.32);
    c->GetPad(15)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth_mu]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth_mu]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth_mu]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth_mu]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth_mu]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth_mu]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth_mu]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarCCoth_mu]->Draw(draw_opt);
    c->Update();
    c->cd(16);
    c->GetPad(16)->SetLogz();
    c->GetPad(16)->SetBottomMargin(0.32);
    c->GetPad(16)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuEBarNC]->Draw(draw_opt);
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
void SuperkMC::NevDebug(SKMCTemplateType_t t, string name)
{
  TFile out(Form("Template_%s.root", name.c_str()),"RECREATE");
  static const int erecobins = 70;
  static const double erecobinedges[erecobins+1] = {
      0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450, 
      0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950, 
      1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  1.300,  1.350,  1.400,  1.450, 
      1.500,  1.550,  1.600,  1.650,  1.700,  1.750,  1.800,  1.850,  1.900,  1.950, 
      2.000,  2.050,  2.100,  2.150,  2.200,  2.250,  2.300,  2.350,  2.400,  2.450, 
      2.500,  2.550,  2.600,  2.650,  2.700,  2.750,  2.800,  2.850,  2.900,  2.950, 
      3.000,  3.100,  3.200,  3.300,  3.400,  3.500,  4.000,  5.000,  7.000, 10.000,
      30.000
  };
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;
    double nsk = 0;
    for(int im = 0; im < knm; im++) {
      SuperkMode_t m = kSKReactionModes[is][im];
      if(m == kSKmUndefined) continue;
      TH2D * mctempl = this -> Template(s,m,t);

      //first write the 2d template, and cout the mode integral
      mctempl->Write();
      cout << name << " " << util::AsString(s) << " " << util::AsString(m) << " " << mctempl->Integral() << endl;
      nsk += mctempl->Integral();

      //fill a 1D ereco spectrum
      TH1D * hrec = new TH1D(Form("erec_int_%s",AsStringCompact_Japan(m).c_str()), Form("erec_int_%s",AsStringCompact_Japan(m).c_str()), erecobins, erecobinedges);
      // loop over Ereco
      for(int ir = 1; ir <= hrec->GetNbinsX(); ir++) {
	//find low and high of ereco for current bin
	double rmin = hrec->GetXaxis()->GetBinLowEdge(ir) + kOneMeV;
	double rmax = hrec->GetXaxis()->GetBinUpEdge(ir)  - kOneMeV;
	//integrate the template
	int ir_min = mctempl->GetXaxis()->FindBin(rmin);
	int ir_max = mctempl->GetXaxis()->FindBin(rmax);
	int it_bins = mctempl->GetNbinsY();
	double nev = mctempl->Integral(ir_min, ir_max, 1, it_bins);
	hrec->SetBinContent(ir,nev);
      }
      hrec->Write();
      delete hrec;

      //fill a 1D ereco spectrum (fine binning)
      TH1D * hrecf = new TH1D(Form("erec_int_%s_fine",AsStringCompact_Japan(m).c_str()), Form("erec_int_%s",AsStringCompact_Japan(m).c_str()), mctempl->GetNbinsX(), mctempl->GetXaxis()->GetXbins()->GetArray());
      // loop over Ereco
      for(int ir = 1; ir <= hrecf->GetNbinsX(); ir++) {
	//find low and high of ereco for current bin
	double rmin = hrecf->GetXaxis()->GetBinLowEdge(ir) + kOneMeV;
	double rmax = hrecf->GetXaxis()->GetBinUpEdge(ir)  - kOneMeV;
	//integrate the template
	int ir_min = mctempl->GetXaxis()->FindBin(rmin);
	int ir_max = mctempl->GetXaxis()->FindBin(rmax);
	int it_bins = mctempl->GetNbinsY();
	double nev = mctempl->Integral(ir_min, ir_max, 1, it_bins);
	hrecf->SetBinContent(ir,nev);
      }
      hrecf->Write();
      delete hrecf;

      //fill a 1D etrue spectrum
      TH1D * htru = new TH1D(Form("enu_int_%s",AsStringCompact_Japan(m).c_str()), Form("enu_int_%s",AsStringCompact_Japan(m).c_str()), util::NumVariableXsEtrueBins(s,m), util::VariableXsEtrueBinEdges(s,m));
      // loop over Etrue
      for(int it = 1; it <= htru->GetNbinsX(); it++) {
	//find low and high of etrue for current bin
	double tmin = htru->GetXaxis()->GetBinLowEdge(it) + kOneMeV;
	double tmax = htru->GetXaxis()->GetBinUpEdge(it)  - kOneMeV;
	//integrate the template
	int it_min = mctempl->GetYaxis()->FindBin(tmin);
	int it_max = mctempl->GetYaxis()->FindBin(tmax);
	int ir_bins = mctempl->GetNbinsX();
	double nev = mctempl->Integral(1, ir_bins, it_min, it_max);
	htru->SetBinContent(it,nev);
      }
      htru->Write();
      delete htru;
    }//mode

    cout << "Total: " << nsk << endl << endl;
  }//sample
  out.Close();
}

