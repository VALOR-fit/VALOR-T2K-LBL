//________________________________________________________________________________________
/*
 Changes since v1.0.0:

 @ 16/03/2011 - CA
   Created from the old SuperkExpectation::SuperkInputs nested class. Removed the numu10a 
   namespace.
 @ 28/03/2011 - CA
   Remove CurrErecoSpectrum() as it doesn't make sense using the combined spectrum
   for all the different samples considered. The default was to return the spectrum for the
   `1-R mu-like' sample but not typing this explicitly could lead to confusion. From now only
   CurrErecoSpectrum(s) and CurrErecoSpectrum(s,m) will be used.
   Similar changes were make to SuperkExpectation and SuperkObservation classes.
 @ 06/06/2011 - CA
   Added PhysHypothesis class.
   Added util::UseEventToFillInTemplate(SuperkMode_t m, int pdg, int reaction_mode).
   Added new (Ereco-dependent) `CCQE ring-counting eff. error' using the Choleslky method.
   Renamed old `CCQE efficiency err' to `CCQE other eff. err' (it now includes all
   sources of errors except ring-counting).
   NumOfEvents_Ereco(), NumOfEvents_Etrue() were modified accordingly to take the new
   Ereco-dependent error into account.
 @ 05/07/2011 - CA
   Modified TweakRingCountingEfficiencyCCQE() to allow specifying whether the forced tweak
   dials are the correlated or de-correlated ones. The de-correlated ones are used if the
   CCQE ring-counting nuisance params are included in the fit. The correlated ones are
   used for certain validation programs. Used the updated PhysHypothesis class.
 @ 20/07/2011 - CA
   Code was re-structured. The systematic variations and the physics hypothesis 
   (eg. oscillations) are no longer applied in a  single-step in all MC templates as 
   the Ereco spectrum gets updated.  All work is now done in UpdateTemplates(). 
   The UpdateErecoSpetra() and UpdateEtrueSpectra are simple loops over MC templates. 
   3-flavour oscillations were plugged-in. At this point matter effects are neglected 
   but they would be included shortly.
 @ 26/07/2011 - CA
   Took newly0-added `nue signal' templates into account. They should be neglected under 
   all but the  3-flavour oscillation hypothesis. Added bool SkipTemplate(SuperkMode_t m) 
   to decide whether a template should be neglected. 
   Added some comments and removed commented-out code.
 @ 29/07/11
   Include matter effects in 3-flavour osc probabilities.
 @ 06/10/11 - CA
   Added `PrintSystematicVariations()'. Tweaked output messages.
 @ 20-25/10/11 - CA,TD
   Substantial structural modifications so that the code can predict both the 1Rmu and 
   1Re distributions. The list of reaction modes kept track of for each sample may be
   different and it is defined in Conventions.h. Implemented 1Re detector systematics.

*/
//________________________________________________________________________________________

#define __show_init_mesg_info__
#define __show_init_mesg_debug__
//#define __show_hypothesis__ 
//#define __show_hypothesis_probabilities__ 
//#define __show_num_of_events_mesg_info__ 
//#define __show_num_of_events_mesg_debug__ 
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
      fMCErecoSpectrum [is][im] = 0;
    }
  }
  
  fEnergyScaleRelErr    = 0.;
  fOthEffCCQERelErr     = 0.;
  fEffCCQERelErrCovarRC = 0;
  fEreco1RCCCQE         = 0.;
  fEreco2RCCCQE         = 0.;      
  fEffCCnonQERelErr     = 0.;
  fEffNCRelErr          = 0.;
  fEffNuECCRelErr       = 0.;
  fEff1ReTSigRelErr     = 0.;
  fEff1ReTBkgRelErr     = 0.;
  

  fEnergyScaleTwkDial   = 0.;
  fOthEffCCQETwkDial    = 0.;
  fRC0EffCCQETwkDial    = 0.;
  fRC1EffCCQETwkDial    = 0.;
  fRC2EffCCQETwkDial    = 0.;
  fEffCCnonQETwkDial    = 0.;
  fEffNCTwkDial         = 0.;
  fEffNuECCTwkDial      = 0.;
  fEff1ReTSigTwkDial    = 0.;
  fEff1ReTBkgTwkDial    = 0.;

  fEffCCQERelErrCovarCholeskyL  = 0;
  
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
  fGenInputs      = 0;
  
  f1RmNuMuEt      = 0;
  f1RmNuMuBarEt   = 0;
  f1RmNuEEt       = 0;
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
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     for(int im = 0; im < knm; im++) {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       if( fMCTemplNom     [s][m] ) delete fMCTemplNom     [s][m];
       if( fMCTemplSys     [s][m] ) delete fMCTemplSys     [s][m];
       if( fMCTemplSysHyp  [s][m] ) delete fMCTemplSysHyp  [s][m];
       if( fMCErecoSpectrum[s][m] ) delete fMCErecoSpectrum[s][m];
     }
  }

  if (f1RmNuMuEt   ) delete f1RmNuMuEt;
  if (f1RmNuMuBarEt) delete f1RmNuMuBarEt;
  if (f1RmNuEEt    ) delete f1RmNuEEt;
}
//________________________________________________________________________________________
bool SuperkMC::Read(TXMLEngine* xml, XMLNodePointer_t root_node)
{
// Read SuperK inputs to disappearance analysis
//
#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> Reading SuperK inputs from input XML file." << endl;
#endif

  //
  // Get data POT
  // 
  fPOT = util::XML2Dbl(xml,root_node,"superk/observation/pot");
#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> " << "POT_{data} = " << fPOT << endl;
#endif
  
  //
  // Get MC trees and MC POT
  //

  fMCPOTNuMu    = util::XML2Dbl  (xml,root_node,"superk/mc/numu/pot");
  fMCPOTNuMuBar = util::XML2Dbl  (xml,root_node,"superk/mc/numubar/pot");
  fMCPOTNuE     = util::XML2Dbl  (xml,root_node,"superk/mc/nuebeam/pot");
  fMCPOTNuESig  = util::XML2Dbl  (xml,root_node,"superk/mc/nuesignal/pot");

  fMCTreeNuMu1Rmu = util::XML2TTree(xml,root_node,"superk/mc/numu/one_mu_ring");
  assert(fMCTreeNuMu1Rmu);
  fMCTreeNuMuBar1Rmu = util::XML2TTree(xml,root_node,"superk/mc/numubar/one_mu_ring");
  assert(fMCTreeNuMuBar1Rmu);
  fMCTreeNuE1Rmu = util::XML2TTree(xml,root_node,"superk/mc/nuebeam/one_mu_ring");
  assert(fMCTreeNuE1Rmu);
  fMCTreeNuESig1Rmu = util::XML2TTree(xml,root_node,"superk/mc/nuesignal/one_mu_ring");
  assert(fMCTreeNuESig1Rmu);

  fMCTreeNuMu1Re = util::XML2TTree(xml,root_node,"superk/mc/numu/one_e_ring");
  assert(fMCTreeNuMu1Re);
  fMCTreeNuMuBar1Re = util::XML2TTree(xml,root_node,"superk/mc/numubar/one_e_ring");
  assert(fMCTreeNuMuBar1Re);
  fMCTreeNuE1Re = util::XML2TTree(xml,root_node,"superk/mc/nuebeam/one_e_ring");
  assert(fMCTreeNuE1Re);
  fMCTreeNuESig1Re = util::XML2TTree(xml,root_node,"superk/mc/nuesignal/one_e_ring");
  assert(fMCTreeNuESig1Re);

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
  
  // create templates 
#ifdef __show_init_mesg_info__
  cout << "Exp/SKMC >> Filling-in SuperK MC templates... " << endl;
#endif
  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     for(int im = 0; im < knm; im++) {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       this->CreateTemplate(s, m);
     }
  }

  //
  // Systematics
  //

  // Get relative error on energy scale
  fEnergyScaleRelErr = util::XML2Dbl(
     xml,root_node,"superk/systematics/energy_scale_relative_error");
  // Get combined relative error on efficiencies
  fOthEffCCQERelErr =
    util::XML2Dbl(xml,root_node,
      "superk/systematics/efficiency_relative_error/one_mu_ring/ccqe/other_error");
  fEffCCQERelErrCovarRC = 
    util::XML2TMatrixD(xml,root_node,
      "superk/systematics/efficiency_relative_error/one_mu_ring/ccqe/rc_error");
  fEreco1RCCCQE =
    util::XML2Dbl(xml,root_node,
      "superk/systematics/efficiency_relative_error/one_mu_ring/ccqe/rc_ereco_high_bin1");
  fEreco2RCCCQE =
    util::XML2Dbl(xml,root_node,
      "superk/systematics/efficiency_relative_error/one_mu_ring/ccqe/rc_ereco_high_bin2");
  fEffCCnonQERelErr   =
    util::XML2Dbl(xml,root_node,
      "superk/systematics/efficiency_relative_error/one_mu_ring/ccnonqe");
  fEffNCRelErr        =
    util::XML2Dbl(xml,root_node,
      "superk/systematics/efficiency_relative_error/one_mu_ring/nc");
  fEffNuECCRelErr     =
    util::XML2Dbl(xml,root_node, 
      "superk/systematics/efficiency_relative_error/one_mu_ring/nuecc");
  fEff1ReTSigRelErr      =
    util::XML2Dbl(xml,root_node, 
      "superk/systematics/efficiency_relative_error/one_e_ring/sig");
  fEff1ReTBkgRelErr      =
    util::XML2Dbl(xml,root_node, 
      "superk/systematics/efficiency_relative_error/one_e_ring/bkg");

  if(!fEffCCQERelErrCovarRC) {
    cerr << "Exp/SK >> No input SuperK CCQE ring-counting error covariance matrix!" << endl;
    exit(1);
  }
      
#ifdef __show_init_mesg_info__      
  cout << "Exp/SKMC >> " << "Error {EnergyScale} : " << fEnergyScaleRelErr  << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency/1RingMu/CCQE/Other} : " << fOthEffCCQERelErr << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency/1RingMu/CCQE/RingCounting} : " << endl;
  fEffCCQERelErrCovarRC->Print();
  cout << "Exp/SKMC >> " << "Ring-counting error matrix is given for the following 3 Ereco bins : ";
  cout << "(0, " << fEreco1RCCCQE << " GeV), "
       << "[" << fEreco1RCCCQE << " GeV, " << fEreco2RCCCQE << " GeV], "
       << "(" << fEreco2RCCCQE << " GeV, oo)" << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency/1RingMu/CCnonQE} : " << fEffCCnonQERelErr  << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency/1RingMu/NC} : "      << fEffNCRelErr       << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency/1RingMu/NuECC} : "   << fEffNuECCRelErr    << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency/1RingE/Sig} : "      << fEff1ReTSigRelErr    << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency/1RingE/Bkg} : "      << fEff1ReTBkgRelErr    << endl;
#endif
  
  // initialise physics hypothesis
  fPhysHypothesis.Init(kHypUndefined);

  //
  // Create histograms
  //

  // Etrue histograms needed for beam norm err calculation
  f1RmNuMuEt    = new TH1D("f1RmNuMuEt",    "Etrue for 1-ring mu-like numu",    200, 0, 10); /// binning expected by Kubo's code
  f1RmNuMuBarEt = new TH1D("f1RmNuMuBarEt", "Etrue for 1-ring mu-like numubar", 200, 0, 10);
  f1RmNuEEt     = new TH1D("f1RmNuMuEEt",   "Etrue for 1-ring mu-like nue",     200, 0, 10);
  f1RmNuMuEt     -> SetDirectory(0);
  f1RmNuMuBarEt  -> SetDirectory(0);
  f1RmNuEEt      -> SetDirectory(0);
    
  // Ereco histograms filled-in by SK spectrum prediction code
  for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;
     fMCErecoSpectrumTot[s] = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
     fMCErecoSpectrumTot[s] -> SetDirectory(0);
     for(int im=0; im < knm; im++) {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;
       fMCErecoSpectrum[s][m] = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
       fMCErecoSpectrum[s][m] -> SetDirectory(0);
     }
  }

#ifdef __show_init_mesg_info__      
  cout << "Exp/SKMC >> Done reading SuperK inputs!" << endl;
#endif
  
  return true;
}
//________________________________________________________________________________________
void SuperkMC::CreateTemplate(SuperkSample_t s, SuperkMode_t m)
{
#ifdef __show_init_mesg_debug__      
  cout << "Exp/SKMC/CreateTemplate >> Creating MC template for: "
       << util::AsString(s) << " / " << util::AsString(m) << endl;
#endif
    
  int is = (int)s;
  int im = (int)m;
    
  fMCTemplNom[is][im] = new TH2D(
     Form("nsk_%d_%d",is,im), util::AsString(m).c_str(),
       kNvErecoBins, kNvErecoBinEdges, kNEtrueBins, kEtrueMin, kEtrueMax);
  fMCTemplNom[is][im]->SetDirectory(0);
  fMCTemplNom[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
  fMCTemplNom[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");

  fMCTemplSys[is][im] = new TH2D(
     Form("nsksyst_%d_%d",is,im), util::AsString(m).c_str(),
       kNvErecoBins, kNvErecoBinEdges, kNEtrueBins, kEtrueMin, kEtrueMax);
  fMCTemplSys[is][im]->SetDirectory(0);
  fMCTemplSys[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
  fMCTemplSys[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");

  fMCTemplSysHyp[is][im] = new TH2D(
     Form("nskcurr_%d_%d",is,im), util::AsString(m).c_str(),
       kNvErecoBins, kNvErecoBinEdges, kNEtrueBins, kEtrueMin, kEtrueMax);
  fMCTemplSysHyp[is][im]->SetDirectory(0);
  fMCTemplSysHyp[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
  fMCTemplSysHyp[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
      
  TTree * mctree = this->MCTree(s,m);
  double  mcpot  = this->MCPOT (s,m);

  if(!mctree) {
     cerr << "Exp/SKMC/CreateTemplate >> Couldn't find event tree for template: " 
          << util::AsString(s) << "/" << util::AsString(m) << endl;
     exit(1);
  }

  double wght = this->AnalyzedSamplePOT() / mcpot;

#ifdef __show_init_mesg_debug__        
  cout << "Exp/SKMC/CreateTemplate >>  - Weight (data POT / MC POT) = " 
       << wght << endl;
#endif
  if(wght <= 0) {
     cerr << "Exp/SKMC/CreateTemplate >> Weight (data POT / MC POT) for " 
          << util::AsString(s) << "/" << util::AsString(m)
          << " is " << wght << endl;
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
     
   //cout << "Exp/SKMC >> Ereco = " << fBrEreco << " GeV, Etrue = " 
   //     << fBrEtrue << " GeV, wght = " << wght << endl;

     fMCTemplNom    [is][im] -> Fill (fBrEreco, fBrEtrue, wght);
     fMCTemplSys    [is][im] -> Fill (fBrEreco, fBrEtrue, wght);
     fMCTemplSysHyp [is][im] -> Fill (fBrEreco, fBrEtrue, wght);

  }//iev

#ifdef __show_init_mesg_debug__        
  cout << "Exp/SKMC/CreateTemplate >>  - MC template was filled with " 
       << fMCTemplNom[is][im]->GetEntries() << " entries." << endl;
#endif
}
//________________________________________________________________________________________
//void SuperkMC::Update(void)
//{
//  this->ApplySystematicsToTemplates();
//  this->ApplyPhysHypothesisToTemplates();
//  this->UpdateErecoSpectra();
//  this->UpdateEtrueSpectra();
//}
//________________________________________________________________________________________
void SuperkMC::UpdateTemplates(void)
{
  this->ApplySystematicsToTemplates();
  this->ApplyPhysHypothesisToTemplates();
}
//________________________________________________________________________________________
void SuperkMC::ApplySystematicsToTemplates(void) 
{
// Apply systematic variation to MC templates

  //
  // Apply a) flux 10a->10d weights, b) SuperK efficiency tweaks (if any),
  // c) SuperK energy scale tweaks (if any), d) flux shape and norm tweaks (if any)
  // and e) cross section tweaks (if any).
  // All systematics are multiplicative factors applied to the number of events
  // and they are either functions of Ereco, or Etrue or have no dependence at all.
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
       TH2D * mctempl = fMCTemplSys[s][m];
       if(!mctempl) continue; 
       // loop over Ereco
       for(int ir = 1; ir <= mctempl->GetXaxis()->GetNbins(); ir++) {
         double Ereco = mctempl->GetXaxis()->GetBinCenter(ir);
         // loop over Etrue
         for(int it = 1; it <= mctempl->GetYaxis()->GetNbins(); it++) {
           double Etrue = mctempl->GetYaxis()->GetBinCenter(it);
           // calculate number of events in this (Ereco,Etrue) bin including all
           // systematic variations
           double nev = this->NumOfEvents(Ereco,Etrue,s,m);
           // update MC template
           mctempl->SetBinContent(ir,it,nev);
         }//it
       }//ir
    }//im
  }//is

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
       TH2D * mctempl = fMCTemplSys[s][m];
       if(!mctempl) continue; 
       // clone
       TH2D * mctempl_no_ereco_scale = new TH2D(*mctempl);
       // get energy scale
       double escale = this->EnergyScale();

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
    double sinsq2theta12 = fPhysHypothesis.CurrValue(kPhypSinsq2Theta12);
    double sinsq2theta13 = fPhysHypothesis.CurrValue(kPhypSinsq2Theta13);
    double sinsq2theta23 = fPhysHypothesis.CurrValue(kPhypSinsq2Theta23);
    double dm12sq        = fPhysHypothesis.CurrValue(kPhypDm12sq);
    double dm23sq        = fPhysHypothesis.CurrValue(kPhypDm23sq);
    double delta_cp      = fPhysHypothesis.CurrValue(kPhypDeltaCP);
#ifdef __show_hypothesis__
    cout << "**  Exp/SKMC >> sin^2(2theta12) = " << sinsq2theta12
         <<   ", sin^2(2theta13) = " << sinsq2theta13
         <<   ", sin^2(2theta23) = " << sinsq2theta23
         <<   ", delta_cp = " << delta_cp
         <<   ", Dm^2_{12} = "   << dm12sq << " eV^2"
         <<   ", |Dm^2_{23}| = " << dm23sq << " eV^2"
         << endl;
#endif
    fNuOsc3Flavour.SetOscParameters(sinsq2theta12,sinsq2theta13,sinsq2theta23,dm12sq,dm23sq,delta_cp);

    // loop over samples
    for(int is = 0; is < kns; is++) {
     SuperkSample_t s = kSKSamples[is];
     if(s == kSKsUndefined) continue;

     assert(fMCTemplSys[s][kSKmNuMuCCQE]);
     TAxis * ereco_axis = fMCTemplSys[s][kSKmNuMuCCQE]->GetXaxis(); 
     TAxis * etrue_axis = fMCTemplSys[s][kSKmNuMuCCQE]->GetYaxis(); 

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
              double nev = fMCTemplSys[s][m] -> GetBinContent(ir,it);
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
              fMCTemplSysHyp[s][m] -> SetBinContent(ir,it, nev_osc);
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
    double sinsq2theta23 = fPhysHypothesis.CurrValue(kPhypSinsq2Theta23);
    double dm23sq        = fPhysHypothesis.CurrValue(kPhypDm23sq);
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
         TH2D * mctempl_src = fMCTemplSys[s][m];
         if(!mctempl_src) continue; 
         TH2D * mctempl_tgt = fMCTemplSysHyp[s][m];
         if(!mctempl_tgt) continue; 
         // loop over Etrue
         for(int it = 1; it <= mctempl_src->GetYaxis()->GetNbins(); it++) {
           double Etrue = mctempl_src->GetYaxis()->GetBinCenter(it);
           double psurv = util::NuMuSurvProbOsc2D(sinsq2theta23, dm23sq, Etrue, m);
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
         TH2D * mctempl_src = fMCTemplSys[s][m];
         if(!mctempl_src) continue; 
         TH2D * mctempl_tgt = fMCTemplSysHyp[s][m];
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
         TH2D * mctempl_src = fMCTemplSys[s][m];
         if(!mctempl_src) continue; 
         TH2D * mctempl_tgt = fMCTemplSysHyp[s][m];
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
double SuperkMC::NumOfEvents(
    double Ereco, double Etrue, SuperkSample_t s, SuperkMode_t m) const
{
// - Apply flux 10a->10d weights.
// - Apply SuperK efficiency tweaks (if any). 
// - Apply SuperK energy scale tweaks (if any).
// - Apply flux tweaks (if any).
// - Apply cross section tweaks (if any).

  TH2D * nsk2d = fMCTemplNom[s][m];
  assert(nsk2d);

  int reco_bin = nsk2d->GetXaxis()->FindBin(Ereco);
  int true_bin = nsk2d->GetYaxis()->FindBin(Etrue);

#ifdef __show_num_of_events_mesg_info__ 
  cout << "Exp/SKMC/NumOfEvents >> "
       << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m) 
       << ", Ereco = " << Ereco << " GeV (bin = " << reco_bin << ")" 
       << ", Etrue = " << Etrue << " GeV (bin = " << true_bin << ")" << endl;
#endif

  // Nominal number of events
  double nev = nsk2d->GetBinContent(reco_bin, true_bin);
        
  //
  // Apply flux weights & systematics = f(Etrue)
  //
  double flux_10d_weight      = 1.; // 10a->10d weight
  double numu_flux_shape_wght = 1.; // systematic variation (flux shape)
  if(fNuFlux) {
     if (util::IsNuMu(m)) {
        flux_10d_weight      = fNuFlux->SuperkNuMuFlux10aTo10dWght(Etrue);
        numu_flux_shape_wght = fNuFlux->SuperkNuMuFluxShapeWeight (Etrue);
     }
     else
     if (util::IsNuMuBar(m)) {
        flux_10d_weight = fNuFlux->SuperkNuMuBarFlux10aTo10dWght(Etrue);
     }
     else
     if (util::IsBeamNuE(m)) {
        flux_10d_weight = fNuFlux->SuperkNuEFlux10aTo10dWght(Etrue);
     } 
     else
     if (util::IsSignalNuE(m)) {
        // signal nue component comes from the numu flux:
        flux_10d_weight      = fNuFlux->SuperkNuMuFlux10aTo10dWght(Etrue);
        numu_flux_shape_wght = fNuFlux->SuperkNuMuFluxShapeWeight (Etrue);
     } 
     else {
        cerr << "Exp/SKMC/NumOfEvents >> Can't handle mode = " << util::AsString(m) << endl;
        exit(1);
     }
  }//flux
#ifdef __show_num_of_events_mesg_debug__
  if(flux_10d_weight != 1) {
    cout << "Exp/SKMC/NumOfEvents >> Applied flux tuning weight = " << flux_10d_weight << endl;
  }
  if(numu_flux_shape_wght != 1) {
    cout << "Exp/SKMC/NumOfEvents >> Applied flux-shape-systematic weight = " << numu_flux_shape_wght << endl;
  }
#endif

  //
  // Apply cross section systematics = f(Etrue)
  //
  double xsec_ccqe_wght  = 1.;
  double xsec_ratio_wght = 1.;
  if(fGenInputs) {
     xsec_ccqe_wght  = fGenInputs->SuperkCCQEXSecWght (Etrue, s, m);
     xsec_ratio_wght = fGenInputs->SuperkXSecRatioWght(Etrue, s, m);
  }
#ifdef __show_num_of_events_mesg_debug__
  if(xsec_ccqe_wght != 1.) {
    cout << "Exp/SKMC/NumOfEvents >> Applied CCQE-xsec-systematic weight = " << xsec_ccqe_wght << endl;
  }
  if(xsec_ratio_wght != 1.) {
    cout << "Exp/SKMC/NumOfEvents >> Applied mode/CCQE-xsec-ratio-systematic weight = " << xsec_ratio_wght << endl;
  }
#endif

  //
  // Apply cross section systematics = f(Ereco)
  //
  double nucl_effect_wght = 1;
  if(fGenInputs) {
     nucl_effect_wght = fGenInputs->SuperkNuclEffectWght(Ereco, s, m);
  }
#ifdef __show_num_of_events_mesg_debug__ 
  if(nucl_effect_wght != 1.) {
    cout << "Exp/SKMC/NumOfEvents >> Applied nucl.-effect-systematic weight = " << nucl_effect_wght << endl;
  }
#endif

  //
  // Plug-in "efficiency" uncertainty due to fiducial volume, reduction chain, OD cut,
  // ring-counting, PID, momentum cut and decay e- cut.
  // The ring-counting error for CCQE has Ereco-dependence.
  //
  double effwght = this->EfficiencyWeight(s,m,Ereco);
#ifdef __show_num_of_events_mesg_debug__
  if(effwght != 1.) {
    cout << "Exp/SKMC/NumOfEvents >> Applied efficiency-systematic weight = " << effwght << endl;
  }
#endif

  //
  // Calculate updated number of events
  //
  double nevtwk = 
     ( nev * 
       flux_10d_weight * 
       numu_flux_shape_wght * 
       xsec_ccqe_wght *
       xsec_ratio_wght *
       nucl_effect_wght *
       effwght
     );

#ifdef __show_num_of_events_mesg_info__ 
  cout << "Exp/SKMC/NumOfEvents >> " << nev << " events (nominal) --> " 
       << nevtwk << " events (current / tweaked)" << endl;
#endif

  return nevtwk;
}
//________________________________________________________________________________________
void SuperkMC::UpdateErecoSpectra(void)
{
  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = (SuperkSample_t) is;

    TH1D * reco_spectrum_tot =  fMCErecoSpectrumTot[is];
    if(!reco_spectrum_tot) continue; 

    // loop over Ereco
    for(int ir = 1; ir <= reco_spectrum_tot->GetNbinsX(); ir++) {
       //double Ereco = mctempl->GetXaxis()->GetBinCenter(ir);
       double nev_tot = 0;

       // loop over modes
       for(int im = 0; im < knm; im++) { 
          SuperkMode_t m = (SuperkMode_t) im;

          // consider template?
          bool skip = this->SkipTemplate(m);
          if(skip) continue;

          TH2D * mctempl = fMCTemplSysHyp[s][m];
          if(!mctempl) continue; 
          TH1D * reco_spectrum =  fMCErecoSpectrum[is][im];
          if(!reco_spectrum) continue;

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
void SuperkMC::UpdateEtrueSpectra(void)
{
//cout << "Exp/SKMC/UpdateEtrueSpectra >> Updating Etrue spectra for beam norm err" << endl;

  TH1D * spectrum = 0;
  TH2D * mctempl1 = 0;
  TH2D * mctempl2 = 0;
  TH2D * mctempl3 = 0;
  TH2D * mctempl4 = 0;

//
//
// NOTE  REQUIRE SAME ETRUE BINNING 
//
//

  // update 1-ring mu-like numu
  spectrum = f1RmNuMuEt;
  mctempl1 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuMuCCQE];
  mctempl2 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi];
  mctempl3 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuMuCCoth_mu];
  mctempl4 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuMuNC];
  for(int it = 1; it <= spectrum->GetNbinsX(); it++) {
//  double Etrue = spectrum->GetBinCenter(ibin);
    double nev = 0;
    for(int ir = 1; ir <= mctempl1->GetXaxis()->GetNbins(); ir++) {
       nev += ( mctempl1->GetBinContent(ir,it)  +
                mctempl2->GetBinContent(ir,it) +
                mctempl3->GetBinContent(ir,it) +
                mctempl4->GetBinContent(ir,it)
              );
    }//ir
#ifdef __show_spectra_update_mesg_debug__
    cout << "Exp/SKMC/UpdateEtrueSpectra >> "
         << "1-ring mu / numu : true bin = " << it << ", nev = " << nev << endl;
#endif
    spectrum->SetBinContent(it,nev);
  }//true bins

  // update 1-ring mu-like numubar
  spectrum = f1RmNuMuBarEt;
  mctempl1 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE];
  mctempl2 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi];
  mctempl3 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth_mu];
  mctempl4 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuMuBarNC];
  for(int it = 1; it <= spectrum->GetNbinsX(); it++) {
//  double Etrue = spectrum->GetBinCenter(ibin);
    double nev = 0;
    for(int ir = 1; ir <= mctempl1->GetXaxis()->GetNbins(); ir++) {
       nev += ( mctempl1->GetBinContent(ir,it)  +
                mctempl2->GetBinContent(ir,it) +
                mctempl3->GetBinContent(ir,it) +
                mctempl4->GetBinContent(ir,it)
              );
    }//ir
#ifdef __show_spectra_update_mesg_debug__
    cout << "Exp/SKMC/UpdateEtrueSpectra >> "
         << "1-ring mu / numubar : true bin = " << it << ", nev = " << nev << endl;
#endif
    spectrum->SetBinContent(it,nev);
  }//true bins


  // update 1-ring mu-like nue
  spectrum = f1RmNuEEt;
  mctempl1 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuECCQE];
  mctempl2 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuECC1pi];
  mctempl3 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuECCoth_mu];
  mctempl4 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuENC];
  for(int it = 1; it <= spectrum->GetNbinsX(); it++) {
//  double Etrue = spectrum->GetBinCenter(ibin);
    double nev = 0;
    for(int ir = 1; ir <= mctempl1->GetXaxis()->GetNbins(); ir++) {
       nev += ( mctempl1->GetBinContent(ir,it)  +
                mctempl2->GetBinContent(ir,it) +
                mctempl3->GetBinContent(ir,it) +
                mctempl4->GetBinContent(ir,it)
              );
    }//ir
#ifdef __show_spectra_update_mesg_debug__
    cout << "Exp/SKMC/UpdateEtrueSpectra() >> "
         << "1-ring mu / nuebkg : true bin = " << it << ", nev = " << nev << endl;
#endif
    spectrum->SetBinContent(it,nev);
  }//true bins

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
void SuperkMC::TweakOtherEfficiencyCCQE(double twkdial) 
{ 
  fOthEffCCQETwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkMC::TweakRingCountingEfficiencyCCQE(
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
  if(!fEffCCQERelErrCovarCholeskyL) 
  {
#ifdef __show_rc_cholesky_mesg_info__
    cout << "Exp/SKMC/TweakRingCountingEffCCQE >> "
         << "Run Cholesky's decomposition method on the CCQE RC error matrix." << endl;
#endif

    TMatrixD ccqe_rc_err_matrix(*fEffCCQERelErrCovarRC);
#ifdef __show_rc_cholesky_mesg_debug__
    cout << "Exp/SKMC/TweakRingCountingEffCCQE >> RC error matrix:" << endl;
    ccqe_rc_err_matrix.Print();
#endif

    TMatrixD L = util::CholeskyDecomposition(ccqe_rc_err_matrix);
    fEffCCQERelErrCovarCholeskyL = new TMatrixD(L);
#ifdef __show_rc_cholesky_mesg_info__
    cout << "Exp/SKMC/TweakRingCountingEffCCQE >> "
         << "Cholesky decomposition of RC error matrix:" << endl;
    fEffCCQERelErrCovarCholeskyL->Print();
#endif
  }


  //
  // Check whether the values of the tweaking dials are forced. 
  // If the forced tweaking dials are supposed to be the correlated ones then we don't
  // need to use the Cholesky decomposition of the  ring-counting error matrix.
  //
  if(force_input_twkdials && corellated) {
    fRC0EffCCQETwkDial = twkdial0;
    fRC1EffCCQETwkDial = twkdial1;
    fRC2EffCCQETwkDial = twkdial2;
    return;
  }

  TMatrixD & L = *fEffCCQERelErrCovarCholeskyL;
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

  TMatrixD & rc_err_matrix = *fEffCCQERelErrCovarRC;
  fRC0EffCCQETwkDial = twkv(0) / TMath::Sqrt(rc_err_matrix(0,0)); 
  fRC1EffCCQETwkDial = twkv(1) / TMath::Sqrt(rc_err_matrix(1,1));
  fRC2EffCCQETwkDial = twkv(2) / TMath::Sqrt(rc_err_matrix(2,2));
}
//________________________________________________________________________________________
void SuperkMC::TweakEfficiencyCCnonQE(double twkdial) 
{ 
  fEffCCnonQETwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkMC::TweakEfficiencyNC(double twkdial) 
{ 
  fEffNCTwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkMC::TweakEfficiencyNuECC(double twkdial) 
{ 
  fEffNuECCTwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkMC::TweakEfficiency1ReTSig(double twkdial) 
{ 
  fEff1ReTSigTwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkMC::TweakEfficiency1ReTBkg(double twkdial) 
{ 
  fEff1ReTBkgTwkDial = twkdial; 
}
//________________________________________________________________________________________
void SuperkMC::PrintSystematicVariations(void) const
{
  cout << "Exp/SKMC >> Twk energy scale = "        << fEnergyScaleTwkDial << endl;
  cout << "Exp/SKMC >> Twk CCQE RC0 efficiency = " << fRC0EffCCQETwkDial  << endl;
  cout << "Exp/SKMC >> Twk CCQE RC1 efficiency = " << fRC1EffCCQETwkDial  << endl;
  cout << "Exp/SKMC >> Twk CCQE RC2 efficiency = " << fRC2EffCCQETwkDial  << endl;
  cout << "Exp/SKMC >> Twk CCQE Oth efficiency = " << fOthEffCCQETwkDial  << endl;
  cout << "Exp/SKMC >> Twk CCnonQE efficiency = "  << fEffCCnonQETwkDial  << endl;
  cout << "Exp/SKMC >> Twk NC efficiency = "       << fEffNCTwkDial       << endl;
  cout << "Exp/SKMC >> Twk nu_e CC efficiency = "  << fEffNuECCTwkDial    << endl;
  cout << "Exp/SKMC >> Twk 1ReT signal efficiency = "  << fEff1ReTSigTwkDial    << endl;
  cout << "Exp/SKMC >> Twk 1ReT background efficiency = "  << fEff1ReTBkgTwkDial    << endl;

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
          double relerr_oth  = fOthEffCCQERelErr;
          double twkdial_oth = fOthEffCCQETwkDial;
          double weight_oth  = (1. + twkdial_oth*relerr_oth);
          // ring-counting efficiency err
          TMatrixD & rc_err_matrix = *fEffCCQERelErrCovarRC;
          double weight_rc0 = 1 + fRC0EffCCQETwkDial * TMath::Sqrt(rc_err_matrix(0,0));
          double weight_rc1 = 1 + fRC1EffCCQETwkDial * TMath::Sqrt(rc_err_matrix(1,1));
          double weight_rc2 = 1 + fRC2EffCCQETwkDial * TMath::Sqrt(rc_err_matrix(2,2));
          double weight_rc  = 1.;
          if      (Ereco < fEreco1RCCCQE) { weight_rc = weight_rc0; }
          else if (Ereco > fEreco2RCCCQE) { weight_rc = weight_rc2; }
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
          double relerr  = fEffCCnonQERelErr;
          double twkdial = fEffCCnonQETwkDial;
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
          double relerr  = fEffNCRelErr;
          double twkdial = fEffNCTwkDial;
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
          double relerr  = fEffNuECCRelErr;
          double twkdial = fEffNuECCTwkDial;
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
	 {
	   double relerr  = fEff1ReTBkgRelErr;
	   double twkdial = fEff1ReTBkgTwkDial;
	   double weight  = (1. + twkdial*relerr);
	   return weight;
	   break;
	 }
       case ( kSKmSigNuECCQE      ) : 
       case ( kSKmSigNuECC1pi     ) : 
       case ( kSKmSigNuECCoth     ) : 
	 {
	   double relerr  = fEff1ReTSigRelErr;
	   double twkdial = fEff1ReTSigTwkDial;
	   double weight  = (1. + twkdial*relerr);
	   return weight;
	   break;
	 }
       default:
	 return 0.;
    }//mode
  }//1-ring e-like sample

  return 0.;
}
//________________________________________________________________________________________
double SuperkMC::EnergyScaleRelErr(void) const
{
// energy scale relative error

  return fEnergyScaleRelErr;
}
//________________________________________________________________________________________
double SuperkMC::EnergyScale(void) const
{
// energy scale

  double twkdial = fEnergyScaleTwkDial;
  double err     = this->EnergyScaleRelErr();
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

  c->Update();
}
//________________________________________________________________________________________
bool SuperkMC::SkipTemplate(t2kana::conv::SuperkMode_t m) 
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


