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

*/
//________________________________________________________________________________________

//#define __debug_level_1__
//#define __debug_level_2__
//#define __debug_level_3__

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
      fMCTemplNom     [is][im] = 0;
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

  fEnergyScaleTwkDial   = 0.;
  fOthEffCCQETwkDial    = 0.;
  fRC0EffCCQETwkDial    = 0.;
  fRC1EffCCQETwkDial    = 0.;
  fRC2EffCCQETwkDial    = 0.;
  fEffCCnonQETwkDial    = 0.;
  fEffNCTwkDial         = 0.;
  fEffNuECCTwkDial      = 0.;

  fEffCCQERelErrCovarCholeskyL  = 0;
  
  fMC1Rmu_numu    = 0;
  fMC1Rmu_numubar = 0;
  fMC1Rmu_nue     = 0;
  fMC1Rmu_nuesig  = 0;

  fMCPOT_numu     = 0;
  fMCPOT_numubar  = 0;
  fMCPOT_nue      = 0;
  fMCPOT_nuesig   = 0;

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
  if(fMC1Rmu_numu   ) delete fMC1Rmu_numu;
  if(fMC1Rmu_numubar) delete fMC1Rmu_numubar;
  if(fMC1Rmu_nue    ) delete fMC1Rmu_nue;
  if(fMC1Rmu_nuesig ) delete fMC1Rmu_nuesig;

  for(int is=0; is<kns; is++) {
    if( fMCErecoSpectrumTot[is] ) delete fMCErecoSpectrumTot[is];
    for(int im=0; im<knm; im++) {
      if( fMCTemplNom     [is][im] ) delete fMCTemplNom     [is][im];
      if( fMCTemplSys     [is][im] ) delete fMCTemplSys     [is][im];
      if( fMCTemplSysHyp  [is][im] ) delete fMCTemplSysHyp  [is][im];
      if( fMCErecoSpectrum[is][im] ) delete fMCErecoSpectrum[is][im];
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
  cout << "Exp/SKMC >> Reading SuperK inputs from input XML file." << endl;

  //
  // Get data POT
  // 
  fPOT = util::XML2Dbl(xml,root_node,"superk/observation/pot");
  cout << "Exp/SKMC >> " << "N_{POT}|data : " << fPOT << endl;
  
  //
  // Get MC trees and MC POT
  //
  fMCPOT_numu = util::XML2Dbl(xml,root_node,"superk/mc/numu/pot");
  fMC1Rmu_numu = util::XML2TTree(xml,root_node,"superk/mc/numu/one_mu_ring");
  cout << "Exp/SKMC >> " << "N_{POT}|(numu MC) : " << fMCPOT_numu << endl;
  assert(fMC1Rmu_numu);
  fMCPOT_numubar = util::XML2Dbl(xml,root_node,"superk/mc/numubar/pot");
  fMC1Rmu_numubar = util::XML2TTree(xml,root_node,"superk/mc/numubar/one_mu_ring");
  cout << "Exp/SKMC >> " << "N_{POT}|(numubar MC) : " << fMCPOT_numubar << endl;
  assert(fMC1Rmu_numubar);
  fMCPOT_nue = util::XML2Dbl(xml,root_node,"superk/mc/nuebeam/pot");
  fMC1Rmu_nue = util::XML2TTree(xml,root_node,"superk/mc/nuebeam/one_mu_ring");
  cout << "Exp/SKMC >> " << "N_{POT}|(beam nue MC) : " << fMCPOT_nue << endl;
  assert(fMC1Rmu_nue);
  fMCPOT_nuesig = util::XML2Dbl(xml,root_node,"superk/mc/nuesignal/pot");
  fMC1Rmu_nuesig = util::XML2TTree(xml,root_node,"superk/mc/nuesignal/one_mu_ring");
  cout << "Exp/SKMC >> " << "N_{POT}|(signal nue MC) : " << fMCPOT_nuesig << endl;
  assert(fMC1Rmu_nuesig);
  
  // create templates 
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuMuCCQE      );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuMuCC1pi     );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuMuCCoth     );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNC            );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuMuBarCCQE   );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuMuBarCC1pi  );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuMuBarCCoth  );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuMuBarNC     );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuECCQE       );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuECC1pi      );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuECCoth      );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmNuENC         );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmSigNuECCQE    );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmSigNuECC1pi   );
  this->CreateTemplate (kSKsFCFV1RingMuLike, kSKmSigNuECCoth   );

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
      
      
  cout << "Exp/SKMC >> " << "Error {EnergyScale} : " << fEnergyScaleRelErr  << endl;

  cout << "Exp/SKMC >> " << "Error {Efficiency/1RingMu/CCQE/Other} : " << fOthEffCCQERelErr << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency/1RingMu/CCQE/RingCounting} : " << endl;
  if(!fEffCCQERelErrCovarRC) {
    cerr << "Exp/SK >> No input SuperK CCQE ring-counting error covariance matrix!" << endl;
    exit(1);
  }
  fEffCCQERelErrCovarRC->Print();
  cout << "Exp/SKMC >> " << "Ring-counting error matrix is given for the following 3 Ereco bins : " << endl;
  cout << "(0, " << fEreco1RCCCQE << " GeV), "
       << "[" << fEreco1RCCCQE << " GeV, " << fEreco2RCCCQE << " GeV], "
       << "(" << fEreco2RCCCQE << " GeV, oo)" << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency/1RingMu/CCnonQE} : " << fEffCCnonQERelErr  << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency/1RingMu/NC} : "      << fEffNCRelErr       << endl;
  cout << "Exp/SKMC >> " << "Error {Efficiency/1RingMu/NuECC} : "   << fEffNuECCRelErr    << endl;
  
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
  fMCErecoSpectrumTot[kSKsFCFV1RingMuLike]                  = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]       = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]      = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]      = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNC]             = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]    = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]   = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]   = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]      = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuECCQE]        = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuECC1pi]       = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuECCoth]       = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuENC]          = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmSigNuECCQE]     = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmSigNuECC1pi]    = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmSigNuECCoth]    = new TH1D("","", kNvErecoBins, kNvErecoBinEdges);
  fMCErecoSpectrumTot[kSKsFCFV1RingMuLike]                  -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]       -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]      -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]      -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNC]             -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]    -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]   -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]   -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]      -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuECCQE]        -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuECC1pi]       -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuECCoth]       -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmNuENC]          -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmSigNuECCQE]     -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmSigNuECC1pi]    -> SetDirectory(0);
  fMCErecoSpectrum[kSKsFCFV1RingMuLike][kSKmSigNuECCoth]    -> SetDirectory(0);

  cout << "Exp/SKMC >> Done reading SuperK inputs!" << endl;
  
  return true;
}
//________________________________________________________________________________________
void SuperkMC::CreateTemplate(SuperkSample_t s, SuperkMode_t m)
{
  cout << "Exp/SKMC >> Filling Ereco-Etrue event template for : "
       << util::AsString(s) << "/" << util::AsString(m) << endl;
    
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
  double  mcpot  = this->MCPOT(s,m);

  double wght = this->AnalyzedSamplePOT() / mcpot;
  
  cout << "Exp/SKMC >> Weight (data POT / MC POT) = " << wght << endl;

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

  cout << "Exp/SKMC >> Filled template with " 
       << fMCTemplNom[is][im]->GetEntries() << " entries." << endl;
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
    // loop over sample modes
    for(int im = 0; im < knm; im++) { 
       // get sample/mode id
       SuperkSample_t s = (SuperkSample_t) is;
       SuperkMode_t   m = (SuperkMode_t  ) im;
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
    // loop over sample modes
    for(int im = 0; im < knm; im++) { 
       // get sample/mode id
       SuperkSample_t s = (SuperkSample_t) is;
       SuperkMode_t   m = (SuperkMode_t  ) im;
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
#ifdef __debug_level_2__
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
     SuperkSample_t s = (SuperkSample_t) is;

     // loop over Etrue
     for(int it = 1; it <= fMCTemplSys[s][kSKmNuMuCCQE]->GetYaxis()->GetNbins(); it++) {

        // calculate oscillation probabilities (incl. matter effects) for the current Etrue value
        double Etrue = fMCTemplSys[s][kSKmNuMuCCQE]->GetYaxis()->GetBinCenter(it);
/*
        ////Due to numerous `Hermitian conjugate of matter mixing matrix may not be its inverse.'
        ////err mesgs, matter effects are not included. Nick tries to understand what is going on
        ////for delta_{cp} around pi.
        double p_numu_numu       = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuMu,    kPdgNuMu,    Etrue);
      //double p_nue_numu        = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuE,     kPdgNuMu,    Etrue);
        double p_numubar_numubar = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuMuBar, kPdgNuMuBar, Etrue);
        double p_nue_nue         = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuE,     kPdgNuE,     Etrue);
        double p_numu_nue        = fNuOsc3Flavour.CalcMatterOscProb(kPdgNuMu,    kPdgNuE,     Etrue);
*/
        double p_numu_numu       = fNuOsc3Flavour.CalcVacuumOscProb(kPdgNuMu,    kPdgNuMu,    Etrue);
      //double p_nue_numu        = fNuOsc3Flavour.CalcVacuumOscProb(kPdgNuE,     kPdgNuMu,    Etrue);
        double p_numubar_numubar = fNuOsc3Flavour.CalcVacuumOscProb(kPdgNuMuBar, kPdgNuMuBar, Etrue);
        double p_nue_nue         = fNuOsc3Flavour.CalcVacuumOscProb(kPdgNuE,     kPdgNuE,     Etrue);
        double p_numu_nue        = fNuOsc3Flavour.CalcVacuumOscProb(kPdgNuMu,    kPdgNuE,     Etrue);

#ifdef __debug_level_3__
        cout << "Exp/SKMC >>  E = " << Etrue << " GeV" << endl;
        cout << "Exp/SKMC >>  * P (numu     -> numu)    = " << p_numu_numu       << endl;
      //cout << "Exp/SKMC >>  * P (nue      -> numu)    = " << p_nue_numu        << endl;
        cout << "Exp/SKMC >>  * P (numubar  -> numubar) = " << p_numubar_numubar << endl;
        cout << "Exp/SKMC >>  * P (nue      -> nue)     = " << p_nue_nue         << endl;
        cout << "Exp/SKMC >>  * P (numu     -> nue)     = " << p_numu_nue        << endl;
#endif

        // loop over Ereco
        for(int ir = 1; ir <= fMCTemplSys[s][kSKmNuMuCCQE]->GetXaxis()->GetNbins(); ir++) {

           // number of events in (Ereco,Etrue) bin (systematics applied)
           double n_numu_ccqe     = fMCTemplSys[s][kSKmNuMuCCQE]     -> GetBinContent(ir,it);
           double n_numu_cc1pi    = fMCTemplSys[s][kSKmNuMuCC1pi]    -> GetBinContent(ir,it);
           double n_numu_ccoth    = fMCTemplSys[s][kSKmNuMuCCoth]    -> GetBinContent(ir,it);
           double n_numu_nc       = fMCTemplSys[s][kSKmNC]           -> GetBinContent(ir,it);
           double n_numubar_ccqe  = fMCTemplSys[s][kSKmNuMuBarCCQE]  -> GetBinContent(ir,it);
           double n_numubar_cc1pi = fMCTemplSys[s][kSKmNuMuBarCC1pi] -> GetBinContent(ir,it);
           double n_numubar_ccoth = fMCTemplSys[s][kSKmNuMuBarCCoth] -> GetBinContent(ir,it);
           double n_numubar_nc    = fMCTemplSys[s][kSKmNuMuBarNC]    -> GetBinContent(ir,it);
           double n_nue_ccqe      = fMCTemplSys[s][kSKmNuECCQE]      -> GetBinContent(ir,it);
           double n_nue_cc1pi     = fMCTemplSys[s][kSKmNuECC1pi]     -> GetBinContent(ir,it);
           double n_nue_ccoth     = fMCTemplSys[s][kSKmNuECCoth]     -> GetBinContent(ir,it);
           double n_nue_nc        = fMCTemplSys[s][kSKmNuENC]        -> GetBinContent(ir,it);
           double n_signue_ccqe   = fMCTemplSys[s][kSKmSigNuECCQE]   -> GetBinContent(ir,it);
           double n_signue_cc1pi  = fMCTemplSys[s][kSKmSigNuECC1pi]  -> GetBinContent(ir,it);
           double n_signue_ccoth  = fMCTemplSys[s][kSKmSigNuECCoth]  -> GetBinContent(ir,it);

           // apply 3-flavour oscillations
           double n_numu_ccqe_osc     = n_numu_ccqe  * p_numu_numu; 
           double n_numu_cc1pi_osc    = n_numu_cc1pi * p_numu_numu; 
           double n_numu_ccoth_osc    = n_numu_ccoth * p_numu_numu; 
           double n_numu_nc_osc       = n_numu_nc;
           double n_numubar_ccqe_osc  = n_numubar_ccqe  * p_numubar_numubar;
           double n_numubar_cc1pi_osc = n_numubar_cc1pi * p_numubar_numubar;
           double n_numubar_ccoth_osc = n_numubar_ccoth * p_numubar_numubar;
           double n_numubar_nc_osc    = n_numubar_nc;
           double n_nue_ccqe_osc      = n_nue_ccqe  * p_nue_nue; 
           double n_nue_cc1pi_osc     = n_nue_cc1pi * p_nue_nue; 
           double n_nue_ccoth_osc     = n_nue_ccoth * p_nue_nue; 
           double n_nue_nc_osc        = n_nue_nc;
           double n_signue_ccqe_osc   = n_signue_ccqe  * p_numu_nue; 
           double n_signue_cc1pi_osc  = n_signue_cc1pi * p_numu_nue; 
           double n_signue_ccoth_osc  = n_signue_ccoth * p_numu_nue; 

           // store oscillated number of events in current (Ereco,Etrue) bin
           fMCTemplSysHyp[s][kSKmNuMuCCQE]     -> SetBinContent(ir,it, n_numu_ccqe_osc    );
           fMCTemplSysHyp[s][kSKmNuMuCC1pi]    -> SetBinContent(ir,it, n_numu_cc1pi_osc   );
           fMCTemplSysHyp[s][kSKmNuMuCCoth]    -> SetBinContent(ir,it, n_numu_ccoth_osc   );
           fMCTemplSysHyp[s][kSKmNC]           -> SetBinContent(ir,it, n_numu_nc_osc      );
           fMCTemplSysHyp[s][kSKmNuMuBarCCQE]  -> SetBinContent(ir,it, n_numubar_ccqe_osc );
           fMCTemplSysHyp[s][kSKmNuMuBarCC1pi] -> SetBinContent(ir,it, n_numubar_cc1pi_osc);
           fMCTemplSysHyp[s][kSKmNuMuBarCCoth] -> SetBinContent(ir,it, n_numubar_ccoth_osc);
           fMCTemplSysHyp[s][kSKmNuMuBarNC]    -> SetBinContent(ir,it, n_numubar_nc_osc   );
           fMCTemplSysHyp[s][kSKmNuECCQE]      -> SetBinContent(ir,it, n_nue_ccqe_osc     );
           fMCTemplSysHyp[s][kSKmNuECC1pi]     -> SetBinContent(ir,it, n_nue_cc1pi_osc    );
           fMCTemplSysHyp[s][kSKmNuECCoth]     -> SetBinContent(ir,it, n_nue_ccoth_osc    );
           fMCTemplSysHyp[s][kSKmNuENC]        -> SetBinContent(ir,it, n_nue_nc_osc       );
           fMCTemplSysHyp[s][kSKmSigNuECCQE]   -> SetBinContent(ir,it, n_signue_ccqe_osc  );
           fMCTemplSysHyp[s][kSKmSigNuECC1pi]  -> SetBinContent(ir,it, n_signue_cc1pi_osc );
           fMCTemplSysHyp[s][kSKmSigNuECCoth]  -> SetBinContent(ir,it, n_signue_ccoth_osc );
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
      SuperkSample_t s = (SuperkSample_t) is;
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
    double sinsq2theta = fPhysHypothesis.CurrValue(kPhypSinsq2Theta23);
    double alpha       = fPhysHypothesis.CurrValue(kPhypAlphaDecay);
    // loop over samples
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = (SuperkSample_t) is;
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
    double sinsq2theta = fPhysHypothesis.CurrValue(kPhypSinsq2Theta23);
    double musq        = fPhysHypothesis.CurrValue(kPhypMusqDecoherence);
    // loop over samples
    for(int is = 0; is < kns; is++) {
      SuperkSample_t s = (SuperkSample_t) is;
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

#ifdef __debug_level_2__ 
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
     switch(m) {
       case (kSKmNuMuCCQE ):
       case (kSKmNuMuCC1pi):
       case (kSKmNuMuCCoth):
       case (kSKmNC       ):
           flux_10d_weight      = fNuFlux->SuperkNuMuFlux10aTo10dWght(Etrue);
           numu_flux_shape_wght = fNuFlux->SuperkNuMuFluxShapeWeight (Etrue);
           break;
       case(kSKmNuMuBarCCQE ):
       case(kSKmNuMuBarCC1pi):
       case(kSKmNuMuBarCCoth):
       case(kSKmNuMuBarNC   ):
           flux_10d_weight = fNuFlux->SuperkNuMuBarFlux10aTo10dWght(Etrue);
           break;
       case(kSKmNuECCQE ):
       case(kSKmNuECC1pi):
       case(kSKmNuECCoth):
       case(kSKmNuENC   ):
           flux_10d_weight = fNuFlux->SuperkNuEFlux10aTo10dWght(Etrue);
           break;
       case(kSKmSigNuECCQE ):
       case(kSKmSigNuECC1pi):
       case(kSKmSigNuECCoth):
           // signal nue component comes from the numu flux:
           flux_10d_weight      = fNuFlux->SuperkNuMuFlux10aTo10dWght(Etrue);
           numu_flux_shape_wght = fNuFlux->SuperkNuMuFluxShapeWeight (Etrue);
           break;
       default:
           return 0;
     }//m switch
  }

  //
  // Apply cross section systematics = f(Etrue)
  //
  double xsec_ccqe_wght  = 1.;
  double xsec_ratio_wght = 1.;
  if(fGenInputs) {
     xsec_ccqe_wght  = fGenInputs->SuperkCCQEXSecWght (Etrue, m);
     xsec_ratio_wght = fGenInputs->SuperkXSecRatioWght(Etrue, m);
  }

  //
  // Apply cross section systematics = f(Ereco)
  //
  double nucl_effect_wght = 1;
  if(fGenInputs) {
     nucl_effect_wght = fGenInputs->SuperkNuclEffectWght(Ereco,m);
  }

#ifdef __debug_level_2__ 
  if(nucl_effect_wght != 1.) {
    cout << "Exp/SKMC/NumOfEvents >> Applied nucl. effect weight = " << nucl_effect_wght << endl;
  }
#endif

  //
  // Plug-in "efficiency" uncertainty due to fiducial volume, reduction chain, OD cut,
  // ring-counting, PID, momentum cut and decay e- cut.
  // The ring-counting error for CCQE has Ereco-dependence.
  //
  double effwght = this->EfficiencyWeight(s,m,Ereco);

#ifdef __debug_level_2__ 
  if(effwght != 1.) {
    cout << "Exp/SKMC/NumOfEvents >> Applied efficiency weight = " << effwght << endl;
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
  mctempl3 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuMuCCoth];
  mctempl4 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNC];
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
#ifdef __debug_level_2__
    cout << "Exp/SKMC/UpdateEtrueSpectra >> "
         << "1-ring mu / numu : true bin = " << ibin << ", nev = " << nev << endl;
#endif
    spectrum->SetBinContent(it,nev);
  }//true bins

  // update 1-ring mu-like numubar
  spectrum = f1RmNuMuBarEt;
  mctempl1 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE];
  mctempl2 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi];
  mctempl3 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth];
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
#ifdef __debug_level_2__
    cout << "Exp/SKMC/UpdateEtrueSpectra >> "
         << "1-ring mu / numubar : true bin = " << ibin << ", nev = " << nev << endl;
#endif
    spectrum->SetBinContent(it,nev);
  }//true bins


  // update 1-ring mu-like nue
  spectrum = f1RmNuEEt;
  mctempl1 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuECCQE];
  mctempl2 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuECC1pi];
  mctempl3 = fMCTemplSysHyp[kSKsFCFV1RingMuLike][kSKmNuECCoth];
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
#ifdef __debug_level_2__
    cout << "Exp/SKMC/UpdateEtrueSpectra() >> "
         << "1-ring mu / nuebkg : true bin = " << ibin << ", nev = " << nev << endl;
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
  assert(s==kSKsFCFV1RingMuLike);

  switch(m) {
     case (kSKmNuMuCCQE ):
     case (kSKmNuMuCC1pi):
     case (kSKmNuMuCCoth):
     case (kSKmNC       ):
        return fMC1Rmu_numu;
        break;
     case(kSKmNuMuBarCCQE ):
     case(kSKmNuMuBarCC1pi):
     case(kSKmNuMuBarCCoth):
     case(kSKmNuMuBarNC   ):
        return fMC1Rmu_numubar;
        break;
     case(kSKmNuECCQE ):
     case(kSKmNuECC1pi):
     case(kSKmNuECCoth):
     case(kSKmNuENC   ):
        return fMC1Rmu_nue;
        break;
     case(kSKmSigNuECCQE ):
     case(kSKmSigNuECC1pi):
     case(kSKmSigNuECCoth):
        return fMC1Rmu_nuesig;
        break;
     default:
        return 0;
  }
  return 0;
}
//________________________________________________________________________________________
double SuperkMC::MCPOT(SuperkSample_t s, SuperkMode_t m)
{
  assert(s==kSKsFCFV1RingMuLike);

  switch(m) {
     case (kSKmNuMuCCQE ):
     case (kSKmNuMuCC1pi):
     case (kSKmNuMuCCoth):
     case (kSKmNC       ):
        return fMCPOT_numu;
        break;
     case(kSKmNuMuBarCCQE ):
     case(kSKmNuMuBarCC1pi):
     case(kSKmNuMuBarCCoth):
     case(kSKmNuMuBarNC   ):
        return fMCPOT_numubar;
        break;
     case(kSKmNuECCQE ):
     case(kSKmNuECC1pi):
     case(kSKmNuECCoth):
     case(kSKmNuENC   ):
        return fMCPOT_nue;
        break;
     case(kSKmSigNuECCQE ):
     case(kSKmSigNuECC1pi):
     case(kSKmSigNuECCoth):
        return fMCPOT_nuesig;
        break;
     default:
        return 0;
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
#ifdef __debug_level_3__
  cout << "Exp/SKMC >> Tweaking CCQE ring-counting efficiency" << endl;
#endif

  //
  // If this is the first time this function is called, then run the Cholesky 
  // decomposition of the ring-counting error matrix (covariance matrix)
  //
  if(!fEffCCQERelErrCovarCholeskyL) {
    cout << "Exp/SKMC/TweakRingCountingEffCCQE >> Running Cholesky's method." << endl;
    TMatrixD ccqe_rc_err_matrix(*fEffCCQERelErrCovarRC);
#ifdef __debug_level_3__
    cout << "Exp/SKMC/TweakRingCountingEffCCQE >> RC error matrix" << endl;
    ccqe_rc_err_matrix.Print();
#endif

    TMatrixD L = util::CholeskyDecomposition(ccqe_rc_err_matrix);
    fEffCCQERelErrCovarCholeskyL = new TMatrixD(L);

    cout << "Exp/SKMC/TweakRingCountingEffCCQE >> "
         << "Cholesky decomposition of RC error matrix" << endl;
    fEffCCQERelErrCovarCholeskyL->Print();
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
       case (kSKmNuMuCCoth)    : 
       case (kSKmNuMuBarCC1pi) : 
       case (kSKmNuMuBarCCoth) : 
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
       case (kSKmNC)        : 
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
       case (kSKmNuECCoth)    : 
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
      fMC1Rmu_numu->Draw("Enu");
      fMC1Rmu_numu->Draw("Enu","Mode==1","same");
      fMC1Rmu_numu->Draw("Enu","Mode==11||Mode==12||Mode==13","same");
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
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->Draw(draw_opt);
    c->cd(4);
    c->GetPad(4)->SetLogz();
    c->GetPad(4)->SetBottomMargin(0.32);
    c->GetPad(4)->SetRightMargin(0.1);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNC]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNC]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNC]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNC]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNC]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNC]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNC]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNC]->Draw(draw_opt);
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
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->Draw(draw_opt);
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
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplNom[kSKsFCFV1RingMuLike][kSKmNuECCoth]->Draw(draw_opt);
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

#ifdef __debug_level_2__ 
  if(skip) {
    cout << "Exp/SKMC >> Skipping template: " << util::AsString(m) << endl;
  }
#endif

  return skip;
}
//________________________________________________________________________________________





//
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
// everything below is old - remove at some point
//


/*
double SuperkMC::NumOfEvents_Ereco(
    double Ereco, SuperkSample_t s, SuperkMode_t m, bool ignore_energy_scale_twk) const
{
// Count events in the given Ereco bin (integrate template over Etrue)
// - Apply oscillations.
// - Apply flux 10a->10d weights.
// - Apply SuperK efficiency tweaks (if any). 
// - Apply SuperK energy scale tweaks (if any).
// - Apply flux tweaks (if any).
// - Apply cross section tweaks (if any).

  TH2D * nsk2d = fMCTemplNom[s][m];
  assert(nsk2d);

  double nevtot = 0;

  int reco_bin = nsk2d->GetXaxis()->FindBin(Ereco);

#ifdef __debug_level_2__ 
  cout << "Exp/SK/NumOfEvents_Ereco() >> "
       << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m) 
       << ", Ereco = " << Ereco << ", bin = " << reco_bin << endl;
#endif

  //
  // loop over Etrue bins
  //
  for(int true_bin = 1; true_bin <= nsk2d->GetYaxis()->GetNbins(); true_bin++) {

      double Etrue  = nsk2d->GetYaxis()->GetBinCenter(true_bin);
      double nev    = nsk2d->GetBinContent(reco_bin, true_bin);
        
      //
      // Oscillation, decay, decoherence surv probability = f(Etrue)
      //
      double prob_surv = this->SurvivalProbability(Etrue, m);

      //
      // Flux weights & systematics = f(Etrue)
      //
      double flux_10d_weight      = 1.; // 10a->10d weight
      double numu_flux_shape_wght = 1.; // systematic variation (flux shape)
      if(fNuFlux) {
         switch(m) {
           case (kSKmNuMuCCQE ):
           case (kSKmNuMuCC1pi):
           case (kSKmNuMuCCoth):
           case (kSKmNC       ):
               flux_10d_weight      = fNuFlux->SuperkNuMuFlux10aTo10dWght(Etrue);
               numu_flux_shape_wght = fNuFlux->SuperkNuMuFluxShapeWeight (Etrue);
               break;
           case(kSKmNuMuBarCCQE ):
           case(kSKmNuMuBarCC1pi):
           case(kSKmNuMuBarCCoth):
           case(kSKmNuMuBarNC   ):
               flux_10d_weight = fNuFlux->SuperkNuMuBarFlux10aTo10dWght(Etrue);
               break;
           case(kSKmNuECCQE ):
           case(kSKmNuECC1pi):
           case(kSKmNuECCoth):
           case(kSKmNuENC   ):
               flux_10d_weight = fNuFlux->SuperkNuEFlux10aTo10dWght(Etrue);
               break;
           default:
               return 0;
         }//m switch
      }

      //
      // Cross section systematics = f(Etrue)
      //
      double xsec_ccqe_wght  = 1.;
      double xsec_ratio_wght = 1.;
      if(fGenInputs) {
         xsec_ccqe_wght  = fGenInputs->SuperkCCQEXSecWght (Etrue, m);
         xsec_ratio_wght = fGenInputs->SuperkXSecRatioWght(Etrue, m);
      }

      nevtot += ( nev * 
                  prob_surv * 
                  flux_10d_weight * 
                  numu_flux_shape_wght * 
                  xsec_ccqe_wght *
                  xsec_ratio_wght
                );

  }//end of Etrue bin loop

  //
  // Cross section systematics = f(Ereco)
  //
  double nucl_effect_wght = 1;
  if(fGenInputs) {
     nucl_effect_wght = fGenInputs->SuperkNuclEffectWght(Ereco,m);
  }
  nevtot *= nucl_effect_wght;

#ifdef __debug_level_2__ 
  if(nucl_effect_wght != 1.) {
    cout << "Exp/SK/NumOfEvents_Ereco() >> Applied nucl. effect weight = " << nucl_effect_wght << endl;
  }
#endif

  //
  // Plug-in "efficiency" uncertainty due to fiducial volume, reduction chain, OD cut,
  // ring-counting, PID, momentum cut and decay e- cut.
  // The ring-counting error for CCQE has Ereco-dependence.
  //
  double effwght = this->EfficiencyWeight(s,m,Ereco);
  nevtot *= effwght; 

#ifdef __debug_level_2__ 
  if(effwght != 1.) {
    cout << "Exp/SK/NumOfEvents_Ereco() >> Applied efficiency weight = " << effwght << endl;
  }
#endif

  //
  // Plug-in energy scale systematic.
  // Let t be the energy scale tweaking dial and sigE the input energy scale err.
  // Estimate the energy shift of the bin edges.
  // If the upper edge of bin i moves higher by dE, and the bin size is EbinSz 
  // then bin i loses dE/EbinSz fraction of its events to bin i+1.
  // But it gains dE'/EbinSz from bin i-1, where dE' is the shift of the lower bin
  // boundary. Vice-versa if the shift is downwards to lower energies.
  //
  if(!ignore_energy_scale_twk) {
     double escale = this->EnergyScale();
     if( TMath::Abs(escale-1) > kASmallNum) {
        double bin_width   = nsk2d->GetXaxis()->GetBinWidth  (reco_bin);
        double bin_max     = nsk2d->GetXaxis()->GetBinUpEdge (reco_bin);
        double bin_min     = nsk2d->GetXaxis()->GetBinLowEdge(reco_bin);
        if(escale > 1) {
           // Losing events to bin i+1 and gaining events from bin i-1
           double lose_to_higher  = (escale - 1) * bin_max / bin_width;
           double gain_from_lower = (escale - 1) * bin_min / bin_width;
           double nev_prev_bin    = 0;
           if(reco_bin > 1) {
             double Ereco_prev_bin = nsk2d->GetXaxis()->GetBinCenter(reco_bin-1);
             nev_prev_bin = this->NumOfEvents_Ereco(Ereco_prev_bin,s,m,true);
           }
           nevtot *= (1-lose_to_higher);
           nevtot += gain_from_lower * nev_prev_bin;
        }//scale>1
        if(escale < 1) {
           // Losing events to bin i-1 and gaining events from bin i+1
           double lose_to_lower    = (1 - escale) * bin_min / bin_width;
           double gain_from_higher = (1 - escale) * bin_max / bin_width;
           double nev_next_bin    = 0;
           if(reco_bin < nsk2d->GetXaxis()->GetNbins()) {
             double Ereco_next_bin = nsk2d->GetXaxis()->GetBinCenter(reco_bin+1);
             nev_next_bin = this->NumOfEvents_Ereco(Ereco_next_bin,s,m,true);
           }
           nevtot *= (1-lose_to_lower);
           nevtot += gain_from_higher * nev_next_bin;
        }//scale<1
     }//energy scale tweaked?
  }//ignore escale

  return nevtot;
}*/
//________________________________________________________________________________________
/*
double SuperkMC::NumOfEvents_Etrue(
    double Etrue, SuperkSample_t s, SuperkMode_t m) const
{
// Count events in the given Ereco bin (integrate template over Ereco)
// - Apply oscillations.
// - Apply flux 10a->10d weights.
// - Apply SuperK efficiency tweaks (if any).
// - Don't apply SuperK energy scale tweaks (if any): Moves events between Ereco bins so 
//   it doesn't change results since we are integrating over Ereco.
// - Apply flux tweaks (if any).
// - Apply cross sections tweak (if any).

  TH2D * nsk2d = fMCTemplNom[s][m];
  assert(nsk2d);

  double nevtot = 0;

  int true_bin = nsk2d->GetYaxis()->FindBin(Etrue);

#ifdef __debug_level_2__ 
  cout << "Exp/SK/NumOfEvents_Ereco() >> "
       << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m) 
       << ", Etrue = " << Etrue << ", bin = " << true_bin << endl;
#endif

  //
  // loop over Ereco bins
  //
  for(int reco_bin = 1; reco_bin <= nsk2d->GetXaxis()->GetNbins(); reco_bin++) {

      double Ereco = nsk2d->GetXaxis()->GetBinCenter(reco_bin);
      double nev   = nsk2d->GetBinContent(reco_bin, true_bin);

      //
      // Cross section systematics = f(Ereco)
      //
      double nucl_effect_wght = 1;
      if(fGenInputs) {
         nucl_effect_wght = fGenInputs->SuperkNuclEffectWght(Ereco,m);
      }
      nev *= nucl_effect_wght;

      //
      // Plug-in "efficiency" uncertainty due to fiducial volume, reduction chain, OD cut,
      // ring-counting, PID, momentum cut and decay e- cut.
      // The ring-counting error for CCQE has Ereco-dependence.
      //
      double effwght = this->EfficiencyWeight(s,m,Ereco);
      nev *= effwght; 

      nevtot += nev;

  }//end of Ereco bin loop

  //
  // Oscillation, decay, decoherence surv probability = f(Etrue)
  //
  double prob_surv = this->SurvivalProbability(Etrue, m);


  //
  // Flux weights & systematics = f(Etrue)
  //
  double flux_10d_weight      = 1.; // 10a->10d weight
  double numu_flux_shape_wght = 1.; // systematic variation (flux shape)
  if(fNuFlux) {
     switch(m) {
        case (kSKmNuMuCCQE ):
        case (kSKmNuMuCC1pi):
        case (kSKmNuMuCCoth):
        case (kSKmNC       ):
            flux_10d_weight      = fNuFlux->SuperkNuMuFlux10aTo10dWght(Etrue);
            numu_flux_shape_wght = fNuFlux->SuperkNuMuFluxShapeWeight (Etrue);
            break;
        case(kSKmNuMuBarCCQE ):
        case(kSKmNuMuBarCC1pi):
        case(kSKmNuMuBarCCoth):
        case(kSKmNuMuBarNC   ):
            flux_10d_weight = fNuFlux->SuperkNuMuBarFlux10aTo10dWght(Etrue);
            break;
        case(kSKmNuECCQE ):
        case(kSKmNuECC1pi):
        case(kSKmNuECCoth):
        case(kSKmNuENC   ):
            flux_10d_weight = fNuFlux->SuperkNuEFlux10aTo10dWght(Etrue);
            break;
        default:
            return 0;
      }//m switch
  }

  //
  // Cross section systematics = f(Etrue)
  //
  double xsec_ccqe_wght  = 1.;
  double xsec_ratio_wght = 1.;
  if(fGenInputs) {
      xsec_ccqe_wght  = fGenInputs->SuperkCCQEXSecWght (Etrue, m);
      xsec_ratio_wght = fGenInputs->SuperkXSecRatioWght(Etrue, m);
  }

  nevtot *= ( prob_surv * 
              flux_10d_weight * 
              numu_flux_shape_wght * 
              xsec_ccqe_wght *
              xsec_ratio_wght
            );

  return nevtot;
}*/
//________________________________________________________________________________________
/*
void SuperkMC::UpdateErecoSpectra(void)
{
  for(int is = 0; is < kns; is++) {
    TH1D * reco_spectrum_tot =  fMCErecoSpectrumTot[is];
    if(!reco_spectrum_tot) continue; /// loop only over samples included in osc fit

    for(int ibin = 1; ibin <= reco_spectrum_tot->GetNbinsX(); ibin++) {
      double Ereco  = reco_spectrum_tot->GetBinCenter(ibin);
      double NevTot = 0;

      for(int im = 0; im < knm; im++) { 
            TH1D * reco_spectrum =  fMCErecoSpectrum[is][im];
            assert(reco_spectrum);
            SuperkSample_t s = (SuperkSample_t) is;
            SuperkMode_t   m = (SuperkMode_t  ) im;
            double Nev = this->NumOfEvents_Ereco(Ereco,s,m);
            reco_spectrum->SetBinContent(ibin,Nev);
            NevTot += Nev;
      }//im
      reco_spectrum_tot->SetBinContent(ibin,NevTot);
    }//bin
 }//is
} */
//________________________________________________________________________________________
/*
void SuperkMC::UpdateEtrueSpectra(void)
{
//cout << "Exp/SK >> Updating Etrue spectra for calculating beam norm err" << endl;

  TH1D * spectrum = 0;

  // update 1-ring mu-like numu
  spectrum = f1RmNuMuEt;
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
    double Etrue = spectrum->GetBinCenter(ibin);
    double nev = 
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuMuCCQE ) +
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuMuCC1pi) +
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuMuCCoth) +
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNC       );
#ifdef __debug_level_2__
    cout << "Exp/SK/UpdateEtrueSpectra() >> "
         << "1-ring mu / numu : bin = " << ibin << ", nev = " << nev << endl;
#endif
    spectrum->SetBinContent(ibin,nev);
  }//bins

  // update 1-ring mu-like numubar
  spectrum = f1RmNuMuBarEt;
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
    double Etrue = spectrum->GetBinCenter(ibin);
    double nev = 
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuMuBarCCQE)  +
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuMuBarCC1pi) +
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuMuBarCCoth) +
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuMuBarNC);
#ifdef __debug_level_2__
    cout << "Exp/SK/UpdateEtrueSpectra() >> "
         << "1-ring mu / numubar : bin = " << ibin << ", nev = " << nev << endl;
#endif
    spectrum->SetBinContent(ibin,nev);
  }//bins

  // update 1-ring mu-like nue
  spectrum = f1RmNuEEt;
  for(int ibin = 1; ibin <= spectrum->GetNbinsX(); ibin++) {
    double Etrue = spectrum->GetBinCenter(ibin);
    double nev = 
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuECCQE ) +
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuECC1pi) +
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuECCoth) +
       this->NumOfEvents_Etrue(Etrue, kSKsFCFV1RingMuLike, kSKmNuENC);
#ifdef __debug_level_2__
    cout << "Exp/SK/UpdateEtrueSpectra() >> "
         << "1-ring mu / nuebkg : bin = " << ibin << ", nev = " << nev << endl;
#endif
    spectrum->SetBinContent(ibin,nev);
  }//bins
} */
//________________________________________________________________________________________
/*
double SuperkMC::SurvivalProbability(double Etrue, SuperkMode_t m) const
{
  if(fPhysHypothesis.Hypothesis() == kHypNuMuDisapp2F)
  {
     double sinsq2theta23 = fPhysHypothesis.CurrValue(kPhypSinsq2Theta23);
     double dm23sq        = fPhysHypothesis.CurrValue(kPhypDm23sq);
     double prob_surv = util::NuMuSurvProbOsc2D(sinsq2theta23, dm23sq, Etrue, m);
     return prob_surv;
  } 
  else
  if(fPhysHypothesis.Hypothesis() == kHypDecay)
  {
     double sinsq2theta = fPhysHypothesis.CurrValue(kPhypSinsq2Theta23);
     double alpha       = fPhysHypothesis.CurrValue(kPhypAlphaDecay);
     double prob_surv   = util::NuMuSurvProbDecay(sinsq2theta, alpha, Etrue, m);
     return prob_surv;
  } 
  else
  if(fPhysHypothesis.Hypothesis() == kHypDecoh)
  {
     double sinsq2theta = fPhysHypothesis.CurrValue(kPhypSinsq2Theta23);
     double musq        = fPhysHypothesis.CurrValue(kPhypMusqDecoherence);
     double prob_surv   = util::NuMuSurvProbDecay(sinsq2theta, musq, Etrue, m);
     return prob_surv;
  }

  return 1.;
}*/
//________________________________________________________________________________________
