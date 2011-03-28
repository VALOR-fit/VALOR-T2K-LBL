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
      fMCTemplates     [is][im] = 0;
      fMCEtrueSpectrum [is][im] = 0;
      fMCErecoSpectrum [is][im] = 0;
    }
  }
  
  fEnergyScaleRelErr  = 0.;
  fEffCCQERelErr      = 0.;
  fEffCCnonQERelErr   = 0.;
  fEffNCRelErr        = 0.;
  fEffNuECCRelErr     = 0.;
         
  fEnergyScaleTwkDial = 0.;
  fEffCCQETwkDial     = 0.;
  fEffCCnonQETwkDial  = 0.;
  fEffNCTwkDial       = 0.;
  fEffNuECCTwkDial    = 0.;
  
  fMCPOT_numu     = 0;
  fMC1Rmu_numu    = 0;
  fMCPOT_numubar  = 0;
  fMC1Rmu_numubar = 0;
  fMCPOT_nue      = 0;
  fMC1Rmu_nue     = 0;
  fBrEreco        = 0;
  fBrEtrue        = 0;
  fBrPdg          = 0;
  fBrGenCode      = 0;
         
  fNuFlux     = 0;
  fGenInputs      = 0;
  
  f1RmNuMuEt      = 0;
  f1RmNuMuBarEt   = 0;
  f1RmNuEEt       = 0;

  fCurrSinsq2Theta23 = 0;
  fCurrDm23sq        = 0;
  fCurrAlphaDecay    = 0;
  fCurrMusqDecoh     = 0;
  
  this -> UseNeutrinoOscillationModel (true);
  this -> UseNeutrinoDecayModel       (false);
  this -> UseNeutrinoDecoherenceModel (false);
}
//________________________________________________________________________________________
void SuperkMC::CleanUp(void)
{
  if(fMC1Rmu_numu   ) delete fMC1Rmu_numu;
  if(fMC1Rmu_numubar) delete fMC1Rmu_numubar;
  if(fMC1Rmu_nue    ) delete fMC1Rmu_nue;

  for(int is=0; is<kns; is++) {
    if( fMCErecoSpectrumTot[is] ) delete fMCErecoSpectrumTot[is];
    for(int im=0; im<knm; im++) {
      if( fMCTemplates[is][im]     ) delete fMCTemplates[is][im];
      if( fMCEtrueSpectrum[is][im] ) delete fMCEtrueSpectrum[is][im];
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
  cout << "Exp/SK >> Reading SuperK inputs from input XML file." << endl;

  //
  // Get data POT
  // 
  fPOT = util::XML2Dbl(xml,root_node,"superk/observation/pot");
  cout << "Exp/SK >> " << "N_{POT}|data : " << fPOT << endl;
  
  //
  // Get MC trees and MC POT
  //
  fMCPOT_numu = util::XML2Dbl(xml,root_node,"superk/mc/numu/pot");
  fMC1Rmu_numu = util::XML2TTree(xml,root_node,"superk/mc/numu/one_mu_ring");
  cout << "Exp/SK >> " << "N_{POT}|(numu MC) : " << fMCPOT_numu << endl;
  assert(fMC1Rmu_numu);
  fMCPOT_numubar = util::XML2Dbl(xml,root_node,"superk/mc/numubar/pot");
  fMC1Rmu_numubar = util::XML2TTree(xml,root_node,"superk/mc/numubar/one_mu_ring");
  cout << "Exp/SK >> " << "N_{POT}|(numubar MC) : " << fMCPOT_numubar << endl;
  assert(fMC1Rmu_numubar);
  fMCPOT_nue = util::XML2Dbl(xml,root_node,"superk/mc/nuebeam/pot");
  fMC1Rmu_nue = util::XML2TTree(xml,root_node,"superk/mc/nuebeam/one_mu_ring");
  cout << "Exp/SK >> " << "N_{POT}|(beam nue MC) : " << fMCPOT_nue << endl;
  assert(fMC1Rmu_nue);
  
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

  //
  // Systematics
  //

  // Get relative error on energy scale
  fEnergyScaleRelErr = util::XML2Dbl(
     xml,root_node,"superk/systematics/energy_scale_relative_error");
  // Get combined relative error on efficiencies
  fEffCCQERelErr      =
    util::XML2Dbl(xml,root_node,
      "superk/systematics/efficiency_relative_error/one_mu_ring/ccqe");
  fEffCCnonQERelErr   =
    util::XML2Dbl(xml,root_node,
      "superk/systematics/efficiency_relative_error/one_mu_ring/ccnonqe");
  fEffNCRelErr        =
    util::XML2Dbl(xml,root_node,
      "superk/systematics/efficiency_relative_error/one_mu_ring/nc");
  fEffNuECCRelErr     =
    util::XML2Dbl(xml,root_node, 
      "superk/systematics/efficiency_relative_error/one_mu_ring/nuecc");
      
      
  cout << "Exp/SK >> " << "RelErr {EnergyScale}                : " << fEnergyScaleRelErr << endl;
  cout << "Exp/SK >> " << "RelErr {Efficiency/1RingMu/CCQE}    : " << fEffCCQERelErr     << endl;
  cout << "Exp/SK >> " << "RelErr {Efficiency/1RingMu/CCnonQE} : " << fEffCCnonQERelErr  << endl;
  cout << "Exp/SK >> " << "RelErr {Efficiency/1RingMu/NC}      : " << fEffNCRelErr       << endl;
  cout << "Exp/SK >> " << "RelErr {Efficiency/1RingMu/NuECC}   : " << fEffNuECCRelErr    << endl;
  
  // Reset osc params
  this->UseOscillationParams(0., 0.);

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

  cout << "Exp/SK >> Done reading SuperK inputs!" << endl;
  
  return true;
}
//________________________________________________________________________________________
void SuperkMC::CreateTemplate(SuperkSample_t s, SuperkMode_t m)
{
  cout << "Exp/SK >> Filling event Ereco-Etrue histogram for : "
       << util::AsString(s) << "/" << util::AsString(m) << endl;
    
  int is = (int)s;
  int im = (int)m;
    
  fMCTemplates[is][im] = new TH2D(
     Form("nsk_%d_%d",is,im), util::AsString(m).c_str(),
       kNvErecoBins, kNvErecoBinEdges, kNEtrueBins, kEtrueMin, kEtrueMax);
  fMCTemplates[is][im]->SetDirectory(0);
  fMCTemplates[is][im]->GetXaxis()->SetTitle("E_{reco} (GeV)");
  fMCTemplates[is][im]->GetYaxis()->SetTitle("E_{true} (GeV)");
      
  fMCEtrueSpectrum[is][im] = new TH1D(
     Form("nsketrue_%d_%d",is,im), "", kNEtrueBins, kEtrueMin, kEtrueMax);
  fMCEtrueSpectrum[is][im]->SetDirectory(0);
  fMCEtrueSpectrum[is][im]->GetXaxis()->SetTitle("E_{true} (GeV)");
  
  TTree * mctree = this->MCTree(s,m);
  double  mcpot  = this->MCPOT(s,m);

  double wght = this->AnalyzedSamplePOT() / mcpot;
  
  cout << "Exp/SK >> Weight (data POT / MC POT) = " << wght << endl;

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

          bool use_it = false;

          if(m == kSKmNuMuCCQE) {
             bool is_numu = (fBrPdg == kPdgNuMu);
             bool is_ccqe = (fBrGenCode == 1);
             use_it = is_numu && is_ccqe;
          } 
          else
          if(m == kSKmNuMuCC1pi) {
             bool is_numu  = (fBrPdg == kPdgNuMu);
             bool is_cc1pi = (fBrGenCode == 11 || 
                              fBrGenCode == 12 || 
                              fBrGenCode == 13);
            use_it = is_numu && is_cc1pi;
          } 
          else
          if(m == kSKmNuMuCCoth) {
             bool is_numu  = (fBrPdg == kPdgNuMu);
             bool is_ccoth = (fBrGenCode == 16 || 
                              fBrGenCode == 17 || 
                              fBrGenCode == 21 ||
                              fBrGenCode == 22 || 
                              fBrGenCode == 23 || 
                              fBrGenCode == 26);
             use_it = is_numu && is_ccoth;
          } 
          else
          if(m == kSKmNC) {
             bool is_numu  = (fBrPdg == kPdgNuMu);
             bool is_nc    = (fBrGenCode == 31 ||
                              fBrGenCode == 32 ||
                              fBrGenCode == 36 ||
                              fBrGenCode == 33 ||
                              fBrGenCode == 34 ||
                              fBrGenCode == 38 ||
                              fBrGenCode == 39 ||
                              fBrGenCode == 41 ||
                              fBrGenCode == 42 ||
                              fBrGenCode == 43 ||
                              fBrGenCode == 46 ||
                              fBrGenCode == 51 ||
                              fBrGenCode == 52);
             use_it = is_numu && is_nc;
          } 
          else
          if(m == kSKmNuMuBarCCQE) {
             bool is_numubar = (fBrPdg == kPdgNuMuBar);
             bool is_ccqe    = (fBrGenCode == -1);
             use_it = is_numubar && is_ccqe;
          } 
          else
          if(m == kSKmNuMuBarCC1pi) {
             bool is_numubar = (fBrPdg == kPdgNuMuBar);
             bool is_cc1pi   = (fBrGenCode == -11 || fBrGenCode == -12 || fBrGenCode == -13);
             use_it = is_numubar && is_cc1pi;
          } 
          else
          if(m == kSKmNuMuBarCCoth) {
             bool is_numubar = (fBrPdg == kPdgNuMuBar);
             bool is_ccqe    = (fBrGenCode == -1);
             bool is_cc1pi   = (fBrGenCode == -11 || fBrGenCode == -12 || fBrGenCode == -13);
             bool is_nc      = (fBrGenCode == -31 ||
                                fBrGenCode == -32 ||
                                fBrGenCode == -36 ||
                                fBrGenCode == -33 ||
                                fBrGenCode == -34 ||
                                fBrGenCode == -38 ||
                                fBrGenCode == -39 ||
                                fBrGenCode == -41 ||
                                fBrGenCode == -42 ||
                                fBrGenCode == -43 ||
                                fBrGenCode == -44 ||
                                fBrGenCode == -45 ||
                                fBrGenCode == -46 ||
                                fBrGenCode == -51 ||
                                fBrGenCode == -52);
             use_it = is_numubar && !(is_nc||is_ccqe||is_cc1pi);
          } 
          else
          if(m == kSKmNuMuBarNC) {
             bool is_numubar  = (fBrPdg == kPdgNuMuBar);
             bool is_nc    = (fBrGenCode == -31 ||
                              fBrGenCode == -32 ||
                              fBrGenCode == -36 ||
                              fBrGenCode == -33 ||
                              fBrGenCode == -34 ||
                              fBrGenCode == -38 ||
                              fBrGenCode == -39 ||
                              fBrGenCode == -41 ||
                              fBrGenCode == -42 ||
                              fBrGenCode == -43 ||
                              fBrGenCode == -44 ||
                              fBrGenCode == -45 ||
                              fBrGenCode == -46 ||
                              fBrGenCode == -51 ||
                              fBrGenCode == -52);
             use_it = is_numubar && is_nc;
          } 
          else
          if(m == kSKmNuECCQE) {
             bool is_nue  = (fBrPdg == kPdgNuE || fBrPdg == kPdgNuEBar);
             bool is_ccqe = (fBrGenCode == 1 || fBrGenCode == -1);
             use_it = is_nue && is_ccqe;
          } 
          else
          if(m == kSKmNuECC1pi) {
             bool is_nue    = (fBrPdg == kPdgNuE || fBrPdg == kPdgNuEBar);
             bool is_cc1pi  = (fBrGenCode ==  11 || 
                               fBrGenCode ==  12 || 
                               fBrGenCode ==  13 ||
                               fBrGenCode == -11 || 
                               fBrGenCode == -12 || 
                               fBrGenCode == -13);
             use_it = is_nue && is_cc1pi;
          }
          else
          if(m == kSKmNuECCoth) {
             bool is_nue    = (fBrPdg == kPdgNuE || fBrPdg == kPdgNuEBar);
             bool is_ccqe   = (fBrGenCode == 1 || 
                               fBrGenCode == -1);
             bool is_cc1pi  = (fBrGenCode ==  11 || 
                               fBrGenCode ==  12 || 
                               fBrGenCode ==  13 ||
                               fBrGenCode == -11 || 
                               fBrGenCode == -12 || 
                               fBrGenCode == -13);
             bool is_nc    =  (fBrGenCode == 31  ||
                               fBrGenCode == 32  ||
                               fBrGenCode == 36  ||
                               fBrGenCode == 33  ||
                               fBrGenCode == 34  ||
                               fBrGenCode == 38  ||
                               fBrGenCode == 39  ||
                               fBrGenCode == 41  ||
                               fBrGenCode == 42  ||
                               fBrGenCode == 43  ||
                               fBrGenCode == 46  ||
                               fBrGenCode == 51  ||
                               fBrGenCode == 52  ||
                               fBrGenCode == -31 ||
                               fBrGenCode == -32 ||
                               fBrGenCode == -36 ||
                               fBrGenCode == -33 ||
                               fBrGenCode == -34 ||
                               fBrGenCode == -38 ||
                               fBrGenCode == -39 ||
                               fBrGenCode == -41 ||
                               fBrGenCode == -42 ||
                               fBrGenCode == -43 ||
                               fBrGenCode == -44 ||
                               fBrGenCode == -45 ||
                               fBrGenCode == -46 ||
                               fBrGenCode == -51 ||
                               fBrGenCode == -52);
             use_it = is_nue && !(is_ccqe||is_cc1pi||is_nc);
          }
          else
          if(m == kSKmNuENC) {
             bool is_nue  = (fBrPdg == kPdgNuE || fBrPdg == kPdgNuEBar);
             bool is_nc    = (fBrGenCode == 31  ||
                              fBrGenCode == 32  ||
                              fBrGenCode == 36  ||
                              fBrGenCode == 33  ||
                              fBrGenCode == 34  ||
                              fBrGenCode == 38  ||
                              fBrGenCode == 39  ||
                              fBrGenCode == 41  ||
                              fBrGenCode == 42  ||
                              fBrGenCode == 43  ||
                              fBrGenCode == 46  ||
                              fBrGenCode == 51  ||
                              fBrGenCode == 52  ||
                              fBrGenCode == -31 ||
                              fBrGenCode == -32 ||
                              fBrGenCode == -36 ||
                              fBrGenCode == -33 ||
                              fBrGenCode == -34 ||
                              fBrGenCode == -38 ||
                              fBrGenCode == -39 ||
                              fBrGenCode == -41 ||
                              fBrGenCode == -42 ||
                              fBrGenCode == -43 ||
                              fBrGenCode == -44 ||
                              fBrGenCode == -45 ||
                              fBrGenCode == -46 ||
                              fBrGenCode == -51 ||
                              fBrGenCode == -52);
             use_it = is_nue && is_nc;
          }
          if(!use_it) continue;

        //cout << "Exp/SK >> Ereco = " << fBrEreco << ", Etrue = " << fBrEtrue << ", wght = " << wght << endl;

          fMCTemplates    [is][im] -> Fill (fBrEreco, fBrEtrue, wght);
          fMCEtrueSpectrum[is][im] -> Fill (fBrEtrue, wght);
  }//iev
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
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->Draw(draw_opt);
    c->cd(2);
    c->GetPad(2)->SetLogz();
    c->GetPad(2)->SetBottomMargin(0.32);
    c->GetPad(2)->SetRightMargin(0.1);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->Draw(draw_opt);
    c->cd(3);
    c->GetPad(3)->SetLogz();
    c->GetPad(3)->SetBottomMargin(0.32);
    c->GetPad(3)->SetRightMargin(0.1);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->Draw(draw_opt);
    c->cd(4);
    c->GetPad(4)->SetLogz();
    c->GetPad(4)->SetBottomMargin(0.32);
    c->GetPad(4)->SetRightMargin(0.1);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNC]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNC]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNC]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNC]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNC]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNC]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNC]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNC]->Draw(draw_opt);
    c->cd(5);
    c->GetPad(5)->SetLogz();
    c->GetPad(5)->SetBottomMargin(0.32);
    c->GetPad(5)->SetRightMargin(0.1);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->Draw(draw_opt);
    c->cd(6);
    c->GetPad(6)->SetLogz();
    c->GetPad(6)->SetBottomMargin(0.32);
    c->GetPad(6)->SetRightMargin(0.1);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->Draw(draw_opt);
    c->cd(7);
    c->GetPad(7)->SetLogz();
    c->GetPad(7)->SetBottomMargin(0.32);
    c->GetPad(7)->SetRightMargin(0.1);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->Draw(draw_opt);
    c->cd(8);
    c->GetPad(8)->SetLogz();
    c->GetPad(8)->SetBottomMargin(0.32);
    c->GetPad(8)->SetRightMargin(0.1);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->Draw(draw_opt);
    c->cd(9);
    c->GetPad(9)->SetLogz();
    c->GetPad(9)->SetBottomMargin(0.32);
    c->GetPad(9)->SetRightMargin(0.1);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCQE]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCQE]->Draw(draw_opt);
    c->cd(10);
    c->GetPad(10)->SetLogz();
    c->GetPad(10)->SetBottomMargin(0.32);
    c->GetPad(10)->SetRightMargin(0.1);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->Draw(draw_opt);
    c->cd(11);
    c->GetPad(11)->SetLogz();
    c->GetPad(11)->SetBottomMargin(0.32);
    c->GetPad(11)->SetRightMargin(0.1);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCoth]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuECCoth]->Draw(draw_opt);
    c->Update();
    c->cd(12);
    c->GetPad(12)->SetLogz();
    c->GetPad(12)->SetBottomMargin(0.32);
    c->GetPad(12)->SetRightMargin(0.1);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuENC]->GetXaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuENC]->GetXaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuENC]->GetYaxis()->SetLabelSize  (0.12);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuENC]->GetYaxis()->SetTitleSize  (0.14);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuENC]->GetYaxis()->SetTitleOffset(0.40);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuENC]->GetXaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuENC]->GetYaxis()->SetRangeUser(emin,emax);
    fMCTemplates[kSKsFCFV1RingMuLike][kSKmNuENC]->Draw(draw_opt);
    c->Update();
  }

  if(opt.Contains("1RingMuLike/EtrueSpectra")) { 
    const char * draw_opt = "";
    c->Divide(2,6);
    c->cd(1);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCQE]->Draw(draw_opt);
    c->cd(2);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCC1pi]->Draw(draw_opt);
    c->cd(3);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuCCoth]->Draw(draw_opt);
    c->cd(4);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNC]->Draw(draw_opt);
    c->cd(5);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCCQE]->Draw(draw_opt);
    c->cd(6);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCC1pi]->Draw(draw_opt);
    c->cd(7);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarCCoth]->Draw(draw_opt);
    c->cd(8);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuMuBarNC]->Draw(draw_opt);
    c->cd(9);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuECCQE]->Draw(draw_opt);
    c->cd(10);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuECC1pi]->Draw(draw_opt);
    c->cd(11);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuECCoth]->Draw(draw_opt);
    c->cd(12);
    fMCEtrueSpectrum[kSKsFCFV1RingMuLike][kSKmNuENC]->Draw(draw_opt);
    c->Update();
  }

  c->Update();
}
//________________________________________________________________________________________
double SuperkMC::NumOfEvents_Ereco(
    double Ereco, SuperkSample_t s, SuperkMode_t m, bool ignore_energy_scale_twk)
{
// Count events
// Apply oscillations, apply flux 10a->10d weights, 
// apply SuperK efficiency tweaks (if any) and SuperK energy scale tweaks (if any),
// apply flux systematic (shape + norm) tweaks.

  TH2D * nsk2d = fMCTemplates[s][m];
  assert(nsk2d);

  double nevtot = 0;

  int reco_bin = nsk2d->GetXaxis()->FindBin(Ereco);

#ifdef __debug_level_2__ 
  cout << "Exp/SK/NumOfEvents_Ereco() >> "
       << "Sample = " << util::AsString(s) << " / Mode = " << util::AsString(m) 
       << ", Ereco = " << Ereco << ", bin = " << reco_bin << endl;
#endif


  for(int true_bin = 1; true_bin <= nsk2d->GetYaxis()->GetNbins(); true_bin++) {

      double Etrue     = nsk2d->GetYaxis()->GetBinCenter(true_bin);
      double nev       = nsk2d->GetBinContent(reco_bin, true_bin);
      double prob_surv = 1.;
      if(fHaveNuOsc){ 
	prob_surv = util::NuMuSurvProbOsc2D(fCurrSinsq2Theta23, fCurrDm23sq, Etrue, m);
      } else
      if(fHaveNuDecay){ 
	prob_surv = util::NuMuSurvProbDecay(fCurrSinsq2Theta23, fCurrAlphaDecay, Etrue, m);
      } else
      if(fHaveNuDecoh){ 
	prob_surv = util::NuMuSurvProbDecoh(fCurrSinsq2Theta23, fCurrMusqDecoh, Etrue, m);
      }

      //
      // Flux weights & systematics
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
  }//Etrue

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
  // Plug-in "efficiency" uncertainty due to fiducial volume, 
  // reduction chain, ring-counting and PID
  //
  double effwght = this->EfficiencyWeight(s,m);
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
}
//________________________________________________________________________________________
double SuperkMC::NumOfEvents_Etrue(
    double Etrue, SuperkSample_t s, SuperkMode_t m)
{
// Count events in Etrue bin.
// + apply oscillations
// + apply flux 10a->10d weights
// + apply SuperK efficiency tweaks
// + apply flux systematic (shape + norm) tweaks

  int is = (int)s;
  int im = (int)m;

  TH1D * nsk1d = fMCEtrueSpectrum[is][im];
  assert(nsk1d);

  int true_bin = nsk1d->GetXaxis()->FindBin(Etrue);
  double nev = nsk1d->GetBinContent(true_bin);

#ifdef __debug_level_2__
  cout << "Exp/SK/NumOfEvents_Etrue() >> Etrue = " << Etrue << ", bin = " << true_bin << endl;
  cout << "Exp/SK/NumOfEvents_Etrue() >> Nev [unweighted] = " << nev  << endl;
#endif

  // 
  // Oscillations
  //
  double prob_surv = 1;
  if(fHaveNuOsc){ 
    prob_surv = util::NuMuSurvProbOsc2D(fCurrSinsq2Theta23, fCurrDm23sq, Etrue, m);
  } else
  if(fHaveNuDecay){ 
    prob_surv = util::NuMuSurvProbDecay(fCurrSinsq2Theta23, fCurrAlphaDecay, Etrue, m);
  } else
  if(fHaveNuDecoh){ 
    prob_surv = util::NuMuSurvProbDecoh(fCurrSinsq2Theta23, fCurrMusqDecoh, Etrue, m);
  }


  nev *= prob_surv;

#ifdef __debug_level_2__
  cout << "Exp/SK/NumOfEvents_Etrue() >> Nev [wght:osc] = " << nev  << endl;
#endif

  //
  // Flux weights & systematics
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
  nev *= (flux_10d_weight*numu_flux_shape_wght);

#ifdef __debug_level_2__
  cout << "Exp/SK/NumOfEvents_Etrue() >> Nev [wght:osc+flux] = " << nev  << endl;
#endif

  //
  // Cross section systematics = f(Etrue)
  //
  double xsec_ccqe_wght  = 1.;
  double xsec_ratio_wght = 1.;
  if(fGenInputs) {
      xsec_ccqe_wght  = fGenInputs->SuperkCCQEXSecWght (Etrue, m);
      xsec_ratio_wght = fGenInputs->SuperkXSecRatioWght(Etrue, m);
  }
  nev *= (xsec_ccqe_wght*xsec_ratio_wght);

#ifdef __debug_level_2__
  cout << "Exp/SK/NumOfEvents_Etrue() >> Nev [wght:osc+flux+xsec] = " << nev  << endl;
#endif

  //
  // Plug-in "efficiency" uncertainty due to fiducial volume, 
  // reduction chain, ring-counting and PID
  //
  double effwght = this->EfficiencyWeight(s,m);
  nev *= effwght; 

#ifdef __debug_level_2__
  cout << "Exp/SK/NumOfEvents_Etrue() >> Nev [wght:osc+flux+xsec+eff] = " << nev  << endl;
#endif

  return nev;
}
//________________________________________________________________________________________
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
}
//________________________________________________________________________________________
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
}
//________________________________________________________________________________________
/*
double SuperkMC::NumOfEventsReadFromTree(
    double Ereco, double ErecoBinSz, SuperkSample_t s, SuperkMode_t m)
{
  double nev = 0;

  TTree * mctree = this->MCTree(s,m);
  double  mcpot  = this->MCPOT(s,m);  

  for(int iev = 0; iev < mctree->GetEntries(); iev++) {
    mctree->GetEntry(iev);

    if(fBrEreco <  Ereco - 0.5*ErecoBinSz) continue;
    if(fBrEreco >= Ereco + 0.5*ErecoBinSz) continue;

    bool use_it = false;

    if(m == kSKmNuMuCCQE) {
      bool is_numu = (fBrPdg == kPdgNuMu);
      bool is_ccqe = (fBrGenCode == 1);
      use_it = is_numu && is_ccqe;
    } 
    else
    if(m == kSKmNuMuCC1pi) {
      bool is_numu  = (fBrPdg == kPdgNuMu);
      bool is_cc1pi = (fBrGenCode == 11 || fBrGenCode == 12 || fBrGenCode == 13);
      use_it = is_numu && is_cc1pi;
    } 
    else
    if(m == kSKmNuMuCCoth) {
      bool is_numu  = (fBrPdg == kPdgNuMu);
      bool is_ccoth = (fBrGenCode == 16 || 
                       fBrGenCode == 17 || 
                       fBrGenCode == 21 ||
                       fBrGenCode == 22 || 
                       fBrGenCode == 23 || 
                       fBrGenCode == 26);
      use_it = is_numu && is_ccoth;
    } 
    else
    if(m == kSKmNC) {
      bool is_numu  = (fBrPdg == kPdgNuMu);
      bool is_nc    = (fBrGenCode == 31 ||
                       fBrGenCode == 32 ||
                       fBrGenCode == 36 ||
                       fBrGenCode == 33 ||
                       fBrGenCode == 34 ||
                       fBrGenCode == 38 ||
                       fBrGenCode == 39 ||
                       fBrGenCode == 41 ||
                       fBrGenCode == 42 ||
                       fBrGenCode == 43 ||
                       fBrGenCode == 46 ||
                       fBrGenCode == 51 ||
                       fBrGenCode == 52);
      use_it = is_numu && is_nc;
    } 
    else
    if(m == kSKmNuMuBarCCQE) {
      bool is_numubar = (fBrPdg == kPdgNuMuBar);
      bool is_ccqe    = (fBrGenCode == -1);
      use_it = is_numubar && is_ccqe;
    } 
    else
    if(m == kSKmNuMuBarCCnonQE) {
      bool is_numubar = (fBrPdg == kPdgNuMuBar);
      bool is_ccqe    = (fBrGenCode == -1);
      bool is_nc    = (fBrGenCode == 31 ||
                       fBrGenCode == 32 ||
                       fBrGenCode == 36 ||
                       fBrGenCode == 33 ||
                       fBrGenCode == 34 ||
                       fBrGenCode == 38 ||
                       fBrGenCode == 39 ||
                       fBrGenCode == 41 ||
                       fBrGenCode == 42 ||
                       fBrGenCode == 43 ||
                       fBrGenCode == 46 ||
                       fBrGenCode == 51 ||
                       fBrGenCode == 52);
      use_it = is_numubar && !(is_nc||is_ccqe);
    } 
    else
    if(m == kSKmNuMuBarNC) {
      bool is_numubar  = (fBrPdg == kPdgNuMuBar);
      bool is_nc    = (fBrGenCode == 31 ||
                       fBrGenCode == 32 ||
                       fBrGenCode == 36 ||
                       fBrGenCode == 33 ||
                       fBrGenCode == 34 ||
                       fBrGenCode == 38 ||
                       fBrGenCode == 39 ||
                       fBrGenCode == 41 ||
                       fBrGenCode == 42 ||
                       fBrGenCode == 43 ||
                       fBrGenCode == 46 ||
                       fBrGenCode == 51 ||
                       fBrGenCode == 52);
      use_it = is_numubar && is_nc;
    } 
    else
    if(m == kSKmNuECCQE) {
        bool is_nue  = (fBrPdg == kPdgNuE || fBrPdg == kPdgNuEBar);
        bool is_ccqe = (fBrGenCode == 1 || fBrGenCode == -1);
        use_it = is_nue && is_ccqe;
     } 
     else
     if(m == kSKmNuECCnonQE) {
         bool is_nue  = (fBrPdg == kPdgNuE || fBrPdg == kPdgNuEBar);
         bool is_ccqe = (fBrGenCode == 1 || fBrGenCode == -1);
         bool is_nc   = (fBrGenCode == 31  ||
                         fBrGenCode == 32  ||
                         fBrGenCode == 36  ||
                         fBrGenCode == 33  ||
                         fBrGenCode == 34  ||
                         fBrGenCode == 38  ||
                         fBrGenCode == 39  ||
                         fBrGenCode == 41  ||
                         fBrGenCode == 42  ||
                         fBrGenCode == 43  ||
                         fBrGenCode == 46  ||
                         fBrGenCode == 51  ||
                         fBrGenCode == 52  ||
                         fBrGenCode == -31 ||
                         fBrGenCode == -32 ||
                         fBrGenCode == -36 ||
                         fBrGenCode == -33 ||
                         fBrGenCode == -34 ||
                         fBrGenCode == -38 ||
                         fBrGenCode == -39 ||
                         fBrGenCode == -41 ||
                         fBrGenCode == -42 ||
                         fBrGenCode == -43 ||
                         fBrGenCode == -44 ||
                         fBrGenCode == -45 ||
                         fBrGenCode == -46 ||
                         fBrGenCode == -51 ||
                         fBrGenCode == -52);
        use_it = is_nue && !(is_ccqe||is_nc);
    }
    else
    if(m == kSKmNuENC) {
        bool is_nue  = (fBrPdg == kPdgNuE || fBrPdg == kPdgNuEBar);
        bool is_nc   = (fBrGenCode == 31  ||
                        fBrGenCode == 32  ||
                        fBrGenCode == 36  ||
                        fBrGenCode == 33  ||
                        fBrGenCode == 34  ||
                        fBrGenCode == 38  ||
                        fBrGenCode == 39  ||
                        fBrGenCode == 41  ||
                        fBrGenCode == 42  ||
                        fBrGenCode == 43  ||
                        fBrGenCode == 46  ||
                        fBrGenCode == 51  ||
                        fBrGenCode == 52  ||
                        fBrGenCode == -31 ||
                        fBrGenCode == -32 ||
                        fBrGenCode == -36 ||
                        fBrGenCode == -33 ||
                        fBrGenCode == -34 ||
                        fBrGenCode == -38 ||
                        fBrGenCode == -39 ||
                        fBrGenCode == -41 ||
                        fBrGenCode == -42 ||
                        fBrGenCode == -43 ||
                        fBrGenCode == -44 ||
                        fBrGenCode == -45 ||
                        fBrGenCode == -46 ||
                        fBrGenCode == -51 ||
                        fBrGenCode == -52);
        use_it = is_nue && is_nc;
    }
    if(!use_it) continue;

    double prob_surv = 1.;
    if(fHaveNuOsc){ 
      prob_surv = util::NuMuSurvProbOsc2D(fCurrSinsq2Theta23, fCurrDm23sq, fBrEtrue, m);
    } else
    if(fHaveNuDecay){ 
      prob_surv = util::NuMuSurvProbDecay(fCurrSinsq2Theta23, fCurrAlphaDecay, fBrEtrue, m);
    } else
    if(fHaveNuDecoh){ 
      prob_surv = util::NuMuSurvProbDecoh(fCurrSinsq2Theta23, fCurrMusqDecoh, fBrEtrue, m);
    }

    nev += prob_surv;

  }//iev

  nev *= (fPOT/mcpot); // MC POT -> DATA POT

  return nev;
}*/
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
     default:
        return 0;
  }
  return 0;
}
//________________________________________________________________________________________
double SuperkMC::EfficiencyRelErr(
    SuperkSample_t s, SuperkMode_t m) const
{
  if(s == kSKsFCFV1RingMuLike) {
    switch(m) {
       case (kSKmNuMuCCQE      ) : return fEffCCQERelErr;    break;
       case (kSKmNuMuCC1pi     ) : return fEffCCnonQERelErr; break;
       case (kSKmNuMuCCoth     ) : return fEffCCnonQERelErr; break;
       case (kSKmNC            ) : return fEffNCRelErr;      break;
       case (kSKmNuMuBarCCQE   ) : return fEffCCQERelErr;    break;
       case (kSKmNuMuBarCC1pi  ) : return fEffCCnonQERelErr; break;
       case (kSKmNuMuBarCCoth  ) : return fEffCCnonQERelErr; break;
       case (kSKmNuMuBarNC     ) : return fEffNCRelErr;      break;
       case (kSKmNuECCQE       ) : return fEffNuECCRelErr;   break;
       case (kSKmNuECC1pi      ) : return fEffNuECCRelErr;   break;
       case (kSKmNuECCoth      ) : return fEffNuECCRelErr;   break;
       case (kSKmNuENC         ) : return fEffNCRelErr;      break;
       default                   : return 0.;                break;
    }
  }

  return 0.;
}
//________________________________________________________________________________________
double SuperkMC::EfficiencyTwkDial(
    SuperkSample_t s, SuperkMode_t m) const
{
  if(s == kSKsFCFV1RingMuLike) {
    switch(m) {
       case (kSKmNuMuCCQE      ) : return fEffCCQETwkDial;    break;
       case (kSKmNuMuCC1pi     ) : return fEffCCnonQETwkDial; break;
       case (kSKmNuMuCCoth     ) : return fEffCCnonQETwkDial; break;
       case (kSKmNC            ) : return fEffNCTwkDial;      break;
       case (kSKmNuMuBarCCQE   ) : return fEffCCQETwkDial;    break;
       case (kSKmNuMuBarCC1pi  ) : return fEffCCnonQETwkDial; break;
       case (kSKmNuMuBarCCoth  ) : return fEffCCnonQETwkDial; break;
       case (kSKmNuMuBarNC     ) : return fEffNCTwkDial;      break;
       case (kSKmNuECCQE       ) : return fEffNuECCTwkDial;   break;
       case (kSKmNuECC1pi      ) : return fEffNuECCTwkDial;   break;
       case (kSKmNuECCoth      ) : return fEffNuECCTwkDial;   break;
       case (kSKmNuENC         ) : return fEffNCTwkDial;      break;
       default                   : return 0.;                 break;
    }
  }

  return 0.;
}
//________________________________________________________________________________________
double SuperkMC::EfficiencyWeight(
    SuperkSample_t s, SuperkMode_t m) const
{
  double relerr  = this->EfficiencyRelErr (s,m);
  double twkdial = this->EfficiencyTwkDial(s,m);
  double weight  = (1. + twkdial*relerr);
  return weight;
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

