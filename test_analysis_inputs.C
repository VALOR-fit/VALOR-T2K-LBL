//
// ------------------------------------------------------------------------------------------------
// Plot numu disappearane analysis inputs
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void test_analysis_inputs()
{
  TString opt = "nd280";

  gSystem->Load("libGpad");  
  gSystem->Load("libPhysics");  
  gSystem->Load("libMatrix");   
  gSystem->Load("libHist");     
  gSystem->Load("libGraf");     
  gSystem->Load("libTree");     
  gSystem->Load("libRIO");   
  gSystem->Load("libXMLIO");
  gSystem->Load("libMinuit");   

  gSystem->Load("$T2KANA/lib/libT2KAnaYr1NuMuOscAnalysis.so");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::yr1numu;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // energy range and binning
  const double Emin  = 0;
  const double Emax  = 4;
  const double dE    = 0.05;
  const int    nbins = (int) ((Emax - Emin)/dE);
  
  if(opt.Contains("Flux")) {
    skexp->GetBeamInputs()->Draw(new TCanvas(), "SKFlux/10a/NuMu");
    skexp->GetBeamInputs()->Draw(new TCanvas(), "SKFlux/10d/NuMu");
    skexp->GetBeamInputs()->Draw(new TCanvas(), "SKFlux/Both10a10d/NuMu");
    skexp->GetBeamInputs()->Draw(new TCanvas(), "SKFlux/Weights/10aTo10d/NuMu");
    skexp->GetBeamInputs()->Draw(new TCanvas(), "SKFlux/CovarianceMatrix/NuMu");
  }//Flux

  if(opt.Contains("SuperK")) {
  //skexp->GetSuperkInputs()->Draw(new TCanvas(), "etrue_spectra");
    skexp->GetSuperkInputs()->Draw(new TCanvas(), "1RingMuLike/ErecoEtrueTemplates");
    skexp->GetSuperkInputs()->Draw(new TCanvas(), "1RingMuLike/EtrueSpectra");
  }//SuperK

  if(opt.Contains("nd280")) {
    cout << "Nobs(CC-like) = " << skexp->GetND280Inputs()->NObsCCLike() << endl;
    cout << "Nexp(CC)      = " << skexp->GetND280Inputs()->NExpCC()  << endl;
    cout << "Nexp(NC-bkg)  = " << skexp->GetND280Inputs()->NExpNCBkg()  << endl;
    cout << "Nexp(all)     = " <<  skexp->GetND280Inputs()->NExpCC() + skexp->GetND280Inputs()->NExpNCBkg() << endl;
    cout << "Obs/Exp       = " << skexp->GetND280Inputs()->CorrectionFactor()  << endl;  
    cout << "Nexp(numu CCQE)     = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuCCQE) << endl;
    cout << "Nexp(numu CCRES)    = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuCCRES) << endl;
    cout << "Nexp(numu CCCOH)    = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuCCCOH) << endl;
    cout << "Nexp(numu CCoth)    = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuCCOth) << endl;
    cout << "Nexp(numu NCRES)    = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuNCRES) << endl;
    cout << "Nexp(numu NCCOH)    = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuNCCOH) << endl;
    cout << "Nexp(numu NCoth)    = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuNCOth) << endl;
    cout << "Nexp(numubar CCQE)  = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuBarCCQE) << endl;
    cout << "Nexp(numubar CCRES) = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuBarCCRES) << endl;
    cout << "Nexp(numubar CCCOH) = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuBarCCCOH) << endl;
    cout << "Nexp(numubar CCoth) = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuBarCCOth) << endl;
    cout << "Nexp(numubar NCRES) = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuBarNCRES) << endl;
    cout << "Nexp(numubar NCCOH) = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuBarNCCOH) << endl;
    cout << "Nexp(numubar NCoth) = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuBarNCOth) << endl;
    cout << "Nexp(nue CCQE)      = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuECCQE) << endl;
    cout << "Nexp(nue CCRES)     = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuECCRES) << endl;
    cout << "Nexp(nue CCCOH)     = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuECCCOH) << endl;
    cout << "Nexp(nue CCoth)     = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuECCOth) << endl;
    cout << "Nexp(nue NCRES)     = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuENCRES) << endl;
    cout << "Nexp(nue NCCOH)     = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuENCCOH) << endl;
    cout << "Nexp(nue NCoth)     = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuENCOth) << endl;

    skexp->GetND280Inputs()->Draw(new TCanvas(), "Saclay/EtrueSpectrum/NuMu/CC");
    skexp->GetND280Inputs()->Draw(new TCanvas(), "Saclay/EtrueSpectrum/NuMu/NC");
    skexp->GetND280Inputs()->Draw(new TCanvas(), "Saclay/EtrueSpectrum/NuMuBar/CC");
    skexp->GetND280Inputs()->Draw(new TCanvas(), "Saclay/EtrueSpectrum/NuMuBar/NC");
    skexp->GetND280Inputs()->Draw(new TCanvas(), "Saclay/EtrueSpectrum/NuE/CC");
    skexp->GetND280Inputs()->Draw(new TCanvas(), "Saclay/EtrueSpectrum/NuE/NC");
    skexp->GetND280Inputs()->Draw(new TCanvas(), "Saclay/Efficiency/NuMu/CC");
    skexp->GetND280Inputs()->Draw(new TCanvas(), "Saclay/Efficiency/NuMu/NC");
    skexp->GetND280Inputs()->Draw(new TCanvas(), "Saclay/Efficiency/NuMuBar/CC");
    skexp->GetND280Inputs()->Draw(new TCanvas(), "Saclay/Efficiency/NuMuBar/NC");
    skexp->GetND280Inputs()->Draw(new TCanvas(), "Saclay/Efficiency/NuE/CC");
    skexp->GetND280Inputs()->Draw(new TCanvas(), "Saclay/Efficiency/NuE/NC");
  }
  if(opt.Contains("NeutrinoGenerator")) {
    //
    //...
    //
  }//NeutrinoGenerator

}

