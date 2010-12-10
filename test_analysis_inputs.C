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
  SuperkObservation * skobs = new SuperkObservation;
  skexp->Read("inputs/analysis_configuration.xml");
  skobs->Read("inputs/analysis_configuration.xml");

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
    skexp->GetSuperkInputs()->Draw(new TCanvas(), "etrue_spectra");
    skexp->GetSuperkInputs()->Draw(new TCanvas(), "1RingMuLike/ErecoEtrueTemplates");
///    skexp->GetSuperkInputs()->Draw(new TCanvas(), "ReconstructionEfficiency");
///    skexp->GetSuperkInputs()->Draw(new TCanvas(), "SelectionEfficiency/OneMuRingTight");
///    skexp->GetSuperkInputs()->Draw(new TCanvas(), "EnergyTransferMatrix/OneMuRingTight");
  }//SuperK


  if(opt.Contains("nd280")) {
    cout << "Nobs(CC-like) = " << skexp->GetND280Inputs()->NObsCCLike() << endl;
    cout << "Nexp(NC-bkg)  = " << skexp->GetND280Inputs()->NExpNCBkg()  << endl;
    cout << "Nexp(CC)      = " << skexp->GetND280Inputs()->NExpCC()     << endl;
    cout << "Nexp(CCQE)    = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuCCQE) << endl;
    cout << "Nexp(CCRES)   = " << skexp->GetND280Inputs()->NumOfEvents(conv::kNDmNuMuCCRES) << endl;
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
//    skexp->GetGeneratorInputs()->Draw(new TCanvas(), "SKSigTotal");
    skexp->GetGeneratorInputs()->Draw(new TCanvas(), "SKSigMode");
//    skexp->GetGeneratorInputs()->Draw(new TCanvas(), "SKFracMode");
//    skexp->GetGeneratorInputs()->Draw(new TCanvas(), "SKErrSigMode");
  }//NeutrinoGenerator

  if(opt.Contains("Observation")) {
    // Draw Ereco spectrum
    skobs->Draw(new TCanvas(), "Ereco/OneMuRingTight");
    // Print-out spectrum values
    double Ereco = Emin;
    while(Ereco < Emax) {
       Ereco+=dE;
       double nev = skobs.NumOfEvents(Ereco);
       cout << "NObs(SK; Ereco = " << Ereco << ") = " << Ereco << endl;
    }
  }//Observation

}

