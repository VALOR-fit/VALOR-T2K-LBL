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
  TString opt = "SuperK";

  using namespace std;
  using namespace t2kana;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");

  // energy range and binning
  const double Emin  = 0;
  const double Emax  = 4;
  const double dE    = 0.05;
  const int    nbins = (int) ((Emax - Emin)/dE);
  
  if(opt.Contains("Flux")) {
    skexp->GetNuFlux()->Draw(new TCanvas(), "SKFlux/10a/NuMu");
    skexp->GetNuFlux()->Draw(new TCanvas(), "SKFlux/10d/NuMu");
    skexp->GetNuFlux()->Draw(new TCanvas(), "SKFlux/Both10a10d/NuMu");
    skexp->GetNuFlux()->Draw(new TCanvas(), "SKFlux/Weights/10aTo10d/NuMu");
    skexp->GetNuFlux()->Draw(new TCanvas(), "SKFlux/CovarianceMatrix/NuMu");
  }//Flux

  if(opt.Contains("SuperK")) {
  //skexp->GetSuperkMC()->Draw(new TCanvas(), "etrue_spectra");
    skexp->GetSuperkMC()->Draw(new TCanvas(), "1RingMuLike/ErecoEtrueTemplates");
    skexp->GetSuperkMC()->Draw(new TCanvas(), "1RingMuLike/EtrueSpectra");
  }//SuperK

  if(opt.Contains("nd280")) {
    cout << "Nobs(CC-like)       = " << skexp->GetND280Measurement()->NObsCCLike() << endl;
    cout << "Nexp(CC)            = " << skexp->GetND280Measurement()->NExpCC()  << endl;
    cout << "Nexp(NC-bkg)        = " << skexp->GetND280Measurement()->NExpNCBkg()  << endl;
    cout << "Nexp(all)           = " << skexp->GetND280Measurement()->NExpCC() + skexp->GetND280Measurement()->NExpNCBkg() << endl;
    cout << "Obs/Exp             = " << skexp->GetND280Measurement()->CorrectionFactor()  << endl;  
    cout << "Nexp(numu CCQE)     = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCQE)     << endl;
    cout << "Nexp(numu CCRES)    = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCRES)    << endl;
    cout << "Nexp(numu CCCOH)    = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCCOH)    << endl;
    cout << "Nexp(numu CCoth)    = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuCCOth)    << endl;
    cout << "Nexp(numu NCRES)    = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuNCRES)    << endl;
    cout << "Nexp(numu NCCOH)    = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuNCCOH)    << endl;
    cout << "Nexp(numu NCoth)    = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuNCOth)    << endl;
    cout << "Nexp(numubar CCQE)  = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCQE)  << endl;
    cout << "Nexp(numubar CCRES) = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCRES) << endl;
    cout << "Nexp(numubar CCCOH) = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCCOH) << endl;
    cout << "Nexp(numubar CCoth) = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarCCOth) << endl;
    cout << "Nexp(numubar NCRES) = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarNCRES) << endl;
    cout << "Nexp(numubar NCCOH) = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarNCCOH) << endl;
    cout << "Nexp(numubar NCoth) = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuMuBarNCOth) << endl;
    cout << "Nexp(nue CCQE)      = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCQE)      << endl;
    cout << "Nexp(nue CCRES)     = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCRES)     << endl;
    cout << "Nexp(nue CCCOH)     = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCCOH)     << endl;
    cout << "Nexp(nue CCoth)     = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuECCOth)     << endl;
    cout << "Nexp(nue NCRES)     = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuENCRES)     << endl;
    cout << "Nexp(nue NCCOH)     = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuENCCOH)     << endl;
    cout << "Nexp(nue NCoth)     = " << skexp->GetND280Measurement()->NumOfEvents(conv::kNDmNuENCOth)     << endl;

    skexp->GetND280Measurement()->Draw(new TCanvas(), "Saclay/EtrueSpectrum/NuMu/CC");
    skexp->GetND280Measurement()->Draw(new TCanvas(), "Saclay/EtrueSpectrum/NuMu/NC");
    skexp->GetND280Measurement()->Draw(new TCanvas(), "Saclay/EtrueSpectrum/NuMuBar/CC");
    skexp->GetND280Measurement()->Draw(new TCanvas(), "Saclay/EtrueSpectrum/NuMuBar/NC");
    skexp->GetND280Measurement()->Draw(new TCanvas(), "Saclay/EtrueSpectrum/NuE/CC");
    skexp->GetND280Measurement()->Draw(new TCanvas(), "Saclay/EtrueSpectrum/NuE/NC");
    skexp->GetND280Measurement()->Draw(new TCanvas(), "Saclay/Efficiency/NuMu/CC");
    skexp->GetND280Measurement()->Draw(new TCanvas(), "Saclay/Efficiency/NuMu/NC");
    skexp->GetND280Measurement()->Draw(new TCanvas(), "Saclay/Efficiency/NuMuBar/CC");
    skexp->GetND280Measurement()->Draw(new TCanvas(), "Saclay/Efficiency/NuMuBar/NC");
    skexp->GetND280Measurement()->Draw(new TCanvas(), "Saclay/Efficiency/NuE/CC");
    skexp->GetND280Measurement()->Draw(new TCanvas(), "Saclay/Efficiency/NuE/NC");
  }
  if(opt.Contains("NeutrinoGenerator")) {
    //
    //...
    //
  }//NeutrinoGenerator

}

