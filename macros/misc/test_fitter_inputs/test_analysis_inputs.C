//
// ------------------------------------------------------------------------------------------------
// Plot numu disappearane analysis inputs
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

//#include "Conventions.h"

void test_analysis_inputs
(
  TString      opt         = "Flux SuperK BANFF NuXSec", //The information to test. Can be combinations of:  SuperK Flux nd280
  const char * config_file = "analysis_configuration.xml"//"2013a_numu_run1+2+3+4_flux11bv3.2_neut13a_banff13v5prefit_sk+fsi.xml"
)
{
  using namespace std;
  using namespace t2kana;

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);

  // energy range and binning
  const double Emin  = 0;
  const double Emax  = 4;
  const double dE    = 0.05;
  const int    nbins = (int) ((Emax - Emin)/dE);
  
  if(opt.Contains("Flux")) {
    skexp->GetNuFlux()->Draw(new TCanvas(), "Flux/SK/NuMu/TuneWghts");
    skexp->GetNuFlux()->Draw(new TCanvas(), "Flux/SK/NuMuBar/TuneWghts");
    skexp->GetNuFlux()->Draw(new TCanvas(), "Flux/SK/NuE/TuneWghts");
    skexp->GetNuFlux()->Draw(new TCanvas(), "Flux/SK/NuEBar/TuneWghts");
    skexp->GetNuFlux()->Draw(new TCanvas(), "Flux/SK/NuMu");
    skexp->GetNuFlux()->Draw(new TCanvas(), "Flux/SK/NuMuBar");
    skexp->GetNuFlux()->Draw(new TCanvas(), "Flux/SK/NuE");
    skexp->GetNuFlux()->Draw(new TCanvas(), "Flux/SK/NuEBar");
  }//Flux

  if(opt.Contains("SuperK")) {
    /* Not currently implemented
    skexp->GetSuperkMC()->Draw(new TCanvas(), "1RingMuLike/ErecoEtrueTemplates");
    */
    skexp->GetSuperkMC()->Draw(new TCanvas(), "EtrueSpectra/NuMu");
  }//SuperK

  if(opt.Contains("NeutrinoGenerator")) {
    //
    //...
    //
  }//NeutrinoGenerator

  if(opt.Contains("BANFF")) {
    util::PrintCompact(skexp->GetBANFF()->BestFitValues(),15,9,"BANFF fit params");
    skexp->GetBANFF()->CovMatrix().Print("f=%8.4f");
    TMatrixDDiag diag(skexp->GetBANFF()->CovMatrix());
    ofstream banfftable;
    banfftable.open("banff_cov_matrix.tex", ofstream::out | ofstream::trunc);
    util::PrintCompactTex(skexp->GetBANFF()->CovMatrix(), 10, 3, banfftable, "banff");
    cout << endl << "SQRT DIAG" << endl;
    for(int i = 0; i < diag.GetNdiags(); i++) {
      cout << std::fixed << std::setprecision(3) << TMath::Sqrt(diag(i)) << endl;
    }
    cout << endl << "Best fit values" << endl;
    util::PrintCompact(skexp->GetBANFF()->BestFitValues(),10,3,"BANFF fit params");
  }//BANFF

  if(opt.Contains("NuXSec")) {
    for(int ip = 0; ip < conv::kNumOfXSecSyst; ip++) {
      conv::PhysParamId_t p = conv::kListOfXSecSyst[ip];
      cout << p << "\t" << util::AsString(p) << endl;
      if(p == conv::kPhypUndefined) continue;
      if(conv::kXSecSystIncludedInBANFF[ip]) continue;
      if(conv::kXSecSystHasWghtSplines[ip]) continue;
      cout << "Relative error for XSec systematic " << util::AsString(p) << " = " << skexp->GetNuXSec()->RelativeError(p) << endl;
    }
  }//NuXSec

  if(opt.Contains("SKDetFSI")) {
   skexp->GetSuperkDetectorFSI()->CovMatrix().Print("f=%8.4f");
   TMatrixDDiag diag(skexp->GetSuperkDetectorFSI()->CovMatrix());
   ofstream skdetfsitable;
   skdetfsitable.open("skdetfsi_cov_matrix.tex", ofstream::out | ofstream::trunc);
   util::PrintCompactTex(skexp->GetSuperkDetectorFSI()->CovMatrix(), 20, 5, skdetfsitable, "skdetfsi");
   cout << fixed << setprecision(3);
   for(int i = 0; i < diag.GetNdiags(); i++) {
     cout << TMath::Sqrt(diag(i)) << endl;
   }
  }//SKDetFSI
  else if(opt.Contains("FSI")) {
   skexp->GetFSI()->CovMatrix().Print();
  }//FSI
}

