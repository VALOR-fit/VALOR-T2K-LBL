//
// ------------------------------------------------------------------------------------------------
// Generate systematics summary.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Conventions.h"

void summarize_systematics
(
  string  file_path,            /* Path to data files produced by bestfit_systematic_study_fitting_job.C */
  double  input_sinsq2theta23,  /* True sinsq2theta23 used in bestfit_systematic_study_fitting_job.C     */
  double  input_dm23sq          /* True dm23sq used in bestfit_systematic_study_fitting_job.C            */
)
{
  using namespace std;
  using namespace t2kana;

  // chain together files for each systematic

  const int nsyst = 14;

  double dm23sq_frac_err        [nsyst];
  double sinsq2theta23_frac_err [nsyst];

  TH1D * h_dm23sq_frac_err        = new TH1D("h_dm23sq_frac_err",        "", nsyst, 0, nsyst);
  TH1D * h_sinsq2theta23_frac_err = new TH1D("h_sinsq2theta23_frac_err", "", nsyst, 0, nsyst);
  TH1D * h_nfits                  = new TH1D("h_nfits",                  "", nsyst, 0, nsyst);

  for(int isyst = 0; isyst < nsyst; isyst++) {

    string syst_label="";
    if (isyst ==  0)  syst_label = "sk_ccqe_other_eff";
    if (isyst ==  1)  syst_label = "sk_ccqe_ring_counting_eff";
    if (isyst ==  2)  syst_label = "sk_ccnonqe_eff";
    if (isyst ==  3)  syst_label = "sk_nc_eff";
    if (isyst ==  4)  syst_label = "sk_nuecc_eff";
    if (isyst ==  5)  syst_label = "sk_energy_scale";
    if (isyst ==  6)  syst_label = "nd280_eff";
    if (isyst ==  7)  syst_label = "flux_norm";
    if (isyst ==  8)  syst_label = "flux_shape";
    if (isyst ==  9)  syst_label = "xsec_ccqe";
    if (isyst == 10)  syst_label = "xsec_cc1pi_to_ccqe";
    if (isyst == 11)  syst_label = "xsec_ccoth_to_ccqe";
    if (isyst == 12)  syst_label = "xsec_nc_to_ccqe";
    if (isyst == 13)  syst_label = "nuclear_effect";

    TChain * chain = new TChain("sysnt");

    TString filename_template = 
       Form("%s/syst_%s-sinsq2theta23_%.2f-dm23sq_%.5feV2_*.root",
            file_path.c_str(), syst_label.c_str(),input_sinsq2theta23,input_dm23sq);

    chain->Add(filename_template.Data());
    cout << "Adding files :  " << filename_template << endl;

    int nfits = chain->GetEntries();
    if(nfits==0) {
	delete chain;
	continue;
    }

    TH1D * hdm23sq        = new TH1D("hdm23sq",        "", 1000, 0.1*input_dm23sq, 3*input_dm23sq);
    TH1D * hsinsq2theta23 = new TH1D("hsinsq2theta23", "",  100, 0., 1.);

    chain->Draw("dm23sq>>hdm23sq",               "", "goff");
    chain->Draw("sinsq2theta23>>hsinsq2theta23", "", "goff");

    dm23sq_frac_err       [isyst] = hdm23sq       ->GetRMS() / input_dm23sq;
    sinsq2theta23_frac_err[isyst] = hsinsq2theta23->GetRMS() / input_sinsq2theta23;

    cout << "* Systematic: " << syst_label << endl;
    cout << "  - number of fits = " << nfits << endl;
    cout << "  - syst_err(dmsq) = "          << 100*dm23sq_frac_err[isyst] << " %" << endl;
    cout << "  - syst_err(sinsq2theta23) = " << 100*sinsq2theta23_frac_err[isyst] << " %" << endl;

    int syst_bin = isyst+1;
    h_dm23sq_frac_err        -> SetBinContent(syst_bin, 100*dm23sq_frac_err[isyst]); 
    h_sinsq2theta23_frac_err -> SetBinContent(syst_bin, 100*sinsq2theta23_frac_err[isyst]);
    h_nfits                  -> SetBinContent(syst_bin, nfits);
    h_dm23sq_frac_err        -> GetXaxis() -> SetBinLabel(syst_bin, syst_label.c_str());
    h_sinsq2theta23_frac_err -> GetXaxis() -> SetBinLabel(syst_bin, syst_label.c_str());
    h_nfits                  -> GetXaxis() -> SetBinLabel(syst_bin, syst_label.c_str());

    delete chain;
    delete hdm23sq;
    delete hsinsq2theta23;

  }//isyst

  gStyle->SetOptStat(0);
//h_dm23sq_frac_err        -> GetXaxis() -> SetTitle("Systematic");
//h_sinsq2theta23_frac_err -> GetXaxis() -> SetTitle("Systematic");
  h_dm23sq_frac_err        -> GetYaxis() -> SetTitle("#delta#Deltam^{2}_{23}/#Deltam^{2}_{23} (%)");
  h_sinsq2theta23_frac_err -> GetYaxis() -> SetTitle("#deltasin^{2}2#theta_{23}/sin^{2}2#theta_{23} (%)");
  TCanvas * c1 = new TCanvas("c1","",10,10,400,400);
  h_dm23sq_frac_err->Draw();
  TCanvas * c2 = new TCanvas("c2","",100,100,500,500);
  h_sinsq2theta23_frac_err->Draw();
  TCanvas * c3 = new TCanvas("c3","",200,200,600,600);
  h_nfits->Draw();

  double dm23sq_combined_frac_err        = 0;
  double sinsq2theta23_combined_frac_err = 0;

  for(int isyst = 0; isyst < nsyst; isyst++) {
      dm23sq_combined_frac_err        += TMath::Power(dm23sq_frac_err[isyst]  ,      2.);
      sinsq2theta23_combined_frac_err += TMath::Power(sinsq2theta23_frac_err[isyst], 2.);
  }
  dm23sq_combined_frac_err        += TMath::Sqrt(dm23sq_combined_frac_err);
  sinsq2theta23_combined_frac_err += TMath::Sqrt(sinsq2theta23_combined_frac_err);

  cout << "* Combined systematics: " << endl;
  cout << "  - syst_err(dmsq) = "          << 100*dm23sq_combined_frac_err << " %" << endl;
  cout << "  - syst_err(sinsq2theta23) = " << 100*sinsq2theta23_combined_frac_err << " %" << endl;
}
