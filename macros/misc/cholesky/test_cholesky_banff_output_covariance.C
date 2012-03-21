//
// ------------------------------------------------------------------------------------------------
// Test BANFF covariance matrix is regenerated after n tweaks
//
// Steve Dennis
// ------------------------------------------------------------------------------------------------
//

void test_cholesky_banff_output_covariance
(
  const char * config_file = "analysis_configuration.xml", /* analysis configuration file */
)
{
  using namespace std;
  using namespace t2kana;
  
  int n_iter = (int)1E5;
  TString * outfile_name = new TString("test_chol_covar.root");
  
  util::SetDefaultStyle();

  const char * config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  cout << "Reading analysis configuration from: " << config << endl;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config);
  
  TMatrixD C = skexp->GetBANFF()->CovMatrix(); 
  
  int ncols = C.GetNcols();
  int nrows = C.GetNrows();

  assert(ncols==nrows);
  
  double * tmparr = new double[nrows];
  
  TMatrixD * cov_matr = new TMatrixD(nrows, nrows);
  cov_matr->Zero();
  
  TVectorD bestfit = skexp->GetBANFF()->BestFitValues();
  
  for (int ii =0 ; ii<n_iter ; ii++) {
    
    skexp->GetBANFF()->Tweak();

    TVectorD curr_vals(nrows);
    
    for (int jj=0; jj < nrows ; jj++) {
      curr_vals[jj] = skexp->GetBANFF()->CurrentValue(conv::kNuisBANFF00+jj);
    }

    for (int jj=0; jj<nrows; jj++) {
      for (int kk=0; kk<nrows; kk++) {
        (*cov_matr)(jj,kk) += (curr_vals(jj)-bestfit(jj)) * (curr_vals(kk)-bestfit(kk));
      }
    }
  }
  
  delete [] tmparr;
  
  *cov_matr *= (1./n_iter);
  
  TCanvas * c1 = new TCanvas("Original");
  c1->cd();
  
  c1->SetTitle("Original BANFF Covariance");

  C.Draw("colz");
  
  cout<< "\033[1;31m==Input Matrix==\033[0m"<<endl;
  C.Print();
  
  TCanvas * c2 = new TCanvas("Recalculated");
  c2->cd();
  c2->SetTitle("Recalculated BANFF Covariance");
  
  cov_matr->Draw("colz");
  
  cout<< "\n\n\033[1;31m==Output Matrix after "<<n_iter<<" tweaks==\033[0m"<<endl;
  
  cov_matr->Print();
  
  if (outfile_name) {
    TFile outfile(*outfile_name, "RECREATE");
    C.Write("BANFF_Cov");
    cov_matr->Write("Recalculated_BANFF_Cov");
    cout<<"Output written to "<<*outfile_name<<endl;
  }
  
  c1->Update();
  c2->Update();
  
}

