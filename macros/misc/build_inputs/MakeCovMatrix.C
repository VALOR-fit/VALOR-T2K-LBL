
//
// ------------------------------------------------------------------------------------------------
// Produce the post-Banff covariance matrix (flux+xsec) with the right format for the fitter
//
// Davide Sgalaberna
// ------------------------------------------------------------------------------------------------
//


#include <iomanip.h>

#define ADDCOVMATR

void MakeCovMatrix( 
		   
		   const char *filename = "$T2KOSCANA/inputs/test/PrelimCovMatrix/flux_covariance_13av1_prelim.root",
		   const char *covname = "total_flux_cov",
		   
		   const char *outfilename = "$T2KOSCANA/outputs/test/flux_prelim_13av1_xsec_postbanff2013.root",
		   

		   const bool CorrMode = false, // true -> correlated FHC-RHC, false -> not correlated

		   // Select the sk sub-matrices:
		   //fhc-only
		   const int first_sk_fhc = 50, 
		   const int last_sk_fhc = 74,
		   // rhc-only
		   const int first_sk_rhc = 75, 
		   const int last_sk_rhc = 99,
		   // fhc+rhc
		   const int first_sk_all = 50, 
		   const int last_sk_all = 99,

		   // Add a piece of covariance matrix
		   const bool AddCovMatr = true,
		   const char *addfilename = "$T2KOSCANA/inputs/official/2013/banff/postfit_banff_2013_joint_v5.root",
		   const char *addcovname = "postfit_cov",
		   // rhc-only
		   const int first_add = 25, 
		   const int last_add = 32,
		   
		   ){
  
  using namespace std;
  using namespace valor;

  TFile *filein = new TFile(filename,"READ");
  
  cout << endl;
  cout << "Take covariance matrix " << covname << " from " << filename << endl;
  cout << endl;
  cout << "Covariance Matrix: ";
  if(CorrMode) 
      cout << " correlated" << endl;
  else 
    cout << " Not correlated" << endl;
  
  TMatrixD *cov = (TMatrixD*)filein->Get(covname);

  cout << endl;
  cout << "FHC mode:" << endl;

  // sk numu
  TAxis *sk_numu_Axis = (TAxis*)filein->Get("sk_numode_numu_bins");
  const TArrayD *sk_numu_bins = (TArrayD*)sk_numu_Axis->GetXbins();
  double sk_numu_bins_N = sk_numu_bins->GetSize();
  cout << "  sk_numu_bins ( " << sk_numu_bins_N - 1 << " ):  ";
  for(int i=0;i<sk_numu_bins_N;i++) cout << sk_numu_bins->At(i) << ", ";
  cout << endl;

  // sk numubar
  TAxis *sk_numubar_Axis = (TAxis*)filein->Get("sk_numode_numub_bins");
  const TArrayD *sk_numubar_bins = (TArrayD*)sk_numubar_Axis->GetXbins();
  double sk_numubar_bins_N = sk_numubar_bins->GetSize();
  cout << "  sk_numubar_bins ( " << sk_numubar_bins_N - 1 << " ):  ";
  for(int i=0;i<sk_numubar_bins_N;i++) cout << sk_numubar_bins->At(i) << ", ";
  cout << endl;

  // sk nue
  TAxis *sk_nue_Axis = (TAxis*)filein->Get("sk_numode_nue_bins");
  const TArrayD *sk_nue_bins = (TArrayD*)sk_nue_Axis->GetXbins();
  double sk_nue_bins_N = sk_nue_bins->GetSize();
  cout << "  sk_nue_bins ( " << sk_nue_bins_N - 1 << " ):  ";
  for(int i=0;i<sk_nue_bins_N;i++) cout << sk_nue_bins->At(i) << ", ";
  cout << endl;

  // sk nuebar
  TAxis *sk_nuebar_Axis = (TAxis*)filein->Get("sk_numode_nueb_bins");
  const TArrayD *sk_nuebar_bins = (TArrayD*)sk_nuebar_Axis->GetXbins();
  double sk_nuebar_bins_N = sk_nuebar_bins->GetSize();
  cout << "  sk_nuebar_bins ( " << sk_nuebar_bins_N - 1 << " ):  ";
  for(int i=0;i<sk_nuebar_bins_N;i++) cout << sk_nuebar_bins->At(i) << ", ";
  cout << endl;

  cout << endl;
  cout << "RHC mode:" << endl;

  // sk numu
  TAxis *sk_numu_Axis = (TAxis*)filein->Get("sk_anumode_numu_bins");
  const TArrayD *sk_numu_bins = (TArrayD*)sk_numu_Axis->GetXbins();
  double sk_numu_bins_N = sk_numu_bins->GetSize();
  cout << "  sk_numu_bins ( " << sk_numu_bins_N - 1 << " ):  ";
  for(int i=0;i<sk_numu_bins_N;i++) cout << sk_numu_bins->At(i) << ", ";
  cout << endl;

  // sk numubar
  TAxis *sk_numubar_Axis = (TAxis*)filein->Get("sk_anumode_numub_bins");
  const TArrayD *sk_numubar_bins = (TArrayD*)sk_numubar_Axis->GetXbins();
  double sk_numubar_bins_N = sk_numubar_bins->GetSize();
  cout << "  sk_numubar_bins ( " << sk_numubar_bins_N - 1 << " ):  ";
  for(int i=0;i<sk_numubar_bins_N;i++) cout << sk_numubar_bins->At(i) << ", ";
  cout << endl;

  // sk nue
  TAxis *sk_nue_Axis = (TAxis*)filein->Get("sk_anumode_nue_bins");
  const TArrayD *sk_nue_bins = (TArrayD*)sk_nue_Axis->GetXbins();
  double sk_nue_bins_N = sk_nue_bins->GetSize();
  cout << "  sk_nue_bins ( " << sk_nue_bins_N - 1 << " ):  ";
  for(int i=0;i<sk_nue_bins_N;i++) cout << sk_nue_bins->At(i) << ", ";
  cout << endl;

  // sk nuebar
  TAxis *sk_nuebar_Axis = (TAxis*)filein->Get("sk_anumode_nueb_bins");
  const TArrayD *sk_nuebar_bins = (TArrayD*)sk_nuebar_Axis->GetXbins();
  double sk_nuebar_bins_N = sk_nuebar_bins->GetSize();
  cout << "  sk_nuebar_bins ( " << sk_nuebar_bins_N - 1 << " ):  ";
  for(int i=0;i<sk_nuebar_bins_N;i++) cout << sk_nuebar_bins->At(i) << ", ";
  cout << endl;

  filein->Close();

  // Take sk sub-matrices

  int firstcol = 50;
  int firstrow = 50;
  int lastcol = cov->GetNcols()-1;
  int lastrow = cov->GetNrows()-1;

  TMatrixD *cov_sk_all = (TMatrixD*) cov->GetSub(first_sk_all,last_sk_all,first_sk_all,last_sk_all); 
  TMatrixD *cov_sk_fhc = (TMatrixD*) cov->GetSub(first_sk_fhc,last_sk_fhc,first_sk_fhc,last_sk_fhc); 
  TMatrixD *cov_sk_rhc = (TMatrixD*) cov->GetSub(first_sk_rhc,last_sk_rhc,first_sk_rhc,last_sk_rhc);


  
  if(CorrMode){
    cout << "Remove correlations" << endl;
    for(int irow=0;irow<cov_sk_all){
      for(){
      }
    }
  }


  // Add matrix 
#ifdef ADDCOVMATR
  
  cout << endl;
  cout << "Add piece of covariance matrix w/o correlations from " << addcovname << endl;
  
  TFile *addfilein = new TFile(addfilename,"READ");
  TMatrixD *cov_sk_add = (TMatrixD*) addfilein->Get("postfit_cov");
  addfilein->Close();

  TMatrixD *cov_sk_add_piece = (TMatrixD*) cov_sk_add->GetSub(first_add,last_add,first_add,last_add);
  
  int nbins_sk_all_final = cov_sk_add_piece->GetNrows() + cov_sk_all->GetNrows();
  int nbins_sk_all = cov_sk_all->GetNrows();
  int nbins_sk_add_piece = cov_sk_add_piece->GetNrows();
  TMatrixD *cov_sk_all_final = new TMatrixD(nbins_sk_all_final,nbins_sk_all_final);

  // Copy flux matrix and additional matrix
  for(int irow=0;irow<nbins_sk_all;irow++)
    for(int icol=0;icol<nbins_sk_all;icol++)
      cov_sk_all_final(irow,icol) = cov_sk_all(irow,icol);
  for(int irow=0;irow<nbins_sk_add_piece;irow++)
    for(int icol=0;icol<nbins_sk_add_piece;icol++)
      cov_sk_all_final(irow+nbins_sk_all,icol+nbins_sk_all) = cov_sk_add_piece(irow,icol);
  

  int nbins_sk_fhc_final = cov_sk_add_piece->GetNrows() + cov_sk_fhc->GetNrows();
  int nbins_sk_fhc = cov_sk_fhc->GetNrows();
  TMatrixD *cov_sk_fhc_final = new TMatrixD(nbins_sk_fhc_final,nbins_sk_fhc_final);
  // Copy flux matrix and additional matrix
  for(int irow=0;irow<nbins_sk_fhc;irow++)
    for(int icol=0;icol<nbins_sk_fhc;icol++)
      cov_sk_fhc_final(irow,icol) = cov_sk_fhc(irow,icol);
  for(int irow=0;irow<nbins_sk_add_piece;irow++)
    for(int icol=0;icol<nbins_sk_add_piece;icol++)
      cov_sk_fhc_final(irow+nbins_sk_fhc,icol+nbins_sk_fhc) = cov_sk_add_piece(irow,icol);


  int nbins_sk_rhc_final = cov_sk_add_piece->GetNrows() + cov_sk_rhc->GetNrows();
  int nbins_sk_rhc = cov_sk_rhc->GetNrows();
  TMatrixD *cov_sk_rhc_final = new TMatrixD(nbins_sk_rhc_final,nbins_sk_rhc_final);
  // Copy flux matrix and additional matrix
  for(int irow=0;irow<nbins_sk_rhc;irow++)
    for(int icol=0;icol<nbins_sk_rhc;icol++)
      cov_sk_rhc_final(irow,icol) = cov_sk_rhc(irow,icol);
  for(int irow=0;irow<nbins_sk_add_piece;irow++)
    for(int icol=0;icol<nbins_sk_add_piece;icol++)
      cov_sk_rhc_final(irow+nbins_sk_rhc,icol+nbins_sk_rhc) = cov_sk_add_piece(irow,icol);
  
#endif 
#ifndef ADDCOVMATR 
    
  TMatrixD *cov_sk_all_final = (TMatrixD*) cov_sk_all->Clone("");
  TMatrixD *cov_sk_fhc_final = (TMatrixD*) cov_sk_fhc->Clone("");
  TMatrixD *cov_sk_rhc_final = (TMatrixD*) cov_sk_rhc->Clone("");
  
#endif

  // Write output covariance matrices in a file
  TFile *outfile = new TFile(outfilename,"RECREATE");
  cov_sk_all_final->Write("prelim_cov_all");
  cov_sk_fhc_final->Write("prelim_cov_fhc");
  cov_sk_rhc_final->Write("prelim_cov_rhc");
  outfile->Close();


#ifdef ADDCOVMATR
  TCanvas *ccov_sk_add_piece = new TCanvas("ccov_sk_add_piece","ccov_sk_add_piece");
  cov_sk_add_piece->Draw("colz");
#endif
  
#ifdef ADDCOVMATR
  TCanvas *ccov_sk_all = new TCanvas("ccov_sk_all","ccov_sk_all");
  cov_sk_all->Draw("colz");
#endif
  TCanvas *ccov_sk_all_final = new TCanvas("ccov_sk_all_final","ccov_sk_all_final");
  cov_sk_all_final->Draw("colz");

#ifdef ADDCOVMATR
  TCanvas *ccov_sk_fhc = new TCanvas("ccov_sk_fhc","ccov_sk_fhc");
  cov_sk_fhc->Draw("colz");
#endif
  TCanvas *ccov_sk_fhc_final = new TCanvas("ccov_sk_fhc_final","ccov_sk_fhc_final");
  cov_sk_fhc_final->Draw("colz");

#ifdef ADDCOVMATR
  TCanvas *ccov_sk_rhc = new TCanvas("ccov_sk_rhc","ccov_sk_rhc");
  cov_sk_rhc->Draw("colz");
#endif
  TCanvas *ccov_sk_rhc_final = new TCanvas("ccov_sk_rhc_final","ccov_sk_rhc_final");
  cov_sk_rhc_final->Draw("colz");
  

 
  //cov_sk_add->Print();
  //cov_sk_add_piece->Print();

   //cov->Print();
  //cov_sk_all->Print();
  //cov_sk_fhc->Print();
  //cov_sk_rhc->Print();

  //cov_sk_all_final->Print();
  //cov_sk_fhc_final->Print();
  //cov_sk_rhc_final->Print();

  
  return;
}
