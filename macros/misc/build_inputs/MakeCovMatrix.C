
//
// ------------------------------------------------------------------------------------------------
// Produce the post-Banff covariance matrix (flux+xsec) with the right format for the fitter
//
// Davide Sgalaberna
// ------------------------------------------------------------------------------------------------
//


#include <iomanip.h>

//

#define SINGLEDIAGMATR // build a single diagonal covariance matrix

#ifndef SINGLEDIAGMATR
#define ADDCOVMATR
#define UNITMATRIX
#endif

void MakeCovMatrix( 
		   
		   const char *filename = "$T2KOSCANA/inputs/official/2014/flux/flux_covariance_13av1_prelim.root",
		   const char *covname = "total_flux_cov",
		   
		   //const char *outfilename = "$T2KOSCANA/outputs/test/cov_flux_prelim_13av1_xsec_postbanff2013.root",
		   //const char *outfilename = "$T2KOSCANA/outputs/test/cov_flux_prelim_13av1_xsec_prelim.root",		   
		   const char *outfilename = "$T2KOSCANA/outputs/test/cov_flux_prelim_diag_2015_skdetfsi.root",		   

		   const bool CorrMode = true, // true -> correlated FHC-RHC, false -> not correlated

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

		   // Add a piece of covariance matrix 1
		   // const char *addfilename = "$T2KOSCANA/inputs/official/2013/banff/postfit_banff_2013_joint_v5.root",
		   // const char *addcovname = "postfit_cov",
		   // const int first_add = 25, // rhc-only 
		   // const int last_add = 32,

		   const char *addfilename = "$T2KOSCANA/inputs/official/2014/xsec/xsec_covariance_2014_v1.root",
		   const char *addcovname = "xsec_cov",
		   const int first_add = 0, // rhc-only
		   const int last_add = 13,
		   
		   //
		   
		   const int NDiagMatr_1 = 12,  
		   const int NDiagMatr_2 = 12,
		    ){
  
  using namespace std;
  using namespace valor;
  
#ifdef SINGLEDIAGMATR

  const int NDiagMatr_all = NDiagMatr_1 + NDiagMatr_2;

  cout << endl;
  cout << "Build a unique diagonal covariance matrix" << endl
       << "Matrix 1 " << NDiagMatr_1 << " x " << NDiagMatr_2 << endl
       << "Matrix 2 " << NDiagMatr_1 << " x " << NDiagMatr_2 << endl
       << "1+2 matrix " << NDiagMatr_all << " x " << NDiagMatr_all
       << endl
       << endl;

  TVectorD *par_diag_1 = new TVectorD(NDiagMatr_1);
  TMatrixD *covdiag_1 = new TMatrixD(NDiagMatr_1,NDiagMatr_1);
  for(int irow=0;irow<covdiag_1->GetNrows();irow++){
    for(int icol=0;icol<covdiag_1->GetNrows();icol++){
      if(irow==icol) covdiag_1(irow,icol) = 0.1;
      else covdiag_1(irow,icol) = 0.;
    }
  }

  TVectorD *par_diag_2 = new TVectorD(NDiagMatr_2);
  TMatrixD *covdiag_2 = new TMatrixD(NDiagMatr_2,NDiagMatr_2);
  for(int irow=0;irow<covdiag_2->GetNrows();irow++){
    for(int icol=0;icol<covdiag_2->GetNrows();icol++){
      if(irow==icol) covdiag_2(irow,icol) = 0.1;
      else covdiag_2(irow,icol) = 0.;
    }
  }

  TVectorD *par_diag_all = new TVectorD(NDiagMatr_all);
  TMatrixD *covdiag_all = new TMatrixD(NDiagMatr_all,NDiagMatr_all);
  for(int irow=0;irow<covdiag_all->GetNrows();irow++){
    for(int icol=0;icol<covdiag_all->GetNrows();icol++){
      if(irow==icol) covdiag_all(irow,icol) = 0.1;
      else covdiag_all(irow,icol) = 0.;
    }
  }

  // Write output covariance matrices in a file
  TFile *outfile = new TFile(outfilename,"RECREATE");
  covdiag_1->Write("cov_diag_app");
  covdiag_2->Write("cov_diag_dis");
  covdiag_all->Write("cov_diag_all");
  //
  par_diag_all->Write("par_diag_all");
  par_diag_1->Write("par_diag_app");
  par_diag_2->Write("par_diag_dis");
  outfile->Close();

  TCanvas *c0 = new TCanvas();
  covdiag_1->Draw("colz");
  TCanvas *c1 = new TCanvas();
  covdiag_2->Draw("colz");
  TCanvas *c2 = new TCanvas();
  covdiag_all->Draw("colz");

  return;

#endif




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


  
  //if(!CorrMode){
  //cout << "Remove correlations" << endl;
  //for(int irow=0;irow<cov_sk_all){
  //for(){
  //}
  //}
  //}
  

  // Add matrix 
#ifdef ADDCOVMATR
  
  cout << endl;
  cout << "///////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "//                                                                           " << endl;
  cout << "// Add piece of covariance matrix w/o correlations from " << addfilename << " " << endl;
  cout << "//                                                                           " << endl;
  cout << "///////////////////////////////////////////////////////////////////////////////" << endl;
  
  TFile *addfilein = new TFile(addfilename,"READ");
  TMatrixD *cov_sk_add = (TMatrixD*) addfilein->Get(addcovname);
  
  //Take parameters
  TVectorD *par_sk_add = (TVectorD*) addfilein->Get("xsec_param_prior");
  TObjArray *parname_sk_add = (TObjArray*) addfilein->Get("paramNames");
  if(parname_sk_add!=0){  
    parname_sk_add->SetOwner(kTRUE);
    int npar = parname_sk_add->GetEntriesFast();
    cout << endl;
    cout << "# of parameters = " << npar << endl;
    cout << endl;
    for(int ipar=0;ipar<npar;ipar++){
      TObjString *parname_objstring = (TObjString*) parname_sk_add->At(ipar);
      cout << parname_objstring->GetName() << endl;
    }
    cout << endl;
  }
  
  addfilein->Close();

  TMatrixD *cov_sk_add_piece = (TMatrixD*) cov_sk_add->GetSub(first_add,last_add,first_add,last_add);
  


  // // Check additional covariance matrix
  // for(int irow=0;irow<cov_sk_add_piece->GetNrows();irow++){
  //   for(int icol=0;icol<cov_sk_add_piece->GetNrows();icol++){
  //     if(cov_sk_add_piece(irow,icol)==1.) cov_sk_add_piece(irow,icol) = 0.1;
  //   }
  // }
  



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


  // Write output parameter vector 

  const int NdimAll = cov_sk_all_final->GetNrows();
  TVectorD *par_sk_all_final = new TVectorD(NdimAll);
  const int NdimFHC = cov_sk_fhc_final->GetNrows();
  TVectorD *par_sk_fhc_final = new TVectorD(NdimFHC);
  const int NdimRHC = cov_sk_rhc_final->GetNrows();
  TVectorD *par_sk_rhc_final = new TVectorD(NdimRHC);

  // Add original parameter values
  const int NdimOrig = cov_sk_all->GetNrows(); 
  for(int ipar=0;ipar<NdimOrig;ipar++) 
    par_sk_all_final(ipar) = 1.;  //set prior values always to 1 for simplicity
#ifdef ADDCOVMATR
  // Add adding parameter values
  const int NdimAdd = cov_sk_add_piece->GetNrows();
  for(int ipar=0;ipar<NdimAdd;ipar++) 
    par_sk_all_final(ipar+NdimOrig) = 1.; //set prior values always to 1 for simplicity
#endif
  
  // Add original parameter values
  const int NdimOrig = cov_sk_fhc->GetNrows(); 
  for(int ipar=0;ipar<NdimOrig;ipar++) 
    par_sk_fhc_final(ipar) = 1.;  //set prior values always to 1 for simplicity
#ifdef ADDCOVMATR
  // Add adding parameter values
  const int NdimAdd = cov_sk_add_piece->GetNrows();
  for(int ipar=0;ipar<NdimAdd;ipar++) 
    par_sk_fhc_final(ipar+NdimOrig) = 1.; //set prior values always to 1 for simplicity
#endif

  // Add original parameter values
  const int NdimOrig = cov_sk_rhc->GetNrows(); 
  for(int ipar=0;ipar<NdimOrig;ipar++) 
    par_sk_rhc_final(ipar) = 1.;  //set prior values always to 1 for simplicity
#ifdef ADDCOVMATR
  // Add adding parameter values
  const int NdimAdd = cov_sk_add_piece->GetNrows();
  for(int ipar=0;ipar<NdimAdd;ipar++) 
    par_sk_rhc_final(ipar+NdimOrig) = 1.; //set prior values always to 1 for simplicity
#endif

  //

#ifdef UNITMATRIX
  for(int irow=0;irow<cov_sk_all_final->GetNrows();irow++){
    for(int icol=0;icol<cov_sk_all_final->GetNrows();icol++){
      if(irow>=50){
	if(irow==icol) cov_sk_all_final(irow,icol) = 0.1;
	else cov_sk_all_final(irow,icol) = 0.;
      }
    }
  }
  for(int irow=0;irow<cov_sk_fhc_final->GetNrows();irow++){
    for(int icol=0;icol<cov_sk_fhc_final->GetNrows();icol++){
      if(irow>=25){
      if(irow==icol) cov_sk_fhc_final(irow,icol) = 0.1;
      else cov_sk_fhc_final(irow,icol) = 0.;
      }
    }
  }
  for(int irow=0;irow<cov_sk_rhc_final->GetNrows();irow++){
    for(int icol=0;icol<cov_sk_rhc_final->GetNrows();icol++){
      if(irow>=25){
	if(irow==icol) cov_sk_rhc_final(irow,icol) = 0.1;
	else cov_sk_rhc_final(irow,icol) = 0.;
      }
    }
  }
#endif
  
  // Write output covariance matrices in a file
  TFile *outfile = new TFile(outfilename,"RECREATE");
  cov_sk_all_final->Write("prelim_cov_all");
  cov_sk_fhc_final->Write("prelim_cov_fhc");
  cov_sk_rhc_final->Write("prelim_cov_rhc");
  //
  par_sk_all_final->Write("prelim_par_all");
  par_sk_fhc_final->Write("prelim_par_fhc");
  par_sk_rhc_final->Write("prelim_par_rhc");
  outfile->Close();


#ifdef ADDCOVMATR
  TCanvas *ccov_sk_add_piece = new TCanvas("ccov_sk_add_piece","ccov_sk_add_piece");
  cov_sk_add_piece->Draw("colz text");
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
  

  // Print out values

  cout << endl;
  cout << "# of selected parameters = " << endl;
  cout << "  - All      --> " << cov_sk_all_final->GetNrows() << endl;
  cout << "  - FHC Only --> " << cov_sk_fhc_final->GetNrows() << endl;
  cout << "  - RHC Only --> " << cov_sk_rhc_final->GetNrows() << endl;
  cout << endl;

  //cov_sk_add->Print();
  //cov_sk_add_piece->Print();

   //cov->Print();
  //cov_sk_all->Print();
  //cov_sk_fhc->Print();
  //cov_sk_rhc->Print();

  //cov_sk_all_final->Print();
  //cov_sk_fhc_final->Print();
  //cov_sk_rhc_final->Print();

  par_sk_all_final->Print();
  par_sk_fhc_final->Print();
  par_sk_rhc_final->Print();
  
  return;
}
