/*
  KDI decided to give us a TMatrixDSym including an SK energy scale error.
  This macro makes a new TMatrixDSym, without the SK energy scale row.
*/

//#define DEBUG_MATRIX

void take_row_from_tmatrixdsym(
			       const char * kdifilename = "inputs/official/2013/skdetfsi/skdet_mat_numuoa_2013a_direct_v5_run1to4_final.root",
			       const char * outfilename = "inputs/official/2013/skdetfsi/skdet_mat_numuoa_2013a_direct_v5_run1to4_final_valor.root")
{

  TFile * f_kdi = new TFile(kdifilename,"read");

  //get kdi matrix
  TMatrixDSym * m_kdi_temp = 0;
  f_kdi->GetObject("covmat_sk_fsi",m_kdi_temp);
  if(!m_kdi_temp) {
    cerr << "KDI matrix not found" << endl;
    exit(1);
  }
  TMatrixDSym m_kdi(*m_kdi_temp);
  cout << "KDI matrix:" << endl;
  m_kdi.Print();

  //do the addition
  const int N = m_kdi.GetNcols()-1;
  TMatrixDSym m_out(N);
  for(int i = 0; i < N; i++ ) {
    for(int j = 0; j < N; j++ ) {
      m_out(i,j) = m_kdi(i,j);
    }//j
  }//i
  cout << "Our VaLOR matrix:" << endl;
  m_out.Print();

  //save
  TFile * fout = new TFile(outfilename,"recreate");  
  m_out.Write("covmat_sk_fsi");
  fout->Close();

  //check inverse
  double det;
  TMatrixD covInv(m_out);
  for (int i = 0; i < 6; i++){
    for (int j = 0; j < 6; j++){
      covInv(i,j) = m_out(i,j);
    }
  }
  covInv.Invert(&det);
  cout << "Determinant = " << det << endl;
  cout << "Inverse covariance matrix:" << endl;
  covInv.Print();

  //check product
  TMatrixD prod = m_out*covInv;
  cout << "Product:" << endl;
  prod.Print();


  TH2D* h = new TH2D("h","SKDET + KDI+SI covariance matrix",6,0,6,6,0,6);
  for (int i = 0; i<6; i++){
    for (int j = 0; j<6; j++){
      h->SetBinContent(i+1,j+1,m_out(i,j));
    }
  }

  TCanvas * c = new TCanvas("c","",20,20,500,500);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  h->GetXaxis()->SetTitle("Parameter");
  h->GetYaxis()->SetTitle("Parameter");
  h->Draw("text");
  c->SaveAs("SKDetFSI.pdf");

}

