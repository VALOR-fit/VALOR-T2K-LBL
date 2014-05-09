/*
take a covariance matrix & scale it -> reduced errors for HK studies

root2k -l -b -q '$T2KOSCANA/macros/misc/test_fitter_inputs/scale_down_covariance_matrix.C(2,"inputs/official/2013/banff/postfit_banff_2013_joint_v5.root","inputs/hk/banff/postfit_banff_2013_joint_v5_scale2.root","postfit_cov")'

root2k -l -b -q '$T2KOSCANA/macros/misc/test_fitter_inputs/scale_down_covariance_matrix.C(10,"inputs/official/2013/skdetfsi/skfsipn_2013_joint_Run4.root","inputs/hk/skdetfsi/skfsipn_2013_joint_Run4_scale10.root","skfsipn_cov")'

*/

//#define DEBUG_MATRIX

void scale_down_covariance_matrix
(
 double       scale       = 2,
 const char * infilename  = "inputs/official/2013/banff/postfit_banff_2013_disappearance_v5.root",
 const char * outfilename = "scaled.root",
 const char * matrixname  = "postfit_cov"
)
{
  using namespace valor;

  //get matrix
  TFile * f = new TFile(infilename,"read");
  TMatrixDSym * mtemp = 0;
  f->GetObject(matrixname,mtemp);
  if(!mtemp) {
    cerr << "Matrix not found" << endl;
    exit(1);
  }
  TMatrixDSym m(*mtemp);
  cout << "matrix:" << endl;
  m.Print();
  const int N = m.GetNcols();

  //convert to correlation matrix
  TMatrixDSym m_corr(N);
  for(int i = 0; i < N; i++ ) {
    for(int j = 0; j < N; j++ ) {
      m_corr(i,j) = m(i,j) / (TMath::Sqrt(m(i,i)) * TMath::Sqrt(m(j,j)));
    }//j
  }//i
  cout << "Correlation matrix:" << endl;
  m_corr.Print();

  //diag matrix with scaled values
  TMatrixDSym m_diag(N);
  for(int i = 0; i < N; i++ ) {
    for(int j = 0; j < N; j++ ) {
      if(i != j)
	m_diag(i,j) = 0;
      else {
	m_diag(i,j) = (TMath::Sqrt(m(i,j)) / scale) * (TMath::Sqrt(m(i,j)) / scale);
	cout << m(i,j) << "\t" << (TMath::Sqrt(m(i,j))) << "\t" << (TMath::Sqrt(m(i,j)) / scale) << endl;
      }
    }//j
  }//i
  cout << "Diagonal matrix:" << endl;
  m_diag.Print();

  //add the new diagonal to the correlation matrix
  TMatrixDSym m_out(N);
  for(int i = 0; i < N; i++ ) {
    for(int j = 0; j < N; j++ ) {
      m_out(i,j) = m_corr(i,j) * TMath::Sqrt(m_diag(i,i)) * TMath::Sqrt(m_diag(j,j));
    }//j
  }//i
  cout << "Matrix scaled by " << scale << ":" << endl;
  m_out.Print();

  //save
  TFile * fout = new TFile(outfilename,"recreate");  
  m_out.Write(matrixname);
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

  //check its decomposable
  TMatrixD L = util::CholeskyDecomposition(m_out);
  cout << "Cholesky decomposition: " << endl;
  L.Print();

}

