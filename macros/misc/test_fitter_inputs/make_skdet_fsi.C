/*
  Make the SK detector + FSI/SI matrix for the numu disapp analysis
  Adapted from Lorena's macro
 */

//#define DEBUG_MATRIX

TMatrixD * TH2D2TMatrixD(TH2D * h);

void make_skdet_fsi(
		    const char * fsifilename = "inputs/official/2013/skdetfsi/FSIandSIErrors2013.root",
		    const char * skdfilename = "inputs/official/2013/skdetfsi/SKDetErrors2013.root",
		    const char * outfilename = "SKDetFSISIErrors2013.root"
){

  TFile * f_fsi = new TFile(fsifilename,"read");
  TFile * f_skd = new TFile(skdfilename,"read");

  //get fsi matrix
  TH2D * h_fsi = 0;
  f_fsi->GetObject("FSICovariance",h_fsi);
  if(!h_fsi) {
    cerr << "FSI matrix not found" << endl;
    exit(1);
  }
  TMatrixD m_fsi = TH2D2TMatrixD(h_fsi);
  cout << "FSI matrix:" << endl;
  m_fsi.Print();

  //get skdet matrix
  TH2D * h_skd = 0;
  f_skd->GetObject("SKDetCovariance",h_skd);
  if(!h_skd) {
    cerr << "SKDet matrix not found" << endl;
    exit(1);
  }
  TMatrixD m_skd = TH2D2TMatrixD(h_skd);
  cout << "SKDet matrix:" << endl;
  m_skd.Print();

  //do the addition
  TMatrixD m_out(m_fsi);
  for(int i = 0; i < 6; i++ ) {
    for(int j = 0; j < 6; j++ ) {
      m_out(i,j) = m_skd(i,j) + m_fsi(i,j);
    }//j
  }//i
  cout << "SKDetFSI matrix:" << endl;
  m_out.Print();

  //save
  TFile * fout = new TFile(outfilename,"recreate");  
  m_out.Write("skdet_fsi_cov");
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


  TH2D* h = new TH2D("h","SKDET + FSI+SI covariance matrix",6,0,6,6,0,6);
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





TMatrixD TH2D2TMatrixD(TH2D * h)
{

  const int nrows = h->GetNbinsX();
  const int ncols = h->GetNbinsY();

  TMatrixD matrix(nrows,ncols);
  for(int i = 1; i <= nrows; i++) {
    for(int j = 1; j <= ncols; j++) {
      matrix(i-1,j-1)=h->GetBinContent(i,j);
#ifdef DEBUG_MATRIX
      cout << i << "\t" << j << "\t" << h->GetBinContent(i,j) << endl;
#endif
    }
  }
#ifdef DEBUG_MATRIX
  matrix.Print();
#endif

  return matrix;
}
