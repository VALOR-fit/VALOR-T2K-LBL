/*
  make the SKDetector matrix, using numbers supplied from Kameda-san in email
  dated 29 May 2013
*/

void make_skdet() {

  TFile * f_out = new TFile("SKDetErrors2013.root","RECREATE");
  const char * name = "SKDetCovariance";
  TH2D * h = new TH2D(name,name, 6,0,6, 6,0,6);

  //input from Kaemda-san
  double diag[6] = {0.0221,0.0362,0.0404,0.0848,1.0000,0.4790};
  double cov[6][6] = {
    { 1.000, 0.820, 0.620, 0.160, 0.000, 0.028 },
    { 0.820, 1.000, 0.760, 0.100, 0.000, 0.017 },
    { 0.620, 0.760, 1.000, 0.090, 0.000, 0.015 },
    { 0.160, 0.100, 0.090, 1.000, 0.000, 0.007 },
    { 0.000, 0.000, 0.000, 0.000, 1.000, 0.000 },
    { 0.028, 0.017, 0.015, 0.007, 0.000, 1.000 }
  }; 

  //set the entries
  for(int i = 1; i <= 6; i++) {
    for(int j = 1; j <= 6; j++) {
      double element = diag[i-1]*diag[j-1]*cov[i-1][j-1];
      h->SetBinContent(i,j, element);
      cout << cov[i-1][j-1] << " * " << diag[i-1] << " * " << diag[j-1]
	   << " =\t" << element << endl;
    }//j
  }//i

  TCanvas * c = new TCanvas();
  h->Draw("text");
  c->SaveAs("SKDet.pdf");

  //h->SetDirectory(0);
  h->Write();
  f_out->Close();
}
