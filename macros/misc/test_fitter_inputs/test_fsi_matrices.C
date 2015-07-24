void test_fsi_matrices() {
  TFile f1("inputs/official/2013/fsi/FSIandSIErrors2013.root");
  TFile f2("inputs/official/2012/fsi/FSIandSIError_ICL_SKDetectorBinning.root");

  TH2D * h = 0;
  f1.GetObject("FSICovariance",h);
  if(!h) {f1.ls(); exit(1);}
  TMatrixTSym<double> * m = 0;
  f2.GetObject("fsi_cov",m);
  if(!m) {f2.ls(); exit(1);}
  TH2D * hcompare = new TH2D(*h);
  
  for(int ix = 1; ix <= h->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= h->GetNbinsX(); iy++) {
      double hcurr = h->GetBinContent(ix,iy);
      double mcurr = m(ix-1,iy-1);
      double comparecurr = (hcurr-mcurr)*2/(hcurr+mcurr);
      cout << hcurr << "\t\t" << mcurr << "\t\t\t" << comparecurr << endl;
      hcompare->SetBinContent(ix,iy, comparecurr);
    }//iy
  }//ix
  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.1);
  hcompare->Draw("colz");
  c->SaveAs("compare.pdf");
}
