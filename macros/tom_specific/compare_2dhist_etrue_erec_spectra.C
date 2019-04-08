//
// Compare an identically-named 2-D histogram from 2 different ROOT files.
// Ratio shown is: histogram in 2nd file / histogram in 1st file
//

#include "/home/hyperk/tdealtry/git-linux/codesnippets/print_out_axis.C"

void compare_2dhist(
  const char * filename1,           /* name of 1st file */
  const char * filename2           /* name of 2nd file */
  //const char * histoname            /* histogram name   */
)

{
  TFile * f1 = new TFile(filename1,"read");
  TFile * f2 = new TFile(filename2,"read");
  f1->ls();
  TH2D * h1, *h2;
  TCanvas * c = new TCanvas("c","",10,10,600,600);
  c->SetRightMargin(0.18);
  c->SaveAs("compare2d.pdf[");
  string samples[4] = {"1Rmu","1Re","RHC1Rmu","RHC1Re"};
  string modes[4] = {"CCQE","CCnQE","NC","TOTAL"};
  TLatex latex;
  //latex.SetTextSize(0.01);
  for(int is = 0; is < 4; is++) {
    for(int im = 0; im < 4; im++) {
      const char * histoname = TString::Format("enu_erec_%s_%s", samples[is].c_str(), modes[im].c_str());
      cout << histoname << endl;
      f1->GetObject(histoname, h1);
      f2->GetObject(histoname, h2);

      h1->SetMarkerStyle(20);
      h1->SetMarkerSize(1);
      
      cout << fixed << setprecision(3) << endl;
      cout << "INTEGRAL 1 " << h1->Integral() << endl;
      cout << "INTEGRAL 2 " << h2->Integral() << endl;

      cout << endl << "X" << endl;
      print_out_axis(h1->GetXaxis());
      print_out_axis(h2->GetXaxis());
      cout << endl << "Y" << endl;
      print_out_axis(h1->GetYaxis());
      print_out_axis(h2->GetYaxis());

      h1->Add(h2, -1);
      h1->GetXaxis()->SetRangeUser(0,4.9);
      h1->GetYaxis()->SetRangeUser(0,4.9);
      double zmax = TMath::Abs(TMath::Max(h1->GetMaximum(), h1->GetMinimum()));
      h1->GetZaxis()->SetRangeUser(-zmax,+zmax);
      h1->Draw("COLZ");

      latex.DrawLatexNDC(0, 0.01, histoname);

      //TH2D * ratio = new TH2D(*h2);
      //ratio->Divide(h1);
      
      //ratio->Draw("colz");
      //gPad->SetLogz();
      c->Update();
      c->SaveAs("compare2d.pdf");

      h1->Add(h2);
      //h1->Divide(h2);
      for(int ix = 1; ix <= h1->GetNbinsX(); ix++) {
	for(int iy = 1; iy <= h1->GetNbinsX(); iy++) {
	  double v1 = h1->GetBinContent(ix, iy);
	  double v2 = h2->GetBinContent(ix, iy);
	  double ratio;
	  if(TMath::Abs(v1) < 1E-16 && TMath::Abs(v2) < 1E-16)
	    ratio = 0;
	  else
	    ratio = (v1 - v2) / (v1 + v2);
	  h1->SetBinContent(ix, iy, ratio);
	}//iy
      }//ix

      zmax = h1->GetMaximum();
      double zmin = h1->GetMinimum();
      h1->GetZaxis()->SetRangeUser(-zmax,+zmax);
      h1->Draw("COLZ");
      latex.DrawLatexNDC(0, 0.01, histoname);
      c->Update();
      c->SaveAs("compare2d.pdf");
    }//im
  }//is
  c->SaveAs("compare2d.pdf]");
}
