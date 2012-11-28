//
// Compare an identically-named 2-D histogram from 2 different ROOT files.
// Ratio shown is: histogram in 2nd file / histogram in 1st file
//

void compare_2dhist(
  const char * filename1,           /* name of 1st file */
  const char * filename2,           /* name of 2nd file */
  const char * histoname            /* histogram name   */
)

{
  TFile * f1 = new TFile(filename1,"read");
  TFile * f2 = new TFile(filename2,"read");

  TH2D * h1 = (TH2D*) f1->Get(histoname);
  TH2D * h2 = (TH2D*) f2->Get(histoname);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);

  TH2D * ratio = new TH2D(*h2);
  ratio->Divide(h1);
      
  TCanvas * c = new TCanvas("c","",10,10,600,600);
  ratio->Draw("colz");
  gPad->SetLogz();
  c->Update();
}
