//
// Compare an identically-named 1-D histogram from 2 different ROOT files.
// Ratio shown is: histogram in 2nd file / histogram in 1st file
//

void compare_1dhist(
  const char * filename1,           /* name of 1st file */
  const char * filename2,           /* name of 2nd file */
  const char * histoname,           /* histogram name   */
  bool zeroOverZeroIsZero = false,  /* use the built in divide where 0/0=0, or use 0/0=1 */
  bool draw_only_ratio = true       /* draw only the ratio? or also plot the 2 histograms overlaid */
)
{
  TFile * f1 = new TFile(filename1,"read");
  TFile * f2 = new TFile(filename2,"read");

  TH1D * h1 = (TH1D*) f1->Get(histoname);
  TH1D * h2 = (TH1D*) f2->Get(histoname);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);

  TH1D * ratio = new TH1D(*h2);
  if (zeroOverZeroIsZero)
    ratio->Divide(h1);
  else {
    Int_t binx;
    Double_t c0,c1,w;
    for (binx=0;binx<=h1->GetNbinsX()+1;binx++) {
      c0  = h2->GetBinContent(binx);
      c1  = h1->GetBinContent(binx);
      if (c1) w = c0/c1;
      else if(!c1 && !c0) w=1;
      else    w = 0;
      ratio->SetBinContent(binx,w);
    }
  }
      
  TCanvas * c = new TCanvas("c","",10,10,600,600);
  if(!draw_only_ratio) {
    c->Divide(1,2);
    c->cd(1);
    h1->Draw("p");
    h2->Draw("same");
    c->cd(2);
  }
  ratio->Draw();
  c->Update();
}
