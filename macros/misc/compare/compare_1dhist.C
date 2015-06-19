//
// Compare an identically-named 1-D histogram from 2 different ROOT files.
// Ratio shown is: histogram in 2nd file / histogram in 1st file
//

void compare_1dhist
(
 const char * filename1,           /* name of 1st file */
 const char * filename2,           /* name of 2nd file */
 string histoname1,                /* histogram name   */
 string histoname2 = "",           /* histogram name (if "", assume same as 1st)  */
 bool zeroOverZeroIsZero = false,  /* use the built in divide where 0/0=0, or use 0/0=1 */
 bool draw_only_ratio = false,     /* draw only the ratio? or also plot the 2 histograms overlaid */
 string tag                        /*if not "", will save in tag.pdf*/
 )
{
  TFile * f1 = new TFile(filename1,"read");
  TFile * f2 = new TFile(filename2,"read");

  if(histoname2 == "")
    histoname2 = histoname1;
  TH1D * h1 = (TH1D*) f1->Get(histoname1.c_str());
  TH1D * h2 = (TH1D*) f2->Get(histoname2.c_str());

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
  if(tag != "")
    c->SaveAs(Form("%s.pdf",tag.c_str()));
}
