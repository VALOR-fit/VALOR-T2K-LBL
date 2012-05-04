//
// Compare all identically-named 1-D histogram from 2 different ROOT files.
// Ratio shown is: histogram in 2nd file / histogram in 1st file
//

void compare_1dhist_all(
  const char * filename1,           /* name of 1st file */
  const char * filename2,           /* name of 2nd file */
  double confidence = 1.0,          /* if a Kolmogorov test returns a value greater than this, assume compatible & don't draw */
  bool zeroOverZeroIsZero = false,  /* use the built in divide where 0/0=0, or use 0/0=1 */
  bool draw_only_ratio = true       /* draw only the ratio? or also plot the 2 histograms overlaid */
)
{
  char *s = new char[5];
  int i = 0;
  TString histoname;

  TFile * f1 = new TFile(filename1,"read");
  TFile * f2 = new TFile(filename2,"read");

  //cycles through histogram names in file 1
  TKey *key;
  TIter next(f1->GetListOfKeys());
  while ((key = (TKey *) next())) {
    histoname = key->GetName();
    //if not skipping, draw histograms
    if(i<=0) {
      cout << "\nViewing comparison for histogram: " << histoname.Data() << endl;
      TH1D * h1 = (TH1D*) f1->Get(histoname.Data());
      TH1D * h2 = (TH1D*) f2->Get(histoname.Data());
      if(!h1||!h2) {cout << "   histogram not found, continuing" << endl; continue;}
      //run kolmoogov test
      double test =  h1->KolmogorovTest(h2);
      cout << "kolmogov test returns: " << test;
      if(test>confidence)
	cout << ", skipping drawing" << endl;
      else {
      	cout << ", drawing histograms" << endl;
	//also get the ratio, with 0/0 = 0 or 0/0 = 1 given in function argument
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
	    if(w>1000000 || w<-1000000)
	      cout << "Point at: " << c1 << " Line at: " << c0 << " Ratio is: " << w << endl;
	  }
	} //zeroOverZeroIsZero?
	//draw
	TCanvas * c = new TCanvas("c","",10,10,600,600);
	c->Divide(1,2);
	c->cd(1);
	h1->Draw("p");
	h2->Draw("same");
	c->cd(2);
	ratio->Draw();
	c->Update();

	gets(s); //waits for a [RETURN] entry on the command line.
	//Can type a number to skip n histograms
	i = atoi(s);
	if(i>0)
	  cout << "skipping " << i << " histograms" << endl;
      }//passes Kolmorogov test?
    } // i<=0
    else
      i-=1;
  }//while there are TKeys
}
