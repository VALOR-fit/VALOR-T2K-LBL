//
// Compare all identically-named 1-D histogram from 2 different ROOT files.
// Ratio shown is: histogram in 2nd file / histogram in 1st file
// If drawn, the original histograms are print with the points from file1,
//  and the line from file2
//
// Tom Dealtry

void compare_1dhist_all(
  const char * filename1,           /* name of 1st file */
  const char * filename2,           /* name of 2nd file */
  const char * tag = "",            // tag for output filename
  string checkstr = "",             /* only check histograms that contain this string in their name */
  string skipstr  = "XXXXXXX",             /* skip checks for histograms with this in their name */
  double confidence = 1.0,          /* if a Kolmogorov test returns a value greater than this, assume compatible & don't draw */
  double ratiotest = 1E-6,          /* maximum ratio needs to be more than this (away from 1) */
  bool zeroOverZeroIsZero = false,  /* use the built in divide where 0/0=0, or use 0/0=1 */
  bool draw_only_ratio = false,     /* draw only the ratio? or also plot the 2 histograms overlaid */
  double xmin = -9999999,           /* minimum of x axis to draw */
  double xmax = +29.99,           /* maximum of x axis to draw */
  bool save_canvas = true,         /* save canvases to a .pdf? Also turns on batch mode*/
  bool dividingTwoRatios = false,   /* special option to use to compare 2 ratio histograms. Needs zerOverZeroIsZero = false */
  bool drawSubtraction = false      /* compare via c2-c1, rather than division */
)
{
  if(save_canvas)
    gROOT->SetBatch(true);
  char *s = new char[50];
  int i = 0;
  TString histoname;

  TFile * f1 = new TFile(filename1,"read");
  TFile * f2 = new TFile(filename2,"read");

  const double testaround = drawSubtraction ? 0 : 1;
  
  //need to make this, otherwise the real canvas isn't made in the loop
  TCanvas * ct = new TCanvas("ct","",10,10,600,600);
  if(!draw_only_ratio)
    ct->Divide(1,2);
  if(save_canvas)
    ct->Print(TString::Format("comparisons_%s.pdf[", tag));

  //cycles through histogram names in file 1
  TKey *key;
  TIter next(f1->GetListOfKeys());
  while ((key = (TKey *) next())) {
    histoname = key->GetName();
    //if not skipping, draw histograms
    if(i<=0) {
      cout << "\nViewing comparison for histogram: " << histoname.Data() << endl;
      if(!histoname.Contains(checkstr.c_str()) || histoname.Contains(skipstr.c_str())) {
	cout << "Skipping due to histogram name" << endl;
	continue;
      }
      TH1D * h1 = (TH1D*) f1->Get(histoname.Data());
      TH1D * h2 = (TH1D*) f2->Get(histoname.Data());
      if(!h1||!h2) {cout << "   histogram not found, continuing" << endl; continue;}
      //run kolmoogov test
      double test =  h1->KolmogorovTest(h2);
      cout << "kolmogov test returns: " << test;
      if(test>confidence) {
	cout << ", skipping drawing" << endl;
	continue;
      }
      cout << ", dividing histograms" << endl;
      //also get the ratio, with 0/0 = 0 or 0/0 = 1 given in function argument
      TH1D * ratio = new TH1D(*h2);
      ratio->GetYaxis()->SetTitle("Division");
      if (zeroOverZeroIsZero)
	ratio->Divide(h1);
      else if(dividingTwoRatios) {
	Int_t binx;
	Double_t c0,c1,w;
	for (binx=0;binx<=h1->GetNbinsX()+1;binx++) {
	  c0  = TMath::Abs(h2->GetBinContent(binx) - 1);
	  c1  = TMath::Abs(h1->GetBinContent(binx) - 1);
	  if (c1) w = c0/c1;
	  else if(!c1 && !c0) w=1;
	  else    w = 0;
	  ratio->SetBinContent(binx,w);
	}
      }
      else if(drawSubtraction) {
	ratio->Add(h1, -1);
	ratio->GetYaxis()->SetTitle("Subtraction");
      }
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
	  if(TMath::Abs(w-1) > 0.000000000001)
	    cout << " Point at: " << c1 << " Line at: " << c0 << " Ratio is: " << w 
		 << " In bin: " << h1->GetBinLowEdge(binx) << "-" << h1->GetBinLowEdge(binx+1) << endl;
	}
      } //zeroOverZeroIsZero?
      //Another check, to see if the ratio is small (e.g. from the migration of a single value from one bin to another)
      if((ratio->GetBinContent(ratio->GetMaximumBin()) < testaround + ratiotest)
	 && (ratio->GetBinContent(ratio->GetMinimumBin()) > testaround - ratiotest)) {
	cout << "skipping due to maximum/minimum ratio test" << endl;
	cout << "Max " << ratio->GetBinContent(ratio->GetMaximumBin()) << "\tMin " << ratio->GetBinContent(ratio->GetMinimumBin()) << endl;
	continue;
      }
      //draw
      TCanvas * c = new TCanvas("c","",10,10,600,600);
      TLegend * leg = new TLegend(0.0,0.0,0.8,0.1);
      leg->SetFillStyle(0);
      leg->AddEntry((TObject*)0,histoname.Data(),"");
      if(!draw_only_ratio) {
	c->Divide(1,2);
	c->cd(1);
	h1->GetXaxis()->SetRangeUser(xmin, xmax);
	h1->SetMarkerSize(2);
	h1->SetMarkerColor(kBlack);
	h2->GetXaxis()->SetRangeUser(xmin, xmax);
	h2->SetLineColor(kRed);
	h1->Draw("p");
	h2->Draw("same");
	c->cd(2);
      }
      ratio->GetXaxis()->SetRangeUser(xmin, xmax);
      ratio->Draw();
      c->Update();
      if(save_canvas) {
	if(!draw_only_ratio) {
	  ct->cd(1);
	  h1->Draw("p");
	  h2->Draw("same");
	  ct->cd(2);
	}
	else
	  ct->cd();
	ratio->Draw();
	ct->Update();
	leg->Draw();
	ct->Print(TString::Format("comparisons_%s.pdf", tag));
      }

      if(!save_canvas)
	gets(s); //waits for a [RETURN] entry on the command line.
      else
	s = "0";
      //Can type a number to skip n histograms
      i = atoi(s);
      if(i>0)
	cout << "skipping " << i << " histograms" << endl;
      //cleanup
      delete c;
      delete leg;
      delete ratio;
    } // i<=0
    else
      i-=1;
  }//while there are TKeys
  if(save_canvas)
    ct->Print(TString::Format("comparisons_%s.pdf]", tag));
}
