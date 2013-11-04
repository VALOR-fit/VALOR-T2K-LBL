int find_first_bin_below(double threshold, TH1D * h) {
  for(int ix = 1; ix <= h->GetNbinsX(); ix++) {
    if(h->GetBinContent(ix) < threshold) {
      int bin = TMath::Max(ix - 1, 1);
      //cout << bin << endl;
      return bin;
    }
  }//ix
}

int find_last_bin_below(double threshold, TH1D * h) {
  for(int ix = h->GetNbinsX(); ix >= 1; ix--) {
    if(h->GetBinContent(ix) < threshold)
      return TMath::Min(ix + 1, h->GetNbinsX());
  }//ix
}


void calculate_sensitvity_p_value
(
 int minexpt = 0,
 int maxexpt = 500,
 double dmix68_l = 0.456, //data limits
 double dmix68_h = 0.572,
 double dmix90_l = 0.428,
 double dmix90_h = 0.599,
 double ddm68_l = 0.00235,
 double ddm68_h = 0.00257,
 double ddm90_l = 0.00229,
 double ddm90_h = 0.00264
)
{

  using namespace t2kana;

  //set up input tree
  TFile f("final.root");
  FitParamList *  fparm = 0;
  int iexpt, idm32sq, imix23;
  TTree * fitnt = (TTree*)f.Get("fitnt");
  TH2D * chisq_surface = (TH2D*)f.Get("chisq_surface");
  fitnt->SetBranchAddress("best_fit_params", &fparm);
  fitnt->SetBranchAddress("iexpt",   &iexpt);
  fitnt->SetBranchAddress("idm32sq", &idm32sq);
  fitnt->SetBranchAddress("imix23",  &imix23);

  //book histograms
  const int nexpt = 510;
  TH2D * h [nexpt] = {0};
  TH1D * hm[nexpt] = {0};
  TH1D * hd[nexpt] = {0};
  TH1D * hm68l = new TH1D("","",
			 chisq_surface->GetNbinsX(), chisq_surface->GetXaxis()->GetBinUpEdge(0), chisq_surface->GetXaxis()->GetBinUpEdge(chisq_surface->GetNbinsX()));
  TH1D * hm90l = new TH1D("","",
			 chisq_surface->GetNbinsX(), chisq_surface->GetXaxis()->GetBinUpEdge(0), chisq_surface->GetXaxis()->GetBinUpEdge(chisq_surface->GetNbinsX()));
  TH1D * hd68l = new TH1D("","",
			 chisq_surface->GetNbinsY(), chisq_surface->GetYaxis()->GetBinUpEdge(0), chisq_surface->GetYaxis()->GetBinUpEdge(chisq_surface->GetNbinsY()));
  TH1D * hd90l = new TH1D("","",
			 chisq_surface->GetNbinsY(), chisq_surface->GetYaxis()->GetBinUpEdge(0), chisq_surface->GetYaxis()->GetBinUpEdge(chisq_surface->GetNbinsY()));
  TH1D * hm68h = new TH1D("","",
			 chisq_surface->GetNbinsX(), chisq_surface->GetXaxis()->GetBinUpEdge(0), chisq_surface->GetXaxis()->GetBinUpEdge(chisq_surface->GetNbinsX()));
  TH1D * hm90h = new TH1D("","",
			 chisq_surface->GetNbinsX(), chisq_surface->GetXaxis()->GetBinUpEdge(0), chisq_surface->GetXaxis()->GetBinUpEdge(chisq_surface->GetNbinsX()));
  TH1D * hd68h = new TH1D("","",
			 chisq_surface->GetNbinsY(), chisq_surface->GetYaxis()->GetBinUpEdge(0), chisq_surface->GetYaxis()->GetBinUpEdge(chisq_surface->GetNbinsY()));
  TH1D * hd90h = new TH1D("","",
			 chisq_surface->GetNbinsY(), chisq_surface->GetYaxis()->GetBinUpEdge(0), chisq_surface->GetYaxis()->GetBinUpEdge(chisq_surface->GetNbinsY()));
  for(int ie = 0; ie < nexpt; ie++) {
    h[ie]  = new TH2D("","",
		      chisq_surface->GetNbinsX(), chisq_surface->GetXaxis()->GetBinUpEdge(0), chisq_surface->GetXaxis()->GetBinUpEdge(chisq_surface->GetNbinsX()),
		      chisq_surface->GetNbinsY(), chisq_surface->GetYaxis()->GetBinUpEdge(0), chisq_surface->GetYaxis()->GetBinUpEdge(chisq_surface->GetNbinsY()));
    hm[ie] = new TH1D("","",
		      chisq_surface->GetNbinsX(), chisq_surface->GetXaxis()->GetBinUpEdge(0), chisq_surface->GetXaxis()->GetBinUpEdge(chisq_surface->GetNbinsX()));
    hd[ie] = new TH1D("","",
		      chisq_surface->GetNbinsY(), chisq_surface->GetYaxis()->GetBinUpEdge(0), chisq_surface->GetYaxis()->GetBinUpEdge(chisq_surface->GetNbinsY()));
    for(int ix = 1; ix <= hm[ie]->GetNbinsX(); ix++)
      hm[ie]->SetBinContent(ix, 9999);
    for(int ix = 1; ix <= hd[ie]->GetNbinsX(); ix++)
      hd[ie]->SetBinContent(ix, 9999);
  }//ie

  //fill histograms (2D & 1D)
  cout << fitnt->GetEntries() << " total entries in tree" << endl;
  for(int i = 0; i < fitnt->GetEntries(); i++) {
    fitnt->GetEntry(i);
    if((i % 10000) == 0)
      cout << i << endl;

    double chisq = fparm->GetBestFitChisq();
    h[iexpt]->SetBinContent(imix23, idm32sq, chisq);
    //cout << chisq << "\t" << hm[iexpt]->GetBinContent(imix23) << "\t" << hd[iexpt]->GetBinContent(idm32sq) << endl;
    if(chisq < hm[iexpt]->GetBinContent(imix23))
      hm[iexpt]->SetBinContent(imix23, chisq);
    if(chisq < hd[iexpt]->GetBinContent(idm32sq))
      hd[iexpt]->SetBinContent(idm32sq, chisq);
  }//tree entries

  //chisq -> dchisq
  for(int ie = 0; ie < nexpt; ie++) {
    //double chisq_min = h[ie]->GetMinimum(); //doesn't work when there are holes in the surface
    double chisq_min = 9999;
    for(int ix = 1; ix <= h[ie]->GetNbinsX(); ix++) 
      for(int iy = 1; iy <= h[ie]->GetNbinsY(); iy++) 
	if(h[ie]->GetBinContent(ix,iy) < chisq_min && h[ie]->GetBinContent(ix,iy) > 1)
	  chisq_min = h[ie]->GetBinContent(ix,iy);
    //cout << ie << "\t" << chisq_min << endl;
    for(int ix = 1; ix <= h[ie]->GetNbinsX(); ix++) {
      for(int iy = 1; iy <= h[ie]->GetNbinsY(); iy++)
	h[ie]->SetBinContent(ix,iy, h[ie]->GetBinContent(ix,iy) - chisq_min);
    }
    for(int ix = 1; ix <= hm[ie]->GetNbinsX(); ix++)
      hm[ie]->SetBinContent(ix, hm[ie]->GetBinContent(ix) - chisq_min);
    for(int ix = 1; ix <= hd[ie]->GetNbinsX(); ix++)
      hd[ie]->SetBinContent(ix, hd[ie]->GetBinContent(ix) - chisq_min);
  }//ie

  //find limits
  //need to be careful with the mix limits - some are double peaked
  // and if we miss the outer edge (due to small chisq surface range)
  // we get stuck on the inner one
  for(int ie = 0; ie < nexpt; ie++) {
    //lo limit
    //dm 90
    int bin = 0; double threshold = 2.71;
    for(int ix = 1; ix <= hd90l->GetNbinsX(); ix++) {
      if(hd[ie]->GetBinContent(ix) < threshold && hd[ie]->GetBinContent(ix) > 1E-8) {
        bin = TMath::Max(ix - 1, 1);
	break;
      }
    }//ix
    hd90l->SetBinContent(bin, hd90l->GetBinContent(bin) + 1);
    //dm 68
    bin = 0; threshold = 1.00;
    for(int ix = 1; ix <= hd68l->GetNbinsX(); ix++) {
      if(hd[ie]->GetBinContent(ix) < threshold && hd[ie]->GetBinContent(ix) > 1E-8) {
        bin = TMath::Max(ix - 1, 1);
	break;
      }
    }//ix
    hd68l->SetBinContent(bin, hd68l->GetBinContent(bin) + 1);
    //mix 90
    bin = 0; threshold = 2.71;
    for(int ix = 1; ix <= hm90l->GetNbinsX(); ix++) {
      if(hm[ie]->GetBinContent(ix) < threshold && hm[ie]->GetBinContent(ix) > 1E-8) {
        bin = TMath::Max(ix - 1, 1);
	break;
      }
    }//ix
    //if(hm[ie]->GetBinContent(1) < threshold)
    //bin = 1;
    hm90l->Fill(hm90l->GetBinCenter(bin));
    //mix 68
    bin = 0; threshold = 1.00;
    for(int ix = 1; ix <= hm68l->GetNbinsX(); ix++) {
      if(hm[ie]->GetBinContent(ix) < threshold && hm[ie]->GetBinContent(ix) > 1E-8) {
        bin = TMath::Max(ix - 1, 1);
	break;
      }
    }
    //if(hm[ie]->GetBinContent(1) < threshold)
    //bin = 1;
    hm68l->Fill(hm68l->GetBinCenter(bin));

    //hi limit
    //dm 90
    bin = 0; threshold = 2.71;
    for(int ix = hd90h->GetNbinsX(); ix >= 1; ix--) {
      if(hd[ie]->GetBinContent(ix) < threshold && hd[ie]->GetBinContent(ix) > 1E-8) {
        bin = TMath::Min(ix + 1, hd[ie]->GetNbinsX());
	break;
      }
    }//ix
    hd90h->SetBinContent(bin, hd90h->GetBinContent(bin) + 1);
    //dm 68
    bin = 0; threshold = 1.00;
    for(int ix = hd68h->GetNbinsX(); ix >= 1; ix--) {
      if(hd[ie]->GetBinContent(ix) < threshold && hd[ie]->GetBinContent(ix) > 1E-8) {
        bin = TMath::Min(ix + 1, hd[ie]->GetNbinsX());
	break;
      }
    }//ix
    hd68h->SetBinContent(bin, hd68h->GetBinContent(bin) + 1);
    //mix 90
    bin = 0; threshold = 2.71;
    for(int ix = hm90h->GetNbinsX() - 1; ix >= 1; ix--) {
      if(hm[ie]->GetBinContent(ix) < threshold && hm[ie]->GetBinContent(ix) > 1E-8) {
        bin = TMath::Min(ix + 1, hm[ie]->GetNbinsX());
	break;
      }
    }//ix
    hm90h->Fill(hm90h->GetBinCenter(bin));
    //mix 68
    bin = 0; threshold = 1.00;
    for(int ix = hm68l->GetNbinsX(); ix >= 1; ix--) {
      //cout << "\t" << ix << "\t" << hm[ie]->GetBinContent(ix) << endl;
      if(hm[ie]->GetBinContent(ix) < threshold && hm[ie]->GetBinContent(ix) > 1E-8) {
        bin = TMath::Min(ix + 1, hm[ie]->GetNbinsX());
	break;
      }
    }
    //cout << "expt: " << ie << "\t" << bin << "\t" << hm68l->GetBinCenter(bin) << endl;
    hm68h->Fill(hm68l->GetBinCenter(bin));
  }//ie

  //save histograms
  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.15);
  c->Divide(1,3); 
  c->SaveAs("test.pdf[");
  TFile fout("test.root", "RECREATE");
  for(int ie = minexpt; ie < maxexpt; ie++) {
    h [ie]->Write(Form("chisq_surface_%d",ie));
    hm[ie]->Write(Form("chisq_profile_mix23_%d",ie));
    hd[ie]->Write(Form("chisq_profile_dm32sq_%d",ie));
    c->cd(1);
    h[ie]->Draw("colz");
    c->cd(2);
    hm[ie]->Draw();
    c->cd(3);
    hd[ie]->Draw();
    c->SaveAs("test.pdf");
  }//ie
  c->SaveAs("test.pdf]");

  TLine * lmix68_l = new TLine(dmix68_l,0,dmix68_l,hm68l->GetMaximum() * 1.05);
  lmix68_l->SetLineColor(kRed);
  TLine * lmix68_h = new TLine(dmix68_h,0,dmix68_h,hm68l->GetMaximum() * 1.05);
  lmix68_h->SetLineColor(kRed);
  double pm68l = hm68l->Integral(0,hm68l->FindBin(dmix68_l)) / hm68l->Integral();
  cout << 1-pm68l << "\t" << "pm68l" << endl;
  double pm68h = hm68h->Integral(hm68h->FindBin(dmix68_h),hm68h->GetNbinsX()) / hm68h->Integral();
  cout << 1-pm68h << "\t" << "pm68h" << endl;
  double pm90l = hm90l->Integral(0,hm90l->FindBin(dmix90_l)) / hm90l->Integral();
  cout << 1-pm90l << "\t" << "pm90l" << endl;
  double pm90h = hm90h->Integral(hm90h->FindBin(dmix90_h),hm90h->GetNbinsX()) / hm90h->Integral();
  cout << 1-pm90h << "\t" << "pm90h" << endl;
  double pd68l = hd68l->Integral(0,hd68l->FindBin(ddm68_l)) / hd68l->Integral();
  cout << 1-pd68l << "\t" << "pd68l" << endl;
  double pd68h = hd68h->Integral(hd68h->FindBin(ddm68_h),hd68h->GetNbinsX()) / hd68h->Integral();
  cout << 1-pd68h << "\t" << "pd68h" << endl;
  double pd90l = hd90l->Integral(0,hd90l->FindBin(ddm90_l)) / hd90l->Integral();
  cout << 1-pd90l << "\t" << "pd90l" << endl;
  double pd90h = hd90h->Integral(hd90h->FindBin(ddm90_h),hd90h->GetNbinsX()) / hd90h->Integral();
  cout << 1-pd90h << "\t" << "pd90h" << endl;

  /*
 double dmix68_l = 0.456,
 double dmix68_h = 0.572,
 double dmix90_l = 0.428,
 double dmix90_h = 0.599,
 double ddm68_l = 0.00235,
 double ddm68_h = 0.00257,
 double ddm90_l = 0.00229,
 double ddm90_h = 0.00264
  */

  TCanvas * c2 = new TCanvas();
  c2->SaveAs("test2.pdf[");
  hm90l->Draw();
  c2->SaveAs("test2.pdf");
  hm68l->Draw();
  lmix68_l->Draw();
  c2->SaveAs("test2.pdf");
  hd90l->Draw();
  c2->SaveAs("test2.pdf");
  hd68l->Draw();
  c2->SaveAs("test2.pdf");
  hm90h->Draw();
  c2->SaveAs("test2.pdf");
  hm68h->Draw();
  c2->SaveAs("test2.pdf");
  hd90h->Draw();
  c2->SaveAs("test2.pdf");
  hd68h->Draw();
  c2->SaveAs("test2.pdf");
  c2->SaveAs("test2.pdf]");

  hm68l->Draw();
  lmix68_l->Draw();
  hm68h->SetLineColor(kBlue);
  hm68h->Draw("same");
  lmix68_h->Draw();
  c2->SaveAs("sensitivity_mix_limits.pdf");

}
