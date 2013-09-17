void plot_1d_dchisq
(
 const char * fname    = "$T2KOSCANA/2012results/1DFC/sintheta/dm/sinthetadm-run_0-expts_3_3-bins_-1_-1.root",
 const char * hname    = "chisqwithsyst_dm23sq_3",
 double       chisq_bf = 56.03,
 const char * xtitle   = "|#Deltam^{2}_{32}|",
 double       ymax     = 8
 )
{
  //get chisq histogram
  TFile f(fname);
  TH1D * htemp = 0;
  f.GetObject(hname,htemp);
  TH1D * h = new TH1D(*htemp);

  //subtract chisq_bf
  for(int i = 1; i <= h->GetNbinsX(); i++) {
    double chisq  = h->GetBinContent(i);
    double center = h->GetBinCenter(i);
    if(TMath::Abs(center - 0.480) < 1E-8
       || TMath::Abs(center - 0.508) < 1E-8
       || TMath::Abs(center - 0.580) < 1E-8
       || TMath::Abs(center - 0.653) < 1E-8
       ) {
      double chisq1 = htemp->GetBinContent(i+1);
      double chisq2 = htemp->GetBinContent(i-1);
      chisq = (chisq1 + chisq2) / 2;
    }
    h->SetBinContent(i, chisq-chisq_bf);
    if(chisq < chisq_bf)
      cout << "dchisq is " << chisq - chisq_bf << " in bin with center " << h->GetBinCenter(i) << endl;
  }

  //canvas
  TCanvas * c = new TCanvas();
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.1);

  //1d dchisq
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle("#Delta#chi^{2}");
  h->GetYaxis()->SetRangeUser(0,ymax);
  h->SetLineWidth(2);
  h->Draw();

  //const dchisq lines
  TLine * l68 = new TLine(h->GetXaxis()->GetXmin(),1.00,h->GetXaxis()->GetXmax(),1.00);
  TLine * l90 = new TLine(h->GetXaxis()->GetXmin(),2.71,h->GetXaxis()->GetXmax(),2.71);
  l68->SetLineStyle(3);
  l90->SetLineStyle(2);
  l68->Draw();
  l90->Draw();

  //legend
  TLegend * leg = new TLegend(0.35,0.7,0.65,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.06);
  leg->AddEntry(h,"sin^{2}#theta_{23} fit","l");
  leg->AddEntry(l68,"68% C.L.","l");
  leg->AddEntry(l90,"90% C.L.","l");
  //leg->Draw();
  c->SaveAs("1D_contour.pdf");
  c->SaveAs("1D_contour.eps");

  //limits
  for(int i = 1; i <= h->GetNbinsX(); i++) {
    if(h->GetBinContent(i) < 2.71) {
      cout << "lower 90% CL " << h->GetBinCenter(i-1) << "\t" << h->GetBinContent(i-1) << endl;
      break;
    }
  }
  for(int i = h->GetNbinsX(); i >= 1; i--) {
    if(h->GetBinContent(i) < 2.71) {
      cout << "upper 90% CL " << h->GetBinCenter(i+1) << "\t" << h->GetBinContent(i+1) << endl;
      break;
    }
  }
  for(int i = 1; i <= h->GetNbinsX(); i++) {
    if(h->GetBinContent(i) < 1.00) {
      cout << "lower 68% CL " << h->GetBinCenter(i-1) << "\t" << h->GetBinContent(i-1) << endl;
      break;
    }
  }
  for(int i = h->GetNbinsX(); i >= 1; i--) {
    if(h->GetBinContent(i) < 1.00) {
      cout << "upper 68% CL " << h->GetBinCenter(i+1) << "\t" << h->GetBinContent(i+1) << endl;
      break;
    }
  }

}

