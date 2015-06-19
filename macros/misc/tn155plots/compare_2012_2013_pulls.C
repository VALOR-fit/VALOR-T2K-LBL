void compare_2012_2013_pulls(
			     const char * f2012 = "2012results/DATA/2all/fitresults_3.root",
			     const char * f2013 = "2013results/DATA/run1-3_nuis3fixdcp_sintheta_normal/combination/run1-3_nuis3fixdcp_sintheta_normal_2_pulls.root",
			     const char * f2013_2 = "2013results/DATA/run1-3_2013_nuis3fixdcp_sintheta_normal/combination/run1-3_2013_nuis3fixdcp_sintheta_normal_pulls.root",
			     double       ymax  = 0.7)
{
  TFile f2(f2012);
  TH1D * h2 = 0;
  f2.GetObject("hsystpull",h2);

  TFile f3(f2013);
  TH1D * h3 = 0;
  f3.GetObject("pulls",h3);

  TFile f4(f2013_2);
  TH1D * h4 = 0;
  f4.GetObject("pulls",h4);

  TH1D * hpull2012   = new TH1D("hpull2012",  "hpull2012",  50,0,50);
  TH1D * hpull2013   = new TH1D("hpull2013",  "hpull2013",  50,0,50);
  TH1D * hpull2013_2 = new TH1D("hpull2013_2","hpull2013_2",50,0,50);
  int ip = 0;

  for(int i = 1; i <= h3->GetNbinsX(); i++) {
    int found = -1;
    TString s3long(h3->GetXaxis()->GetBinLabel(i));
    TString s3 = s3long.Strip(TString::kBoth,' ');//get rid of spaces
    for(int j = 1; j <= h2->GetNbinsX(); j++) {
      TString s2long(h2->GetXaxis()->GetBinLabel(j));
      TString s2 = s2long.Strip(TString::kBoth,' ');//get rid of spaces
      if(s2.Contains(s3)) {
	found = j;
	break;
      }
    }//j
    if(found <= 0)
      cout << "not found matching entries for " << s3 << endl;
    else {
      ip++;
      if(ip % 2) {
	hpull2012->  GetXaxis()->SetBinLabel(ip,s3.Data());
	hpull2013->  GetXaxis()->SetBinLabel(ip,s3.Data());
	hpull2013_2->GetXaxis()->SetBinLabel(ip,s3.Data());
      }
      else {
	hpull2012->  GetXaxis()->SetBinLabel(ip,Form("%s                    ",s3.Data()));
	hpull2013->  GetXaxis()->SetBinLabel(ip,Form("%s",s3.Data()));
	hpull2013_2->GetXaxis()->SetBinLabel(ip,Form("%s",s3.Data()));
      }
      hpull2012->  SetBinContent(ip, h2->GetBinContent(found));
      hpull2013->  SetBinContent(ip, h3->GetBinContent(i));
      hpull2013_2->SetBinContent(ip, h4->GetBinContent(i));
      cout << s3 << "\t" << h2->GetBinContent(found) << "\t" << h3->GetBinContent(i) << "\t" << h4->GetBinContent(i) << endl;
    }
  }//i

  TCanvas * c = new TCanvas();
  c->SetGridx(1);
  c->SetBottomMargin(0.4);
  c->cd();

  hpull2012->SetLineColor(kRed);
  hpull2012->SetMarkerColor(kRed);
  hpull2012->GetYaxis()->SetRangeUser(-ymax,ymax);
  hpull2012->GetXaxis()->LabelsOption("v");
  hpull2012->GetXaxis()->SetRange(1,ip);
  hpull2012-> GetYaxis()->SetTitle("Pull = (bestfit - input) / #sigma_{fit}");
  hpull2013->SetMarkerStyle(26);
  hpull2013_2->SetMarkerStyle(25);
  hpull2013_2->SetMarkerColor(kBlue);

  TLegend * l = new TLegend(0.0,0.0,0.2,0.2);
  l->SetFillColor(kWhite);
  l->SetBorderSize(0);
  l->SetTextSize(0.045);
  l->AddEntry(hpull2012,  "TN141 fit (2012 data)","P");
  l->AddEntry(hpull2013,  "TN155 fit (2012 data)","P");
  l->AddEntry(hpull2013_2,"TN155 fit (2013 data)","P");

  hpull2012->Draw("P");
  hpull2013->Draw("Psame");
  hpull2013_2->Draw("Psame");
  TLine * lpull = new TLine(0,0,hpull2012->GetXaxis()->GetLast(),0);
  lpull->Draw();
  l->Draw();
  c->SaveAs("pulls_comparison.pdf");
  c->SaveAs("pulls_comparison.eps");
}

