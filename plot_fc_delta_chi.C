//
// Plot Feldman-Cousins \delta\chi^{2}
//
// C.Andreopoulos
//

{
  const double sinsq2theta_min = 0.5;
  const double sinsq2theta_max = 1.0;
  const double dchisq_2p_68    = 2.30; /// canonical 68% CL \delta\chi^{2} for 2 parameters
  const double dchisq_2p_90    = 4.61; /// canonical 90% CL \delta\chi^{2} for 2 parameters

  // Input files & histograms

  TFile * fstatonly = new TFile("DeltaChi2MCStat.root",        "read");
  TFile * fstatsyst = new TFile("DeltaChi2MCStatPlusSyst.root","read");

  TH2D * hdchisq68_statonly = (TH2D*) fstatonly -> Get("hDeltaChi2Surface68"); 
  TH2D * hdchisq90_statonly = (TH2D*) fstatonly -> Get("hDeltaChi2Surface90"); 
  TH2D * hdchisq68_statsyst = (TH2D*) fstatsyst -> Get("hDeltaChi2Surface68"); 
  TH2D * hdchisq90_statsyst = (TH2D*) fstatsyst -> Get("hDeltaChi2Surface90"); 

  hdchisq68_statonly -> Smooth();
  hdchisq90_statonly -> Smooth();
  hdchisq68_statsyst -> Smooth();
  hdchisq90_statsyst -> Smooth();

  gStyle->SetOptStat(0);

  //
  // Plot \delta\chi^{2}_{68} and \delta\chi^{2}_{90} both for 'stat' and 'stat+syst' cases.
  //

  hdchisq68_statonly  -> GetXaxis() -> SetTitle("sin^{2}2#theta_{23}");
  hdchisq90_statonly  -> GetXaxis() -> SetTitle("sin^{2}2#theta_{23}");
  hdchisq68_statsyst  -> GetXaxis() -> SetTitle("sin^{2}2#theta_{23}");
  hdchisq90_statsyst  -> GetXaxis() -> SetTitle("sin^{2}2#theta_{23}");
  hdchisq68_statonly  -> GetYaxis() -> SetTitle("|#Deltam^{2}_{23}| (eV^{2}/c^{4}");
  hdchisq90_statonly  -> GetYaxis() -> SetTitle("|#Deltam^{2}_{23}| (eV^{2}/c^{4}");
  hdchisq68_statsyst  -> GetYaxis() -> SetTitle("|#Deltam^{2}_{23}| (eV^{2}/c^{4}");
  hdchisq90_statsyst  -> GetYaxis() -> SetTitle("|#Deltam^{2}_{23}| (eV^{2}/c^{4}");
  hdchisq68_statonly  -> GetYaxis() -> SetTitleOffset(1.6);
  hdchisq90_statonly  -> GetYaxis() -> SetTitleOffset(1.6);
  hdchisq68_statsyst  -> GetYaxis() -> SetTitleOffset(1.6);
  hdchisq90_statsyst  -> GetYaxis() -> SetTitleOffset(1.6);

  TCanvas * c1 = new TCanvas("c1","",10,20,500,500);
  c1->Divide(2,2);
  c1->cd(1);
  c1->GetPad(1)->SetLeftMargin(0.3);
  c1->GetPad(1)->SetRightMargin(0.15);
  hdchisq68_statonly->Draw("colz");
  c1->cd(2);
  c1->GetPad(2)->SetLeftMargin(0.15);
  c1->GetPad(2)->SetRightMargin(0.3);
  hdchisq90_statonly->Draw("colz");
  c1->cd(3);
  c1->GetPad(3)->SetLeftMargin(0.3);
  c1->GetPad(3)->SetRightMargin(0.15);
  hdchisq68_statsyst->Draw("colz");
  c1->cd(4);
  c1->GetPad(4)->SetLeftMargin(0.15);
  c1->GetPad(4)->SetRightMargin(0.3);
  hdchisq90_statsyst->Draw("colz");
  c1->Update();


  //
  // Plot \delta\chi^{2}_{68} and \delta\chi^{2}_{90} vs sinsq2theta for 3 dmsq.
  //

  double dmsq1 = 2.0E-3;
  double dmsq2 = 2.5E-3;
  double dmsq3 = 3.0E-3;

  int n = hdchisq90_statonly->GetXaxis()->GetNbins();

  double * sinsq      = new double[n];
  double * dchisq1_90 = new double[n];
  double * dchisq2_90 = new double[n];
  double * dchisq3_90 = new double[n];
  double * dchisq1_68 = new double[n];
  double * dchisq2_68 = new double[n];
  double * dchisq3_68 = new double[n];

  int idmsq1 = hdchisq90_statonly->GetYaxis()->FindBin(dmsq1);
  for(int i=0; i<n; i++) {
    dchisq1_90[i] = hdchisq90_statonly->GetBinContent(i+1,idmsq1);
    dchisq1_68[i] = hdchisq68_statonly->GetBinContent(i+1,idmsq1);
    sinsq[i]      = hdchisq90_statonly->GetXaxis()->GetBinCenter(i+1);
  }
  int idmsq2 = hdchisq90_statonly->GetYaxis()->FindBin(dmsq2);
  for(int i=0; i<n; i++) {
    dchisq2_90[i] = hdchisq90_statonly->GetBinContent(i+1,idmsq2);
    dchisq2_68[i] = hdchisq68_statonly->GetBinContent(i+1,idmsq2);
  }
  int idmsq3 = hdchisq90_statonly->GetYaxis()->FindBin(dmsq3);
  for(int i=0; i<n; i++) {
    dchisq3_90[i] = hdchisq90_statonly->GetBinContent(i+1,idmsq3);
    dchisq3_68[i] = hdchisq68_statonly->GetBinContent(i+1,idmsq3);
  }


  TCanvas * c2 = new TCanvas("c2","",10,20,500,500);
  c2->Divide(2,1);
  c2->cd(1);
  TGraph * gr1_68 = new TGraph(n,sinsq,dchisq1_68);
  TGraph * gr2_68 = new TGraph(n,sinsq,dchisq2_68);
  TGraph * gr3_68 = new TGraph(n,sinsq,dchisq3_68);
  gr1_68->SetLineColor(kBlue);
  gr2_68->SetLineColor(kGreen);
  gr3_68->SetLineColor(kRed);
  TH1F * hframe_68 = (TH1F*)c2->GetPad(1)->DrawFrame(sinsq2theta_min, 0., sinsq2theta_max, 4.);
  hframe_68->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  hframe_68->GetYaxis()->SetTitle("#Delta#chi^{2}_{68}");
  gr1_68->Draw("c");
  gr2_68->Draw("c");
  gr3_68->Draw("c");
  TLegend * legend_68 = new TLegend(0.5, 0.5, 0.9, 0.9);
  legend_68->SetFillColor(0);
  legend_68->SetBorderSize(0);
  legend_68->AddEntry(gr1_68, Form("#Delta#chi^{2}_{68}, #Deltam^{2}_{23} = %.5f eV^{2}/c^{4}", dmsq1), "L");
  legend_68->AddEntry(gr2_68, Form("#Delta#chi^{2}_{68}, #Deltam^{2}_{23} = %.5f eV^{2}/c^{4}", dmsq2), "L");
  legend_68->AddEntry(gr3_68, Form("#Delta#chi^{2}_{68}, #Deltam^{2}_{23} = %.5f eV^{2}/c^{4}", dmsq3), "L");
  legend_68->Draw();
  TLine * lc_68 = new TLine(sinsq2theta_min, dchisq_2p_68, sinsq2theta_max, dchisq_2p_68); /// canonical \delta\chi^{2}
  lc_68->Draw();
  c2->cd(2);
  TGraph * gr1_90 = new TGraph(n,sinsq,dchisq1_90);
  TGraph * gr2_90 = new TGraph(n,sinsq,dchisq2_90);
  TGraph * gr3_90 = new TGraph(n,sinsq,dchisq3_90);
  gr1_90->SetLineColor(kBlue);
  gr2_90->SetLineColor(kGreen);
  gr3_90->SetLineColor(kRed);
  TH1F * hframe_90 = (TH1F*)c2->GetPad(1)->DrawFrame(sinsq2theta_min, 0., sinsq2theta_max, 8.);
  hframe_90->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  hframe_90->GetYaxis()->SetTitle("#Delta#chi^{2}_{68}");
  gr1_90->Draw("c");
  gr2_90->Draw("c");
  gr3_90->Draw("c");
  TLegend * legend_90 = new TLegend(0.5, 0.5, 0.9, 0.9);
  legend_90->SetFillColor(0);
  legend_90->SetBorderSize(0);
  legend_90->AddEntry(gr1_90, Form("#Delta#chi^{2}_{90}, #Deltam^{2}_{23} = %.5f eV^{2}/c^{2}", dmsq1), "L");
  legend_90->AddEntry(gr2_90, Form("#Delta#chi^{2}_{90}, #Deltam^{2}_{23} = %.5f eV^{2}/c^{2}", dmsq2), "L");
  legend_90->AddEntry(gr3_90, Form("#Delta#chi^{2}_{90}, #Deltam^{2}_{23} = %.5f eV^{2}/c^{2}", dmsq3), "L");
  legend_90->Draw();
  TLine * lc_90 = new TLine(sinsq2theta_min, dchisq_2p_90, sinsq2theta_max, dchisq_2p_90); /// canonical \delta\chi^{2}
  lc_90->Draw();
}

