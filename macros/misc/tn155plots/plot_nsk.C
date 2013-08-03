{
  using namespace t2kana;

  TFile f_tun("2013results/SPECTRA/nsk_plots/banfftuned_normal/nsk_vs_osc_params.root");
  TFile f_flux("2013results/SPECTRA/nsk_plots/fluxtuned_normal/nsk_vs_osc_params.root");
  TFile f_nom("2013results/SPECTRA/nsk_plots/nominal_normal/nsk_vs_osc_params.root");
  TFile f_inv("2013results/SPECTRA/nsk_plots/banfftuned_inverted/nsk_vs_osc_params.root");
  TFile f_t23("2013results/SPECTRA/nsk_plots/thrownon23_normal/nsk_vs_osc_params.root");

  //histograms
  TH2D * h_tun = 0;
  f_tun.GetObject("nsk_vs_osc_params",h_tun);
  TH2D * h_flux = 0;
  f_flux.GetObject("nsk_vs_osc_params",h_flux);
  TH2D * h_nom = 0;
  f_nom.GetObject("nsk_vs_osc_params",h_nom);
  TH2D * h_inv = 0;
  f_inv.GetObject("nsk_vs_osc_params",h_inv);
  TH2D * h_t23 = 0;
  f_t23.GetObject("nsk_vs_osc_params",h_t23);

  //ratios
  TH2D * h_ratio_tunflux = new TH2D(*h_tun);
  TH2D * h_ratio_nomflux = new TH2D(*h_tun);
  TH2D * h_ratio_invtun = new TH2D(*h_tun);
  TH2D * h_ratio_t23tun = new TH2D(*h_tun);
  for(int i = 1; i <= h_tun->GetNbinsX(); i++) {
    for(int j = 1; j <= h_tun->GetNbinsY(); j++) {
      double n_tun = h_tun->GetBinContent(i,j);
      double n_flux = h_flux->GetBinContent(i,j);
      double n_nom = h_nom->GetBinContent(i,j);
      double n_inv = h_inv->GetBinContent(i,j);
      double n_t23 = h_t23->GetBinContent(i,j);
      h_ratio_tunflux->SetBinContent(i,j, (n_tun-n_flux)/n_flux);
      h_ratio_nomflux->SetBinContent(i,j, (n_flux-n_nom)/n_nom);
      h_ratio_invtun->SetBinContent(i,j, (n_inv-n_tun)/n_tun);
      h_ratio_t23tun->SetBinContent(i,j, (n_t23-n_tun)/n_tun);
    }
  }


  TCanvas * c = new TCanvas();
  c->SetLeftMargin(0.15);
  c->SetRightMargin(0.15);
  h_ratio_tunflux->GetYaxis()->SetTitleOffset(1.3);
  h_ratio_nomflux->GetYaxis()->SetTitleOffset(1.3);
  h_tun->GetYaxis()->SetTitleOffset(1.3);
  h_ratio_invtun->GetYaxis()->SetTitleOffset(1.3);
  h_ratio_t23tun->GetYaxis()->SetTitleOffset(1.3);

  h_ratio_tunflux->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");;
  h_ratio_nomflux->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");;
  h_tun->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");;
  h_ratio_invtun->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");;
  h_ratio_t23tun->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");;

  /*
  //minos
  TFile fext("$T2KOSCANA/inputs/contours/t2k_2013paper_likelihoodratio_contour_and_deltachisq.root");
  TGraph * gext = 0;
  fext.GetObject("g_contour_90_o2",gext);
  gext->SetLineColor(kBlack);
  gext->SetLineWidth(2);
  */
  //t2k run1-3
  TFile fext("$T2KOSCANA/2012results/DATA/sintheta23/combine/combination/fitresults_3.root");
  TH2D * hext = 0;
  fext.GetObject("hdchisq",hext);
  hext->SetLineColor(kBlack);
  hext->SetLineWidth(2);
  hext->SetContour(1);
  hext->SetContourLevel(0,4.61);

  //draw ratios
  h_ratio_tunflux->Draw("colz");
  hext->Draw("cont3same");
  //gext->Draw("l");
  c->SaveAs("nsk_diff_tunflux.eps");
  h_ratio_nomflux->Draw("colz");
  hext->Draw("cont3same");
  //gext->Draw("l");
  c->SaveAs("nsk_diff_nomflux.eps");
  h_ratio_invtun->Draw("colz");
  hext->Draw("cont3same");
  //gext->Draw("l");
  c->SaveAs("nsk_diff_inv.eps");
  h_ratio_t23tun->Draw("colz");
  hext->Draw("cont3same");
  //gext->Draw("l");
  c->SaveAs("nsk_diff_non23throw.eps");

  //draw normal nsk
  h_tun->Draw("colz");
  hext->Draw("cont3same");
  //gext->Draw("l");
  c->SaveAs("nsk_nominalbanff.pdf");
  c->SaveAs("nsk_nominalbanff.eps");

  double nmax = 0;
  double btunmin = 9999, btunmax = -9999; //banff tuned min/max ratios
  double ftunmin = 9999, ftunmax = -9999; //flux tuned min/max ratios
  for(int ix = 1; ix <= hext->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= hext->GetNbinsX(); iy++) {
      double dchisq = hext->GetBinContent(ix,iy);
      if(dchisq <= 4.61 && dchisq > 0) {
	double mix = hext->GetXaxis()->GetBinCenter(ix);
	double dm2 = hext->GetYaxis()->GetBinCenter(iy);
	double ncurr = h_tun->GetBinContent(h_tun->FindBin(mix,dm2));
	double btuncurr = h_ratio_tunflux->GetBinContent(h_ratio_tunflux->FindBin(mix,dm2));
	double ftuncurr = h_ratio_nomflux->GetBinContent(h_ratio_nomflux->FindBin(mix,dm2));
	if(ncurr > nmax) {
	  nmax = ncurr;
	  cout << "Banff tuned events        \t" << mix << "\t" << dm2 << "\t" << nmax << endl;
	}
	if(btuncurr < btunmin) {
	  btunmin = btuncurr;
	  cout << "Minimum flux->Banff ratio \t" << mix << "\t" << dm2 << "\t" << btunmin << endl;
	}
	if(btuncurr > btunmax) {
	  btunmax = btuncurr;
	  cout << "Maximum flux->Banff ratio \t" << mix << "\t" << dm2 << "\t" << btunmax << endl;
	}
	if(ftuncurr < ftunmin) {
	  ftunmin = ftuncurr;
	  cout << "Minimum nom->flux ratio   \t" << mix << "\t" << dm2 << "\t" << ftunmin << endl;
	}
	if(ftuncurr > ftunmax) {
	  ftunmax = ftuncurr;
	  cout << "Maximum nom->flux ratio   \t" << mix << "\t" << dm2 << "\t" << ftunmax << endl;
	}
      }
    }//iy
  }//ix
  cout << "Banff tuned events        \t" << nmax << endl;
  cout << "Minimum flux->Banff ratio \t" << btunmin << endl;
  cout << "Maximum flux->Banff ratio \t" << btunmax << endl;
  cout << "Minimum nom->flux ratio   \t" << ftunmin << endl;
  cout << "Maximum nom->flux ratio   \t" << ftunmax << endl;

}


/*
do the following in a ROOT session to look at Nsk at the edge of the T2K contour (to check the above calculation):



  TFile f_tun("2013results/SPECTRA/nsk_plots/banfftuned_normal/nsk_vs_osc_params.root");
  gStyle->SetPaintTextFormat("3.0f");
  nsk_vs_osc_params->Draw("text");
  TFile fext("$T2KOSCANA/2012results/DATA/sintheta23/combine/combination/fitresults_3.root");
  TH2D * hext = 0;
  fext.GetObject("hdchisq",hext);
  hext->SetLineColor(kBlack);
  hext->SetLineWidth(2);
  hext->SetContour(1);
  hext->SetContourLevel(0,4.61);
  hext->Draw("cont3same");


*/
