void plot_sensitivity_sinsqtheta(
		      const char * filename  = "2013results/SENSITIVITY/sinsqtheta/run1-4abf/combination/final.root",
		      double       in_mix    = 0.514,
		      double       in_dm     = 0.002463
)
{
  TFile fsens(filename);
  TH2D * h_combchisq = 0;
  fsens.GetObject("chisq_surface",h_combchisq);

  //make the combined chisq surface a an average chisq surface (average)
  TH2D * h_avchisq = new TH2D(*h_combchisq);
  for(int ix = 1; ix <= h_combchisq->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= h_combchisq->GetNbinsX(); iy++) {
      double chisq_combine = h_combchisq->GetBinContent(ix,iy);
      double chisq_average = (chisq_combine < 700) ? chisq_combine / (double)5 : 75;
      /*
      //hack for a missing job
      if(chisq_combine < 1)
	chisq_average = 80;
      else if(chisq_combine < 100)
	chisq_average = chisq_combine / (double)1;
      else if(chisq_combine < 200)
	chisq_average = chisq_combine / (double)2;
      */
      h_avchisq->SetBinContent(ix,iy, chisq_average);
    }//iy
  }//ix

  //find the minimum
  double chisq_min = h_avchisq->GetMinimum();
  int minbin = h_avchisq->GetMinimumBin();
  int minx, miny, minz;
  h_avchisq->GetBinXYZ(minbin, minx, miny, minz);
  cout << endl << endl << "min chisq " << chisq_min << " at " << h_avchisq->GetXaxis()->GetBinCenter(minx) << " " << h_avchisq->GetYaxis()->GetBinCenter(miny) << endl;

  //find the value of the input
  double chisq_input = h_avchisq->GetBinContent(h_avchisq->FindBin(in_mix,in_dm));
  cout << chisq_input << endl;
  chisq_min = chisq_input;

  //make the average chisq surface a dchisq surface (subtract min)
  TH2D * h_sens_90 = new TH2D(*h_combchisq);
  for(int ix = 1; ix <= h_combchisq->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= h_combchisq->GetNbinsX(); iy++) {
      double chisq_av = h_avchisq->GetBinContent(ix,iy);
      double dchisq   = chisq_av - chisq_min;
      h_sens_90->SetBinContent(ix, iy, dchisq);
    }//iy
  }//ix

  //set up the sensitivity histogram for contour drawing
  h_sens_90->SetContour(1);
  h_sens_90->SetContourLevel(0,4.61);
  h_sens_90->SetLineStyle(1);
  TH2D * h_sens_68 = new TH2D(*h_sens_90);
  h_sens_68->SetContourLevel(0,2.30);
  h_sens_68->SetLineStyle(2);

  //Input marker
  TMarker * m_input = new TMarker(in_mix,in_dm,29);
  m_input->SetMarkerSize(1.8);

  //
  //Draw
  //
  // param space
  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.16);
  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.0;
  double mixing23_max   = 1.0;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1E-3;
  double dm23sq_max     = 6E-3;
  double dm23sq_step    = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + TMath::CeilNint ((mixing23_max - mixing23_min) / mixing23_step);
  int dm23sq_nbins   = 1 + TMath::CeilNint ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  );
  TH2D * osc_param_space = new TH2D("osc_param_space", "",
       mixing23_nbins, 
       mixing23_min - 0.5*mixing23_step, 
       mixing23_max + 0.5*mixing23_step,
       dm23sq_nbins, 
       dm23sq_min - 0.5*dm23sq_step, 
       dm23sq_max + 0.5*dm23sq_step);
  osc_param_space->GetXaxis()->SetTitle("sin^{2}(#theta_{23})"); // adjust if mixing param is sinsqtheta rather than sinsq2theta...
  osc_param_space->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  osc_param_space->GetYaxis()->SetTitleOffset(1.4);

  //run 1-4 result
  TFile fdata1234a("2013results/DATA/run1-4_nuis3fixdcp_sintheta_normal/combination/run1-4_nuis3fixdcp_sintheta_normal-run_0-expt_6.root");
  TH2D * h_data1234a_90 = 0, * h_data1234a_68 = 0;
  fdata1234a.GetObject("dchisq_cont68",h_data1234a_68);
  fdata1234a.GetObject("dchisq_cont90",h_data1234a_90);

  //style
  h_sens_90->SetLineWidth(2);
  h_sens_90->SetLineColor(kRed);
  h_sens_68->SetLineWidth(2);
  h_sens_68->SetLineColor(kRed);
  h_data1234a_90->SetLineWidth(2);
  h_data1234a_68->SetLineStyle(2);
  h_data1234a_68->SetLineWidth(2);

  //draw
  osc_param_space->GetXaxis()->SetRangeUser(0.35,0.65);
  osc_param_space->GetYaxis()->SetRangeUser(0.0022,0.0033);
  osc_param_space->Draw();
  h_sens_90->Draw("cont3same");
  h_sens_68->Draw("cont3same");
  h_data1234a_68->Draw("cont3same");
  h_data1234a_90->Draw("cont3same");
  m_input->Draw();
  //legend
  TLegend * leg = new TLegend(0.18, 0.65, 0.45, 0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(132);
  leg->SetTextSize(0.05);
  leg->AddEntry(h_sens_68,  "Run 1-4a 68% C.L. sensitivity","l");
  leg->AddEntry(h_sens_90,  "Run 1-4a 90% C.L. sensitivity","l");
  leg->AddEntry(h_data1234a_68,  "Run 1-4a 68% C.L. result","l");
  leg->AddEntry(h_data1234a_90,  "Run 1-4a 90% C.L. result","l");
  leg->AddEntry(m_input,   "Run 1-4a best fit","p");
  leg->Draw();
  //hiding device
  TLegend * leg2 = new TLegend(0.86, 0.53, 0.91, 0.61);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->Draw();
  c->SaveAs("sensitivity.pdf");

}
