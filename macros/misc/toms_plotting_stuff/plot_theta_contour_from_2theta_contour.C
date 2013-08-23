//plot contours in sinsqtheta23-dm32sq space
//
//takes 3 inputs
// 2 sinsq2theta23 contours (one from each octant
// 1 sinsqthet23 contour
// 1 contour from numu-only part of joint analysis
//sinsq2theta23 contours are converted to sinsqtheta23 contours

void plot_theta_contour_from_2theta_contour
(
 const char * filename1 = "2012results/NEWRESULTS/datafit/dchisq_2all.root",
 const char * filename2 = "2012results/DATA/octant2_23only/combination/fitresults_3.root",
 const char * filename_thetafit = "2012results/DATA/sintheta23/combine/combination/fitresults_3.root",
 const char * filename_thetafit_DM = "2012results/DATA/sintheta23_DM/contours_disap.root",
 bool draw_DM = false
)
{
  TFile f1(filename1);
  TH2D * hsinsq2theta1 = 0;
  f1.GetObject("hdchisq",hsinsq2theta1);
  TH2D * hsinsq2theta1_68 = new TH2D(*hsinsq2theta1);
  TH2D * hsinsq2theta1_90 = new TH2D(*hsinsq2theta1);
  hsinsq2theta1_68->SetContour(1);
  hsinsq2theta1_68->SetContourLevel(0,2.30);
  hsinsq2theta1_68->SetLineStyle(3);
  hsinsq2theta1_90->SetContour(1);
  hsinsq2theta1_90->SetContourLevel(0,4.61);
  hsinsq2theta1_90->SetLineStyle(1);

  TFile f2(filename2);
  TH2D * hsinsq2theta2 = 0;
  f2.GetObject("hdchisq",hsinsq2theta2);
  TH2D * hsinsq2theta2_68 = new TH2D(*hsinsq2theta2);
  TH2D * hsinsq2theta2_90 = new TH2D(*hsinsq2theta2);
  hsinsq2theta2_68->SetContour(1);
  hsinsq2theta2_68->SetContourLevel(0,2.30);
  hsinsq2theta2_68->SetLineStyle(3);
  hsinsq2theta2_90->SetContour(1);
  hsinsq2theta2_90->SetContourLevel(0,4.61);
  hsinsq2theta2_90->SetLineStyle(1);

  TFile f3(filename_thetafit);
  TH2D * hsinsqtheta3 = 0;
  f3.GetObject("hdchisq",hsinsqtheta3);
  TH2D * hsinsqtheta3_68 = new TH2D(*hsinsqtheta3);
  TH2D * hsinsqtheta3_90 = new TH2D(*hsinsqtheta3);
  hsinsqtheta3_68->SetContour(1);
  hsinsqtheta3_68->SetContourLevel(0,2.30);
  hsinsqtheta3_68->SetLineStyle(3);
  hsinsqtheta3_90->SetContour(1);
  hsinsqtheta3_90->SetContourLevel(0,4.61);
  hsinsqtheta3_90->SetLineStyle(1);
  hsinsqtheta3_68->SetLineColor(kBlue);
  hsinsqtheta3_90->SetLineColor(kBlue);

  TFile f4(filename_thetafit_DM);
  TH2D * hsinsqtheta4_DM_68 = 0, * hsinsqtheta4_DM_90 = 0;
  f4.GetObject("h_theta23_dm23sq_disap",   hsinsqtheta4_DM_90);
  f4.GetObject("h_theta23_dm23sq_disap_68",hsinsqtheta4_DM_68);
  hsinsqtheta4_DM_68->SetLineStyle(3);
  hsinsqtheta4_DM_90->SetLineStyle(1);
  hsinsqtheta4_DM_68->SetLineColor(kYellow+1);
  hsinsqtheta4_DM_90->SetLineColor(kYellow+1);

  //convert from DMsq to dm32sq
  // get y axis array
  int nxbins_DM = hsinsqtheta4_DM_68->GetNbinsX();
  int nybins_DM = hsinsqtheta4_DM_68->GetNbinsY();
  double xbins_dm32sq[100] = {0};
  double ybins_dm32sq[100] = {0};
  for(int i = 0; i <= nybins_DM; i++) {
    ybins_dm32sq[i] = hsinsqtheta4_DM_68->GetYaxis()->GetBinUpEdge(i) - 3.75E-5; //subtract dm12sq/2
  }
  for(int i = 0; i <= nxbins_DM; i++) {
    xbins_dm32sq[i] = hsinsqtheta4_DM_68->GetXaxis()->GetBinUpEdge(i);
  }
  //make new histograms
  TH2D * hsinsqtheta4_68 = new TH2D("","",nxbins_DM,xbins_dm32sq,nybins_DM,ybins_dm32sq);
  TH2D * hsinsqtheta4_90 = new TH2D("","",nxbins_DM,xbins_dm32sq,nybins_DM,ybins_dm32sq);
  for(int ix = 1; ix <= nxbins_DM; ix++) {
    for(int iy = 1; iy <= nybins_DM; iy++) {
      hsinsqtheta4_68->SetBinContent(ix,iy, hsinsqtheta4_DM_68->GetBinContent(ix,iy));
      hsinsqtheta4_90->SetBinContent(ix,iy, hsinsqtheta4_DM_90->GetBinContent(ix,iy));
    }
  }
  hsinsqtheta4_68->SetContour(1);
  hsinsqtheta4_68->SetContourLevel(0,0);
  hsinsqtheta4_90->SetContour(1);
  hsinsqtheta4_90->SetContourLevel(0,0);
  hsinsqtheta4_68->SetLineStyle(3);
  hsinsqtheta4_90->SetLineStyle(1);
  hsinsqtheta4_68->SetLineColor(kYellow+1);
  hsinsqtheta4_90->SetLineColor(kYellow+1);

  

  //make an array of sinsqtheta values from sinsq2theta axis
  const int nbinsx1  = hsinsq2theta1->GetNbinsX();
  const int nbinsy1  = hsinsq2theta1->GetNbinsY();
  cout << nbinsx1 << "\t" << nbinsy1 << " bins" << endl;
  const double binminy1 = hsinsq2theta1->GetYaxis()->GetXmin();
  const double binmaxy1 = hsinsq2theta1->GetYaxis()->GetXmax();
  double axis1o1[100] = {0};
  double axis1o2[100] = {0};
  for(int ix = 1; ix <= nbinsx1; ix++) {
    //double sinsq2theta23 = hsinsq2theta1->GetXaxis()->GetBinCenter(ix);
    double sinsq2theta23 = hsinsq2theta1->GetXaxis()->GetBinLowEdge(ix);
    //sinsq2theta23 = TMath::Min(sinsq2theta23, 1.00);
    double s23sq_octant1 = (0.5 * (1 - sqrt(1 - sinsq2theta23)));
    double s23sq_octant2 = 1 - (0.5 * (1 - sqrt(1 - sinsq2theta23)));
    //cout << sinsq2theta23 << "\t" << s23sq_octant1 << endl;
    axis1o1[ix-1] = s23sq_octant1;
    //cout << "octant1 bin edge " << ix-1 << "\t" << s23sq_octant1 << endl;
    axis1o2[nbinsx1-ix+1] = s23sq_octant2;
    //cout << "octant2 bin edge " << nbinsx1-ix+1 << "\t" << s23sq_octant2 << endl;
  }
  axis1o1[nbinsx1] = 1 - axis1o1[nbinsx1-1];
  //cout << "octant1 bin edge " << nbinsx1 << "\t" << axis1o1[nbinsx1] << endl;
  axis1o2[0] = 1 - axis1o2[1];
  //cout << "octant2 bin edge " << 0 << "\t" << axis1o2[0] << endl;
  for(int ix = 0; ix <= nbinsx1; ix++) 
    cout << "bin edges for octant 1,2: " << "\t" << axis1o1[ix] << "\t" << axis1o2[ix] << endl;

  //make an array of sinsqtheta values from sinsq2theta axis
  const int nbinsx2  = hsinsq2theta2->GetNbinsX();
  const int nbinsy2  = hsinsq2theta2->GetNbinsY();
  cout << nbinsx2 << "\t" << nbinsy2 << " bins" << endl;
  const double binminy2 = hsinsq2theta2->GetYaxis()->GetXmin();
  const double binmaxy2 = hsinsq2theta2->GetYaxis()->GetXmax();
  double axis2o1[100] = {0};
  double axis2o2[100] = {0};
  for(int ix = 1; ix <= nbinsx2; ix++) {
    //double sinsq2theta23 = hsinsq2theta2->GetXaxis()->GetBinCenter(ix);
    double sinsq2theta23 = hsinsq2theta2->GetXaxis()->GetBinLowEdge(ix);
    //sinsq2theta23 = TMath::Min(sinsq2theta23, 1.00);
    double s23sq_octant1 = (0.5 * (1 - sqrt(1 - sinsq2theta23)));
    double s23sq_octant2 = 1 - (0.5 * (1 - sqrt(1 - sinsq2theta23)));
    //cout << sinsq2theta23 << "\t" << s23sq_octant1 << endl;
    axis2o1[ix-1] = s23sq_octant1;
    //cout << "octant1 bin edge " << ix-1 << "\t" << s23sq_octant1 << endl;
    axis2o2[nbinsx2-ix+1] = s23sq_octant2;
    //cout << "octant2 bin edge " << nbinsx2-ix+1 << "\t" << s23sq_octant2 << endl;
  }
  axis2o1[nbinsx2] = 1 - axis2o1[nbinsx2-1];
  //cout << "octant1 bin edge " << nbinsx2 << "\t" << axis2o1[nbinsx1] << endl;
  axis2o2[0] = 1 - axis2o2[1];
  //cout << "octant2 bin edge " << 0 << "\t" << axis2o2[0] << endl;
  for(int ix = 0; ix <= nbinsx2; ix++) 
    cout << "bin edges for octant 1,2: " << "\t" << axis2o1[ix] << "\t" << axis2o2[ix] << endl;

  //make a histogram with the sinsqtheta x axis octant 1 & fill it
  TH2D * hsinsqtheta1o1_68 = new TH2D("","", nbinsx1, axis1o1, nbinsy1, binminy1, binmaxy1);
  TH2D * hsinsqtheta1o1_90 = new TH2D("","", nbinsx1, axis1o1, nbinsy1, binminy1, binmaxy1);
  hsinsqtheta1o1_68->SetContour(1);
  hsinsqtheta1o1_68->SetContourLevel(0,2.30);
  hsinsqtheta1o1_68->SetLineStyle(3);
  hsinsqtheta1o1_90->SetContour(1);
  hsinsqtheta1o1_90->SetContourLevel(0,4.61);
  hsinsqtheta1o1_90->SetLineStyle(1);
  hsinsqtheta1o1_68->SetLineColor(kBlack);
  hsinsqtheta1o1_90->SetLineColor(kBlack);
  for(int ix = 1; ix <= nbinsx1; ix++) {
    for(int iy = 1; iy <= nbinsy1; iy++) {
      double dchisq = hsinsq2theta1->GetBinContent(ix,iy);
      hsinsqtheta1o1_68->SetBinContent(ix,iy, dchisq);
      hsinsqtheta1o1_90->SetBinContent(ix,iy, dchisq);
    }//iy
  }//ix

  //make a histogram with the sinsqtheta x axis octant 2 & fill it
  TH2D * hsinsqtheta1o2_68 = new TH2D("","", nbinsx1, axis1o2, nbinsy1, binminy1, binmaxy1);
  TH2D * hsinsqtheta1o2_90 = new TH2D("","", nbinsx1, axis1o2, nbinsy1, binminy1, binmaxy1);
  hsinsqtheta1o2_68->SetContour(1);
  hsinsqtheta1o2_68->SetContourLevel(0,2.30);
  hsinsqtheta1o2_68->SetLineStyle(3);
  hsinsqtheta1o2_90->SetContour(1);
  hsinsqtheta1o2_90->SetContourLevel(0,4.61);
  hsinsqtheta1o2_90->SetLineStyle(1);
  hsinsqtheta1o2_68->SetLineColor(kRed);
  hsinsqtheta1o2_90->SetLineColor(kRed);
  for(int ix = 1; ix <= nbinsx1; ix++) {
    for(int iy = 1; iy <= nbinsy1; iy++) {
      double dchisq = hsinsq2theta1->GetBinContent(ix,iy);
      hsinsqtheta1o2_68->SetBinContent(nbinsx1-ix+1,iy, dchisq);
      hsinsqtheta1o2_90->SetBinContent(nbinsx1-ix+1,iy, dchisq);
    }//iy
  }//ix

  //make a histogram with the sinsqtheta x axis octant 1 & fill it
  TH2D * hsinsqtheta2o1_68 = new TH2D("","", nbinsx2, axis2o1, nbinsy2, binminy2, binmaxy2);
  TH2D * hsinsqtheta2o1_90 = new TH2D("","", nbinsx2, axis2o1, nbinsy2, binminy2, binmaxy2);
  hsinsqtheta2o1_68->SetContour(1);
  hsinsqtheta2o1_68->SetContourLevel(0,2.30);
  hsinsqtheta2o1_68->SetLineStyle(3);
  hsinsqtheta2o1_90->SetContour(1);
  hsinsqtheta2o1_90->SetContourLevel(0,4.61);
  hsinsqtheta2o1_90->SetLineStyle(1);
  hsinsqtheta2o1_68->SetLineColor(kRed);
  hsinsqtheta2o1_90->SetLineColor(kRed);
  for(int ix = 1; ix <= nbinsx2; ix++) {
    for(int iy = 1; iy <= nbinsy2; iy++) {
      double dchisq = hsinsq2theta2->GetBinContent(ix,iy);
      hsinsqtheta2o1_68->SetBinContent(ix,iy, dchisq);
      hsinsqtheta2o1_90->SetBinContent(ix,iy, dchisq);
    }//iy
  }//ix

  //make a histogram with the sinsqtheta x axis octant 2 & fill it
  TH2D * hsinsqtheta2o2_68 = new TH2D("","", nbinsx2, axis2o2, nbinsy2, binminy2, binmaxy2);
  TH2D * hsinsqtheta2o2_90 = new TH2D("","", nbinsx2, axis2o2, nbinsy2, binminy2, binmaxy2);
  hsinsqtheta2o2_68->SetContour(1);
  hsinsqtheta2o2_68->SetContourLevel(0,2.30);
  hsinsqtheta2o2_68->SetLineStyle(3);
  hsinsqtheta2o2_90->SetContour(1);
  hsinsqtheta2o2_90->SetContourLevel(0,4.61);
  hsinsqtheta2o2_90->SetLineStyle(1);
  hsinsqtheta2o2_68->SetLineColor(kBlack);
  hsinsqtheta2o2_90->SetLineColor(kBlack);
  for(int ix = 1; ix <= nbinsx2; ix++) {
    for(int iy = 1; iy <= nbinsy2; iy++) {
      double dchisq = hsinsq2theta2->GetBinContent(ix,iy);
      hsinsqtheta2o2_68->SetBinContent(nbinsx2-ix+1,iy, dchisq);
      hsinsqtheta2o2_90->SetBinContent(nbinsx2-ix+1,iy, dchisq);
      //cout << ix << "\t" << iy << "\t" << dchisq << endl
      //<< hsinsqtheta2o2_68->GetXaxis()->GetBinCenter(ix) << "\t"
      //<< hsinsqtheta2o2_68->GetYaxis()->GetBinCenter(iy) << endl;
    }//iy
    cout << ix << "\t" << nbinsx2-ix+1 <<  "\t" << hsinsqtheta2o2_68->GetXaxis()->GetBinCenter(ix) << endl;
  }//ix

  //best fit markers
  double bestfit_mixo1 = 1.000;
  double bestfit_mixo2 = 0.999;
  double bestfit_dmo1  = 0.002443;
  double bestfit_dmo2  = 0.002441;
  double bestfit_sinsqtheta1o1 = (0.5 * (1 - sqrt(1 - bestfit_mixo1)));
  double bestfit_sinsqtheta1o2 = 1 - (0.5 * (1 - sqrt(1 - bestfit_mixo1)));
  double bestfit_sinsqtheta2o1 = (0.5 * (1 - sqrt(1 - bestfit_mixo2)));
  double bestfit_sinsqtheta2o2 = 1 - (0.5 * (1 - sqrt(1 - bestfit_mixo2)));
  TMarker bf1o1(bestfit_sinsqtheta1o1, bestfit_dmo1, 25);
  TMarker bf1o2(bestfit_sinsqtheta1o2, bestfit_dmo1, 29);
  TMarker bf2o1(bestfit_sinsqtheta2o1, bestfit_dmo2, 29);
  TMarker bf2o2(bestfit_sinsqtheta2o2, bestfit_dmo2, 29);
  TMarker bftheta(0.5139, 0.002441, 29);
  bf2o1.SetMarkerColor(kRed);
  bf2o2.SetMarkerColor(kRed);
  bftheta.SetMarkerColor(kBlue);

  //osc param space
  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.35;
  double mixing23_max   = 0.68;
  double mixing23_step  = 0.001;
  double dm23sq_min     = 2.1E-3;
  double dm23sq_max     = 3.1E-3;
  double dm23sq_step    = 5E-5;
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
  //osc_param_space->GetXaxis()->SetTitleOffset(1.2);
  osc_param_space->GetYaxis()->SetTitleOffset(1.4);

  //draw
  TCanvas * c = new TCanvas();
  c->SetLeftMargin(0.18);
  osc_param_space->Draw();
  hsinsqtheta1o1_90->Draw("cont3same");
  hsinsqtheta1o1_68->Draw("cont3same");
  //hsinsqtheta1o2_90->Draw("cont3same");
  //hsinsqtheta1o2_68->Draw("cont3same");
  //hsinsqtheta2o1_90->Draw("cont3same");
  //hsinsqtheta2o1_68->Draw("cont3same");
  hsinsqtheta2o2_90->SetLineColor(kRed);
  hsinsqtheta2o2_68->SetLineColor(kRed);
  hsinsqtheta2o2_90->Draw("cont3same");
  hsinsqtheta2o2_68->Draw("cont3same");
  hsinsqtheta3_90->Draw("cont3same");
  hsinsqtheta3_68->Draw("cont3same");
  if(draw_DM) {
    hsinsqtheta4_90->Draw("cont3same");
    hsinsqtheta4_68->Draw("cont3same");
  }
  bf1o1.Draw();
  //bf1o2.Draw();
  //bf2o1.Draw();
  bf2o2.Draw();
  bftheta.Draw();
  bf1o1.Print();
  bf2o2.Print();
  TLine l1(0.5,0.0021,0.5,0.0031);
  l1.SetLineStyle(2);
  l1.Draw();
  TLegend * leg1 = new TLegend(0.25,0.75,0.55,0.9);
  leg1->SetFillStyle(0);
  leg1->SetLineColor(kWhite);
  leg1->AddEntry(hsinsqtheta1o1_90, "Octant 1", "l");
  leg1->AddEntry(hsinsqtheta2o2_90, "Octant 2", "l");
  leg1->AddEntry(hsinsqtheta3_90, "sin^{2}#theta_{23}", "l");
  if(draw_DM)
    leg1->AddEntry(hsinsqtheta3_90, "sin^{2}#theta_{23}-DMSQ", "l");
  leg1->Draw();
  c->SaveAs("contour_compare_octants_sinsqtheta.pdf");
  c->SaveAs("contour_compare_octants_sinsqtheta.eps");


  //draw only conversion contours
  osc_param_space->Draw();
  hsinsqtheta3_90->SetLineColor(kBlack);
  hsinsqtheta3_68->SetLineColor(kBlack);
  hsinsqtheta1o1_90->SetLineColor(kBlack);
  hsinsqtheta1o1_68->SetLineColor(kBlack);
  hsinsqtheta2o2_90->SetLineColor(kBlack);
  hsinsqtheta2o2_68->SetLineColor(kBlack);
  bftheta.SetMarkerColor(kBlack);
  hsinsqtheta1o1_90->Draw("cont3same");
  hsinsqtheta2o2_90->Draw("cont3same");
  hsinsqtheta1o1_68->Draw("cont3same");
  hsinsqtheta2o2_68->Draw("cont3same");
  bftheta.Draw();
  delete leg1;
  leg1 = new TLegend(0.25,0.75,0.55,0.9);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.06);
  leg1->AddEntry(hsinsqtheta3_68, "T2K Run 1-3 68% C.L.", "l");
  leg1->AddEntry(hsinsqtheta3_90, "T2K Run 1-3 90% C.L.", "l");
  leg1->AddEntry(&bftheta, "T2K Run 1-3 best fit", "p");
  leg1->Draw();
  c->SaveAs("contour_2012_run123_sinsqtheta_conversion.pdf");

  //draw only special fit contour
  osc_param_space->Draw();
  hsinsqtheta3_90->Draw("cont3same");
  hsinsqtheta3_68->Draw("cont3same");
  bftheta.Draw();
  delete leg1;
  leg1 = new TLegend(0.25,0.75,0.55,0.9);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.06);
  leg1->AddEntry(hsinsqtheta3_68, "T2K Run 1-3 68% C.L.", "l");
  leg1->AddEntry(hsinsqtheta3_90, "T2K Run 1-3 90% C.L.", "l");
  leg1->AddEntry(&bftheta, "T2K Run 1-3 best fit", "p");
  leg1->Draw();
  c->SaveAs("contour_2012_run123_sinsqtheta_fit.pdf");
}
