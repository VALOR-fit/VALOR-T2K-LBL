void plot_theta_contour_from_2theta_contour
(
 const char * filename1 = "2012results/NEWRESULTS/datafit/dchisq_2all.root",
 const char * filename2 = "2012results/DATA/octant2_23only/combination/fitresults_3.root"
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

  //make an array of sinsqtheta values from sinsq2theta axis
  const int nbinsx1  = hsinsq2theta1->GetNbinsX();
  const int nbinsy1  = hsinsq2theta1->GetNbinsY();
  cout << nbinsx1 << "\t" << nbinsy1 << " bins" << endl;
  const double binminy1 = hsinsq2theta1->GetYaxis()->GetXmin();
  const double binmaxy1 = hsinsq2theta1->GetYaxis()->GetXmax();
  double axis1o1[100] = {0};
  double axis1o2[100] = {0};
  for(int ix = 1; ix <= nbinsx1; ix++) {
    double sinsq2theta23 = hsinsq2theta1->GetXaxis()->GetBinCenter(ix);
    double s23sq_octant1 = (0.5 * (1 - sqrt(1 - sinsq2theta23)));
    double s23sq_octant2 = 1 - (0.5 * (1 - sqrt(1 - sinsq2theta23)));
    //cout << sinsq2theta23 << "\t" << s23sq_octant1 << endl;
    axis1o1[ix] = s23sq_octant1;
    axis1o2[nbinsx1-ix] = s23sq_octant2;
    cout << "bin " << nbinsx1-ix << "\t" << s23sq_octant2 << endl;
  }
  axis1o2[nbinsx1] = axis1o2[nbinsx1-1];
  for(int ix = 0; ix < nbinsx1+1; ix++) {
    cout << ix << "\t" << axis1o2[ix] << "\t" << axis1o1[ix] << endl;
  }

  //make an array of sinsqtheta values from sinsq2theta axis
  const int nbinsx2  = hsinsq2theta2->GetNbinsX();
  const int nbinsy2  = hsinsq2theta2->GetNbinsY();
  cout << nbinsx2 << "\t" << nbinsy2 << " bins" << endl;
  const double binminy2 = hsinsq2theta2->GetYaxis()->GetXmin();
  const double binmaxy2 = hsinsq2theta2->GetYaxis()->GetXmax();
  double axis2o1[100] = {0};
  double axis2o2[100] = {0};
  for(int ix = 1; ix <= nbinsx2; ix++) {
    double sinsq2theta23 = hsinsq2theta2->GetXaxis()->GetBinCenter(ix);
    double s23sq_octant1 = (0.5 * (1 - sqrt(1 - sinsq2theta23)));
    double s23sq_octant2 = 1 - (0.5 * (1 - sqrt(1 - sinsq2theta23)));
    //cout << sinsq2theta23 << "\t" << s23sq_octant1 << endl;
    axis2o1[ix] = s23sq_octant1;
    axis2o2[nbinsx2-ix] = s23sq_octant2;
    cout << "bin " << nbinsx2-ix << "\t" << s23sq_octant2 << endl;
  }
  axis2o2[nbinsx2] = axis2o2[nbinsx2-1];
  for(int ix = 0; ix < nbinsx2+1; ix++) {
    cout << ix << "\t" << axis2o2[ix] << "\t" << axis2o1[ix] << endl;
  }

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
      hsinsqtheta1o2_68->SetBinContent(nbinsx1-ix,iy, dchisq);
      hsinsqtheta1o2_90->SetBinContent(nbinsx1-ix,iy, dchisq);
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
      hsinsqtheta2o2_68->SetBinContent(nbinsx2-ix,iy, dchisq);
      hsinsqtheta2o2_90->SetBinContent(nbinsx2-ix,iy, dchisq);
      cout << ix << "\t" << iy << "\t" << dchisq << endl
	   << hsinsqtheta2o2_68->GetXaxis()->GetBinCenter(ix) << "\t"
	   << hsinsqtheta2o2_68->GetYaxis()->GetBinCenter(iy) << endl;
    }//iy
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
  TMarker bf1o1(bestfit_sinsqtheta1o1, bestfit_dmo1, 29);
  TMarker bf1o2(bestfit_sinsqtheta1o2, bestfit_dmo1, 29);
  TMarker bf2o1(bestfit_sinsqtheta2o1, bestfit_dmo2, 29);
  TMarker bf2o2(bestfit_sinsqtheta2o2, bestfit_dmo2, 29);
  bf2o1.SetMarkerColor(kRed);
  //bf2o2.SetMarkerColor(kRed);

  //osc param space
  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.15;
  double mixing23_max   = 0.95;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 2E-3;
  double dm23sq_max     = 3E-3;
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
  osc_param_space->GetXaxis()->SetTitle("sin^{2}(2#theta_{23})"); // adjust if mixing param is sinsqtheta rather than sinsq2theta...
  osc_param_space->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
  //osc_param_space->GetXaxis()->SetTitleOffset(1.2);
  osc_param_space->GetYaxis()->SetTitleOffset(1.4);

  //draw
  TCanvas * c = new TCanvas();
  c->SetLeftMargin(0.18);
  osc_param_space->Draw();
  hsinsqtheta1o1_90->Draw("cont3same");
  hsinsqtheta1o1_68->Draw("cont3same");
  hsinsqtheta1o2_90->Draw("cont3same");
  hsinsqtheta1o2_68->Draw("cont3same");
  hsinsqtheta2o1_90->Draw("cont3same");
  hsinsqtheta2o1_68->Draw("cont3same");
  hsinsqtheta2o2_90->Draw("cont3same");
  hsinsqtheta2o2_68->Draw("cont3same");
  bf1o1.Draw();
  bf1o2.Draw();
  bf2o1.Draw();
  bf2o2.Draw();
  bf1o1.Print();
  bf2o2.Print();
  TLine l1(0.5,0.002,0.5,0.003);
  l1.Draw();
  c->SaveAs("contour_compare_octants_sinsqtheta.pdf");
}
