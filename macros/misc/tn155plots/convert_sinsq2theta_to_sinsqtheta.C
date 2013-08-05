void convert_sinsq2theta_to_sinsqtheta
(
 const char * outfilename = "nuis3fixdcp_normal_sinsq2theta.root",
 const char * infilename  = "2012results/NEWRESULTS/datafit/dchisq_2all.root",
 const char * incontour   = "hdchisq"
 )
{
  cout << "doesn't work yet" << endl;
  return;

  TFile fin(infilename);
  TH2D * hsinsq2theta = 0;
  fin.GetObject(incontour,hsinsq2theta);


  //make an array of sinsqtheta values from sinsq2theta axis
  const int nbinsx1  = hsinsq2theta->GetNbinsX();
  const int nbinsy1  = hsinsq2theta->GetNbinsY();
  cout << nbinsx1 << "\t" << nbinsy1 << " bins" << endl;
  const double binminy1 = hsinsq2theta->GetYaxis()->GetXmin();
  const double binmaxy1 = hsinsq2theta->GetYaxis()->GetXmax();
  double axis1o1[100] = {0};
  double axis1o2[100] = {0};
  for(int ix = 1; ix <= nbinsx1; ix++) {
    //double sinsq2theta23 = hsinsq2theta->GetXaxis()->GetBinCenter(ix);
    double sinsq2theta23 = hsinsq2theta->GetXaxis()->GetBinLowEdge(ix);
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


}
