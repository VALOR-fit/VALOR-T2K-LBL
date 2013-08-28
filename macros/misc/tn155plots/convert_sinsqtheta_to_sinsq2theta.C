void convert_sinsqtheta_to_sinsq2theta
(
 const char * outfilename = "nuis3fixdcp_normal_sinsq2theta.root",
 const char * infilename  = "2013results/DATA/run1-4_nuis3fixdcp_sintheta_normal/combination/run1-4_nuis3fixdcp_sintheta_normal-run_0-expt_6.root",
 const char * incontour   = "dchisq",
 bool         is_constdchisq = true, // if set, will apply the 2.30,4.61 values to output contourlevel. otherwise, apply 0 (the FC value)
 int          force_octant   = 0     // if 0 will marginalise, otherwise it's 1 or 2
 )
{
  TFile fin(infilename);
  TH2D * hsinsqtheta = 0;
  fin.GetObject(incontour,hsinsqtheta);


  //make an array of sinsqtheta values from sinsq2theta axis
  const int nbinsx  = hsinsqtheta->GetNbinsX();
  const int nbinsy  = hsinsqtheta->GetNbinsY();
  cout << nbinsx << "\t" << nbinsy << " bins" << endl;
  const double binminy = hsinsqtheta->GetYaxis()->GetXmin();
  const double binmaxy = hsinsqtheta->GetYaxis()->GetXmax();
  double axis[1000] = {0};
  int nsinsq2thetabins = -1;
  for(int ix = 1; ix <= nbinsx; ix++) {
    //double s23sq = hsinsqtheta->GetXaxis()->GetBinCenter(ix);
    double s23sq = hsinsqtheta->GetXaxis()->GetBinLowEdge(ix);
    if(s23sq > 0.5)
      break;
    double sinsq2theta23 = 1 - ((1 - (s23sq * 2)) * (1 - (s23sq * 2)));
    axis[ix-1] = sinsq2theta23;
    //cout << s23sq << " goes to new bin edge: " << ix-1 << "\t" << sinsq2theta23 << endl;
    nsinsq2thetabins++;
  }

  axis[nsinsq2thetabins+1] = 1 + TMath::Abs(1 - axis[nsinsq2thetabins]);
  cout << nsinsq2thetabins << " sinsq2theta bins:" << endl;
  for(int ix = 0; ix <= nsinsq2thetabins+1; ix++) 
    cout << "bin edges: " << "\t" << axis[ix] << endl;

  //make a histogram with the sinsq2theta x axis & fill it
  TH2D * hsinsq2theta_68 = new TH2D("hdchisq_68","", nsinsq2thetabins, axis, nbinsy, binminy, binmaxy);
  TH2D * hsinsq2theta_90 = new TH2D("hdchisq_90","", nsinsq2thetabins, axis, nbinsy, binminy, binmaxy);
  hsinsq2theta_68->SetContour(1);
  if(is_constdchisq)
    hsinsq2theta_68->SetContourLevel(0,2.30);
  else
    hsinsq2theta_68->SetContourLevel(0,0.00);
  hsinsq2theta_68->SetLineStyle(3);
  hsinsq2theta_90->SetContour(1);
  if(is_constdchisq)
    hsinsq2theta_90->SetContourLevel(0,4.61);
  else
    hsinsq2theta_90->SetContourLevel(0,0.00);
  hsinsq2theta_90->SetLineStyle(1);
  hsinsq2theta_68->SetLineColor(kBlack);
  hsinsq2theta_90->SetLineColor(kBlack);

  if(force_octant == 0) {
    for(int ix = 1; ix <= nbinsx; ix++) {
      double s23sq_curr = hsinsqtheta->GetXaxis()->GetBinCenter(ix);
      double s23sq_mirr = 0.5;
      if(s23sq_curr < 0.5)
	s23sq_mirr = 0.5 + (TMath::Abs(0.5 - s23sq_curr));
      else
	s23sq_mirr = 0.5 - (TMath::Abs(0.5 - s23sq_curr));
      int ix_alt = hsinsqtheta->GetXaxis()->FindBin(s23sq_mirr);
      //cout << ix << "\t" << ix_alt << endl; continue;
      for(int iy = 1; iy <= nbinsy; iy++) {
	double dchisq     = hsinsqtheta->GetBinContent(ix,    iy);
	double dchisq_alt = hsinsqtheta->GetBinContent(ix_alt,iy);
	//cout << dchisq << "\t" << dchisq_alt << endl;
	double dchisq_min = TMath::Min(dchisq, dchisq_alt);
	hsinsq2theta_68->SetBinContent(ix,iy, dchisq_min);
	hsinsq2theta_90->SetBinContent(ix,iy, dchisq_min);
      }//iy
    }//ix
  }//force_octant == 0
  else if(force_octant == 1) {
    for(int ix = 1; ix <= nbinsx; ix++) {
      double s23sq_curr = hsinsqtheta->GetXaxis()->GetBinCenter(ix);
      if(s23sq_curr > 0.5)
	continue;
      //cout << ix << "\t" << ix_alt << endl; continue;
      for(int iy = 1; iy <= nbinsy; iy++) {
	double dchisq     = hsinsqtheta->GetBinContent(ix,    iy);
	hsinsq2theta_68->SetBinContent(ix,iy, dchisq);
	hsinsq2theta_90->SetBinContent(ix,iy, dchisq);
      }//iy
    }//ix
  }//force_octant == 1
  else if(force_octant == 2) {
    for(int ix = 1; ix <= nbinsx; ix++) {
      double s23sq_curr = hsinsqtheta->GetXaxis()->GetBinCenter(ix);
      if(s23sq_curr < 0.5)
	continue;
      //cout << ix << "\t" << ix_alt << endl; continue;
      for(int iy = 1; iy <= nbinsy; iy++) {
	double dchisq     = hsinsqtheta->GetBinContent(ix,    iy);
	hsinsq2theta_68->SetBinContent(TMath::Abs(ix-nbinsx)+1,iy, dchisq);
	hsinsq2theta_90->SetBinContent(TMath::Abs(ix-nbinsx)+1,iy, dchisq);
      }//iy
      cout << TMath::Abs(ix-nbinsx)+1 << "\t" << s23sq_curr << endl;
    }//ix
  }//force_octant == 2
  else {
    cerr << "Unknown force_octant value: " << force_octant << endl;
    exit(1);
  }

  TFile fout(outfilename,"RECREATE");
  fout.cd();
  hsinsq2theta_68->Write();
  hsinsq2theta_90->Write();
  fout.Close();
}


/* 
check it works!

root2k nuis3fixdcp_normal_sinsq2theta.root 2013results/DATA/run1-4_nuis3fixdcp_sin2theta_octant2_normal/combination/run1-4_nuis3fixdcp_sin2theta_octant2_normal-run_0-expt_6.root 
_file0->cd();
hdchisq_90->Draw("cont3");
hdchisq_68->Draw("cont3same");
_file1->cd();
dchisq_cont68->SetLineColor(kRed);
dchisq_cont68->SetLineStyle(3);
dchisq_cont90->SetLineColor(kRed);
dchisq_cont90->SetLineStyle(3);
dchisq_cont68->Draw("cont3same");
dchisq_cont90->Draw("cont3same");

*/
