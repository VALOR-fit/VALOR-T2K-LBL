void oned_fc_contour_from_twod_dchisq_surfaces(
  const char * dchisq_crit_filename = "2012results/FC/delta_chisq_critical_2all.root",
  const char * dchisq_filename      = "2012results/DATA/2all/fitresults_3.root",
  const char * dchisq_crit_histname = "hDeltaChi2Surface90",
  const char * dchisq_histname      = "hdchisq",
  double       const_dchisq_level   = 4.61
)
{
  //load dchisq_crit surface
  TFile f_dchisq_crit(dchisq_crit_filename);
  TH2D * h_dchisq_crit = 0;
  f_dchisq_crit.GetObject(dchisq_crit_histname,h_dchisq_crit);
  if(!h_dchisq_crit) {
    cerr << "dchisq_crit histogram not found. Exiting..." << endl;
    exit(1);
  }
  //load dchisq surface
  TFile f_dchisq(dchisq_filename);
  TH2D * h_dchisq = 0;
  f_dchisq.GetObject(dchisq_histname,h_dchisq);
  if(!h_dchisq) {
    cerr << "dchisq histogram not found. Exiting..." << endl;
    exit(1);
  }
  //initialise the 1D histograms
  TH1D * h_dchisq_x       = new TH1D("","", h_dchisq->GetNbinsX(), h_dchisq->GetXaxis()->GetXmin(), h_dchisq->GetXaxis()->GetXmax());
  TH1D * h_contour_x      = new TH1D("","", h_dchisq->GetNbinsX(), h_dchisq->GetXaxis()->GetXmin(), h_dchisq->GetXaxis()->GetXmax());
  TH1D * h_contour2_x     = new TH1D("","", h_dchisq->GetNbinsX(), h_dchisq->GetXaxis()->GetXmin(), h_dchisq->GetXaxis()->GetXmax());
  TH1D * h_constcontour_x = new TH1D("","", h_dchisq->GetNbinsX(), h_dchisq->GetXaxis()->GetXmin(), h_dchisq->GetXaxis()->GetXmax());
  TH1D * h_crit_x         = new TH1D("","", h_dchisq->GetNbinsX(), h_dchisq->GetXaxis()->GetXmin(), h_dchisq->GetXaxis()->GetXmax());
  TH1D * h_crit2_x        = new TH1D("","", h_dchisq->GetNbinsX(), h_dchisq->GetXaxis()->GetXmin(), h_dchisq->GetXaxis()->GetXmax());
  h_dchisq_x      ->GetXaxis()->SetRangeUser(0.9,1);
  h_contour_x     ->GetXaxis()->SetRangeUser(0.9,1);
  h_contour2_x    ->GetXaxis()->SetRangeUser(0.9,1);
  h_constcontour_x->GetXaxis()->SetRangeUser(0.9,1);
  h_crit_x        ->GetXaxis()->SetRangeUser(0.9,1);
  h_crit2_x       ->GetXaxis()->SetRangeUser(0.9,1);
  //loop over x axis
  for(int ix = 1; ix <= h_dchisq->GetNbinsX(); ix++) {
    double min_dchisq = 99999;
    double crit_dchisq = 99999;
    //margnialise over y in dchisq & find the corresponding dchisq_crit
    for(int iy = 1; iy <= h_dchisq->GetNbinsY(); iy++) {
      double current_dchisq = h_dchisq->GetBinContent(ix,iy);
      if(current_dchisq < min_dchisq) {
	min_dchisq = current_dchisq;
	int i_crit = h_dchisq_crit->FindBin(h_dchisq->GetXaxis()->GetBinCenter(ix), h_dchisq->GetYaxis()->GetBinCenter(iy)); 
	crit_dchisq = h_dchisq_crit->GetBinContent(i_crit);
      }
    }//iy
    //try an alternative way of getting the dchisq_crit (marganilise over that surface too)
    double crit_dchisq2 = 99999;
    for(int iy = 1; iy <= h_dchisq_crit->GetNbinsY(); iy++) {
      double current_dchisq_crit = h_dchisq_crit->GetBinContent(h_dchisq_crit->GetXaxis()->FindBin(h_dchisq->GetXaxis()->GetBinCenter(ix)),iy);
      //cout << current_dchisq_crit << endl;
      if((current_dchisq_crit < crit_dchisq2) && (current_dchisq_crit > 0.1)) {
	crit_dchisq2 = current_dchisq_crit;
      }
    }//iy

    //fill the histograms
    h_dchisq_x->SetBinContent(ix, min_dchisq);
    h_crit_x->SetBinContent(ix, crit_dchisq);
    h_crit2_x->SetBinContent(ix, crit_dchisq2);
    h_contour_x->SetBinContent(ix, min_dchisq - crit_dchisq);
    h_contour2_x->SetBinContent(ix, min_dchisq - crit_dchisq2);
    h_constcontour_x->SetBinContent(ix, min_dchisq - const_dchisq_level);

    //write out bins contained with the contour
    if((min_dchisq - crit_dchisq) < 0)
      cout << h_contour_x->GetBinCenter(ix) << "\t is contained within the contour" << endl;
    if((min_dchisq - crit_dchisq2) < 0 && crit_dchisq2 < 99998)
      cout << h_contour_x->GetBinCenter(ix) << "\t is contained within the alternative contour " 
	   << min_dchisq << "\t" << crit_dchisq2 << "\t" << min_dchisq - crit_dchisq2 << endl;
    if((min_dchisq - const_dchisq_level) < 0)
      cout << h_contour_x->GetBinCenter(ix) << "\t is contained within the constant contour" << endl;
  }//ix
  //Draw
  TCanvas * c = new TCanvas();
  TLine l0(0,0,1,0);
  h_dchisq_x->Draw();
  c->SaveAs("1d_dchisq.pdf");
  h_contour_x->Draw();
  l0.Draw();
  c->SaveAs("1d_contour.pdf");
  h_contour2_x->Draw();
  l0.Draw();
  c->SaveAs("1d_contour2.pdf");
  h_constcontour_x->Draw();
  l0.Draw();
  c->SaveAs("1d_constcontour.pdf");
  h_crit_x->Draw();
  c->SaveAs("1d_crit.pdf");
  h_crit2_x->Draw();
  c->SaveAs("1d_crit2.pdf");
}
