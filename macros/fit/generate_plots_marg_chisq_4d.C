
void generate_plots_marg_chisq_4d(){

  double kpi = 3.14159;
  double dcp_min     = -1 * kpi;
  double dcp_max     = kpi;
  double dcp_step    = 2*kpi/90.;
  double mix13_min   = 0.0;
  double mix13_max   = 0.12;
  double mix13_step  = 0.0012;
  double mix23_min   = 0.3;
  double mix23_max   = 0.7;
  double mix23_step  = 0.005;
  double dm23sq_min  = 1.0E-3;
  double dm23sq_max  = 5.0E-3;
  double dm23sq_step = 1E-4;

  // Make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin                                                 
  int dcp_nbins    = 1 + TMath::Nint( (dcp_max    - dcp_min   ) / dcp_step    );
  int mix13_nbins  = 1 + TMath::Nint( (mix13_max  - mix13_min ) / mix13_step  );
  int mix23_nbins  = 1 + TMath::Nint( (mix23_max  - mix23_min ) / mix23_step  );
  int dm23sq_nbins = 1 + TMath::Nint( (dm23sq_max - dm23sq_min) / dm23sq_step );

  //input file
  TFile *_file0 = TFile::Open("marginalised_chisq_4d.root");
  gStyle->SetFillColor(kWhite);
  gStyle->SetOptStat(0);

  TCanvas * c = new TCanvas("dm23sq", "Marginalised chisq for dm23sq", 1200, 700);
  dm23sq_marg_chisq->Draw("geoff");
  double min = dm23sq_marg_chisq->GetMinimum();
  for(int ibin = 1; ibin <= dm23sq_nbins; ibin++) {
    double content = dm23sq_marg_chisq->GetBinContent(ibin);
    dm23sq_marg_chisq->SetBinContent(ibin, content - min);
  }
  dm23sq_marg_chisqmu->Draw("geoff");
  min = dm23sq_marg_chisqmu->GetMinimum();
  for(int ibin = 1; ibin <= dm23sq_nbins; ibin++) {
    double content = dm23sq_marg_chisqmu->GetBinContent(ibin);
    dm23sq_marg_chisqmu->SetBinContent(ibin, content - min);
  }
  dm23sq_marg_chisqe->Draw("geoff");
  min = dm23sq_marg_chisqe->GetMinimum();
  for(int ibin = 1; ibin <= dm23sq_nbins; ibin++) {
    double content = dm23sq_marg_chisqe->GetBinContent(ibin);
    dm23sq_marg_chisqe->SetBinContent(ibin, content - min);
  }
  dm23sq_marg_chisq->Draw("c");
  dm23sq_marg_chisqmu->SetLineColor(kBlue);
  dm23sq_marg_chisqmu->Draw("csame");
  dm23sq_marg_chisqe->SetLineColor(kRed);
  dm23sq_marg_chisqe->Draw("csame");

  TLegend * legend = new TLegend (0.7, 0.8, 1.0, 1.0);
  legend -> AddEntry ( dm23sq_marg_chisq,"total" , "l");
  legend -> AddEntry ( dm23sq_marg_chisqmu,"1Rmu" , "l");
  legend -> AddEntry ( dm23sq_marg_chisqe,"1Re" , "l");
  legend->Draw("same");

  TCanvas * c1 = new TCanvas("dcp", "Marginalised chisq for dcp", 1200, 700);
  dcp_marg_chisq->Draw("geoff");
  double min = dcp_marg_chisq->GetMinimum();
  for(int ibin = 1; ibin <= dcp_nbins; ibin++) {
    double content = dcp_marg_chisq->GetBinContent(ibin);
    dcp_marg_chisq->SetBinContent(ibin, content - min);
  }
  dcp_marg_chisqmu->Draw("geoff");
  min = dcp_marg_chisqmu->GetMinimum();
  for(int ibin = 1; ibin <= dcp_nbins; ibin++) {
    double content = dcp_marg_chisqmu->GetBinContent(ibin);
    dcp_marg_chisqmu->SetBinContent(ibin, content - min);
  }
  dcp_marg_chisqe->Draw("geoff");
  min = dcp_marg_chisqe->GetMinimum();
  for(int ibin = 1; ibin <= dcp_nbins; ibin++) {
    double content = dcp_marg_chisqe->GetBinContent(ibin);
    dcp_marg_chisqe->SetBinContent(ibin, content - min);
  }

  dcp_marg_chisq->Draw("c");
  dcp_marg_chisqmu->SetLineColor(kBlue);
  dcp_marg_chisqmu->Draw("csame");
  dcp_marg_chisqe->SetLineColor(kRed);
  dcp_marg_chisqe->Draw("csame");

  TLegend * legend1 = new TLegend (0.7, 0.8, 1.0, 1.0);
  legend1 -> AddEntry ( dcp_marg_chisq,"total" , "l");
  legend1 -> AddEntry ( dcp_marg_chisqmu,"1Rmu" , "l");
  legend1 -> AddEntry ( dcp_marg_chisqe,"1Re" , "l");
  legend1->Draw("same");


  TCanvas * c2 = new TCanvas("mix23", "Marginalised chisq for mix23", 1200, 700);
  mix23_marg_chisq->Draw("geoff");
  min = mix23_marg_chisq->GetMinimum();
  for(int ibin = 1; ibin <= mix23_nbins; ibin++) {
    double content = mix23_marg_chisq->GetBinContent(ibin);
    mix23_marg_chisq->SetBinContent(ibin, content - min);
  }
  mix23_marg_chisqmu->Draw("geoff");
  min = mix23_marg_chisqmu->GetMinimum();
  for(int ibin = 1; ibin <= mix23_nbins; ibin++) {
    double content = mix23_marg_chisqmu->GetBinContent(ibin);
    mix23_marg_chisqmu->SetBinContent(ibin, content - min);
  }
  mix23_marg_chisqe->Draw("geoff");
  min = mix23_marg_chisqe->GetMinimum();
  for(int ibin = 1; ibin <= mix23_nbins; ibin++) {
    double content = mix23_marg_chisqe->GetBinContent(ibin);
    mix23_marg_chisqe->SetBinContent(ibin, content - min);
  }

  mix23_marg_chisq->Draw("c");
  mix23_marg_chisqmu->SetLineColor(kBlue);
  mix23_marg_chisqmu->Draw("csame");
  mix23_marg_chisqe->SetLineColor(kRed);
  mix23_marg_chisqe->Draw("csame");

  TLegend * legend2 = new TLegend (0.7, 0.8, 1.0, 1.0);
  legend2 -> AddEntry ( mix23_marg_chisq,"total" , "l");
  legend2 -> AddEntry ( mix23_marg_chisqmu,"1Rmu" , "l");
  legend2 -> AddEntry ( mix23_marg_chisqe,"1Re" , "l");
  legend2->Draw("same");

  TCanvas * c3 = new TCanvas("mix13", "Marginalised chisq for mix13", 1200, 700);
  mix13_marg_chisq->Draw("geoff");
  min = mix13_marg_chisq->GetMinimum();
  for(int ibin = 1; ibin <= mix13_nbins; ibin++) {
    double content = mix13_marg_chisq->GetBinContent(ibin);
    mix13_marg_chisq->SetBinContent(ibin, content - min);
  }
  mix13_marg_chisqmu->Draw("geoff");
  min = mix13_marg_chisqmu->GetMinimum();
  for(int ibin = 1; ibin <= mix13_nbins; ibin++) {
    double content = mix13_marg_chisqmu->GetBinContent(ibin);
    mix13_marg_chisqmu->SetBinContent(ibin, content - min);
  }
  mix13_marg_chisqe->Draw("geoff");
  min = mix13_marg_chisqe->GetMinimum();
  for(int ibin = 1; ibin <= mix13_nbins; ibin++) {
    double content = mix13_marg_chisqe->GetBinContent(ibin);
    mix13_marg_chisqe->SetBinContent(ibin, content - min);
  }

  mix13_marg_chisq->Draw("c");
  mix13_marg_chisqmu->SetLineColor(kBlue);
  mix13_marg_chisqmu->Draw("csame");
  mix13_marg_chisqe->SetLineColor(kRed);
  mix13_marg_chisqe->Draw("csame");

  TLegend * legend3 = new TLegend (0.7, 0.8, 1.0, 1.0);
  legend3 -> AddEntry ( mix13_marg_chisq,"total" , "l");
  legend3 -> AddEntry ( mix13_marg_chisqmu,"1Rmu" , "l");
  legend3 -> AddEntry ( mix13_marg_chisqe,"1Re" , "l");
  legend3->Draw("same");

}
