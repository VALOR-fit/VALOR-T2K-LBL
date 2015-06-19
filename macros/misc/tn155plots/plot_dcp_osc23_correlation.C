//plot correlation between Input DELTACP & Fitted 23-param
//uses output of mock_data_fit_pulls_new.C
//
//plot used in replies to comments v1 for tn155

void plot_dcp_osc23_correlation
(
 const char * filename = "valor_fit_bias_study_fixdcp_sys_0.37.root",
 int mode = 1 //1: sinsqtheta, 2: dm32sq
)
{
  TFile f(filename);
  TH2D * h = 0;
  if(mode == 1)
    f.GetObject("h_dcp_sinsqtheta23",h);
  else if(mode == 2)
    f.GetObject("h_dcp_dm32sq",h);


  h->RebinX(20);
  if(mode == 1)
    h->RebinY(10);
  else if(mode == 2)
    h->RebinY(5);

  h->GetXaxis()->SetTitle("Input #delta_{CP}");
  if(mode == 1)
    h->GetYaxis()->SetTitle("Fitted sin^{2}#theta_{23} resiudal");
  else if(mode == 2)
    h->GetYaxis()->SetTitle("Fitted |#Delta m^{2}_{32}| residual");

  h->GetXaxis()->SetRangeUser(-3.5,+3.5);
  h->GetYaxis()->SetRangeUser(+0.2,+0.8);

  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.1);
  h->Draw("colz");
  c->SaveAs("corr.pdf");
  c->SaveAs("corr.eps");
}
