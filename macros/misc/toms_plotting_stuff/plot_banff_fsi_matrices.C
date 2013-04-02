void plot_banff_fsi_matrices
(
 const char * banfffile = "inputs/official/2012/banff/postfit_banff_v8_osc_marg_MAQE_EBfix.root",
 const char * fsifile   = "inputs/official/2012/fsi/FSIandSIError_ICL_SKDetectorBinning.root",
 bool         summary   = false
 )
{
  TFile fb(banfffile);
  TFile ff(fsifile);

  TCanvas * c = new TCanvas();
  if(!summary)
    c->SetRightMargin(0.15);

  // FSI+SI

  TMatrixDSym * vfsi = 0;
  ff.GetObject("fsi_cov",vfsi);
  const TMatrixDSym & vvfsi = *vfsi;
  const int sizefsi = vvfsi.GetNcols();
  TH2D * fsi = new TH2D("","",sizefsi, 0, sizefsi,sizefsi, 0, sizefsi);
  TH1D * fsi_diag_sqrt = new TH1D("","",sizefsi, 0, sizefsi);
  for(int i = 1; i <= sizefsi; i++) {
    fsi_diag_sqrt->SetBinContent(i,TMath::Sqrt(vvfsi[i-1][i-1]));
    for(int j = 1; j <= sizefsi; j++) {
      fsi->SetBinContent(i,j,vvfsi[i-1][j-1]);
    }
  }

  TAxis * fsi_xaxis = fsi->GetXaxis();
  TAxis * fsi_yaxis = fsi->GetYaxis();
  TAxis * fsi_diag_sqrt_xaxis = fsi_diag_sqrt->GetXaxis();
  fsi_diag_sqrt->GetYaxis()->SetTitle("1#sigma fractional error");
/*
  fsi_xaxis->SetBinLabel(1,"CCQE E0");
  fsi_xaxis->SetBinLabel(2,"CCQE E1");
  fsi_xaxis->SetBinLabel(3,"CCQE E2");
  fsi_xaxis->SetBinLabel(4,"CCnQE");
  fsi_xaxis->SetBinLabel(5,"CC#nu_e");
  fsi_xaxis->SetBinLabel(6,"NC");
*/
  for(int i = 0; i< sizefsi; i++) {
    fsi_xaxis->SetBinLabel(i+1,Form("f_{%d; t, r}^{fsi}", i));
    fsi_yaxis->SetBinLabel(i+1,Form("f_{%d; t, r}^{fsi}", i));
    fsi_diag_sqrt_xaxis->SetBinLabel(i+1,Form("f_{%d; t, r}^{fsi}", i));
  }
  
  fsi->Draw("colz");
  c->SaveAs("fsi.eps");

  fsi_diag_sqrt->SetFillColor(4);
  fsi_diag_sqrt->SetBarWidth(0.4);
  fsi_diag_sqrt->SetBarOffset(0.1);
  fsi_diag_sqrt->SetMinimum(0);
  fsi_diag_sqrt->SetMaximum(0.15);
  fsi_diag_sqrt->Draw("bar");
  c->SaveAs("fsi_diag_sqrt.eps");

  // BANFF

  TMatrixDSym * vbanff     = 0;
  TMatrixDSym * vbanff_pre = 0;
  fb.GetObject("postfit_cov",vbanff);
  fb.GetObject("prefit_cov", vbanff_pre);
  const TMatrixDSym & vvbanff     = *vbanff;
  const TMatrixDSym & vvbanff_pre = *vbanff_pre;
  const int sizebanff = vvbanff.GetNcols();
  TH2D * banff = new TH2D("","",sizebanff, 0, sizebanff,sizebanff, 0, sizebanff);
  TH1D * banff_diag_sqrt     = new TH1D("","",sizebanff, 0, sizebanff);
  TH1D * banff_diag_sqrt_pre = new TH1D("","",sizebanff, 0, sizebanff);
  for(int i = 1; i <= sizebanff; i++) {
    banff_diag_sqrt    ->SetBinContent(i,TMath::Sqrt(vvbanff    [i-1][i-1]));
    banff_diag_sqrt_pre->SetBinContent(i,TMath::Sqrt(vvbanff_pre[i-1][i-1]));
    for(int j = 1; j <= sizebanff; j++) {
      banff->SetBinContent(i,j,vvbanff[i-1][j-1]);
    }
  }

  TAxis * banff_xaxis = banff->GetXaxis();
  TAxis * banff_yaxis = banff->GetYaxis();
  TAxis * banff_diag_sqrt_xaxis = banff_diag_sqrt->GetXaxis();
  banff_diag_sqrt->GetYaxis()->SetTitle("1#sigma fractional error");
/*
  banff_xaxis->SetBinLabel(1,"#nu_{#mu} flux");
  banff_xaxis->SetBinLabel(12,"#bar{#nu_{#mu}} flux");
  banff_xaxis->SetBinLabel(17,"M_{A} QE");
  banff_xaxis->SetBinLabel(18,"M_{A} RES");
  banff_xaxis->SetBinLabel(19,"CCQE");
  banff_xaxis->SetBinLabel(22,"CC1#pi");
*/
  banff_xaxis->LabelsOption("v");
  for(int i = 0; i< sizebanff; i++) {
    if(i%2 == 0) {
      banff_xaxis->SetBinLabel(i+1,Form("f_{%d; t, r}^{banff}", i));
      banff_yaxis->SetBinLabel(i+1,Form("f_{%d; t, r}^{banff}", i));
      if(summary)
	banff_diag_sqrt_xaxis->SetBinLabel(i+1,Form("", i));
      else
	banff_diag_sqrt_xaxis->SetBinLabel(i+1,Form("f_{%d; t, r}^{banff}", i));
    }
  }
  
  banff->Draw("colz");
  c->SaveAs("banff.eps");

  banff_diag_sqrt->SetFillColor(4);
  banff_diag_sqrt->SetBarWidth(0.4);
  banff_diag_sqrt->SetBarOffset(0.1);
  banff_diag_sqrt->SetMinimum(0);
  banff_diag_sqrt->SetMaximum(0.5);
  banff_diag_sqrt->Draw("bar");
  banff_diag_sqrt_pre->SetFillColor(38);
  banff_diag_sqrt_pre->SetBarWidth(0.4);
  banff_diag_sqrt_pre->SetBarOffset(0.5);
  banff_diag_sqrt_pre->SetStats(0);
  banff_diag_sqrt_pre->Draw("bar same");

  TLatex l1(0.2,-0.05,"#uparrow #nu_{#mu} flux #rightarrow");
  l1.SetTextSize(0.05);
  TLatex l2(11.2,-0.05,"#uparrow #bar{#nu_{#mu}} flux #rightarrow");
  l2.SetTextSize(0.05);
  TLatex l3(16.2,-0.05,"#uparrow cross section #rightarrow");
  l3.SetTextSize(0.05);
  
  TLegend * leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->SetFillColor(kRed);
  leg->AddEntry(banff_diag_sqrt_pre, "Pre ND280 fit", "f");
  leg->AddEntry(banff_diag_sqrt, "Post ND280 fit", "f");

  if(summary) {
    l1.Draw();
    l2.Draw();
    l3.Draw();
    leg->Draw(); 
  }

  c->Update();

  c->SaveAs("banff_diag_sqrt.eps");
}
