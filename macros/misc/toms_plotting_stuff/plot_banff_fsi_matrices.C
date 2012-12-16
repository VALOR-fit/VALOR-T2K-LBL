void plot_banff_fsi_matrices
(

 )
{
  TFile fb("inputs/official/2012/banff/postfit_banff_v8_osc_marg_MAQE_EBfix.root");
  TFile ff("inputs/official/2012/fsi/FSIandSIError_ICL_SKDetectorBinning.root");

  TCanvas * c = new TCanvas();
  c->SetRightMargin(0.15);

  TMatrixDSym * vfsi = 0;
  ff.GetObject("fsi_cov",vfsi);
  const TMatrixDSym & vvfsi = *vfsi;
  const int sizefsi = vvfsi.GetNcols();
  TH2D * fsi = new TH2D("","",sizefsi, 0, sizefsi,sizefsi, 0, sizefsi);
  for(int i = 1; i <= sizefsi; i++) {
    for(int j = 1; j <= sizefsi; j++) {
      fsi->SetBinContent(i,j,vvfsi[i-1][j-1]);
    }
  }

  TAxis * fsi_xaxis = fsi->GetXaxis();
  fsi_xaxis->SetBinLabel(1,"CCQE E0");
  fsi_xaxis->SetBinLabel(2,"CCQE E1");
  fsi_xaxis->SetBinLabel(3,"CCQE E2");
  fsi_xaxis->SetBinLabel(4,"CCnQE");
  fsi_xaxis->SetBinLabel(5,"CC#nu_e");
  fsi_xaxis->SetBinLabel(6,"NC");
  TAxis * fsi_yaxis = fsi->GetYaxis();
  for(int i = 0; i< sizefsi; i++) {
    fsi_yaxis->SetBinLabel(i+1,Form("f_{%d; t, r}^{fsi}", i));
  }
  
  fsi->Draw("colz");
  c->SaveAs("fsi.eps");




  TMatrixDSym * vbanff = 0;
  fb.GetObject("postfit_cov",vbanff);
  const TMatrixDSym & vvbanff = *vbanff;
  const int sizebanff = vvbanff.GetNcols();
  TH2D * banff = new TH2D("","",sizebanff, 0, sizebanff,sizebanff, 0, sizebanff);
  for(int i = 1; i <= sizebanff; i++) {
    for(int j = 1; j <= sizebanff; j++) {
      banff->SetBinContent(i,j,vvbanff[i-1][j-1]);
    }
  }

  TAxis * banff_xaxis = banff->GetXaxis();
  TAxis * banff_yaxis = banff->GetYaxis();
  banff_xaxis->SetBinLabel(1,"#nu_{#mu} flux");
  banff_xaxis->SetBinLabel(12,"#bar{#nu_{#mu}} flux");
  banff_xaxis->SetBinLabel(17,"M_{A} QE");
  banff_xaxis->SetBinLabel(18,"M_{A} RES");
  banff_xaxis->SetBinLabel(19,"CCQE");
  banff_xaxis->SetBinLabel(22,"CC1#pi");
  banff_xaxis->LabelsOption("v");
  for(int i = 0; i< sizebanff; i++) {
    if(i%2 == 0)
      banff_yaxis->SetBinLabel(i+1,Form("f_{%d; t, r}^{banff}", i));
  }
  
  banff->Draw("colz");
  c->SaveAs("banff.eps");
}
