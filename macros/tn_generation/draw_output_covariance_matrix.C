// Plots the output covariance matrix from the fit
// Lorena Escudero, March 2016

void draw_output_covariance_matrix(
  const char * input_file,
  int entry = 0)
{
  const int n_tot_flux = 52;
  const int n_int = 44;
  const int n_det = 0;
  
  const int n_params = n_tot_flux+n_int+n_det;
  
  TFile * tf = TFile::Open(input_file);
  
  valor::PhysParamList * true_params     = 0; // input true params
  valor::FitParamList *  best_fit_params = 0; // best-fit params

  TTree * tree = (TTree*)tf->Get("valor_dune_nd_ntp");
  tree->SetBranchAddress("true_params",        &true_params);
  tree->SetBranchAddress("best_fit_params",    &best_fit_params);
  
  assert(tree->GetEntries()>entry);
  tree->GetEntry(entry);
  
  TMatrixD cov(*(best_fit_params->GetCovMatrix()));
  
  TH2D * h = new TH2D("h","Fractional covariance matrix", n_params,0.,n_params,n_params,0,n_params);

  for(int i=0; i<n_params; i++){
    for(int j=0; j<n_params; j++){
      double value = cov(i,j);
      h->SetBinContent(i+1,j+1,value);
    }
  }

  TCanvas *c0 = new TCanvas("c0","c0",500,500);
  gStyle->SetOptStat(0);
  cov.Draw("colz");

  TCanvas *c1 = new TCanvas("c1","c1",750,750);
  double xy_margin = 0.174;
  c1->SetLogz();
  c1->SetLeftMargin(xy_margin);
  c1->SetRightMargin(0.15);
  c1->SetBottomMargin(xy_margin);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
 
  for (int iaxis = 0; iaxis < 2 ; iaxis++) {
    TAxis * axis = (iaxis ? h->GetYaxis() : h->GetXaxis() );
    axis->SetBinLabel(10 ,"FHC #nu_{#mu}");
    axis->SetBinLabel(30 ,"FHC #bar{#nu}_{#mu}");
    axis->SetBinLabel(42 ,"FHC #nu_{e}");
    axis->SetBinLabel(49 ,"FHC #bar{#nu}_{e}");
    //~ axis->SetBinLabel(63 ,"RHC #nu_{#mu}");
    //~ axis->SetBinLabel(82 ,"RHC #bar{#nu}_{#mu}");
    //~ axis->SetBinLabel(95 ,"RHC #nu_{e}");
    //~ axis->SetBinLabel(101,"RHC #bar{#nu}_{e}");

    axis->SetBinLabel(4+n_tot_flux,"CCQE");
    //~ axis->SetBinLabel(4+n_tot_flux,"#bar{#nu}_{#mu} CCQE");
    axis->SetBinLabel(8+n_tot_flux,"MEC         ");
    //~ axis->SetBinLabel(8+n_tot_flux,"#bar{#nu}_{#mu} MEC");
    axis->SetBinLabel(12+n_tot_flux,"CC1#pi^{#pm}");
    //~ axis->SetBinLabel(12+n_tot_flux,"#bar{#nu}_{#mu} CC1#pi^{#pm}");
    axis->SetBinLabel(18+n_tot_flux,"CC1#pi^{0}");
    //~ axis->SetBinLabel(18+n_tot_flux,"#bar{#nu}_{#mu} CC1#pi^{0}");
    axis->SetBinLabel(21+n_tot_flux,"CC2#pi         ");
    axis->SetBinLabel(26+n_tot_flux,"DIS");
    //~ axis->SetBinLabel(26+n_tot_flux,"#bar{#nu}_{#mu} DIS");
    axis->SetBinLabel(30+n_tot_flux,"CC coh       ");
    axis->SetBinLabel(32+n_tot_flux,"NC");
    //~ axis->SetBinLabel(33+n_tot_flux,"CC #nu_{e}/#nu_{#mu}");
    //~ axis->SetBinLabel(34+n_tot_flux,"MFP #pi");
    //~ axis->SetBinLabel(35+n_tot_flux,"MFP N");
    axis->SetBinLabel(36+n_tot_flux,"#pi FSI");
    //~ axis->SetBinLabel(37+n_tot_flux,"FrInel #pi");
    //~ axis->SetBinLabel(38+n_tot_flux,"FrAbs #pi");
    //~ axis->SetBinLabel(39+n_tot_flux,"FrPiProd #pi");
    axis->SetBinLabel(42+n_tot_flux,"N FSI");
    //~ axis->SetBinLabel(41+n_tot_flux,"FrInel N");
    //~ axis->SetBinLabel(42+n_tot_flux,"FrAbs N");
    //~ axis->SetBinLabel(n_params,"FrPiProd N");
    //~ axis->LabelsOption("v");

    axis->SetTickLength(0.);
  }

  h->GetYaxis()->LabelsOption("v");



  //  h->GetYaxis()->SetLabelSize("0.03");
  h->Draw("colz");
  h->GetZaxis()->SetRangeUser(10E-10,10E-6);
  //  h->Draw("colz");
  
  TLine * lines[100] = {};

  lines[0] = new TLine(n_tot_flux,0,n_tot_flux,n_params);
  lines[1] = new TLine(0,n_tot_flux,n_params,n_tot_flux);

  lines[2] = new TLine(6+n_tot_flux,0,6+n_tot_flux,n_params);
  lines[3] = new TLine(0,6+n_tot_flux,n_params,6+n_tot_flux);

  lines[4] = new TLine(8+n_tot_flux,0,8+n_tot_flux,n_params);
  lines[5] = new TLine(0,8+n_tot_flux,n_params,8+n_tot_flux);

  lines[6] = new TLine(19,0,19,n_params);
  lines[7] = new TLine(0,19,n_params,19);

  lines[34] = new TLine(14+n_tot_flux,0,14+n_tot_flux,n_params);
  lines[35] = new TLine(0,14+n_tot_flux,n_params,14+n_tot_flux);

  lines[8] = new TLine(20+n_tot_flux,0,20+n_tot_flux,n_params);
  lines[9] = new TLine(0,20+n_tot_flux,n_params,20+n_tot_flux);

  lines[10] = new TLine(22+n_tot_flux,0,22+n_tot_flux,n_params);
  lines[11] = new TLine(0,22+n_tot_flux,n_params,22+n_tot_flux);

  lines[12] = new TLine(38,0,38,n_params);
  lines[13] = new TLine(0,38,n_params,38);

  lines[14] = new TLine(28+n_tot_flux,0,28+n_tot_flux,n_params);
  lines[15] = new TLine(0,28+n_tot_flux,n_params,28+n_tot_flux);

  lines[16] = new TLine(30+n_tot_flux,0,30+n_tot_flux,n_params);
  lines[17] = new TLine(0,30+n_tot_flux,n_params,30+n_tot_flux);

  lines[18] = new TLine(32+n_tot_flux,0,32+n_tot_flux,n_params);
  lines[19] = new TLine(0,32+n_tot_flux,n_params,32+n_tot_flux);

  lines[20] = new TLine(45,0,45,n_params);
  lines[21] = new TLine(0,45,n_params,45);

  //~ lines[22] = new TLine(35+n_tot_flux,0,35+n_tot_flux,n_params);
  //~ lines[23] = new TLine(0,35+n_tot_flux,n_params,35+n_tot_flux);

  lines[24] = new TLine(52,0,52,n_params);
  lines[25] = new TLine(0,52,n_params,52);

  lines[26] = new TLine(71,0,71,n_params);
  lines[27] = new TLine(0,71,n_params,71);

  lines[28] = new TLine(90,0,90,n_params);
  lines[29] = new TLine(0,90,n_params,90);

  //~ lines[30] = new TLine(39+n_tot_flux,0,39+n_tot_flux,n_params);
  //~ lines[31] = new TLine(0,39+n_tot_flux,n_params,39+n_tot_flux);

  lines[32] = new TLine(97,0,97,n_params);
  lines[33] = new TLine(0,97,n_params,97);

/*
  lines[36] = new TLine(17,0,17,43);
  lines[37] = new TLine(0,17,43,17);

  lines[38] = new TLine(32,0,32,43);
  lines[39] = new TLine(0,32,43,32);

  lines[40] = new TLine(41,0,41,43);
  lines[41] = new TLine(0,41,43,41);

  lines[44] = new TLine(42,0,42,43);
  lines[45] = new TLine(0,42,43,42);
  */
  //~ lines[42] = new TLine(33+n_tot_flux,0,33+n_tot_flux,n_params);
  //~ lines[43] = new TLine(0,33+n_tot_flux,n_params,33+n_tot_flux);
  
  for (int il = 0 ; il < 100 ; il++) {
    if (lines[il]) {
      // Don't draw any lines off the edges
      if (lines[il]->GetX1() <= n_params && lines[il]->GetY1() <= n_params &&
          lines[il]->GetX2() <= n_params && lines[il]->GetY2() <= n_params) {
        lines[il]->SetLineWidth(2);
        lines[il]->Draw("same");
      }
    }
  }
  
  c1->Print("cov.eps");
  
  TMatrixD corr = valor::util::CovarianceMatrixToCorrelationMatrix(cov);
  TCanvas * corr_canv = new TCanvas("corrcanv","corrcanv",750,750);
  corr_canv->SetLeftMargin(xy_margin);
  corr_canv->SetRightMargin(0.15);
  corr_canv->SetBottomMargin(xy_margin);
  TH2D * corr_hist = new TH2D(*h);
  for(int i=0; i<n_params; i++){
    for(int j=0; j<n_params; j++){
      double value = corr(i,j);
      corr_hist->SetBinContent(i+1,j+1,value);
    }
  }
  const int nshades = 100;
  int red_blue_palette = valor::util::CreateRedToBlueColorPalette(nshades);
  corr_hist->SetContour(nshades);
  
  corr_hist->Draw("colz");
  corr_hist->GetZaxis()->SetRangeUser(-1.,1.);
  
  for (int il = 0 ; il < 100 ; il++) {
    if (lines[il]) {
      // Don't draw any lines off the edges
      if (lines[il]->GetX1() <= n_params && lines[il]->GetY1() <= n_params &&
          lines[il]->GetX2() <= n_params && lines[il]->GetY2() <= n_params) {
        lines[il]->Draw("same");
      }
    }
  }
  
  corr_canv->Print("corr.eps");

}
