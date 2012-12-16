// run from $T2KOSCANA/2012/SENSITIVITIES

void plot_minos_sensitivity_studies() {

  const int npoints = 15;

  //input points
  TMarker * in[npoints+1] = {0};
  in[1] = new TMarker(1.000, 0.00209, 29); in[1]->SetMarkerSize(1.8); in[1]->SetMarkerColor(1);
  in[2] = new TMarker(1.000, 0.00229, 29); in[2]->SetMarkerSize(1.8); in[2]->SetMarkerColor(1);
  in[3] = new TMarker(1.000, 0.00239, 29); in[3]->SetMarkerSize(1.8); in[3]->SetMarkerColor(1);
  in[4] = new TMarker(1.000, 0.00248, 29); in[4]->SetMarkerSize(1.8); in[4]->SetMarkerColor(1);
  in[5] = new TMarker(1.000, 0.00266, 29); in[5]->SetMarkerSize(1.8); in[5]->SetMarkerColor(1);
  in[6] = new TMarker(0.957, 0.00209, 29); in[6]->SetMarkerSize(1.8); in[6]->SetMarkerColor(1);
  in[7] = new TMarker(0.957, 0.00229, 29); in[7]->SetMarkerSize(1.8); in[7]->SetMarkerColor(1);
  in[8] = new TMarker(0.957, 0.00239, 29); in[8]->SetMarkerSize(1.8); in[8]->SetMarkerColor(1);
  in[9] = new TMarker(0.957, 0.00248, 29); in[9]->SetMarkerSize(1.8); in[9]->SetMarkerColor(1);
  in[10] = new TMarker(0.957, 0.00266, 29); in[10]->SetMarkerSize(1.8); in[10]->SetMarkerColor(1);
  in[11] = new TMarker(0.877, 0.00209, 29); in[11]->SetMarkerSize(1.8); in[11]->SetMarkerColor(1);
  in[12] = new TMarker(0.877, 0.00229, 29); in[12]->SetMarkerSize(1.8); in[12]->SetMarkerColor(1);
  in[13] = new TMarker(0.877, 0.00239, 29); in[13]->SetMarkerSize(1.8); in[13]->SetMarkerColor(1);
  in[14] = new TMarker(0.877, 0.00248, 29); in[14]->SetMarkerSize(1.8); in[14]->SetMarkerColor(1);
  in[15] = new TMarker(0.877, 0.00266, 29); in[15]->SetMarkerSize(1.8); in[15]->SetMarkerColor(1);

  //files
  TFile * f[npoints+1] = {0};
  f[1] = new TFile("minos_awayfrombest_28/1/sensitivity.root");
  f[2] = new TFile("minos_awayfrombest_28/2/sensitivity.root");
  f[3] = new TFile("minos_awayfrombest_28/3/sensitivity.root");
  f[4] = new TFile("minos_awayfrombest_28/4/sensitivity.root");
  f[5] = new TFile("minos_awayfrombest_28/5/sensitivity.root");
  f[6] = new TFile("minos_awayfrombest_28/6/sensitivity.root");
  f[7] = new TFile("minos_awayfrombest_28/7/sensitivity.root");
  f[8] = new TFile("minos_28/sensitivity.root");
  f[9] = new TFile("minos_awayfrombest_28/9/sensitivity.root");
  /*f[10] = new TFile("minos_awayfrombest_28/10/sensitivity.root");
  f[11] = new TFile("minos_awayfrombest_28/11/sensitivity.root");
  f[12] = new TFile("minos_awayfrombest_28/12/sensitivity.root");
  f[13] = new TFile("minos_awayfrombest_28/13/sensitivity.root");
  f[14] = new TFile("minos_awayfrombest_28/14/sensitivity.root");
  f[15] = new TFile("minos_awayfrombest_28/15/sensitivity.root");
  */

  //contours
  TH2D * hdchisq68[npoints+1] = {0}, * hdchisq90[npoints+1] = {0};
  for(int i = 1; i < npoints+1; i++) {
    if(!f[i] || f[i]->IsZombie())
      continue;
    f[i]->GetObject("dchisq_for_cont1_68",hdchisq68[i]);
    f[i]->GetObject("dchisq_for_cont1_90",hdchisq90[i]);
    if(hdchisq68[i])
      hdchisq68[i]->SetLineStyle(2);
  }

  // param space
  // dm23sq, mixing23 range and step
  double mixing23_min   = 0.6;
  double mixing23_max   = 1.;
  double mixing23_step  = 0.005;
  double dm23sq_min     = 1.0E-3;
  double dm23sq_max     = 4.0E-3;
  double dm23sq_step    = 1E-4;
  // make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin
  int mixing23_nbins = 1 + (int) ((mixing23_max - mixing23_min) / mixing23_step);
  int dm23sq_nbins   = 1 + (int) ((dm23sq_max   - dm23sq_min  ) / dm23sq_step  );
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
  osc_param_space->GetYaxis()->SetTitleOffset(1.28);

  //setup canvases
  TCanvas * c[5] = {0};
  for(int ic = 0; ic < 5; ic++) {
    c[ic] = new TCanvas();
    c[ic]->SetLeftMargin(0.15);
    c[ic]->cd();
    osc_param_space->Draw();
  }

  //actually draw
  for(int i = 1; i < npoints+1; i++) {
    //correct canvas
    c[i%5]->cd();
    //set colors
    Color_t col = kBlack;
    if(i <= 5) {
      col = kBlack;
    }
    else if(i <= 10) {
      col = kRed;
    }
    else if(i <= 15) {
      col = kBlue;
    }
    if(hdchisq68[i]) {
      hdchisq68[i]->SetLineColor(col);
      hdchisq68[i]->Draw("cont3same");
    }
    if(hdchisq90[i]) {
      hdchisq90[i]->SetLineColor(col);
      hdchisq90[i]->Draw("cont3same");
    }
    in[i]->SetMarkerColor(col);
    in[i]->Draw("Psame");
  }

  for(int ic = 0; ic < 5; ic++) {
    c[ic]->SaveAs(Form("sensitivity_%d.eps",ic));
  }

}
