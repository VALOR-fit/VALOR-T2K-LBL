void set_contour_style(TGraph *contour, Int_t color, Int_t width, Int_t style);

void plot_ext_contours_and_sensitivity_grid(int mode = 1 /* 1 - draw sensitivity. 2 - draw star plot */ )
{
  if(mode != 1 && mode != 2) {
    cout << "Invalid mode" << endl;
    return;
  }

  double sinsq2theta23_min   = 0.7; // 0.0;
  double sinsq2theta23_max   = 1.0;
  double sinsq2theta23_step  = 0.01;
  double dm23sq_min          = 1.0E-3; //1E-3;
  double dm23sq_max          = 5.0E-3; //6E-3;
  double dm23sq_step         = 1E-4;


  TCanvas * c = new TCanvas("c","",20,20,400,400);
  c->SetLeftMargin(0.22);

  TH2F * hframe = 0;
  if(mode == 1)
    hframe = (TH2F*)c->DrawFrame(0.8,0.0015,1,0.0040);
  else if(mode == 2)
    hframe = (TH2F*)c->DrawFrame(0.75,0.0015,1,0.0040);
  hframe->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
  hframe->GetYaxis()->SetTitle("|#Deltam_{32}^{2}| (eV^{2}/c^{4})");
//hframe->GetXaxis()->SetLimits(0.,1.002);
//hframe->GetYaxis()->SetLimits(0.001,0.004);
  hframe->GetXaxis()->SetTitleOffset(1.2);
  hframe->GetYaxis()->SetTitleOffset(1.8);
  hframe->GetXaxis()->CenterTitle();
  hframe->GetYaxis()->CenterTitle();
  hframe->Draw();
  hframe->Draw();

  const char * extcontourfile = Form("%s/inputs/contours/contours_nont2k.root", gSystem->Getenv("T2KOSCANA"));
  TFile f(extcontourfile,"read");
  TGraph * minos2012_90cl = new TGraph(* (TGraph*)f.Get("minos2012_combined_90CL") );
  TGraph * skle2012_90cl  = new TGraph(* (TGraph*)f.Get("sk2012_le2nu_90CL")       );
  TGraph * skzen2012_90cl = new TGraph(* (TGraph*)f.Get("sk2012_zenith3nu_90CL")   );
  minos2012_90cl->SetName("minos2012_90cl");
  skle2012_90cl->SetName("skle2012_90cl");
  skzen2012_90cl->SetName("skzen2012_90cl");
  set_contour_style(minos2012_90cl, kBlack,     1, kSolid);
  set_contour_style(skle2012_90cl,  kMagenta-3, 1, kSolid);
  set_contour_style(skzen2012_90cl, kGreen+2,   1, kSolid);
  minos2012_90cl->Draw("l");
  skle2012_90cl ->Draw("l");
  skzen2012_90cl->Draw("l");
  f.Close();

  TLegend * legend = new TLegend(0.25,0.72,0.9,0.9);
  legend->SetFillStyle(0);
  legend->SetTextFont(42);
  
  legend->AddEntry("minos2012_90cl", "MINOS 2012 (10.71#times10^{20} POT + 37.9 kton*yr), 90% CL", "l");
  legend->AddEntry("skzen2012_90cl", "SK1-4 zenith angle, 90% CL", "l");
  legend->AddEntry("skle2012_90cl",  "SK1-4 L/E, 90% CL", "l");
  legend->SetBorderSize(0);
  legend->Draw();

  double minos_best_fit = 0.00239;
  double minos_plus_sig = 0.00009;
  double minos_mins_sig = 0.00010;
  
  const int nsinsq2theta23 = 3;
  const int nnsig = 5;

  double sinsq2theta23[nsinsq2theta23] = { 0.877, 0.957, 1.000 };
  double nsig[nnsig] = { 0, +1, -1, +3, -3 };

  TMarker * m = 0;

  if(mode == 1) {
    for(int i=0; i<nsinsq2theta23; i++) {
      for(int j=0; j<nnsig; j++) {
	double dmsq = minos_best_fit;
        if(nsig[j] > 0) { dmsq += (nsig[j]*minos_plus_sig); }
        if(nsig[j] < 0) { dmsq += (nsig[j]*minos_mins_sig); }
	m = new TMarker(sinsq2theta23[i], dmsq, 20);
        m->SetMarkerSize(1);
        m->Draw();
      }
    }
  }

  else if(mode == 2) {
    const int nstar = 14;
    double starmix[nstar] = {0.8, 0.8, 0.8, 0.9, 0.9, 0.9, 0.95, 0.95, 0.95, 1.0, 1.0, 1.0, 1.0, 1.0};
    double stardm[nstar] = {0.0020, 0.0024, 0.0030, 0.0022, 0.0024, 0.0026, 0.0022, 0.0024, 0.0026, 0.0020, 0.0022, 0.0024, 0.0026, 0.0030};
    for(int i = 0; i < nstar; i++) {
      m = new TMarker(starmix[i], stardm[i], 20);
      m->SetMarkerSize(1);
      m->Draw();
    }
  }

}

void set_contour_style(TGraph *contour, Int_t color, Int_t width, Int_t style)
{
  contour->SetLineWidth(width);
  contour->SetLineColor(color);
  contour->SetLineStyle(style);
}

