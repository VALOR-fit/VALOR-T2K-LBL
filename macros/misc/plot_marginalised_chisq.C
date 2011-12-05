//------------------------------------------------------------------------------------------------
// Plot marginalised \Delta\chi^{2} distributions as function of |Dm^{2}_{23}| or sin^2(2theta23).
// This is the macro used for Figs.25,26,35,36,37,38,39,40 of T2K-TN-64
//
// C.Andreopoulos
//------------------------------------------------------------------------------------------------


TH1D * marginalise_chisq(TH2D * chisqsurf, Option_t * option, int max_bin=-1);

void plot_marginalised_chisq(
  const char * label1,          /* data-set label */
  const char * fitfile1,        /* fitter output file (contains \chi^{2} surface) */
  const char * label2    = 0,   /* ... can superimpose \Delta\chi^{2} for a 2nd data-set */
  const char * fitfile2  = 0,
  const char * label3    = 0,   /* ... can superimpose \Delta\chi^{2} for a 3rd data-set */
  const char * fitfile3  = 0,
  const char * label4    = 0,   /* ... can superimpose \Delta\chi^{2} for a 4th data-set */
  const char * fitfile4  = 0
)

{
  const int nmax = 4;

//int max_dmsq_bin = 51;
  int max_dmsq_bin = -1;

  const char * label   [nmax] = { label1,   label2,   label3,   label4   };
  const char * fitfile [nmax] = { fitfile1, fitfile2, fitfile3, fitfile4 };
  int          lwidth  [nmax] = { 2,        2,        2,        2        };
  int          lstyle  [nmax] = { kSolid,   kDashed,  6,        kDotted  };

  TFile *  f [nmax];
  TH2D  *  chisqsurf [nmax];
  TH1D  *  mrgchisq_dm23sq [nmax];
  TH1D  *  mrgchisq_sinsq2q23 [nmax];

  for(int i=0; i<nmax; i++) {
    f[i] = 0;
    chisqsurf[i] = 0;
    mrgchisq_dm23sq[i]    = 0;
    mrgchisq_sinsq2q23[i] = 0;
    if(fitfile[i]) {
       f[i] = new TFile(fitfile[i],"read");
       TH2D * h2d = (TH2D*) f[i]->Get("chisq_surface");
       chisqsurf[i] = new TH2D(*h2d);
       mrgchisq_dm23sq   [i] = marginalise_chisq(chisqsurf[i], "marginalize_mixing23 + subtract_chisqmin",max_dmsq_bin);
       mrgchisq_sinsq2q23[i] = marginalise_chisq(chisqsurf[i], "marginalize_dm23sq + subtract_chisqmin");
       mrgchisq_dm23sq   [i] -> SetLineWidth (lwidth[i]);
       mrgchisq_sinsq2q23[i] -> SetLineWidth (lwidth[i]);
       mrgchisq_dm23sq   [i] -> SetLineStyle (lstyle[i]);
       mrgchisq_sinsq2q23[i] -> SetLineStyle (lstyle[i]);
    }
  }

  TCanvas * c1 = new TCanvas("c1","",10,10,500,500);
  TH2D * c1frame = (TH2D*)c1->DrawFrame(0.001,0,0.030,70);
  c1frame->GetXaxis()->SetTitle("|#Deltam^{2}_{23}| (eV^{2}/c^{4})");
  c1frame->GetYaxis()->SetTitle("#Delta#chi^{2}");
  for(int i=0; i<nmax; i++) {
    if(mrgchisq_dm23sq[i]) {
       mrgchisq_dm23sq[i]->Draw("csame");
    }
  }
  TLegend * legend1 = new TLegend(0.5, 0.5, 0.9, 0.9);
  legend1->SetBorderSize(0);
  legend1->SetFillColor(0);
  for(int i=0; i<nmax; i++) {
    if(mrgchisq_dm23sq[i]) {
      legend1->AddEntry(mrgchisq_dm23sq[i], label[i], "L");
    }
  }
  legend1->Draw();


  TCanvas * c2 = new TCanvas("c2","",110,110,600,600);
  TH2D * c2frame = (TH2D*)c2->DrawFrame(0.0,0,1.0,70);
  c2frame->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
  c2frame->GetYaxis()->SetTitle("#Delta#chi^{2}");
  for(int i=0; i<nmax; i++) {
    if(mrgchisq_sinsq2q23[i]) {
       mrgchisq_sinsq2q23[i]->Draw("csame");
    }
  }
  TLegend * legend2 = new TLegend(0.5, 0.5, 0.9, 0.9);
  legend2->SetBorderSize(0);
  legend2->SetFillColor(0);
  for(int i=0; i<nmax; i++) {
    if(mrgchisq_sinsq2q23[i]) {
      legend2->AddEntry(mrgchisq_sinsq2q23[i], label[i], "L");
    }
  }
  legend2->Draw();
}
//------------------------------------------------------------------------------------------------
TH1D * marginalise_chisq(TH2D * chisqsurf, Option_t * option, int max_bin)
{
  TH1D * mchisq = 0;

  TString opt = option;

  int    mixing23_nbins = chisqsurf->GetXaxis()->GetNbins();
  double mixing23_min   = chisqsurf->GetXaxis()->GetXmin();
  double mixing23_max   = chisqsurf->GetXaxis()->GetXmax();
  int    dm23sq_nbins        = chisqsurf->GetYaxis()->GetNbins();
  double dm23sq_min          = chisqsurf->GetYaxis()->GetXmin();
  double dm23sq_max          = chisqsurf->GetYaxis()->GetXmax();

  //
  // pick minimum chisq along marginalised parameter
  //
  if(opt.Contains("marginalize_dm23sq"))
  {
     //
     // integrate out Dm23sq
     //
     mchisq = new TH1D("","", mixing23_nbins, mixing23_min, mixing23_max);
     for(int ibin = 1; ibin <= mixing23_nbins; ibin++) {
        double chisq = 999999999;        
        for(int jbin = 1; jbin <= dm23sq_nbins; jbin++) {
            chisq = TMath::Min(chisq, chisqsurf->GetBinContent(ibin,jbin));
        }//jbin
        mchisq->SetBinContent(ibin, chisq);
     }//ibin   
  }
  else
  if(opt.Contains("marginalize_mixing23"))
  {
     //
     // integrate out mixing23
     //
     mchisq = new TH1D("","", dm23sq_nbins, dm23sq_min, dm23sq_max);
     for(int jbin = 1; jbin <= dm23sq_nbins; jbin++) {
        if(max_bin>0 && jbin>max_bin) { mchisq->SetBinContent(jbin, 9999999); continue; }
        double chisq = 999999999;
        for(int ibin = 1; ibin <= mixing23_nbins; ibin++) {
            chisq = TMath::Min(chisq, chisqsurf->GetBinContent(ibin,jbin));
        }//ibin
        mchisq->SetBinContent(jbin, chisq);
     }//jbin
  }


  if(opt.Contains("subtract_chisqmin")) {
    double chisqmin = mchisq->GetMinimum();
     for(int ibin = 1; ibin <= mchisq->GetNbinsX(); ibin++) {
       double chisq = mchisq->GetBinContent(ibin);
       chisq -= chisqmin;
       chisq = TMath::Max(0., chisq);
       mchisq->SetBinContent(ibin, chisq);
     }
  }

  return mchisq;

}
//------------------------------------------------------------------------------------------------
