//
// ------------------------------------------------------------------------------------------------
// Plot the \delta\chi^2 surface from the input fit file.
//
// C.Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void plot_deltachisq_surface(
  const char * fitfile
)
{
  TFile * ffit  = new TFile(fitfile, "read");
  TH2D * hchisqsurf = new TH2D(*chisq_surface);

  double chisqmin = 999999;
  double best_sinsq2theta = 0;
  double best_dmsq        = 0;
  for(int i = 1; i <= hchisqsurf->GetNbinsX(); i++) {
     for(int j = 1; j <= hchisqsurf->GetNbinsY(); j++) {
          double chisq = hchisqsurf->GetBinContent(i,j);
          if(chisq < chisqmin) {
              chisqmin = chisq;
              best_sinsq2theta = hchisqsurf->GetXaxis()->GetBinCenter(i);
              best_dmsq        = hchisqsurf->GetYaxis()->GetBinCenter(j);
          }
     }//i
  }//j
  for(int i = 1; i <= hchisqsurf->GetNbinsX(); i++) {
     for(int j = 1; j <= hchisqsurf->GetNbinsY(); j++) {
          hchisqsurf->SetBinContent(i,j, hchisqsurf->GetBinContent(i,j) - chisqmin);
      }//i
  }//j
  TMarker * best_fit_point = new TMarker(best_sinsq2theta, best_dmsq, 29);
  best_fit_point->SetMarkerSize(1.4);

  hchisqsurf->Draw("colz");
  best_fit_point->Draw();
}

