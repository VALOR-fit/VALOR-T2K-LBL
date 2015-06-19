//the dchisq surface for the nuis4 fit is not completely smooth.
//this macro is an attempt to correlate failures with things we check (fit time, fit attempts, status codes)

void look_at_failures_in_dchisq_surface(const char * filename = "final.root")
{
  using namespace valor;

  TFile f(filename);
  TTree * fitnt = (TTree*)f.Get("fitnt");
  FitParamList *  fparm = 0;
  fitnt->SetBranchAddress("best_fit_params", &fparm);

  TH2D * h = 0;
  f.GetObject("chisq_surface",h);
  TH2D * h2 = new TH2D(*h);

  for(int ix = 1; ix <= h2->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= h2->GetNbinsY(); iy++) {
      h2->SetBinContent(ix,iy, 0);
    }//iy
  }//ix

  for(int i = 0; i < fitnt->GetEntries(); i++) {
    fitnt->GetEntry(i);
    double mix = fparm->BestFitValue(conv::kPhypMixing23);
    double dm  = fparm->BestFitValue(conv::kPhypDm23sq);
    int    bin = h2->FindBin(mix, dm);
    double chisq = fparm->GetBestFitChisq();
    if(
       //fparm->GetFitTime() > 250
       //||
       //fparm->GetFitTime() > 3
       //||
       fparm->GetStatusCode("migrad_min") != 0
       ||
       fparm->GetStatusCode("migrad_err") != 3
	)
      continue;
    h2->SetBinContent(bin, chisq);
  }//tree entries

  TCanvas * c = new TCanvas();
  h2->Draw("colz");
  c->SaveAs("test1.pdf");
  h->SetLineColor(kRed);
  h->Draw("colz");
  c->SaveAs("test2.pdf");
  h->Draw("cont3");
  c->SaveAs("test3.pdf");
  
}
