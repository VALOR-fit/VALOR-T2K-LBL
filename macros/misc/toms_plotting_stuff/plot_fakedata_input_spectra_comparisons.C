void plot_fakedata_input_spectra_comparisons()
{
  TFile * fin[6] = {0};
  TH1D * hin[6] = {0};
  TFile * f2a[6] = {0};
  TH1D * h2a[6] = {0};
  TFile * f28[6] = {0};
  TH1D * h28[6] = {0};
  TH1D * hob[6] = {0};
  TCanvas * c = new TCanvas();


  for(int i = 0; i < 6; i++) {
    hin[i] = 0;
    h2a[i] = 0;
    h28[i] = 0;
    hob[i] = 0;
    fin[i] = new TFile(Form("$T2KOSCANA/2012results/DATASPECTRA/FAKEDATAINPUTS/%d/input_spectra_%d.root",i,i));
    fin[i]->GetObject("sk_1Rmutotal",hin[i]);
    f2a[i] = new TFile(Form("$T2KOSCANA/2012results/DATA/FAKE2all/fitresults_%d.root"));
    f2a[i]->GetObject("bf_spectrum",h2a[i]);
    f28[i] = new TFile(Form("$T2KOSCANA/2012results/DATA/bigger_surface/FAKE28/fitresults_%d.root"));
    f28[i]->GetObject("bf_spectrum",h28[i]);
    f28[i]->GetObject("obs_spectrum",hob[i]);

    TLegend * leg = new TLegend(0.75,0.7,0.95,0.9);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    hob[i]->GetXaxis()->SetTitle("E_{reco}");
    hob[i]->GetXaxis()->SetRangeUser(0,5.9999999);

    h28[i]->SetLineColor(kRed);
    h2a[i]->SetLineColor(kBlue);

    hob[i]->Draw("PE");
    h28[i]->Draw("same");
    hin[i]->Draw("same");
    h2a[i]->Draw("same");

    leg->AddEntry(hob[i], "Fake data", "p");
    leg->AddEntry(hin[i], "Input", "l");
    leg->AddEntry(h28[i], "2+8 best fit", "l");
    leg->AddEntry(h2a[i], "2+N best fit", "l");
    leg->Draw();

    c->SaveAs(Form("fakedata_spectra_comparison_%d.eps",i));

    delete leg;
  }
  

}
