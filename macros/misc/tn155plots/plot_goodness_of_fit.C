void plot_goodness_of_fit
(
 const char * tag,
 const char * filename = "goodness-of-fit-123good28-run_3.root",
 const int rebin = 8,
 const double overrule_pvalue = -1 /* when histogram has been hadd-ed (since hadd doesn't work for vectors) */
)
{
  TFile f(filename);
  f.ls();
  TH1D * h = 0;
  f.GetObject("gof",h);
  TLine * l = 0;
  f.GetObject("line", l);
  TVectorD * v = 0;
  f.GetObject("goodness_params", v);

  h->Rebin(rebin);
  h->GetXaxis()->SetRangeUser(0,20);
  h->GetXaxis()->SetTitle("#chi^{2}_{gof}");
  l->SetY2(h->GetMaximum()*1.05);

  const TVectorD & vec = *v;
  TPaveText pt(13,10,19,18,"tr");
  pt.SetTextSize(0.06);
  pt.AddText(Form("#chi^{2}_{data} = %.2f", vec[1]));
  if(overrule_pvalue < 0)
    pt.AddText(Form("p-value = %.2f", vec[0]));
  else
    pt.AddText(Form("p-value = %.2f", overrule_pvalue));
  pt.SetShadowColor(0);
  pt.SetFillColor(0);
  pt.SetFillStyle(0);

  TCanvas * c = new TCanvas();
  h->Draw();
  l->Draw("same");
  pt.Draw("");
  c->SaveAs(Form("goodness_%s.eps",tag));

  //for PAC
  c->SaveAs(Form("goodness_%s.png",tag));
  c->SaveAs(Form("goodness_%s.jpg",tag));
  c->SaveAs(Form("goodness_%s.pdf",tag));
  c->SaveAs(Form("goodness_%s.C",tag));
}
