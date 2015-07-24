void plot_goodness_of_fit
(
 const int id = 3,
 const char * filename = "goodness-of-fit-123good28-run_3.root",
 const int rebin = 8,
 const double overrule_pvalue = -1 /* when histogram has been hadd-ed (since hadd doesn't work for vectors) */
)
{
  TFile f(filename);
  TH1D * h = 0;
  f.GetObject(Form("chisq_%d",id),h);
  TLine * l = 0;
  f.GetObject(Form("line_%d",id), l);
  TVectorD * v = 0;
  f.GetObject("goodness_params", v);

  h->Rebin(rebin);
  h->GetXaxis()->SetRangeUser(0,20);
  h->GetXaxis()->SetTitle("#chi^{2}_{gof}");
  l->SetY2(h->GetMaximum()*1.05);

  const TVectorD & vec = *v;
  TPaveText pt(13,10,19,16,"tr");
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
  c->SaveAs(Form("goodness_%d.eps",id));

  //for PAC
  c->SaveAs(Form("goodness_%d.png",id));
  c->SaveAs(Form("goodness_%d.jpg",id));
  c->SaveAs(Form("goodness_%d.pdf",id));
}
