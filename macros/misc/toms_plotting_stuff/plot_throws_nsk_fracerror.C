void plot_throws_nsk_fracerror()
{
  TFile f("sk_ereco_spectrum-twk_sys_throw_sinsq2theta[12_0.857,13_0.098,23_1.000],dcp_0.000,dmsq[12_0.00007500,23_0.002400]eV2.root");
  const int nfits = 4;
  TH1D * h[nfits] = {0};
  TH1D * h2[nfits] = {0};
  string name[nfits];
  
  //get objects
  f.GetObject("BANFF_correlatedprefit_nsk",h[0]);
  f.GetObject("BANFF_correlated_nsk;1",h[1]);
  f.GetObject("unconstrained_xsec_nsk",h[2]);
  f.GetObject("SKDet_FSI_nsk",h[3]);

  f.GetObject("BANFF_correlatedprefit_pluserr",h2[0]);
  f.GetObject("BANFF_correlated_pluserr;1",h2[1]);
  f.GetObject("unconstrained_xsec_pluserr",h2[2]);
  f.GetObject("SKDet_FSI_pluserr",h2[3]);

  name[0] = "BANFF pre";
  name[1] = "BANFF post";
  name[2] = "Unconstrained XSec";
  name[3] = "SK detector + FSI + SI";

  //setup drawing
  TCanvas * c = new TCanvas();
  TCanvas * c2 = new TCanvas();
  TLegend * l = new TLegend(0.6,0.7,0.95,0.9);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  TLegend * l2 = new TLegend(0.6,0.7,0.95,0.9);
  l2->SetFillStyle(0);
  l2->SetBorderSize(0);

  //setup colors, titles, binning, legend
  for(int i = 0; i < nfits; i++) {
    Color_t col = kBlack; 
    switch(i) {
    case 0:
      col = kBlack; break;
    case 1:
      col = kRed; break;
    case 2:
      col = kBlue; break;
    case 3:
      col = kGreen+2; break;
    default:
      col = kBlack; break;
    }
    h[i]->SetLineColor(col);
    h2[i]->SetLineColor(col);
    h[i]->Rebin(5);
    l->AddEntry(h[i],name[i].c_str(), "l");
    l2->AddEntry(h2[i],name[i].c_str(), "l");
    h[i]->GetXaxis()->SetRangeUser(0,120);
    h[i]->GetXaxis()->SetTitle("N_{SK}");
    h2[i]->GetXaxis()->SetTitle("E_{reco} / GeV");
    h2[i]->GetYaxis()->SetTitle("Fractional error");
    h2[i]->GetXaxis()->SetRangeUser(0,5.99);
  }

  //actually draw (need to do in a different order)
  c->cd();
  h[1]->Draw();
  h[0]->Draw("same");
  h[2]->Draw("same");
  h[3]->Draw("same");
  l->Draw();
  c->SaveAs("nsk.png");
  c->SaveAs("nsk.jpg");
  c->SaveAs("nsk.eps");
  c->SaveAs("nsk.pdf");

  c2->cd();
  h2[3]->Draw();
  h2[0]->Draw("same");
  h2[1]->Draw("same");
  h2[2]->Draw("same");
  l2->Draw();
  c2->SaveAs("fracerror.png");
  c2->SaveAs("fracerror.jpg");
  c2->SaveAs("fracerror.eps");
  c2->SaveAs("fracerror.pdf");
}
