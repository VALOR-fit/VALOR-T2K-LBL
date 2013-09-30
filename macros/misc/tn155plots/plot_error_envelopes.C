void plot_error_envelopes(
			  const char * filename  = "$T2KOSCANA/2013results/SPECTRA/error_envelope/100000/sk_ereco_spectrum_oscillated.root",
			  double       max_ereco = 7
)
{
  TFile f(filename);
  
  TCanvas * c = new TCanvas();
 
  //draw modes errors
  SKDet_FSI_pluserr->GetYaxis()->SetTitle("1 #sigma error");
  SKDet_FSI_pluserr->GetXaxis()->SetRangeUser(0, max_ereco - 0.001);
  SKDet_FSI_pluserr->SetLineColor(kBlack);
  SKDet_FSI_pluserr->Draw();
  BANFF_pluserr->SetLineColor(kBlue);
  BANFF_pluserr->Draw("same");
  unconstrained_xsec_pluserr->SetLineColor(kRed);
  unconstrained_xsec_pluserr->Draw("same");
  Osc_param_pluserr->SetLineColor(kGreen+2);
  Osc_param_pluserr->Draw("same");

  //legend
  TLegend * leg = new TLegend(0.50,0.67,0.95,0.93);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(132);
  leg->AddEntry(SKDet_FSI_pluserr, "SuperK detector + FSI + SI","l");
  leg->AddEntry(BANFF_pluserr, "BANFF (postfit)","l");
  leg->AddEntry(unconstrained_xsec_pluserr, "Uncorrelated XSec","l");
  leg->AddEntry(Osc_param_pluserr, "#it{Non-23-sector} oscillation parameters","l");
  leg->Draw();

  c->SaveAs("error_envelopes_modes.pdf");
  c->SaveAs("error_envelopes_modes.eps");

  //draw all with envelopes
  TCanvas * c2 = new TCanvas();
  c2->Divide(1,2);

  c2->cd(1);
  all_envpluserr->SetLineStyle(2);
  all_envpluserr->GetXaxis()->SetRangeUser(0, max_ereco - 0.001);
  all_envpluserr->Draw();
  all_envminuserr->SetLineStyle(2);
  all_envminuserr->Draw("same");
  nominal->Draw("same");

  c2->cd(2);
  all_pluserr->GetYaxis()->SetTitle("1 #sigma error");
  all_pluserr->GetXaxis()->SetRangeUser(0, max_ereco - 0.001);
  all_pluserr->Draw();

  c2->SaveAs("error_envelopes_all.pdf");
  c2->SaveAs("error_envelopes_all.eps");

}
