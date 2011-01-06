//
// ----------------------------------------------------------------------------------------------
// Plot (Nsk_NEUT - Nsk_GENIE) / (Nsk_NEUT + Nsk_GENIE)/2 Nsk as functions of sinsq2theta, dmsq.
// Input files are the ones produced by test_nsk_twk_osc.C.
//
// Costas Andreopoulos
// ----------------------------------------------------------------------------------------------
//

void plot_delta_nsk(const char * genie_file, const char * neut_file)
{
  TFile * fg = new TFile(genie_file,"read");
  TH2D * hg = (TH2D*)fg->Get("nsk_vs_osc_params");

  TFile * fn = new TFile(neut_file,"read");
  TH2D * hn = (TH2D*)fn->Get("nsk_vs_osc_params");

  TH2D * hdiff = new TH2D(*hg);

  for(int i = 1; i <= hdiff->GetXaxis()->GetNbins(); i++) {
    for(int j = 1; j <= hdiff->GetYaxis()->GetNbins(); j++) {
	double nskg = hg->GetBinContent(i,j);
	double nskn = hn->GetBinContent(i,j);
	double diff = 100 * 0.5 * (nskn-nskg) / (nskn+nskg);
	hdiff->SetBinContent(i,j,diff);
    }
  }

  hdiff->Draw("colz");
}
