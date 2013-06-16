// Find the energy resolution at the T2K beam peak energy
// the question came up for the Run1-3 PRL
//
// Input is from SuperkMC run with __save_templates_debug__
void energy_resolution(bool all = false)
{
  TFile f("sk_spectra_banffpst_noosc.root");
  const int nmodes = 28;
  TH2D * templ = 0;
  for(int i = 0; i < nmodes; i++) {
    if(!all)
      if(i != 0 && i != 6) continue; //not numu(bar) CCQE
    TH2D * temp = 0;
    f.GetObject(Form("nsk_tun_0_%d",i),temp);
    if(!templ) {
      templ = new TH2D(*temp);
    }
    else {
      for(int ix = 1; ix <= templ->GetNbinsX(); ix++) {
	for(int iy = 1; iy <= templ->GetNbinsY(); iy++) {
	  double nev_new = temp-> GetBinContent(ix,iy);
	  double nev_cur = templ->GetBinContent(ix,iy);
	  templ->SetBinContent(ix, iy, nev_new+nev_cur);
	}//iy
      }//ix
    }//templ
  }//nmodes

  int lowbin = 18, highbin = 19;
  for(int i = lowbin; i <= highbin + 1; i++)
    cout << templ->GetYaxis()->GetBinLowEdge(i) << endl;

  TH1D * spectrum = templ->ProjectionX("spectrum",lowbin,highbin);

  cout << "Resolution: " << spectrum->GetRMS() << endl;

  TFile fout("eres.root","RECREATE");
  templ->Write();
  spectrum->Write();

  gStyle->SetOptStat(1111);

  TCanvas * c = new TCanvas();
  spectrum->GetXaxis()->SetRangeUser(0,1.999);
  //spectrum->GetYaxis()->SetRangeUser(0,0.05);
  spectrum->Draw();
  if(all)
    c->SaveAs("spectrum_all.pdf");
  else
    c->SaveAs("spectrum_ccqe.pdf");

}
