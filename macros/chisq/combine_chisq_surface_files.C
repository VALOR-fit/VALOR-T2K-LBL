//
// ------------------------------------------------------------------------------------------------
// Take results of the chisq surface fit result, and combine them into a single surface
// (hadd doesn't work if there are some failed jobs)
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

void min_bin(TH2D * h);

void combine_chisq_surface_files(
				 const char * filename_tree = "final.root",
				 const char * filename_hist = "final.root",
				 const char * filename_out  = "test.root",
				 int verbose                = 0
)
{
  using namespace t2kana;

  TChain * fitnt = new TChain("fitnt");
  fitnt->Add(filename_tree);
  FitParamList *  fparm = 0;
  int idm32sq, imix23, iexp;
  fitnt->SetBranchAddress("best_fit_params", &fparm);
  fitnt->SetBranchAddress("idm32sq", &idm32sq);
  fitnt->SetBranchAddress("imix23",  &imix23);
  fitnt->SetBranchAddress("iexpt",   &iexp);

  TCanvas * c = new TCanvas("c","",20,20,500,500);

  TFile f(filename_hist,"READ");
  TH2D * chisq_surface = 0;
  f.GetObject("chisq_surface",chisq_surface);
  TH2D * hchisq  = new TH2D(*chisq_surface);
  //TH2D * nevents = new TH2D(*chisq_surface);
  TH2I * nevents = new TH2I("nevents","nevents",
			    chisq_surface->GetNbinsX(), chisq_surface->GetXaxis()->GetXbins()->GetArray(),
			    chisq_surface->GetNbinsY(), chisq_surface->GetYaxis()->GetXbins()->GetArray());
  TH1D * chisq = new TH1D("chisq","chisq",10000,0,1000);
  const int nhist = 1000;
  TH2D * hchisqn[nhist];
  for(int ih = 0; ih < nhist; ih++) {
    hchisqn[ih] = new TH2D(*chisq_surface);
    for(int ix = 1; ix <= hchisq->GetNbinsX(); ix++) {
      for(int iy = 1; iy <= hchisq->GetNbinsY(); iy++) {
	hchisqn[ih]->SetBinContent(ix,iy, 0);
      }//iy
    }//ix
  }//ih

  for(int iev=0; iev<fitnt->GetEntries(); iev++) {
    fitnt->GetEntry(iev);

    if(verbose > 0 || (iev % 10000 == 0)) {
      cout << "*******************************************************" << endl;
      cout << "Entry: # " << iev << endl;
      cout << "*******************************************************" << endl;
      if(verbose > 1) {
	cout << "Fit parameter list: " << endl;
	cout << *fparm;
	cout << endl;
      }
    }//verbosity

    double this_chisq = fparm->GetBestFitChisq();
    if(this_chisq > 200 || this_chisq < 0) {
      cout << this_chisq << endl;
      continue;
    }
    chisq->Fill(this_chisq);

    double current = hchisq->GetBinContent(imix23,idm32sq);
    hchisq ->SetBinContent(imix23,idm32sq, current + this_chisq);
    int current_ev = nevents->GetBinContent(imix23,idm32sq);
    nevents->SetBinContent(imix23,idm32sq, current_ev + 1);

    if(iexp <= 10)
    cout << iexp << "\t" << imix23 << "\t" << idm32sq << "\t" << this_chisq << endl;
    hchisqn[iexp]->SetBinContent(ix,iy, this_chisq);
  }//iev

  //average
  TH2D * hchisq_av = new TH2D(*hchisq);
  for(int ix = 1; ix <= hchisq->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= hchisq->GetNbinsY(); iy++) {
      double this_chisq    = hchisq->GetBinContent(ix,iy);
      double this_events   = (double)nevents->GetBinContent(ix,iy);
      double this_chisq_av = this_chisq / this_events;
      if(this_chisq_av < 40)
	cout << this_chisq << "\t" << nevents->GetBinContent(ix,iy) << "\t" << this_chisq_av << endl;
      hchisq_av->SetBinContent(ix,iy, this_chisq_av);
    }//iy
  }//ix

  //subtract min
  TH2D * hdchisq = new TH2D(*hchisq_av);
  double chisq_min = hchisq_av->GetMinimum();
  min_bin(hchisq_av);
  for(int ix = 1; ix <= hchisq->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= hchisq->GetNbinsY(); iy++) {
      double this_chisq = hchisq_av->GetBinContent(ix,iy);
      hdchisq->SetBinContent(ix,iy, this_chisq - chisq_min);
    }//iy
  }//ix


  TFile fout(filename_out,"RECREATE");
  chisq->Write();
  hchisq->Write("hchisq");
  hchisq_av->Write("hchisq_av");
  double crit[2] = {2.30,4.61};
  hdchisq->SetContour(2,crit);
  hdchisq->Write("hdchisq");
  nevents->Write("nevents");
  for(int ih = 0; ih < nhist; ih++) {
    if(hchisqn[ih]->GetEntries() <= 0)
      continue;
    min_bin(hchisqn[ih]);
    hchisqn[ih]->Write(Form("hdchisq_%d",ih));
  }

  TCanvas * c = new TCanvas();
  hchisqn[10]->Draw("colz");
  c->SaveAs("test.pdf");
}

void min_bin(TH2D * h) {
  double chisq_min = h->GetMaximum();
  int minbin = h->GetMinimumBin();
  int minx, miny, minz;
  h->GetBinXYZ(minbin, minx, miny, minz);
  cout << "average chisq " << chisq_min << " at " << h->GetXaxis()->GetBinCenter(minx) << " " << h->GetYaxis()->GetBinCenter(miny) << endl;

}
