void dcp_sensitivity
(
 const int    npoints   = 25,
 const char * filenames = "./*.root", 
 const char * tag       = ""
)
{
  using namespace std;
  using namespace std::string;
  using namespace t2kana;
  using namespace t2kana::conv;

  //set up chain
  TChain *ch = new TChain("fitnt");
  ch->Add(filenames);
  FitParamList *fparm = 0;
  PhysParamList *iparm = 0;
  FitParamList *fcp0 = 0;
  TObjArray *dat_array = 0;
  TObjArray * pdf_array = 0;
  ch->SetBranchAddress("best_fit_params", &fparm);
  ch->SetBranchAddress("input_params",    &iparm);
  ch->SetBranchAddress("data",            &dat_array);
  ch->SetBranchAddress("best_fit_pdf",    &pdf_array);
  ch->SetBranchAddress("best_fit_params_dcp0", &fcp0);

  //book vectors
  vector<double> chisq_bf_distribution  [npoints];
  vector<double> chisq_dcp0_distribution[npoints];
  vector<double> dchisq_distribution    [npoints];
  double boundary = TMath::Pi() + (TMath::Pi()/(npoints - 1.));
  if(!(npoints % 2))
    boundary = TMath::Pi() + (TMath::Pi()/(npoints));
  TH1D * h = new TH1D("","", npoints, -boundary, +boundary);
  h->GetXaxis()->SetTitle("#delta_{CP}");
  h->GetYaxis()->SetTitle("#sqrt{#Delta#chi^{2}}");

  //loop over chain
  int nfits = ch->GetEntries();
  cout << nfits << " total fits in the TChain" << endl;
  for(int i = 0; i < nfits ; i++) {
    ch->GetEntry(i);

    //get values
    double cpin   = iparm->Value(kPhypDeltaCP);
    double cpfit  = fparm->BestFitValue(kPhypDeltaCP);
    double chi    = fparm->GetBestFitChisq();
    double chicp0 = fcp0->GetBestFitChisq();

    //find bin
    int bin = (int)h->GetXaxis()->FindBin(cpin) - 1;
    if(TMath::Abs(cpin - h->GetXaxis()->GetBinCenter(bin+1)) > 1E-4)
      cout << cpin << "\t" << h->GetXaxis()->GetBinCenter(bin+1) << endl;

    //fill vectors
    chisq_bf_distribution  [bin].push_back(chi);
    chisq_dcp0_distribution[bin].push_back(chicp0);
    dchisq_distribution    [bin].push_back(chicp0 - chi);
  }//entries (i)

  //fill histogram
  for(int ip = 0; ip < npoints; ip++) {
    double total_dchisq = 0;
    for(int iv = 0; iv < dchisq_distribution[ip].size(); iv++) {
      total_dchisq += dchisq_distribution[ip][iv];
    }
    if(total_dchisq)
      total_dchisq /= (double)dchisq_distribution[ip].size();
    h->SetBinContent(ip+1, TMath::Sqrt(total_dchisq));
    cout << ip << "\tAverage dchisq is " << total_dchisq << "\tbased on " << dchisq_distribution[ip].size() << " fits" << endl;
  }//points (ip)

  //draw
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas();
  h->Draw("PC");
  c->SaveAs(Form("dcp_sensitivity_%s.pdf",tag));
  c->SaveAs(Form("dcp_sensitivity_%s.eps",tag));
}
