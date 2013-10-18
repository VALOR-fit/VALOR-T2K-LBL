//
// ------------------------------------------------------------------------------------------------
// The nuis4 disappearance fit doesn't work properly in the region where there is a minimum in the chisq at the boundary (±pi)
// This causes the chisq surface to be not as good as the nuis3fixdcp surface
// This macro:
//  takes the result of dcp_fit_correlation.C 1D scans in deltacp
//  reorganises the output pdf, and also adds the best fit point to the plot, & whether it's SIMPLEX or MINUIT
//  can apply various cuts such as only print pdf for specific dm32sq value or when the fit & the scan minima are different (hardcoded currently)
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

void min_bin(TH2D * h);

void dcp_fit_correlation_draw_scans(
			 const char * filename_tree = "run1-4_nuis4_sintheta_normal_ox/combination/final.root",
			 const char * filename_hist = "data_dcp_scan_normal_hierarchy.root",
			 const char * filename_out  = "test2.root",
			 int verbose                = 0
)
{
  using namespace t2kana;
  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsqTheta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  TChain * fitnt = new TChain("fitnt");
  fitnt->Add(filename_tree);
  FitParamList *  fparm = 0;
  int idm32sq, imix23, iexp;
  fitnt->SetBranchAddress("best_fit_params", &fparm);
  fitnt->SetBranchAddress("idm32sq", &idm32sq);
  fitnt->SetBranchAddress("imix23",  &imix23);
  fitnt->SetBranchAddress("iexpt",   &iexp);

  TFile f(filename_hist,"READ");

  int nbinsdm = 0, nbinsmix = 0;
  double mixmin = 1, mixmax = 0, dmmin = 1, dmmax = 0;
  for(int iev=0; iev<fitnt->GetEntries(); iev++) {
    fitnt->GetEntry(iev);
    nbinsdm  = TMath::Max(nbinsdm, idm32sq);
    nbinsmix = TMath::Max(nbinsmix,imix23);
    mixmin   = TMath::Min(mixmin,  fparm->BestFitValue(conv::kPhypMixing23));
    mixmax   = TMath::Max(mixmax,  fparm->BestFitValue(conv::kPhypMixing23));
    dmmin    = TMath::Min(dmmin,   fparm->BestFitValue(conv::kPhypDm23sq));
    dmmax    = TMath::Max(dmmax,   fparm->BestFitValue(conv::kPhypDm23sq));
  }
  const int nbins1 = nbinsdm + 1;
  const int nbins2 = nbinsmix + 1;
  TH1D * hdcp[nbins1][nbins2] = {0};
  TCanvas * c[nbins1][nbins2] = {0};
  TCanvas * ctemp = new TCanvas();
  ctemp->SaveAs("dcp_profiles_draw.pdf[");
 
  cout << "Starting loop over " << fitnt->GetEntries() << " entries" << endl;
  for(int iev=0; iev<fitnt->GetEntries(); iev++) {
    fitnt->GetEntry(iev);
    hdcp[idm32sq][imix23] = 0;
    f.GetObject(Form("mix_%d_dm_%d",imix23,idm32sq), hdcp[idm32sq][imix23]);
    if(!hdcp[idm32sq][imix23])
      continue;

    //if(TMath::Abs(hdcp[idm32sq][imix23]->GetMinimum() - fparm->GetBestFitChisq()) < 1E-3) continue;

    //if(fparm->GetFitTries() != 99) continue;

    c[idm32sq][imix23] = new TCanvas();
    hdcp[idm32sq][imix23]->Draw();
    TLegend * l = new TLegend(0,0,0.3,0.1);
    l->AddEntry(hdcp[idm32sq][imix23], hdcp[idm32sq][imix23]->GetTitle(), "l");
    l->Draw();
    TLegend * leg = new TLegend(0.3,0,0.8,0.1);
    leg->AddEntry(hdcp[idm32sq][imix23], Form("#chi^{2}_{bf}=%f, tries=%d, #delta_{CP}=%f",fparm->GetBestFitChisq(),fparm->GetFitTries(),fparm->BestFitValue(conv::kPhypDeltaCP)), "l");
    leg->SetTextSize(0.03);
    leg->Draw();
    TLine *l1 = new TLine(-conv::kpi,fparm->GetBestFitChisq(),+conv::kpi,fparm->GetBestFitChisq());
    l1->SetLineColor(kRed);
    l1->Draw();
    TLine *l2 = new TLine(fparm->BestFitValue(conv::kPhypDeltaCP),0,fparm->BestFitValue(conv::kPhypDeltaCP),100);
    l2->SetLineColor(kRed);
    l2->Draw();
  }//iev

  //for(int idm32sq = 0; idm32sq < nbinsdm; idm32sq++)
    idm32sq = 9;//this is the transition bin 2.5E-3 in the chisq surface i was studying
    for(int imix23 = 0; imix23 < nbinsmix; imix23++)
      if(c[idm32sq][imix23])
	c[idm32sq][imix23]->SaveAs("dcp_profiles_draw.pdf");


  ctemp->SaveAs("dcp_profiles_draw.pdf]");

}

