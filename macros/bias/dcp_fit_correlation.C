//
// ------------------------------------------------------------------------------------------------
// The nuis4 disappearance fit doesn't work properly in the region where there is a minimum in the chisq at the boundary (±pi)
// This causes the chisq surface to be not as good as the nuis3fixdcp surface
// This macro:
//  prints out & histograms the best fit deltacp for each point on the oscillation surface
//  has a mode to create 1D scans of deltacp at every point on the surface (time consuming)
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

void min_bin(TH2D * h);

void dcp_fit_correlation(
			 const int    npibins       = 100,
			 const char * filename_tree = "run1-4_nuis4_sintheta_normal_ox/combination/final.root",
			 //const char * filename_tree = "run1-4_nuis3fixdcp_sintheta_normal_ox/combination/final.root",
			 const char * filename_hist = "final.root",
			 const char * filename_out  = "test.root",
			 int verbose                = 0,
			 bool scan_dcp              = false
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

  TCanvas * c = new TCanvas("c","",20,20,500,500);

  /*
  TFile f(filename_hist,"READ");
  TH2D * chisq_surface = 0;
  f.GetObject("chisq_surface",chisq_surface);
  TH2D * hchisq  = new TH2D(*chisq_surface);
  TH2I * nevents = new TH2I("nevents","nevents",
			    chisq_surface->GetNbinsX(), chisq_surface->GetXaxis()->GetXbins()->GetArray(),
			    chisq_surface->GetNbinsY(), chisq_surface->GetYaxis()->GetXbins()->GetArray());
  */

  const char * config_file             = "analysis_configuration.xml";
  TString config = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);
  SuperkExpectation * skexp = new SuperkExpectation;
  if(scan_dcp)
    skexp->Read(config.Data());

  const char * xmlfile = "data.xml";
  const int iexperiment = 8;
  bool using_xml_file = true;
  TString xmlfile2 = Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile);
  SuperkObservation * skobs = new SuperkObservation;
  if(scan_dcp)
    skobs->Read(xmlfile2.Data(),iexperiment,using_xml_file);

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
  TH2D * hdcpbf = new TH2D("hdcpbf","",nbinsmix,mixmin,mixmax,nbinsdm,dmmin,dmmax);

  TCanvas * c;
  if(scan_dcp) {
    c = new TCanvas();
    c->SaveAs("dcp_profiles.pdf[");
  }
  TFile fout;
  if(scan_dcp)
    fout.Open(filename_out,"RECREATE");
  
  cout << "Starting loop over " << fitnt->GetEntries() << " entries" << endl;
  for(int iev=0; iev<fitnt->GetEntries(); iev++) {
    fitnt->GetEntry(iev);
    if(!hdcp[idm32sq][imix23] && scan_dcp)
      hdcp[idm32sq][imix23] = new TH1D("",
				       Form("mix %f, dm %f;#delta_{CP}",fparm->BestFitValue(conv::kPhypMixing23),fparm->BestFitValue(conv::kPhypDm23sq)),
				       npibins,-conv::kpi,+conv::kpi);
    if(verbose > 0 || (iev % 10000 == 0)) {
      cout << "*******************************************************" << endl;
      cout << "Entry: # " << iev << endl;
      cout << "*******************************************************" << endl;
      if(verbose > 1) {
	cout << "Fit parameter list: " << endl;
	cout << *fparm;
	cout << endl;
      }
    }//verbosity<< this_chisq << endl;

    cout << fparm->BestFitValue(conv::kPhypDm23sq) << "\t" << fparm->BestFitValue(conv::kPhypDeltaCP) << endl;
    if(TMath::Abs(fparm->BestFitValue(conv::kPhypDeltaCP) < 1E-4) || TMath::Abs(fparm->BestFitValue(conv::kPhypDeltaCP) - conv::kpi) < 1E-4 || TMath::Abs(fparm->BestFitValue(conv::kPhypDeltaCP) + conv::kpi) < 1E-4)
      hdcpbf->SetBinContent(imix23,idm32sq,1);
    hdcpbf->SetBinContent(imix23,idm32sq,fparm->BestFitValue(conv::kPhypDeltaCP));

    if(scan_dcp) {
      fparm->CopyBestFitToCurrValues();

      FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);

      //cout << PhysHypothesis(*fparm) << endl;
      cout << "\t<<<<<  " << fparm->GetBestFitChisq() << endl;
      for(int ix = 1; ix <= hdcp[idm32sq][imix23]->GetNbinsX(); ix++) {
	double deltacp = hdcp[idm32sq][imix23]->GetBinCenter(ix);
	fparm->SetInitValue(conv::kPhypDeltaCP, deltacp);
	skexp->Generate(*fparm);      
	double chisq = fit_inputs->Chisquare();//("poisson",conv::kSKsFCFV1RingMuLike);
	hdcp[idm32sq][imix23]->Fill(deltacp, chisq);
	cout << chisq << "\t";
      }
      cout << endl;

      //save & draw
      hdcp[idm32sq][imix23]->Write(Form("mix_%d_dm_%d",imix23,idm32sq));
      hdcp[idm32sq][imix23]->Draw();
      TLegend * l = new TLegend(0,0,0.3,0.1);
      l->AddEntry(hdcp[idm32sq][imix23], hdcp[idm32sq][imix23]->GetTitle(), "l");
      l->Draw();
      c->SaveAs("dcp_profiles.pdf");
      delete fit_inputs;
    }//scan_dcp
  }//iev

  if(scan_dcp)
    c->SaveAs("dcp_profiles.pdf]");

  c->SetRightMargin(0.15);
  hdcpbf->Draw("colz");
  c->SaveAs("hdcp.pdf");
}

