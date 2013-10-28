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
			 //const char * filename_tree = "run1-4_nuis4_sintheta_normal/combination/final.root",
			 const char * filename_tree = "run1-4_nuis4_sintheta_inverted31/combination/final.root",
			 const char * filename_hist = "final.root",
			 const char * filename_out  = "dcp_fit_correlation_temp.root",
			 int verbose                = 0,
			 bool scan_dcp              = false,
			 bool scan_mix13            = false,
			 bool scan_mix12            = false,
			 bool scan_dm21sq           = false
)
{
  using namespace t2kana;
  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsqTheta;
  //Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;
  //Opt::Instance()->gOscDmsqScheme        = conv::kNuDmsq12andAbs23;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhInverted;
  Opt::Instance()->gOscDmsqScheme        = conv::kNuDmsq12andAbs31;

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
  if(scan_dcp || scan_mix13 || scan_mix12 || scan_dm21sq)
    skexp->Read(config.Data());

  const char * xmlfile = "data.xml";
  const int iexperiment = 8;
  bool using_xml_file = true;
  TString xmlfile2 = Form("%s/inputs/%s", gSystem->Getenv("T2KOSCANA"), xmlfile);
  SuperkObservation * skobs = new SuperkObservation;
  if(scan_dcp || scan_mix13 || scan_mix12 || scan_dm21sq)
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
  TH2D * hmix13bf = new TH2D("hmix13bf","",nbinsmix,mixmin,mixmax,nbinsdm,dmmin,dmmax);
  TH2D * hmix12bf = new TH2D("hmix13bf","",nbinsmix,mixmin,mixmax,nbinsdm,dmmin,dmmax);
  TH2D * hdm21sqbf = new TH2D("hmix13bf","",nbinsmix,mixmin,mixmax,nbinsdm,dmmin,dmmax);

  if(scan_dcp)
    c->SaveAs("dcp_profiles.pdf[");
  if(scan_mix13)
    c->SaveAs("mix13_profiles.pdf[");
  if(scan_mix12)
    c->SaveAs("mix12_profiles.pdf[");
  if(scan_dm21sq)
    c->SaveAs("dm21sq_profiles.pdf[");

  TFile fout;
  if(scan_dcp || scan_mix13 || scan_mix12 || scan_dm21sq)
    fout.Open(filename_out,"RECREATE");
  
  cout << "Starting loop over " << fitnt->GetEntries() << " entries" << endl;
  for(int iev=0; iev<fitnt->GetEntries(); iev++) {
    fitnt->GetEntry(iev);
    //if(fparm->GetFitTries() == 1) continue;
    //if(fparm->GetFitTries() != 99) continue;
 
    if(!hdcp[idm32sq][imix23] && scan_dcp)
      hdcp[idm32sq][imix23] = new TH1D("",
				       Form("mix %f, dm %f;#delta_{CP}",fparm->BestFitValue(conv::kPhypMixing23),fparm->BestFitValue(conv::kPhypDm23sq)),
				       npibins,-conv::kpi,+conv::kpi);
    if(!hmix13[idm32sq][imix23] && scan_mix13)
      hmix13[idm32sq][imix23] = new TH1D("",
					 Form("mix %f, dm %f;sin^{2}#theta_{13}",fparm->BestFitValue(conv::kPhypMixing23),fparm->BestFitValue(conv::kPhypDm23sq)),
					 npibins,kGlobFitSinsqTheta13*(1-(kGlobFitErrSinsqTheta13*7)),kGlobFitSinsqTheta13*(1+(kGlobFitErrSinsqTheta13*7))); 
    if(!hmix12[idm32sq][imix23] && scan_mix12)
      hmix12[idm32sq][imix23] = new TH1D("",
				       Form("mix %f, dm %f;sin^{2}#theta_{12}",fparm->BestFitValue(conv::kPhypMixing23),fparm->BestFitValue(conv::kPhypDm23sq)),
				       npibins,kGlobFitSinsqTheta12*(1-(kGlobFitErrSinsqTheta12*7)),kGlobFitSinsqTheta12*(1+(kGlobFitErrSinsqTheta12*7))); 
    if(!hdm21sq[idm32sq][imix23] && scan_dm21sq)
      hdm21sq[idm32sq][imix23] = new TH1D("",
				       Form("mix %f, dm %f;#Deltam^{2}_{21}",fparm->BestFitValue(conv::kPhypMixing23),fparm->BestFitValue(conv::kPhypDm23sq)),
				       npibins,kGlobFitDm23sq*(1-(kGlobFitErrDm23sq*7)),kGlobFitDm23sq*(1+(kGlobFitErrDm23sq*7))); 
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

    cout << fparm->GetStatusCode("migrad_min") << "\t" << fparm->GetStatusCode("migrad_err") << "\t" << fparm->GetFitTries() << "\t" << fparm->BestFitValue(conv::kPhypDm23sq) << "\t" << fparm->BestFitValue(conv::kPhypMixing23) << "\t" << fparm->BestFitValue(conv::kPhypDeltaCP) << endl;
    //if(TMath::Abs(fparm->BestFitValue(conv::kPhypDeltaCP) < 1E-4) || TMath::Abs(fparm->BestFitValue(conv::kPhypDeltaCP) - conv::kpi) < 1E-4 || TMath::Abs(fparm->BestFitValue(conv::kPhypDeltaCP) + conv::kpi) < 1E-4)
    //if(TMath::Abs(fparm->BestFitValue(conv::kPhypDeltaCP)) > conv::kpi)
    //hdcpbf->SetBinContent(imix23,idm32sq,1);
    hdcpbf->SetBinContent(imix23,idm32sq,fparm->BestFitValue(conv::kPhypDeltaCP));
    hmix13bf->SetBinContent(imix23,idm32sq,fparm->BestFitValue(conv::kPhypMixin13));
    hmix12bf->SetBinContent(imix23,idm32sq,fparm->BestFitValue(conv::kPhypMixing12));
    hdm21sqbf->SetBinContent(imix23,idm32sq,fparm->BestFitValue(conv::kPhypDm12sq));

    if(scan_dcp || scan_mix13 || scan_mix12 || scan_dm21sq) {
      //if(fparm->BestFitValue(conv::kPhypDm23sq) < 2.41E-3 || fparm->BestFitValue(conv::kPhypDm23sq) > 2.71E-3)
      //if(fparm->BestFitValue(conv::kPhypDm23sq) < 2.99E-3)
      //continue;
      //if(fparm->BestFitValue(conv::kPhypMixing23) < 0.399 || fparm->BestFitValue(conv::kPhypMixing23) > 0.451)
      //continue;
      fparm->CopyBestFitToCurrValues();

      FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);

      //cout << PhysHypothesis(*fparm) << endl;
      cout << "\t<<<<<  " << fparm->GetBestFitChisq() << endl;


      if(scan_dcp) {
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
      }//scan_dcp


      ///need to copy above & remember to reset!




      delete fit_inputs;
    }//scan_dcp
  }//iev

  if(scan_dcp)
    c->SaveAs("dcp_profiles.pdf]");
  if(scan_mix13)
    c->SaveAs("mix13_profiles.pdf]");
  if(scan_mix12)
    c->SaveAs("mix12_profiles.pdf]");
  if(scan_dm21sq)
    c->SaveAs("dm21sq_profiles.pdf]");

  c->SetRightMargin(0.15);
  hdcpbf->Draw("colz");
  c->SaveAs("hdcp.pdf");
  hmix13bf->Draw("colz");
  c->SaveAs("hmix13.pdf");
  hmix12bf->Draw("colz");
  c->SaveAs("hmix12.pdf");
  hdm21sqbf->Draw("colz");
  c->SaveAs("hdm21sq.pdf");
}

