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
			 const int    npibins       = 50,
			 //const char * filename_tree = "run1-4_nuis4_sintheta_normal/combination/final.root",
			 const char * filename_tree = "run1-4_nuis4_sintheta_inverted31/combination/final.root",
			 const char * filename_out  = "dcp_fit_correlation_temp.root",
			 int verbose                = 0,
			 bool scan_dcp              = true,
			 bool scan_mix13            = true,
			 bool scan_mix12            = true,
			 bool scan_dm21sq           = true
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
  TH1D * hdcp   [nbins1][nbins2] = {0};
  TH1D * hmix13 [nbins1][nbins2] = {0};
  TH1D * hmix12 [nbins1][nbins2] = {0};
  TH1D * hdm21sq[nbins1][nbins2] = {0};
  TH2D * hdcpbf    = new TH2D("hdcpbf",   "",nbinsmix,mixmin,mixmax,nbinsdm,dmmin,dmmax);
  TH2D * hmix13bf  = new TH2D("hmix13bf", "",nbinsmix,mixmin,mixmax,nbinsdm,dmmin,dmmax);
  TH2D * hmix12bf  = new TH2D("hmix12bf", "",nbinsmix,mixmin,mixmax,nbinsdm,dmmin,dmmax);
  TH2D * hdm21sqbf = new TH2D("hdm21sqbf","",nbinsmix,mixmin,mixmax,nbinsdm,dmmin,dmmax);

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
					 npibins,conv::kGlobFitSinsqTheta13*(1-(conv::kGlobFitErrSinsqTheta13*7)),conv::kGlobFitSinsqTheta13*(1+(conv::kGlobFitErrSinsqTheta13*7))); 
    if(!hmix12[idm32sq][imix23] && scan_mix12)
      hmix12[idm32sq][imix23] = new TH1D("",
				       Form("mix %f, dm %f;sin^{2}#theta_{12}",fparm->BestFitValue(conv::kPhypMixing23),fparm->BestFitValue(conv::kPhypDm23sq)),
				       npibins,conv::kGlobFitSinsqTheta12*(1-(conv::kGlobFitErrSinsqTheta12*7)),conv::kGlobFitSinsqTheta12*(1+(conv::kGlobFitErrSinsqTheta12*7))); 
    if(!hdm21sq[idm32sq][imix23] && scan_dm21sq)
      hdm21sq[idm32sq][imix23] = new TH1D("",
				       Form("mix %f, dm %f;#Deltam^{2}_{21}",fparm->BestFitValue(conv::kPhypMixing23),fparm->BestFitValue(conv::kPhypDm23sq)),
				       npibins,conv::kGlobFitDm12sq*(1-(conv::kGlobFitErrDm12sq*7)),conv::kGlobFitDm23sq*(1+(conv::kGlobFitErrDm12sq*7))); 
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
    hmix13bf->SetBinContent(imix23,idm32sq,fparm->BestFitValue(conv::kPhypMixing13));
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
	  double dcp = hdcp[idm32sq][imix23]->GetBinCenter(ix);
	  fparm->SetInitValue(conv::kPhypDeltaCP, dcp);
	  skexp->Generate(*fparm);      
	  double chisq = fit_inputs->Chisquare();//("poisson",conv::kSKsFCFV1RingMuLike);
	  hdcp[idm32sq][imix23]->Fill(dcp, chisq);
	  cout << chisq << "\t";
	}
	cout << endl;
	//save & draw
	hmix13[idm32sq][imix23]->GetYaxis()->SetTitle("#chi^{2}");
#line 0
	hmix13[idm32sq][imix23]->GetXaxis()->SetTitle(util::AsString(conv::kPhypDeltaCP).c_str());
	hdcp[idm32sq][imix23]->Write(Form("dcp_mix_%d_dm_%d",imix23,idm32sq));
	hdcp[idm32sq][imix23]->Draw();
	TLegend * l = new TLegend(0,0,0.3,0.1);
	l->AddEntry(hdcp[idm32sq][imix23], hdcp[idm32sq][imix23]->GetTitle(), "l");
	l->Draw();
	c->SaveAs("dcp_profiles.pdf");
      }//scan_dcp

      fparm->CopyBestFitToCurrValues();
      if(scan_mix13) {
	for(int ix = 1; ix <= hmix13[idm32sq][imix23]->GetNbinsX(); ix++) {
	  double mix13 = hmix13[idm32sq][imix23]->GetBinCenter(ix);
	  fparm->SetInitValue(conv::kPhypMixing13, mix13);
	  skexp->Generate(*fparm);      
	  double chisq = fit_inputs->Chisquare();//("poisson",conv::kSKsFCFV1RingMuLike);
	  hmix13[idm32sq][imix23]->Fill(mix13, chisq);
	  cout << chisq << "\t";
	}
	cout << endl;
	//save & draw
	hmix13[idm32sq][imix23]->GetYaxis()->SetTitle("#chi^{2}");
	hmix13[idm32sq][imix23]->GetXaxis()->SetTitle(util::AsString(conv::kPhypMixing13).c_str());
	hmix13[idm32sq][imix23]->Write(Form("mix13_mix_%d_dm_%d",imix23,idm32sq));
	hmix13[idm32sq][imix23]->Draw();
	TLegend * l = new TLegend(0,0,0.3,0.1);
	l->AddEntry(hmix13[idm32sq][imix23], hmix13[idm32sq][imix23]->GetTitle(), "l");
	l->Draw();
	c->SaveAs("mix13_profiles.pdf");
      }//scan_mix13

      fparm->CopyBestFitToCurrValues();
      if(scan_mix12) {
	for(int ix = 1; ix <= hmix12[idm32sq][imix23]->GetNbinsX(); ix++) {
	  double mix12 = hmix12[idm32sq][imix23]->GetBinCenter(ix);
	  fparm->SetInitValue(conv::kPhypMixing12, mix12);
	  skexp->Generate(*fparm);      
	  double chisq = fit_inputs->Chisquare();//("poisson",conv::kSKsFCFV1RingMuLike);
	  hmix12[idm32sq][imix23]->Fill(mix12, chisq);
	  cout << chisq << "\t";
	}
	cout << endl;
	//save & draw
	hmix12[idm32sq][imix23]->GetYaxis()->SetTitle("#chi^{2}");
	hmix12[idm32sq][imix23]->GetXaxis()->SetTitle(util::AsString(conv::kPhypMixing12).c_str());
	hmix12[idm32sq][imix23]->Write(Form("mix12_mix_%d_dm_%d",imix23,idm32sq));
	hmix12[idm32sq][imix23]->Draw();
	TLegend * l = new TLegend(0,0,0.3,0.1);
	l->AddEntry(hmix12[idm32sq][imix23], hmix12[idm32sq][imix23]->GetTitle(), "l");
	l->Draw();
	c->SaveAs("mix12_profiles.pdf");
      }//scan_mix12

      fparm->CopyBestFitToCurrValues();
      if(scan_dm21sq) {
	for(int ix = 1; ix <= hdm21sq[idm32sq][imix23]->GetNbinsX(); ix++) {
	  double dm21sq = hdm21sq[idm32sq][imix23]->GetBinCenter(ix);
	  fparm->SetInitValue(conv::kPhypDm12sq, dm21sq);
	  skexp->Generate(*fparm);      
	  double chisq = fit_inputs->Chisquare();//("poisson",conv::kSKsFCFV1RingMuLike);
	  hdm21sq[idm32sq][imix23]->Fill(dm21sq, chisq);
	  cout << chisq << "\t";
	}
	cout << endl;
	//save & draw
	hdm21sq[idm32sq][imix23]->GetYaxis()->SetTitle("#chi^{2}");
	hdm21sq[idm32sq][imix23]->GetXaxis()->SetTitle(util::AsString(conv::kPhypDm12sq).c_str());
	hdm21sq[idm32sq][imix23]->Write(Form("dm21sq_mix_%d_dm_%d",imix23,idm32sq));
	hdm21sq[idm32sq][imix23]->Draw();
	TLegend * l = new TLegend(0,0,0.3,0.1);
	l->AddEntry(hdm21sq[idm32sq][imix23], hdm21sq[idm32sq][imix23]->GetTitle(), "l");
	l->Draw();
	c->SaveAs("dm21sq_profiles.pdf");
      }//scan_dm21sq

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

