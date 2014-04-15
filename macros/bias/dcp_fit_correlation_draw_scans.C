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
				    //const char * filename_tree = "run1-4_nuis4_sintheta_normal/combination/final.root",
				    const char * filename_tree = "run1-4_nuis4_sintheta_inverted31/combination/final.root",
				    const char * filename_hist = "dcp_fit_correlation_inverted31.root",
				    const char * filename_out  = "dcp_fit_correlation_inverted31_scans.root",
				    int verbose                = 0,
				    bool scan_dcp = true,
				    bool scan_energies = true
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

  TFile f(filename_hist,"READ");
  TFile fout;
  if(scan_energies)
    fout.Open(filename_out,"RECREATE");

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
  const int nmodes = 6; //numu->numu etc
  const int nenergies = conv::kNEtrueBins1Rmu;
  TCanvas * ctemp = new TCanvas();
  ctemp->SaveAs("dcp_profiles_draw.pdf[");
  ctemp->SaveAs("dcp_profiles_probabilities_unphysical.pdf[");
 
  NuOsc3Flavour  fNuOsc3Flavour;

  cout << "Starting loop over " << fitnt->GetEntries() << " entries" << endl;
  int i_energy_scans = 0;
  for(int iev=0; iev<fitnt->GetEntries(); iev++) {
    fitnt->GetEntry(iev);
    hdcp[idm32sq][imix23] = 0;
    f.GetObject(Form("mix_%d_dm_%d",imix23,idm32sq), hdcp[idm32sq][imix23]);
    if(!hdcp[idm32sq][imix23]) continue;

    if(TMath::Abs(hdcp[idm32sq][imix23]->GetMinimum() - fparm->GetBestFitChisq()) < 1E-3) continue;

    //if(fparm->GetFitTries() != 99) continue;


    //if(fparm->BestFitValue(conv::kPhypDm23sq) < 2.29E-3) continue;
    //if(fparm->BestFitValue(conv::kPhypDm23sq) > 2.61E-3) continue;
    //if(fparm->BestFitValue(conv::kPhypMixing23) > 0.503) { continue; }
    //if(fparm->BestFitValue(conv::kPhypMixing23) < 0.497) { continue; }
 
    cout << endl << endl << endl << *fparm << endl;

    double chisq = 0;
    if(scan_dcp) {
      fparm->CopyBestFitToCurrValues();
      FitInputs * fit_inputs = new FitInputs(skobs,skexp,fparm);
      chisq = fit_inputs->Chisquare();//("poisson",conv::kSKsFCFV1RingMuLike);
    }//scan_dcp

    if(scan_energies) {
      double mixing12 = fparm->BestFitValue(conv::kPhypMixing12);
      double mixing13 = fparm->BestFitValue(conv::kPhypMixing13);
      double mixing23 = fparm->BestFitValue(conv::kPhypMixing23);
      double dmsq     = fparm->BestFitValue(conv::kPhypDm12sq);
      double Dmsq     = fparm->BestFitValue(conv::kPhypDm23sq);

      TCanvas * cprob[nmodes][nenergies] = {0};
      TH1D * hprob[nmodes][nenergies];
      for(int im = 0; im < nmodes; im++) {
	string modename;
	switch(im) {
	case 0: modename = "numu_to_numu"; break;
	case 1: modename = "numubar_to_numubar"; break;
	case 2: modename = "nue_to_nue"; break;
	case 3: modename = "nuebar_to_nuebar"; break;
	case 4: modename = "numu_to_nue"; break;
	case 5: modename = "nue_to_numu"; break;
	default: modename = ""; break;
	}//switch im
	ctemp->SaveAs(Form("dcp_profiles_probabilities_%s_%d.pdf[",modename.c_str(),i_energy_scans));
	for(int ie = 0; ie < nenergies; ie++) {
	  //hprob[im][ie] = new TH1D("",Form(";;%s",modename.c_str()),100,fparm->BestFitValue(conv::kPhypCPV    ),fparm->BestFitValue(conv::kPhypCPV    )+conv::kpi);
	  hprob[im][ie] = new TH1D("",Form(";;%s",modename.c_str()),200,-conv::kpi,conv::kpi);
	  double Etrue = (util::EtrueBinEdges(conv::kSKsFCFV1RingMuLike,ie-1) + util::EtrueBinEdges(conv::kSKsFCFV1RingMuLike,ie) ) / 2;
	  for(int ip = 1; ip <= hprob[im][ie]->GetNbinsX(); ip++) {
	    double delta_cp = hprob[im][ie]->GetBinLowEdge(ip);      
	    fNuOsc3Flavour.SetOscParametersSinSqTAll(mixing12,mixing13,mixing23,dmsq,Dmsq,delta_cp);
	    double prob = -1;
	    switch(im) {
	    case 0:
	      prob = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuMu,    conv::kPdgNuMu,    Etrue);
	      break;
	    case 1:
	      prob = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuMuBar, conv::kPdgNuMuBar, Etrue);
	      break;
	    case 2:
	      prob = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuE,     conv::kPdgNuE,     Etrue);
	      break;
	    case 3:
	      prob = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuEBar,  conv::kPdgNuEBar,  Etrue);
	      break;
	    case 4:
	      prob = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuMu,    conv::kPdgNuE,     Etrue);
	      break;
	    case 5:
	      prob = fNuOsc3Flavour.CalcMatterOscProb(conv::kPdgNuE,     conv::kPdgNuMu,    Etrue);
	      break;
	    default:
	      prob = -1;
	      break;
	    }
	    hprob[im][ie]->SetBinContent(ip,prob);
	  }//ip
	  //hprob[im][ie]->SetDirectory(0);
	  //fout.cd();
	  hprob[im][ie]->Write(Form("dcp_profiles_probabilities_%s_%d_%f",modename.c_str(),i_energy_scans,Etrue));
	  double minimum = hprob[im][ie]->GetMinimum();
	  double maximum = hprob[im][ie]->GetMaximum();
	  double extrema = TMath::Abs(maximum - minimum);
	  //extrema = 0.001;//5;
	  //apply a cut - ones with 
	  bool draw = false;
	  for(int ip = 1; ip < hprob[im][ie]->GetNbinsX(); ip++) {
	    double now  = hprob[im][ie]->GetBinContent(ip);
	    double next = hprob[im][ie]->GetBinContent(ip+1);
	    if(TMath::Abs(now - next) > extrema*0.5) {
	      draw = true;
	      break;
	    }
	  }//ip
	  //draw = true;
	  //if unphysical, draw it.
	  if(minimum < 0 || maximum > 1 || draw) {
	    ctemp->cd();
	    hprob[im][ie]->Draw();
	    TLegend * legg = new TLegend(0,0,1,0.1);
	    legg->SetTextSize(0.025);
	    legg->AddEntry((TObject*)0,Form("E:%f,#Delta_{32}=%f,#theta_{23}=%f,#theta_{12}=%f,#theta_{13}=%f,#Delta_{21}=%f,event:",Etrue,Dmsq,mixing23,mixing12,mixing13,dmsq,iev),"l");
	    legg->Draw();
	    ctemp->SaveAs("dcp_profiles_probabilities_unphysical.pdf");
	    delete legg;
	  }
	  if(draw) {
	    //draw everything that passes the cut
	    cprob[im][ie] = new TCanvas();
	    hprob[im][ie]->Draw();
	    TLegend * legg = new TLegend(0,0,1,0.1);
	    legg->SetTextSize(0.03);
	    legg->AddEntry((TObject*)0,Form("E:%f,#Delta_{32}=%f,#theta_{23}=%f,#theta_{12}=%f,#theta_{13}=%f,#Delta_{21}=%f",Etrue,Dmsq,mixing23,mixing12,mixing13,dmsq),"l");
	    legg->Draw();
	    cprob[im][ie]->SaveAs(Form("dcp_profiles_probabilities_%s_%d.pdf",modename.c_str(),i_energy_scans));
	    delete legg;
	  }
	  delete cprob[im][ie];
	  delete hprob[im][ie];
	}//ie
	ctemp->SaveAs(Form("dcp_profiles_probabilities_%s_%d.pdf]",modename.c_str(),i_energy_scans));
      }//im
      i_energy_scans++;
    }//scan energies?

    c[idm32sq][imix23] = new TCanvas();
    hdcp[idm32sq][imix23]->Draw();
    TLegend * l = new TLegend(0,0,0.3,0.1);
    l->AddEntry((TObject*)0, hdcp[idm32sq][imix23]->GetTitle(), "l");
    l->Draw();
    TLegend * leg = new TLegend(0.3,0,1.0,0.1);
    leg->AddEntry((TObject*)0, Form("#chi^{2}_{bf}=%f, tries=%d, #delta_{CP}=%f [%d,%d], #chi^{2}=%f", fparm->GetBestFitChisq(), fparm->GetFitTries(), fparm->BestFitValue(conv::kPhypCPV    ), imix23, idm32sq, chisq), "l");
    leg->SetTextSize(0.03);
    leg->Draw();
    TLine *l1 = new TLine(-conv::kpi,fparm->GetBestFitChisq(),+conv::kpi,fparm->GetBestFitChisq());
    l1->SetLineColor(kRed);
    l1->Draw();
    TLine *l2 = new TLine(fparm->BestFitValue(conv::kPhypCPV    ),0,fparm->BestFitValue(conv::kPhypCPV    ),200);
    l2->SetLineColor(kRed);
    l2->Draw();
  }//iev

  for(int idm32sq = 0; idm32sq < nbinsdm; idm32sq++)
  //idm32sq = 9;//this is the transition bin 2.5E-3 in the chisq surface i was studying
    for(int imix23 = 0; imix23 < nbinsmix; imix23++)
      if(c[idm32sq][imix23])
	c[idm32sq][imix23]->SaveAs("dcp_profiles_draw.pdf");


  ctemp->SaveAs("dcp_profiles_draw.pdf]");
  ctemp->SaveAs("dcp_profiles_probabilities_unphysical.pdf]");
}

