//
// ------------------------------------------------------------------------------------------------
// Plot predicted SuperK reconstructed neutrino energy spectrum for each of the cross section
// systematics that use splines (those both in the BANFF fit, and outside of it).
// Do this with unoscillated, untweaked spectra, so as to check that our implementation of
// how to distribute weights from big bins which the splines are created in, to our small
// 50MeV bins is correct.
//
// Tom Dealtry
// ------------------------------------------------------------------------------------------------
//

void test_spline_implementation
(
  bool         draw               = false,                          /* draw & save the canvas? */
  const char * config_file        = "analysis_configuration.xml"    /* */
)
{
  using namespace std;
  using namespace t2kana;

  // Analysis configuration
  const char * config  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config_file);

  // Global options
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  // Set physics hypothesis - no oscillations!
  PhysHypothesis ph("2f_numu_disap");
  ph.SetCurrValue(conv::kPhypMixing23,  0);
  ph.SetCurrValue(conv::kPhypDm23sq,    0);
  
  // Print options
  cout << "\n[Job inputs]" << endl;
  cout << "- Reading analysis configuration from: " << config << endl;
  cout << *Opt::Instance() << endl;

  // Create & config superk spectrum predictor
  SuperkExpectation skexp;
  skexp.Read(config);

  // Get the BANFF covariance matrix (for "1-sigma" BANFF tweaks)
  TMatrixD cov = skexp.GetBANFF()->CovMatrix();

  // Tweak range & histogram drawing options
  const int ntwkdials = 5;
  double twkdial          [ntwkdials] = { -2.0,    -1.0,    0.0,       +1.0,    +2.0    };
  const char * twkdialstr [ntwkdials] = { "m2",    "minus", "nominal", "plus",  "p2"    };
  int hwidth              [ntwkdials] = { 1,       1,       2,         1,       1       };
  int hcolor              [ntwkdials] = { kBlue,   kBlue,   kBlack,    kRed,    kRed    };
  int hstyle              [ntwkdials] = { kDotted, kDashed, kSolid,    kDashed, kDotted };
 
  // Histograms & canvases (a few too many than I actually need) & output file
  TH1D * hEreco   [conv::kNumOfXSecSyst][conv::knm][ntwkdials];
  TCanvas * c[conv::kNumOfXSecSyst][conv::knm];
  TLegend * legend = 0;
  TFile f("test_spline_implementation.root","recreate");

  // Integer representation of 1Rmu
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(conv::kSKsFCFV1RingMuLike);

  // Get the tweaks
  for(int ip = 0; ip < conv::kNumOfXSecSyst; ip++) {
    //get the parameter to tweak
    conv::PhysParamId_t p = conv::kListOfXSecSyst[ip];
    if(p == conv::kPhypUndefined) continue;
    if(!conv::kXSecSystHasWghtSplines[ip]) continue;
    string pstrc = util::AsStringCompact(p);
    string pstr  = util::AsString(p);
    cout << " *** " << pstr << endl;
    //get the banff parameter number (if it is a banff param)
    int bindex = -1;
    if(conv::kXSecSystIncludedInBANFF[ip]) {
      switch(p) {
        case ( conv::kNuisBANFF22 ) : bindex = 22; break;
        case ( conv::kNuisBANFF23 ) : bindex = 23; break;
        default : cout << "unknown banff parameter " << pstr << endl; return;
      }
    }//banff parameter?

    for (int imode = 0; imode < conv::knm; imode++) {
      //get the mode to look at
      conv::SuperkMode_t m = conv::kSKReactionModes[is][imode];
      if( m == conv::kSKmUndefined ) continue;
      string mstr = util::AsString(m);
      string mstrc = util::AsStringCompact(m);
      cout << " *** " << mstr << endl;
      if(draw) {
	c[ip][imode] = new TCanvas(Form("c%d_%d",ip,imode),"",20,20,500,500);
	c[ip][imode]->cd(1);
	legend = new TLegend(0.3, 0.7, 0.85, 0.9);
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	legend->SetHeader(pstr.c_str());
      }//draw?

      //get the MC template (etrue-ereco)
      TH2D * mctempl = skexp.GetSuperkMC()->Template(conv::kSKsFCFV1RingMuLike,m,conv::kSKMCtNom);

      //loop over dials
      for(int idial = 0; idial < ntwkdials; idial++) {
	//get the spectrum histogram (as a template for the binning)
	hEreco[ip][imode][idial] = new TH1D(*skexp.CurrPredictedErecoSpectrum(conv::kSKsFCFV1RingMuLike,m));
	hEreco[ip][imode][idial]->SetDirectory(0);
	// loop over Ereco
	for(int ir = 1; ir <= hEreco[ip][imode][idial]->GetNbinsX(); ir++) {
	  double nev = 0;
	  double Ereco = mctempl->GetXaxis()->GetBinCenter(ir);
          // loop over Etrue
	  for(int it = 1; it <= mctempl->GetYaxis()->GetNbins(); it++) {
	    double Etrue = mctempl->GetYaxis()->GetBinCenter(it);
	    //apply the systematic?
	    //if(!skexp.GetNuXSec()->ApplySyst(p,s,m,Ereco,Etrue)) continue; //private method
	    bool apply = conv::kXsApplyXSecSystToSKMode[is][imode][ip] &&
	      (Ereco >= conv::kXSecSystErecoCut[ip][0] && Ereco < conv::kXSecSystErecoCut[ip][1]) &&
	      (Etrue >= conv::kXSecSystEtrueCut[ip][0] && Etrue < conv::kXSecSystEtrueCut[ip][1]);
	    if(!apply) continue;
	    //get the weight
	    TSpline3 * wght_spl = skexp.GetNuXSec()->FindWghtSpline(p,conv::kSKsFCFV1RingMuLike,m,Ereco,Etrue);
	    assert(wght_spl);
	    double relerr = 0;
	    if(conv::kXSecSystIncludedInBANFF[ip]) {
	      util::PrintCompact(cov);
	      cout << bindex << endl;
	      skexp.GetBANFF()->CovMatrix().Print();
	      cout << skexp.GetBANFF()->CovMatrix()(bindex,bindex) << endl;
	      relerr = skexp.GetBANFF()->CovMatrix()(bindex,bindex);
	    }
	    else
	      relerr = skexp.GetNuXSec()->RelativeError(p);
	    double wght = wght_spl->Eval(1 + (twkdial[idial] * relerr)) / wght_spl->Eval(1);
	    wght = TMath::Max(0.,wght);
	    //add the weighted number of events in the et-er bin to the nev for the er bin
	    nev += mctempl->GetBinContent(ir,it) * wght;
          }//it
          hEreco[ip][imode][idial]->SetBinContent(ir,nev);
	}//ir
	const char * label = Form("tweak dial = %.2f", twkdial[idial]);
	if(draw) {
	  hEreco[ip][imode][idial]->SetLineWidth(hwidth[idial]);
	  hEreco[ip][imode][idial]->SetLineColor(hcolor[idial]);
	  hEreco[ip][imode][idial]->SetLineStyle(hstyle[idial]);
	  const char * opt   = (idial==0) ? "" : "same";
	  if(idial==0) {
	    hEreco[ip][imode][idial]->GetXaxis()->SetTitle("E_{reco} (GeV)");
	    hEreco[ip][imode][idial]->GetYaxis()->SetTitle("Number of events");
	  }
	  hEreco[ip][imode][idial]->Draw(opt);
	  legend->AddEntry(hEreco[ip][imode][idial], label, "L");
	}//draw?
	cout << " * " << label << " : NumEvents = " << hEreco[ip][imode][idial]->Integral() << endl;
	// Write-out in file
	const char * title = Form("tweakEreco_1Rmu_%s_%s_%s", mstrc.c_str(), pstrc.c_str(), twkdialstr[idial]);
	cout << title << endl;
	hEreco[ip][imode][idial]->Write(title);
      }//idial
      //last bit of drawing
      if(draw) {
	legend->Draw();
	c[ip][imode]->Update();
	delete legend;
      }//draw?
    }//imode
  }//ip

  f.Close();
}

