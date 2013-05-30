using namespace std;
using namespace t2kana;
using namespace t2kana::util;

void test_valor_kdi_splines(
			    bool   compare_with_old_valor = false,
			    double check                  = -1
) {
  TFile * fval         = new TFile("$T2KOSCANA/inputs/official/2013/xsec/XSecSplinesCoarse2013_coarseEreco.root");
  TFile * fkdi_numu    = new TFile("$T2KOSCANA/inputs/official/2013/xsec/spl_numuoa_sk13a_flx11bv32_numu.root");
  TFile * fkdi_numubar = new TFile("$T2KOSCANA/inputs/official/2013/xsec/spl_numuoa_sk13a_flx11bv32_numubar.root");
  TFile * fkdi_nue     = new TFile("$T2KOSCANA/inputs/official/2013/xsec/spl_numuoa_sk13a_flx11bv32_nue.root");
  TFile * fkdi_nuebar  = new TFile("$T2KOSCANA/inputs/official/2013/xsec/spl_numuoa_sk13a_flx11bv32_nuebar.root");
  TFile * fkdi_signue  = new TFile("$T2KOSCANA/inputs/official/2013/xsec/spl_numuoa_sk13a_flx11bv32_signalnue.root");
  TFile * fkdi         = 0;
  if(compare_with_old_valor)
    fkdi = new TFile("$T2KOSCANA/inputs/official/2012/xsec/XSecSplinesCoarse2012_coarseEreco.root");

  TCanvas can;

  stringstream problem_free;

  //sample
  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  // loop over modes for current sample
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    if(!compare_with_old_valor) {
      //find the KDI spline file
      if(util::IsNuMu(m))
	fkdi = fkdi_numu;
      else if(util::IsNuMuBar(m))
	fkdi = fkdi_numubar;
      else if(util::IsNuE(m))
	fkdi = fkdi_nue;
      else if(util::IsNuEBar(m))
	fkdi = fkdi_nuebar;
      else if(util::IsSignalNuE(m))
	fkdi = fkdi_signue;
    }
    vector<string> mode = util::Tokens(util::AsStringCompact_Japan(m), "_");
    //loop over systematics
    for(int ip = 0; ip < conv::kNumOfXSecSyst; ip++) {
      conv::PhysParamId_t p = conv::kListOfXSecSyst[ip];    
      if(p == conv::kPhypUndefined) continue;
      // do we need weight splines for the current systematic parameter?
      bool need_splines = conv::kXSecSystHasWghtSplines[ip];
      if(!need_splines) continue;
      // is the current systematic parameter going to be applied to current sample/mode?
      bool is_applied = conv::kXsApplyXSecSystToSKMode[is][im][ip];
      if(!is_applied) continue;

      //open the output pdf
      TString filename, filenameo, filenamec;
      filenameo.Form("%s_%s.pdf[",util::AsStringCompact(p).c_str(),util::AsStringCompact(m).c_str());
      filename.Form ("%s_%s.pdf", util::AsStringCompact(p).c_str(),util::AsStringCompact(m).c_str());
      filenamec.Form("%s_%s.pdf]",util::AsStringCompact(p).c_str(),util::AsStringCompact(m).c_str());
      can.Print(filenameo.Data(),"pdf");

      //loop over energy bins
      for(int ir = 0; ir < conv::kXsNErecoBins; ir++) {
	double Ereco = conv::kXsErecoBinEdges[ir] + conv::kOneMeV;
	for(int it = 0; it < conv::kXsNEtrueBins; it++) {
	  double Etrue = conv::kXsEtrueBinEdges[is][im][it] + conv::kOneMeV;
	  
	  //now finally can get the splines!
	  ostringstream kdi_splname;
	  kdi_splname << "spl_"
		      << util::AsStringCompact_Japan(p) << "_"
		      << mode[1] << "_"
		      << "enu"  << it << "_" 
		      << "erec" << ir;
	  ostringstream val_splname;
	  val_splname << "spl_"
		      << util::AsStringCompact(s) << "_"
		      << util::AsStringCompact(m) << "_"
		      << util::AsStringCompact(p) << "_"
		      << ir << "_" 
		      << it;

	  // Copy splines from input file
	  TSpline3 * val_spl = (TSpline3 *) fval->Get(val_splname.str().c_str());
	  TSpline3 * kdi_spl = 0;
	  if(compare_with_old_valor)
	    //use val_splname
	    kdi_spl = (TSpline3 *) fkdi->Get(val_splname.str().c_str());
	  else {
	    //use kit_splname
	    TSpline3 * kdi_spl_temp = (TSpline3 *) fkdi->Get(kdi_splname.str().c_str());
	    
	    const int n = kdi_spl_temp->GetNp();
	    double x[n] = {0}, y[n] = 0;
	    for(int k = 0; k < n; k++) {
	      kdi_spl_temp->GetKnot(k, x[k], y[k]);
	      //x[k] += 1;
	      //y[k] += 1;
	    }//k
	    kdi_spl = new TSpline3("",x,y,n);
	  }//compare_with_old_valor

	  if((!val_spl) || (!kdi_spl)) continue;

	  bool problem = false;
	  for(int ibin = -3; ibin <= 3; ibin++) {
	    if(TMath::Abs(val_spl->Eval(ibin) - kdi_spl->Eval(ibin)) > check) {
	      kdi_spl->SetLineColor(kRed);
	      kdi_spl->SetMarkerColor(kRed);
	      kdi_spl->SetMarkerSize(1);
	      kdi_spl->SetMarkerStyle(8);
	      val_spl->SetMarkerSize(1);
	      val_spl->SetMarkerStyle(8);
	      val_spl->Draw("LP");
	      kdi_spl->Draw("LPsame");
	      TPaveText * leg1 = new TPaveText (0,0,1.,0.05,"NDC");
	      leg1->SetBorderSize(0);
	      leg1->SetFillStyle(0);
	      leg1->AddText(Form("%s %s",val_splname.str().c_str(), val_spl->GetTitle()));
	      leg1->Draw();
	      TPaveText * leg2 = new TPaveText (0,0.05,1.,0.1,"NDC");
	      leg2->SetBorderSize(0);
	      leg2->SetFillStyle(0);
	      if(compare_with_old_valor)
		leg2->AddText(Form("%s %s",val_splname.str().c_str(), kdi_spl->GetTitle()));
	      else
		leg2->AddText(Form("%s %s",kdi_splname.str().c_str(), kdi_spl->GetTitle()));
	      leg2->SetTextColor(kRed);
	      leg2->Draw();
	      can.Print(filename.Data(),"pdf");
	      problem = true;
	      break;
	    }
	  }//ibin
	  if(!problem)
	    problem_free << val_splname.str() << "_" << val_splname.str() << endl;
	  delete kdi_spl;
	}//it
      }//ir
      //close the output pdf
      can.Print(filenamec.Data(),"pdf");
    }//ip
  }//im
  cout << problem_free.str() << endl;

  ofstream outfile ("spl_compare.txt");
  long size;
  problem_free.seekg(0,ifstream::end);
  size=problem_free.tellg();
  problem_free.seekg(0);
  outfile.write (problem_free.str().c_str(),size);
  outfile.close();
}
