using namespace std;
using namespace t2kana;
using namespace t2kana::util;

void test_valor_kdi_splines
(
   bool   compare_with_old_valor = false,
   double check                  = -1
) 
{
  // compare response funcs here:
  TFile * fval = new TFile("$T2KOSCANA/inputs/official/2013/xsec/XSecSplinesCoarse2013_coarseEreco.root");

  // with reference response funcs here:
  TFile * fref = 0;

  // locations of reference response funcs (old VaLOR, or KDI)
  TFile * fval_old     = new TFile("$T2KOSCANA/inputs/official/2012/xsec/XSecSplinesCoarse2012_coarseEreco.root");
  TFile * fkdi_numu    = new TFile("$T2KOSCANA/inputs/official/2013/xsec/spl_numuoa_sk13a_flx11bv32_numu.root");
  TFile * fkdi_numubar = new TFile("$T2KOSCANA/inputs/official/2013/xsec/spl_numuoa_sk13a_flx11bv32_numubar.root");
  TFile * fkdi_nue     = new TFile("$T2KOSCANA/inputs/official/2013/xsec/spl_numuoa_sk13a_flx11bv32_nue.root");
  TFile * fkdi_nuebar  = new TFile("$T2KOSCANA/inputs/official/2013/xsec/spl_numuoa_sk13a_flx11bv32_nuebar.root");
  TFile * fkdi_signue  = new TFile("$T2KOSCANA/inputs/official/2013/xsec/spl_numuoa_sk13a_flx11bv32_signalnue.root");

  TCanvas * can = new TCanvas("can");
  can->SetTopMargin(0.13);
  can->SetBottomMargin(0.12);

  stringstream problem_free;

  //print number and names of cross-section systematics considered by this macro
  cout<<"conv::kNumOfXSecSyst = "<<conv::kNumOfXSecSyst<<endl;
  for(int ip = 0; ip < conv::kNumOfXSecSyst; ip++) {
    conv::PhysParamId_t p = conv::kListOfXSecSyst[ip];
    cout<<util::AsString(p) << endl;
  }

  //sample
  conv::SuperkSample_t s = conv::kSKsFCFV1RingMuLike;
  int is = util::ArrIdx::Instance()->SuperkSampleIdx(s);
  // loop over modes for current sample
  for(int im = 0; im < conv::knm; im++) {
    conv::SuperkMode_t m = conv::kSKReactionModes[is][im];
    if(m == conv::kSKmUndefined) continue;
    cout << "Current sample/mode: " << util::AsString(s) << "/" << util::AsString(m) << endl;
    if(!compare_with_old_valor) {
      //find the KDI spline file
      if(util::IsNuMu(m))
	fref = fkdi_numu;
      else if(util::IsNuMuBar(m))
	fref = fkdi_numubar;
      else if(util::IsNuE(m))
	fref = fkdi_nue;
      else if(util::IsNuEBar(m))
	fref = fkdi_nuebar;
      else if(util::IsSignalNuE(m))
	fref = fkdi_signue;
    } 
    else {
      fref = fval_old;
    }
    assert(fref!=0);
    vector<string> mode = util::Tokens(util::AsStringCompact_KDI(m), "_");
    //loop over systematics
    for(int ip = 0; ip < conv::kNumOfXSecSyst; ip++) {
      conv::PhysParamId_t p = conv::kListOfXSecSyst[ip];    
      cout << "** Comparing response functions for cross-section systematic: " << util::AsString(p) << endl;
      if(p == conv::kPhypUndefined) continue;
      // do we need weight splines for the current systematic parameter?
      bool need_splines = conv::kXSecSystHasWghtSplines[ip];
      if(!need_splines) {
         cout << "We don't use response functions for systematic: " << util::AsString(p) << endl;
         continue;
      }
      // is the current systematic parameter going to be applied to current sample/mode?
      bool is_applied = conv::kXsApplyXSecSystToSKMode[is][im][ip];
      if(!is_applied) {
         cout << "Systematic: " << util::AsString(p) 
              << " is not applied to current sample/mode: " 
              << util::AsString(s) << "/" << util::AsString(m) << endl;
         continue;
      }

      TString filename;
      filename.Form ("%s_%s.ps", util::AsStringCompact(p).c_str(),util::AsStringCompact(m).c_str());
      TPostScript * ps = new TPostScript(filename.Data(), 111);

      //loop over energy bins
      for(int ir = 0; ir < conv::kXsNErecoBins1Rmu; ir++) {
	double Ereco = conv::kXsErecoBinEdges1Rmu[ir] + conv::kOneMeV;
	for(int it = 0; it < conv::kXsNEtrueBins; it++) {
	  double Etrue = conv::kXsEtrueBinEdges[is][im][it] + conv::kOneMeV;	  
	  //now finally can get the splines!
	  ostringstream kdi_splname;
	  kdi_splname << "spl_"
		      << util::AsStringCompact_KDI(p) << "_"
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

	  // Get VaLOR spline
          cout << "Getting VaLOR response function: " << val_splname.str() << endl;
          //ip == 0 is MaQE, ip == 1 is MaRes
	  TSpline3 * val_spl = 0;
	  if(ip == 0 || ip == 1) {
            TSpline3 * val_spl_temp = (TSpline3 *) fval->Get(val_splname.str().c_str());
            if(val_spl_temp) {
	      const int nval = val_spl_temp->GetNp();
	      double * xval = new double[nval];
	      double * yval = new double[nval];
	      for(int k = 0; k < nval; k++) {
	        val_spl_temp->GetKnot(k, xval[k], yval[k]);
                //convert horizontal axis from units of x_Ma to units of sigma for VaLOR MaQE and MaRes splines
	        if(ip == 0)
	          xval[k] = (xval[k] * conv::kNomMAQE - conv::kNomMAQE) / conv::kXSecSystFracErr[ip];
                else if(ip == 1)
                  xval[k] = (xval[k] * conv::kNomMARES - conv::kNomMARES) / conv::kXSecSystFracErr[ip];
	      }//k
	      val_spl = new TSpline3("",xval,yval,nval);
	      val_spl->SetTitle(val_spl_temp->GetTitle());
	    }//val_spl_temp
	    else
              cerr << "No " << val_splname.str() << " response func" << endl; 
	  }//ip
          else {
	    val_spl = (TSpline3 *) fval->Get(val_splname.str().c_str());
	    if(!val_spl) {
	      cerr << "No " << val_splname.str() << " response func" << endl;
	    }
	  }
          // Get reference spline
	  TSpline3 * ref_spl = 0;
	  if(compare_with_old_valor) {
	    //use val_splname
            cout << "Getting VaLOR old response function: " << val_splname.str() << endl;
	    ref_spl = (TSpline3 *) fref->Get(val_splname.str().c_str());
	    if(!ref_spl) {
	      cerr << "No " << val_splname.str() << " response func" << endl;
	      }
	  }
	  else {
	    //use kdi_splname
            cout << "Getting KDI response function: " << kdi_splname.str() << endl;
	    TSpline3 * ref_spl_temp = (TSpline3 *) fref->Get(kdi_splname.str().c_str());
            if(ref_spl_temp) {
	      const int n = ref_spl_temp->GetNp();
	      double * x = new double[n];
	      double * y = new double[n];
	      for(int k = 0; k < n; k++) {
		ref_spl_temp->GetKnot(k, x[k], y[k]);
		//the horizontal axes of KDI splines are in units of fractional change in the parameter
		//divide by fractional error in parameter to convert horizontal axes to units of sigma
		x[k] = x[k] / conv::kXSecSystFracErr[ip];
	      }//k
	      ref_spl = new TSpline3("",x,y,n);
	      ref_spl->SetTitle(fref->GetTitle());
//            delete [] x;
//            delete [] y;
	    }//ref_spl_temp
	    else
	      cerr << "No " << kdi_splname.str() << " response func" << endl;
	  }//compare_with_old_valor

	  if((!val_spl) || (!ref_spl)) continue;

	  bool problem = false;
	  for(int ibin = -3; ibin <= 3; ibin++) {
	    if(TMath::Abs(val_spl->Eval(ibin) - ref_spl->Eval(ibin)) > check) {
              // draw
	      ref_spl->SetLineColor(kRed);
	      ref_spl->SetMarkerColor(kRed);
	      ref_spl->SetMarkerSize(1);
	      ref_spl->SetMarkerStyle(8); 
	      val_spl->SetMarkerSize(1);
	      val_spl->SetMarkerStyle(8);
	      val_spl->Draw("LP");
	      ref_spl->Draw("LPsame");
	      TPaveText * leg1 = new TPaveText (0,0.94,1.,1.,"NDC");
	      leg1->SetBorderSize(0);
	      leg1->SetFillStyle(0);
	      leg1->AddText(Form("VaLOR: %s %s",val_splname.str().c_str(), val_spl->GetTitle()));
	      leg1->Draw();
	      TPaveText * leg2 = new TPaveText (0,0.88,1.,0.94,"NDC");
	      leg2->SetBorderSize(0);
	      leg2->SetFillStyle(0);
	      if(compare_with_old_valor)
		leg2->AddText(Form("VaLOR: %s %s",val_splname.str().c_str(), ref_spl->GetTitle()));
	      else
		leg2->AddText(Form("KDI: %s %s",kdi_splname.str().c_str(), ref_spl->GetTitle()));
	      leg2->SetTextColor(kRed);
	      leg2->Draw();
              can->Update(); 
	      problem = true;
	      break;
	    }
	  }//ibin
	  if(!problem)
	    problem_free << val_splname.str() << "_" << val_splname.str() << endl;
	  delete val_spl;
	  delete ref_spl;
	}//it
      }//ir

      can->Clear();
      can->Update();
      ps->Close();
      delete ps;

      // convert to pdf and erase original ps
      gSystem->Exec(Form("ps2pdf14 %s",filename.Data()));
      gSystem->Exec(Form("rm %s",filename.Data()));

    }//ip
  }//im
  cout << problem_free.str() << endl;

  delete can;

  ofstream outfile ("spl_compare.txt");
  long size;
  problem_free.seekg(0,ifstream::end);
  size=problem_free.tellg();
  problem_free.seekg(0);
  outfile.write (problem_free.str().c_str(),size);
  outfile.close();
}
