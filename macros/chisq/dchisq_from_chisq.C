//
// Convert a chisq surface into a dchisq surface
// chisq_min comes either from argument, or the minimum of the surface
//
void dchisq_from_chisq
(
 const char * chisq_surf_filename,
 const char * chisq_surf_histname,
 double       chisq_bestfit,           //if negative, subtract the minimum
 const char * outfilename,
 const char * outhistname = "hdchisq",
 double       chisq_scale = 1          //multiply the CHISQ surface by this
)
{
  //read in histogram
  TFile f_in(chisq_surf_filename);
  TH2D * h = 0;
  f_in.GetObject(chisq_surf_histname,h);
  if(!h) { cerr << "Can't find histogram " << chisq_surf_histname << endl; exit(1); }
  
  //get minimum
  if(chisq_bestfit < 0) {
    for(int i = 1; i <= h->GetNbinsX(); i++) {
      chisq_bestfit = 9999 * chisq_scale;
      for(int j = 1; j <= h->GetNbinsY(); j++) {
	if(TMath::Abs(h->GetBinContent(i,j)) < 1E-5) {
	  h->SetBinContent(i,j, 10000 * chisq_scale);
	  cout << h->GetBinContent(i,j)  << "\t" << i << "\t" << j<< endl;
	}
	//if(h->GetBinContent(i,j) < chisq_bestfit)
	//chisq_bestfit = h->GetBinContent(i,j);
      }//j
    }//i
    chisq_bestfit = h->GetBinContent(h->GetMinimumBin());
  }

  //turn into DCHISQ
  for(int i = 1; i <= h->GetNbinsX(); i++) {
    for(int j = 1; j <= h->GetNbinsY(); j++) {
      h->SetBinContent(i,j, (h->GetBinContent(i,j) * chisq_scale) - chisq_bestfit);
    }
  }

  //print out the min (to make sure)
  double min_dchisq = 1;
  for(int i = 1; i <= h->GetNbinsX(); i++) {
    for(int j = 1; j <= h->GetNbinsY(); j++) {
      if(h->GetBinContent(i,j) < min_dchisq) {
	min_dchisq = h->GetBinContent(i,j);
	cout << h->GetXaxis()->GetBinCenter(i) << "\t" << h->GetYaxis()->GetBinCenter(j) << "\t" << h->GetBinContent(i,j) << "\t" << i << "\t" << j << endl;
      }
    }
  }

  //save
  TFile f_out(outfilename,"RECREATE");
  h->Write(outhistname);
  h->SetContour(2);
  h->SetContourLevel(0,2.30);
  h->SetContourLevel(1,4.61);
  h->Write(Form("%s_cont3",outhistname));
  f_out.Close();
}
