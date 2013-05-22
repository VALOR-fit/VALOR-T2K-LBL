//-----------------------------------------------------------------------------------------------------
//Use the outputs of "2d_nuisance_fit.C" and make ChiSqSurface histograms
// ..................................................................
//
// Jelena Ilic
// ---------------------------------------------------------------------------------------    

void makeChiSqSurface_nonFC(const char * input_data_directory)  
{
      gROOT->Reset();

      // Define dm23sq, sinsq2theta23 range and step  
      double sinsq2theta23_min   = 0.0;
      double sinsq2theta23_max   = 1.0;
      double sinsq2theta23_step  = 0.01;
      double dm23sq_min          = 1.0E-3;
      double dm23sq_max          = 6.0E-3;
      double dm23sq_step         = 1.0E-4;

      // Define grid so that min and max are at the bin centre 
      // and not at the low/high edge of the first/last bin 
      const int sinsq2theta23_nbins = 1 + TMath::CeilNint ((sinsq2theta23_max - sinsq2theta23_min) / sinsq2theta23_step);
      const int dm23sq_nbins        = 1 + TMath::CeilNint ((dm23sq_max        - dm23sq_min       ) / dm23sq_step       );

      double sinsq2theta23_grid_min = sinsq2theta23_min - 0.5*sinsq2theta23_step;
      double sinsq2theta23_grid_max = sinsq2theta23_max + 0.5*sinsq2theta23_step;
      double dm23sq_grid_min        = dm23sq_min        - 0.5*dm23sq_step;
      double dm23sq_grid_max        = dm23sq_max        + 0.5*dm23sq_step;


      // Book some histograms (chiSqSurface) 

      //chiSqSurface
      TH2D* hChi2SurfaceStat = new TH2D("hChi2SurfaceStat", "",
                  sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
      TH2D* hChi2SurfaceSyst = new TH2D("hChi2SurfaceSyst", "",
                  sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);
      //----
      //chiSqSurface - min
      TH2D* hChi2SurfaceStat_min = new TH2D("hChi2SurfaceStat_min", "",
                  sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);

      TH2D* hChi2SurfaceSyst_min = new TH2D("hChi2SurfaceSyst_min", "",
                  sinsq2theta23_nbins, sinsq2theta23_grid_min, sinsq2theta23_grid_max, dm23sq_nbins, dm23sq_grid_min, dm23sq_grid_max);

      //==========================

      // Load the root files
      TChain chain("clnt");
      chain.Add(Form("%s/*root", input_data_directory));


      //Fill in the chiSqSurface histograms
      for(int isinsq2theta23_bin = 1; isinsq2theta23_bin <= sinsq2theta23_nbins; isinsq2theta23_bin++) {


            for(int idm23sq_bin = 1; idm23sq_bin <= dm23sq_nbins; idm23sq_bin++){ 


                  double sinsq2theta23_low   = hChi2SurfaceStat->GetXaxis()->GetBinLowEdge(isinsq2theta23_bin);
                  double sinsq2theta23_high  = hChi2SurfaceStat->GetXaxis()->GetBinUpEdge(isinsq2theta23_bin);
                  double sinsq2theta23_center = hChi2SurfaceStat->GetXaxis()->GetBinCenter(isinsq2theta23_bin);
                  double dm23sq_low          = hChi2SurfaceStat->GetYaxis()->GetBinLowEdge(idm23sq_bin);
                  double dm23sq_high         = hChi2SurfaceStat->GetYaxis()->GetBinUpEdge(idm23sq_bin);
                  double dm23sq_center = hChi2SurfaceStat->GetYaxis()->GetBinCenter(idm23sq_bin);

                  //cut 
                  const char * cut = Form("input_sinsq2theta23 >= %f && input_sinsq2theta23 <= %f && input_dm23sq >= %f && input_dm23sq <= %f",
                              sinsq2theta23_low, sinsq2theta23_high, dm23sq_low, dm23sq_high);
                  cout << "Looking for fits in 2-D bin defined by: " << cut << endl;
                  //======================== 

                  double chi2   = 0.0;
                  double avChi2 =0.0;
                  double avChi2Stat =0.0;
                  double tmpv =-1.0;

                  //stat+syst chi2
                  chain.Draw("chisq_sys", cut, "goff");
                  int nfits = chain.GetSelectedRows();
                  if(nfits != 0) {
                        cout << "Found " << nfits << " entries in this bin" << endl;
                  }

                  if(nfits > 0) {
                        for(int j = 0; j < nfits; j++) {
                              tmpv = chain.GetV1()[j];
                              chi2=tmpv;
                        }
                        avChi2=chi2;// should be only one fit per bin!
                        hChi2SurfaceSyst->SetBinContent(isinsq2theta23_bin, idm23sq_bin, avChi2);    
                  }else{
                        hChi2SurfaceSyst->SetBinContent(isinsq2theta23_bin, idm23sq_bin, 0.0);    
                  }

                  //stat only chi2
                  chain.Draw("chisq_stat", cut, "goff");
                  int nfits = chain.GetSelectedRows();
                  if (nfits != 0) {
                        cout << "Found " << nfits << " entries in this bin" << endl;
                  }
                  chi2 =0.0;
                  tmpv =-1.0;

                  if(nfits > 0) {
                        for(int j = 0; j < 1; j++) {
                              tmpv = chain.GetV1()[j];
                              chi2=tmpv;
                        }
                        avChi2Stat=chi2;
                        hChi2SurfaceStat->SetBinContent(isinsq2theta23_bin, idm23sq_bin, avChi2Stat);    
                  }else{
                        hChi2SurfaceStat->SetBinContent(isinsq2theta23_bin, idm23sq_bin, 0.0);    
                  }


            }// Dm2 loop 

      }//sinSq2theta loop


      //-----------------------------------------------
      //Fill in the DeltaChiSqSurface histograms
      double chisq_min                    = 999999.0;
      double sinsq2theta23_bestfit        = -1;
      int    sinsq2theta23_bestfit_bin    = -10;
      double dm23sq_bestfit               = -1;
      int    dm23sq_bestfit_bin           = -10;


      //find the minimum in hChi2SurfaceStat
      for(int i = 1; i <= hChi2SurfaceStat->GetNbinsX(); i++) {
            for(int j = 1; j <= hChi2SurfaceStat->GetNbinsY(); j++) {
                  double chisq = hChi2SurfaceStat->GetBinContent(i,j);
                  if(chisq==0.0){
                        ;
                  }else{
                        if(chisq < chisq_min) {
                              chisq_min = chisq;
                              sinsq2theta23_bestfit = hChi2SurfaceStat->GetXaxis()->GetBinCenter(i);
                              dm23sq_bestfit        = hChi2SurfaceStat->GetYaxis()->GetBinCenter(j);
                              dm23sq_bestfit_bin = j;
                              sinsq2theta23_bestfit_bin = i;
                        }
                  }
            }//j loop
      }//i loop

      //fill in the hChi2SurfaceStat_min 
      for(int i = 1; i <= hChi2SurfaceStat_min->GetNbinsX(); i++) {
            for(int j = 1; j <= hChi2SurfaceStat_min->GetNbinsY(); j++) {
                  hChi2SurfaceStat_min->SetBinContent(i,j, chisq_min);
            }//j loop
      }//i loop

      cout << "chi\2_{min} = " << chisq_min << endl;
      cout << "dm23sq_{bestfit} = " << dm23sq_bestfit << "      ";
      cout << "dm23sq_{bestfit bin} = " << dm23sq_bestfit_bin << endl;
      cout << "sinsq2theta_{bestfit} = " << sinsq2theta23_bestfit << "      ";
      cout << "sinsq2theta_{bestfit bin} = " << sinsq2theta23_bestfit_bin << endl;


      chisq_min =999991.0; 

      //find the minimum in hChi2SurfaceSyst
      for(int i = 1; i <= hChi2SurfaceSyst->GetNbinsX(); i++) {
            for(int j = 1; j <= hChi2SurfaceSyst->GetNbinsY(); j++) {
                  double chisq = hChi2SurfaceSyst->GetBinContent(i,j);
                  if(chisq==0.0){
                        ;
                  }else{
                        if(chisq < chisq_min) {
                              chisq_min = chisq;
                              sinsq2theta23_bestfit = hChi2SurfaceSyst->GetXaxis()->GetBinCenter(i);
                              dm23sq_bestfit        = hChi2SurfaceSyst->GetYaxis()->GetBinCenter(j);
                              dm23sq_bestfit_bin = j;
                              sinsq2theta23_bestfit_bin = i;
                        }
                  }
            }//j loop
      }//i loop

      cout << "chi\2_{min stat + sys} = " << chisq_min << endl;
      cout << "dm23sq_{bestfit stat + sys} = " << dm23sq_bestfit << "      ";
      cout << "dm23sq_{bestfit bin stat + sys} = " << dm23sq_bestfit_bin << endl;
      cout << "sinsq2theta_{bestfit stat + sys} = " << sinsq2theta23_bestfit << "      ";
      cout << "sinsq2theta_{bestfit bin stat + sys} = " << sinsq2theta23_bestfit_bin << endl;

      //fill in the hChi2SurfaceSyst_min 
      for(int i = 1; i <= hChi2SurfaceSyst_min->GetNbinsX(); i++) {
            for(int j = 1; j <= hChi2SurfaceSyst_min->GetNbinsY(); j++) {
                  hChi2SurfaceSyst_min->SetBinContent(i,j, chisq_min);
            }//j loop
      }//i loop

      //=========================================================== 
      //write out the histos
      TFile * fout = new TFile("ChisqSurface.root", "RECREATE");
      hChi2SurfaceStat->Write();
      hChi2SurfaceSyst->Write();

      hChi2SurfaceStat_min->Write();
      hChi2SurfaceSyst_min->Write();
      fout->Close();


}
