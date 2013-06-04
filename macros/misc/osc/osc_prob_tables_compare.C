
void osc_prob_tables_compare()
{
  //Notice: make sure these ranges coincide with input root files 
  double kpi = 3.14159;
  double mixing13_min   = 0.0;
  double mixing13_max   = 0.4;
  double mixing13_step  = 0.004;
  double deltacp_min    = -1 * kpi;
  double deltacp_max    = kpi;
  double deltacp_step   = 2*kpi/100.;
  double mixing23_min   = 0.8;
  double mixing23_max   = 1.0;
  double mixing23_step  = 0.002;
  double dm23sq_min     = 1.5E-3;
  double dm23sq_max     = 3.5E-3;
  double dm23sq_step    = 2.0E-5;

  // Make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin                               
  int deltacp_nbins   = 1 + TMath::Nint( (deltacp_max    - deltacp_min ) / deltacp_step   );
  int mixing13_nbins  = 1 + TMath::Nint( (mixing13_max  - mixing13_min ) / mixing13_step  );
  int mixing23_nbins  = 1 + TMath::Nint( (mixing23_max  - mixing23_min ) / mixing23_step  );
  int dm23sq_nbins    = 1 + TMath::Nint( (dm23sq_max - dm23sq_min)       / dm23sq_step    );

  //cout << "deltacp_nbins=" << deltacp_nbins << endl;
  //cout << "dm23sq_nbins=" << dm23sq_nbins << endl;
  //cout << "mixing13_nbins=" << mixing13_nbins << endl;
  //cout << "mixing23_nbins=" << mixing23_nbins << endl;

  //energy
  double energy_max  = 10.000;
  double energy_min  = 0.250;
  int energy_nbins = 500;
  double energy_step = (energy_max - energy_min)/(energy_nbins-1);

  //define arrays to store results of reading
  double * nuetonue_dcp_first = new double[energy_nbins][deltacp_nbins];
  double * nuetonumu_dcp_first = new double[energy_nbins][deltacp_nbins];
  double * numutonue_dcp_first = new double[energy_nbins][deltacp_nbins];
  double * numutonumu_dcp_first = new double[energy_nbins][deltacp_nbins];

  double * nuetonue_dcp_second = new double[energy_nbins][deltacp_nbins];
  double * nuetonumu_dcp_second = new double[energy_nbins][deltacp_nbins];
  double * numutonue_dcp_second = new double[energy_nbins][deltacp_nbins];
  double * numutonumu_dcp_second = new double[energy_nbins][deltacp_nbins];

  double * nuetonue_mix13_first = new double[energy_nbins][mixing13_nbins];
  double * nuetonumu_mix13_first = new double[energy_nbins][mixing13_nbins];
  double * numutonue_mix13_first = new double[energy_nbins][mixing13_nbins];
  double * numutonumu_mix13_first = new double[energy_nbins][mixing13_nbins];

  double * nuetonue_mix13_second = new double[energy_nbins][mixing13_nbins];
  double * nuetonumu_mix13_second = new double[energy_nbins][mixing13_nbins];
  double * numutonue_mix13_second = new double[energy_nbins][mixing13_nbins];
  double * numutonumu_mix13_second = new double[energy_nbins][mixing13_nbins];


  //histograms to store the difference
  TH2D *NuEToNuE3f_diff_dcp         =  new TH2D("NuEToNuE3f_diff_dcp","3 Flavour P_{#nu_{e}#rightarrow#nu_{e}}",
                                           energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
                                           deltacp_nbins, deltacp_min - 0.5*deltacp_step, deltacp_max + 0.5*deltacp_step);
  TH2D *NuEToNuMu3f_diff_dcp        =  new TH2D("NuEToNuMu3f_diff_dcp","3 Flavour P_{#nu_{e}#rightarrow#nu_{#mu}}",
                                           energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
                                           deltacp_nbins, deltacp_min - 0.5*deltacp_step, deltacp_max + 0.5*deltacp_step);
  TH2D *NuMuToNuE3f_diff_dcp        =  new TH2D("NuMuToNuE3f_diff_dcp","3 Flavour P_{#nu_{#mu}#rightarrow#nu_{e}}",
                                           energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
                                           deltacp_nbins, deltacp_min - 0.5*deltacp_step, deltacp_max + 0.5*deltacp_step);
  TH2D *NuMuToNuMu3f_diff_dcp       =  new TH2D("NuMuToNuMu3f_diff_dcp","3 Flavour P_{#nu_{#mu}#rightarrow#nu_{#mu}}",
                                           energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
                                           deltacp_nbins, deltacp_min - 0.5*deltacp_step, deltacp_max + 0.5*deltacp_step);
  //mix13 vs energy surfaces                 
  TH2D *NuEToNuE3f_diff_mix13       =  new TH2D("NuEToNuE3f_diff_mix13","3 Flavour P_{#nu_{e}#rightarrow#nu_{e}}",
                                           energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
                                           mixing13_nbins, mixing13_min - 0.5*mixing13_step, mixing13_max + 0.5*mixing13_step);
  TH2D *NuEToNuMu3f_diff_mix13      =  new TH2D("NuEToNuMu3f_diff_mix13","3 Flavour P_{#nu_{e}#rightarrow#nu_{#mu}}",
                                           energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
                                           mixing13_nbins, mixing13_min - 0.5*mixing13_step, mixing13_max + 0.5*mixing13_step);
  TH2D *NuMuToNuE3f_diff_mix13      =  new TH2D("NuMuToNuE3f_diff_mix13","3 Flavour P_{#nu_{#mu}#rightarrow#nu_{e}}",
                                           energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
                                           mixing13_nbins, mixing13_min - 0.5*mixing13_step, mixing13_max + 0.5*mixing13_step);
  TH2D *NuMuToNuMu3f_diff_mix13     =  new TH2D("NuMuToNuMu3f_diff_mix13","3 Flavour P_{#nu_{#mu}#rightarrow#nu_{#mu}}",
                                           energy_nbins,  energy_min -0.5*energy_step , energy_max + 0.5*energy_step,
                                           mixing13_nbins, mixing13_min - 0.5*mixing13_step, mixing13_max + 0.5*mixing13_step);


  //delta_cp vs sinsq2theta13 surfaces to store the difference
                                                               
  TH2D *NuEToNuE3f_diff             =  new TH2D("NuEToNuE3f_diff","P_{#nu_{e}#rightarrow#nu_{e}}",
                                           mixing13_nbins, mixing13_min - 0.5*mixing13_step, mixing13_max + 0.5*mixing13_step,
                                           deltacp_nbins, deltacp_min - 0.5*deltacp_step, deltacp_max + 0.5*deltacp_step);

  TH2D *NuEBarToNuEBar3f_diff       =  new TH2D("NuEBarToNuEBar3f_diff","P_{#bar{#nu_{e}}#rightarrow#bar{#nu_{e}}}",
                                           mixing13_nbins, mixing13_min - 0.5*mixing13_step, mixing13_max + 0.5*mixing13_step,
                                           deltacp_nbins, deltacp_min - 0.5*deltacp_step, deltacp_max + 0.5*deltacp_step);

  TH2D *NuMuToNuE3f_diff            =  new TH2D("NuMuToNuE3f_diff","P_{#nu_{#mu}#rightarrow#nu_{e}}",
                                           mixing13_nbins, mixing13_min - 0.5*mixing13_step, mixing13_max + 0.5*mixing13_step,
                                           deltacp_nbins, deltacp_min - 0.5*deltacp_step, deltacp_max + 0.5*deltacp_step);

  //dm23sq vs sinsq2theta23 surfaces                                                                                                                                                                          
  TH2D *NuMuToNuMu3f_diff           =  new TH2D("NuMuToNuMu3f_diff","P_{#nu_{#mu}#rightarrow#nu_{#mu}}",
                                           mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step,
                                           dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);

  TH2D *NuMuBarToNuMuBar3f_diff     =  new TH2D("NuMuBarToNuMuBar3f_diff","P_{#bar{#nu_{#mu}}#rightarrow#bar{#nu_{#mu}}}",
                                           mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step,
                                           dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);

  TH2D *NuEToNuMu3f_diff            =  new TH2D("NuEToNuMu3f_diff","P_{#nu_{e}#rightarrow#nu_{#mu}}",
                                           mixing23_nbins, mixing23_min - 0.5*mixing23_step, mixing23_max + 0.5*mixing23_step,
                                           dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);




  //open files
  //these two files contain the oscillation tables to be compared
  //the first one is generated with macros/misc/generate_osc_tables.C
  //the second one is the output file from Prob3++
  //Be sure the ranges of osc parameters are the same in both
  TFile *_file0 = TFile::Open("OurOscProbTable_inverted.root"); // with our code
  TFile *_file1 = TFile::Open("OscProbTables_inverted.root");   // Prob3++
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  /*
  for(int jbin = 1; jbin <= deltacp_nbins; jbin++) {
    double dcp = deltacp_min + (jbin-1)*deltacp_step;
    for (int ibin = 1; ibin<= energy_nbins; ibin++){
      double energy = energy_min + (ibin-1)*energy_step;
      _file0->cd();
      nuetonue_dcp_first[ibin][jbin]   = NuEToNuE3f_dcp->GetBinContent(ibin,jbin);
      nuetonumu_dcp_first[ibin][jbin]  = NuEToNuMu3f_dcp->GetBinContent(ibin,jbin);
      numutonue_dcp_first[ibin][jbin]  = NuMuToNuE3f_dcp->GetBinContent(ibin,jbin);
      numutonumu_dcp_first[ibin][jbin] = NuMuToNuMu3f_dcp->GetBinContent(ibin,jbin);

      _file1->cd();
      nuetonue_dcp_second[ibin][jbin]   = NuEToNuE3f_dcp->GetBinContent(ibin,jbin);
      nuetonumu_dcp_second[ibin][jbin]  = NuEToNuMu3f_dcp->GetBinContent(ibin,jbin);
      numutonue_dcp_second[ibin][jbin]  = NuMuToNuE3f_dcp->GetBinContent(ibin,jbin);
      numutonumu_dcp_second[ibin][jbin] = NuMuToNuMu3f_dcp->GetBinContent(ibin,jbin);

      // fill the histograms with the difference
      NuEToNuE3f_diff_dcp  ->Fill(energy , dcp  , (nuetonue_dcp_first[ibin][jbin]- nuetonue_dcp_second[ibin][jbin]));
      NuEToNuMu3f_diff_dcp ->Fill(energy , dcp  , (nuetonumu_dcp_first[ibin][jbin]- nuetonumu_dcp_second[ibin][jbin]));
      NuMuToNuE3f_diff_dcp ->Fill(energy , dcp  , (numutonue_dcp_first[ibin][jbin]- numutonue_dcp_second[ibin][jbin]));
      NuMuToNuMu3f_diff_dcp->Fill(energy , dcp  , (numutonumu_dcp_first[ibin][jbin]- numutonumu_dcp_second[ibin][jbin]));
    }
  }

  for(int jbin = 1; jbin <= mixing13_nbins; jbin++) {
    double mix13 = mixing13_min + (jbin-1)*mixing13_step;
    for (int ibin = 1; ibin<= energy_nbins; ibin++){
      double energy = energy_min + (ibin-1)*energy_step;
      _file0->cd();
      nuetonue_mix13_first[ibin][jbin]   = NuEToNuE3f_mix13->GetBinContent(ibin,jbin);
      nuetonumu_mix13_first[ibin][jbin]  = NuEToNuMu3f_mix13->GetBinContent(ibin,jbin);
      numutonue_mix13_first[ibin][jbin]  = NuMuToNuE3f_mix13->GetBinContent(ibin,jbin);
      numutonumu_mix13_first[ibin][jbin] = NuMuToNuMu3f_mix13->GetBinContent(ibin,jbin);

      _file1->cd();
      nuetonue_mix13_second[ibin][jbin]   = NuEToNuE3f_mix13->GetBinContent(ibin,jbin);
      nuetonumu_mix13_second[ibin][jbin]  = NuEToNuMu3f_mix13->GetBinContent(ibin,jbin);
      numutonue_mix13_second[ibin][jbin]  = NuMuToNuE3f_mix13->GetBinContent(ibin,jbin);
      numutonumu_mix13_second[ibin][jbin] = NuMuToNuMu3f_mix13->GetBinContent(ibin,jbin);

      // fill the histograms with the difference
      NuEToNuE3f_diff_mix13  ->Fill(energy , mix13  , (nuetonue_mix13_first[ibin][jbin]- nuetonue_mix13_second[ibin][jbin]));
      NuEToNuMu3f_diff_mix13 ->Fill(energy , mix13  , (nuetonumu_mix13_first[ibin][jbin]- nuetonumu_mix13_second[ibin][jbin]));
      NuMuToNuE3f_diff_mix13 ->Fill(energy , mix13  , (numutonue_mix13_first[ibin][jbin]- numutonue_mix13_second[ibin][jbin]));
      NuMuToNuMu3f_diff_mix13->Fill(energy , mix13  , (numutonumu_mix13_first[ibin][jbin]- numutonumu_mix13_second[ibin][jbin]));
    }
  }

  */


  for(int jbin = 1; jbin <= deltacp_nbins; jbin++) {
    double dcp = deltacp_min + (jbin-1)*deltacp_step;
    for (int ibin = 1; ibin<= mixing13_nbins; ibin++){
      double mix13 = mixing13_min + (ibin-1)*mixing13_step;

      _file0->cd();
      double nuetonue_first         = NuEToNuE3f->GetBinContent(ibin,jbin);
      double nuebartonuebar_first   = NuEBarToNuEBar3f->GetBinContent(ibin,jbin);
      double numutonue_first        = NuMuToNuE3f->GetBinContent(ibin,jbin);
      
      _file1->cd();
      double nuetonue_second         = NuEToNuE3f->GetBinContent(ibin,jbin);
      double nuebartonuebar_second   = NuEBarToNuEBar3f->GetBinContent(ibin,jbin);
      double numutonue_second        = NuMuToNuE3f->GetBinContent(ibin,jbin);
      
      // fill the histograms with the difference
      // it will be fractional difference, so difference divided by reference value
      // it will take our result as reference value

      NuEToNuE3f_diff  ->Fill(mix13 , dcp  , (nuetonue_first - nuetonue_second)/nuetonue_first );
      NuEBarToNuEBar3f_diff  ->Fill(mix13 , dcp  , (nuebartonuebar_first - nuebartonuebar_second)/nuebartonuebar_first );
      NuMuToNuE3f_diff ->Fill(mix13 , dcp  , (numutonue_first - numutonue_second)/numutonue_first );
      
    }
  }

  for(int jbin = 1; jbin <= dm23sq_nbins; jbin++) {
    double dm23sq = dm23sq_min + (jbin-1)*dm23sq_step;
    for (int ibin = 1; ibin<= mixing23_nbins; ibin++){
      double mix23 = mixing23_min + (ibin-1)*mixing23_step;

      _file0->cd();
      double nuetonumu_first        = NuEToNuMu3f->GetBinContent(ibin,jbin);
      double numutonumu_first       = NuMuToNuMu3f->GetBinContent(ibin,jbin);
      double numubartonumubar_first = NuMuBarToNuMuBar3f->GetBinContent(ibin,jbin);

      _file1->cd();
      double nuetonumu_second        = NuEToNuMu3f->GetBinContent(ibin,jbin);
      double numutonumu_second       = NuMuToNuMu3f->GetBinContent(ibin,jbin);
      double numubartonumubar_second = NuMuBarToNuMuBar3f->GetBinContent(ibin,jbin);


      // fill the histograms with the difference
      // it will be fractional difference, so difference divided by reference value
      // it will take our result as reference value

      NuEToNuMu3f_diff ->Fill(mix23 , dm23sq  , (nuetonumu_first - nuetonumu_second)/nuetonumu_first );
      NuMuToNuMu3f_diff->Fill(mix23 , dm23sq  , (numutonumu_first - numutonumu_second)/numutonumu_first );
      NuMuBarToNuMuBar3f_diff->Fill(mix23 , dm23sq  , (numubartonumubar_first - numubartonumubar_second)/numubartonumubar_first );
      
    }
  }

  NuEToNuE3f_diff      ->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  NuEBarToNuEBar3f_diff->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");
  NuMuToNuE3f_diff     ->GetXaxis()->SetTitle("sin^{2}2#theta_{13}");

  NuEToNuE3f_diff      ->GetYaxis()->SetTitle("#delta_{CP}");
  NuEBarToNuEBar3f_diff->GetYaxis()->SetTitle("#delta_{CP}");
  NuMuToNuE3f_diff     ->GetYaxis()->SetTitle("#delta_{CP}");

  NuEToNuE3f_diff      ->GetZaxis()->SetTitle("FRACTIONAL DIFFERENCE");
  NuEBarToNuEBar3f_diff->GetZaxis()->SetTitle("FRACTIONAL DIFFERENCE");
  NuMuToNuE3f_diff     ->GetZaxis()->SetTitle("FRACTIONAL DIFFERENCE");

  NuMuToNuMu3f_diff          ->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
  NuMuBarToNuMuBar3f_diff    ->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");
  NuEToNuMu3f_diff           ->GetXaxis()->SetTitle("sin^{2}2#theta_{23}");

  NuMuToNuMu3f_diff          ->GetYaxis()->SetTitle("|#Deltam^{2}_{23}|");
  NuMuBarToNuMuBar3f_diff    ->GetYaxis()->SetTitle("|#Deltam^{2}_{23}|");
  NuEToNuMu3f_diff           ->GetYaxis()->SetTitle("|#Deltam^{2}_{23}|");

  NuMuToNuMu3f_diff          ->GetZaxis()->SetTitle("FRACTIONAL DIFFERENCE");
  NuMuBarToNuMuBar3f_diff    ->GetZaxis()->SetTitle("FRACTIONAL DIFFERENCE");
  NuEToNuMu3f_diff           ->GetZaxis()->SetTitle("FRACTIONAL DIFFERENCE");


  TFile f("diff_osc_prob_tables_inverted.root", "recreate");
  f.cd();

  NuEToNuE3f_diff     ->Write();
  NuEBarToNuEBar3f_diff     ->Write();
  NuEToNuMu3f_diff    ->Write();
  NuMuToNuE3f_diff    ->Write();
  NuMuToNuMu3f_diff   ->Write();
  NuMuBarToNuMuBar3f_diff   ->Write();

  f.Write();
  f.Close();
}
