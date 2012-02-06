//
// --------------------------------------------------------------------------------------------------------
// Build the marginalised chisq histograms from the merged output of build_chisq_surface_4d.C jobs.
// Each such job would have compute \chi^{2}, \chi^{2}_e and \chi^{2}_mu for part of the oscillation parameter space.
//
// Inputs
//  - inputname    : input file name 
//  - outputfile : name of output ROOT file with merged \chi^{2} surfaces
//              
// C.Andreopoulos, L.Escudero
// --------------------------------------------------------------------------------------------------------
//
 
void build_marginalised_chisq_4d(
     int dm23sq_min_bin=-1,
     int dm23sq_max_bin=999999,
     const char * inputname = "merged_chisq_4d_surface.root",
     const char * outputfile = "marginalised_chisq_4d.root")
{
  //
  // Define grid
  //
    
  using namespace t2kana;//just for kpi

  double dcp_min     = -1 * conv::kpi;
  double dcp_max     = conv::kpi;
  double dcp_step    = 2*conv::kpi/90.;
  double mix13_min   = 0.0;
  double mix13_fix   = 0.02;//fixed double chooz value
  double mix13_max   = 0.12;
  double mix13_step  = 0.0012;
  double mix23_min   = 0.3;
  double mix23_fix   = 0.43; //fixed value with numu disap best fit
  double mix23_max   = 0.7;
  double mix23_step  = 0.005;
  double dm23sq_min  = 1.0E-3;
  double dm23sq_fix  = 2.4E-3;//fixed value
  double dm23sq_max  = 5.0E-3;
  double dm23sq_step = 1E-4;

  // Make sure that min and max are at the bin centre and not at the low/high edge of the first/last bin                                              
  int dcp_nbins    = 1 + TMath::Nint( (dcp_max    - dcp_min   ) / dcp_step    );
  int mix13_nbins  = 1 + TMath::Nint( (mix13_max  - mix13_min ) / mix13_step  );
  int mix23_nbins  = 1 + TMath::Nint( (mix23_max  - mix23_min ) / mix23_step  );
  int dm23sq_nbins = 1 + TMath::Nint( (dm23sq_max - dm23sq_min) / dm23sq_step );
  int nentries = mix13_nbins * mix23_nbins ;

  //1D histograms for 1D marginalised
  TH1D *  dm23sq_marg_chisq    = new TH1D("dm23sq_marg_chisq", "Marginalised #chi^{2} for |#Deltam^{2}_{23}|",
					  dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);

  TH1D *  dm23sq_marg_chisqmu  = new TH1D("dm23sq_marg_chisqmu", "Marginalised #chi^{2} 1Rmu for |#Deltam^{2}_{23}|",
					  dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);

  TH1D *  dm23sq_marg_chisqe   = new TH1D("dm23sq_marg_chisqe", "Marginalised #chi^{2} 1Re for |#Deltam^{2}_{23}|",
					  dm23sq_nbins, dm23sq_min - 0.5*dm23sq_step, dm23sq_max + 0.5*dm23sq_step);


  TH1D *  dcp_marg_chisq       = new TH1D("dcp_marg_chisq", "Marginalised #chi^{2} for #delta_{CP}",
					  dcp_nbins, dcp_min - 0.5*dcp_step, dcp_max + 0.5*dcp_step);

  TH1D *  dcp_marg_chisqmu     = new TH1D("dcp_marg_chisqmu", "Marginalised #chi^{2} 1Rmu for #delta_{CP}",
					  dcp_nbins, dcp_min - 0.5*dcp_step, dcp_max + 0.5*dcp_step);

  TH1D *  dcp_marg_chisqe      = new TH1D("dcp_marg_chisqe", "Marginalised #chi^{2} 1Re for #delta_{CP}",
					  dcp_nbins, dcp_min - 0.5*dcp_step, dcp_max + 0.5*dcp_step);


  TH1D *  mix13_marg_chisq     = new TH1D("mix13_marg_chisq", "Marginalised #chi^{2} for sin^{2}#theta_{13}",
					  mix13_nbins, mix13_min - 0.5*mix13_step, mix13_max + 0.5*mix13_step);

  TH1D *  mix13_marg_chisqmu   = new TH1D("mix13_marg_chisqmu", "Marginalised #chi^{2} 1Rmu for sin^{2}#theta_{13}",
					  mix13_nbins, mix13_min - 0.5*mix13_step, mix13_max + 0.5*mix13_step);

  TH1D *  mix13_marg_chisqe    = new TH1D("mix13_marg_chisqe", "Marginalised #chi^{2} 1Re for sin^{2}#theta_{13}",
					  mix13_nbins, mix13_min - 0.5*mix13_step, mix13_max + 0.5*mix13_step);


  TH1D *  mix23_marg_chisq     = new TH1D("mix23_marg_chisq", "Marginalised #chi^{2} for sin^{2}#theta_{23}",
					  mix23_nbins, mix23_min - 0.5*mix23_step, mix23_max + 0.5*mix23_step);

  TH1D *  mix23_marg_chisqmu   = new TH1D("mix23_marg_chisqmu", "Marginalised #chi^{2} 1Rmu for sin^{2}#theta_{23}",
					  mix23_nbins, mix23_min - 0.5*mix23_step, mix23_max + 0.5*mix23_step);

  TH1D *  mix23_marg_chisqe    = new TH1D("mix23_marg_chisqe", "Marginalised #chi^{2} 1Re for sin^{2}#theta_{23}",
					  mix23_nbins, mix23_min - 0.5*mix23_step, mix23_max + 0.5*mix23_step);


  //write titles in 1D histograms
  dm23sq_marg_chisq      ->GetXaxis()->SetTitle("|#Deltam^{2}_{23}| (eV^{2}/c^{4})");
  dm23sq_marg_chisqmu    ->GetXaxis()->SetTitle("|#Deltam^{2}_{23}| (eV^{2}/c^{4})");
  dm23sq_marg_chisqe     ->GetXaxis()->SetTitle("|#Deltam^{2}_{23}| (eV^{2}/c^{4})");

  dcp_marg_chisq         ->GetXaxis()->SetTitle("#delta_{CP}");
  dcp_marg_chisqmu       ->GetXaxis()->SetTitle("#delta_{CP}");
  dcp_marg_chisqe        ->GetXaxis()->SetTitle("#delta_{CP}");


  mix13_marg_chisq       ->GetXaxis()->SetTitle("sin^{2}#theta_{13}");
  mix13_marg_chisqmu     ->GetXaxis()->SetTitle("sin^{2}#theta_{13}");
  mix13_marg_chisqe      ->GetXaxis()->SetTitle("sin^{2}#theta_{13}");

  mix23_marg_chisq       ->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  mix23_marg_chisqmu     ->GetXaxis()->SetTitle("sin^{2}#theta_{23}");
  mix23_marg_chisqe      ->GetXaxis()->SetTitle("sin^{2}#theta_{23}");

  // Load input file and ntuple
  TFile* infile = TFile::Open(inputname, "read"); // open the file
  TNtupleD *ntuple = (TNtupleD*)infile->Get("ntchisq");
  
  //loop over dcp bins
  for(int ibin = 1; ibin <= dcp_nbins; ibin++) {
    double dcp = dcp_min + (ibin-1)*dcp_step;
    //    const char* cut = Form("chisq!=0 && chisq1Rmu !=0 && chisq1Re !=0 && abs(dcp-%f)<%f && abs(dm23sq-%f)<%f && abs(mix23-%f)<%f && abs(mix13-%f)<%f", dcp , 1.E-3, dm23sq_fix , 1.E-5, mix23_fix , 1.E-4, mix13_fix, 1.E-04 );
    const char* cut = Form("chisq!=0 && chisq1Rmu !=0 && chisq1Re !=0 && abs(dcp-%f)<%f", dcp , 1.E-3);
    
    //chisq
    ntuple->Draw("chisq", cut, "goff");
    double * chisqarray = ntuple->GetV1();
    Int_t idx[999999];
    TMath::Sort(ntuple->GetSelectedRows(), chisqarray, idx, false);
    dcp_marg_chisq->SetBinContent(ibin,chisqarray[idx[0]]);
    
    //chisq 1Rmu
    ntuple->Draw("chisq1Rmu", cut, "goff");
    double * chisqarray_mu = ntuple->GetV1();
    Int_t idx_mu[999999];
    TMath::Sort(ntuple->GetSelectedRows(), chisqarray_mu, idx_mu, false);
    //    dcp_marg_chisqmu->SetBinContent(ibin,chisqarray_mu[idx_mu[0]]); //minimum for chisqmu
    dcp_marg_chisqmu->SetBinContent(ibin,chisqarray_mu[idx[0]]); //minimum for total chisq
    

    //chisq 1Re
    ntuple->Draw("chisq1Re", cut, "goff");
    double * chisqarray_e = ntuple->GetV1();
    Int_t idx_e[999999];
    TMath::Sort(ntuple->GetSelectedRows(), chisqarray_e, idx_e, false);
    //    dcp_marg_chisqe->SetBinContent(ibin,chisqarray_e[idx_e[0]]); //minimum for chisqe
    dcp_marg_chisqe->SetBinContent(ibin,chisqarray_e[idx[0]]); //minimum for total chisq

  } // end loop ibin

  
  //loop for dm23sq bins
  for(int ibin = 1; ibin <= dm23sq_nbins; ibin++) {
    double dm23sq = dm23sq_min + (ibin-1)*dm23sq_step;
    const char* cut = Form("chisq!=0 && chisq1Rmu !=0 && chisq1Re !=0 && abs(dm23sq-%f)<%f", dm23sq , 1.E-5);

    //chisq
    ntuple->Draw("chisq", cut, "goff");
    double * chisqarray = ntuple->GetV1();
    Int_t idx[999999];
    TMath::Sort(ntuple->GetSelectedRows(), chisqarray, idx, false);
    dm23sq_marg_chisq->SetBinContent(ibin,chisqarray[idx[0]]);

    //chisq 1Rmu
    ntuple->Draw("chisq1Rmu", cut, "goff");
    double * chisqarray_mu = ntuple->GetV1();
    Int_t idx_mu[999999];
    TMath::Sort(ntuple->GetSelectedRows(), chisqarray_mu, idx_mu, false);
    //    dm23sq_marg_chisqmu->SetBinContent(ibin,chisqarray_mu[idx_mu[0]]); //minimum for chisqmu
    dm23sq_marg_chisqmu->SetBinContent(ibin,chisqarray_mu[idx[0]]);//minimum for total chisq

    //chisq 1Re
    ntuple->Draw("chisq1Re", cut, "goff");
    double * chisqarray_e = ntuple->GetV1();
    Int_t idx_e[999999];
    TMath::Sort(ntuple->GetSelectedRows(), chisqarray_e, idx_e, false);
    //    dm23sq_marg_chisqe->SetBinContent(ibin,chisqarray_e[idx_e[0]]); //minimum for chisqe
    dm23sq_marg_chisqe->SetBinContent(ibin,chisqarray_e[idx[0]]); //minimum for total chisq

  } // end loop ibin
  
  
  //loop for mix23 bins
  for(int ibin = 1; ibin <= mix23_nbins; ibin++) {
    double mix23 = mix23_min + (ibin-1)*mix23_step;
    //    const char* cut = Form("chisq!=0 && chisq1Rmu !=0 && chisq1Re !=0 && abs(mix23-%f)<%f && abs(mix13-%f)<%f", mix23 , 1.E-4, mix13_fix, 1.E-04);
    const char* cut = Form("chisq!=0 && chisq1Rmu !=0 && chisq1Re !=0 && abs(mix23-%f)<%f", mix23 , 1.E-4);
    
    //chisq
    ntuple->Draw("chisq", cut, "goff");
    double * chisqarray = ntuple->GetV1();
    Int_t idx[999999];
    TMath::Sort(ntuple->GetSelectedRows(), chisqarray, idx, false);
    mix23_marg_chisq->SetBinContent(ibin,chisqarray[idx[0]]);
    
    //chisq 1Rmu
    ntuple->Draw("chisq1Rmu", cut, "goff");
    double * chisqarray_mu = ntuple->GetV1();
    Int_t idx_mu[999999];
    TMath::Sort(ntuple->GetSelectedRows(), chisqarray_mu, idx_mu, false);
    //    mix23_marg_chisqmu->SetBinContent(ibin,chisqarray_mu[idx_mu[0]]); //minimum for chisqmu
    mix23_marg_chisqmu->SetBinContent(ibin,chisqarray_mu[idx[0]]);//minimum for total chisq
    
    //chisq 1Re
    ntuple->Draw("chisq1Re", cut, "goff");
    double * chisqarray_e = ntuple->GetV1();
    Int_t idx_e[999999];
    TMath::Sort(ntuple->GetSelectedRows(), chisqarray_e, idx_e, false);
    //    mix23_marg_chisqe->SetBinContent(ibin,chisqarray_e[idx_e[0]]); //minimum for chisqe
    mix23_marg_chisqe->SetBinContent(ibin,chisqarray_e[idx[0]]); //minimum for total chisq

  } // end loop ibin

  
  //loop for mix13 bins
  for(int ibin = 1; ibin <= mix13_nbins; ibin++) {
    double mix13 = mix13_min + (ibin-1)*mix13_step;
    //    const char* cut = Form("chisq!=0 && chisq1Rmu !=0 && chisq1Re !=0 && abs(mix13-%f)<%f && abs(mix23-%f)<%f", mix13 , 1.E-4, mix23_fix, 1.E-4);
    const char* cut = Form("chisq!=0 && chisq1Rmu !=0 && chisq1Re !=0 && abs(mix13-%f)<%f", mix13 , 1.E-4);
    
    //chisq
    ntuple->Draw("chisq", cut, "goff");
    double * chisqarray = ntuple->GetV1();
    Int_t idx[999999];
    TMath::Sort(ntuple->GetSelectedRows(), chisqarray, idx, false);
    mix13_marg_chisq->SetBinContent(ibin,chisqarray[idx[0]]);
    
     
    //chisq 1Rmu
    ntuple->Draw("chisq1Rmu", cut, "goff");
    double * chisqarray_mu = ntuple->GetV1();
    Int_t idx_mu[999999];
    TMath::Sort(ntuple->GetSelectedRows(), chisqarray_mu, idx_mu, false);
    //    mix13_marg_chisqmu->SetBinContent(ibin,chisqarray_mu[idx_mu[0]]); //minimum for chisqmu
    mix13_marg_chisqmu->SetBinContent(ibin,chisqarray_mu[idx[0]]); //minimum for total chisq

    //chisq 1Re
    ntuple->Draw("chisq1Re", cut, "goff");
    double * chisqarray_e = ntuple->GetV1();
    Int_t idx_e[999999];
    TMath::Sort(ntuple->GetSelectedRows(), chisqarray_e, idx_e, false);
    //    mix13_marg_chisqe->SetBinContent(ibin,chisqarray_e[idx_e[0]]); //minimum for chisqe
    mix13_marg_chisqe->SetBinContent(ibin,chisqarray_e[idx[0]]); //minimum for total chisq

  } // end loop ibin


  //
  // Save histograms
  //
  
  cout << "Saving histograms to " << outputfile << endl;

  TFile * fout = new TFile(outputfile, "RECREATE");
  dm23sq_marg_chisq  ->Write();
  dm23sq_marg_chisqmu->Write();
  dm23sq_marg_chisqe ->Write();
  dcp_marg_chisq     ->Write();
  dcp_marg_chisqmu   ->Write();
  dcp_marg_chisqe    ->Write();
  mix13_marg_chisq   ->Write();
  mix13_marg_chisqmu ->Write();
  mix13_marg_chisqe  ->Write();
  mix23_marg_chisq   ->Write();
  mix23_marg_chisqmu ->Write();
  mix23_marg_chisqe  ->Write();
  
  fout->Close();
  
}
