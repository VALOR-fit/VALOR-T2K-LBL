void draw_average_spectrum(
			   const char * name = "joint_residuals_new_fix_point0_stats_marg_octant",
			   const char * folder = "",
			   const char * treename = "fit_tree",
			   bool average_data = true /* if true, average data spectra. otherwise, average best-fit spectra */
)
{

  using namespace std;
  using namespace valor;

  //set up the tchain
  TChain ch(treename);
  if(folder == "") {
    TString filename = Form("%s.root",name);
    ch.Add(filename);
  }
  else {
    TString filenames = Form("%s/*.root",folder);
    ch.Add(filenames);
  }
  FitParamList *  fparm = 0;
  PhysParamList * iparm = 0;
  TObjArray * dat_array = 0;
  TObjArray * pdf_array = 0;
  fitnt->SetBranchAddress("best_fit_params", &fparm);
  fitnt->SetBranchAddress("input_params",    &iparm);
  fitnt->SetBranchAddress("data",            &dat_array);
  fitnt->SetBranchAddress("best_fit_pdf",    &pdf_array);


  const char * outfile0 = Form("%s_average.root",name);
  TFile f0(outfile0,"recreate");

  gStyle->SetOptStat(1100);
  int N = ch.GetEntries();

  TH1D * spectrum_peak           = new TH1D("spectrum_peak", " Average spectrum in the max disappearance peak",
				       util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * spectrum_away           = new TH1D("spectrum_away", " Average spectrum away from the max disappearance peak",
					util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * spectrum_all            = new TH1D("spectrum_all",  " Average spectrum for all fits",
					util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));
 
  const int ntoys = ch.GetEntries();
  int ntoys_peak = 0;
  int ntoys_away = 0;
  double array_peak[util::NumErecoBins (conv::kSKsFCFV1RingMuLike)] = {0};
  double array_away[util::NumErecoBins (conv::kSKsFCFV1RingMuLike)] = {0};
  double array_all [util::NumErecoBins (conv::kSKsFCFV1RingMuLike)] = {0};
  for(int ie = 0; ie < ntoys; ie++) {
    if(ie % 1000 == 0)
      cout << "Getting entry " << ie << endl;
    ch.GetEntry(ie);
    TH1D * data1Rmu    = (TH1D*) dat_array->FindObject("obs_1Rmu");
    TH1D * bestfit1Rmu = (TH1D*) pdf_array->FindObject("pdf_1Rmu");
    TH1D * hist = 0;
    if(average_data)
      hist = data1Rmu;
    else
      hist = bestfit1Rmu;
    for(int ib = 0; ib < util::NumErecoBins (conv::kSKsFCFV1RingMuLike); ib++) {
      if(TMath::Abs(fparm->BestFitValue(conv::kPhypMixing23)-iparm->Value(conv::kPhypMixing23)) < 0.02) {
	array_peak[ib] += hist->GetBinContent(ib);
	ntoys_peak++;
      }
      else {
	array_away[ib] += hist->GetBinContent(ib);
	ntoys_away++;
      }
      array_all[ib] += hist->GetBinContent(ib);
    }//ib
  }//ie

  //average & fill histograms
  for(int ib = 0; ib < util::NumErecoBins (conv::kSKsFCFV1RingMuLike); ib++) {
    spectrum_peak->SetBinContent(ib+1, array_peak[ib] / ntoys_peak);
    spectrum_away->SetBinContent(ib+1, array_away[ib] / ntoys_away);
    spectrum_all ->SetBinContent(ib+1, array_all [ib] / ntoys);
  }

  //draw
  TCanvas * c_1 = new TCanvas("c_1","",20,20,1000,500);
  spectrum_peak->SetLineColor(kRed);
  spectrum_peak->GetXaxis()->SetRangeUser(0.0,7.0);
  spectrum_peak->Draw("");
  TCanvas * c_2 = new TCanvas("c_2","",20,20,1000,500);
  spectrum_away->SetLineColor(kBlue);
  spectrum_away->GetXaxis()->SetRangeUser(0.0,7.0);
  spectrum_away->Draw("");
  TCanvas * c_total = new TCanvas("c_total","",20,20,1000,500);
  spectrum_all->GetXaxis()->SetRangeUser(0.0,7.0);
  spectrum_all->Draw("");

  //save
  f0.cd();
  spectrum_peak->Write();
  spectrum_away->Write();
  spectrum_all->Write();
}
