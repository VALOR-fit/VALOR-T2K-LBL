void draw_average_spectrum(
  const char * name = "joint_residuals_new_fix_point0_stats_marg_octant"
)
{

  using namespace std;
  using namespace valor;

  const char * filename = Form("%s.root",name);
  TFile * f = new TFile(filename,"read");
  TTree * tree ;
  f->GetObject("fit_tree",tree);

  const char * outfile0 = Form("%s_average.root",name);
  TFile f0(outfile0,"recreate");

  gStyle->SetOptStat(1100);
  int N = tree->GetEntries();

  TH1D * spectrum_1           = new TH1D("spectrum_1", " Average spectrum ",
				       util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * spectrum_2           = new TH1D("spectrum_2", " Average spectrum ",
					util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));

  TH1D * spectrum_total       = new TH1D("spectrum_total", " Average spectrum ",
					util::NumErecoBins (conv::kSKsFCFV1RingMuLike), util::ErecoBinEdges(conv::kSKsFCFV1RingMuLike));
 

  //now draw osc params     
  for ( int i = 0; i<util::NumErecoBins (conv::kSKsFCFV1RingMuLike); i++){
    if(i<10) tree->Draw(Form("bin0%d",i),"abs(bestfit1_mixing23-input_mixing23)<0.02","");
    else  tree->Draw(Form("bin%d",i),"abs(bestfit1_mixing23-input_mixing23)<0.02","");
    int ndraw = tree->GetSelectedRows();
    double average = 0.;
    double *v = tree->GetV1();
    for (int j = 0; j<ndraw; j++){
      average +=  v[j];
    }
    average /= ndraw;
    spectrum_1->SetBinContent(i+1,average);
  }

  for ( int i = 0; i<util::NumErecoBins (conv::kSKsFCFV1RingMuLike); i++){
    if(i<10) tree->Draw(Form("bin0%d",i),"abs(bestfit1_mixing23-input_mixing23)>0.02","");
    else  tree->Draw(Form("bin%d",i),"abs(bestfit1_mixing23-input_mixing23)>0.02","");
    int ndraw = tree->GetSelectedRows();
    double average = 0.;
    double *v = tree->GetV1();
    for (int j = 0; j<ndraw; j++){
      average +=  v[j];
    }
    average /= ndraw;
    spectrum_2->SetBinContent(i+1,average);
  }

  for ( int i = 0; i<util::NumErecoBins (conv::kSKsFCFV1RingMuLike); i++){
    if(i<10) tree->Draw(Form("bin0%d",i),"","");
    else  tree->Draw(Form("bin%d",i),"","");
    int ndraw = tree->GetSelectedRows();
    double average = 0.;
    double *v = tree->GetV1();
    for (int j = 0; j<ndraw; j++){
      average +=  v[j];
    }
    average /= ndraw;
    spectrum_total->SetBinContent(i+1,average);
  }

  TCanvas * c_1 = new TCanvas("c_1","",20,20,1000,500);
  spectrum_1->SetLineColor(kRed);
  spectrum_1->GetXaxis()->SetRangeUser(0.0,7.0);
  spectrum_1->Draw("");
  TCanvas * c_2 = new TCanvas("c_2","",20,20,1000,500);
  spectrum_2->SetLineColor(kBlue);
  spectrum_2->GetXaxis()->SetRangeUser(0.0,7.0);
  spectrum_2->Draw("");
  TCanvas * c_total = new TCanvas("c_total","",20,20,1000,500);
  spectrum_total->GetXaxis()->SetRangeUser(0.0,7.0);
  spectrum_total->Draw("");

  f0.cd();
  spectrum_1->Write();
  spectrum_2->Write();
  spectrum_total->Write();
}
