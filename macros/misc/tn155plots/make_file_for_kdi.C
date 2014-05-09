//collect the various files we have & give KDI only the information we need for data comparison plots

void make_file_for_kdi(
		       const char * tag = "run1-4full_skdetmathfix_run1-4abanff_nuis4_sintheta_inverted31",
		       int run = 0,
		       int expt = 8,
		       bool smooth = true,
		       const char * tag2 = "run1-4full_nuis4_sintheta_inverted31"
)
{
  using namespace valor;

  if(tag2 == "")
    tag2 = tag;

  TFile fpull(Form("%s_pulls.root",tag));
  TFile fall(Form("%s-run_%d-expt_%d.root",tag,run,expt));
  TFile flimit(Form("$T2KOSCANA/2013results/DATA/PROFILES/%s/combination/final.root",tag2));

  TH1D * pulls = 0;
  fpull.GetObject("pulls",pulls);
  TH1D * pdf_1Rmu = 0;
  fall.GetObject("pdf_1Rmu",pdf_1Rmu);
  TH2D * dchisq = 0;
  fall.GetObject("dchisq",dchisq);
  TH2D * dchisq_cont68 = 0;
  fall.GetObject("dchisq_cont68",dchisq_cont68);
  TH2D * dchisq_cont90 = 0;
  fall.GetObject("dchisq_cont90",dchisq_cont90);
  FitParamList * fparm = 0;
  fall.GetObject("fparm",fparm);
  TH1D * profile_mix = 0;
  flimit.GetObject("chisq_mixing23",profile_mix);
  double chisq_min = profile_mix->GetMinimum();
  for(int i = 1; i <= profile_mix->GetNbinsX(); i++)
    profile_mix->SetBinContent(i, profile_mix->GetBinContent(i) - chisq_min);
  TH1D * profile_dm = 0;
  flimit.GetObject("chisq_dm32sq",profile_dm);
  double chisq_min = profile_dm->GetMinimum();
  for(int i = 1; i <= profile_dm->GetNbinsX(); i++)
    profile_dm->SetBinContent(i, profile_dm->GetBinContent(i) - chisq_min);


  TMarker bf(fparm->BestFitValue(conv::kPhypMixing23),fparm->BestFitValue(conv::kPhypDm23sq),29);

  if(smooth) {
    dchisq->Smooth();
    dchisq_cont68->Smooth();
    dchisq_cont90->Smooth();
  }

  TFile fout(Form("%s_hist.root",tag),"RECREATE");
  fout.cd();
  pulls->Write("pulls");
  pdf_1Rmu->Write("pdf_1Rmu");
  dchisq->Write("dchisq");
  dchisq_cont68->Write("dchisq_cont68");
  dchisq_cont90->Write("dchisq_cont90");
  bf.Write("bestfit");
  profile_mix->Write("profile_mixing23");
  profile_dm->Write("profile_dm32sw");

}
