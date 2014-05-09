void get_starplot_staterror(
			    const char * folder    = "./"
)
{
  using namespace valor;
  Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta;

  TChain c("fitnt");
  c.Add(Form("%s/stat*root",folder));
  cout << c.GetEntries() << " entries in trees, from location " << Form("%s/stat*",folder) << endl;
  c.SetBranchStatus("*",0);
  c.SetBranchStatus("best_fit_params",1);
  c.SetBranchStatus("input_params",1);
  PhysParamList * iparm = 0;
  FitParamList  * fparm = 0;
  c.SetBranchAddress("best_fit_params", &fparm);
  c.SetBranchAddress("input_params",    &iparm);

  double sd_mixing23 = 0;
  double sd_dm32sq = 0;
  TH1D * bfmix = new TH1D("","",100,0.2,0.8);
  for(int i=0; i<c.GetEntries(); i++) {
    c.GetEntry(i);
    double in_dm = iparm->Value(conv::kPhypDm23sq);
    double bf_dm = fparm->BestFitValue(conv::kPhypDm23sq);
    sd_dm32sq += (bf_dm - in_dm) * (bf_dm - in_dm);
    //need to be more careful with sinsqtheta23 - it has a mirror point
    double in_mix    = iparm->Value(conv::kPhypMixing23);
    double in13      = iparm->Value(conv::kPhypMixing13);
    double centre    = 1 / (2 * (1 - in13));
    double in_mirror = in_mix;
    if     (in_mix < centre)
      in_mirror = centre + TMath::Abs(centre - in_mix);
    else if(in_mix > centre)
      in_mirror = centre - TMath::Abs(centre - in_mix);
    double bf_mix = fparm->BestFitValue(conv::kPhypMixing23);
    double devsq_in     = (bf_mix - in_mix) * (bf_mix - in_mix);
    double devsq_mirror = (bf_mix - in_mirror) * (bf_mix - in_mirror);
    //cout << in_mix << "\t" << in_mirror << endl;
    sd_mixing23 += TMath::Min(devsq_in,devsq_mirror);
    //sd_mixing23 += devsq_in;
    bfmix->Fill(bf_mix);
  }
  sd_mixing23 = sqrt(sd_mixing23 / (c.GetEntries() - 1));
  sd_dm32sq = sqrt(sd_dm32sq / (c.GetEntries() - 1));

  cout << sd_dm32sq << endl;
  cout << sd_mixing23 << endl;
  cout << "lalalala\t" << bfmix->GetRMS() << "\t" << in_mix << endl;
  bfmix->GetXaxis()->SetRangeUser(0.2,0.51);
  cout << "lalala\t" << bfmix->GetRMS() << "\t" << in_mix << endl;

  ofstream stat;
  stat.open (Form("%s/stat.data",folder));
  stat << sd_mixing23 << endl;
  stat << sd_dm32sq << endl;
  stat.close();

  TCanvas * can = new TCanvas();
  bfmix->Draw();
  can->SaveAs(Form("%s/bf_%f.pdf",folder,in_mix));
}


/*
root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_-99_-99_-99_0.0022")'
root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_-99_-99_-99_-99")'
root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_-99_-99_-99_0.0026")'
root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.95_-99_-99_0.0022")'
root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.95_-99_-99_-99")'
root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.95_-99_-99_0.0026")'
root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.9_-99_-99_0.0022")'
root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.9_-99_-99_-99")'
root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.9_-99_-99_0.0026")'

root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.8_-99_-99_0.002")'
root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.8_-99_-99_0.003")'
root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.8_-99_-99_0.0024")'
root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_1.0_-99_-99_0.002")'
root2k -q '$T2KOSCANA/macros/syst_studies/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_1.0_-99_-99_0.003")'




cp nuosc3f_numudisapp_fix4/pulls_-99_-99_-99_-99_-99_0.0022/stat.data  bestfit_point_shift_1Rmu_-99_-99_-99_-99_-99_0.0022/stat.data
cp nuosc3f_numudisapp_fix4/pulls_-99_-99_-99_-99_-99_-99/stat.data     bestfit_point_shift_1Rmu_-99_-99_-99_-99_-99_-99/stat.data
cp nuosc3f_numudisapp_fix4/pulls_-99_-99_-99_-99_-99_0.0026/stat.data  bestfit_point_shift_1Rmu_-99_-99_-99_-99_-99_0.0026/stat.data
cp nuosc3f_numudisapp_fix4/pulls_-99_-99_0.95_-99_-99_0.0022/stat.data bestfit_point_shift_1Rmu_-99_-99_0.95_-99_-99_0.0022/stat.data
cp nuosc3f_numudisapp_fix4/pulls_-99_-99_0.95_-99_-99_-99/stat.data    bestfit_point_shift_1Rmu_-99_-99_0.95_-99_-99_-99/stat.data
cp nuosc3f_numudisapp_fix4/pulls_-99_-99_0.95_-99_-99_0.0026/stat.data bestfit_point_shift_1Rmu_-99_-99_0.95_-99_-99_0.0026/stat.data
cp nuosc3f_numudisapp_fix4/pulls_-99_-99_0.9_-99_-99_0.0022/stat.data  bestfit_point_shift_1Rmu_-99_-99_0.9_-99_-99_0.0022/stat.data
cp nuosc3f_numudisapp_fix4/pulls_-99_-99_0.9_-99_-99_-99/stat.data     bestfit_point_shift_1Rmu_-99_-99_0.9_-99_-99_-99/stat.data
cp nuosc3f_numudisapp_fix4/pulls_-99_-99_0.9_-99_-99_0.0026/stat.data  bestfit_point_shift_1Rmu_-99_-99_0.9_-99_-99_0.0026/stat.data
cp nuosc3f_numudisapp_fix4/pulls_-99_-99_0.8_-99_-99_0.002/stat.data   bestfit_point_shift_1Rmu_-99_-99_0.8_-99_-99_0.002/stat.data
cp nuosc3f_numudisapp_fix4/pulls_-99_-99_0.8_-99_-99_0.0024/stat.data  bestfit_point_shift_1Rmu_-99_-99_0.8_-99_-99_0.0024/stat.data
cp nuosc3f_numudisapp_fix4/pulls_-99_-99_0.8_-99_-99_0.003/stat.data   bestfit_point_shift_1Rmu_-99_-99_0.8_-99_-99_0.003/stat.data
cp nuosc3f_numudisapp_fix4/pulls_-99_-99_1.0_-99_-99_0.002/stat.data   bestfit_point_shift_1Rmu_-99_-99_1.0_-99_-99_0.002/stat.data
cp nuosc3f_numudisapp_fix4/pulls_-99_-99_1.0_-99_-99_0.003/stat.data   bestfit_point_shift_1Rmu_-99_-99_1.0_-99_-99_0.003/stat.data

 */
