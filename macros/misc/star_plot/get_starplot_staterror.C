void get_starplot_staterror(const char * folder = "./") {

  TChain c("nfitnt");
  c.Add(Form("%s/starstat-run_*",folder));
  cout << c.GetEntries() << " entries in trees, from location " << Form("%s/starstat-run_*",folder) << endl;

  double input_mixing23, bestfit_mixing23;
  c.SetBranchAddress("input_mixing23",&input_mixing23);
  c.SetBranchAddress("bestfit_mixing23",&bestfit_mixing23);
  double sd_mixing23 = 0;
  for(int i=0; i<c.GetEntries(); i++) {
    c.GetEntry(i);
    sd_mixing23 += (bestfit_mixing23 - input_mixing23) * (bestfit_mixing23 - input_mixing23);
    //cout << input_mixing23 << endl;
  }
  sd_mixing23 = sqrt(sd_mixing23 / (c.GetEntries() - 1));

  double input_dm23sq, bestfit_dm23sq;
  c.SetBranchAddress("input_dm23sq",&input_dm23sq);
  c.SetBranchAddress("bestfit_dm23sq",&bestfit_dm23sq);
  double sd_dm23sq = 0;
  for(int i=0; i<c.GetEntries(); i++) {
    c.GetEntry(i);
    sd_dm23sq += (bestfit_dm23sq - input_dm23sq) * (bestfit_dm23sq - input_dm23sq);
    //cout << input_dm23sq << "\t" << bestfit_dm23sq << endl;
  }
  sd_dm23sq = sqrt(sd_dm23sq / (c.GetEntries() - 1));
  cout << sd_dm23sq << "\t" << input_dm23sq << endl;
  cout << sd_mixing23 << "\t" << input_mixing23 << endl;

  ofstream stat;
  stat.open (Form("%s/stat.data",folder));
  stat << sd_mixing23 << endl;
  stat << sd_dm23sq << endl;
  stat.close();
}


/*
root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_-99_-99_-99_0.0022")'
root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_-99_-99_-99_-99")'
root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_-99_-99_-99_0.0026")'
root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.95_-99_-99_0.0022")'
root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.95_-99_-99_-99")'
root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.95_-99_-99_0.0026")'
root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.9_-99_-99_0.0022")'
root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.9_-99_-99_-99")'
root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.9_-99_-99_0.0026")'

root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.8_-99_-99_0.002")'
root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.8_-99_-99_0.003")'
root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_0.8_-99_-99_0.0024")'
root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_1.0_-99_-99_0.002")'
root -q '$T2KOSCANA/macros/misc/star_plot/get_starplot_staterror.C("nuosc3f_numudisapp_fix4/pulls_-99_-99_1.0_-99_-99_0.003")'




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
