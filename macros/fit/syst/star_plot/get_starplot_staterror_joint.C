void get_starplot_staterror_joint(const char * folder = "./") {

  TChain c("fit_tree");
  c.Add(Form("%s/pulls_theta_2333.root",folder));
  cout << c.GetEntries() << " entries in trees, from location " << Form("%s/pulls_theta_2333.root",folder) << endl;

  int fitstatus, migradfitstatus;
  c.SetBranchAddress("fitstatus",&fitstatus);
  c.SetBranchAddress("migradfitstatus",&migradfitstatus);

  double input_mixing23, bestfit_mixing23;
  c.SetBranchAddress("input_mixing23",&input_mixing23);
  c.SetBranchAddress("bestfit_mixing23",&bestfit_mixing23);
  double sd_mixing23 = 0;
  int entries = 0;
  for(int i=0; i<c.GetEntries(); i++) {
    c.GetEntry(i);
    if ((fitstatus != 3 ) || (migradfitstatus == 4)) continue;
    sd_mixing23 += (bestfit_mixing23 - input_mixing23) * (bestfit_mixing23 - input_mixing23);
    //cout << input_mixing23 << endl;
    entries ++ ;
  }
  sd_mixing23 = sqrt(sd_mixing23 / (entries - 1));

  double input_mixing13, bestfit_mixing13;
  c.SetBranchAddress("input_mixing13",&input_mixing13);
  c.SetBranchAddress("bestfit_mixing13",&bestfit_mixing13);
  double sd_mixing13 = 0;
  for(int i=0; i<c.GetEntries(); i++) {
    c.GetEntry(i);
    if ((fitstatus != 3 ) || (migradfitstatus == 4)) continue;
    sd_mixing13 += (bestfit_mixing13 - input_mixing13) * (bestfit_mixing13 - input_mixing13);
    //cout << input_mixing13 << endl;
  }
  sd_mixing13 = sqrt(sd_mixing13 / (entries - 1));

  double input_dm23sq, bestfit_dm23sq;
  c.SetBranchAddress("input_dm23sq",&input_dm23sq);
  c.SetBranchAddress("bestfit_dm23sq",&bestfit_dm23sq);
  double sd_dm23sq = 0;
  for(int i=0; i<c.GetEntries(); i++) {
    c.GetEntry(i);
    if ((fitstatus != 3 ) || (migradfitstatus == 4)) continue;
    sd_dm23sq += (bestfit_dm23sq - input_dm23sq) * (bestfit_dm23sq - input_dm23sq);
    //cout << input_dm23sq << "\t" << bestfit_dm23sq << endl;
  }
  sd_dm23sq = sqrt(sd_dm23sq / (entries - 1));

  double input_deltacp, bestfit_deltacp;
  c.SetBranchAddress("input_deltacp",&input_deltacp);
  c.SetBranchAddress("bestfit_deltacp",&bestfit_deltacp);
  double sd_deltacp = 0;
  for(int i=0; i<c.GetEntries(); i++) {
    c.GetEntry(i);
    if ((fitstatus != 3 ) || (migradfitstatus == 4)) continue;
    sd_deltacp += (bestfit_deltacp - input_deltacp) * (bestfit_deltacp - input_deltacp);
    //cout << input_deltacp << "\t" << bestfit_deltacp << endl;
  }
  sd_deltacp = sqrt(sd_deltacp / (entries - 1));


  cout << "entries = " << entries << endl;
  cout << sd_dm23sq << "\t" << input_dm23sq << endl;
  cout << sd_mixing23 << "\t" << input_mixing23 << endl;

  cout << sd_mixing13 << "\t" << input_mixing13 << endl;
  cout << sd_deltacp << "\t" << input_deltacp << endl;

  ofstream stat;
  stat.open (Form("%s/stat.data",folder));
  stat << sd_mixing23 << endl;
  stat << sd_dm23sq << endl;
  stat << sd_mixing13 << endl;
  stat << sd_deltacp << endl;
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
