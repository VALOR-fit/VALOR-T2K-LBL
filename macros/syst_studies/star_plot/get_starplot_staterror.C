void get_starplot_staterror(
  const char * folder    = "./",
  const char * par1      = "mixing23",
  const char * par2      = "dmsq23"
)
{
  using namespace valor;
  using namespace valor::conv;
  Opt::Instance()->gOscMixingParamChoice == conv::kNuMixSinsqTheta;
  
  PhysParamId_t p1 = util::PhysParamIdFromString(par1);
  PhysParamId_t p2 = util::PhysParamIdFromString(par2);
  
  assert(p1!=kPhypUndefined);
  assert(p2!=kPhypUndefined);
  
  string pattern = Form("%s/star_plot_stat*root",folder);
  TChain c("fitnt");
  c.Add(pattern.c_str());
  cout << c.GetEntries() << " entries in trees, from location " << pattern << endl;
  c.SetBranchStatus("*",0);
  c.SetBranchStatus("best_fit_params",1);
  c.SetBranchStatus("input_params",1);
  PhysParamList * iparm = 0;
  FitParamList  * fparm = 0;
  c.SetBranchAddress("best_fit_params", &fparm);
  c.SetBranchAddress("input_params",    &iparm);
  
  c.GetEntry(0);

  double sd_p1 = 0;
  double sd_p2 = 0;
  
  assert(fparm->Exists(p1));
  assert(fparm->Exists(p2));
  
  TH1D * h_bf_p1 = new TH1D("","",100,fparm->MinValue(p1),fparm->MaxValue(p1));
  TH1D * h_bf_p2 = new TH1D("","",100,fparm->MinValue(p2),fparm->MaxValue(p2));
  
  for(int i=0; i<c.GetEntries(); i++) {
    c.GetEntry(i);
    
    double in_p1=-999,in_p2=-999,bf_p1=-999,bf_p2=-999,devsq_in_p1=-999,devsq_in_p2=-999;
    double in_mirror_p1    = 1E100;
    double in_mirror_p2    = 1E100;
    double devsq_mirror_p1 = 1E100;
    double devsq_mirror_p2 = 1E100;
    
    in_p1 = iparm->Value(p1);
    bf_p1 = fparm->BestFitValue(p1);
    in_p2 = iparm->Value(p2);
    bf_p2 = fparm->BestFitValue(p2);
    
    if (p1 == conv::kPhypMixing23 || p1 == conv::kPhypMixing23Bar) { // handle the mirror point if necessary
      double in13 = iparm->Value(conv::kPhypMixing13);
      double centre    = 1 / (2 * (1 - in13));
      double in_p1     = iparm->Value(p1);
      if     (in_p1 < centre) { in_mirror_p1 = centre + TMath::Abs(centre - in_p1);}
      else if(in_p1 > centre) { in_mirror_p1 = centre - TMath::Abs(centre - in_p1);}
      else {in_mirror_p1 = in_p1;}
      
      devsq_in_p1     = (bf_p1 - in_p1       ) * (bf_p1 - in_p1       );
      devsq_mirror_p1 = (bf_p1 - in_mirror_p1) * (bf_p1 - in_mirror_p1);
    }
    else {
      devsq_in_p1 = pow((bf_p2 - in_p2),2);
    }
    
    if (p2 == conv::kPhypMixing23 || p2 == conv::kPhypMixing23Bar) { // handle the mirror point if necessary
      double in13 = iparm->Value(conv::kPhypMixing13);
      double centre    = 1 / (2 * (1 - in13));
      double in_p2     = iparm->Value(p2);
      if     (in_p2 < centre) { in_mirror_p2 = centre + TMath::Abs(centre - in_p2);}
      else if(in_p2 > centre) { in_mirror_p2 = centre - TMath::Abs(centre - in_p2);}
      else {in_mirror_p2 = in_p2;}
      
      devsq_in_p2     = (bf_p2 - in_p2       ) * (bf_p2 - in_p2       );
      devsq_mirror_p2 = (bf_p2 - in_mirror_p2) * (bf_p2 - in_mirror_p2);
    }
    else {
      devsq_in_p2 = pow((bf_p2 - in_p2),2);
    }
    
    //cout << in_p1 << "\t" << in_mirror << endl;
    sd_p1 += TMath::Min(devsq_in_p1,devsq_mirror_p1);
    sd_p2 += TMath::Min(devsq_in_p2,devsq_mirror_p2);
    //sd_p1 += devsq_in;
    h_bf_p1->Fill(bf_p1);
    h_bf_p2->Fill(bf_p2);
  }
  sd_p1 = sqrt(sd_p1 / (c.GetEntries() - 1));
  sd_p2 = sqrt(sd_p2 / (c.GetEntries() - 1));
  
  cout << sd_p1 << endl;
  cout << sd_p2 << endl;
  cout << "P1 RMS:\t" << h_bf_p1->GetRMS() << "\t" << in_p1 << "\t"<<util::AsString(p1)<<endl;
  cout << "P2 RMS:\t" << h_bf_p2->GetRMS() << "\t" << in_p2 << "\t"<<util::AsString(p2)<<endl;

  ofstream stat;
  stat.open (Form("%s/stat.data",folder));
  stat << sd_p1 << endl;
  stat << sd_p2 << endl;
  stat.close();

  TCanvas * can_p1 = new TCanvas();
  h_bf_p1->Draw();
  can_p1->SaveAs(Form("%s/bf_%s_%lf.pdf",folder,util::AsStringSimple(p1).c_str(),in_p1));
  TCanvas * can_p2 = new TCanvas();
  h_bf_p2->Draw();
  can_p2->SaveAs(Form("%s/bf_%s_%lf.pdf",folder,util::AsStringSimple(p2).c_str(),in_p2));
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
