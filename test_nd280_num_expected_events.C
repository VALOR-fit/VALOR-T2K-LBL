{
  TFile f("./inputs/official/nd280/saclay/neut/10d_v1/ND280NuEnergyHistoNewFV.root", "read");

  TH1F * h01 = (TH1F*)f.Get("EnergyNuMu_true_QE");
  TH1F * h02 = (TH1F*)f.Get("EnergyNuMu_true_res");
  TH1F * h03 = (TH1F*)f.Get("EnergyNuMu_true_oth");
  TH1F * h04 = (TH1F*)f.Get("EnergyNuMu_true_coh");
  TH1F * h05 = (TH1F*)f.Get("EnergyNuMu_true_ncel");
  TH1F * h06 = (TH1F*)f.Get("EnergyNuMu_true_ncres");
  TH1F * h07 = (TH1F*)f.Get("EnergyNuMu_true_ncoth");
  TH1F * h08 = (TH1F*)f.Get("EnergyNuMu_true_nccoh");
  TH1F * h09 = (TH1F*)f.Get("EnergyANuMu_true_QE");
  TH1F * h10 = (TH1F*)f.Get("EnergyANuMu_true_res");
  TH1F * h11 = (TH1F*)f.Get("EnergyANuMu_true_oth");
  TH1F * h12 = (TH1F*)f.Get("EnergyANuMu_true_coh");
  TH1F * h13 = (TH1F*)f.Get("EnergyANuMu_true_ncel");
  TH1F * h14 = (TH1F*)f.Get("EnergyANuMu_true_ncres");
  TH1F * h15 = (TH1F*)f.Get("EnergyANuMu_true_ncoth");
  TH1F * h16 = (TH1F*)f.Get("EnergyANuMu_true_nccoh");
  TH1F * h17 = (TH1F*)f.Get("EnergyNuE_true_QE");
  TH1F * h18 = (TH1F*)f.Get("EnergyNuE_true_res");
  TH1F * h19 = (TH1F*)f.Get("EnergyNuE_true_oth");
  TH1F * h20 = (TH1F*)f.Get("EnergyNuE_true_coh");
  TH1F * h21 = (TH1F*)f.Get("EnergyNuE_true_ncel");
  TH1F * h22 = (TH1F*)f.Get("EnergyNuE_true_ncres");
  TH1F * h23 = (TH1F*)f.Get("EnergyNuE_true_ncoth");
  TH1F * h24 = (TH1F*)f.Get("EnergyNuE_true_nccoh");

  double nd280_data_pot = 2.92E19;
  double nd280_mc_pot   = 49.85E19;
  double scale_factor   = nd280_data_pot/nd280_mc_pot;

  cout << "MC -> Data scaling factor = " << scale_factor << endl;

  h01 -> Scale (scale_factor);
  h02 -> Scale (scale_factor);
  h03 -> Scale (scale_factor);
  h04 -> Scale (scale_factor);
  h05 -> Scale (scale_factor);
  h06 -> Scale (scale_factor);
  h07 -> Scale (scale_factor);
  h08 -> Scale (scale_factor);
  h09 -> Scale (scale_factor);
  h10 -> Scale (scale_factor);
  h11 -> Scale (scale_factor);
  h12 -> Scale (scale_factor);
  h13 -> Scale (scale_factor);
  h14 -> Scale (scale_factor);
  h15 -> Scale (scale_factor);
  h16 -> Scale (scale_factor);
  h17 -> Scale (scale_factor);
  h18 -> Scale (scale_factor);
  h19 -> Scale (scale_factor);
  h20 -> Scale (scale_factor);
  h21 -> Scale (scale_factor);
  h22 -> Scale (scale_factor);
  h23 -> Scale (scale_factor);
  h24 -> Scale (scale_factor);

  double nexp = 
   	  h01 -> Integral () +
	  h02 -> Integral () +
	  h03 -> Integral () +
	  h04 -> Integral () +
	  h05 -> Integral () +
	  h06 -> Integral () +
	  h07 -> Integral () +
	  h08 -> Integral () +
	  h09 -> Integral () +
	  h10 -> Integral () +
	  h11 -> Integral () +
	  h12 -> Integral () +
	  h13 -> Integral () +
	  h14 -> Integral () +
	  h15 -> Integral () +
	  h16 -> Integral () +
	  h17 -> Integral () +
	  h18 -> Integral () +
	  h19 -> Integral () +
	  h20 -> Integral () +
	  h21 -> Integral () +
	  h22 -> Integral () +
	  h23 -> Integral () +
	  h24 -> Integral ();

  cout << "N_{exp} " << nexp << endl;
}
