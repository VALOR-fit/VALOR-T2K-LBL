//
// ------------------------------------------------------------------------------------------------
// Test histograms provided by the ND280 / Saclay group
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

void test_nd280_histograms
(
   string generator = "neut",     /*  "neut" or "genie" */
   int    flux10dvrs = 3          /*   2 --> 10d_v2, 3 --> 10d_v3 */
)

{
  double nd280_data_pot     = 2.88E19;
  double nd280_mc_pot       = 0;
  double nd280_mc_pot_neut  = 4.98 * (1.025*0.982) * 1E20; // NEUT
  double nd280_mc_pot_genie = 4.94  *(1.025*0.982) * 1E20; // GENIE 

  string filename = "";

  if(generator == "neut" && flux10dvrs == 2) {
     filename = "inputs/official/2010/nd280/saclay/10d_v2/ND280NuEnergyHisto_NEUT.root";
     nd280_mc_pot = nd280_mc_pot_neut;
  }
  if(generator == "neut" && flux10dvrs == 3) {
     filename = "inputs/official/2010/nd280/saclay/10d_v3/ND280NuEnergyHisto_NEUT10d_v3.root";
     nd280_mc_pot = nd280_mc_pot_neut;
  }
  if(generator == "genie" && flux10dvrs == 2) {
     filename = "inputs/official/2010/nd280/saclay/10d_v2/ND280NuEnergyHisto_GENIE.root";
     nd280_mc_pot = nd280_mc_pot_genie;
  }
  if(generator == "genie" && flux10dvrs == 3) {
     filename = "inputs/official/2010/nd280/saclay/10d_v3/ND280NuEnergyHisto_GENIE10d_v3.root";
     nd280_mc_pot = nd280_mc_pot_genie;
  }

  if(filename == "") {
     cout << "Can't find nd280 histograms" << endl;
     return;
  }

  string topdir   = gSystem->Getenv("T2KOSCANA");
  filename = topdir + "/" + filename;

  TFile f(filename.c_str(),  "read");

  TH1F * h01 = (TH1F*)f.Get("EnergyNuMu_true_QE");
  TH1F * h02 = (TH1F*)f.Get("EnergyNuMu_true_res");
  TH1F * h03 = (TH1F*)f.Get("EnergyNuMu_true_oth");
  TH1F * h04 = (TH1F*)f.Get("EnergyNuMu_true_coh");
  TH1F * h05 = (TH1F*)f.Get("EnergyNuMu_true_ncres");
  TH1F * h06 = (TH1F*)f.Get("EnergyNuMu_true_ncoth");
  TH1F * h07 = (TH1F*)f.Get("EnergyNuMu_true_nccoh");
  TH1F * h08 = (TH1F*)f.Get("EnergyANuMu_true_QE");
  TH1F * h09 = (TH1F*)f.Get("EnergyANuMu_true_res");
  TH1F * h10 = (TH1F*)f.Get("EnergyANuMu_true_oth");
  TH1F * h11 = (TH1F*)f.Get("EnergyANuMu_true_coh");
  TH1F * h12 = (TH1F*)f.Get("EnergyANuMu_true_ncres");
  TH1F * h13 = (TH1F*)f.Get("EnergyANuMu_true_ncoth");
  TH1F * h14 = (TH1F*)f.Get("EnergyANuMu_true_nccoh");
  TH1F * h15 = (TH1F*)f.Get("EnergyNuE_true_QE");
  TH1F * h16 = (TH1F*)f.Get("EnergyNuE_true_res");
  TH1F * h17 = (TH1F*)f.Get("EnergyNuE_true_oth");
  TH1F * h18 = (TH1F*)f.Get("EnergyNuE_true_coh");
  TH1F * h19 = (TH1F*)f.Get("EnergyNuE_true_ncres");
  TH1F * h20 = (TH1F*)f.Get("EnergyNuE_true_ncoth");
  TH1F * h21 = (TH1F*)f.Get("EnergyNuE_true_nccoh");

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
	  h21 -> Integral ();

  cout << "N_{exp} = " << nexp << endl;

  double nexp_qe = 
   	  h01 -> Integral () +
   	  h08 -> Integral () +
   	  h15 -> Integral ();

  cout << "N_{exp;QE} = " << nexp_qe << endl;

}
