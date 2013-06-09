/*
Takes the output of the star plot (in tree form), finds the maximum deviation from nominal (for a given parameter,
 as a fraction of the stat error, from a variety of points)

Tom
 */

void max_deviation_table(const char * filename = "readfile.root")
{
  using namespace std;
  using namespace t2kana;

  /*TFile f(filename);
  TTree * tree = 0;
  f.GetObject("tree",tree);*/
  TChain tree("tree");
  cout << tree.Add(filename) << " files" << endl;
  cout << tree.GetEntries() << " entries" << endl;
  double maxmix = 0;
  double maxdm = 0;

  double var_mix23P = 0, var_mix23M = 0, var_inmix23 = 0, var_stmix23 = 0;
  double var_dm23P = 0, var_dm23M = 0, var_indm23 = 0, var_stdm23 = 0;
  char var_name[100];
  tree.SetBranchAddress("mix23P",&var_mix23P);
  tree.SetBranchAddress("mix23M",&var_mix23M);
  tree.SetBranchAddress("inmix23",&var_inmix23);
  tree.SetBranchAddress("stmix23",&var_stmix23);
  tree.SetBranchAddress("dm23P",&var_dm23P);
  tree.SetBranchAddress("dm23M",&var_dm23M);
  tree.SetBranchAddress("indm23",&var_indm23);
  tree.SetBranchAddress("stdm23",&var_stdm23);
  tree.SetBranchAddress("name",var_name);

  cout << fixed;
  cout.precision(4);

  vector<string> systematics;
  systematics.push_back("sk_energy_scale");
  /*
  systematics.push_back("sk_1Rmu_ccqe_other_eff");
  systematics.push_back("sk_1Rmu_ccqe_ring_counting_eff_0");
  systematics.push_back("sk_1Rmu_ccqe_ring_counting_eff_1");
  systematics.push_back("sk_1Rmu_ccqe_ring_counting_eff_2");
  systematics.push_back("sk_1Rmu_ccnonqe_eff");
  systematics.push_back("sk_1Rmu_nc_eff");
  systematics.push_back("sk_1Rmu_nuecc_eff");
  systematics.push_back("fsi_000");
  systematics.push_back("fsi_001");
  systematics.push_back("fsi_002");
  systematics.push_back("fsi_003");
  systematics.push_back("fsi_004");
  systematics.push_back("fsi_005");
  */
  systematics.push_back("skdetfsi_000");
  systematics.push_back("skdetfsi_001");
  systematics.push_back("skdetfsi_002");
  systematics.push_back("skdetfsi_003");
  systematics.push_back("skdetfsi_004");
  systematics.push_back("skdetfsi_005");
  systematics.push_back("banff_00");
  systematics.push_back("banff_01");
  systematics.push_back("banff_02");
  systematics.push_back("banff_03");
  systematics.push_back("banff_04");
  systematics.push_back("banff_05");
  systematics.push_back("banff_06");
  systematics.push_back("banff_07");
  systematics.push_back("banff_08");
  systematics.push_back("banff_09");
  systematics.push_back("banff_10");
  systematics.push_back("banff_11");
  systematics.push_back("banff_12");
  systematics.push_back("banff_13");
  systematics.push_back("banff_14");
  systematics.push_back("banff_15");
  systematics.push_back("banff_maqe");
  systematics.push_back("banff_mares");
  systematics.push_back("banff_ccqeE0");
  systematics.push_back("banff_ccqeE1");
  systematics.push_back("banff_ccqeE2");
  systematics.push_back("banff_cc1piE0");
  systematics.push_back("banff_cc1piE1");
  systematics.push_back("xsec_pf");
  systematics.push_back("xsec_bindE");
  systematics.push_back("xsec_sf");
  systematics.push_back("xsec_pilessdelta");
  systematics.push_back("xsec_shape_ccoth");
  systematics.push_back("xsec_wshape");
  systematics.push_back("xsec_norm_cccoh");
  systematics.push_back("xsec_norm_nc1piC");
  systematics.push_back("xsec_norm_ncoth");
  systematics.push_back("xsec_norm_nuetonumu");
  systematics.push_back("xsec_norm_antinutonu");


  for(int isys = 0; isys < systematics.size(); isys++) {
    string name = systematics[isys];

    TString cutstr = Form("name==\"%s\"", name.c_str());
    const char * cut = cutstr.Data();

    maxmix = 0; maxdm = 0;

    for(int ie = 0; ie < tree.GetEntries(); ie++) {
      tree.GetEntry(ie);
      if(var_name != name) continue;
      
      maxmix = TMath::Max(maxmix, TMath::Max( TMath::Abs((var_mix23P-var_inmix23)/var_stmix23), TMath::Abs((var_mix23M-var_inmix23)/var_stmix23) ));
      maxdm = TMath::Max(maxdm, TMath::Max( TMath::Abs((var_dm23P-var_indm23)/var_stdm23), TMath::Abs((var_dm23M-var_indm23)/var_stdm23) ));

    }
    if(maxmix > 0.25 || maxdm > 0.25)
      cout << "{\\color{red}";

    cout << util::AsStringPaperTex(util::PhysParamIdFromString(name));
    if(maxmix > 0.25 || maxdm > 0.25)
      cout << "}";

    cout << "\t & " << maxmix << "\t & " << maxdm << " \\\\" << endl;

  }

}
