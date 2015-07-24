void weight_selection_trees()

{
  using namespace std;
  using namespace valor;

  // Global options                                                                                                                                                                                     
  Opt::Instance()->gOscMixingParamChoice = conv::kNuMixSinsq2Theta;
  Opt::Instance()->gMassHierarchyAssumed = conv::kNmhNormal;

  //files
  const char * filename1 = "inputs/official/2012/superk/mc/neut/11c/superk_1ReT_numu.root";
  TFile * f1 = new TFile(filename1,"read");
  const char * filename2 = "inputs/official/2012/superk/mc/neut/11c/superk_1ReT_numubar.root";
  TFile * f2 = new TFile(filename2,"read");
  const char * filename3 = "inputs/official/2012/superk/mc/neut/11c/superk_1ReT_nue.root";
  TFile * f3 = new TFile(filename3,"read");
  const char * filename4 = "inputs/official/2012/superk/mc/neut/11c/superk_1ReT_nuebar.root";
  TFile * f4 = new TFile(filename4,"read");
  const char * filename5 = "inputs/official/2012/superk/mc/neut/11c/superk_1ReT_nuesig.root";
  TFile * f5 = new TFile(filename5,"read");



  //config file
  const char * config    = "analysis_configuration.xml";
  const char * config_file  = Form("%s/config/%s", gSystem->Getenv("T2KOSCANA"), config);

  // 
  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read(config_file);

  //read tree
  string treename1 = "selection_tree";
  TTree * t1 = (TTree*) f1->Get(treename1.c_str());
  TTree * t2 = (TTree*) f2->Get(treename1.c_str());
  TTree * t3 = (TTree*) f3->Get(treename1.c_str());
  TTree * t4 = (TTree*) f4->Get(treename1.c_str());
  TTree * t5 = (TTree*) f5->Get(treename1.c_str());

  int mode1, mode2, mode3, mode4, mode5; 
  int ipnu1, ipnu2, ipnu3, ipnu4, ipnu5;
  double erec1, erec2, erec3, erec4, erec5;
  double pnu1, pnu2, pnu3, pnu4, pnu5;

  t1->SetBranchAddress("Mode",&mode1);
  t1->SetBranchAddress("Pdg",&ipnu1);
  t1->SetBranchAddress("EnuReco",&erec1);
  t1->SetBranchAddress("Enu",&pnu1);

  t2->SetBranchAddress("Mode",&mode2);
  t2->SetBranchAddress("Pdg",&ipnu2);
  t2->SetBranchAddress("EnuReco",&erec2);
  t2->SetBranchAddress("Enu",&pnu2);

  t3->SetBranchAddress("Mode",&mode3);
  t3->SetBranchAddress("Pdg",&ipnu3);
  t3->SetBranchAddress("EnuReco",&erec3);
  t3->SetBranchAddress("Enu",&pnu3);

  t4->SetBranchAddress("Mode",&mode4);
  t4->SetBranchAddress("Pdg",&ipnu4);
  t4->SetBranchAddress("EnuReco",&erec4);
  t4->SetBranchAddress("Enu",&pnu4);

  t5->SetBranchAddress("Mode",&mode5);
  t5->SetBranchAddress("Pdg",&ipnu5);
  t5->SetBranchAddress("EnuReco",&erec5);
  t5->SetBranchAddress("Enu",&pnu5);

  const char * nameout = "sk11c_1ReT_flx11bv32fine_10acuts.root";
  TFile * fout = new TFile(nameout,"recreate");

  TTree * tout1 = new TTree("tnumu","Summary information from superk ntuple");
  TTree * tout2 = new TTree("tnumubar","Summary information from superk ntuple");
  TTree * tout3 = new TTree("tnue","Summary information from superk ntuple");
  TTree * tout4 = new TTree("tnuebar","Summary information from superk ntuple");
  TTree * tout5 = new TTree("tnuesig","Summary information from superk ntuple");

  int ipnu_out1, ipnu_out2, ipnu_out3, ipnu_out4, ipnu_out5;
  float pnu_out1, pnu_out2, pnu_out3, pnu_out4, pnu_out5;
  int mode_out1, mode_out2, mode_out3, mode_out4, mode_out5;
  float erec_out1, erec_out2, erec_out3, erec_out4, erec_out5;
  double wgtflx1, wgtflx2, wgtflx3, wgtflx4, wgtflx5;

  tout1->Branch("mode", &mode_out1, "mode/I");
  tout1->Branch("enu", &pnu_out1, "pnu/F");
  tout1->Branch("erec", &erec_out1, "erec/F");
  tout1->Branch("wgt", &wgtflx1, "wgt/D");

  tout2->Branch("mode", &mode_out2, "mode/I");
  tout2->Branch("enu", &pnu_out2, "pnu/F");
  tout2->Branch("erec", &erec_out2, "erec/F");
  tout2->Branch("wgt", &wgtflx2, "wgt/D");

  tout3->Branch("mode", &mode_out3, "mode/I");
  tout3->Branch("enu", &pnu_out3, "pnu/F");
  tout3->Branch("erec", &erec_out3, "erec/F");
  tout3->Branch("wgt", &wgtflx3, "wgt/D");

  tout4->Branch("mode", &mode_out4, "mode/I");
  tout4->Branch("enu", &pnu_out4, "pnu/F");
  tout4->Branch("erec", &erec_out4, "erec/F");
  tout4->Branch("wgt", &wgtflx4, "wgt/D");

  tout5->Branch("mode", &mode_out5, "mode/I");
  tout5->Branch("enu", &pnu_out5, "pnu/F");
  tout5->Branch("erec", &erec_out5, "erec/F");
  tout5->Branch("wgt", &wgtflx5, "wgt/D");



  //calculate weights
  for(int i=0; i<t1->GetEntries(); i++) {

    t1->GetEntry(i);
    int i_ipnu = ipnu1;
    double i_pnu = pnu1;
    wgtflx1 = skexp->GetNuFlux()->SuperkFluxExtTuningWght(i_ipnu,i_pnu);

    pnu_out1 = (float)pnu1;
    ipnu_out1 = ipnu1;
    erec_out1 = (float) erec1*1000.;
    mode_out1 = (int) mode1;

    tout1->Fill();
  }

  for(int i=0; i<t2->GetEntries(); i++) {

    t2->GetEntry(i);
    int i_ipnu = ipnu2;
    double i_pnu = pnu2;
    wgtflx2 = skexp->GetNuFlux()->SuperkFluxExtTuningWght(i_ipnu,i_pnu);

    pnu_out2 = (float)pnu2;
    ipnu_out2 = ipnu2;
    erec_out2 = (float) erec2*1000.;
    mode_out2 = (int) mode2;

    tout2->Fill();
  }

  for(int i=0; i<t3->GetEntries(); i++) {

    t3->GetEntry(i);
    int i_ipnu = ipnu3;
    double i_pnu = pnu3;
    wgtflx3 = skexp->GetNuFlux()->SuperkFluxExtTuningWght(i_ipnu,i_pnu);

    pnu_out3 = (float)pnu3;
    ipnu_out3 = ipnu3;
    erec_out3 = (float) erec3*1000.;
    mode_out3 = (int) mode3;

    tout3->Fill();
  }

  for(int i=0; i<t4->GetEntries(); i++) {

    t4->GetEntry(i);
    int i_ipnu = ipnu4;
    double i_pnu = pnu4;
    wgtflx4 = skexp->GetNuFlux()->SuperkFluxExtTuningWght(i_ipnu,i_pnu);

    pnu_out4 = (float)pnu4;
    ipnu_out4 = ipnu4;
    erec_out4 = (float) erec4*1000.;
    mode_out4 = (int) mode4;

    tout4->Fill();
  }

  for(int i=0; i<t5->GetEntries(); i++) {

    t5->GetEntry(i);
    int i_ipnu = ipnu5;
    double i_pnu = pnu5;
    //for nue signal, use numu flux
    wgtflx5 = skexp->GetNuFlux()->SuperkFluxExtTuningWght(conv::kPdgNuMu,i_pnu);

    pnu_out5 = (float)pnu5;
    ipnu_out5 = ipnu5;
    erec_out5 = (float) erec5*1000.;
    mode_out5 = (int) mode5;

    tout5->Fill();
  }

  //  cout << "filename2=" << filename2 << endl;
  fout->cd();
  tout1->Write();
  tout2->Write();
  tout3->Write();
  tout4->Write();
  tout5->Write();

}
