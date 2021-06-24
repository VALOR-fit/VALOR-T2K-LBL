//
// Convert the data minituple from the SK group into a VALOR data DST format to be read by GSSFD::FillFromDST
//
// NOTE: Be aware that the sample, beam and reac mode for each event is stored as an int (in the output tree),
// so if any of these enumerations change, the MC file will have to be regenerated
//
// Francis Bench <F.Bench \at liverpool.ac.uk>
// November 2019


#include <vector>
#include <utility>


//___________________________________________________________________________
// Hard-coded info

// T2K beam direction
double dir_nu[3] = {0.669764, -0.742179, 0.024223};

// SK selection codes
int sk_code_1re     = 11;
int sk_code_1rmu    = 13;
int sk_code_1re1dcy = 14;

const char* valor_minidst_name = "valor_t2k_data_mini_dst";

bool print_low_level_info = false;

//___________________________________________________________________________
// Prototypes
bool   EventIsMuLike   (int iclass);
bool   EventIsELike    (int iclass);
bool   EventIsFHC      (int horn_current);
double CalcThetaDegrees(float dir[][7][3], int iclass);


//___________________________________________________________________________
void convert_t2k_data_minituple_into_root_fitqun_valor_dst
(
 const char* output_tag     = "SK_data_Run1-9d_valor_minidst",
 const char* input_file     = "$VALOR_EXPT/data/official/2018/data/t2ksk.14c.fqv4r0.full_run1-9.data.minituple_cuts_applied.root",
 const char* input_treename = "minituple"
 )
{
  using namespace valor;
  using namespace valor::conv;

  // At this point, we don't have an MC, so can't load in the config files to read the beam and sample lists, so hard-code them here
  std::vector<BeamConf_t> beams = {kbJPARC_FHC, kbJPARC_RHC};
  std::map<BeamConf_t, std::vector<Sample_t>> samples;
  std::vector<Sample_t> vs = {ks1Rmu, ks1Re, ksMultiRe};
  samples.insert(std::make_pair(kbJPARC_FHC, vs));
  vs = {ks1Rmu, ks1Re};
  samples.insert(std::make_pair(kbJPARC_RHC, vs));

  // Read in the input SK tree
  TFile * fin = new TFile(input_file,"OPEN");
  assert(fin);
  TTree * tree = (TTree*) fin->Get(input_treename);
  assert(tree);

  // Set up the input tree branches
  double erec;          // Neutrino reconstructed energy
  int    iclass;        // SK selection code. 11 is single-ring e-like, 13 is single-ring mu-like, 14 is 1Re1dcy (i.e. e CC1pi-like)
  float  dir[2][7][3];  // Array containing the direction of the outgoing lepton
  double horn_current;  // FHC (>0) or RHC (<0)
  tree->SetBranchAddress("erec"        , &erec);
  tree->SetBranchAddress("iclass"      , &iclass);
  tree->SetBranchAddress("horn_current", &horn_current);
  tree->SetBranchAddress("fq1rdir"     , &dir);

  // Set up the output file and VALOR miniDST tree
  TString out_file_name = TString::Format("%s.root", output_tag);
  TFile* fout = new TFile(out_file_name, "RECREATE");

  int    valor_minidst_beam       = -1;  // see Conventions/BeamConf.h
  int    valor_minidst_sample     = -1;  // see Conventions/Sample.h
  double valor_minidst_enu_reco   = -1.; // reco neutrino energy, in GeV
  double valor_minidst_theta_reco = -1.; // reco angle between the outgoing lepton and the beam direction
  double valor_minidst_wght       = -1.; // misc weight

  TTree * valor_minidst = new TTree(valor_minidst_name, "Mini-DST containing VALOR/T2K data");

  valor_minidst->Branch("beam"      , &valor_minidst_beam        ,"beam/I"      );
  valor_minidst->Branch("sample"    , &valor_minidst_sample      ,"sample/I"    );
  valor_minidst->Branch("enu_reco"  , &valor_minidst_enu_reco    ,"enu_reco/D"  );
  valor_minidst->Branch("Theta_reco", &valor_minidst_theta_reco  ,"Theta_reco/D");
  valor_minidst->Branch("wght"      , &valor_minidst_wght        ,"wght/D"      );

  // Also set up a vector for each sample, containing the (Erec, theta) values for each event
  std::map<conv::BeamConf_t, std::map<conv::Sample_t, std::vector<std::map<conv::KineVar_t, double>>>> data_kinematics;
  for(conv::BeamConf_t b : beams) {
    for(conv::Sample_t s : samples.at(b)) {
      data_kinematics[b][s] = std::vector<std::map<conv::KineVar_t, double>>();
    }
  }

  // Loop over the input tree to determine the kinematics, sample, beam and weights of each event and fill the output tree
  int n_fhc_mu=0,n_rhc_mu=0,n_fhc_e=0,n_rhc_e=0,n_multie=0;
  int nentries = tree->GetEntriesFast();
  int nev_survived_cuts = 0;
  for(int i = 0; i < nentries; i++) {
    tree->GetEntry(i);

    // Apply unit conversions
    erec /= 1000.;  // In the input tree, Erec is in Mev, but VALOR uses GeV

    // Calculate theta
    double costheta = 0.0;
    if(EventIsELike(iclass)) {
       costheta = dir_nu[0] * dir[0][1][0] + dir_nu[1] * dir[0][1][1] + dir_nu[2] * dir[0][1][2];
    }
    else if(EventIsMuLike(iclass)) {
       costheta = dir_nu[0] * dir[0][2][0] + dir_nu[1] * dir[0][2][1] + dir_nu[2] * dir[0][2][2];
    }
    double LepTheta = TMath::ACos(costheta);          // radians
    double theta    = LepTheta / TMath::Pi() * 180.;  // degrees

    // Apply cuts
    if(erec < 0.) { continue; }
    if(EventIsELike(iclass) && erec > 1.25) { continue; }
    if(theta < 0.) { continue; }
    nev_survived_cuts++;

    // Determine the VALOR SK sample enum and count the total number of events in each sample
    int selection = -1;
    if(iclass == sk_code_1re) {
      selection = (int) conv::ks1Re;
      if(EventIsFHC(horn_current)){ n_fhc_e++; }
      else                        { n_rhc_e++; }
    }
    else if(iclass == sk_code_1re1dcy) {
      selection = (int) conv::ks1Re1dcy;
      if(EventIsFHC(horn_current)){ n_multie++; }
      else {
        std::cerr << "ERROR: An RHC 1Re1dcy event. This shouldn't happen!" << std::endl;
        exit(1);
      }
    }
    else if(iclass == sk_code_1rmu) {
      selection = (int) conv::ks1Rmu;
      if(EventIsFHC(horn_current)){ n_fhc_mu++; }
      else                        { n_rhc_mu++; }
    }
    if(selection < 0) {
      std::cerr << "ERROR: Selection for event #" << i << " was undefined." << std::endl;
      exit(1);
    }

    // Determine the VALOR beam enum
    int beam = -1;
    if(EventIsFHC(horn_current)) { beam = (int) conv::kbJPARC_FHC; }
    else                         { beam = (int) conv::kbJPARC_RHC; }
    if(beam < 0) {
      std::cerr << "ERROR: Beam for event #" << i << " was undefined." << std::endl;
      exit(1);
    }

    // Fill the vector of kinematic values
    std::map<conv::KineVar_t, double> kinematic_map;
    kinematic_map[conv::kvEreco]     = erec;
    kinematic_map[conv::kvThetareco] = theta;
    data_kinematics.at((conv::BeamConf_t)beam).at((conv::Sample_t)selection).push_back(kinematic_map);

    // Save event in VALOR format to the output VALOR mini-DST
    valor_minidst_beam       = beam;
    valor_minidst_sample     = selection;
    valor_minidst_enu_reco   = erec;
    valor_minidst_theta_reco = theta;
    valor_minidst_wght       = 1.; // The data has no weights, so set the weight branch to 1
    valor_minidst->Fill();

    if(i < 10 || print_low_level_info){
      std::cout << "\nInformation copied to VALOR T2K data Mini-DST: "
                << "\n Beam          = " << util::AsString(util::BeamConfFromInt(valor_minidst_beam))
                << "\n Sample        = " << util::AsString(util::SampleFromInt(valor_minidst_sample))
                << "\n Ev (reco)     = " << valor_minidst_enu_reco << " GeV"
                << "\n theta (reco)  = " << valor_minidst_theta_reco
                << "\n misc weight   = " << valor_minidst_wght
                << std::endl;
    }
  } // i event in SK tree

  // Print # of events
  std::cout << "\n\nEvent rates:" << endl;
  std::cout << " FHC 1Rmu: " << n_fhc_mu << std::endl;
  std::cout << " FHC 1Re:  " << n_fhc_e  << std::endl;
  std::cout << " RHC 1Rmu: " << n_rhc_mu << std::endl;
  std::cout << " RHC 1Re:  " << n_rhc_e  << std::endl;
  std::cout << " 1Re1dcy:  " << n_multie << std::endl;
  
  fin->Close();
  
  // Write the output tree
  std::cout << "Writing the output VALOR data mini-DST to the file \'" << out_file_name.Data() << "\'" << std::endl;
  valor_minidst->Write(valor_minidst_name);

  // Write a 2D TGraph for each sample with the kinematics of each event
  conv::Detector_t d = conv::kdSuperK;
  for(conv::BeamConf_t b : beams) {
    for(conv::Sample_t s : samples.at(b)) {
      std::vector <std::map<conv::KineVar_t, double>> v = data_kinematics.at(b).at(s);
      const int n_data = v.size();
      if(n_data == 0) { continue; }

      TGraph* g = new TGraph();
      TString graph_name = TString::Format("g_kinematics_data_%s_%s_%s", util::AsString(d).c_str(),
                                           util::AsString(b).c_str(), util::AsString(s).c_str());
      g->SetName(graph_name.Data());

      TString graph_title = TString::Format("%s;%s;%s", graph_name.Data(),
                                            util::AsString(conv::kvEreco).c_str(),
                                            util::AsString(conv::kvThetareco).c_str());
      g->SetTitle(graph_title.Data());

      unsigned int i = 0;
      for(std::map<conv::KineVar_t, double> event_kinematics : v) {
        g->SetPoint(i, event_kinematics.at(conv::kvEreco), event_kinematics.at(conv::kvThetareco));
        i++;
      } // loop over events

      g->Write(graph_name.Data());
      delete g;
    }
  }

  fout->Close();
}
//___________________________________________________________________________
bool EventIsMuLike(int iclass)
{
  return iclass == sk_code_1rmu;
}
//___________________________________________________________________________
bool EventIsELike(int iclass)
{
  return (iclass == sk_code_1re || iclass == sk_code_1re1dcy);
}
//___________________________________________________________________________
bool EventIsFHC(int horn_current)
{
  if (horn_current == 0) {
    std::cerr << "ERROR: got horn current code = 0. This shouldn't happen!" << std::endl;
    exit(1);
  }
  return horn_current > 0;
}
//___________________________________________________________________________
// Takes the scalar product between the beam and lepton directions
double CalcThetaDegrees(float dir[][7][3], int iclass)
{
  double costheta = 0.0;
  if(EventIsELike(iclass)) {
     costheta = dir_nu[0] * dir[0][1][0] + dir_nu[1] * dir[0][1][1] + dir_nu[2] * dir[0][1][2];
  }
  else if(EventIsMuLike(iclass)) {
     costheta = dir_nu[0] * dir[0][2][0] + dir_nu[1] * dir[0][2][1] + dir_nu[2] * dir[0][2][2];
  }
  else {
    std::cerr << "ERROR: Event was neither mu-like nor e-like. Shouldn't be here!" << std::endl;
    exit(1);
  }
  
  if( costheta<-1.00001 || costheta>+1.00001 ){
    std::cerr << std::endl;
    std::cerr << std::setprecision(10);
    std::cerr << "costheta = " << costheta << " --> outside allowed range" << std::endl;
    std::cerr << std::endl;
    exit(1);
  }
  double LepTheta        = TMath::ACos(costheta);         // radians
  double LepThetaDegrees = LepTheta / TMath::Pi() * 180.; // degrees
  if(LepTheta<0.){
    std::cerr << std::endl;
    std::cerr << "Theta lepton is negative: " << LepThetaDegrees << " degrees" << std::endl;
    std::cerr << std::endl;
    exit(1);
  }
  return LepThetaDegrees;
}
//___________________________________________________________________________
