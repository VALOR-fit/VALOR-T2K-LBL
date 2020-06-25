//____________________________________________________________________________
/*!

  \app       valor_dune-make_caf_minidst

  \brief     Converts input DUNE FD analysis DSTs to the VALOR mini-DST format.

  \syntax    valor_dune-make_nd_minidst 
  --input  input_analysis_file.root
  --output output_file_containing_valor_minidst.root
  --format format_of_the_input_file <hpgartpc, lartpc, fgt>
  [--disable-checks]

  \author    Steve Dennis <s.r.dennis \at liverpool.ac.uk>
             University of Liverpool

  \created   August 2016

  \copyright (c) 2010-2016, The VALOR Fitting Group, All Rights Reserved

*/
//____________________________________________________________________________

#include <string>
#include <iostream>
#include <sstream>
#include <map>
#include <utility>
#include <algorithm>

//
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TPRegexp.h>

// valor/dune
#include "valor/Conventions/BeamConf.h"
#include "valor/Conventions/Constants.h"
#include "valor/Conventions/Sample.h"
#include "valor/Conventions/ReacMode.h"
#include "valor/Conventions/NuclearTarget.h"
#include "valor/Messenger/Messenger.h"
#include "valor/Utilities/Env.h"
#include "valor/CmdLine/CmdLnArgList.h"
#include "valor/Utilities/UtlNumStat.h"
#include "valor/Utilities/UtlStyle.h"
#include "valor/Utilities/UtlStr.h"
#include "valor/Utilities/UtlSystem.h"
#include "valor/Utilities/Rnd.h"



#include "valor/GENIEUtilities/GENIEIncludes.h"

using std::string;
using std::ostringstream;

static const string gTreeName     = "caf";
static const string gMetaTreeName = "meta";

void GetCommandLineArgs (int /*argc*/, char ** /*argv*/);

static const int gPrintEveryNEvents = 50000;

//___________________________________________________________________
int main(int argc, char ** argv)
{
  
  
  // Cuts from Cal
  std::map<std::string,double> cut_upper_limit_numu;
  std::map<std::string,double> cut_upper_limit_nue;
  std::map<std::string,double> cut_lower_limit_numu;
  std::map<std::string,double> cut_lower_limit_nue;

  cut_lower_limit_numu["cvnnumu"] = 0.50;
  cut_lower_limit_nue ["cvnnue" ] = 0.85;
  
  cut_upper_limit_nue ["cvnnumu"] = 0.50;
  cut_upper_limit_numu["cvnnue" ] = 0.85;

  cut_lower_limit_numu["vtx_x"] = -310. ;
  cut_upper_limit_numu["vtx_x"] =  310. ;
  cut_lower_limit_numu["vtx_y"] = -550. ;
  cut_upper_limit_numu["vtx_y"] =  550. ;
  cut_lower_limit_numu["vtx_z"] =  50.  ;
  cut_upper_limit_numu["vtx_z"] =  1244.;

  cut_lower_limit_nue ["vtx_x"] = -310. ;
  cut_upper_limit_nue ["vtx_x"] =  310. ;
  cut_lower_limit_nue ["vtx_y"] = -550. ;
  cut_upper_limit_nue ["vtx_y"] =  550. ;
  cut_lower_limit_nue ["vtx_z"] =  50.  ;
  cut_upper_limit_nue ["vtx_z"] =  1244.;

  using namespace valor;
  
  
 
  CmdLnArgList arg;
  arg.Read(argc,argv,"$VALOR_EXPT/config/apps/make_caf_minidst.arg");
  LOG("valor_dune", pNOTICE) << arg;
  
  std::string thresholds = util::CompleteConfigPath(arg.GetString("mesg-thresholds"));
  Messenger::Instance()->Configure(thresholds);
  
  string opt_output_name       = arg.GetString("output");
  string opt_input_name        = arg.GetString("input");
  bool   opt_disable_checks    = arg.GetBool("disable-checks");
  bool   disable_extra_cuts    = arg.GetBool("disable-extra-cuts");
  bool   split_sample_by_res   = arg.GetBool("split-sample-by-res");
  double nc_factor             = arg.GetDouble("nc-factor");
   
  int    valor_minidst_sample       = -1;
  int    valor_minidst_beam         = -1;
  int    valor_minidst_reac_mode    = -1;  // see Conventions/ReacMode.h 
  int    valor_minidst_target       = -1;  // nuclear target -  GENIE PDG convention
  double valor_minidst_enu_true     = -1.; // true neutrino energy, in GeV
  double valor_minidst_Q2_true      = -1.; // true neutrino energy, in GeV
  double valor_minidst_enu_reco     = -1.; // reco neutrino energy, in GeV
  double valor_minidst_y_reco       = -1.; // inelasticity
  double valor_minidst_wght_pot     = -1.; // POT weight 
  double valor_minidst_wght         = -1.; // misc weight 
  double valor_minidst_wght_smear   = -1.; // energy smearing weight 
  
  TFile fout(opt_output_name.c_str(),"recreate");
  fout.cd();

  TTree * valor_minidst = 
    new TTree("valor_dune_caf_minidst","Mini-DST input for VALOR/DUNE Fit");

  valor_minidst->Branch("beam"      , &valor_minidst_beam        ,"beam/I"       );
  valor_minidst->Branch("sample"    , &valor_minidst_sample      ,"sample/I"     );
  valor_minidst->Branch("reac_mode" , &valor_minidst_reac_mode   ,"reac_mode/I"  );
  valor_minidst->Branch("target"    , &valor_minidst_target      ,"target/I"     );
  valor_minidst->Branch("enu_true"  , &valor_minidst_enu_true    ,"enu_true/D"   );
  valor_minidst->Branch("Q2_true"   , &valor_minidst_Q2_true     ,"q2_true/D"    );
  valor_minidst->Branch("enu_reco"  , &valor_minidst_enu_reco    ,"enu_reco/D"   );
  valor_minidst->Branch("y_reco"    , &valor_minidst_y_reco      ,"y_reco/D"     );
  valor_minidst->Branch("wght_pot"  , &valor_minidst_wght_pot    ,"wght_pot/D"   );
  valor_minidst->Branch("wght"      , &valor_minidst_wght        ,"wght/D"       );
  valor_minidst->Branch("wght_smear", &valor_minidst_wght_smear  ,"wght_smear/D" );
  
  
  double valor_minidst_cvnnumu  = -1.; 
  double valor_minidst_cvnnue   = -1.; 
  valor_minidst->Branch("cvnnumu"   , &valor_minidst_cvnnumu    ,"cvnnumu/D"       );
  valor_minidst->Branch("cvnnue"    , &valor_minidst_cvnnue     ,"cvnnue/D"       );
  
  // First get the POT metadata
  TChain * metadata_tree = new TChain(gMetaTreeName.c_str());
  metadata_tree->Add(opt_input_name.c_str());
  if(!metadata_tree || !metadata_tree->GetEntries()) {
    VLOG("make_caf_minidst", pFATAL) 
      << "Unable to read " << gMetaTreeName
      <<" from "<< opt_input_name;
    valor::util::Die(-1);
  }
  double br_pot = 0.;
  metadata_tree->SetBranchAddress("pot",&br_pot);
  std::vector<double> pot_per_file;
  int last_file_number = -1;
  for (Long64_t ientry = 0 ; ientry < metadata_tree->GetEntries() ; ientry++) {
    metadata_tree->GetEntry(ientry);
    int this_file_number = metadata_tree->GetTreeNumber();
    if (this_file_number != last_file_number) {
      pot_per_file.push_back(0.);
      last_file_number = this_file_number;
    }
    pot_per_file.back() += br_pot;
  }
  VLOG("make_can_minidst",pNOTICE) << "Reading "<<pot_per_file.size()<<" files";
  for (size_t ifile = 0 ; ifile < pot_per_file.size() ; ifile++) {
    VLOG("make_can_minidst",pNOTICE) << "Input analysis tree in file "<<ifile<<" contains "<<pot_per_file[ifile]<<" POT";
  }
  
  // Then the analysis tree
  TChain * analysis_tree = new TChain(gTreeName.c_str());
  analysis_tree->Add(opt_input_name.c_str());
  if(!analysis_tree || !analysis_tree->GetEntries()) {
    VLOG("make_caf_minidst", pFATAL) 
      << "Unable to read " << gTreeName
      <<" from "<< opt_input_name;
    valor::util::Die(-1);
  }

  Long64_t nev = analysis_tree->GetEntries();
  VLOG("make_caf_minidst", pNOTICE) 
    << "Input analysis tree has " << nev << ((nev==1) ? " entry." : " entries.");
  analysis_tree->GetEntry(0);
  
  int run = -1;
  double enu_true         =  -999.;
  double enu_reco_nue     =  -999.;
  double enu_reco_numu    =  -999.;
  double had_en_reco_nue  =  -999.;
  double had_en_reco_numu =  -999.;
  double Q2_true          =  -999.;
  
  int init_pdg  = -1;
  int final_pdg  = -1;
  
  int is_fhc   = -1;
  int mode     = -1;
  int nin      = -1;
  int nip      = -1;
  int nipip    = -1;
  int nipim    = -1;
  int nipi0    = -1;
  int nikp     = -1;
  int nikm     = -1;
  int nik0     = -1;
  int niem     = -1;
  int niother  = -1;
  
  int tgt = -1;
  
  int iscc   = -1;
  
  TBranch * branchptr = 0;
  vector<TBranch*> active_branches;
  
  // Hack for "isFHC" while that's broken in CAFs
  TString filename(opt_input_name.c_str());
  if      (filename.Contains("FHC")) {is_fhc = 1;}
  else if (filename.Contains("RHC")) {is_fhc = 0;}
  else {
    std::cerr<<"Could not work out if beam is FHC or RHC from filename "<<opt_input_name<<"  - exiting."<<std::endl;
    util::Die(-1);
  }
    
  //~ analysis_tree->SetBranchAddress("isFHC"        ,  &is_fhc          , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("isCC"         ,  &iscc            , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("mode"         ,  &mode            , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("Ev"           ,  &enu_true        , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("Ev_reco_nue"  ,  &enu_reco_nue    , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("Ev_reco_numu" ,  &enu_reco_numu   , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("RecoHadEnNue" ,  &had_en_reco_nue , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("RecoHadEnNumu",  &had_en_reco_numu, &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("Q2"           ,  &Q2_true         , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("nuPDGunosc"   ,  &init_pdg        , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("nuPDG"        ,  &final_pdg       , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("nP"           ,  &nip             , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("nN"           ,  &nin             , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("nipip"        ,  &nipip           , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("nipim"        ,  &nipim           , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("nipi0"        ,  &nipi0           , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("nikp"         ,  &nikp            , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("nikm"         ,  &nikm            , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("nik0"         ,  &nik0            , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("niem"         ,  &niem            , &branchptr); active_branches.push_back(branchptr);
  analysis_tree->SetBranchAddress("niother"      ,  &niother         , &branchptr); active_branches.push_back(branchptr);
  //~ analysis_tree->SetBranchAddress("tgt"          ,  &tgt              , &branchptr); active_branches.push_back(branchptr);
  
  // Branches for the cut variables
  std::set<std::string> cut_branch_names;
  for (const auto & it : cut_lower_limit_numu) { cut_branch_names.insert(it.first);  }
  for (const auto & it : cut_upper_limit_numu) { cut_branch_names.insert(it.first);  }
  for (const auto & it : cut_lower_limit_nue ) { cut_branch_names.insert(it.first);  }
  for (const auto & it : cut_upper_limit_nue ) { cut_branch_names.insert(it.first);  }
  
  std::map<std::string,double> cut_branches;
  for (const std::string & branchname : cut_branch_names) {
    double & branch_ref = cut_branches[branchname];
    analysis_tree->SetBranchAddress(branchname.c_str(), &branch_ref, &branchptr); active_branches.push_back(branchptr);
  }
  
  // Get list of systematics in tree
  std::vector<std::string> syst_list;
  TPRegexp syst_pattern("(\\w+)_cvwgt");
  TObjArray * branch_list = analysis_tree->GetListOfBranches();
  for (size_t ib = 0 ; ib < branch_list->GetEntries(); ib++) {
    const char * branchname = branch_list->At(ib)->GetName();
    TObjArray * matches = syst_pattern.MatchS(branchname);
    if (matches->GetEntries()==2) {
      syst_list.push_back(std::string(dynamic_cast<TObjString*>(matches->At(1))->GetString().Data()));
    }
    delete matches;
  }
  
  std::map<std::string,double> cvwght_branches;
  for (const std::string & syst_name : syst_list) {
    double & branch_ref = cvwght_branches[syst_name];
    std::string branchname(Form("%s_cvwgt",syst_name.c_str()));
    analysis_tree->SetBranchAddress(branchname.c_str(), &branch_ref, &branchptr); active_branches.push_back(branchptr);
    // And the VALOR branch for the CVWght
    std::string valor_branch_name (Form("autowght_%s",syst_name.c_str()));
    valor_minidst->Branch(valor_branch_name.c_str(),&branch_ref,Form("%s/D",valor_branch_name.c_str()));
  }
  
  // Temporarily set target to argon
  tgt = conv::kntArgon;
  
  // Why the hell doesn't this work? Probably indicate of a bug elsewhere.
  //~ analysis_tree->SetBranchStatus("*",0);
  
  // Only read the branches we're using and their sub-branches
  for (unsigned int ib = 0 ; ib < active_branches.size(); ib++) {
    active_branches[ib]->SetStatus(1);
    TObjArray * bar = active_branches[ib]->GetListOfBranches();
    for (int ie = 0 ; ie < bar->GetEntries(); ie++) {
      TBranch * b = dynamic_cast<TBranch*>(bar->At(ie));
      assert(b);
      b->SetStatus(1);
    }
  }
  
  //
  // event loop
  //

  Long64_t nev_accept = 0;
  
  std::vector<int> nev_accept_by_mode;
  for (int im = 0 ; im < valor::conv::kNReacModes; im++) {
    nev_accept_by_mode.push_back(0);
  }

  for(Long64_t iev = 0; iev < nev; iev++) {

    if (gPrintEveryNEvents && !(iev%gPrintEveryNEvents)) {
      VLOG("make_caf_minidst", pNOTICE) << "** Reading event: " << iev;
    }
    else {
      VLOG("make_caf_minidst", pDEBUG) << "** Reading event: " << iev;
    }

    analysis_tree->GetEntry(iev);
    
    VLOG("make_caf_minidst",pDEBUG) << "Read event: "
      << "Etrue            = "<<enu_true  << " GeV; "
      << "Ereco_nue        = "<<enu_reco_nue  << " GeV; "
      << "Ereco_numu       = "<<enu_reco_numu  << " GeV; "
      << "RecoHadEnNue     = "<<had_en_reco_nue  << ";"
      << "RecoHadEnNuMu    = "<<had_en_reco_numu  << ";";
    
    // Construct the reaction mode
    conv::ReacType_t t = conv::kmtUndefined;
    if (iscc) {
      if      (mode == genie::kScQuasiElastic      ) { t = conv::kmtCCQE;    }
      else if (mode == genie::kScMEC               ) { t = conv::kmtCCMEC;   }
      else if (mode == genie::kScCoherent          ) { t = conv::kmtCCcoh;   }
      else if (mode == genie::kScInverseMuDecay    ) { t = conv::kmtCCIMD;   }
      else if (mode == genie::kScNuElectronElastic ) { t = conv::kmtCCELScat;}
      
      else { // start counting pions
        if      (nikp || nik0 || nikm || niother) { t = conv::kmtCCoth;}
        else if (nipip+nipim == 1 && nipi0 == 0)  { t = conv::kmtCC1piC;}
        else if (nipip+nipim == 0 && nipi0 == 1)  { t = conv::kmtCC1pi0;}
        else if (nipip+nipim == 2 && nipi0 == 0)  { t = conv::kmtCC2piC;}
        else if (nipip+nipim == 0 && nipi0 == 2)  { t = conv::kmtCC2pi0;}
        else if (nipip+nipim == 1 && nipi0 == 1)  { t = conv::kmtCCpiCpi0;}
        else                                      { t = conv::kmtCCoth;}
      }
    }
    else {
      if      (mode == genie::kScQuasiElastic      ) { t = conv::kmtNCEL;    }
    //else if (mode == genie::kScMEC               ) { t = conv::kmtNCMEC;   }
      else if (mode == genie::kScCoherent          ) { t = conv::kmtNCcoh;   }
    //else if (mode == genie::kScNuElectronElastic ) { t = conv::kmtNCELScat;}
      
      else { // start counting pions
        if      (nikp || nik0 || nikm || niother) { t = conv::kmtNCoth;}
        else if (nipip+nipim == 1 && nipi0 == 0)  { t = conv::kmtNC1piC;}
        else if (nipip+nipim == 0 && nipi0 == 1)  { t = conv::kmtNC1pi0;}
      //else if (nipip+nipim == 2 && nipi0 == 0)  { t = conv::kmtNC2piC;}
      //else if (nipip+nipim == 0 && nipi0 == 2)  { t = conv::kmtNC2pi0;}
      //else if (nipip+nipim == 1 && nipi0 == 1)  { t = conv::kmtNCpiCpi0;}
        else                                      { t = conv::kmtNCoth;}
      }
    }
    // Apparently we're setting iscc to true for these events? Who knows
    //~ else if (isNuEl) {
      //~ t = conv::kmtELScat;
    //~ }
  
    //~ std::cerr<<util::AsString(t)<<"\t"<<init_pdg<<"\t"<<final_pdg<<std::endl;
    
    conv::NuFlavour_t init_flav  = valor::util::FlavourFromPDG(init_pdg );
    conv::NuFlavour_t final_flav = valor::util::FlavourFromPDG(final_pdg);
    
    // If it's NC, we'll ignore oscillated events (to avoid double counting)
    if (!iscc) {
      if (final_flav != init_flav) {continue;}
    }
    
    // Only have the combination mode for final state nues
    if (t == conv::kmtCCELScat && (final_flav == conv::kNuE || final_flav == conv::kNuEBar)) {
      t = conv::kmtELScat;
    }
    // NCel is not in play yet
    if (t == conv::kmtNCEL) {
      if (final_flav == conv::kNuE || final_flav == conv::kNuEBar) {
        t = conv::kmtELScat;
      }
      else {
        t = conv::kmtNCoth;
      }
    }
    
    // Get the real valor mode
    conv::ReacMode_t  reac_mode  = valor::util::ReacTypeToMode(t,init_flav,final_flav);
    
    //~ std::cerr<<util::AsString(t)<<"\t"<<util::AsString(init_flav)<<"\t"<<util::AsString(final_flav)<<std::endl;
    
    if(!opt_disable_checks) {
      //
      // ignore events that can not be used in the VALOR analysis
      //
      if(reac_mode == valor::conv::kmUndefined) {
        VLOG("make_caf_minidst", pWARN) << "- Rejecting event [bad mode]\n"
        << util::AsString(init_flav)<<" -> "<<util::AsString(final_flav)<<"\t"<<util::AsString(t);
        continue;
      }
    }
    
    conv::Sample_t sample = conv::ksUndefined;
    
    // Do the cuts
    sample = conv::ksNuELikeCChigh;
    for (const auto & cut_upper_limit : cut_upper_limit_nue) {
      const std::string & branch_name = cut_upper_limit.first;
      if (cut_branches[branch_name] > cut_upper_limit.second) {
        sample = conv::ksUndefined;
      }
    }
    for (const auto & cut_lower_limit : cut_lower_limit_nue) {
      const std::string & branch_name = cut_lower_limit.first;
      if (cut_branches[branch_name] < cut_lower_limit.second) {
        sample = conv::ksUndefined;
      }
    }
    if (sample==conv::ksUndefined) {
      sample = conv::ksNuMuLikeCChigh;
      for (const auto & cut_upper_limit : cut_upper_limit_numu) {
        const std::string & branch_name = cut_upper_limit.first;
        if (cut_branches[branch_name] > cut_upper_limit.second) {
          sample = conv::ksUndefined;
        }
      }
      for (const auto & cut_lower_limit : cut_lower_limit_numu) {
        const std::string & branch_name = cut_lower_limit.first;
        if (cut_branches[branch_name] < cut_lower_limit.second) {
          sample = conv::ksUndefined;
        }
      }
    }
    
    double enu_reco = (util::IsELike(sample) ? enu_reco_nue : enu_reco_numu);
    double y_reco   = (util::IsELike(sample) ? had_en_reco_nue/enu_reco : had_en_reco_numu / enu_reco);
    if(!opt_disable_checks) {
      
      //
      // ignore events that can not be used in the VALOR analysis
      //
      if(enu_reco <= 0.) {
        VLOG("make_caf_minidst", pNOTICE) << "- Rejecting event [Ev (reco) <= 0]: "<<enu_reco;
        continue;
      }
      
      if(sample == conv::ksUndefined) {
        VLOG("make_caf_minidst", pINFO) << "- Rejecting event not in one of our samples";
        continue;
      }

      VLOG("make_caf_minidst", pINFO) << "+ Accepting event";
    }
      
    nev_accept++;
    nev_accept_by_mode[reac_mode]++;

      //~ VLOG("make_caf_minidst", pDEBUG) 
        //~ << "GENIE MC truth info: " << event;
        
    double wght = 1.;
    // Write in the thing to apply the cvwghts here
    
    valor_minidst_beam      = (is_fhc ? conv::kbLBNF_FHC : conv::kbLBNF_RHC);
    valor_minidst_sample    = sample; 
    valor_minidst_reac_mode = reac_mode;
    valor_minidst_enu_true  = enu_true; 
    valor_minidst_enu_reco  = enu_reco; 
    valor_minidst_y_reco    = y_reco;
    valor_minidst_wght      = wght;
    valor_minidst_wght_pot  = 1./pot_per_file[analysis_tree->GetTreeNumber()];
    valor_minidst_wght_smear= 1.; 
    valor_minidst_Q2_true   = Q2_true;
      
    valor_minidst_target    = tgt;
    
    valor_minidst_cvnnumu = cut_branches.at("cvnnumu");
    valor_minidst_cvnnue  = cut_branches.at("cvnnue" );
      
    VLOG("make_caf_minidst", pINFO) 
      << "\n Information copied to VALOR DUNE CAF Mini-DST: "
      << "\n Beam          = " << valor_minidst_beam   
      << "\n Sample        = " << valor_minidst_sample    
      << "\n Reaction mode = " << valor_minidst_reac_mode 
      << "\n Target        = " << valor_minidst_target
      << "\n Q2 (true)     = " << valor_minidst_Q2_true << " GeV^2"
      << "\n Ev (true)     = " << valor_minidst_enu_true << " GeV"
      << "\n Ev (reco)     = " << valor_minidst_enu_reco << " GeV" 
      << "\n y (reco)      = " << valor_minidst_y_reco    
      << "\n POT weight    = " << valor_minidst_wght_pot
      << "\n Smear weight  = " << valor_minidst_wght_smear;
    
    valor_minidst->Fill();
  } // event loop
  
  for (unsigned int im = 0 ; im < nev_accept_by_mode.size(); im++) {
    valor::conv::ReacMode_t m = (valor::conv::ReacMode_t)im;
    VLOG("make_caf_minidst", pNOTICE) << "Selected "<<nev_accept_by_mode[im]<<" events of mode "<<valor::util::AsString(m);
  }
  
  double f = 100. * (double)nev_accept / (double)nev;
  VLOG("make_caf_minidst", pNOTICE) <<
    "Used " << f << "% (" << nev_accept << "/" << nev 
    << ") of input events.";
    
  //~ VLOG("make_caf_minidst", pNOTICE) <<
    //~ "Total POT read: "<<total_pot_read;

  // Saving VALOR mini-FST to output file 
  VLOG("make_caf_minidst", pNOTICE) 
    << "Saving VALOR/DUNE mini-DST to: " << opt_output_name;
  valor_minidst->Write();  
  fout.Close();

  VLOG("make_caf_minidst", pNOTICE)  << "Done!";

  return 0;
}
//___________________________________________________________________
