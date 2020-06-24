#include <iostream>
#include <cassert>
#include <vector>

#include "TChain.h"
#include "TBranch.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"

// GENIE
#include "EVGCore/EventRecord.h"
#include "Ntuple/NtpMCEventRecord.h"

// NDTF
#include "AnalysisOutput/Base/TaskForceAna.h"

using namespace std;



// Chris's FGT target function, declared inside here to avoid having to split the loading of the libraries into a different file
//___________________________________________________________________________________________________________________________
int FGTwgt( ndtf::TaskForceAna & ana )
  {

    // For two types of interactions we will use the STT events
    bool neutrino_electron_scattering = ( ana.sample_reco == 25 );
    bool nu_e_CC = ( ana.sample_reco >= 17 && ana.sample_reco <= 20 ); // 4 different nu_e CC modes

    // Figure out if the reco interaction vertex is in one of the filled target modules
    bool targetModule = ( (ana.recoVtxZ > 260.608 && ana.recoVtxZ < 261.309) || // Ca target (no enclosure simulated)
                          (ana.recoVtxZ > 267.508 && ana.recoVtxZ < 272.792) || // Filled Ar 1 (including C tube)
                          (ana.recoVtxZ > 278.991 && ana.recoVtxZ < 284.275) || // Filled Ar 2 (including C tube)
                          (ana.recoVtxZ > 290.474 && ana.recoVtxZ < 295.759) ); // Filled Ar 3 (including C tube)

    // or in the empty target module
    bool emptyTargetModule = (ana.recoVtxZ > 301.917 && ana.recoVtxZ < 307.202); // C tube filled with air

    // or in the STT; we will not use events just downstream of the target region, or in the tubes between target modules
    bool STT = (ana.recoVtxZ > 315.25); // downstream boundary already enforced in FGT code

    // true target nucleus A
    int tgtA = (ana.target%10000)/10; // 10LZZZAAAI is nucleus, %10000 pulls off AAAI, integer division by 10 strips last digit

    // If it's nu-e scattering or nu_e CC, then keep the events iff they are in the STT
    if( neutrino_electron_scattering || nu_e_CC ) {
      if( STT ) return 1;
      else return 0; // don't use event in VALOR
    }

    // For nu_mu CC, and NC of any flavor: keep only target module events
    // For 2+ tracks, we get vertex resoution ~100 microns, assume perfect and keep true 40Ar and 40Ca only, with no subtraction
    // For 1 track, we only know the event happened somewhere in the target module
    // subtract the empty module, with a factor of 3 because there are 3 filled modules to subtract
    // this neglects z-dependent acceptance effects, if there are any
    // NC pi0 we probably don't have good enough direction resolution on the photons to get a great vertex
    // but let's assume we can get it to within the target module, and treat it like a 1-track event

    // isFound == true for reconstructed particles, but photons and pi0s will not have tracks
    int ntracks = 0;
    int npi0 = 0;
    std::vector<ndtf::FinalStateParticle> & recoParticles = ana.recoFSParticles;
    for( int ip = 0 ; ip < recoParticles.size() ; ip++ ) {
      ndtf::FinalStateParticle & p(recoParticles.at(ip));
      if( /*p.isFound && */ p.pdg != 22 && p.pdg != 111 ) ntracks++;
      if( /*p.isFound && */ p.pdg == 111 ) npi0++;
    }

    if( ntracks >= 2 ) {
      if( targetModule && tgtA == 40 ) return 1;
      else return 0; // don't use event in VALOR
    }

    if( ntracks == 1 || npi0 >= 1 ) {
      if( targetModule ) return 1;
      else if( emptyTargetModule ) return -3;
      else return 0;
    }

    // if there are particles, but no tracks and no pi0s then we have single gammas and we will have no vertex resolution so throw them out
    if( recoParticles.size() ) return 0;

    // Also can have 0 reconstructed particles
    // These are events that are not reconstructed, and they will have an undefined sample code
    // We want to keep them to study how often the various event categories fail to be reconstructed at all
    // But if we only keep Ca or Ar target events, then we only want to know how often this happens for true target events
    if( targetModule && tgtA == 40 ) return 1;
    else return 0;

  }
//___________________________________________________________________________________________________________________________


//~ void fgt_fixer(
  //~ TString input_file,
  //~ TString output_file = "dune_dst.root",
  //~ bool is_fgt = false
 //~ )
//~ {
int main(int argc, char *argv[])
{
  TString input_file  = argv[1];
  TString output_file = argv[2];
  bool is_fgt = true;
  
  const char * treename = "dune_dst";

  TChain chain(treename);

  chain.Add(input_file);
  assert(chain.GetEntries());
  chain.GetEntry(0);
  
  int sample=-2, reactionmode=-2;
  float Ev=-9E9,Ev_reco=-9E9,y_reco=-9E9,POTWeight=0,TGTWeight=-99;
  
  TBranch * branchptr = NULL;
  ndtf::TaskForceAna * ana = NULL;
  genie::NtpMCEventRecord * genie_mc_truth = NULL;
  
  chain.SetBranchAddress("genie_mc_truth",&genie_mc_truth,&branchptr);

  chain.SetBranchStatus("genie_mc_truth",1);
  chain.SetBranchStatus("ana",1);
  chain.SetBranchStatus("xyzzy",0);
  
  TFile * f = new TFile(output_file.Data(),"RECREATE");
  f->cd();
  
  TTree * out_tree = new TTree(treename,treename);
  
  chain.SetBranchAddress("ana",&ana,&branchptr);
  
  out_tree->Branch("genie_mc_truth","genie::NtpMCEventRecord", &genie_mc_truth, 32000, 1);
  out_tree->Branch("sample"        ,&sample        ,"sample/I"      );
  out_tree->Branch("reactionmode"  ,&reactionmode  ,"reactionmode/I");
  out_tree->Branch("Ev"            ,&Ev            ,"Ev/F"          );
  out_tree->Branch("Ev_reco"       ,&Ev_reco       ,"Ev_reco/F"     );
  out_tree->Branch("y_reco"        ,&y_reco        ,"y_reco/F"      );
  out_tree->Branch("POTWeight"     ,&POTWeight     ,"POTWeight/F"   );
  out_tree->Branch("TGTWeight"     ,&TGTWeight     ,"TGTWeight/F"   );
  
  int sample_reco = -3;
  out_tree->Branch("sample_reco"        ,&sample_reco        ,"sample_reco/I"      );
  
  //~ out_tree->SetBranchStatus("ana",0);


  for (int i = 0 ; i < chain.GetEntries(); i++) {
    chain.GetEntry(i);
    sample       = ana->sample;
    reactionmode = ana->reactionmode;
    Ev           = ana->Ev;
    Ev_reco      = ana->Ev_reco;
    y_reco       = ana->y_reco;
    POTWeight    = ana->POTWeight;
    sample_reco  = ana->sample_reco;
    
    if (is_fgt) {TGTWeight = FGTwgt(*ana);}
    else {TGTWeight = 1.;}
    
    out_tree->Fill();
    if (!(i%5000)) {
      cout<<"Event #"<<i<<" of "<<chain.GetEntries()<<endl;
      //~ out_tree->AutoSave();
    }
    genie_mc_truth->Clear();
  }
  
  cout<<"Last sample: "<<sample<<endl;
    
  out_tree->Write();
  cout<<"Written "<<out_tree->GetEntries()<<" events"<<endl;
  f->Close();
  
  return 0;
}
