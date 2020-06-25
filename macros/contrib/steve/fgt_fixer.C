//~ #include <iostream>
//~ #include <cassert>
//~ 
//~ #include "Cintex/Cintex.h"
//~ #include "TChain.h"
//~ #include "TFile.h"
//~ #include "TString.h"
//~ #include "TSystem.h"
//~ #include "TTree.h"
//~ #include "TROOT.h"

using namespace std;

void loadlib(TString lib)
{
  gSystem->ExpandPathName(lib);
  cout<<lib<<" load status: "<<gSystem->Load(lib.Data())<<endl;
}

void fgt_fixer(
  TString input_file,
  TString output_file = "dune_dst.root"
 )
{
  
  const char * treename = "dune_dst";
  
  loadlib("libCintex");
  ROOT::Cintex::Cintex::Enable();
  
  TString genie_libs = "$GENIE/src/scripts/gcint/genie.C";
  gSystem->ExpandPathName(genie_libs);
  gROOT->ProcessLine(".x "+genie_libs);
  
  loadlib("$DUNENDTF_LIB/libAnalysisOutputBase.so"      );
  loadlib("$DUNENDTF_LIB/libAnalysisOutput_Base_map.so" );
  loadlib("$DUNENDTF_LIB/libAnalysisOutputART.so"       );
  loadlib("$DUNENDTF_LIB/libAnalysisOutput_Base_dict.so");
  loadlib("$DUNENDTF_LIB/libTaskForceOutputer_module.so");

  TChain chain(treename);

  chain.Add(input_file);
  assert(chain.GetEntries());
  chain.GetEntry(0);
  
  int sample=-2, reactionmode=-2;
  float Ev=-9E9,Ev_reco=-9E9,y_reco=-9E9,POTWeight=0;
  
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
    out_tree->Fill();
    if (!(i%5000)) {
      cout<<"Event #"<<i<<" of "<<chain.GetEntries()<<endl;
      out_tree->AutoSave();
    }
    genie_mc_truth->Clear();
  }
  
  cout<<"Last sample: "<<sample<<endl;
    
  out_tree->Write();
  cout<<"Written "<<out_tree->GetEntries()<<" events"<<endl;
  f->Close();
  
}
