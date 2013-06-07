/*
  This is a bit of a hack...
Takes as input something similar to a normal star*.data file (i.e. list with name, deviations for +sigma, deviations for -1 sigma, inputs)
but slightly modified (with STAT appended on each line) (also currently only expects 23-parameters)

Outputs the data in root file
 */

void text_to_root
(
 double inmix23, //input
 double indm23,  //input
 double stmix23, //stat error
 double stdm23,  //stat error
 const char * filename = "all.data"
 )
{
  TFile * f = new TFile(Form("readfile_%f_%f.root",inmix23,indm23),"RECREATE");

  TTree * tree = new TTree("tree","tree");

  //tree->ReadFile("all.data","name/C:mix23P/D:dm23P/D:mix23M/D:dm23M/D:inmix23/D:indm23/D:stmix23/D:stdm23/D");
  

  tree->ReadFile(filename,"name/C:mix12P/D:mix13P/D:mix23P/D:dcpP/D:dm12P/D:dm23P/D:mix12M/D:mix13M/D:mix23M/D:dcpM/D:dm12M/D:dm23M/D");

  TBranch *statmix = tree->Branch("stmix23", &stmix23, "stmix23/D");
  TBranch *statdm = tree->Branch("stdm23", &stdm23, "stdm23/D");
  TBranch *inputmix = tree->Branch("inmix23", &inmix23, "inmix23/D");
  TBranch *inputdm = tree->Branch("indm23", &indm23, "indm23/D");

  for(int i = 0; i < tree->GetEntries(); i++) {
    statmix->Fill();
    statdm->Fill();
    inputmix->Fill();
    inputdm->Fill();
  }

  tree->Write("tree",TObject::kOverwrite);

}
