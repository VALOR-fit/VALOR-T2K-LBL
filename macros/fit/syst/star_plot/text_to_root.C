/*
  This is a bit of a hack...
Takes as input something similar to a normal star*.data file (i.e. list with name, deviations for +sigma, deviations for -1 sigma, inputs)
but slightly modified (with STAT appended on each line) (also currently only expects 23-parameters)

Outputs the data in root file
 */

void text_to_root()
{
  TFile * f = new TFile("readfile.root","RECREATE");

  TTree * tree = new TTree("tree","tree");

  tree->ReadFile("all.data","name/C:mix23P/D:dm23P/D:mix23M/D:dm23M/D:inmix23/D:indm23/D:stmix23/D:stdm23/D");

  tree->Write();

}
