//
// ------------------------------------------------------------------------------------------------
// Generate / store toy-MC experiments.
// The data can be stored in an XML and/or ROOT format.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include "Build.h"

void WriteXMLfromROOTfile
(
 const char * tag                     = "fakedata_201411_F",                        /* */
 
 //char *path_mulike = "$T2KOSCANA/inputs/fakedata/201411.fakedata/minituples/t2ksk.fakedt.minituple.setA_p250kA.0p681e21_pot.mask.root",	
 //char *path_mulike = "$T2KOSCANA/inputs/fakedata/201411.fakedata/minituples/t2ksk.fakedt.minituple.setB_p250kA.0p681e21_pot.mask.root",	
 //char *path_mulike = "$T2KOSCANA/inputs/fakedata/201411.fakedata/minituples/t2ksk.fakedt.minituple.setC_p250kA.0p681e21_pot.mask.root",	
 //char *path_mulike = "$T2KOSCANA/inputs/fakedata/201411.fakedata/minituples/t2ksk.fakedt.minituple.setD_p250kA.0p681e21_pot.mask.root",	
 //char *path_mulike = "$T2KOSCANA/inputs/fakedata/201411.fakedata/minituples/t2ksk.fakedt.minituple.setE_p250kA.0p681e21_pot.mask.root",	
 char *path_mulike = "$T2KOSCANA/inputs/fakedata/201411.fakedata/minituples/t2ksk.fakedt.minituple.setF_p250kA.0p681e21_pot.mask.root",

 char *path_elike,	

 //char *path_mulike_rhc = "$T2KOSCANA/inputs/fakedata/201411.fakedata/minituples/t2ksk.fakedt.minituple.setA_m250kA.0p051e21_pot.mask.root",	
 //char *path_mulike_rhc = "$T2KOSCANA/inputs/fakedata/201411.fakedata/minituples/t2ksk.fakedt.minituple.setB_m250kA.0p051e21_pot.mask.root",	
 //char *path_mulike_rhc = "$T2KOSCANA/inputs/fakedata/201411.fakedata/minituples/t2ksk.fakedt.minituple.setC_m250kA.0p051e21_pot.mask.root",	
 //char *path_mulike_rhc = "$T2KOSCANA/inputs/fakedata/201411.fakedata/minituples/t2ksk.fakedt.minituple.setD_m250kA.0p051e21_pot.mask.root",	
 //char *path_mulike_rhc = "$T2KOSCANA/inputs/fakedata/201411.fakedata/minituples/t2ksk.fakedt.minituple.setE_m250kA.0p051e21_pot.mask.root",
 char *path_mulike_rhc = "$T2KOSCANA/inputs/fakedata/201411.fakedata/minituples/t2ksk.fakedt.minituple.setF_m250kA.0p051e21_pot.mask.root",

 char *path_elike_rhc,

 char *treename = "minituple"
 )
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;
  
  ofstream fxml;
  ofstream ftxt;
  ostringstream xml_filename;
  xml_filename << tag << ".xml";
  fxml.open(xml_filename.str().c_str());
  
  fxml << endl;
  fxml << "<!-- ";
  fxml << "fakedata set ";
  fxml << tag;
  fxml << " -->" << endl;
  fxml << endl;

  fxml << "<?xml version=\"1.0\"?> " << endl;
  fxml << "<analysis_inputs>       " << endl;
  fxml << "  <superk>              " << endl;
  fxml << "    <observation>       " << endl;  

  // Dataset info
  int iexperiment = 0;
  fxml << endl;
  fxml << "      " << Form("<data_set_%d>",iexperiment) << endl;


  // mu-like fakedata set  
  fxml << "	<one_mu_ring>" << endl;
  fxml << "          <ereco_values>" << endl;

  if(path_mulike){
    TFile *fin_mulike = new TFile(path_mulike,"OPEN");
    if(fin_mulike){
      TTree *tree_mulike = (TTree*) fin_mulike->Get(treename);
      double erec;
      tree_mulike->SetBranchAddress("erec",&erec);
      int nentries = tree_mulike->GetEntriesFast(); 
      for(int i=0;i<nentries;i++){
	tree_mulike->GetEntry(i);     
	if(erec <= 0.) continue;
	fxml << erec << ", ";
	//cout << erec << ", ";
      }
    }
  }
  fxml << endl;
  fxml << "          </ereco_values> " << endl;
  fxml << "        </one_e_ring>    " << endl;


  // e-like fakedata set  
  fxml << "	<one_e_ring>" << endl;
  fxml << "          <ereco_values>" << endl;
  if(path_elike){
    TFile *fin_elike = new TFile(path_elike,"OPEN");
    if(fin_elike){
      TTree *tree_elike = (TTree*) fin_elike->Get(treename);
      double erec;
      tree_elike->SetBranchAddress("erec",&erec);
      int nentries = tree_elike->GetEntriesFast(); 
      for(int i=0;i<nentries;i++){
	tree_elike->GetEntry(i);     
	if(erec <= 0.) continue;
	fxml << erec << ", ";
      }
    }
  }
  fxml << endl;
  fxml << "          </ereco_values> " << endl;
  fxml << "        </one_e_ring>    " << endl;


  // mu-like-rhc fakedata set  
  fxml << "	<one_mu_ring_rhc>" << endl;
  fxml << "          <ereco_values>" << endl;
  if(path_mulike_rhc){
    TFile *fin_mulike_rhc = new TFile(path_mulike_rhc,"OPEN");
    if(fin_mulike_rhc)){
      TTree *tree_mulike_rhc = (TTree*) fin_mulike_rhc->Get(treename);
      double erec;
      tree_mulike_rhc->SetBranchAddress("erec",&erec);
      int nentries = tree_mulike_rhc->GetEntriesFast(); 
      for(int i=0;i<nentries;i++){
	tree_mulike_rhc->GetEntry(i);     
	if(erec <= 0.) continue;
	fxml << erec << ", ";
      }
    }
  }
  fxml << endl;
  fxml << "          </ereco_values> " << endl;
  fxml << "        </one_mu_ring_rhc>    " << endl;


  // e-like-rhc fakedata set  
  fxml << "	<one_e_ring_rhc>" << endl;
  fxml << "          <ereco_values>" << endl;
  if(path_elike_rhc){
    TFile *fin_elike_rhc = new TFile(path_elike_rhc,"OPEN");
    if(fin_elike_rhc){
      TTree *tree_elike_rhc = (TTree*) fin_elike_rhc->Get(treename);
      double erec;
      tree_elike_rhc->SetBranchAddress("erec",&erec);
      int nentries = tree_elike_rhc->GetEntriesFast(); 
      for(int i=0;i<nentries;i++){
	tree_elike_rhc->GetEntry(i);     
	if(erec <= 0.) continue;
	fxml << erec << ", ";
      }
    }
  }
  fxml << endl;
  fxml << "          </ereco_values> " << endl;
  fxml << "        </one_e_ring_rhc>    " << endl;

  
// 
  
  fxml << "      " << Form("</data_set_%d>\n",iexperiment) << endl;
  fxml << "    </observation>       " << endl;
  fxml << "  </superk>              " << endl;
  fxml << "</analysis_inputs>       ";// << endl;
  fxml.close();
  
  return;
}





