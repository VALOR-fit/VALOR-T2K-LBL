
#include "Build.h"

void ParseDataMinitupleIntoXML
(
 const char * tag                     = "skdata",                        /* */
 //const char * input_file              = "$T2KOSCANA/inputs/t2ksk.14c.fqv4r0.run1-6.data.minituple.root",
 const char * input_file              = "$T2KOSCANA/inputs/data_run1to6/t2ksk.14c.fqv4r0.run1-6.mr63.data.minituple.root",
 
 char *treename = "minituple"
 )
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;
  
  
  const int good_spill_value_for_fhc = 1;
  const int good_spill_value_for_rhc = -1;
  
  const int entries_per_line = 10;
  
  TFile * fin = new TFile(input_file,"OPEN");
  assert(fin);
  TTree * tree = (TTree*) fin->Get(treename);
  
  double amom;
  double amomm;
  double amome;
  double erec;
  int iclass;
  int good_spill_flag;
  
  tree->SetBranchAddress("erec",&erec);
  tree->SetBranchAddress("amom",&amom);
  tree->SetBranchAddress("amomm",&amomm);
  tree->SetBranchAddress("amome",&amome);
  tree->SetBranchAddress("iclass",&iclass);
  tree->SetBranchAddress("good_spill_flag",&good_spill_flag);
  int nentries = tree->GetEntriesFast(); 
  
  ofstream fxml;
  ofstream ftxt;
  ostringstream xml_filename;
  xml_filename << tag << ".xml";
  fxml.open(xml_filename.str().c_str());

  fxml << "<?xml version=\"1.0\"?> " << endl;
  fxml << "<analysis_inputs>       " << endl;
  fxml << "  <superk>              " << endl;
  fxml << "    <observation>       " << endl;  

  // Dataset info
  int iexperiment = 0;
  fxml << endl;
  fxml << "      " << Form("<data_set_%d>",iexperiment) << endl;
  
  int n_fhc_mu=0,n_rhc_mu=0,n_fhc_e=0,n_rhc_e=0;

  // mu-like fakedata set  
  fxml << "        <one_mu_ring>" << endl;
  fxml << "          <ereco_values>" << endl;
  
  cout << endl;
  cout << endl;
  cout << "Mu-like FHC" << endl;
  cout << endl;
  


  for(int i=0;i<nentries;i++){
    tree->GetEntry(i);     
    if(erec <= 0.) continue;
    

    
    if(iclass != 3) {continue;} // iclass = 3 --> mu-like
    if(good_spill_flag != good_spill_value_for_fhc) {continue;}
    fxml << erec / 1000. << ", ";
    cout << erec / 1000. << ", ";
    //cout << erec / 1000. << ", ";
    n_fhc_mu++;
    if (n_fhc_mu%entries_per_line==0) {fxml<<"\n";}
  }
  fxml << endl;
  fxml << "          </ereco_values> " << endl;
  fxml << "        </one_mu_ring>    " << endl;


  // e-like fakedata set  
  fxml << "        <one_e_ring>" << endl;
  fxml << "          <ereco_values>" << endl;
  
  cout << endl;
  cout << endl;
  cout << "E-like FHC" << endl;
  cout << endl;
  
  for(int i=0;i<nentries;i++){
    tree->GetEntry(i);  

    if(erec <= 0.) continue;



    if(iclass != 1) continue;  // iclass = 1 --> e-like
    if(good_spill_flag != good_spill_value_for_fhc) {continue;}
    fxml << erec / 1000. << ", ";
    cout << erec / 1000. << ", ";
    n_fhc_e++;
    if (n_fhc_e%entries_per_line==0) {fxml<<"\n";}
  }
    
  fxml << endl;
  fxml << "          </ereco_values> " << endl;
  fxml << "        </one_e_ring>    " << endl;
  
  // mu-like-rhc fakedata set  
  fxml << "        <one_mu_ring_rhc>" << endl;
  fxml << "          <ereco_values>" << endl;
  
  cout << endl;
  cout << endl;
  cout << "Mu-like RHC" << endl;
  cout << endl;
  
  for(int i=0;i<nentries;i++){
    tree->GetEntry(i);     
    if(erec <= 0.) {continue;}
    
    if(iclass != 3) {continue;}  // iclass = 3 --> mu-like
    if(good_spill_flag != good_spill_value_for_rhc) {continue;}
    
    fxml << erec / 1000. << ", ";
    cout << erec / 1000. << ", ";
    n_rhc_mu++;
    if (n_rhc_mu%entries_per_line==0) {fxml<<"\n";}
  }
      
  fxml << endl;
  fxml << "          </ereco_values> " << endl; 
  fxml << "        </one_mu_ring_rhc>    " << endl;


  // e-like-rhc fakedata set  
  fxml << "        <one_e_ring_rhc>" << endl;
  fxml << "          <ereco_values>" << endl;
  
  cout << endl;
  cout << endl;
  cout << "E-like RHC" << endl;
  cout << endl;

  for(int i=0;i<nentries;i++){
    tree->GetEntry(i);     
    if(erec <= 0.) continue;



    if(iclass != 1) continue;  // iclass = 1 --> e-like
    if(good_spill_flag != good_spill_value_for_rhc) {continue;}
    
    fxml << erec / 1000. << ", ";
    cout << erec / 1000. << ", ";
    // cout << amom[2] / 1000. << ", ";
    // cout << amome[2] / 1000. << ", ";
    // cout << amomm[2] / 1000. << ", ";
    //1cout << endl;
    n_rhc_e++;
    if (n_rhc_e%entries_per_line==0) {fxml<<"\n";}
  }
  
      
  fxml << endl;
  fxml << "          </ereco_values> " << endl;
  fxml << "        </one_e_ring_rhc>    " << endl;

  
// 
  
  fxml << "      " << Form("</data_set_%d>\n",iexperiment) << endl;
  fxml << "    </observation>       " << endl;
  fxml << "  </superk>              " << endl;
  fxml << "</analysis_inputs>       " << endl;
  fxml.close();
  
  cout <<"\n\nEvent rates:"<<endl;
  cout <<" FHC 1Rmu: "<<n_fhc_mu<<endl;
  cout <<" FHC 1Re:  "<<n_fhc_e <<endl;
  cout <<" RHC 1Rmu: "<<n_rhc_mu<<endl;
  cout <<" RHC 1Re:  "<<n_rhc_e <<endl;
  
  return;
}





