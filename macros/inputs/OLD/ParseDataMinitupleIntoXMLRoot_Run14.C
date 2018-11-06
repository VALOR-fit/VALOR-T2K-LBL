
#include "Build.h"
#include "AnalysisOptions.h"

double CalcThetaDegrees(float dir[10][3]);

void ParseDataMinitupleIntoXMLRoot_Run14
(
 const char * tag                     = "SKdata_Run14",                        /* */
 const char * input_file              = "$T2KOSCANA/inputs/data_run1to4/Run1-4Data.root",
 
 char *treename = "mtuple"
 )
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;
    
  const int entries_per_line = 10;

  const int nereco_1rmu = util::NumErecoBins(conv::kSKsFCFV1RingMuLike);
  const int nq2_1rmu = util::NumQ2Bins(conv::kSKsFCFV1RingMuLike);
  TH2D *mock_data_1Rmu_0 = new TH2D("mock_data_1Rmu_0","mock_data_1Rmu_0",nq2_1rmu,util::Q2BinEdges(kSKsFCFV1RingMuLike),nereco_1rmu,util::ErecoBinEdges(kSKsFCFV1RingMuLike));
  
  const int nereco_rhc1rmu = util::NumErecoBins(conv::kSKsRHCFCFV1RingMuLike);
  const int nq2_rhc1rmu = util::NumQ2Bins(conv::kSKsRHCFCFV1RingMuLike);
  TH2D *mock_data_RHC1Rmu_0 = new TH2D("mock_data_RHC1Rmu_0","mock_data_RHC1Rmu_0",nq2_rhc1rmu,util::Q2BinEdges(kSKsRHCFCFV1RingMuLike),nereco_rhc1rmu,util::ErecoBinEdges(kSKsRHCFCFV1RingMuLike));
 
  const int nereco_1re = util::NumErecoBins(conv::kSKsFCFV1RingELike);
  const int nq2_1re = util::NumQ2Bins(conv::kSKsFCFV1RingELike);
  TH2D *mock_data_1Re_0 = new TH2D("mock_data_1Re_0","mock_data_1Re_0",nq2_1re,util::Q2BinEdges(kSKsFCFV1RingELike),nereco_1re,util::ErecoBinEdges(kSKsFCFV1RingELike));
  
  const int nereco_rhc1re = util::NumErecoBins(conv::kSKsRHCFCFV1RingELike);
  const int nq2_rhc1re = util::NumQ2Bins(conv::kSKsRHCFCFV1RingELike);
  TH2D *mock_data_RHC1Re_0 = new TH2D("mock_data_RHC1Re_0","mock_data_RHC1Re_0",nq2_rhc1re,util::Q2BinEdges(kSKsRHCFCFV1RingELike),nereco_rhc1re,util::ErecoBinEdges(kSKsRHCFCFV1RingELike));

  TFile * fin = new TFile(input_file,"OPEN");
  assert(fin);
  TTree * tree = (TTree*) fin->Get(treename);
  
  float amom;
  float amomm;
  float amome;
  float erec;
  int mode;
  float dir[10][3];
  
  tree->SetBranchAddress("erec",&erec);
  tree->SetBranchAddress("amom",&amom);
  tree->SetBranchAddress("amomm",&amomm);
  tree->SetBranchAddress("amome",&amome);
  tree->SetBranchAddress("mode",&mode);
  tree->SetBranchAddress("dir",&dir);
  int nentries = tree->GetEntriesFast(); 
  
  ofstream fxml;
  ostringstream xml_filename;
  xml_filename << tag << ".xml";
  fxml.open(xml_filename.str().c_str());

  ofstream ftxt;
  ostringstream txt_filename;
  txt_filename << tag << ".txt";
  ftxt.open(txt_filename.str().c_str());
 
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
  //fxml << "          <ereco_theta_values>" << endl;
  
  cout << endl;
  cout << endl;
  cout << "Mu-like FHC" << endl;
  cout << endl;

  int ievtot=0;      
  for(int i=0;i<nentries;i++){
    tree->GetEntry(i);     

    if(erec <= 0.) continue;
    erec /= 1000.;

    double theta = CalcThetaDegrees(dir);
    if(theta <= 0.) continue;
    
    if(mode != 0) {continue;} // In mtuple, chose mode==0 for 1Rmu and mode==1 for 1Re

    cout << "(" << erec << "," << theta << ") ";

    // Fill histogram
    int iflag = mock_data_1Rmu_0->Fill(theta,erec);
    if(iflag<=0){
      cerr << endl;
      cerr << "Flag of histo->Fill() = " << iflag << endl;
      cerr << "erec = " << erec << " - theta = " << theta << endl;
      cerr << endl;
      exit(1);
    }
    
    // Fill xml file
    // #ifdef __USE_2D_TEMPLATES__
    //     fxml << "            <ev" << ievtot << ">"
    // 	 << erec << "," << theta
    // 	 << "</ev" << ievtot << ">"
    // 	 << endl;
    //     cout << "            <ev" << ievtot << ">"
    // 	 << erec << "," << theta
    // 	 << "</ev" << ievtot << ">"
    // 	 << endl;
    // #else
    fxml << erec << ",";
    ievtot++;
    if (ievtot%entries_per_line==0) {fxml<<"\n";}
    //#endif
  }
  // #ifdef __USE_2D_TEMPLATES__
  //   fxml << "            <nevents> " << ievtot << "</nevents>" << endl;
  // #endif
    
  n_fhc_mu = ievtot;
  
  fxml << endl;
  //fxml << "          </ereco_theta_values> " << endl;
  fxml << "          </ereco_values> " << endl;
  fxml << "        </one_mu_ring>    " << endl;

  // e-like fakedata set  
  fxml << "        <one_e_ring>" << endl;
  //fxml << "          <ereco_theta_values>" << endl;
  fxml << "          <ereco_values>" << endl;
  
  cout << endl;
  cout << endl;
  cout << "E-like FHC" << endl;
  cout << endl;
  
  ievtot=0;
  for(int i=0;i<nentries;i++){
    tree->GetEntry(i);  

    if(erec <= 0.) continue;
    erec /= 1000.;

    double theta = CalcThetaDegrees(dir);
    if(theta <= 0.) continue;
    
    if(mode != 1) continue;  // mode = 1 --> e-like
    
    cout << "(" << erec << "," << theta << ") ";

    // Fill histogram
    int iflag = mock_data_1Re_0->Fill(theta,erec);
    if(iflag<=0){
      cerr << endl;
      cerr << "Flag of histo->Fill() = " << iflag << endl;
      cerr << "erec = " << erec << " - theta = " << theta << endl;
      cerr << endl;
      exit(1);
    }

    // Fill xml file
    // #ifdef __USE_2D_TEMPLATES__
    //     fxml << "            <ev" << ievtot << ">"
    // 	 << erec << "," << theta
    // 	 << "</ev" << ievtot << ">"
    // 	 << endl;
    //     cout << "            <ev" << ievtot << ">"
    // 	 << erec << "," << theta
    // 	 << "</ev" << ievtot << ">"
    // 	 << endl;
    // #else
    fxml << erec << ",";
    ievtot++;
    if (ievtot%entries_per_line==0) {fxml<<"\n";}
    //#endif
  }
  // #ifdef __USE_2D_TEMPLATES__
  //   fxml << "            <nevents> " << ievtot << "</nevents>" << endl;
  // #endif
    
  n_fhc_e = ievtot;

  fxml << endl;
  fxml << "          </ereco_values> " << endl;
  //fxml << "          </ereco_theta_values> " << endl;
  fxml << "        </one_e_ring>    " << endl;
  
  // mu-like-rhc fakedata set  
  fxml << "        <one_mu_ring_rhc>" << endl;
  fxml << "          <ereco_values>" << endl;
  //fxml << "          <ereco_theta_values>" << endl;
  fxml << endl;
  fxml << "          </ereco_values> " << endl; 
  //fxml << "          </ereco_theta_values> " << endl; 
  fxml << "        </one_mu_ring_rhc>    " << endl;

  // e-like-rhc fakedata set  
  fxml << "        <one_e_ring_rhc>" << endl;
  fxml << "          <ereco_values>" << endl;
  //fxml << "          <ereco_theta_values>" << endl;
  fxml << endl;
  fxml << "          </ereco_values> " << endl;
  //fxml << "          </ereco_theta_values> " << endl;
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
  
  fin->Close();

  TCanvas *c_1Rmu = new TCanvas("c_1Rmu","c_1Rmu");
  mock_data_1Rmu_0->Draw("colz");

  TCanvas *c_RHC1Rmu = new TCanvas("c_RHC1Rmu","c_RHC1Rmu");
  mock_data_RHC1Rmu_0->Draw("colz");

  TCanvas *c_1Re = new TCanvas("c_1Re","c_1Re");
  mock_data_1Re_0->Draw("colz");

  TCanvas *c_RHC1Re = new TCanvas("c_RHC1Re","c_RHC1Re");
  mock_data_RHC1Re_0->Draw("colz");

  // Store histograms

  char name[200];
  sprintf(name,"%s.root",tag);
  TFile *fout = new TFile(name,"RECREATE");
  mock_data_1Rmu_0->Write();
  mock_data_RHC1Rmu_0->Write();
  mock_data_1Re_0->Write();
  mock_data_RHC1Re_0->Write();
  fout->Close();

  return;
}





double CalcThetaDegrees(float dir[10][3]){
  //cout << dir[0][0] << " " << dir[0][1] << " " << dir[0][2] << endl;
  double dir_nu[3] = {0.669764, -0.742179, 0.024223}; // T2K beam direction
  double costheta = dir_nu[0] * dir[0][0] + dir_nu[1] * dir[0][1] + dir_nu[2] * dir[0][2];
  if( costheta<-1.00001 || costheta>+1.00001 ){
    cerr << endl;
    cerr << std::setprecision(10);
    cerr << "SuperkMC::FillTemplateHistogram" << endl;
    cerr << "costheta = " << costheta << " --> outside allowed range" << endl;
    cerr << endl;
    exit(1);
  }    
  double LepTheta        = TMath::ACos(costheta); // radiant    
  double LepThetaDegrees = LepTheta / TMath::Pi() * 180.; // degrees   
  if(LepTheta<0.){
    cerr << endl;
    cerr << "Theta lepton is negative: " << LepThetaDegrees << " degrees" << endl;
    cerr << endl;
    exit(1);
  }
  return LepThetaDegrees;
}
