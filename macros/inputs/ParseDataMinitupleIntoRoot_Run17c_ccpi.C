
#include "Build.h"
#include "AnalysisOptions.h"

double CalcThetaDegreesData(float dir[10][3]);

void ParseDataMinitupleIntoRoot_Run17c_ccpi
(
 const char * tag                     = "SKdata_Run17c_ccpi_1D",                        /* */
 const char * input_file              = "$T2KOSCANA/inputs/data_run1to7c/t2ksk.14c.fqv4r0.run1-7.mr68.data.minituple.root",
 //const char * input_file              = "t2ksk.14c.fqv4r0.run1-6.mr63.data.minituple.root",
 
 char *treename = "minituple"
 )
{
  using namespace std;
  using namespace valor;
  using namespace valor::conv;
  
  
  const int good_spill_value_for_fhc = 1;
  const int good_spill_value_for_fhc_205kA = 2; // FHC run3 with 205 kA
  const int good_spill_value_for_rhc = -1;
  
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

  const int nereco_multire = util::NumErecoBins(conv::kSKsFCFVMultiRingELike);
  const int nq2_multire = util::NumQ2Bins(conv::kSKsFCFVMultiRingELike);
  TH2D *mock_data_MultiRe_0 = new TH2D("mock_data_MultiRe_0","mock_data_MultiRe_0",nq2_multire,util::Q2BinEdges(kSKsFCFVMultiRingELike),nereco_multire,util::ErecoBinEdges(kSKsFCFVMultiRingELike));

  TFile * fin = new TFile(input_file,"OPEN");
  assert(fin);
  TTree * tree = (TTree*) fin->Get(treename);
  
  float amom;
  float amomm;
  float amome;
  double erec;
  int iclass;
  int good_spill_flag;
  float dir[10][3];
  
  tree->SetBranchAddress("erec",&erec);
  tree->SetBranchAddress("amom",&amom);
  tree->SetBranchAddress("amomm",&amomm);
  tree->SetBranchAddress("amome",&amome);
  tree->SetBranchAddress("iclass",&iclass);
  tree->SetBranchAddress("good_spill_flag",&good_spill_flag);
  tree->SetBranchAddress("dir",&dir);
  int nentries = tree->GetEntriesFast();
  
  
  // Dataset info
  int iexperiment = 0;
  
  int n_fhc_mu=0,n_rhc_mu=0,n_fhc_e=0,n_rhc_e=0,n_multie=0;
  
  SuperkMC skmc;
  
  for(int i=0;i<nentries;i++){
    tree->GetEntry(i);
    
    if(erec < 0.) continue;
    erec /= 1000.;
    
    //double theta = CalcThetaDegreesData(dir);
    //cout << "theta1 = " << theta << " - ";
    
    double theta = CalcThetaDegrees(dir);
    //cout << "theta2 = " << theta << endl;
    
    if(theta < 0.) continue;
    
    if(iclass == 1){ // e-like
      
      if(erec <= 1.25){
        
        if((good_spill_flag == good_spill_value_for_fhc) || (good_spill_flag == good_spill_value_for_fhc_205kA)){ // fhc
          mock_data_1Re_0->Fill(theta,erec);
          n_fhc_e++;
          //cout << "Erec = " << erec << " - Theta = " << theta << endl;
        }
        else if(good_spill_flag == good_spill_value_for_rhc){ // rhc
          mock_data_RHC1Re_0->Fill(theta,erec);
          n_rhc_e++;
          //cout << "Erec = " << erec << " - Theta = " << theta << endl;
        }
        
      } // erec <= 1.25
      
      
    }
    else if(iclass == 4 && erec <= 1.25){ //MultiRe
      if((good_spill_flag == good_spill_value_for_fhc) || (good_spill_flag == good_spill_value_for_fhc_205kA)){ // fhc
        mock_data_MultiRe_0->Fill(theta,erec);
        n_multie++;
        //cout << "Erec = " << erec << " - Theta = " << theta << endl;
      }
    }
    else if(iclass == 3){ // mu-like
                          //if(good_spill_flag == good_spill_value_for_fhc){ // fhc
      if((good_spill_flag == good_spill_value_for_fhc) || (good_spill_flag == good_spill_value_for_fhc_205kA)){ // fhc
        mock_data_1Rmu_0->Fill(theta,erec);
        n_fhc_mu++;
        //cout << "Erec = " << erec << " - Theta = " << theta << endl;
        
      }
      else if(good_spill_flag == good_spill_value_for_rhc){ // rhc
        mock_data_RHC1Rmu_0->Fill(theta,erec);
        n_rhc_mu++;
        cout << "Erec = " << erec << " - Theta = " << theta << endl;
      }
    }
    
  }
  
  
  
  
  
  
  
  
  
  // Print # of events
  
  cout <<"\n\nEvent rates:"<<endl;
  cout <<" FHC 1Rmu: "<<n_fhc_mu<<endl;
  cout <<" FHC 1Re:  "<<n_fhc_e <<endl;
  cout <<" RHC 1Rmu: "<<n_rhc_mu<<endl;
  cout <<" RHC 1Re:  "<<n_rhc_e <<endl;
  cout <<" MultiRe:  "<<n_multie <<endl;
  
  fin->Close();
  
  TCanvas *c_1Rmu = new TCanvas("c_1Rmu","c_1Rmu");
  mock_data_1Rmu_0->Draw("colz");
  
  TCanvas *c_RHC1Rmu = new TCanvas("c_RHC1Rmu","c_RHC1Rmu");
  mock_data_RHC1Rmu_0->Draw("colz");
  
  TCanvas *c_1Re = new TCanvas("c_1Re","c_1Re");
  mock_data_1Re_0->Draw("colz");
  
  TCanvas *c_RHC1Re = new TCanvas("c_RHC1Re","c_RHC1Re");
  mock_data_RHC1Re_0->Draw("colz");
  
  TCanvas *c_MultiRe = new TCanvas("c_MultiRe","c_MultiRe");
  mock_data_MultiRe_0->Draw("colz");
  
  
  // Store histograms
  
  char name[200];
  sprintf(name,"%s.root",tag);
  TFile *fout = new TFile(name,"RECREATE");
  mock_data_1Rmu_0->Write();
  mock_data_RHC1Rmu_0->Write();
  mock_data_1Re_0->Write();
  mock_data_RHC1Re_0->Write();
  mock_data_MultiRe_0->Write();
  fout->Close();
  
  return;
}






double CalcThetaDegrees(float dir[10][3]){
  
  //double dir_nu[3] = {0.669764, -0.742179, 0.024223}; // T2K beam direction
  //double costheta = dir_nu[0] * dir[0][0] + dir_nu[1] * dir[0][1] + dir_nu[2] * dir[0][2];
  double costheta = kBeamDir[0] * dir[0][0] + kBeamDir[1] * dir[0][1] + kBeamDir[2] * dir[0][2];
  
  if( costheta<-1.00001 || costheta>+1.00001 ){
    cerr << endl;
    cerr << std::setprecision(10);
    cerr << "SuperkMC::FillTemplateHistogram" << endl;
    cerr << "costheta = " << costheta << " --> outside allowed range" << endl;
    cerr << endl;
    util::Die(1); return 0;
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



